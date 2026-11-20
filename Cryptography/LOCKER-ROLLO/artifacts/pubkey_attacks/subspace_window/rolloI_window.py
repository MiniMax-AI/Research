#!/usr/bin/env python3
"""
rolloI_window.py -- Exact arithmetic for the 2020-04-21 ROLLO-I specification and
the "Window Lemma" key-recovery attack on known-window weak keys.

Spec-authentic polynomials (ROLLO spec v2.0, 2020-04-21, Table 2 p.18 and Table 3 p.20):
  ROLLO-I-128 : n=83, m=67, r=7, d=8 ; P = X^83+X^7+X^4+X^2+1 ; g = X^67+X^5+X^2+X+1
  ROLLO-I-192 : n=97, m=79, r=8, d=8 ; P = X^97+X^6+1          ; g = X^79+X^9+1
  ROLLO-I-256 : n=113,m=97, r=9, d=9 ; P = X^113+X^9+1         ; g = X^97+X^6+1

Model:  K = GF(2^m) = GF(2)[T]/(g(T)),  L = GF(2^n) = GF(2)[X]/(P(X)),
        R = K[X]/(P(X)) = GF(2^(mn))   (gcd(m,n)=1 so P stays irreducible over K),
        R ~= K (x)_{F2} L as F2-algebras.
  x in R is stored as [x_0,...,x_{n-1}], x_i in K;  Supp_{F2}(x) = {x_i != 0};
  Sigma(x) = span_{F2} Supp(x) subseteq K;  rank(x) = dim Sigma(x).

CENTRAL IDENTITY (Window Lemma).  Supp(x) subseteq F  <=>  x in V_F := F (x)_{F2} L,
and V_F is an L-SUBSPACE of R of L-dimension dim_{F2} F.  Multiplication by h is
L-linear.  Hence for any D-dim F2-subspace W of K,

    { x in R : Supp(x) subseteq W and Supp(xh) subseteq W } = ker M_W,
    M_W in Mat_{m-D, D}(L),   computable in O((m-D) D^2 n^2) bit operations.

Moreover this kernel is SCALE INVARIANT:  ker M_{gamma W} = gamma . ker M_W for
gamma in K^* (because R is a field, gamma commutes with h).  So ONE window W
covers the whole family {gamma W : gamma in K^*}.
"""
import random

# ============================================================ GF(2^k) ========
class GF2k:
    """GF(2^k) = GF(2)[z]/(z^k + low), elements = python ints, bit j <-> z^j."""
    __slots__ = ("k", "low", "mod", "mask", "_red")

    def __init__(self, k, low):
        self.k = k
        self.low = low                 # low has bit j set iff z^j appears, j<k
        self.mod = (1 << k) | low
        self.mask = (1 << k) - 1
        # reduction table: for the top byte positions
        self._red = None

    def mul(self, a, b):
        r = 0
        k, low, mask = self.k, self.low, self.mask
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> k:
                a = (a & mask) ^ low
        return r

    def sq(self, a):
        return self.mul(a, a)

    def pow(self, a, e):
        r, base = 1, a
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.mul(base, base)
            e >>= 1
        return r

    def inv(self, a):
        assert a, "inverse of 0"
        # extended Euclid on GF(2)[z]
        u, v = a, self.mod
        g1, g2 = 1, 0
        while u != 1:
            j = u.bit_length() - v.bit_length()
            if j < 0:
                u, v, g1, g2, j = v, u, g2, g1, -j
            u ^= v << j
            g1 ^= g2 << j
        k, low, mask = self.k, self.low, self.mask
        while g1 >> k:
            hi = g1 >> k
            g1 = (g1 & mask) ^ self._mul_low(hi, low, k)
        return g1

    @staticmethod
    def _mul_low(a, b, k):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
        return r


def gf2_poly_mulmod(a, b, mod, k):
    r = 0
    mask = (1 << k) - 1
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> k:
            a = (a & mask) ^ (mod & mask)
    return r


def is_irreducible(k, low):
    mod = (1 << k) | low
    if not (low & 1):
        return False
    if gf2_poly_mulmod(2, 2, mod, k) is None:
        return False
    # x^(2^k) == x (mod mod)
    v = 2
    for _ in range(k):
        v = gf2_poly_mulmod(v, v, mod, k)
    if v != 2:
        return False
    # factor k
    mm, fac = k, []
    f = 2
    while f * f <= mm:
        while mm % f == 0:
            fac.append(f); mm //= f
        f += 1
    if mm > 1:
        fac.append(mm)
    for q in set(fac):
        w = 2
        for _ in range(k // q):
            w = gf2_poly_mulmod(w, w, mod, k)
        g = gf2_gcd(w ^ 2, mod)
        if g != 1:
            return False
    return True


def gf2_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length() and a:
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def find_irreducible(k, rng):
    for _ in range(200000):
        low = rng.randrange(1, 1 << k) | 1
        if is_irreducible(k, low):
            return low
    raise RuntimeError("no irreducible polynomial found for k=%d" % k)


# ===================================================== R = K[X]/(P) =========
class RolloI:
    """ROLLO-I instance. R-elements are python lists of n GF(2^m) ints."""

    SPEC = {
        #  name              m   n   r   d   g_low(deg m field poly)  P_low(deg n modulus)
        "ROLLO-I-128": (67, 83, 7, 8, (1 << 5) | (1 << 2) | (1 << 1) | 1,
                        (1 << 7) | (1 << 4) | (1 << 2) | 1),
        "ROLLO-I-192": (79, 97, 8, 8, (1 << 9) | 1, (1 << 6) | 1),
        "ROLLO-I-256": (97, 113, 9, 9, (1 << 6) | 1, (1 << 9) | 1),
    }

    def __init__(self, m=None, n=None, d=None, r=None, g_low=None, P_low=None,
                 seed=0, name=None):
        if name is not None:
            m, n, r, d, g_low, P_low = self.SPEC[name]
        rng = random.Random(seed)
        self.name = name or "custom"
        self.m, self.n, self.d, self.r = m, n, d, r
        self.rng = rng
        self.K = GF2k(m, g_low if g_low is not None else find_irreducible(m, rng))
        self.P_low = P_low if P_low is not None else find_irreducible(n, rng)
        self.L = GF2k(n, self.P_low)
        self.Pbits = [j for j in range(n) if (self.P_low >> j) & 1]

    # ---- R arithmetic ----------------------------------------------------
    def rmul(self, a, b):
        n, K = self.n, self.K
        tmp = [0] * (2 * n - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        tmp[i + j] ^= K.mul(ai, bj)
        for deg in range(2 * n - 2, n - 1, -1):
            c = tmp[deg]
            if c:
                sh = deg - n
                for j in self.Pbits:
                    tmp[sh + j] ^= c
        return tmp[:n]

    def rscal_L(self, lam, a):
        """multiply R-element a by lambda in L, where L is embedded in R as
        {0,1}-coefficient X-polynomials:  lam * a = rmul(lam_as_R_element, a)."""
        n = self.n
        lamR = [1 if (lam >> i) & 1 else 0 for i in range(n)]
        return self.rmul(lamR, a)

    def rinv(self, a):
        K = self.K
        def trim(p):
            while len(p) > 1 and p[-1] == 0:
                p.pop()
            return p
        def pmul(p, q):
            out = [0] * (len(p) + len(q) - 1)
            for i, ai in enumerate(p):
                if ai:
                    for j, bj in enumerate(q):
                        if bj:
                            out[i + j] ^= K.mul(ai, bj)
            return trim(out)
        def padd(p, q):
            return trim([(p[i] if i < len(p) else 0) ^ (q[i] if i < len(q) else 0)
                         for i in range(max(len(p), len(q)))])
        def pdivmod(num, den):
            num, den = trim(num[:]), trim(den[:])
            if len(num) < len(den):
                return [0], num
            q = [0] * (len(num) - len(den) + 1)
            il = K.inv(den[-1])
            while len(num) >= len(den) and num != [0]:
                sh = len(num) - len(den)
                sc = K.mul(num[-1], il)
                q[sh] ^= sc
                for i, c in enumerate(den):
                    if c:
                        num[i + sh] ^= K.mul(sc, c)
                num = trim(num)
            return trim(q), num
        mod = [0] * (self.n + 1)
        mod[self.n] = 1
        for j in self.Pbits:
            mod[j] = 1
        old_r, rr = trim(mod[:]), trim(a[:])
        old_t, t = [0], [1]
        while rr != [0]:
            q, rem = pdivmod(old_r, rr)
            old_r, rr = rr, rem
            old_t, t = t, padd(old_t, pmul(q, t))
        if len(old_r) != 1 or old_r[0] == 0:
            return None
        iv = K.inv(old_r[0])
        out = [K.mul(c, iv) for c in old_t] + [0] * self.n
        return out[:self.n]

    # ---- sampling --------------------------------------------------------
    def rand_subspace(self, dim):
        m = self.m
        while True:
            bas = gf2_rref([self.rng.randrange(1, 1 << m) for _ in range(dim)], m)
            if len(bas) == dim:
                return bas

    def rand_dense(self, basis):
        dim, n = len(basis), self.n
        while True:
            v = [0] * n
            for i in range(n):
                msk = self.rng.randrange(1 << dim)
                acc = 0
                for k in range(dim):
                    if (msk >> k) & 1:
                        acc ^= basis[k]
                v[i] = acc
            if gf2_rank(v, self.m) == dim:
                return v

    def keygen(self, F=None):
        """ROLLO-I KeyGen (spec Fig.2): F full-rank dim d, (x,y) in S_d^{2n}, h = x^-1 y."""
        F = F if F is not None else self.rand_subspace(self.d)
        while True:
            x = self.rand_dense(F)
            xi = self.rinv(x)
            if xi is not None:
                y = self.rand_dense(F)
                if gf2_rank(x + y, self.m) == self.d:      # dim Supp(x,y) = d exactly
                    return F, x, y, self.rmul(xi, y)


# ================================================= GF(2) linear algebra =====
def gf2_rref(rows, nbits):
    rows = [r for r in rows]
    rk = 0
    for col in range(nbits):
        piv = -1
        for i in range(rk, len(rows)):
            if (rows[i] >> col) & 1:
                piv = i
                break
        if piv < 0:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rk]
        rk += 1
    return rows[:rk]


def gf2_rank(rows, nbits):
    return len(gf2_rref(rows, nbits))


def gf2_mat_inverse(rows, n):
    A = list(rows)
    Iv = [1 << i for i in range(n)]
    for col in range(n):
        piv = -1
        for i in range(col, n):
            if (A[i] >> col) & 1:
                piv = i
                break
        if piv < 0:
            raise ValueError("singular")
        A[col], A[piv] = A[piv], A[col]
        Iv[col], Iv[piv] = Iv[piv], Iv[col]
        for i in range(n):
            if i != col and ((A[i] >> col) & 1):
                A[i] ^= A[col]
                Iv[i] ^= Iv[col]
    return Iv


def parity(x):
    return bin(x).count("1") & 1


def dual_basis(wbas, m):
    """Given an F2-basis wbas[0..m-1] of GF(2)^m (ints, bit j = coefficient of T^j),
    return the dual basis: dual[k] satisfies  parity(dual[k] & wbas[l]) = delta_{kl}.
    Then coord_k(v) = parity(dual[k] & v)."""
    inv = gf2_mat_inverse(wbas, m)          # rows: A^{-1}, A[k] = wbas[k]
    dual = [0] * m                          # dual[k] = column k of A^{-1}
    for r in range(m):
        ir = inv[r]
        for k in range(m):
            if (ir >> k) & 1:
                dual[k] |= 1 << r
    return dual


# ================================================== the Window Lemma ========
def extend_basis(W, m):
    """Extend independent K-elements W to a full F2-basis of K=GF(2)^m,
    keeping the first len(W) entries EXACTLY equal to gf2_rref(W) (a basis of W)."""
    bas = gf2_rref(W[:], m)
    assert len(bas) == len(W), "W is not linearly independent"
    red = list(bas)
    for j in range(m):
        if len(red) == m:
            break
        v = 1 << j
        for b in red:
            lb = b & -b
            if v & lb:
                v ^= b
        if v:
            red.append(v)
    assert len(red) == m, "could not extend to a full basis"
    return red


def window_matrix(inst, h, W):
    """M_W in Mat_{m-D, D}(L):
         ker M_W = { lam in L^D : x = sum_j lam_j w_j has Supp(x), Supp(xh) subseteq W }.
       Entry (k,j) = sum_{i<n} [ w_{D+k}-coordinate of (w_j * h_i) in K ] * X^i  in L.
       (Note w_j*h needs NO reduction mod P: deg_X < n.)"""
    m, n, K = inst.m, inst.n, inst.K
    D = len(W)
    wbas = extend_basis(W, m)
    dual = dual_basis(wbas, m)          # coord_k(v) = parity(dual[k] & v)
    M = [[0] * D for _ in range(m - D)]
    for j in range(D):
        wj = wbas[j]
        for i in range(n):
            if h[i] == 0:
                continue
            prod = K.mul(wj, h[i])
            for k in range(m - D):
                if parity(prod & dual[D + k]):
                    M[k][j] ^= (1 << i)
    return M, wbas


# ============================================ GF(2^n) linear algebra ========
def L_reduce(M, nrows, ncols, L):
    A = [row[:] for row in M]
    piv = []
    r = 0
    for c in range(ncols):
        sel = -1
        for i in range(r, nrows):
            if A[i][c]:
                sel = i
                break
        if sel < 0:
            continue
        A[r], A[sel] = A[sel], A[r]
        iv = L.inv(A[r][c])
        A[r] = [L.mul(iv, v) for v in A[r]]
        for i in range(nrows):
            if i != r and A[i][c]:
                f = A[i][c]
                Ai = A[i]
                Ar = A[r]
                for k in range(ncols):
                    if Ar[k]:
                        Ai[k] ^= L.mul(f, Ar[k])
        piv.append(c)
        r += 1
        if r == nrows:
            break
    return A, piv, r


def L_rank(M, nrows, ncols, L):
    return L_reduce(M, nrows, ncols, L)[2]


def L_nullspace(M, nrows, ncols, L):
    A, piv, r = L_reduce(M, nrows, ncols, L)
    free = [c for c in range(ncols) if c not in piv]
    ker = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for i in range(len(piv) - 1, -1, -1):
            s = 0
            for k in range(piv[i] + 1, ncols):
                if A[i][k] and v[k]:
                    s ^= L.mul(A[i][k], v[k])
            v[piv[i]] = s
        ker.append(v)
    return ker


def kernel_to_R(inst, v, wbas):
    """lam in L^D  ->  x = sum_j lam_j w_j  in R."""
    n, D = inst.n, len(v)
    x = [0] * n
    for j in range(D):
        lam = v[j]
        if not lam:
            continue
        for i in range(n):
            if (lam >> i) & 1:
                x[i] ^= wbas[j]
    return x
