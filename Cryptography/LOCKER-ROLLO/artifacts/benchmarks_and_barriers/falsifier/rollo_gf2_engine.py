#!/usr/bin/env python3
"""
rollo_gf2_engine.py -- EXACT GF(2) engine for the ROLLO/LOCKER Transposed-Syndrome-Map (TSM)
system, over TRUE GF(2^n)/GF(2^m)/GF(2^(mn)) with the REAL partial Frobenius rho.

Mathematical model (verified against parameterized_tsm/tsm.cpp):
  K = GF(2)[T]/(g) = GF(2^m),  L = GF(2)[X]/(P) = GF(2^n),  gcd(m,n)=1
  R = K[X]/P = GF(2)[T,X]/(g,P) =: Q = L[T]/g          (the SAME ring, two bases)
  Phi : R -> Q is the n x m bit-matrix transpose (a ring isomorphism).
  rho : Q -> Q,  rho(sum_t gamma_t T^t) = sum_t gamma_t^2 T^t
        (= Frobenius of L extended by fixing T; it FIXES K pointwise).
  In the R-basis rho acts as  X -> X^2, T -> T.

  ROLLO-I ciphertext:  h = x^{-1} y,  u = e1 + e2 h,  Supp(e1),Supp(e2) subset E,
  dim_{GF(2)} E = r.  Writing e_i = sum_k b_k a_k with a_k in L (F_2-coefficients!),
      Phi(u) = sum_k b~_k (alpha_k + H_0 beta_k),   alpha_k,beta_k in L.
  Applying rho^q:
      U_q = sum_k b~_k ( rho^q(alpha_k) + H_q rho^q(beta_k) ),  H_q = rho^q(H_0), U_q = rho^q(U_0).

  ==> THE P CONJUGATE EQUATIONS ARE rho-IMAGES OF THE q=0 EQUATION.  <==  (Lemma REDUND)
"""
import itertools, random
from math import gcd

# ---------------------------------------------------------------- GF(2)[x] polys as ints
def pdeg(a): return a.bit_length() - 1

def pmul(a, b):
    r = 0
    while b:
        lb = b & -b
        r ^= a << (lb.bit_length() - 1)
        b ^= lb
    return r

def pmod(a, m):
    dm = pdeg(m)
    while a.bit_length() - 1 >= dm and a:
        a ^= m << (pdeg(a) - dm)
    return a

def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a

def is_irred(f):
    d = pdeg(f)
    if d <= 0: return False
    if d == 1: return True
    # x^(2^d) == x (mod f) and gcd(x^(2^(d/p)) - x, f)=1 for all primes p | d
    x = 2
    for p in _primes(d):
        y = x
        for _ in range(d // p):
            y = pmod(pmul(y, y), f)
        if pgcd(y ^ 2, f) != 1: return False
    y = x
    for _ in range(d):
        y = pmod(pmul(y, y), f)
    return y == 2

def _primes(n):
    out, d, m = [], 2, n
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: out.append(m)
    return out

def irred_poly(d, seed=0):
    """Smallest irreducible poly of degree d (with the leading and constant terms)."""
    for c in range(1, 1 << (d - 1)):
        f = (1 << d) | (c << 1) | 1
        if is_irred(f): return f
    raise RuntimeError

# ---------------------------------------------------------------- small fields
class F2k:
    """GF(2^k) as ints, modulus `mod` (bit-encoded irreducible poly of degree k)."""
    def __init__(self, k, mod=None):
        self.k = k
        self.mod = mod if mod is not None else irred_poly(k)
        assert pdeg(self.mod) == k
        self.mask = (1 << k) - 1
        # log / antilog tables w.r.t. a generator (search for one)
        self._build_tables()

    def mul(self, a, b):
        r = 0
        while b:
            lb = b & -b
            r ^= a << (lb.bit_length() - 1)
            b ^= lb
        return pmod(r, self.mod)

    def sq(self, a): return self.mul(a, a)

    def _build_tables(self):
        n = (1 << self.k) - 1
        for g in range(2, 1 << self.k):
            seen, x, ok = [0] * (n + 1), 1, True
            for i in range(n):
                if seen[x]: ok = False; break
                seen[x] = 1
                x = self.mul(x, g)
            if ok and x == 1:
                self.gen = g
                self.expg = [0] * (2 * n)
                self.logg = [0] * (1 << self.k)
                x = 1
                for i in range(n):
                    self.expg[i] = x
                    self.expg[i + n] = x
                    self.logg[x] = i
                    x = self.mul(x, g)
                return
        raise RuntimeError("no generator")

    def fmul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.expg[self.logg[a] + self.logg[b]]

    def finv(self, a):
        assert a
        n = (1 << self.k) - 1
        return self.expg[(n - self.logg[a]) % n]

    def fdiv(self, a, b):
        if a == 0: return 0
        n = (1 << self.k) - 1
        return self.expg[(self.logg[a] - self.logg[b]) % n]

    def fpow(self, a, e):
        r, b = 1, a
        while e:
            if e & 1: r = self.fmul(r, b)
            b = self.fmul(b, b)
            e >>= 1
        return r

# ---------------------------------------------------------------- Q = L[T]/g  (elements: list of m L-elems)
class RingQ:
    """Q = L[T]/(g(T)); elements are tuples of length m with entries in L (ints < 2^n).
       Also provides the R-view (n coefficients in K) and the transpose Phi."""
    def __init__(self, m, n, gmod=None, pmod=None):
        assert gcd(m, n) == 1, "need gcd(m,n)=1 so that R = Q = GF(2^(mn))"
        self.m, self.n = m, n
        self.K = F2k(m, gmod)
        self.L = F2k(n, pmod)

    # ---- Q arithmetic (list of m L-coeffs, index = T-degree)
    def qzero(self): return [0] * self.m
    def qadd(self, a, b): return [x ^ y for x, y in zip(a, b)]
    def qmul(self, a, b):
        m, L, gm = self.m, self.L, self.K.mod
        tmp = [0] * (2 * m - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj: tmp[i + j] ^= L.fmul(ai, bj)
        for d in range(2 * m - 2, m - 1, -1):
            c = tmp[d]
            if c:
                sh = d - m
                for t in range(m):
                    if (gm >> t) & 1:
                        tmp[sh + t] ^= c
        return tmp[:m]

    def rho(self, a):  # coefficientwise squaring in L
        L = self.L
        return [L.fmul(x, x) for x in a]

    def K_to_Q(self, z):
        """Embed K (int, m bits) into Q: z = sum_t z_t T^t, z_t in {0,1} subset L."""
        return [(z >> t) & 1 for t in range(self.m)]

    # ---- R arithmetic (list of n K-coeffs, index = X-degree)
    def rmul(self, a, b):
        n, K = self.n, self.K
        tmp = [0] * (2 * n - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj: tmp[i + j] ^= K.mul(ai, bj)
        pm = self.L.mod
        for d in range(2 * n - 2, n - 1, -1):
            c = tmp[d]
            if c:
                sh = d - n
                for t in range(n):
                    if (pm >> t) & 1: tmp[sh + t] ^= c
        return tmp[:n]

    def rinv(self, a):
        """inverse in R = GF(2^(mn)) by Fermat: a^(2^(mn)-2)."""
        e = (1 << (self.m * self.n)) - 2
        r, b = [1] + [0] * (self.n - 1), a
        while e:
            if e & 1: r = self.rmul(r, b)
            b = self.rmul(b, b)
            e >>= 1
        return r

    def rho_R(self, a):
        """rho in the R-basis: X -> X^2, coefficients (in K) fixed."""
        out = [0] * self.n
        for j, aj in enumerate(a):
            if aj:
                d = (2 * j) % self.n
                # X^(2j) mod P : need reduction; do it by repeated squaring of X
                out = [x ^ y for x, y in zip(out, self._Xpow(2 * j, aj))]
        return out

    def _Xpow_cache(self):
        if not hasattr(self, "_xpow"):
            self._xpow = {}
        return self._xpow

    def _Xpow(self, e, coef):
        """coef * X^e  in R."""
        cache = self._Xpow_cache()
        if e not in cache:
            base = [0] * self.n; base[1 % self.n] = 1
            r = [1] + [0] * (self.n - 1)
            b, ee = base, e
            while ee:
                if ee & 1: r = self.rmul(r, b)
                b = self.rmul(b, b)
                ee >>= 1
            cache[e] = r
        v = cache[e]
        K = self.K
        return [K.mul(coef, x) for x in v]

    def Phi(self, r):
        """R (n K-coeffs) -> Q (m L-coeffs): bit-matrix transpose."""
        out = [0] * self.m
        for j, rj in enumerate(r):
            xj = 1 << j          # X^j in L
            t = 0
            while rj:
                if rj & 1: out[t] ^= xj
                rj >>= 1; t += 1
        return out

    def Phi_inv(self, q):
        out = [0] * self.n
        for t, qt in enumerate(q):
            j = 0
            while qt:
                if qt & 1: out[j] ^= (1 << t)
                qt >>= 1; j += 1
        return out

# ---------------------------------------------------------------- linear algebra over L
def L_rank(rows, L, ncols):
    """rank over L of a list of rows (each a list of L-elements)."""
    rows = [r[:] for r in rows]
    rk, piv = 0, []
    for c in range(ncols):
        p = -1
        for i in range(rk, len(rows)):
            if rows[i][c]: p = i; break
        if p < 0: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        iv = L.finv(rows[rk][c])
        rows[rk] = [L.fmul(x, iv) for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]
                rows[i] = [x ^ L.fmul(f, y) for x, y in zip(rows[i], rows[rk])]
        piv.append(c); rk += 1
        if rk == len(rows): break
    return rk

def L_in_span(vecs, target, L, ncols):
    """is `target` in the L-span of `vecs` (all rows of length ncols)?"""
    return L_rank(list(vecs) + [target], L, ncols) == L_rank(list(vecs), L, ncols)

# ---------------------------------------------------------------- GF(2) bit-matrix rank (int rows)
def gf2_rank(rows, ncols):
    rows = [r for r in rows]
    rk = 0
    for c in range(ncols):
        p = -1
        for i in range(rk, len(rows)):
            if (rows[i] >> c) & 1: p = i; break
        if p < 0: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        pr = rows[rk]
        for i in range(len(rows)):
            if i != rk and ((rows[i] >> c) & 1): rows[i] ^= pr
        rk += 1
        if rk == len(rows): break
    return rk

# ---------------------------------------------------------------- subspaces of GF(2)^m
def rref_subspaces(m, r):
    """Yield every r-dimensional GF(2) subspace of GF(2)^m as a basis of r ints (RREF)."""
    for pivots in itertools.combinations(range(m), r):
        ps = set(pivots)
        nonpiv = [c for c in range(m) if c not in ps]
        free = [[c for c in nonpiv if c > p] for p in pivots]
        for bits in itertools.product(*[range(1 << len(f)) for f in free]):
            rows = []
            for i, p in enumerate(pivots):
                v = 1 << p
                for k, c in enumerate(free[i]):
                    if (bits[i] >> k) & 1: v |= 1 << c
                rows.append(v)
            yield rows

def count_subspaces(m, r):
    num, den = 1, 1
    for i in range(r):
        num *= (1 << (m - i)) - 1
        den *= (1 << (r - i)) - 1
    return num // den
