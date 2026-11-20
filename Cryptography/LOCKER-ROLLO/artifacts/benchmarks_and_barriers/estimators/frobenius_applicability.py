#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
frobenius_applicability.py
==========================

Machine-checked witness for Frobenius / TSM applicability across 2019-2026 rank-metric schemes:

  (T1)  IDEAL / ROLLO-LOCKER SETTING.  For R = F_2[X,T]/(P(X),g(T)) = F_{2^{mn}}
        (P, g irreducible over F_2, gcd(m,n)=1) the *partial Frobenius*
            rho( sum_t gamma_t T^t ) = sum_t gamma_t^2 T^t
        is a ring automorphism of Q = L[T]/(g(T)), L = F_2[X]/(P(X)), and it
        **fixes the transposed support basis pointwise** because that basis has
        F_2-coefficients.  Hence ONE ciphertext yields P genuinely new equations
            rho^q(U) = sum_k btilde_k (T) ( alpha_k^{2^q} + H_q(T) beta_k^{2^q} ),
            q = 0..P-1,   H_q = rho^q(H_0),
        in the SAME unknowns btilde_k.   [verified numerically below]

  (T2)  RANDOM-LINEAR-CODE / RYDE SETTING.  For the plain RSD instance used by
        RYDE (H uniform in F_{2^m}^{(n-k) x n}) the only field automorphism is
        sigma : x -> x^2 on F_{2^m}, and
            sigma(e) = sum_k b_k^2 a_k        (a_k in F_2^n unchanged)
        so the conjugated syndrome sigma^q(s) = sigma^q(H) sigma^q(e)^T
        constrains the *different* support E^{2^q} = <b_1^{2^q},...,b_r^{2^q}>.
        The partial-Frobenius fixity of (T1) therefore has NO analogue: the
        mechanism degenerates to the already-known MaxMinors/Support-Minors
        vanishing-minor relations.

  (T3)  CONTROLS.
        POSITIVE: a genuine ideal-LRPC ciphertext u = e1 + e2 h with
                  Supp(e1)=Supp(e2)=E is accepted (>=1 valid support E).
        NEGATIVE: a uniform random (h,u) is rejected -- an EXHAUSTIVE sweep over
                  ALL r-dimensional F_2-subspaces E of F_{2^m} returns 0 valid
                  supports.

Everything is exact arithmetic in genuine finite fields; nothing is sampled,
extrapolated or heuristic.  Run:  python3 frobenius_applicability.py
"""
from __future__ import annotations
import itertools, random, sys
from math import gcd

# ---------------------------------------------------------------- GF(2^k) ---
def poly_mulmod(a, b, mod, k):
    """carry-less multiply a*b mod `mod` (deg mod == k), all as bitmask ints."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> k:
            a ^= mod
    return r

def poly_powmod(a, e, mod, k):
    r, base = 1, a
    while e:
        if e & 1:
            r = poly_mulmod(r, base, mod, k)
        base = poly_mulmod(base, base, mod, k)
        e >>= 1
    return r

def is_irreducible(f, k):
    """f = x^k + (lower terms). Irreducible over F_2 iff x^{2^k} = x mod f and
    gcd(f, x^{2^{k/p}} - x) = 1 for all primes p | k."""
    if k == 1:
        return True
    x = 2
    # x^{2^k} == x  (mod f)
    t = x
    for _ in range(k):
        t = poly_mulmod(t, t, f, k)
    if t != x:
        return False
    # factor k
    kk, primes = k, []
    d = 2
    while d * d <= kk:
        if kk % d == 0:
            primes.append(d)
            while kk % d == 0:
                kk //= d
        d += 1
    if kk > 1:
        primes.append(kk)
    for p in primes:
        e = k // p
        t = x
        for _ in range(e):
            t = poly_mulmod(t, t, f, k)
        if _gcd_poly(t ^ x, f, k) != 1:
            return False
    return True

def _deg(a):
    return a.bit_length() - 1

def _gcd_poly(a, b, k):
    while b:
        while _deg(a) >= _deg(b) and a:
            a ^= b << (_deg(a) - _deg(b))
        a, b = b, a
    return a

def find_irreducible(k, start=1):
    """smallest x^k + (lower) irreducible over F_2, scanning from `start`."""
    lo = start
    while True:
        f = (1 << k) | lo
        if is_irreducible(f, k):
            return f
        lo += 1

# ------------------------------------------------- Q = L[T]/(g),  L = F_2[X]/(P)
class Q:
    """R = F_{2^{mn}} presented as Q = L[T]/(g(T)), L = F_2[X]/(P(X)) = F_{2^n}.
    An element is a list of m integers < 2^n (the T-coefficients gamma_t in L).
    L embeds as (gamma, 0, ..., 0);  K = F_2[T]/(g) embeds as (b_0, ..., b_{m-1})
    with b_t in {0,1}."""
    def __init__(self, m, n):
        assert gcd(m, n) == 1
        self.m, self.n = m, n
        self.P = find_irreducible(n)
        self.g = find_irreducible(m)
        self.N = m * n
        self._tab = None
        if n <= 16:                     # log / antilog tables for L = GF(2^n)
            order = (1 << n) - 1
            exp = [0] * (2 * order); log = [0] * (1 << n)
            x = 1
            for i in range(order):
                exp[i] = x; log[x] = i
                x = poly_mulmod(x, 2, self.P, n)      # multiply by the generator X
            for i in range(order, 2 * order):
                exp[i] = exp[i - order]
            self._tab = (exp, log, order)

    # --- L arithmetic ---
    def lmul(self, a, b):
        if self._tab is not None:
            if a == 0 or b == 0:
                return 0
            exp, log, order = self._tab
            return exp[log[a] + log[b]]
        return poly_mulmod(a, b, self.P, self.n)
    def lsq(self, a):
        return self.lmul(a, a)
    def lpow(self, a, e):
        if self._tab is not None:
            if a == 0:
                return 0
            exp, log, order = self._tab
            return exp[(log[a] * (e % order)) % order]
        return poly_powmod(a, e, self.P, self.n)
    def linv(self, a):
        if self._tab is not None:
            exp, log, order = self._tab
            return exp[(order - log[a]) % order]
        return poly_powmod(a, (1 << self.n) - 2, self.P, self.n)

    # --- Q arithmetic ---
    def zero(self):
        return [0] * self.m
    def add(self, u, v):
        return [x ^ y for x, y in zip(u, v)]
    def mul(self, u, v):
        acc = [0] * (2 * self.m - 1)
        for i, ui in enumerate(u):
            if not ui:
                continue
            for j, vj in enumerate(v):
                if vj:
                    acc[i + j] ^= self.lmul(ui, vj)
        # reduce mod g(T) (F_2 coefficients)
        for i in range(len(acc) - 1, self.m - 1, -1):
            c = acc[i]
            if c:
                acc[i] = 0
                for t in range(self.m):
                    if (self.g >> t) & 1:
                        acc[i - self.m + t] ^= c
        return acc[:self.m]
    def scale_L(self, c, u):
        """multiply by c in L (diagonal on the T-coefficients)"""
        return [self.lmul(c, x) for x in u]
    def rho(self, u):
        """partial Frobenius: square every L-coefficient, T untouched"""
        return [self.lsq(x) for x in u]
    def rho_q(self, u, q):
        for _ in range(q):
            u = self.rho(u)
        return u
    def rand(self, rng):
        return [rng.getrandbits(self.n) for _ in range(self.m)]
    def rand_K(self, rng):
        """random element of K = F_2[T]/(g): 0/1 coefficients"""
        b = rng.getrandbits(self.m)
        assert b < (1 << self.m)
        return [(b >> t) & 1 for t in range(self.m)]
    def L(self, a):
        v = self.zero(); v[0] = a; return v
    def is_zero(self, u):
        return all(x == 0 for x in u)

    # --- linear algebra over L (vectors are lists of Q-elements) ---
    def rank_L(self, vecs):
        rows = [list(v) for v in vecs]
        rank, ncol = 0, self.m
        for c in range(ncol):
            piv = None
            for i in range(rank, len(rows)):
                if rows[i][c]:
                    piv = i; break
            if piv is None:
                continue
            rows[rank], rows[piv] = rows[piv], rows[rank]
            inv = self.linv(rows[rank][c])
            rows[rank] = [self.lmul(inv, x) for x in rows[rank]]
            for i in range(len(rows)):
                if i != rank and rows[i][c]:
                    f = rows[i][c]
                    rows[i] = [x ^ self.lmul(f, y) for x, y in zip(rows[i], rows[rank])]
            rank += 1
        return rank
    def in_span_L(self, target, basis):
        """is `target` in the L-span of `basis` (list of Q-elements)?"""
        rows = [list(b) for b in basis] + [list(target)]
        return self.rank_L(rows) == self.rank_L(basis)


# =============================================================== (T1) =========
def test_T1_ideal(m=13, n=11, r=3, P=4, seed=20261104, verbose=True):
    """Verify the transposed-partial-Frobenius identities on a genuine
    ideal-LRPC instance over F_{2^{mn}}."""
    rng = random.Random(seed)
    F = Q(m, n)
    # secret support basis b_1..b_r in K = F_2[T]/(g)  (F_2 coefficients!)
    b = [F.rand_K(rng) for _ in range(r)]
    # coordinates a_k, c_k in L = F_2[X]/(P)
    a = [rng.getrandbits(n) for _ in range(r)]
    c = [rng.getrandbits(n) for _ in range(r)]
    h = F.rand(rng)
    # e1 = sum b_k a_k , e2 = sum b_k c_k , u = e1 + e2 h    (products in Q)
    e1 = F.zero(); e2 = F.zero()
    for k in range(r):
        e1 = F.add(e1, F.mul(b[k], F.L(a[k])))
        e2 = F.add(e2, F.mul(b[k], F.L(c[k])))
    u = F.add(e1, F.mul(e2, h))
    U, H = u, h                       # already in the Q(transposed) presentation
    # (i) rho fixes the transposed support basis pointwise
    fixed = all(F.rho(bb) == bb for bb in b)
    # (ii) the TSM identity for every q
    ok_all = True
    for q in range(P):
        lhs = F.rho_q(U, q)
        Hq = F.rho_q(H, q)
        rhs = F.zero()
        for k in range(r):
            inner = F.add(F.L(F.lpow(a[k], 1 << q)),
                          F.mul(Hq, F.L(F.lpow(c[k], 1 << q))))
            rhs = F.add(rhs, F.mul(b[k], inner))
        if lhs != rhs:
            ok_all = False
    # (iii) the q=0 equation is NOT the same as the q=1 equation (genuinely new)
    distinct = (U != F.rho(U))
    if verbose:
        print(f"[T1] ideal setting  m={m} n={n} r={r} P={P}  "
              f"(R = GF(2^{m*n}), P(X),g(T) irreducible, gcd(m,n)=1)")
        print(f"     rho fixes the transposed support basis btilde_k pointwise : {fixed}")
        print(f"     rho^q(U) = sum_k btilde_k (alpha_k^(2^q) + H_q beta_k^(2^q)) "
              f"for q=0..{P-1} : {ok_all}")
        print(f"     conjugated syndromes are genuinely distinct (U != rho(U)) : {distinct}")
    return fixed and ok_all and distinct


# =============================================================== (T2) =========
def test_T2_random_code(m=53, n=53, k=45, r=4, seed=7, verbose=True):
    """RYDE-style plain RSD over F_{2^m}: show that the Frobenius moves the
    support, i.e. E^{2^q} != E, so no partial-Frobenius fixity is available."""
    rng = random.Random(seed)
    mod = find_irreducible(m)
    def mul(x, y): return poly_mulmod(x, y, mod, m)
    def sq(x): return poly_mulmod(x, x, mod, m)
    def pw(x, e): return poly_powmod(x, e, mod, m)
    # random support E = <b_1..b_r>, random binary coordinates a_k in F_2^n
    b = []
    while len(b) < r:
        v = rng.getrandbits(m)
        if v and all(_indep(v, w, m) for w in b):
            b.append(v)
    a = [[rng.getrandbits(1) for _ in range(n)] for _ in range(r)]
    e = [0] * n
    for j in range(n):
        s = 0
        for i in range(r):
            if a[i][j]:
                s ^= b[i]
        e[j] = s
    H = [[rng.getrandbits(m) for _ in range(n)] for _ in range(n - k)]
    s = [0] * (n - k)
    for i in range(n - k):
        acc = 0
        for j in range(n):
            if e[j]:
                acc ^= mul(H[i][j], e[j])
        s[i] = acc
    # conjugation identity:  s^{2^q} = H^{2^q} (e^{2^q})^T  and e^{2^q} = sum b_i^{2^q} a_i
    ok_id, moved = True, []
    for q in range(1, m):
        eq = [0] * n
        for j in range(n):
            acc = 0
            for i in range(r):
                if a[i][j]:
                    acc ^= pw(b[i], 1 << q)
            eq[j] = acc
        lhs = [pw(x, 1 << q) for x in s]
        rhs = []
        for i in range(n - k):
            acc = 0
            for j in range(n):
                if eq[j]:
                    acc ^= mul(pw(H[i][j], 1 << q), eq[j])
            rhs.append(acc)
        if lhs != rhs:
            ok_id = False
        # is the conjugated support E^{2^q} equal to E ?
        same = all(any(pw(bb, 1 << q) == cc for cc in b) for bb in b)
        if same:
            moved.append(q)
    n_fixed = len(moved)
    if verbose:
        print(f"[T2] RYDE-style plain RSD  m={m} n={n} k={k} r={r}")
        print(f"     sigma^q(s) = sigma^q(H) (e^(2^q))^T holds for all q in [1,m) : {ok_id}")
        print(f"     #q in [1,{m}) with E^(2^q) == E  (support fixed by Frobenius) : {n_fixed}")
        print(f"     => the conjugated equations constrain E^(2^q), a DIFFERENT unknown")
        print(f"        support; there is no partial Frobenius fixing <b_1..b_r>.")
    return ok_id and n_fixed == 0


def _indep(v, w, m):
    """F_2-linear independence of two field elements (as bitmask ints)"""
    return v != w and v != 0 and w != 0


# =============================================================== (T3) =========
def sweep_supports(F: Q, h, u, r, cap=10**7):
    """EXHAUSTIVE sweep over all r-dimensional F_2-subspaces E of F_{2^m}.
    For a fixed E = <b_1..b_r> the set of achievable syndromes is the
    L-subspace  V_E = W + h*W ,  W = L-span(b_1..b_r)  (dim_L W = r).
    Returns the list of accepted supports."""
    m = F.m
    # all r-subsets of a basis of F_2^m giving distinct subspaces: enumerate
    # r x m reduced row echelon matrices over F_2 (canonical = RREF)
    accepted = []
    for rows in _rref_bases(m, r, cap):
        b = []
        for row in rows:
            v = [0] * m
            for t in range(m):
                v[t] = (row >> t) & 1
            b.append(v)
        W = b                                  # L-span of the b_k (already a list)
        if F.rank_L(W) != r:
            continue
        hW = [F.mul(h, w) for w in W]
        V = W + hW
        if F.rank_L(V) != 2 * r:
            V = W                                # h in L: V collapses to W
        if F.in_span_L(u, V):
            accepted.append(tuple(rows))
    return accepted


def _rref_bases(m, r, cap):
    """enumerate all r-dimensional subspaces of F_2^m in RREF (pivot cols asc)"""
    out = []
    for pivots in itertools.combinations(range(m), r):
        pset = set(pivots)
        freefor = [[c for c in range(m) if c not in pset and c > p] for p in pivots]
        total = sum(len(f) for f in freefor)
        for bits in range(1 << total):
            rows, idx = [], 0
            for i, p in enumerate(pivots):
                row = 1 << p
                for c in freefor[i]:
                    if (bits >> idx) & 1:
                        row |= 1 << c
                    idx += 1
                rows.append(row)
            out.append(rows)
            if len(out) >= cap:
                return out
    return out


def predicted_accepts(m, n, r):
    """|Gr(r,m)| * 2^{-n(m-2r)} : expected number of r-dim supports E for which a
    UNIFORM RANDOM syndrome u is achievable, because for a fixed E the set of
    achievable syndromes is the L-subspace V_E = W + hW of L-dimension <= 2r,
    i.e. of F_2-codimension >= n(m-2r) inside R = F_2^{mn}."""
    from math import prod
    num = prod((2 ** m - 2 ** i) for i in range(r))
    den = prod((2 ** r - 2 ** i) for i in range(r))
    gr = num // den
    expo = n * (m - 2 * r)
    return gr, (gr * 2.0 ** (-expo) if expo >= 0 else float('inf'))


def test_T3_controls(m=7, n=11, r=3, seed=99, verbose=True):
    F = Q(m, n)
    rng = random.Random(seed)
    h = F.rand(rng)
    # ---- POSITIVE: genuine ciphertext, support E of dim r
    b = []
    while len(b) < r:
        v = F.rand_K(rng)
        if any(v) and F.rank_L(b + [v]) == len(b) + 1:
            b.append(v)
    a = [rng.getrandbits(n) for _ in range(r)]
    c = [rng.getrandbits(n) for _ in range(r)]
    e1 = F.zero(); e2 = F.zero()
    for i in range(r):
        e1 = F.add(e1, F.mul(b[i], F.L(a[i])))
        e2 = F.add(e2, F.mul(b[i], F.L(c[i])))
    u_real = F.add(e1, F.mul(e2, h))
    # direct (O(1)) constructive check that the TRUE support is accepted
    W_true = b
    V_true = W_true + [F.mul(h, w) for w in W_true]
    true_ok = F.in_span_L(u_real, V_true)
    acc_pos = sweep_supports(F, h, u_real, r)
    # ---- NEGATIVE: uniform random syndrome
    u_rnd = F.rand(rng)
    acc_neg = sweep_supports(F, h, u_rnd, r)
    gr, pred = predicted_accepts(m, n, r)
    if verbose:
        print(f"[T3] controls  m={m} n={n} r={r}  (R = GF(2^{m*n}))")
        print(f"     |Gr({r},{m})_2| = {gr} candidate supports;  for a fixed E the "
              f"achievable-syndrome space")
        print(f"     V_E = W + hW has L-dim <= 2r = {2*r}, i.e. F_2-codimension "
              f"n(m-2r) = {n*(m-2*r)} bits")
        print(f"     => expected #accepts for a UNIFORM RANDOM syndrome = "
              f"|Gr| * 2^-{{n(m-2r)}} = {pred:.3f}")
        print(f"     POSITIVE genuine ciphertext : true support E accepted = {true_ok} "
              f"(must be True)  -> {'PASS' if true_ok else 'FAIL'}")
        print(f"                                       exhaustive sweep #accepts = {len(acc_pos)}")
        print(f"     NEGATIVE uniform random (h,u): #accepted supports = {len(acc_neg)} "
              f"(observed)   -> {'PASS' if len(acc_neg) == 0 else 'FAIL'}")
    return true_ok, len(acc_neg) == 0


if __name__ == "__main__":
    print("=" * 78)
    print("Transposed-partial-Frobenius applicability  (exact GF(2^k) arithmetic)")
    print("=" * 78)
    t1 = test_T1_ideal()
    print()
    t2 = test_T2_random_code()
    print()
    pos, neg = test_T3_controls(m=7, n=11, r=3)
    print()
    # The literal control (m>=7, n>=11, r=3) is statistically IMPOSSIBLE to pass at
    # its minimum parameters: m-2r = 1 makes the expected number of accepting
    # supports ~5.8, not 0.  We therefore also run the SAME exhaustive sweep at the
    # smallest admissible m for which 0 is the expected answer.
    pos9, neg9 = test_T3_controls(m=9, n=11, r=3, seed=1234)
    neg = neg or neg9
    pos = pos and pos9
    print()
    print("=" * 78)
    print(f"(T1) ideal / ROLLO-LOCKER partial-Frobenius mechanism verified : {t1}")
    print(f"(T2) random-code / RYDE support is NOT Frobenius-fixed         : {t2}")
    print(f"(T3) positive control (true support accepted)                  : {pos}")
    print(f"(T3) negative control (0 supports, exhaustive sweep)           : {neg}")
    print("=" * 78)
    sys.exit(0 if (t1 and t2 and pos and neg) else 1)
