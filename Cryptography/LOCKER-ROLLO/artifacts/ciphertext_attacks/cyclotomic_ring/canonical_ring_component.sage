#!/usr/bin/env sage
# canonical_ring_component.sage
#
# Portability of the partial-Frobenius / MaxMinors support-recovery attack to
# the canonical ROLLO/LOCKER ring modulus Z^n - 1.
#
# Running the linearisation directly over the whole quotient ring
# R = F[Z]/(Z^n - 1) produces a Frobenius system with trivial kernel; this
# script isolates why and verifies the cyclotomic-component projection fix.
#
# Structure.  For n prime and 2 a primitive root mod n, over GF(2)
#       Z^n - 1 = (Z+1) * Phi_n ,   Phi_n irreducible of degree n-1,
# so    R = F[Z]/(Z^n-1)  ~=  F  x  F[Z]/(Phi_n)   (a product, NOT a field).
# The partial Frobenius rho = id_F (x) sigma acts as the IDENTITY on the first
# component (F2 has trivial Frobenius) and as the F-Frobenius on the second.
#
# Consequence for the MaxMinors linearisation:
#   * component 1 (Z -> 1): every row of the MinRank matrix is identical, so
#     every (r+1)x(r+1) minor vanishes IDENTICALLY -- the component carries no
#     equation at all, and the "Frobenius" relation p(S+1,T+1) = p(S,T)^2
#     degenerates to the LINEAR relation p(S+1,T+1) = p(S,T).  That is exactly
#     the observed failure.
#   * component 2 (mod Phi_n): a genuine MinRank instance over the field
#     F (x) GF(2^{n-1}), identical in shape to the irreducible-P build.
#
# Fix under test: run the linear algebra in component 2 only, i.e. reduce every
# ring element modulo Phi_n first.  Then read the support E off the recovered
# e2's F-coordinates.  Note pi(e2) has F-coordinates e2_i + e2_{n-1}, which all
# lie in E and generically still span it -- so E is recovered with no ambiguity.
#
# Run:  sage canonical_ring_component.sage

from itertools import combinations
import sys

# ------------------------------------------------------------------ parameters
# small but faithful: n prime, 2 primitive mod n (ord_11(2) = 10), Phi_11 irreducible
n = 11
r = 3
p = 6
m = 31                      # gcd(31,10)=1 so F (x) GF(2^10) is a field
assert (2 ** 10) % n == 1 and (2 ** 5) % n == n - 1     # ord_n(2) = n-1

def log(*a):
    print(*a); sys.stdout.flush()

# ------------------------------------------------------------------- the rings
ZF = PolynomialRing(GF(2), 'Z')
Zv = ZF.gen()
Zn1 = Zv ** n + 1                                   # Z^n - 1  (= Z^n + 1 in char 2)
Phi = sum(Zv ** i for i in range(n))                # (Z^n - 1)/(Z+1) = 1+Z+...+Z^{n-1}
assert Zn1 == (Zv + 1) * Phi
assert Phi.is_irreducible()

F = GF(2 ** m, 'T')                                 # support field
Lp = GF(2 ** (n - 1), 'z', modulus=Phi)             # component-2 base field
z = Lp.gen()

S = PolynomialRing(Lp, 'X')
fq = S([Lp(c) for c in F.modulus().list()])         # same modulus as F, over Lp
assert fq.is_irreducible()
Q = S.quotient(fq, 'X')                             # Q = Lp[T]/(f) = F (x) Lp
X = Q.gen()


# ---- field element <-> integer bit-vector helpers (GF(2^k) as F2[T]/(f)) ----
def _mk(K):
    g = K.gen(); d = K.degree()
    def from_int(v):
        e = K(0); t = K(1)
        for j in range(d):
            if (v >> j) & 1: e += t
            t *= g
        return e
    def to_int(x):
        v = 0
        for j, c in enumerate(x.polynomial().list()):
            if c: v |= 1 << j
        return v
    return from_int, to_int

# ------------------------------------------------------------------- tau maps
# R' = F[Z]/(Phi)  <->  Q = Lp[T]/(f):   sum_i c_i Z^i  <->  sum_j (sum_i c_{ij} z^i) X^j
# with c_i = sum_j c_{ij} T^j, c_{ij} in GF(2).
Ffrom, Fto = _mk(F)
Lpfrom, Lpto = _mk(Lp)

def tau(coeffs):
    """coeffs: list of n-1 elements of F  ->  element of Q"""
    out = [Lp(0)] * m
    for i, c in enumerate(coeffs):
        bits = Fto(c)
        zi = z ** i
        for j in range(m):
            if (bits >> j) & 1:
                out[j] += zi
    return Q(out)

def tau_inv(y):
    """element of Q -> list of n-1 elements of F"""
    yl = y.lift().list()
    yl += [Lp(0)] * (m - len(yl))
    cols = [0] * (n - 1)
    for j in range(m):
        b = Lpto(yl[j])
        for i in range(n - 1):
            if (b >> i) & 1:
                cols[i] ^^= (1 << j)
    return [Ffrom(c) for c in cols]

def rho(y):
    """partial Frobenius id_F (x) sigma_{Lp}: square every Lp-coefficient."""
    yl = y.lift().list()
    yl += [Lp(0)] * (m - len(yl))
    return Q([c ** 2 for c in yl])

# --------------------------------------------------------------- the instance
# canonical LOCKER instance over R = F[Z]/(Z^n-1): e1, e2 have all F-coordinates
# in a common r-dimensional GF(2)-subspace E of F; h uniform; u = e1 + e2*h.
def make_instance(seed):
    set_random_seed(seed)
    def rowE(v):
        return [GF(2)((v >> j) & 1) for j in range(m)]
    def rspace(vs):
        return matrix(GF(2), len(vs), m, [rowE(v) for v in vs]).row_space()
    a = [Ffrom(randint(1, 2 ** m - 1)) for _ in range(r)]
    while rspace([Fto(x) for x in a]).dimension() < r:
        a = [Ffrom(randint(1, 2 ** m - 1)) for _ in range(r)]
    E = rspace([Fto(x) for x in a])

    def rand_E():
        cs = [randint(0, 1) for _ in range(r)]
        v = 0
        for k in range(r):
            if cs[k]:
                v ^^= Fto(a[k])
        return Ffrom(v)

    e1c = [rand_E() for _ in range(n)]
    e2c = [rand_E() for _ in range(n)]
    hc = [Ffrom(randint(0, 2 ** m - 1)) for _ in range(n)]

    RF = PolynomialRing(F, 'Z')
    Rq = RF.quotient(RF(Zn1))
    e1 = Rq(RF(e1c)); e2 = Rq(RF(e2c)); h = Rq(RF(hc))
    u = e1 + e2 * h
    return E, a, e1c, e2c, hc, u.lift().list()

def project(coeffs):
    """R = F[Z]/(Z^n-1) -> R' = F[Z]/(Phi): reduce mod Phi."""
    RF = PolynomialRing(F, 'Z')
    return (RF(coeffs) % RF(Phi.change_ring(F))).list()

# ------------------------------------------------------- MaxMinors on R' = F(x)Lp
Ts = list(combinations(range(p), r))
cols = [(T, s) for T in Ts for s in range(2 ** r)]
cidx = {c: i for i, c in enumerate(cols)}

def attack(hq, uq, verbose=True):
    """hq, uq: lists of p elements of Q (h', u' and their rho-powers).
    Returns the Lp-kernel dimension and a basis of the Plucker space."""
    H = [hq[0]]; U = [uq[0]]
    for q in range(1, p):
        H.append(rho(H[-1])); U.append(rho(U[-1]))
    if verbose:
        assert H[1] == hq[1] and U[1] == uq[1]

    # subset products of H
    prod = {0: Q(1)}
    for mask in range(1, 1 << p):
        b = mask & (-mask); i = b.bit_length() - 1
        prod[mask] = prod[mask ^^ b] * H[i]

    rows = []
    for Sset in combinations(range(p), r + 1):
        # coefficient (in Q) of each Plucker coordinate
        coefs = {}
        for jj in range(r + 1):
            j = Sset[jj]
            T = tuple(x for x in Sset if x != j)
            for s in range(2 ** r):
                hm = 0
                for i in range(r):
                    if (s >> i) & 1:
                        hm |= 1 << T[i]
                coefs[(T, s)] = coefs.get((T, s), Q(0)) + U[j] * prod[hm]
        # one Q-equation -> m equations over Lp (coefficient of X^0..X^{m-1})
        for deg in range(m):
            row = [Lp(0)] * len(cols)
            for c, v in coefs.items():
                vl = v.lift().list(); vl += [Lp(0)] * (m - len(vl))
                row[cidx[c]] = vl[deg]
            rows.append(row)
    Mx = matrix(Lp, rows)
    K = Mx.right_kernel()
    return Mx, K

def recover_support(K, hq, uq):
    """from the 1-dimensional Plucker kernel: generator -> coefficients -> e2'."""
    H = [hq[0]]; U = [uq[0]]
    for q in range(1, p):
        H.append(rho(H[-1])); U.append(rho(U[-1]))
    vec = K.basis()[0]
    P = {c: vec[cidx[c]] for c in cols}
    chart = tuple((i, 1) for i in range(r))
    cm = P[(tuple(range(r)), (1 << r) - 1)]
    assert cm != 0
    ci = cm ** -1
    G = [[Lp(0)] * (2 * p) for _ in range(r)]
    for i in range(r):
        G[i][2 * i + 1] = Lp(1)
    def plucker(sel):
        T = tuple(sorted(q for q, s in sel))
        s = 0
        for i, q in enumerate(T):
            for qq, ss in sel:
                if qq == q:
                    s |= ss << i
        return P[(T, s)]
    for pair in range(r, p):
        for side in (0, 1):
            for row in range(r):
                sel = [(i, 1) for i in range(r) if i != row] + [(pair, side)]
                G[row][2 * pair + side] = plucker(sel) * ci
    for i in range(r):
        G[i][2 * i] = plucker([(j, 1 if j != i else 0) for j in range(r)]) * ci
    # columns (i,0) with i < r duplicate a chart pair, so use a two-element
    # chart replacement (same route as kemsolve.cpp) to fill rows != i
    for i in range(r):
        pair = side = -1
        for pp in range(r, p):
            for ss in (0, 1):
                if G[i][2 * pp + ss] != 0:
                    pair, side = pp, ss; break
            if pair >= 0: break
        assert pair >= 0, "no nonsingular alternate chart"
        for row in range(r):
            if row == i: continue
            sel = [(j, 1) for j in range(r) if j != i and j != row] + [(i, 0), (pair, side)]
            d = plucker(sel) * ci
            G[row][2 * i] = (d + G[i][2 * i] * G[row][2 * pair + side]) / G[i][2 * pair + side]
    # solve  U_q = sum_k c_k (G[k][2q] + H_q G[k][2q+1])   over Q, r unknowns
    Amat = matrix(Q, p, r, [Q(G[k][2 * q]) + H[q] * Q(G[k][2 * q + 1])
                            for q in range(p) for k in range(r)])
    b = vector(Q, U)
    sol = Amat.solve_right(b)
    e2q = sum(sol[k] * Q(G[k][1]) for k in range(r))
    return e2q

# ------------------------------------------------------------------- the test
log("=" * 78)
log("canonical ring R = F[Z]/(Z^%d-1),  F = GF(2^%d),  r = %d,  p = %d" % (n, m, r, p))
log("Z^%d-1 = (Z+1) * Phi,  Phi irreducible of degree %d" % (n, n - 1))
log("N_u = C(p,r)*2^r = %d,  N_e = C(p-1,r)*m = %d  -> overdetermined"
    % (len(cols), len(list(combinations(range(p), r + 1))) * m))
log("=" * 78)

ok = 0
for seed in range(1, 21):
    E, a, e1c, e2c, hc, uc = make_instance(seed)

    # ---- (C) why the naive attack fails: the Z->1 component is vacuous
    # rho is the identity there, so all rows of the MinRank matrix coincide and
    # every (r+1)x(r+1) minor vanishes.  Check numerically: evaluate every
    # equation coefficient at Z = 1.
    def ev1(coeffs):
        return sum(coeffs)
    H1 = [ev1(hc)] * p
    U1 = [ev1(uc)] * p
    prod1 = {0: F(1)}
    for mask in range(1, 1 << p):
        b = mask & (-mask); i = b.bit_length() - 1
        prod1[mask] = prod1[mask ^^ b] * H1[i]
    vac = True
    for Sset in combinations(range(p), r + 1):
        acc = F(0)
        for jj in range(r + 1):
            j = Sset[jj]; T = tuple(x for x in Sset if x != j)
            for s in range(2 ** r):
                hm = 0
                for i in range(r):
                    if (s >> i) & 1:
                        hm |= 1 << T[i]
                acc += U1[j] * prod1[hm]
        if acc != 0:
            vac = False
    # ---- (B) the fix: run the attack in component 2 (mod Phi)
    e2p = project(e2c)
    hq = [tau(project(hc))]
    uq = [tau(project(uc))]
    hq.append(rho(hq[0])); uq.append(rho(uq[0]))
    Mx, K = attack(hq, uq)
    dim = K.dimension()
    line = "seed %d: naive Z->1 component vacuous = %s | component-2 kernel dim = %d" % (seed, vac, dim)
    if dim == 1:
        e2rec = recover_support(K, hq, uq)
        cf = tau_inv(e2rec)
        Er = matrix(GF(2), len(cf), m,
                    [[GF(2)((Fto(x) >> j) & 1) for j in range(m)] for x in cf]).row_space()
        good = (Er == E)
        ok += int(good)
        line += " | support recovered = %s (dim %d)" % (good, Er.dimension())
    log(line)

log("-" * 78)
log("component-2 attack recovered the support on %d/20 canonical-ring instances" % ok)
