#!/usr/bin/env sage
# canonical_ring_frobenius.sage
#
# Second half of the canonical-ring portability test (see
# canonical_ring_component.sage for the structural argument).
#
# The first script used an OVER-determined parameter set, where the Lp-kernel of
# the MaxMinors linearisation is already 1-dimensional.  The LOCKER-I
# parameters (m = 71, r = 5, p = 9) are UNDER-determined by exactly
#       dim = C(p-1,r) * (2^r*p - m*(p-r)) / (p-r)
# and need the Frobenius relations  P[T+1][s] = P[T][s]^2  to cut the kernel to a
# single line.  This script repeats the canonical-ring test in that regime, so
# the WHOLE pipeline (linearisation -> Frobenius completion -> generator ->
# coefficients -> support) is exercised on R = F[Z]/(Z^n - 1).
#
# Parameters: n = 11, r = 3, p = 7, m = 13 -> dim = 20*(56-52)/4 = 20 > 1.
# (m must be coprime to n-1 so that F (x) GF(2^{n-1}) is a field, and the number
#  of shift relations C(p-1,r)*2^r = 160 must comfortably exceed the deficiency,
#  exactly as 1792 >> 56 does at the LOCKER-I parameters m=71, r=5, p=9.)
#
# Run:  sage canonical_ring_frobenius.sage

from itertools import combinations
import sys

n, r, p, m = 11, 3, 7, 13
assert (2 ** 10) % n == 1 and (2 ** 5) % n == n - 1

def log(*a):
    print(*a); sys.stdout.flush()

ZF = PolynomialRing(GF(2), 'Z'); Zv = ZF.gen()
Zn1 = Zv ** n + 1
Phi = sum(Zv ** i for i in range(n))
assert Zn1 == (Zv + 1) * Phi and Phi.is_irreducible()

F = GF(2 ** m, 'T')
Lp = GF(2 ** (n - 1), 'z', modulus=Phi); z = Lp.gen()
S = PolynomialRing(Lp, 'X')
fq = S([Lp(c) for c in F.modulus().list()])
assert fq.is_irreducible()
Q = S.quotient(fq, 'X')

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

Ffrom, Fto = _mk(F)
Lpfrom, Lpto = _mk(Lp)
dL = Lp.degree()

def tau(coeffs):
    out = [Lp(0)] * m
    for i, c in enumerate(coeffs):
        bits = Fto(c); zi = z ** i
        for j in range(m):
            if (bits >> j) & 1: out[j] += zi
    return Q(out)

def tau_inv(y):
    yl = y.lift().list(); yl += [Lp(0)] * (m - len(yl))
    cols = [0] * (n - 1)
    for j in range(m):
        b = Lpto(yl[j])
        for i in range(n - 1):
            if (b >> i) & 1: cols[i] ^^= (1 << j)
    return [Ffrom(c) for c in cols]

def rho(y):
    yl = y.lift().list(); yl += [Lp(0)] * (m - len(yl))
    return Q([c ** 2 for c in yl])

def project(coeffs):
    RF = PolynomialRing(F, 'Z')
    return (RF(coeffs) % RF(Phi.change_ring(F))).list()

def make_instance(seed):
    set_random_seed(seed)
    rspace = lambda vs: matrix(GF(2), len(vs), m,
        [[GF(2)((v >> j) & 1) for j in range(m)] for v in vs]).row_space()
    a = [Ffrom(randint(1, 2 ** m - 1)) for _ in range(r)]
    while rspace([Fto(x) for x in a]).dimension() < r:
        a = [Ffrom(randint(1, 2 ** m - 1)) for _ in range(r)]
    E = rspace([Fto(x) for x in a])
    def rand_E():
        v = 0
        for k in range(r):
            if randint(0, 1): v ^^= Fto(a[k])
        return Ffrom(v)
    e1c = [rand_E() for _ in range(n)]
    e2c = [rand_E() for _ in range(n)]
    hc = [Ffrom(randint(0, 2 ** m - 1)) for _ in range(n)]
    RF = PolynomialRing(F, 'Z'); Rq = RF.quotient(RF(Zn1))
    u = Rq(RF(e1c)) + Rq(RF(e2c)) * Rq(RF(hc))
    return E, hc, u.lift().list()

Ts = list(combinations(range(p), r))
cols = [(T, s) for T in Ts for s in range(2 ** r)]
cidx = {c: i for i, c in enumerate(cols)}

def attack(hq):
    H = [hq]; U = [uq0]
    for q in range(1, p):
        H.append(rho(H[-1])); U.append(rho(U[-1]))
    prod = {0: Q(1)}
    for mask in range(1, 1 << p):
        b = mask & (-mask); i = b.bit_length() - 1
        prod[mask] = prod[mask ^^ b] * H[i]
    rows = []
    for Sset in combinations(range(p), r + 1):
        coefs = {}
        for jj in range(r + 1):
            j = Sset[jj]; T = tuple(x for x in Sset if x != j)
            for s in range(2 ** r):
                hm = 0
                for i in range(r):
                    if (s >> i) & 1: hm |= 1 << T[i]
                coefs[(T, s)] = coefs.get((T, s), Q(0)) + U[j] * prod[hm]
        for deg in range(m):
            row = [Lp(0)] * len(cols)
            for c, v in coefs.items():
                vl = v.lift().list(); vl += [Lp(0)] * (m - len(vl))
                row[cidx[c]] = vl[deg]
            rows.append(row)
    Mx = matrix(Lp, rows)
    return H, U, Mx.right_kernel()

def frobenius(K):
    """impose P[T+1][s] = P[T][s]^2 on the Lp-kernel; returns candidate vectors."""
    B = K.basis(); dim = len(B)
    ebasis = [z ** i for i in range(dL)]
    rows = []
    nrel = 0
    for T in combinations(range(p - 1), r):
        T1 = tuple(x + 1 for x in T)
        if T1[-1] > p - 1: continue
        for s in range(2 ** r):
            nrel += 1
            a = [B[f][cidx[(T1, s)]] for f in range(dim)]
            b = [B[f][cidx[(T, s)]] ** 2 for f in range(dim)]
            # one Lp-equation -> dL F2 rows; the entry at (output bit i,
            # unknown lambda_{f,j}) is bit i of  z^j*a_f + z^(2j)*c_f
            tab = [[Lpto(ebasis[j] * a[f] + (ebasis[j] ** 2) * b[f])
                    for j in range(dL)] for f in range(dim)]
            for i in range(dL):
                row = [0] * (dim * dL)
                for f in range(dim):
                    for j in range(dL):
                        if (tab[f][j] >> i) & 1: row[f * dL + j] = 1
                rows.append(row)
    M2 = matrix(GF(2), rows)
    K2 = M2.right_kernel()
    out = []
    for lam_bits in K2.basis():
        vec = [Lp(0)] * len(cols)
        for f in range(dim):
            lam = Lpfrom(sum(int(lam_bits[f * dL + i]) << i for i in range(dL)))
            if lam == 0: continue
            for j in range(len(cols)):
                vec[j] += lam * B[f][j]
        out.append(vec)
    return nrel, K2.dimension(), out

def recover(H, U, vec):
    P = {c: vec[cidx[c]] for c in cols}
    cm = P[(tuple(range(r)), (1 << r) - 1)]
    if cm == 0: return None
    ci = cm ** -1
    G = [[Lp(0)] * (2 * p) for _ in range(r)]
    for i in range(r): G[i][2 * i + 1] = Lp(1)
    def plucker(sel):
        T = tuple(sorted(q for q, s in sel)); s = 0
        for i, q in enumerate(T):
            for qq, ss in sel:
                if qq == q: s |= ss << i
        return P[(T, s)]
    for pair in range(r, p):
        for side in (0, 1):
            for row in range(r):
                G[row][2 * pair + side] = plucker(
                    [(i, 1) for i in range(r) if i != row] + [(pair, side)]) * ci
    for i in range(r):
        G[i][2 * i] = plucker([(j, 1 if j != i else 0) for j in range(r)]) * ci
    for i in range(r):
        pair = side = -1
        for pp in range(r, p):
            for ss in (0, 1):
                if G[i][2 * pp + ss] != 0: pair, side = pp, ss; break
            if pair >= 0: break
        if pair < 0: return None
        for row in range(r):
            if row == i: continue
            sel = [(j, 1) for j in range(r) if j != i and j != row] + [(i, 0), (pair, side)]
            d = plucker(sel) * ci
            G[row][2 * i] = (d + G[i][2 * i] * G[row][2 * pair + side]) / G[i][2 * pair + side]
    Amat = matrix(Q, p, r, [Q(G[k][2 * q]) + H[q] * Q(G[k][2 * q + 1])
                            for q in range(p) for k in range(r)])
    try:
        sol = Amat.solve_right(vector(Q, U))
    except ValueError:
        return None
    return tau_inv(sum(sol[k] * Q(G[k][1]) for k in range(r)))

log("=" * 78)
log("canonical ring R = F[Z]/(Z^%d-1), F = GF(2^%d), r = %d, p = %d" % (n, m, r, p))
log("N_u = %d, N_e = %d  -> UNDER-determined, Frobenius stage required"
    % (len(cols), len(list(combinations(range(p), r + 1))) * m))
log("predicted Lp-kernel dim = C(p-1,r)*(2^r*p - m*(p-r))/(p-r) = %d"
    % (binomial(p - 1, r) * (2 ** r * p - m * (p - r)) // (p - r)))
log("=" * 78)

ok = 0
for seed in range(1, 21):
    E, hc, uc = make_instance(seed)
    hq0 = tau(project(hc)); uq0 = tau(project(uc))
    globals()['uq0'] = uq0
    H, U, K = attack(hq0)
    nrel, kdim, cands = frobenius(K)
    got = None
    for vec in cands:
        cf = recover(H, U, vec)
        if cf is None: continue
        Er = matrix(GF(2), len(cf), m,
                    [[GF(2)((Fto(x) >> j) & 1) for j in range(m)] for x in cf]).row_space()
        if Er == E: got = Er
    ok += int(got is not None)
    log("seed %d: Lp-kernel dim = %d | Frobenius relations = %d -> binary kernel dim = %d"
        " | support recovered = %s" % (seed, K.dimension(), nrel, kdim, got is not None))

log("-" * 78)
log("FULL pipeline (linearisation + Frobenius completion) recovered the support on"
    " %d/20 canonical-ring instances" % ok)
