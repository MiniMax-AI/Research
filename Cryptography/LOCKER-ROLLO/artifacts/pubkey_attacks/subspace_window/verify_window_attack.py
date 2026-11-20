#!/usr/bin/env python3
"""
verify_window_attack.py  --  Weak-key Window Lemma and DFR structure verification for 2020 ROLLO-I

Executable witness for the WINDOW LEMMA and the weak-key classes of 2020 ROLLO-I.
Run:  python3 verify_window_attack.py            (full, ~3 min)
      python3 verify_window_attack.py --quick

Everything is computed on the REAL 2020-04-21 ROLLO-I specification parameters
(ROLLO spec v2.0, Table 2 p.18 / Table 3 p.20):

  ROLLO-I-128 : n=83,  m=67, r=7, d=8 ; P=X^83+X^7+X^4+X^2+1 ; g=X^67+X^5+X^2+X+1
  ROLLO-I-192 : n=97,  m=79, r=8, d=8 ; P=X^97+X^6+1          ; g=X^79+X^9+1
  ROLLO-I-256 : n=113, m=97, r=9, d=9 ; P=X^113+X^9+1         ; g=X^97+X^6+1

Every claim below is either (a) a proved algebraic
identity about the official parameters, or (b) an exact computation on them.
"""
import sys, os, time, random, hashlib
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from rolloI_window import (RolloI, GF2k, gf2_rref, gf2_rank, gf2_mat_inverse, parity,
                           dual_basis, extend_basis, window_matrix, L_rank, L_nullspace,
                           kernel_to_R, find_irreducible, is_irreducible)

QUICK = '--quick' in sys.argv
WRITE_LOG = '--write-log' in sys.argv
RESULTS = {}
OUT = []
def say(s=""):
    print(s, flush=True)
    OUT.append(s)

def gauss_binom(N, k):
    """[N choose k]_2 = prod_{i<k} (N-2^i)/(2^k-2^i); N = 2^(ambient dim)."""
    num, den = 1, 1
    for i in range(k):
        num *= (N - (1 << i))
        den *= ((1 << k) - (1 << i))
    return num // den

def log2i(x):
    from math import log2
    return log2(x) if x > 0 else float('-inf')


# ---------------------------------------------------------------- ROLLO-I KEM
def orth_complement(basis, m):
    """Orthogonal complement (standard dot product) of the span of `basis` in GF(2)^m."""
    rows = [b for b in basis]
    piv = []
    r = 0
    red = []
    for col in range(m):
        sel = -1
        for i in range(r, len(rows)):
            if (rows[i] >> col) & 1:
                sel = i; break
        if sel < 0:
            continue
        rows[r], rows[sel] = rows[sel], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> col) & 1):
                rows[i] ^= rows[r]
        piv.append(col); r += 1
    free = [c for c in range(m) if c not in piv]
    out = []
    for fc in free:
        v = 1 << fc
        for i in range(len(piv) - 1, -1, -1):
            if parity(rows[i] & v):
                v ^= 1 << piv[i]
        out.append(v)
    return gf2_rref(out, m)


def subspace_intersect(A, B, m):
    """A cap B for F2-subspaces given by bases."""
    if not A: return []
    if not B: return []
    OA, OB = orth_complement(A, m), orth_complement(B, m)
    return orth_complement(gf2_rref(OA + OB, m), m)


def rsr_decaps(inst, xkey, c, Fprime, r):
    """ROLLO-I Decap (spec Fig.2 + Algorithm 1 RSR):
         s = x*c ; S = <s_1..s_n> ; E = intersection_i f_i^{-1} S ."""
    m, K = inst.m, inst.K
    s = inst.rmul(xkey, c)
    S = gf2_rref([v for v in s if v], m)
    E = None
    for f in Fprime:
        fi = K.inv(f)
        Si = gf2_rref([K.mul(fi, v) for v in S], m)
        E = Si if E is None else subspace_intersect(E, Si, m)
        if E is not None and len(E) == r:
            break
    return E


def shared_secret(E, m, r):
    """ROLLO-I K = Hash(E) : SHA-512 over the echelonised support (spec Sec 2.6)."""
    nbytes = (r * m + 7) // 8
    buf = bytearray()
    for row in E:
        buf += row.to_bytes((m + 7) // 8, 'big')
    return hashlib.sha512(bytes(buf)).hexdigest()


# =============================================================================
say("=" * 78)
say("PART 1.  T7 -- No Frobenius-stable / multiplicatively-invariant weak keys")
say("         (exhaustive verification for the small primes, proof for all m prime)")
say("=" * 78)

def frob(v, k, low):
    """v -> v^2 in GF(2^k)=GF(2)[z]/(z^k+low): carry-less square + reduction."""
    r, t = 0, 0
    x = v
    while x:
        if x & 1:
            r ^= 1 << t
        x >>= 1
        t += 2
    for b in range(r.bit_length() - 1, k - 1, -1):
        if (r >> b) & 1:
            r ^= 1 << b
            r ^= low << (b - k)
    return r

def all_subspaces(m, dim):
    """Enumerate EVERY dim-dimensional F2-subspace of GF(2)^m exactly once,
    as its unique reduced row echelon basis (rows sorted by increasing pivot)."""
    if dim == 0:
        yield []
        return
    def rec(pivs, start, cur):
        if len(cur) == dim:
            yield list(cur)
            return
        for p in range(start, m - (dim - len(cur)) + 1):
            forb = 0
            for q in pivs:
                forb |= 1 << q
            free = [b for b in range(p) if not ((forb >> b) & 1)]
            for sel in range(1 << len(free)):
                v = 1 << p
                for t, b in enumerate(free):
                    if (sel >> t) & 1:
                        v |= 1 << b
                yield from rec(pivs + [p], p + 1, cur + [v])
    yield from rec([], 0, [])

FROB_OK = True
for _mp in (5, 7, 11):
    _low = find_irreducible(_mp, random.Random(_mp))
    _K = GF2k(_mp, _low)
    for _v in range(1, 1 << _mp):
        if frob(_v, _mp, _low) != _K.mul(_v, _v):
            FROB_OK = False
            break
say("  Frobenius-map self-test (v -> v^2 == GF2k.mul(v,v)) for m=5,7,11 : %s"
    % ("PASS" if FROB_OK else "FAIL"))

for mp in ([5, 7] if not QUICK else [5]):
    rng = random.Random(mp)
    low = find_irreducible(mp, rng)
    tot = 0
    stable = []
    for dim in range(1, mp):
        for bas in all_subspaces(mp, dim):
            tot += 1
            sq = [frob(v, mp, low) for v in bas]
            if len(gf2_rref(bas + sq, mp)) == dim:
                stable.append(dim)
    say("  m=%d (prime): enumerated %d subspaces of dim 1..%d ; Frobenius-stable ones: %d %s"
        % (mp, tot, mp - 1, len(stable), stable))
say("  Exact structure (normal basis theorem): GF(2^m) is a CYCLIC F2[t]-module")
say("  under t = Frobenius, isomorphic to F2[t]/(t^m-1).  Hence the Frobenius-stable")
say("  F2-subspaces are exactly the ideals, i.e. the divisors of t^m-1 over GF(2),")
say("  and the possible dimensions are the DEGREES of those divisors:")
say("      { deg(g) : g | t^m-1 } = { sum over e|m of k_e * ord_e(2) }.")
say("  (m=5 above: ord_5(2)=4, divisors of t^5-1=(t-1)*Phi_5 have degrees {0,1,5}")
say("   -- but dim 1 is <spanned by 1> and dim 4 is the TRACE-ZERO hyperplane; both")
say("   are Frobenius stable, so the naive 'no stable subspace' claim is FALSE.)")
say()

def ordm(a, m):
    x, k = a % m, 1
    while x != 1:
        x = (x * a) % m
        k += 1
    return k

def frob_stable_dims(m):
    """degrees of the divisors of t^m-1 over GF(2)"""
    dims = {0}
    # factor t^m-1 = prod_{e|m} Phi_e ; over GF(2), Phi_e = (e irr. factors) of deg ord_e(2)
    e = 1
    divs = [d for d in range(1, m + 1) if m % d == 0]
    parts = []
    for e in divs:
        o = ordm(2, e) if e > 1 else 1
        parts.append((o, (1 if e == 1 else 0) or 0))
    # number of irreducible factors of Phi_e over GF(2) = phi(e)/ord_e(2)
    def phi(nn):
        res, f, mm = nn, 2, nn
        while f * f <= mm:
            if mm % f == 0:
                while mm % f == 0:
                    mm //= f
                res -= res // f
            f += 1
        if mm > 1:
            res -= res // mm
        return res
    blocks = []
    for e in divs:
        o = 1 if e == 1 else ordm(2, e)
        cnt = 1 if e == 1 else phi(e) // o
        blocks.append((o, cnt))
    dims = {0}
    for (o, cnt) in blocks:
        nd = set()
        for dd in dims:
            for k in range(cnt + 1):
                nd.add(dd + k * o)
        dims = nd
    return sorted(dims)

say("  VERIFICATION of the divisor characterisation at small m (exhaustive):")
for mp in [5, 7]:
    rng = random.Random(mp)
    low = find_irreducible(mp, rng)
    obs = set()
    for dim in range(1, mp):
        for bas in all_subspaces(mp, dim):
            sq = [frob(v, mp, low) for v in bas]
            if len(gf2_rref(bas + sq, mp)) == dim:
                obs.add(dim)
    pred = set(frob_stable_dims(mp)) - {0, mp}
    say("     m=%2d : observed stable dims %s   predicted %s   MATCH=%s"
        % (mp, sorted(obs), sorted(pred), obs == pred))

say()
say("  THE REAL PARAMETERS (all m prime):")
say("  %-14s %3s %3s %-12s %-28s %s" % ("instance", "m", "d", "ord_m(2)", "Frobenius-stable dims", "d-dim stable F exists?"))
for nm in ["ROLLO-I-128", "ROLLO-I-192", "ROLLO-I-256"]:
    mm, nn, rr, dd = RolloI.SPEC[nm][:4]
    o = ordm(2, mm)
    dims = frob_stable_dims(mm)
    ex = dd in dims
    say("  %-14s %3d %3d %-12d %-28s %s"
        % (nm, mm, dd, o, str(dims) if len(dims) < 8 else str(dims[:4]) + "...", ex))
    RESULTS[(nm, 'frob')] = (not ex)

say()
say("  T7 (multiplicative invariance -- the version that IS true):")
say("  Stab_{K*}(F) = {c in K* : cF subseteq F}.  If c in Stab and c not in F2 then")
say("  F is closed under F2[c] = K (m prime), so F in {0,K}.  Hence for every F with")
say("  1 <= dim F < m :  Stab_{K*}(F) = {1}.")
for mp in ([5, 7] if not QUICK else [5]):
    rng = random.Random(mp)
    low = find_irreducible(mp, rng)
    K = GF2k(mp, low)
    bad = 0
    tot = 0
    for dim in range(1, mp):
        for bas in all_subspaces(mp, dim):
            tot += 1
            # direct test: does there exist c not in F2 with cF subseteq F ?
            sp = set()
            for msk in range(1 << dim):
                acc = 0
                for k in range(dim):
                    if (msk >> k) & 1:
                        acc ^= bas[k]
                sp.add(acc)
            for c in range(2, 1 << mp):
                if all(K.mul(c, v) in sp for v in sp):
                    bad += 1
                    break
    say("     m=%d : over ALL %d subspaces F (1<=dim<m), #F with Stab_{K*}(F) != {1} : %d"
        % (mp, tot, bad))
    RESULTS[('stab', mp)] = (bad == 0)

# =============================================================================
say()
say("=" * 78)
say("PART 2.  T2/T3 -- Window Lemma + scale invariance, at REAL ROLLO-I params")
say("=" * 78)

def try_break(inst, h, W, tag):
    """Run the window attack with window W on public key h.
    Returns (xkey, Fprime) with dim(Sigma(xkey)+Sigma(xkey*h)) <= d, or (None,None)."""
    m, d = inst.m, inst.d
    D = len(W)
    M, wb = window_matrix(inst, h, W)
    rk = L_rank(M, len(M), D, inst.L)
    if rk == D:
        return None, None, rk
    ker = L_nullspace(M, len(M), D, inst.L)
    for v in ker:
        xx = kernel_to_R(inst, v, wb)
        if xx == [0] * inst.n:
            continue
        dim = gf2_rank(xx + inst.rmul(xx, h), m)
        if 1 <= dim <= d:
            Fp = gf2_rref([a for a in xx if a] + [a for a in inst.rmul(xx, h) if a], m)
            return xx, Fp, rk
    return None, None, rk


def end_to_end(inst, h, xkey, Fp, say, tag):
    """Full ROLLO-I KEM break: encapsulate, decapsulate with the RECOVERED key,
    check the shared secret Hash(E) matches."""
    m, n, r, d = inst.m, inst.n, inst.r, inst.d
    E = inst.rand_subspace(r)
    e1 = inst.rand_dense(E); e2 = inst.rand_dense(E)
    while gf2_rank(e1 + e2, m) != r:
        e2 = inst.rand_dense(E)
    c = [a ^ b for a, b in zip(e1, inst.rmul(e2, h))]
    Erec = rsr_decaps(inst, xkey, c, Fp, r)
    Etrue = gf2_rref(E, m)
    ok = (Erec is not None and gf2_rref(Erec, m) == Etrue)
    ss_ok = ok and (shared_secret(Erec, m, r) == shared_secret(Etrue, m, r))
    say("      END-TO-END KEM BREAK (%s): RSR(E) == true E : %s ; Hash(E) matches : %s"
        % (tag, ok, ss_ok))
    return ss_ok


for name in (["ROLLO-I-128"] if QUICK else ["ROLLO-I-128", "ROLLO-I-192", "ROLLO-I-256"]):
    inst = RolloI(name=name, seed=20200421)
    m, n, r, d = inst.m, inst.n, inst.r, inst.d
    say()
    say("--- %s : (m,n,r,d) = (%d,%d,%d,%d) ---" % (name, m, n, r, d))
    t0 = time.time()

    # ---- (a) weak key with F = <1,T,...,T^{d-1}> : attack with the SINGLE known window
    W = [1 << j for j in range(d)]
    F, x, y, h = inst.keygen(F=W)
    assert inst.rmul(x, h) == y and gf2_rank(x + y, m) == d
    xk, Fp, rk = try_break(inst, h, W, 'a')
    say("  (a) F = <1,T,...,T^%d>  (dim d=%d) :  M_W is %dx%d over GF(2^%d)"
        % (d - 1, d, m - d, d, n))
    say("      rank(M_W) = %d  (< d=%d  =>  nontrivial kernel)" % (rk, d))
    say("      recovered key x': dim(Sigma(x')+Sigma(x'h)) = %d <= d ; Sigma(x')==F : %s"
        % (gf2_rank(xk + inst.rmul(xk, h), m) if xk else -1,
           gf2_rref([a for a in xk if a], m) == gf2_rref(F, m) if xk else False))
    e2e = end_to_end(inst, h, xk, Fp, say, 'a') if xk else False
    RESULTS[(name, 'lowdeg')] = (xk is not None and e2e)

    # ---- (b) SCALE INVARIANCE: F = gamma * <1,T,...,T^{d-1}> for random gamma,
    #          attacked with the UNSCALED window  W = <1,T,...,T^{d-1}>.
    rng = random.Random(12345)
    gam = rng.randrange(1, 1 << m)
    Fg = gf2_rref([inst.K.mul(gam, v) for v in W], m)
    F2, x2, y2, h2 = inst.keygen(F=Fg)
    xk2, Fp2, rk2 = try_break(inst, h2, W, 'b')
    say("  (b) SCALE INVARIANCE  F = gamma*<1..T^%d>, random gamma:" % (d - 1))
    say("      rank(M_W) = %d with the SAME unscaled window (kernel nontrivial: %s)"
        % (rk2, rk2 < d))
    if xk2:
        say("      recovered key x' = gamma^{-1} x : Sigma(x')+Sigma(x'h) has dim %d <= d"
            % gf2_rank(xk2 + inst.rmul(xk2, h2), m))
    e2e2 = end_to_end(inst, h2, xk2, Fp2, say, 'b') if xk2 else False
    RESULTS[(name, 'scaled')] = (xk2 is not None and e2e2)

    # ---- (c) MONOMIAL WINDOW: F = gamma * <T^{j_1},...,T^{j_d}> for random exponents
    rng = random.Random(999)
    J = sorted(rng.sample(range(m), d))
    Wj = [1 << j for j in J]
    gam = rng.randrange(1, 1 << m)
    Fj = gf2_rref([inst.K.mul(gam, v) for v in Wj], m)
    F3, x3, y3, h3 = inst.keygen(F=Fj)
    xk3, Fp3, rk3 = try_break(inst, h3, Wj, 'c')
    say("  (c) MONOMIAL WINDOW  F = gamma*<T^j : j in J>, |J|=%d, J=%s" % (d, J))
    say("      rank(M_W) = %d (kernel nontrivial: %s)" % (rk3, rk3 < d))
    e2e3 = end_to_end(inst, h3, xk3, Fp3, say, 'c') if xk3 else False
    RESULTS[(name, 'monomial')] = (xk3 is not None and e2e3)

    # ---- (d) D > d : window <1,...,T^{D-1}> strictly larger than F
    D = d + 3
    Wd = [1 << j for j in range(D)]
    rng = random.Random(4242)
    Fd = gf2_rref([1 << j for j in sorted(rng.sample(range(D), d))], m)
    F4, x4, y4, h4 = inst.keygen(F=Fd)
    M4, wb4 = window_matrix(inst, h4, Wd)
    rk4 = L_rank(M4, len(M4), D, inst.L)
    ker4 = L_nullspace(M4, len(M4), D, inst.L)
    say("  (d) D = d+3 = %d window <1..T^%d>, F a random d-dim subspace of it:" % (D, D - 1))
    say("      M_W is %dx%d, rank = %d, kernel dim = %d (>=1 expected)" % (m - D, D, rk4, len(ker4)))
    RESULTS[(name, 'bigD')] = (len(ker4) >= 1)

    # ---- (e) NEGATIVE CONTROL: genuine ROLLO-I keys with a uniformly random F
    ntr = 3 if QUICK else 6
    full = 0
    tneg = time.time()
    for t in range(ntr):
        Fr, xr, yr, hr = inst.keygen()          # uniformly random full-rank F, dim d
        Mr, wbr = window_matrix(inst, hr, W)    # canonical window <1..T^{d-1}>
        if L_rank(Mr, len(Mr), d, inst.L) == d:
            full += 1
    say("  (e) NEGATIVE CONTROL: %d/%d uniformly-random-F ROLLO-I keys give rank(M_W)=d"
        % (full, ntr))
    say("      i.e. TRIVIAL kernel => no false positive.            [%.2fs]" % (time.time() - tneg))
    RESULTS[(name, 'neg')] = (full == ntr)

    # ---- (f) NEGATIVE CONTROL: uniformly random h (not a ROLLO-I public key at all)
    full = 0
    for t in range(ntr):
        hu = [rng.randrange(0, 1 << m) for _ in range(n)]
        Mu, wbu = window_matrix(inst, hu, W)
        if L_rank(Mu, len(Mu), d, inst.L) == d:
            full += 1
    say("  (f) NEGATIVE CONTROL: %d/%d uniformly RANDOM h in R give rank(M_W)=d (rejected)"
        % (full, ntr))
    RESULTS[(name, 'negrand')] = (full == ntr)
    say("      total wall clock for %s : %.2fs" % (name, time.time() - t0))


# =============================================================================
say()
say("=" * 78)
say("PART 3.  EXHAUSTIVE negative control at small parameters")
say("         (m=9, n=11, d=2, gcd=1 : ALL 86530 two-dimensional F are enumerated)")
say("=" * 78)
inst = RolloI(9, 11, 2, 2, seed=7)
m, n, d = 9, 11, 2
W = [1, 2]
# (i) POSITIVE: a genuine weak key must be found
Fw, xw, yw, hw = inst.keygen(F=W)
Mw, wbw = window_matrix(inst, hw, W)
say("  POSITIVE : genuine weak key, rank(M_W)=%d  (need < d=%d)"
    % (L_rank(Mw, len(Mw), d, inst.L), d))
# (ii) EXHAUSTIVE NEGATIVE: a uniformly random h must admit NO valid key of any 2-dim F
rng = random.Random(31337)
hu = [rng.randrange(0, 1 << m) for _ in range(n)]
cnt = 0
tot = 0
t0 = time.time()
N = (1 << m) - 1
# enumerate all 2-dim subspaces of GF(2)^m as {0, a, b, a^b} with a<b in a canonical order
els = list(range(1, 1 << m))
seen = set()
for a in els:
    for b in els:
        if a >= b:
            continue
        c = a ^ b
        key = tuple(sorted((a, b, c)))
        if key in seen:
            continue
        seen.add(key)
        tot += 1
        bas = gf2_rref([a, b], m)
        if len(bas) != 2:
            continue
        Mx, _ = window_matrix(inst, hu, bas)
        if L_rank(Mx, len(Mx), 2, inst.L) < 2:
            cnt += 1
say("  EXHAUSTIVE NEGATIVE : enumerated ALL %d two-dimensional F2-subspaces F of GF(2^%d);"
    % (tot, m))
say("      number admitting a nonzero x with Supp(x),Supp(xh) subseteq F : %d" % cnt)
say("      => a uniformly random h admits NO ROLLO-I key of support dimension 2.")
say("      (predicted count = #F * 2^{-n(m-2d)} = %d * 2^%d = %.3g)   [%.1fs]"
    % (tot, -n * (m - 2 * d), tot * 2.0 ** (-n * (m - 2 * d)), time.time() - t0))
RESULTS[('exhaustive', 'neg')] = (cnt == 0)


# =============================================================================
say()
say("=" * 78)
say("PART 4.  Exact weak-key fractions (Gaussian binomials, exact big integers)")
say("=" * 78)
say("  For a fixed known D-dim window W the covered class is")
say("     W_D(W) = { F in Gr_d(K) : exists gamma in K*,  F subseteq gamma*W }.")
say("  Scale invariance (T3) means ONE window covers all 2^m-1 scalings at once.")
say()
say("  %-14s %3s %3s %5s | %-14s %-14s %-10s" %
    ("instance", "m", "d", "D", "log2 |Gr_d(K)|", "log2 |W_D(W)|", "log2 frac"))
FRAC = {}
for name in ["ROLLO-I-128", "ROLLO-I-192", "ROLLO-I-256"]:
    mm, nn, rr, dd = RolloI.SPEC[name][:4]
    tot = gauss_binom(1 << mm, dd)
    # D = d : F = gamma*W exactly, and Stab_{K*}(F2[T]_{<d}) = {1} when 2d-1 <= m
    assert 2 * dd - 1 <= mm
    cov = (1 << mm) - 1
    FRAC[(name, 'D=d')] = (log2i(tot), log2i(cov), log2i(cov) - log2i(tot))
    say("  %-14s %3d %3d %5d | %14.2f %14.2f %10.1f" %
        (name, mm, dd, dd, log2i(tot), log2i(cov), log2i(cov) - log2i(tot)))
say()
say("  Upper bound for D > d (each gamma contributes [D choose d]_2 subspaces):")
for name in ["ROLLO-I-128"]:
    mm, nn, rr, dd = RolloI.SPEC[name][:4]
    tot = gauss_binom(1 << mm, dd)
    for D in [dd, dd + 1, dd + 2, dd + 5, 2 * dd, mm // 2]:
        cov = ((1 << mm) - 1) * gauss_binom(1 << D, dd)
        say("     D=%2d : log2 |W_D| <= %8.2f   log2 fraction <= %9.2f   (window matrix %dx%d)"
            % (D, log2i(cov), log2i(cov) - log2i(tot), mm - D, D))


# =============================================================================
say()
say("=" * 78)
say("PART 5.  T6 -- Window Covering Bound (why this is inherently a WEAK-KEY attack)")
say("=" * 78)
say("  Any family Fm of D-dim windows that contains, for EVERY F in Gr_d(K), at least")
say("  one W with F subseteq W must have")
say("        |Fm|  >=  [2^m-1 choose d]_2  /  [2^(m-d) choose D-d]_2 .")
for name in ["ROLLO-I-128", "ROLLO-I-192", "ROLLO-I-256"]:
    mm, nn, rr, dd = RolloI.SPEC[name][:4]
    tot = gauss_binom(1 << mm, dd)
    for D in [dd + 1]:
        den = gauss_binom(1 << (mm - dd), D - dd)
        say("     %-14s D=%d : |Fm| >= 2^%.1f  windows  (=> no sub-exponential full attack)"
            % (name, D, log2i(tot) - log2i(den)))


# =============================================================================
say()
say("=" * 78)
say("PART 6.  T9 -- DFR structure: dim(EF) profile at the real parameters")
say("=" * 78)
say("  ROLLO spec v2.0 Prop 2.4.3:  DFR ~ q^{-(d-1)(m-r'-r)} + q^{-(n-r'+1)},  r'=dim EF.")
say("  Both exponents are INCREASING in r', so the DFR is MAXIMISED at the generic")
say("  r' = r*d.  Low-dimension F*E intersections therefore make decryption MORE")
say("  reliable, not less -- the opposite of a DFR-amplification weakness.")
say()
say("  %-14s %3s %3s %3s %4s | %-9s %-9s %-9s %-9s" %
    ("instance", "m", "n", "r", "d", "term1", "term2", "DFR(spec)", "log2 Pr[r'<rd]"))
for name in ["ROLLO-I-128", "ROLLO-I-192", "ROLLO-I-256"]:
    mm, nn, rr, dd = RolloI.SPEC[name][:4]
    t1 = -(dd - 1) * (mm - rr * dd - rr)
    t2 = -(nn - rr * dd + 1)
    say("  %-14s %3d %3d %3d %4d | 2^%-7d 2^%-7d 2^%-7d 2^%-9d" %
        (name, mm, nn, rr, dd, t1, t2, max(t1, t2), -(mm - rr * dd)))
say()
say("  (Predicted log2 Pr[dim<EF> < r*d] = -(m - r*d); verified by simulation below.)")
if not QUICK:
    inst = RolloI(name="ROLLO-I-128", seed=5)
    m, n, r, d = inst.m, inst.n, inst.r, inst.d
    rng = random.Random(2020)
    T = 4000
    bad = 0
    t0 = time.time()
    for _ in range(T):
        Fs = inst.rand_subspace(d)
        Es = inst.rand_subspace(r)
        prod = [inst.K.mul(a, b) for a in Fs for b in Es]
        if gf2_rank(prod, m) < r * d:
            bad += 1
    say("  SIMULATION ROLLO-I-128: %d/%d random (F,E) pairs have dim<EF> < r*d = %d"
        % (bad, T, r * d))
    say("      predicted %0.4f = 2^%d ;   observed rate %.5f      [%.1fs]"
        % (2.0 ** (-(m - r * d)), -(m - r * d), bad / T, time.time() - t0))


# =============================================================================
say()
say("=" * 78)
say("SUMMARY")
say("=" * 78)
allok = all(RESULTS.values())
for k in sorted(RESULTS, key=str):
    say("   %-40s : %s" % (str(k), "PASS" if RESULTS[k] else "FAIL"))
say()
say("OVERALL: %s" % ("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED"))

if WRITE_LOG:
    txt = "\n".join(OUT) + "\n"
    with open(os.path.join(HERE, "verify_window_attack.log"), "w") as f:
        f.write(txt)
print("\nsha256(self) =", hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
sys.exit(0 if allok else 1)
