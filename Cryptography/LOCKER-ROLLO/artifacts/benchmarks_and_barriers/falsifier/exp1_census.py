#!/usr/bin/env python3
"""
exp1_census.py -- EXACT exhaustive census of the TSM solution set over TRUE GF(2^(mn)).

Target degree-1 solvability identity:
   "each conjugate q imposes m-2r > 0 constraints, so the solution is already
    algebraically unique over L once P >= ceil(r(m-r)/(m-2r))"

We enumerate EVERY r-dimensional GF(2)-subspace E' of K = GF(2^m) and decide EXACTLY
whether the TSM system with P conjugates admits a solution for that E'.  We do this for
  (A) the TRUE system      : H_q = rho^q(H_0), U_q = rho^q(U_0)                [real ROLLO]
  (B) the PROBE/RANDOM model: H_q independent uniform, U_q built from independent
                              uniform alpha^(q), beta^(q)   [degree-2 bi-monomial baseline]
and print the exact count as a function of P.
"""
import sys, os, json, random, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollo_gf2_engine import RingQ, rref_subspaces, count_subspaces, F2k

# ---------------------------------------------------------------- instance generation
def rand_subspace(K, dim, rng):
    while True:
        rows, piv = [], []
        for _ in range(dim):
            while True:
                v = rng.randrange(1, K.mask + 1)
                t = v & -v
                if t not in piv: break
            piv.append(v & -v)
            rows.append(v)
        # gaussian elim
        basis, rk = [], 0
        rs = rows[:]
        for c in range(K.k):
            p = -1
            for i in range(rk, len(rs)):
                if (rs[i] >> c) & 1: p = i; break
            if p < 0: continue
            rs[rk], rs[p] = rs[p], rs[rk]
            for i in range(len(rs)):
                if i != rk and ((rs[i] >> c) & 1): rs[i] ^= rs[rk]
            basis.append(rs[rk]); rk += 1
        if rk == dim: return basis

def sample_supported(basis, n, rng, K):
    """random R-element whose n coefficients all lie in span(basis), of full rank."""
    while True:
        v = []
        for _ in range(n):
            msk = rng.randrange(1 << len(basis))
            x = 0
            for k in range(len(basis)):
                if (msk >> k) & 1: x ^= basis[k]
            v.append(x)
        # full GF(2) rank of the coefficients
        rs, rk = v[:], 0
        for c in range(K.k):
            p = -1
            for i in range(rk, len(rs)):
                if (rs[i] >> c) & 1: p = i; break
            if p < 0: continue
            rs[rk], rs[p] = rs[p], rs[rk]
            for i in range(len(rs)):
                if i != rk and ((rs[i] >> c) & 1): rs[i] ^= rs[rk]
            rk += 1
        if rk == len(basis): return v

def canon(basis, k):
    """canonical RREF basis (sorted list of ints) of the GF(2)-span of `basis`."""
    rs = list(basis); out = []; rk = 0
    for c in range(k):
        p = -1
        for i in range(rk, len(rs)):
            if (rs[i] >> c) & 1: p = i; break
        if p < 0: continue
        rs[rk], rs[p] = rs[p], rs[rk]
        for i in range(len(rs)):
            if i != rk and ((rs[i] >> c) & 1): rs[i] ^= rs[rk]
        rk += 1
    return tuple(rs[:rk])   # RREF order = increasing pivot column (NOT sorted by value!)

def make_instance(m, n, r, d, seed):
    rng = random.Random(seed)
    R = RingQ(m, n)
    K = R.K
    F = rand_subspace(K, d, rng)
    while True:
        x = sample_supported(F, n, rng, K)
        if any(x): break
    y = sample_supported(F, n, rng, K)
    xi = R.rinv(x)
    h = R.rmul(xi, y)
    E = rand_subspace(K, r, rng)
    e1 = sample_supported(E, n, rng, K)
    e2 = sample_supported(E, n, rng, K)
    u = [a ^ b for a, b in zip(e1, R.rmul(e2, h))]
    return R, dict(F=F, x=x, y=y, h=h, E=E, e1=e1, e2=e2, u=u)

# ---------------------------------------------------------------- L-linear algebra helpers
def span_contains(vecs, target, L, ncols):
    """True iff target is in the L-span of vecs (list of rows, length ncols)."""
    rows = [v[:] for v in vecs]
    rk = 0
    pivots = []
    for c in range(ncols):
        p = -1
        for i in range(rk, len(rows)):
            if rows[i][c]: p = i; break
        if p < 0: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        iv = L.finv(rows[rk][c])
        if iv != 1: rows[rk] = [L.fmul(x, iv) for x in rows[rk]]
        pr = rows[rk]
        for i in range(rk + 1, len(rows)):
            if rows[i][c]:
                f = rows[i][c]
                rows[i] = [x ^ L.fmul(f, y) for x, y in zip(rows[i], pr)]
        pivots.append((c, pr)); rk += 1
        if rk == len(rows): break
    t = target[:]
    for c, pr in pivots:
        if t[c]:
            f = t[c]
            t = [x ^ L.fmul(f, y) for x, y in zip(t, pr)]
    return not any(t), rk

def rank_L(vecs, L, ncols):
    rows = [v[:] for v in vecs]
    rk = 0
    for c in range(ncols):
        p = -1
        for i in range(rk, len(rows)):
            if rows[i][c]: p = i; break
        if p < 0: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        iv = L.finv(rows[rk][c])
        rows[rk] = [L.fmul(x, iv) for x in rows[rk]]
        for i in range(rk + 1, len(rows)):
            if rows[i][c]:
                f = rows[i][c]
                rows[i] = [x ^ L.fmul(f, y) for x, y in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows): break
    return rk

# ---------------------------------------------------------------- census
def census(R, blocks, r, list_all=False):
    """blocks = list of (H_q, U_q) pairs (Q-elements).  Return list of E' passing ALL blocks."""
    m, L = R.m, R.L
    HT = []
    for (H, U) in blocks:
        tab = [None] * (1 << m)
        for z in range(1 << m):
            tab[z] = R.qmul(H, R.K_to_Q(z))
        HT.append(tab)
    sols = []
    for basis in rref_subspaces(m, r):
        ok = True
        for bi, (H, U) in enumerate(blocks):
            tab = HT[bi]
            vecs = []
            for b in basis:
                vecs.append(R.K_to_Q(b)); vecs.append(tab[b])
            c, _ = span_contains(vecs, U, L, m)
            if not c: ok = False; break
        if ok: sols.append(tuple(basis))
    return sols

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--ns", type=str, default="5,9,11,13")
    ap.add_argument("--P-max", type=int, default=14)
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--json-out", type=str, default="")
    args = ap.parse_args()

    print(f"### EXACT TSM CENSUS  m={args.m} r={args.r} d={args.d}  (2r={2*args.r}, 2^r={2**args.r})")
    print(f"### predicted formula  P0 = ceil(r(m-r)/(m-2r)) = {-(-args.r*(args.m-args.r)//(args.m-2*args.r)) if args.m>2*args.r else 'N/A (m<=2r)'}")
    print(f"### #r-dim subspaces of GF(2^{args.m}) = {count_subspaces(args.m,args.r)}")
    out = {"m": args.m, "r": args.r, "d": args.d, "rows": []}
    for n in [int(x) for x in args.ns.split(",")]:
        for t in range(args.trials):
            seed = args.seed + 1000 * n + t
            R, inst = make_instance(args.m, n, args.r, args.d, seed)
            H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
            Etrue = canon(inst["E"], args.m)
            # sanity: Phi is a ring iso, rho is the partial Frobenius
            assert R.Phi(R.rmul(inst["e2"], inst["h"])) == R.qmul(R.Phi(inst["e2"]), H0)
            assert R.rho(H0) == R.Phi(R.rho_R(inst["h"])), "rho/Phi compatibility"
            # ---- (A) TRUE system, P = 1..P-max
            Hs, Us = [H0], [U0]
            for q in range(1, args.P_max + 1):
                Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
            counts_true = {}
            for P in range(1, args.P_max + 1):
                if P > n: break
                sols = census(R, list(zip(Hs[:P], Us[:P])), args.r)
                counts_true[P] = len(sols)
                if P == 1:
                    assert any(canon(s, args.m) == Etrue for s in sols), "true E missing!"
            # ---- (B) RANDOM (probe) model
            rng = random.Random(seed ^ 0xabcdef)
            Hr, Ur = [], []
            for q in range(args.P_max + 1):
                Hq = [rng.randrange(1 << n) for _ in range(args.m)]
                aq = [rng.randrange(1 << n) for _ in range(args.r)]
                bq = [rng.randrange(1 << n) for _ in range(args.r)]
                # U_q = sum_k b~_k (alpha_k + H_q beta_k)   [PROBE model: alpha,beta independent per q]
                Uq = [0] * args.m
                for k in range(args.r):
                    t1 = R.K_to_Q(inst["E"][k])
                    inner = [aq[k]] + [0] * (args.m - 1)
                    t2 = R.qmul(Hq, [bq[k]] + [0] * (args.m - 1))
                    Uq = [a ^ b ^ c for a, b, c in zip(Uq, R.qmul(t1, inner), R.qmul(t1, t2))]
                Hr.append(Hq); Ur.append(Uq)
            counts_rand = {}
            for P in range(1, args.P_max + 1):
                if P > n: break
                sols = census(R, list(zip(Hr[:P], Ur[:P])), args.r)
                counts_rand[P] = len(sols)
            pred = args.r * (args.m - args.r) - n * (args.m - 2 * args.r)
            print(f"n={n:3d} trial={t}  predicted #spurious ~ 2^({args.r}(m-{args.r}) - n(m-2r)) = 2^{pred}"
                  f"   TRUE counts by P: {counts_true}")
            print(f"                    {'':18s}RANDOM(probe) counts by P: {counts_rand}")
            out["rows"].append(dict(n=n, trial=t, seed=seed, pred_exponent=pred,
                                    true_counts=counts_true, random_counts=counts_rand))
    if args.json_out:
        json.dump(out, open(args.json_out, "w"), indent=1)

if __name__ == "__main__":
    main()
