#!/usr/bin/env python3
"""
exp3_pubkey_census.py -- EXACT exhaustive census for the HOMOGENEOUS PUBLIC-KEY problem
    h = x^{-1} y,  Supp(x), Supp(y) subset F,  dim_{GF(2)} F = d.
Transposed (TSM) form:   sum_k f~_k (alpha_k H_0 + beta_k) = 0,  alpha_k, beta_k in L.
Writing W_0 = <1,H_0>_L (2-dimensional over L), the condition is
    dim_L( F~L + H_0 F~L ) < 2d          (i.e. the 2d vectors f~_k, H_0 f~_k are L-DEPENDENT).
By Lemma REDUND the q-th conjugate condition  dim_L(F~L + H_q F~L) < 2d  is the rho^q-image of
the q=0 one, hence adds nothing: the solution set is P-independent.

We enumerate EVERY d-dimensional GF(2)-subspace F' of K and count exactly.
"""
import sys, os, json, random, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollo_gf2_engine import RingQ, rref_subspaces, count_subspaces
from exp1_census import make_instance, canon, rank_L

def census_pubkey(R, Hs, d):
    """count d-dim F' with dim_L(F'L + H_q F'L) < 2d for EVERY q."""
    m, L = R.m, R.L
    HT = []
    for H in Hs:
        tab = [None] * (1 << m)
        for z in range(1 << m):
            tab[z] = R.qmul(H, R.K_to_Q(z))
        HT.append(tab)
    sols = []
    for basis in rref_subspaces(m, d):
        ok = True
        for tab in HT:
            vecs = []
            for b in basis:
                vecs.append(R.K_to_Q(b)); vecs.append(tab[b])
            if rank_L(vecs, L, m) == 2 * d: ok = False; break
        if ok: sols.append(tuple(basis))
    return sols

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--ns", type=str, default="5,7,9,11,13")
    ap.add_argument("--P-max", type=int, default=14)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--seed", type=int, default=555)
    ap.add_argument("--json-out", type=str, default="")
    args = ap.parse_args()
    print(f"### PUBLIC-KEY (homogeneous) CENSUS  m={args.m} d={args.d}  (2d={2*args.d}, 2^d={2**args.d})")
    print(f"### #d-dim subspaces of GF(2^{args.m}) = {count_subspaces(args.m,args.d)}")
    print(f"### predicted #solutions ~ C(m,d)_2 * 2^(-n(m-2d+1))")
    out = {"m": args.m, "d": args.d, "rows": []}
    for n in [int(x) for x in args.ns.split(",")]:
        for t in range(args.trials):
            seed = args.seed + 1000 * n + t
            R, inst = make_instance(args.m, n, args.r, args.d, seed)
            H0 = R.Phi(inst["h"])
            Ftrue = canon(inst["F"], args.m)
            Hs = [H0]
            for q in range(1, args.P_max + 1):
                Hs.append(R.rho(Hs[-1]))
            counts, orbits = {}, {}
            for P in range(1, min(args.P_max, n) + 1):
                sols = census_pubkey(R, Hs[:P], args.d)
                counts[P] = len(sols)
                if P == 1:
                    assert any(canon(s, args.m) == Ftrue for s in sols), "true F missing!"
                    # K^*-orbit structure: cF is a solution whenever F is
                    base = set(canon(s, args.m) for s in sols)
                    orb = set()
                    for s in base:
                        for c in range(1, 1 << args.m):
                            orb.add(canon([R.K.mul(c, x) for x in s], args.m))
                    orbits[P] = (len(base), len(orb), len(base) * ((1 << args.m) - 1) == len(orb))
            import math
            pred = math.log2(count_subspaces(args.m, args.d)) - n * (args.m - 2 * args.d + 1)
            print(f"n={n:3d} trial={t}  log2 pred #extra = {pred:+.2f}  counts by P: {counts}")
            print(f"          K^*-orbit check (|sols|, |K^*-orbit of true F|, closed): {orbits.get(1)}")
            out["rows"].append(dict(n=n, trial=t, seed=seed, log2_pred=pred, counts=counts, orbit=orbits.get(1)))
    if args.json_out: json.dump(out, open(args.json_out, "w"), indent=1)

if __name__ == "__main__":
    main()
