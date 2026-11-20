#!/usr/bin/env python3
"""
exp4_negative_control.py -- NEGATIVE CONTROL (uniform random syndrome rejection).

A uniform random syndrome u in R (NOT of the form e1 + e2*h with Supp(e1),Supp(e2) in a
rank-r subspace E) must be REJECTED: the exact census must return 0 valid rank-r supports E.

We also re-run the POSITIVE side (a genuine ROLLO-I instance) in the same process, so the
checker is validated in both directions.
"""
import sys, os, random, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollo_gf2_engine import RingQ, count_subspaces
from exp1_census import make_instance, canon, census

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--n", type=int, default=11)
    ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--seed", type=int, default=8080)
    args = ap.parse_args()
    print(f"### NEGATIVE CONTROL  m={args.m} n={args.n} r={args.r}: random (h,u) must give 0 valid E")
    print(f"### (exhaustive over all {count_subspaces(args.m,args.r)} rank-{args.r} subspaces of GF(2^{args.m}))")
    allok = True
    for t in range(args.trials):
        seed = args.seed + t
        # ---- POSITIVE side: genuine instance, census must be non-empty and contain the true E
        R, inst = make_instance(args.m, args.n, args.r, args.d, seed)
        H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
        sols = census(R, [(H0, U0)], args.r)
        Et = canon(inst["E"], args.m)
        pos_ok = any(canon(s, args.m) == Et for s in sols)
        # ---- NEGATIVE side: uniform random (h,u)
        rng = random.Random(seed ^ 0xBADF00D)
        h = [rng.randrange(1 << args.m) for _ in range(args.n)]
        u = [rng.randrange(1 << args.m) for _ in range(args.n)]
        Hr, Ur = R.Phi(h), R.Phi(u)
        nsol = len(census(R, [(Hr, Ur)], args.r))
        neg_ok = (nsol == 0)
        allok &= pos_ok and neg_ok
        print(f"trial {t}: POSITIVE genuine instance -> {len(sols)} valid E (contains true E: {pos_ok}) | "
              f"NEGATIVE random (h,u) -> {nsol} valid E (expected 0: {neg_ok})")
    print("NEGATIVE CONTROL PASSED" if allok else "NEGATIVE CONTROL FAILED")
    return 0 if allok else 1

if __name__ == "__main__":
    sys.exit(main())
