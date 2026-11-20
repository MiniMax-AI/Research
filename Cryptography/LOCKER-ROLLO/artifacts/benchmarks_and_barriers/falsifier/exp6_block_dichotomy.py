#!/usr/bin/env python3
"""
exp6_block_dichotomy.py -- STRUCTURAL DICHOTOMY of the TSM support-minor (Plucker) blocks.

THEOREM (block dichotomy).  In the TRUE Frobenius system (H_q = rho^q(H_0), U_q = rho^q(U_0)),
after the Frobenius relations  p_{(Q+1,mask)} = sigma(p_{(Q,mask)})  are imposed, the Laplace
equation blocks  E_{S'},  S' in C([P], r+1), 0 in S',  split EXACTLY as

    E_{S'} is INDEPENDENT   <=>   {0,1} subset S'          (there are C(P-2,r-1) of them)
    E_{S'} is REDUNDANT     <=>   1 not in S'              (there are C(P-2,r)   of them)

and every independent block contributes EXACTLY m*n bits of rank.  Consequently

    rank_{GF(2)}[ Frobenius ]              = n * 2^r * C(P-1, r)
    rank_{GF(2)}[ Laplace (+) Frobenius ]  = n * ( 2^r C(P-1,r) + m C(P-2,r-1) )
    nullity                                = n * ( 2^r C(P-1,r-1) - m C(P-2,r-1) )

Mechanism: writing  gamma_k = alpha_k + beta_k H_0  one has
C_{q,k} = rho^q(gamma_k) and U_q = rho^q(U_0), i.e. the q-th TSM equation is literally the
rho^q-image of the q=0 one; and modulo the Frobenius relations  E_{S'+1} = rho(E_{S'}).
Only the blocks that simultaneously see q = 0 and q = 1 escape this collapse.
"""
import sys, os, itertools, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exp1_census import make_instance
import exp2_rank_truth_table as E2

def probe(m, n, r, P, seed=2024):
    R, inst = make_instance(m, n, r, 3, seed)
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Hs, Us = [H0], [U0]
    for q in range(1, P): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    POW, BLK = E2.mul_tables(R); SIG = E2.sigma_blk(R)
    fives = list(itertools.combinations(range(P), r)); fidx = {t: i for i, t in enumerate(fives)}
    ncols = len(fives) * (1 << r); ncols_gf2 = ncols * n
    def block(sel):
        A = np.zeros((m, ncols), dtype=np.uint32)
        for omit in range(r + 1):
            pairs = tuple(x for i, x in enumerate(sel) if i != omit); blk = fidx[pairs]
            for sm in range(1 << r):
                coef = Us[sel[omit]]
                for k in range(r):
                    if (sm >> k) & 1: coef = R.qmul(coef, Hs[pairs[k]])
                for d in range(m): A[d, blk * (1 << r) + sm] = coef[d]
        return E2.flatten(A, BLK, n)
    frows = []
    for Q in itertools.combinations(range(P - 1), r):
        Q1 = tuple(x + 1 for x in Q)
        if Q1 not in fidx: continue
        c0 = fidx[Q] * (1 << r); c1 = fidx[Q1] * (1 << r)
        for sm in range(1 << r):
            for i in range(n):
                v = (1 << i) << ((c1 + sm) * n); v ^= SIG[i] << ((c0 + sm) * n); frows.append(v)
    sixes = [s for s in itertools.combinations(range(P), r + 1) if s[0] == 0]
    cur = list(frows); rk = E2.gf2_rank(cur, ncols_gf2, "f"); contrib = {}
    for s in sixes:
        cur = cur + block(s); nrk = E2.gf2_rank(cur, ncols_gf2, "i"); contrib[s] = nrk - rk; rk = nrk
    return sixes, contrib, rk, m * n

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfgs", type=str,
        default="5,7,3,4 5,7,3,5 5,7,3,6 6,7,3,5 6,7,3,6 7,5,3,5 7,5,3,6 7,5,3,7 7,13,3,5 7,13,3,6 "
                "7,13,3,7 9,5,4,6 9,5,4,7 11,5,4,7 13,7,4,7 9,5,5,7 11,5,5,7 17,5,5,7 9,11,6,8")
    a = ap.parse_args()
    print(f"{'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'#indep':>7}{'C(P-2,r-1)':>11} | {'indep<=>1inS':>12} | "
          f"{'each=mn':>8} | {'rank':>7}{'formula':>9} {'ok':>3}")
    allok = True
    for cfg in a.cfgs.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        sixes, contrib, rk, mn = probe(m, n, r, P)
        indep = [s for s in sixes if contrib[s] > 0]
        dich = all(((1 in s) == (contrib[s] > 0)) for s in sixes)
        each = all((contrib[s] == mn) for s in indep)
        pred = n * ((1 << r) * math.comb(P - 1, r) + m * math.comb(P - 2, r - 1))
        ok = (len(indep) == math.comb(P - 2, r - 1)) and dich and each and rk == pred
        allok &= ok
        print(f"{m:>3}{n:>3}{r:>3}{P:>3} | {len(indep):>7}{math.comb(P-2,r-1):>11} | {str(dich):>12} | "
              f"{str(each):>8} | {rk:>7}{pred:>9} {'OK' if ok else 'NO':>3}")
    print("STRUCTURAL DICHOTOMY VERIFIED ON ALL CONFIGS:", allok)
    sys.exit(0 if allok else 1)
