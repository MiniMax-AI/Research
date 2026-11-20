#!/usr/bin/env python3
"""
exp5_sigma_equivariance.py -- verification of the structural mechanism behind the
Frobenius redundancy of the TSM support-minor (Plucker) linearisation.

PROPOSITION (sigma-equivariance).  Let  S' subset [P], |S'| = r+1, 0 in S', P-1 not in S',
and let E_{S'} be the Laplace equation block
      E_{S'} :  sum_t sum_{mu in {0,1}^r}  U_{s_t} * ( prod_{i: mu_i=1} H_{q_i} ) * p_{(Q^(t),mu)} = 0
with Q^(t) = S' \ {s_t} = {q_1<...<q_r}.  Then, coefficientwise application of rho
(= squaring in L, i.e. rho(U_q)=U_{q+1}, rho(H_q)=H_{q+1}) maps the coefficient of E_{S'}
sitting at column (Q,mu) to the coefficient of E_{S'+1} sitting at column (Q+1,mu).

Because  P-1 not in S'  implies  Q^(t) subset [P-2]  for every t, the Frobenius relation
p_{(Q+1,mu)} = sigma(p_{(Q,mu)}) is available for every coordinate occurring in E_{S'}.
Hence, MODULO THE FROBENIUS RELATIONS,  E_{S'+1} = rho(E_{S'}) :  every Laplace block is the
rho-image of the block one step to its left.  The C(P-1,r) blocks therefore form a handful of
sigma-chains instead of C(P-1,r) independent constraints -- this is the "extra Frobenius
syzygy" that the random-H_q model does not possess.
"""
import sys, os, itertools, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exp1_census import make_instance

def block(R, Hs, Us, sel, fidx, ncols):
    m = R.m
    A = np.zeros((m, ncols), dtype=np.uint32)
    for omit in range(len(sel)):
        pairs = tuple(x for i, x in enumerate(sel) if i != omit)
        blk = fidx[pairs]
        for sm in range(1 << len(pairs)):
            coef = Us[sel[omit]]
            for k in range(len(pairs)):
                if (sm >> k) & 1: coef = R.qmul(coef, Hs[pairs[k]])
            for d in range(m): A[d, blk * (1 << len(pairs)) + sm] = coef[d]
    return A

def check(m, n, r, P, seed=2024):
    R, inst = make_instance(m, n, r, 3, seed)
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Hs, Us = [H0], [U0]
    for q in range(1, P): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    fives = list(itertools.combinations(range(P), r)); fidx = {t: i for i, t in enumerate(fives)}
    ncols = len(fives) * (1 << r)
    ok = tot = 0
    for S in itertools.combinations(range(P), r + 1):
        if S[0] != 0 or (P - 1) in S: continue
        A, A1 = block(R, Hs, Us, S, fidx, ncols), block(R, Hs, Us, tuple(x + 1 for x in S), fidx, ncols)
        tot += 1; good = True
        for omit in range(r + 1):
            pairs = tuple(x for i, x in enumerate(S) if i != omit)
            blk, blk1 = fidx[pairs], fidx[tuple(x + 1 for x in pairs)]
            for sm in range(1 << r):
                for d in range(m):
                    a = int(A[d, blk * (1 << r) + sm])
                    if R.L.mul(a, a) != int(A1[d, blk1 * (1 << r) + sm]): good = False
        ok += int(good)
    return ok, tot

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfgs", type=str, default="7,13,3,6 7,13,3,7 9,5,4,7 11,5,5,7 13,7,4,7 5,7,3,6")
    a = ap.parse_args()
    allok = True
    for cfg in a.cfgs.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        ok, tot = check(m, n, r, P)
        allok &= (ok == tot)
        print(f"(m,n,r,P)=({m},{n},{r},{P}): rho(coeff of E_S' at (Q,mu)) == coeff of E_(S'+1) at (Q+1,mu)"
              f"  for {ok}/{tot} admissible S'  -> {'VERIFIED' if ok==tot else 'FAILED'}")
    print("PROPOSITION (sigma-equivariance) VERIFIED" if allok else "PROPOSITION FAILED")
    sys.exit(0 if allok else 1)
