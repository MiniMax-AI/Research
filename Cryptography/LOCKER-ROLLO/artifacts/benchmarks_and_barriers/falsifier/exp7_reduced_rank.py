#!/usr/bin/env python3
"""
exp7_reduced_rank.py -- INDEPENDENT CROSS-VALIDATION of the nullity law.

exp2_rank_truth_table.py computes  rank_{GF(2)}[ Laplace-rows (+) Frobenius-rows ]  by stacking
the two row families and calling M4RI.  This script instead SUBSTITUTES the Frobenius relations
explicitly -- every coordinate p_{(Q,mask)} is rewritten as sigma^{-h}(x_{ter(Q),mask}) with
h = (P-1) - max(Q) and ter(Q) = Q + h  (so ter(Q) always contains P-1) -- and then ranks the
resulting system in the FREE coordinates only.  This is a completely different code path
(different matrix, different column space, no Frobenius rows at all), so agreement is a strong
check against implementation bugs.

Claim checked:   rank_{GF(2)}[ Laplace restricted to the Frobenius quotient ]
                     =  m * n * C(P-2, r-1)
and hence  nullity = n * ( 2^r C(P-1,r-1) - m C(P-2,r-1) ).
"""
import sys, os, itertools, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exp1_census import make_instance
import exp2_rank_truth_table as E2

def analyse(m, n, r, P, seed=2024):
    R, inst = make_instance(m, n, r, 3, seed); L = R.L
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Hs, Us = [H0], [U0]
    for q in range(1, P): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    fives = list(itertools.combinations(range(P), r)); fidx = {t: i for i, t in enumerate(fives)}
    ncols = len(fives) * (1 << r)
    def sig_mat(e):
        e %= n; cols = []
        for j in range(n):
            v = 1 << j
            for _ in range(e): v = L.mul(v, v)
            cols.append(v)
        return cols
    free = [(Q, sm) for Q in fives if (P - 1) in Q for sm in range(1 << r)]
    fidx2 = {c: i for i, c in enumerate(free)}
    allc = [(Q, sm) for Q in fives for sm in range(1 << r)]
    def ter(c):
        Q, sm = c; h = (P - 1) - max(Q); return (tuple(x + h for x in Q), sm), h
    nfree = len(free) * n
    TER = {c: ter(c) for c in allc}; SIG = {h: sig_mat(-h) for h in range(P)}
    def reduce_rows(fullrow):
        out = [0] * n
        for ci, c in enumerate(allc):
            a = int(fullrow[ci])
            if not a: continue
            t, h = TER[c]; ti = fidx2[t]; S = SIG[h]
            for i in range(n):
                blk = 0
                for j in range(n):
                    if (L.mul(a, S[j]) >> i) & 1: blk |= 1 << (ti * n + j)
                out[i] ^= blk
        return out
    sixes = [s for s in itertools.combinations(range(P), r + 1) if s[0] == 0]
    rows = []
    for sel in sixes:
        A = np.zeros((m, ncols), dtype=np.uint32)
        for omit in range(r + 1):
            pairs = tuple(x for i, x in enumerate(sel) if i != omit); blk = fidx[pairs]
            for sm in range(1 << r):
                coef = Us[sel[omit]]
                for k in range(r):
                    if (sm >> k) & 1: coef = R.qmul(coef, Hs[pairs[k]])
                for d in range(m): A[d, blk * (1 << r) + sm] = coef[d]
        for d in range(m): rows.extend(reduce_rows(A[d]))
    rk = E2.gf2_rank_py(rows, nfree)
    return len(rows), nfree, rk, m * n

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfgs", type=str, default="5,7,3,5 6,7,3,5 7,5,3,5 7,5,3,6 7,13,3,5 7,13,3,6 9,5,4,6 9,5,4,7 13,7,4,7 11,5,5,7")
    a = ap.parse_args()
    print(f"{'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'rows':>6}{'freecols':>9} | {'rank(subst)':>12}{'m*n*C(P-2,r-1)':>16} {'ok':>3}")
    allok = True
    for cfg in a.cfgs.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        nr, nfree, rk, mn = analyse(m, n, r, P)
        pred = mn * math.comb(P - 2, r - 1)
        ok = (rk == pred); allok &= ok
        print(f"{m:>3}{n:>3}{r:>3}{P:>3} | {nr:>6}{nfree:>9} | {rk:>12}{pred:>16} {'OK' if ok else 'NO':>3}")
    print("INDEPENDENT CROSS-VALIDATION OF THE RANK LAW:", "PASSED" if allok else "FAILED")
    sys.exit(0 if allok else 1)
