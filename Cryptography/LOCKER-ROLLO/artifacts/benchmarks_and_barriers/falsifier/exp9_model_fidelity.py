#!/usr/bin/env python3
"""
exp9_model_fidelity.py -- two checks that the linear-algebra model used here
faithfully represents the reference solver parameterized_tsm/tsm.cpp.

CHECK A (no missing equations).  We build the Laplace stage from the (r+1)-subsets
S' that CONTAIN 0 (that is what tsm.cpp does: `for (auto &c : sixes_all) if (c[0]==0) ...`).
Does using ALL C(P,r+1) subsets add rank?  Measured answer: NO, exactly zero gain.

CHECK B (same Frobenius stage).  tsm.cpp's Frobenius stage imposes the L-LINEAR condition
        sum_f x_f * ( b_f[theta c] + sigma(b_f[c]) ) = 0            for every relation c,
on the L-kernel coordinates x of the Stage-1 system, whereas the TRUE Frobenius relation is
the GF(2)-linear
        sum_f x_f b_f[theta c]  =  sigma( sum_f x_f b_f[c] )  =  sum_f sigma(x_f) sigma(b_f[c]).
(The two agree iff x is sigma-fixed.)  Measured answer: in every m <= 2^r configuration tested
the two give IDENTICAL nullities, so the model represents the reference solver in the regime
m <= 2^r.
"""
import sys, os, itertools, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exp1_census import make_instance
import exp2_rank_truth_table as E2

def _laplace(R, Hs, Us, P, r, allsub):
    m = R.m
    fives = list(itertools.combinations(range(P), r)); fidx = {t: i for i, t in enumerate(fives)}
    ncols = len(fives) * (1 << r)
    sixes = [s for s in itertools.combinations(range(P), r + 1) if (s[0] == 0 or allsub)]
    A = np.zeros((len(sixes) * m, ncols), dtype=np.uint32)
    for eq, sel in enumerate(sixes):
        for omit in range(r + 1):
            pairs = tuple(x for i, x in enumerate(sel) if i != omit); blk = fidx[pairs]
            for sm in range(1 << r):
                coef = Us[sel[omit]]
                for k in range(r):
                    if (sm >> k) & 1: coef = R.qmul(coef, Hs[pairs[k]])
                for d in range(m): A[eq * m + d, blk * (1 << r) + sm] = coef[d]
    return fives, fidx, ncols, A

def check_A(m, n, r, P, seed=2024):
    R, inst = make_instance(m, n, r, 3, seed)
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Hs, Us = [H0], [U0]
    for q in range(1, P): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    POW, BLK = E2.mul_tables(R); SIG = E2.sigma_blk(R)
    out = []
    for allsub in (False, True):
        fives, fidx, ncols, A = _laplace(R, Hs, Us, P, r, allsub)
        frows = []
        for Q in itertools.combinations(range(P - 1), r):
            Q1 = tuple(x + 1 for x in Q)
            if Q1 not in fidx: continue
            c0 = fidx[Q] * (1 << r); c1 = fidx[Q1] * (1 << r)
            for sm in range(1 << r):
                for i in range(n):
                    v = (1 << i) << ((c1 + sm) * n); v ^= SIG[i] << ((c0 + sm) * n); frows.append(v)
        rows = E2.flatten(A, BLK, n) + frows
        out.append(ncols * n - E2.gf2_rank_py(rows, ncols * n))
    return out

def check_B(m, n, r, P, seed=2024):
    """nullity under tsm.cpp's L-linear Frobenius condition (Stage-1 L-kernel first)."""
    R, inst = make_instance(m, n, r, 3, seed); L = R.L
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Hs, Us = [H0], [U0]
    for q in range(1, P): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    fives, fidx, ncols, A = _laplace(R, Hs, Us, P, r, False)
    M = [list(map(int, A[i])) for i in range(A.shape[0])]
    piv = []; rk = 0
    for c in range(ncols):
        p = -1
        for i in range(rk, len(M)):
            if M[i][c]: p = i; break
        if p < 0: continue
        M[rk], M[p] = M[p], M[rk]
        iv = L.finv(M[rk][c]); M[rk] = [L.fmul(x, iv) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [x ^ L.fmul(f, y) for x, y in zip(M[i], M[rk])]
        piv.append(c); rk += 1
        if rk == len(M): break
    D = ncols - rk
    free = [c for c in range(ncols) if c not in set(piv)]
    basis = []
    for fc in free:
        v = [0] * ncols; v[fc] = 1
        for i in range(len(piv) - 1, -1, -1):
            c = piv[i]; s = 0
            for j in range(c + 1, ncols):
                if v[j] and M[i][j]: s ^= L.fmul(M[i][j], v[j])
            v[c] = s
        basis.append(v)
    rel = []
    for ci, Q in enumerate(fives):
        if Q[-1] == P - 1: continue
        si = fidx[tuple(x + 1 for x in Q)]
        for sm in range(1 << r): rel.append((ci * (1 << r) + sm, si * (1 << r) + sm))
    rows = []
    for (c0, c1) in rel:
        for i in range(n):
            row = 0
            for f in range(D):
                a = basis[f][c1]; b = L.fmul(basis[f][c0], basis[f][c0]); v = a ^ b
                if not v: continue
                for j in range(n):
                    if (L.mul(v, 1 << j) >> i) & 1: row |= 1 << (f * n + j)
            rows.append(row)
    return D, D * n - E2.gf2_rank_py(rows, D * n)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=str, default="7,13,3,6 7,13,3,7 6,7,3,6 17,19,3,6 11,13,3,7")
    ap.add_argument("--B", type=str, default="7,13,3,5 7,13,3,6 7,13,3,7 6,7,3,6")
    a = ap.parse_args()
    print("=== CHECK A: does using ALL C(P,r+1) subsets add rank? ===")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'only-0 null':>12}{'ALL null':>10}{'gain':>7}")
    okA = True
    for cfg in a.A.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        n0, n1 = check_A(m, n, r, P)
        okA &= (n0 == n1)
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {n0:>12}{n1:>10}{n0-n1:>7}")
    print(f"  CHECK A {'PASSED (no missing equations)' if okA else 'FAILED'}")
    print("\n=== CHECK B: tsm.cpp's L-linear Frobenius condition vs the true GF(2)-linear one ===")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'m<=2^r':>7}{'Stage1 D':>10}{'tsm-frob null':>15}{'pred null':>13}{'same':>6}")
    okB = True
    for cfg in a.B.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        D, nT = check_B(m, n, r, P)
        rep = n * ((1 << r) * math.comb(P - 1, r - 1) - m * math.comb(P - 2, r - 1))
        same = (nT == rep); okB &= same
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {str(m<=2**r):>7}{D:>10}{nT:>15}{rep:>13}{str(same):>6}")
    print(f"  CHECK B {'PASSED (identical nullities for m<=2^r)' if okB else 'FAILED'}")
    sys.exit(0 if (okA and okB) else 1)
