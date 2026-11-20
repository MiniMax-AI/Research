#!/usr/bin/env python3
"""
exp10_degree2_obstruction.py -- the Degree-2 (Plucker) stage is underdetermined by construction.

After the Degree-1 stage the solution is only known modulo a GF(2)-kernel of dimension
    D_mask = n [ 2^r C(P-1,r-1) - m C(P-2,r-1) ]      (mask parameterisation, parameterized_tsm/tsm.cpp)
    D_simp = m n C(P-2,r-2)                           ((SIMP) parameterisation, exp8)
The only remaining algebraic constraints are the Plucker relations of the relevant Grassmannian
Gr(r,N), N = 2P (mask) resp. N = P (SIMP); their number is
    #rel = C(N, r-1) * C(N-r+1, r+1),        and this is 0 unless N >= 2r.
Linearising them in the degree-2 monomials of the kernel coordinates needs C(D+1,2) monomials.

This script (a) tabulates the forced counts for the official parameter sets and (b) verifies the
counting at toy scale by actually building the Degree-2 system on the TRUE (SIMP) kernel.
"""
import sys, os, itertools, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exp1_census import make_instance
from rollo_gf2_engine import RingQ
from exp8_simple_linearization import qinv, qadd
RingQ.qinv = qinv

SETS = [("ROLLO-I-128 (2020)", 67, 83, 7, 8), ("ROLLO-I-192 (2020)", 79, 97, 8, 8),
        ("ROLLO-I-256 (2020)", 97, 113, 9, 9), ("ROLLO-II-128 (2020)", 83, 189, 7, 8),
        ("ROLLO-II-192 (2020)", 97, 193, 8, 8), ("ROLLO-II-256 (2020)", 97, 211, 8, 9),
        ("new2ROLLO-I-128", 73, 83, 7, 8), ("new2ROLLO-I-192", 89, 97, 8, 8),
        ("new2ROLLO-I-256", 103, 113, 9, 9), ("ROLLO-I-256 (2019)", 113, 67, 7, 7),
        ("ROLLO-II-256 (2019)", 127, 157, 7, 9), ("LOCKER-I (2017)", 71, 83, 5, 7)]

def fmt(x):
    if x == 0: return "0"
    e = int(math.log10(x)); return f"{x/10**e:.1f}e{e}"

def table():
    print("DEGREE-2 (Plucker) LINEARISATION at the planned P")
    print("  mask parameterisation Gr(r,2P):  #rel = C(2P,r-1)*C(2P-r+1,r+1)")
    print("  SIMP parameterisation Gr(r,P) :  #rel = C(P,r-1)*C(P-r+1,r+1)")
    print(f"{'parameter set':<21}{'m':>4}{'n':>5}{'r':>3}{'P':>3} | {'D_mask(bits)':>13}{'#monomials':>12}"
          f"{'#rel Gr(r,2P)':>14}{'deficit':>10} | {'#rel Gr(r,P)':>13}")
    for (nm, m, n, r, P) in SETS:
        D = max(n * ((1 << r) * math.comb(P - 1, r - 1) - m * math.comb(P - 2, r - 1)), 1)
        mono = D * (D + 1) // 2
        rm = math.comb(2 * P, r - 1) * math.comb(2 * P - r + 1, r + 1) if 2 * P - r + 1 >= r + 1 else 0
        rs = math.comb(P, r - 1) * math.comb(P - r + 1, r + 1) if P - r + 1 >= r + 1 else 0
        print(f"{nm:<21}{m:>4}{n:>5}{r:>3}{P:>3} | {fmt(D):>13}{fmt(mono):>12}{fmt(rm):>14}"
              f"{(f'{mono/rm:.0f}x' if rm else 'inf'):>10} | {fmt(rs):>13}")
    print("\nSmallest P for which ANY Plucker relation exists (SIMP Gr(r,P) needs P >= 2r):")
    for (nm, m, n, r, P) in SETS[:6]:
        print(f"  {nm:<21} r={r}: P >= 2r = {2*r}")

def toy(m, n, r, P, seed=2024):
    """actually build the Degree-2 system on the TRUE (SIMP) kernel and rank it"""
    import exp8_simple_linearization as S
    R, inst = make_instance(m, n, r, 3, seed)
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    Us = [U0]
    for q in range(1, P): Us.append(R.rho(Us[-1]))
    fives = list(itertools.combinations(range(P), r))
    free = [Q for Q in fives if (P - 1) in Q]; fidx = {Q: i for i, Q in enumerate(free)}
    mn = m * n
    TER = {Q: ((P - 1) - max(Q), tuple(x + ((P - 1) - max(Q)) for x in Q)) for Q in fives}
    RH = {}
    for h in range(P):
        e = (-h) % n; lst = []
        for j in range(mn):
            t, i = divmod(j, n); v = [0] * m; v[t] = 1 << i
            for _ in range(e): v = R.rho(v)
            lst.append(v)
        RH[h] = lst
    memo = {}
    def blockmat(h, s):
        if (h, s) in memo: return memo[(h, s)]
        coef = Us[s]; lst = RH[h]; out = [0] * mn
        for j in range(mn):
            prod = R.qmul(coef, lst[j])
            for b in range(mn):
                t, i = divmod(b, n)
                if (prod[t] >> i) & 1: out[b] |= 1 << j
        memo[(h, s)] = out; return out
    nfree = len(free) * mn
    rows = []
    for S_ in itertools.combinations(range(P), r + 1):
        if S_[0] != 0: continue
        acc = [0] * mn
        for t in range(r + 1):
            Q = tuple(x for i, x in enumerate(S_) if i != t)
            h, T = TER[Q]; ti = fidx[T]; M = blockmat(h, S_[t])
            for b in range(mn): acc[b] ^= M[b] << (ti * mn)
        rows.extend(acc)
    Mx = [x for x in rows]; piv = []; rk = 0
    for c in range(nfree):
        p = -1
        for i in range(rk, len(Mx)):
            if (Mx[i] >> c) & 1: p = i; break
        if p < 0: continue
        Mx[rk], Mx[p] = Mx[p], Mx[rk]
        for i in range(len(Mx)):
            if i != rk and ((Mx[i] >> c) & 1): Mx[i] ^= Mx[rk]
        piv.append(c); rk += 1
        if rk == len(Mx): break
    D = nfree - rk
    fr = [c for c in range(nfree) if c not in set(piv)]
    basis = []
    for fc in fr:
        v = 1 << fc
        for i in range(len(piv) - 1, -1, -1):
            if (Mx[i] >> fc) & 1: v |= 1 << piv[i]
        basis.append(v)
    # linear forms: bit b of p_Q  as a D-bit form
    forms = {}
    for Q in fives:
        h, T = TER[Q]; ti = fidx[T]
        for b in range(mn):
            f = 0
            for j in range(mn):
                t2, i2 = divmod(j, n); v = [0] * m; v[t2] = 1 << i2
                for _ in range((-h) % n): v = R.rho(v)
                tb, ib = divmod(b, n)
                if (v[tb] >> ib) & 1:
                    for d in range(D):
                        if (basis[d] >> (ti * mn + j)) & 1: f ^= 1 << d
            forms[(Q, b)] = f
    nm = D * (D + 1) // 2
    idx = {}; c = 0
    for i in range(D):
        for j in range(i, D): idx[(i, j)] = c; c += 1
    Is = list(itertools.combinations(range(P), r - 1)); Js = list(itertools.combinations(range(P), r + 1))
    qrows = []
    for I in Is:
        for J in Js:
            if set(I) & set(J): continue
            acc = 0
            for jk in J:
                A = tuple(sorted(I + (jk,))); B = tuple(x for x in J if x != jk)
                for b1 in range(mn):
                    f1 = forms[(A, b1)]
                    if not f1: continue
                    for b2 in range(mn):
                        f2 = forms[(B, b2)]
                        if not f2: continue
                        a1 = np.unpackbits(np.frombuffer(f1.to_bytes((D + 7) // 8, "little"), dtype=np.uint8), bitorder="little")[:D]
                        a2 = np.unpackbits(np.frombuffer(f2.to_bytes((D + 7) // 8, "little"), dtype=np.uint8), bitorder="little")[:D]
                        pr = np.outer(a1, a2) ^ np.outer(a2, a1)
                        ii, jj = np.nonzero(np.triu(pr, 1))
                        for a, b in zip(ii, jj): acc ^= 1 << idx[(int(a), int(b))]
                        dg = np.nonzero(a1 & a2)[0]
                        for a in dg: acc ^= 1 << idx[(int(a), int(a))]
            if acc: qrows.append(acc)
    import exp2_rank_truth_table as E2
    rk2 = E2.gf2_rank_py(qrows, nm) if qrows else 0
    return D, nm, len(qrows), rk2, nm - rk2

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy", type=str, default="3,5,2,4 5,7,2,5")
    a = ap.parse_args()
    table()
    print("\n=== toy-scale verification: actually build the Degree-2 system on the TRUE (SIMP) kernel ===")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'D':>4}{'#monomials':>11}{'#relations':>11}{'rank':>7}{'NULLITY':>8}")
    for cfg in a.toy.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        D, nm, nr, rk, nl = toy(m, n, r, P)
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {D:>4}{nm:>11}{nr:>11}{rk:>7}{nl:>8}")
    print("  => the Degree-2 linearisation is underdetermined already at toy size.")
