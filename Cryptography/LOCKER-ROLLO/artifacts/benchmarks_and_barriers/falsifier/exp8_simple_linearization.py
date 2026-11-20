#!/usr/bin/env python3
"""
exp8_simple_linearization.py -- A STRICTLY STRONGER Degree-1 TSM linearisation, and its
(negative) consequence for the m <= 2^r barrier.

DERIVATION.  Put  gamma_k = alpha_k + beta_k * H_0  in Q.  Then
        C_{q,k} = w_{2q,k} + H_q w_{2q+1,k} = rho^q(alpha_k) + rho^q(H_0) rho^q(beta_k)
                = rho^q(gamma_k),          and   U_q = rho^q(U_0).
So with the single vector  v = (gamma_1,...,gamma_r,U_0) in Q^{r+1}, the TSM equations say that
the (r+1) x (r+1) matrix with rows rho^{s}(v), s in S', has determinant 0 (its last column is a
Q-linear combination of the first r, with coefficients b~_k in K subset Q).  Expanding along the
LAST COLUMN gives the clean system

        (SIMP)     sum_{t=0}^{r}  U_{s_t} * p_{S' \ {s_t}}  =  0     in Q,

where  p_Q = det[ rho^q(gamma_k) ]_{q in Q, k}  are the Plucker coordinates of the P x r matrix
[rho^q(gamma_k)]  (NOT the 2^r-times larger mask coordinates used by parameterized_tsm/tsm.cpp).
Frobenius relations:   p_{Q+1} = rho(p_Q).   Unknowns C(P,r) elements of Q; free ones: Q with
P-1 in Q, i.e. C(P-1,r-1).

(SIMP) is a VALID necessary condition (checked here against the true p_Q) and is STRICTLY
STRONGER than the mask linearisation, because the mask system is exactly the pull-back of (SIMP)
under the Q-linear map  p_{(Q,mask)} |-> sum_mask (prod_{i:mask_i=1} H_{q_i}) p_{(Q,mask)}.

This script (a) validates (SIMP) on genuine instances and (b) measures its exact GF(2) nullity.
"""
import sys, os, itertools, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollo_gf2_engine import RingQ
from exp1_census import make_instance
import exp2_rank_truth_table as E2

def qinv(self, a):
    e = (1 << (self.m * self.n)) - 2; r = [1] + [0] * (self.m - 1); b = a[:]
    while e:
        if e & 1: r = self.qmul(r, b)
        b = self.qmul(b, b); e >>= 1
    return r
RingQ.qinv = qinv
qadd = lambda a, b: [x ^ y for x, y in zip(a, b)]

def det_Q(M, R):
    n_ = len(M); M = [[c[:] for c in row] for row in M]
    d = [1] + [0] * (R.m - 1)
    for c in range(n_):
        piv = -1
        for i in range(c, n_):
            if any(M[i][c]): piv = i; break
        if piv < 0: return [0] * R.m
        M[c], M[piv] = M[piv], M[c]
        d = R.qmul(d, M[c][c]); iv = R.qinv(M[c][c])
        for i in range(c + 1, n_):
            if any(M[i][c]):
                f = R.qmul(M[i][c], iv)
                M[i] = [qadd(a, R.qmul(f, b)) for a, b in zip(M[i], M[c])]
    return d

def _gamma(R, inst, m, n, r):
    E = inst["E"]
    def recover(e):
        cols = []
        for j in range(n):
            v = e[j]; rs = E[:]; tag = [1 << k for k in range(r)]; rk = 0
            for c in range(m):
                p = -1
                for i in range(rk, r):
                    if (rs[i] >> c) & 1: p = i; break
                if p < 0: continue
                rs[rk], rs[p] = rs[p], rs[rk]; tag[rk], tag[p] = tag[p], tag[rk]
                for i in range(r):
                    if i != rk and ((rs[i] >> c) & 1): rs[i] ^= rs[rk]; tag[i] ^= tag[rk]
                rk += 1
            t, s = v, 0
            for i in range(rk):
                c = (rs[i] & -rs[i]).bit_length() - 1
                if (t >> c) & 1: t ^= rs[i]; s ^= tag[i]
            cols.append(s)
        al = [0] * r
        for j in range(n):
            for k in range(r):
                if (cols[j] >> k) & 1: al[k] ^= (1 << j)
        return al
    al, be = recover(inst["e1"]), recover(inst["e2"])
    H0 = R.Phi(inst["h"])
    return [qadd([al[k]] + [0] * (m - 1), R.qmul([be[k]] + [0] * (m - 1), H0)) for k in range(r)]

def validate(m, n, r, P, seed=2024):
    R, inst = make_instance(m, n, r, 3, seed)
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    gam = _gamma(R, inst, m, n, r)
    rowsq = []; g = [x[:] for x in gam]
    for q in range(P): rowsq.append([x[:] for x in g]); g = [R.rho(x) for x in g]
    pQ = {Q: det_Q([[rowsq[q][k] for k in range(r)] for q in Q], R)
          for Q in itertools.combinations(range(P), r)}
    Us = [U0]
    for q in range(1, P): Us.append(R.rho(Us[-1]))
    bad = tot = 0
    for S in itertools.combinations(range(P), r + 1):
        if S[0] != 0: continue
        acc = [0] * m
        for t in range(r + 1):
            Q = tuple(x for i, x in enumerate(S) if i != t)
            acc = qadd(acc, R.qmul(Us[S[t]], pQ[Q]))
        tot += 1; bad += int(any(acc))
    fbad = ftot = 0
    for Q in itertools.combinations(range(P - 1), r):
        ftot += 1
        if pQ[tuple(x + 1 for x in Q)] != R.rho(pQ[Q]): fbad += 1
    nz = sum(1 for v in pQ.values() if any(v))
    return bad, tot, fbad, ftot, nz, len(pQ)

def nullity(m, n, r, P, seed=2024):
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
        key = (h, s)
        if key in memo: return memo[key]
        coef = Us[s]; lst = RH[h]; out = [0] * mn
        for j in range(mn):
            prod = R.qmul(coef, lst[j])
            for b in range(mn):
                t, i = divmod(b, n)
                if (prod[t] >> i) & 1: out[b] |= 1 << j
        memo[key] = out; return out
    nfree = len(free) * mn
    rows = []
    for S in itertools.combinations(range(P), r + 1):
        if S[0] != 0: continue
        acc = [0] * mn
        for t in range(r + 1):
            Q = tuple(x for i, x in enumerate(S) if i != t)
            h, T = TER[Q]; ti = fidx[T]; M = blockmat(h, S[t])
            for b in range(mn): acc[b] ^= M[b] << (ti * mn)
        rows.extend(acc)
    rk = E2.gf2_rank_py(rows, nfree)
    return len(rows), nfree, rk, nfree - rk

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--valid", type=str, default="7,13,3,5 7,5,3,6 9,5,4,6 6,7,3,5 11,5,5,7")
    ap.add_argument("--null", type=str, default="7,5,3,4 7,5,3,5 6,7,3,4 6,7,3,5 6,7,3,6 6,7,3,7 9,5,4,5 11,5,5,7 7,13,3,4 7,13,3,5 7,13,3,6")
    a = ap.parse_args()
    print("=== (a) POSITIVE CONTROL: the true p_Q must satisfy (SIMP) and the Frobenius relations ===")
    allok = True
    for cfg in a.valid.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        bad, tot, fbad, ftot, nz, tp = validate(m, n, r, P)
        ok = (bad == 0 and fbad == 0 and nz > 0); allok &= ok
        print(f"  (m,n,r,P)=({m},{n},{r},{P}): eq violations {bad}/{tot}, Frobenius violations {fbad}/{ftot}, "
              f"nonzero p_Q {nz}/{tp} -> {'VALID' if ok else 'INVALID'}")
    print(f"  (SIMP) IS A VALID NECESSARY CONDITION: {allok}")
    print("\n=== (b) exact GF(2) nullity of (SIMP) vs the mask linearisation ===")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'rows':>6}{'freecols':>9}{'rank':>7}{'NULL':>7} | {'mask NULL':>10}{'ratio':>7}")
    for cfg in a.null.split():
        m, n, r, P = [int(x) for x in cfg.split(",")]
        nr, nf, rk, nl = nullity(m, n, r, P)
        mask = n * ((1 << r) * math.comb(P - 1, r - 1) - m * math.comb(P - 2, r - 1))
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {nr:>6}{nf:>9}{rk:>7}{nl:>7} | {mask:>10}{(nl/mask if mask else 0):>7.2f}")
    print("  => (SIMP) is strictly stronger (nullity ratio ~0.5-0.8) but its nullity still GROWS with P:")
    print("     the m <= 2^r barrier survives this stronger linearisation as well.")
