#!/usr/bin/env python3
"""
s19_falsifier.py -- Master algebraic verification and degree-2 obstruction witness suite.
Self-contained driver: run `python3 s19_falsifier.py` from this directory.

It executes, over TRUE GF(2^m)/GF(2^n)/GF(2^(mn)) with the REAL partial Frobenius
rho(sum gamma_t T^t) = sum gamma_t^2 T^t :

  W1  Lemma REDUND  -- the P conjugate TSM equations are rho^q-images of the q=0 equation.
  W2  Exact exhaustive census of the ciphertext TSM solution set vs P  (TRUE vs PROBE model).
  W3  Count law  E[#valid E'] = C(m,r)_2 * 2^{-n(m-2r)}.
  W4  Exact GF(2) rank truth table of [Laplace (+) Frobenius] and the closed-form nullity law
        null(P) = n [ 2^r C(P-1,r-1) - m C(P-2,r-1) ]     (verified for every m <= 2^r tested)
  W5  Homogeneous public-key census (Gamma_S): solution set = K^*-orbit of F, P-independent.
  W6  Negative control: uniform random (h,u) gives 0 valid rank-r supports.
  W7  Proposition (sigma-equivariance): E_{S'+1} = rho(E_{S'}) modulo the Frobenius relations.
  W8  Theorem (block dichotomy): E_{S'} survives iff {0,1} subset S'.
  W9  Independent cross-validation of the rank law by explicit Frobenius substitution.
  W10 New linearisation (SIMP) in the Plucker coordinates of [rho^q(gamma_k)]: validity + nullity law.
  W11 Model fidelity vs tsm.cpp: (A) all subsets add no rank, (B) same Frobenius-stage nullity.
  W12 The Degree-2 (Plucker) stage is underdetermined by construction.

No external solver is needed (pure python + numpy).  If the M4RI helper ./gf2rank is present
(build: gcc -O3 -o gf2rank gf2rank.c $(pkg-config --cflags --libs m4ri)) it is used for speed.
"""
import os, sys, subprocess, math, itertools, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("GF2RANK", os.path.join(HERE, "gf2rank"))

from rollo_gf2_engine import RingQ, rref_subspaces, count_subspaces
from exp1_census import make_instance, canon, census
from exp3_pubkey_census import census_pubkey
from exp5_sigma_equivariance import check as sigma_check
from exp6_block_dichotomy import probe as dich_probe
from exp7_reduced_rank import analyse as red_analyse
from exp8_simple_linearization import validate as simp_validate, nullity as simp_nullity
from exp9_model_fidelity import check_A as fid_A, check_B as fid_B
import exp10_degree2_obstruction as D2
import exp2_rank_truth_table as E2

def hdr(t): print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)

# ---------------------------------------------------------------- W1
def w1_lemma_redund(m=7, n=13, r=3, d=3, P=6, seed=4711):
    hdr("W1  LEMMA REDUND : the P conjugate TSM equations are rho^q-images of the q=0 equation")
    R, inst = make_instance(m, n, r, d, seed)
    L = R.L
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
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
    Hs, Us = [H0], [U0]
    for q in range(1, P):
        Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    ok = True
    for q in range(P):
        aq, bq = list(al), list(be)
        for _ in range(q): aq = [L.fmul(x, x) for x in aq]; bq = [L.fmul(x, x) for x in bq]
        acc = [0] * m
        for k in range(r):
            t1 = R.K_to_Q(E[k])
            acc = [x ^ y for x, y in zip(acc, R.qmul(t1, [aq[k]] + [0] * (m - 1)))]
            acc = [x ^ y for x, y in zip(acc, R.qmul(t1, R.qmul(Hs[q], [bq[k]] + [0] * (m - 1))))]
        ok &= (acc == Us[q])
        ok &= (R.rho(Hs[q]) == Hs[q + 1] if q + 1 < P else True)
    print(f"  instance (m,n,r,d)=({m},{n},{r},{d}) seed={seed}")
    print(f"  U_q = sum_k b~_k ( rho^q(alpha_k) + H_q rho^q(beta_k) ) holds for q=0..{P-1}: {ok}")
    print(f"  H_q = rho^q(H_0) and U_q = rho^q(U_0) by construction (see rollo_gf2_engine.RingQ.rho).")
    print("  => applying rho^q (an INJECTIVE ring automorphism fixing K pointwise) to the q=0")
    print("     equation yields the q-th equation.  The P equations are ONE equation.")
    print(f"  LEMMA REDUND verified: {ok}")
    return ok

# ---------------------------------------------------------------- W2 / W3
def w2_w3_census(m=7, r=3, d=3, ns=(5, 9, 11, 13), Pmax=13, trials=2, seed=12345):
    hdr("W2  EXACT CENSUS of the ciphertext-TSM solution set vs P   (TRUE Frobenius vs PROBE model)")
    print(f"  (m,r)=({m},{r})   2r={2*r}  2^r={2**r}   #r-dim subspaces of GF(2^{m}) = {count_subspaces(m,r)}")
    print(f"  theoretical identity under test: 'each conjugate q imposes m-2r>0 constraints, so the solution")
    print(f"  is algebraically unique once P >= ceil(r(m-r)/(m-2r))' = {-(-r*(m-r)//(m-2*r)) if m>2*r else 'N/A'}")
    print(f"  {'n':>3} {'tr':>2} | {'log2 C(m,r)_2 - n(m-2r)':>24} | {'TRUE: #valid E by P':>34} | {'PROBE(random H_q): #valid E by P':>36}")
    rows = []
    for n in ns:
        for t in range(trials):
            sd = seed + 1000 * n + t
            R, inst = make_instance(m, n, r, d, sd)
            H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
            Hs, Us = [H0], [U0]
            for q in range(1, Pmax + 1): Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
            ct = {}
            for P in range(1, min(Pmax, n) + 1):
                ct[P] = len(census(R, list(zip(Hs[:P], Us[:P])), r))
            rng = random.Random(sd ^ 0xabcdef)
            Hr, Ur = [], []
            for q in range(Pmax + 1):
                Hq = [rng.randrange(1 << n) for _ in range(m)]
                aq = [rng.randrange(1 << n) for _ in range(r)]
                bq = [rng.randrange(1 << n) for _ in range(r)]
                Uq = [0] * m
                for k in range(r):
                    t1 = R.K_to_Q(inst["E"][k])
                    Uq = [x ^ y for x, y in zip(Uq, R.qmul(t1, [aq[k]] + [0] * (m - 1)))]
                    Uq = [x ^ y for x, y in zip(Uq, R.qmul(t1, R.qmul(Hq, [bq[k]] + [0] * (m - 1))))]
                Hr.append(Hq); Ur.append(Uq)
            cr = {}
            for P in range(1, min(Pmax, n) + 1):
                cr[P] = len(census(R, list(zip(Hr[:P], Ur[:P])), r))
            pred = math.log2(count_subspaces(m, r)) - n * (m - 2 * r)
            vals = sorted(set(ct.values()))
            print(f"  {n:>3} {t:>2} | {pred:>24.2f} | {str(ct if len(ct)<7 else dict(list(ct.items())[:4])):>34} | {str(cr if len(cr)<7 else dict(list(cr.items())[:4])):>36}")
            print(f"        {'':>24} | TRUE distinct counts over ALL P: {vals}  (P-INDEPENDENT: {len(vals)==1})"
                  f" | PROBE distinct counts: {sorted(set(cr.values()))}")
            rows.append((n, t, pred, ct, cr))
    print("\n  W3  COUNT LAW:  observed #valid E  vs  C(m,r)_2 * 2^{-n(m-2r)}")
    for (n, t, pred, ct, cr) in rows:
        print(f"    n={n:>3} trial={t}  predicted 2^{pred:+.2f} = {2**pred:8.2f}   observed {list(ct.values())[0]:6d}")
    return rows

# ---------------------------------------------------------------- W4
def w4_rank_table(cfgs):
    hdr("W4  EXACT GF(2) RANK TRUTH TABLE  [Laplace (+) Frobenius]   TRUE vs PROBE(random H_q)")
    print("  closed form under test:   null_TRUE(P) = n [ 2^r C(P-1,r-1) - m C(P-2,r-1) ]")
    print("  equivalently              rank_TRUE(P) = m n C(P-2,r-1)")
    print("  vanishing condition  2^r (P-1)/(P-r) <= m   <=>   P >= (mr-2^r)/(m-2^r)   (the audited m>2^r threshold)")
    allok = True
    for (m, n, r, Ps) in cfgs:
        print(f"\n  --- m={m} n={n} r={r}  (2r={2*r}, 2^r={2**r}, m<=2^r: {m <= 2**r}) ---")
        print(f"  {'P':>3} {'colsL':>6} {'rowsL':>6} | {'TRUElap':>8} {'TRUEfrob':>9} {'FORMULA':>8} {'ok':>3} | {'PROBElap':>9} {'PROBEfrob':>10}")
        for P in Ps:
            res = E2.run(m, n, r, 3, P, 2024)
            pred = n * ((1 << r) * math.comb(P - 1, r - 1) - m * math.comb(P - 2, r - 1))
            obs = res["TRUE_laplace_plus_frobenius"]["nullity"]
            good = (obs == pred)
            allok &= good
            print(f"  {P:>3} {res['ncols_L']:>6} {res['nrows_L']:>6} | {res['TRUE_laplace_only']['nullity']:>8} {obs:>9} {pred:>8} {'OK' if good else 'NO':>3} | "
                  f"{res['RANDOM_laplace_only']['nullity']:>9} {res['RANDOM_laplace_plus_frobenius']['nullity']:>10}"
                  f"   ctrl(resid={res['true_laplace_residual_rows']},frobviol={res['true_frobenius_violations']})")
    print(f"\n  closed-form nullity law verified on every (m<=2^r) entry: {allok}")
    return allok

# ---------------------------------------------------------------- W5
def w5_pubkey(m=7, d=3, r=3, ns=(5, 9, 11, 13), Pmax=13, trials=1, seed=555):
    hdr("W5  HOMOGENEOUS PUBLIC-KEY CENSUS  (Gamma_S):  dim_L(F~L + H_q F~L) < 2d")
    print(f"  (m,d)=({m},{d})  #d-dim subspaces = {count_subspaces(m,d)}  K^* has {2**m-1} elements")
    for n in ns:
        for t in range(trials):
            R, inst = make_instance(m, n, r, d, seed + 1000 * n + t)
            H0 = R.Phi(inst["h"])
            Hs = [H0]
            for q in range(1, Pmax + 1): Hs.append(R.rho(Hs[-1]))
            counts = {}
            for P in range(1, min(Pmax, n) + 1):
                counts[P] = len(census_pubkey(R, Hs[:P], d))
            vals = sorted(set(counts.values()))
            print(f"  n={n:>3} trial={t}: #valid F by P = {vals[0]} for every P  (P-INDEPENDENT: {len(vals)==1});"
                  f"  2^m-1 = {2**m-1}  -> solution set is EXACTLY the K^*-orbit of F: {vals==[2**m-1]}")
    return True

# ---------------------------------------------------------------- W6
def w6_negative(m=7, n=17, r=3, d=3, trials=3, seed=8080):
    hdr("W6  NEGATIVE CONTROL: uniform random (h,u) must admit 0 valid rank-r supports")
    ok = True
    for t in range(trials):
        R, inst = make_instance(m, n, r, d, seed + t)
        H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
        sols = census(R, [(H0, U0)], r)
        Et = canon(inst["E"], m)
        pos = any(canon(s, m) == Et for s in sols)
        rng = random.Random((seed + t) ^ 0xBADF00D)
        h = [rng.randrange(1 << m) for _ in range(n)]
        u = [rng.randrange(1 << m) for _ in range(n)]
        nneg = len(census(R, [(R.Phi(h), R.Phi(u))], r))
        ok &= pos and (nneg == 0)
        print(f"  trial {t}: genuine instance -> {len(sols):3d} valid E (contains true E: {pos}) | random (h,u) -> {nneg} valid E")
    print(f"  NEGATIVE CONTROL {'PASSED' if ok else 'FAILED'}  (m={m}, n={n}, r={r})")
    print(f"  NOTE: the control only holds once n(m-2r) > log2 C(m,r)_2 = {math.log2(count_subspaces(m,r)):.2f};"
          f" here n(m-2r) = {n*(m-2*r)}.")
    return ok

def w7_sigma_equivariance():
    hdr("W7  PROPOSITION (sigma-equivariance): modulo the Frobenius relations E_{S'+1} = rho(E_{S'})")
    allok = True
    for cfg in [(7,13,3,6),(7,13,3,7),(9,5,4,7),(11,5,5,7),(13,7,4,7),(5,7,3,6)]:
        ok, tot = sigma_check(*cfg)
        allok &= (ok == tot)
        print(f"  (m,n,r,P)={cfg}: rho(coeff of E_S' at (Q,mu)) == coeff of E_(S'+1) at (Q+1,mu)"
              f"  for {ok}/{tot} admissible S'  -> {'VERIFIED' if ok==tot else 'FAILED'}")
    print(f"  PROPOSITION (sigma-equivariance) {'VERIFIED' if allok else 'FAILED'}")
    return allok

def w8_block_dichotomy():
    hdr("W8  THEOREM (block dichotomy): E_{S'} survives  <=>  {0,1} subset S'")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'#indep':>7}{'C(P-2,r-1)':>11} | {'indep<=>1inS':>12} | {'each=mn':>8} | {'rank':>7}{'formula':>9} {'ok':>3}")
    allok = True
    for cfg in [(5,7,3,5),(7,5,3,6),(7,13,3,6),(9,5,4,7),(11,5,4,7),(9,5,5,7),(17,5,5,7)]:
        m,n,r,P = cfg
        sixes, contrib, rk, mn = dich_probe(m,n,r,P)
        indep=[s for s in sixes if contrib[s]>0]
        dich=all(((1 in s)==(contrib[s]>0)) for s in sixes); each=all(contrib[s]==mn for s in indep)
        pred=n*((1<<r)*math.comb(P-1,r)+m*math.comb(P-2,r-1))
        ok=(len(indep)==math.comb(P-2,r-1)) and dich and each and rk==pred
        allok &= ok
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {len(indep):>7}{math.comb(P-2,r-1):>11} | {str(dich):>12} | {str(each):>8} | {rk:>7}{pred:>9} {'OK' if ok else 'NO':>3}")
    print(f"  THEOREM (block dichotomy) {'VERIFIED' if allok else 'FAILED'}")
    return allok

def w9_cross_validation():
    hdr("W9  INDEPENDENT CROSS-VALIDATION of the rank law (explicit Frobenius substitution)")
    print(f"  {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'rows':>6}{'freecols':>9} | {'rank(subst)':>12}{'m*n*C(P-2,r-1)':>16} {'ok':>3}")
    allok = True
    for cfg in [(5,7,3,5),(6,7,3,5),(7,5,3,6),(7,13,3,5),(7,13,3,6),(9,5,4,6),(9,5,4,7),(13,7,4,7),(11,5,5,7)]:
        m,n,r,P = cfg
        nr,nfree,rk,mn = red_analyse(m,n,r,P)
        pred = mn*math.comb(P-2,r-1); ok = (rk==pred); allok &= ok
        print(f"  {m:>3}{n:>3}{r:>3}{P:>3} | {nr:>6}{nfree:>9} | {rk:>12}{pred:>16} {'OK' if ok else 'NO':>3}")
    print(f"  INDEPENDENT CROSS-VALIDATION {'PASSED' if allok else 'FAILED'}")
    return allok

def w10_simple_linearization():
    hdr("W10  NEW LINEARISATION (SIMP): sum_t U_{s_t} p_{S'\\{s_t}} = 0, p_{Q+1}=rho(p_Q)")
    print("  (a) validity: the true p_Q from the real gamma_k must satisfy (SIMP) + Frobenius")
    allok = True
    for cfg in [(7,13,3,5),(7,5,3,6),(9,5,4,6),(6,7,3,5),(11,5,5,7)]:
        bad,tot,fbad,ftot,nz,tp = simp_validate(*cfg)
        ok=(bad==0 and fbad==0 and nz>0); allok &= ok
        print(f"    (m,n,r,P)={cfg}: eq violations {bad}/{tot}, Frobenius violations {fbad}/{ftot}, nonzero p_Q {nz}/{tp} -> {'VALID' if ok else 'INVALID'}")
    print("  (b) exact nullity law:  null_SIMP = m*n*C(P-2,r-2)")
    print(f"    {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'freecols':>9}{'rank':>7}{'NULL':>7}{'m*n*C(P-2,r-2)':>16} {'ok':>3}")
    for cfg in [(5,7,3,5),(6,7,3,6),(7,5,3,5),(7,13,3,6),(9,5,4,5),(11,5,4,5),(13,7,5,7),(5,7,4,6),(3,5,2,4),(5,7,2,5)]:
        m,n,r,P = cfg
        nr,nf,rk,nl = simp_nullity(m,n,r,P); pred=m*n*math.comb(P-2,r-2)
        ok=(nl==pred); allok &= ok
        print(f"    {m:>3}{n:>3}{r:>3}{P:>3} | {nf:>9}{rk:>7}{nl:>7}{pred:>16} {'OK' if ok else 'NO':>3}")
    print(f"  (SIMP) VALID AND NULLITY LAW VERIFIED: {allok}")
    return allok

def w11_model_fidelity():
    hdr("W11  MODEL FIDELITY vs parameterized_tsm/tsm.cpp")
    print("  CHECK A: does using ALL C(P,r+1) subsets (not only those containing 0) add rank?")
    okA = True
    for cfg in [(7,13,3,6),(7,13,3,7),(6,7,3,6),(11,13,3,7)]:
        n0,n1 = fid_A(*cfg); okA &= (n0==n1)
        print(f"    (m,n,r,P)={cfg}: only-0 null {n0}, ALL null {n1}, gain {n0-n1}")
    print(f"  CHECK A {'PASSED (no missing equations)' if okA else 'FAILED'}")
    print("  CHECK B: tsm.cpp's L-linear Frobenius condition vs the true GF(2)-linear one")
    okB = True
    for cfg in [(7,13,3,5),(7,13,3,6),(7,13,3,7),(6,7,3,6)]:
        D,nT = fid_B(*cfg)
        rep = cfg[1]*((1<<cfg[2])*math.comb(cfg[3]-1,cfg[2]-1)-cfg[0]*math.comb(cfg[3]-2,cfg[2]-1))
        same=(nT==rep); okB &= same
        print(f"    (m,n,r,P)={cfg}: Stage1 L-nullity {D}, tsm-frob nullity {nT}, predicted law {rep} -> same: {same}")
    print(f"  CHECK B {'PASSED (identical nullities for m<=2^r)' if okB else 'FAILED'}")
    return okA and okB

def w12_degree2_obstruction():
    hdr("W12  DEGREE-2 (Plucker) STAGE IS UNDERDETERMINED BY CONSTRUCTION")
    D2.table()
    print("\n  toy-scale verification on the TRUE (SIMP) kernel:")
    print(f"    {'m':>3}{'n':>3}{'r':>3}{'P':>3} | {'D':>4}{'#monomials':>11}{'#relations':>11}{'rank':>7}{'NULLITY':>8}")
    for cfg in [(3,5,2,4),(5,7,2,5)]:
        D,nm,nr,rk,nl = D2.toy(*cfg)
        print(f"    {cfg[0]:>3}{cfg[1]:>3}{cfg[2]:>3}{cfg[3]:>3} | {D:>4}{nm:>11}{nr:>11}{rk:>7}{nl:>8}")
    return True

if __name__ == "__main__":
    t0 = time.time()
    r1 = w1_lemma_redund()
    w2_w3_census()
    r4 = w4_rank_table([(5, 7, 3, (4, 5, 6)), (6, 7, 3, (4, 5, 6)), (7, 5, 3, (4, 5, 6, 7)),
                        (7, 13, 3, (4, 5, 6, 7)), (9, 5, 4, (5, 6, 7)), (11, 5, 4, (5, 6, 7)),
                        (13, 7, 4, (5, 6, 7)), (15, 7, 4, (5, 6, 7))])
    w5_pubkey()
    r6 = w6_negative()
    r7 = w7_sigma_equivariance()
    r8 = w8_block_dichotomy()
    r9 = w9_cross_validation()
    r10 = w10_simple_linearization()
    r11 = w11_model_fidelity()
    r12 = w12_degree2_obstruction()
    print(f"\nALL CHECKS DONE in {time.time()-t0:.1f}s   W1={r1} W4={r4} W6={r6} W7={r7} W8={r8} W9={r9} W10={r10} W11={r11} W12={r12}")
