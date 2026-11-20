#!/usr/bin/env python3
"""Verify the *reformulation theorem* that underlies the TSM attack.

Setting.  R = GF(2^m)[X]/(f), deg f = n, gcd(m,n)=1, so R is the field GF(2^(mn))
with subfields K = GF(2^m) (the constants) and Lsub = GF(2^n) = GF(2)[X]/(f).
The ideal-LRPC (ROLLO-I / LOCKER-I) encapsulation publishes h and
    u = e1 + e2 h,   supp(e1) = supp(e2) = E,  dim_{GF(2)} E = r.

THEOREM (reformulation).  Put W = E * GF(2^n) = span_{GF(2^n)}(E) inside R.  Then
  (T1)  dim_{GF(2^n)} W = r                                   (E is GF(2)-independent)
  (T2)  W  cap  GF(2^m)  =  E                                 (exactly, as GF(2)-spaces)
  (T3)  u in W + h W  =  span_{GF(2^n)}(E  u  hE), whose dimension is 2r generically
  (T4)  conversely, u in W + hW with W of the form E*GF(2^n) recovers E = W cap GF(2^m)
  (T5)  every Frobenius conjugate equation  U_q = sum_k E~_k sigma^q(c_k)  is a
        consequence of the single equation q = 0 together with E~_k in GF(2^m)
        (this is what makes the P-fold "free" system consistent, and what the
         Frobenius pinning stage re-imposes after linearisation).

(T1)+(T2) are exactly where gcd(m,n)=1 is used: they say the tensor map
GF(2^n) (x) GF(2^m) -> R is an isomorphism of GF(2)-spaces, so no GF(2^n)-combination
of GF(2^m)-elements can collapse.

usage: verify_reformulation.sage --m 71 --n 83 --r 5 [--seed 3] [--trials 3]
"""
import argparse, itertools, json
from sage.all import GF, PolynomialRing, matrix, set_random_seed, gcd


def build(m, n, seed):
    set_random_seed(seed)
    P = PolynomialRing(GF(2), "x")
    g = P.irreducible_element(m)
    f = P.irreducible_element(n)
    K = GF(2 ** m, "t", modulus=g)     # GF(2^m)  (the "E" side)
    L = GF(2 ** n, "a", modulus=f)     # GF(2^n)  (the "coefficient" side)
    Q = PolynomialRing(L, "T").quotient(K.modulus().change_ring(L), "z")  # GF(2^n)[T]/(g)
    return P, g, f, K, L, Q


def phi(coeffs, K, L, Q, m, n):
    """R = GF(2^m)[X]/(f)  ->  Q = GF(2^n)[T]/(g):  Phi(a) = sum_b (sum_j (a_j)_b X^j) T^b."""
    bits = [[(int(c.to_integer()) >> b) & 1 for b in range(m)] for c in coeffs]
    return Q([L([bits[j][b] for j in range(n)]) for b in range(m)])


def qrank(vecs, L, Q, MDEG):
    """rank over GF(2^n) of a list of Q-elements (rows = T-coefficients)."""
    if not vecs:
        return 0
    rows = []
    for v in vecs:
        cs = list(v.lift().list())
        cs += [L(0)] * (MDEG - len(cs))
        rows.append([L(c) for c in cs])
    return matrix(L, rows).rank()


def subfield_intersection(basis, K, L, Q, m, n):
    """W cap GF(2^m) as a GF(2)-space, given a GF(2^n)-basis of W (list of Q elts)."""
    d = len(basis)
    # x = sum lambda_k basis_k in GF(2^m)  <=>  x^{2^m} = x  (K is the fixed field of sigma^m)
    # sigma^m(lambda) = lambda^{2^m}.  Unknowns: lambda in GF(2^n)^d  (d*n GF(2)-bits).
    rows = []
    for i in range(d):
        for j in range(m):                      # T-coefficient index of the equation
            pass
    # Solve over GF(2): sum_k (lambda_k^{2^m} - lambda_k) * basis_k = 0 in Q.
    # Represent lambda_k by its n GF(2)-bits; lambda_k^{2^m} is GF(2)-linear in them.
    # Build the linear map  lambda -> sum_k lambda_k^{2^m} basis_k  and  lambda -> sum lambda_k basis_k
    B = []
    for v in basis:
        cs = list(v.lift().list())
        cs += [L(0)] * (m - len(cs))
        B.append([L(c) for c in cs])                          # d x m over GF(2^n)
    # columns of the unknown space: bit (k, t) = t-th bit of lambda_k in the L-basis
    colsA, colsB = [], []
    for k in range(d):
        for t in range(n):
            e = L.gen() ** t
            em = e ** (2 ** m)                                 # Frobenius 2^m on GF(2^n)
            colsA.append([B[k][j] * e for j in range(m)])      # lambda_k basis_k
            colsB.append([B[k][j] * em for j in range(m)])     # lambda_k^{2^m} basis_k
    rows = []
    for j in range(m):
        for bit in range(n):
            row = []
            for c in range(len(colsA)):
                a = int((colsA[c][j] + colsB[c][j]).to_integer())
                row.append((a >> bit) & 1)
            rows.append(row)
    Mx = matrix(GF(2), rows)
    ker = Mx.right_kernel().basis()
    # map kernel vectors back to elements of Q and keep those in GF(2^m)
    out = []
    for v in ker:
        x = Q(0)
        idx = 0
        for k in range(d):
            lam = L([int(v[idx + t]) for t in range(n)])
            idx += n
            x += Q(lam) * basis[k]
        out.append(x)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=71)
    ap.add_argument("--n", type=int, default=83)
    ap.add_argument("--r", type=int, default=5)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--pairs", type=int, default=9)
    args = ap.parse_args()
    m, n, r = args.m, args.n, args.r
    P, g, f, K, L, Q = build(m, n, args.seed)
    report = {"m": m, "n": n, "r": r, "gcd": int(gcd(m, n)), "trials": []}

    def to_Q_from_R(coeffs):
        return phi(coeffs, K, L, Q, m, n)

    def rand_bin():
        return [GF(2).random_element() for _ in range(n)]

    def lift(elem_R):                      # Q-element -> list of n K-coefficients
        pass

    for trial in range(args.trials):
        while True:
            E = [K.random_element() for _ in range(r)]
            if matrix(GF(2), r, m, lambda i, j: (E[i].to_integer() >> j) & 1).rank() == r:
                break
        a = [rand_bin() for _ in range(r)]
        b = [rand_bin() for _ in range(r)]

        def build_err(vec):
            coeffs = [K(0)] * n
            for k in range(r):
                for j in range(n):
                    if vec[k][j]:
                        coeffs[j] += E[k]
            return coeffs

        e1c, e2c = build_err(a), build_err(b)
        hc = [K.random_element() for _ in range(n)]
        # u = e1 + e2 * h  in R = K[X]/(f)
        fexp = [int(e) for e in f.exponents()]        # x^n + sum of lower terms
        low = [e for e in fexp if e < n]

        def mulmod(p, q):
            acc = [K(0)] * (2 * n - 1)
            for i in range(n):
                if not p[i]:
                    continue
                for j in range(n):
                    if q[j]:
                        acc[i + j] += p[i] * q[j]
            for d in range(2 * n - 2, n - 1, -1):      # x^d = sum_low x^{d-n+low}
                c = acc[d]
                if c:
                    for e in low:
                        acc[d - n + e] += c
                    acc[d] = K(0)
            return acc[:n]
        e2h = mulmod(e2c, hc)
        uc = [e1c[j] + e2h[j] for j in range(n)]

        Eq = [to_Q_from_R([E[k]] + [K(0)] * (n - 1)) for k in range(r)]   # E~_k
        hq = to_Q_from_R(hc)
        uq = to_Q_from_R(uc)
        ak = [to_Q_from_R([K(a[k][j]) for j in range(n)]) for k in range(r)]
        bk = [to_Q_from_R([K(b[k][j]) for j in range(n)]) for k in range(r)]
        e1q = to_Q_from_R(e1c)
        e2q = to_Q_from_R(e2c)

        # --- basic identities -------------------------------------------------
        T0 = all(qrank([Eq[k]], L, Q, m) == 1 for k in range(r))
        model_ok = (uq == e1q + e2q * hq) and \
                   all(e1q == sum(Q(ak[k]) * Eq[k] for k in range(r)) for _ in [0]) and \
                   all(e2q == sum(Q(bk[k]) * Eq[k] for k in range(r)) for _ in [0])

        # --- (T1) dim_{GF(2^n)} W = r ----------------------------------------
        dimW = qrank(Eq, L, Q, m)
        # --- (T3) u in W + hW, dim(W + hW) -----------------------------------
        span = Eq + [hq * Eq[k] for k in range(r)]
        dimV = qrank(span, L, Q, m)
        u_in_V = qrank(span + [uq], L, Q, m) == dimV
        # --- (T2) W cap GF(2^m) = E ------------------------------------------
        inter = subfield_intersection(Eq, K, L, Q, m, n)
        # express each intersection element in the T-basis and read off its bits
        inter_bits = []
        for x in inter:
            cs = list(x.lift().list()) + [L(0)] * (m - len(x.lift().list()))
            inter_bits.append(tuple(int(c.to_integer()) for c in cs))
        # E as bit vectors in the same (T-basis) coordinates
        E_bits = []
        for k in range(r):
            v = Eq[k].lift().list() + [L(0)] * (m - len(Eq[k].lift().list()))
            E_bits.append(tuple(int(c.to_integer()) for c in v))
        ME = matrix(GF(2), [list(v) for v in E_bits])
        MI = matrix(GF(2), [list(v) for v in inter_bits]) if inter_bits else matrix(GF(2), 0, n * m)
        same = (MI.row_space() == ME.row_space())

        # --- (T5) all P Frobenius equations follow from q=0 -------------------
        P_ = args.pairs
        sigma = lambda z: Q([c ** 2 for c in z.lift().list()])
        conj_ok = True
        Hq, Uq = hq, uq
        for q in range(P_):
            lhs = sum(Eq[k] * (Q(ak[k]) ** (2 ** q) + Hq * (Q(bk[k]) ** (2 ** q))
                               ) for k in range(r))
            if lhs != Uq:
                conj_ok = False
            Uq, Hq = sigma(Uq), sigma(Hq)
        # and the "one equation" form: sigma^q(u) is forced
        one_eq = all(sigma(e1q) == sum(sigma(Q(ak[k])) * Eq[k] for k in range(r)) for _ in [0])

        # --- (T4) a random r-dim GF(2^n)-space almost never works -------------
        good = 0
        for _ in range(20):
            Wp = [Q([L.random_element() for _ in range(m)]) for _ in range(r)]
            if qrank(Wp, L, Q, m) != r:
                continue
            sp = Wp + [hq * w for w in Wp]
            if qrank(sp + [uq], L, Q, m) == qrank(sp, L, Q, m):
                good += 1
        report["trials"].append({
            "trial": trial,
            "model_ok": bool(model_ok),
            "dim_W_over_GF2n": int(dimW),                 # T1: expect r
            "dim_W_plus_hW": int(dimV),                   # T3: expect 2r
            "u_in_W_plus_hW": bool(u_in_V),               # T3
            "dim_W_cap_GF2m": int(MI.nrows()),            # T2: expect r
            "W_cap_GF2m_equals_E": bool(same),            # T2
            "frobenius_conjugates_hold": bool(conj_ok),   # T5
            "random_r_space_hits": good,                  # T4: expect 0/20
        })
        print(json.dumps(report["trials"][-1], separators=(",", ":")), flush=True)

    ok = all(t["dim_W_over_GF2n"] == r and t["dim_W_plus_hW"] == 2 * r
             and t["u_in_W_plus_hW"] and t["dim_W_cap_GF2m"] == r
             and t["W_cap_GF2m_equals_E"] and t["frobenius_conjugates_hold"]
             and t["random_r_space_hits"] == 0 and t["model_ok"]
             for t in report["trials"])
    report["all_pass"] = bool(ok)
    print(json.dumps({"ALL_PASS": bool(ok)}, separators=(",", ":")))


if __name__ == "__main__":
    main()
