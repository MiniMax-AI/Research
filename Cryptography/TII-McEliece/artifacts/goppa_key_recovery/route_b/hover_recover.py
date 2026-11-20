"""Stage 5 of Route B (Frobenius-reduced HOVER): Goppa polynomial recovery and
O(n) support completion from Segre-Veronese fibers.

Setting
-------
D = published binary code (dim n - mt),  C = D^perp (dim mt, binary),
W = span_Fq(C) (dim mt).  For the true key (alpha, y = 1/g(alpha)) the Frobenius orbit
of the GRS code fills W:

    W = sum_{j=0}^{m-1} GRS_t(alpha, y)^{[2^j]} .                            (W1)

Consequences used here.

(W2)  For *any* PGL_2 image beta = sigma(alpha), sigma = [[a,b],[c,d]], there is a
      multiplier w, unique up to scale, with GRS_t(beta,w) = GRS_t(alpha,y) subset W:
      w_i = lam * y_i * (c alpha_i + d)^{t-1}.  It is found by linear algebra, since
      GRS_t(beta,w) subset W  <=>  beta^u (*) w in W for all u < t.

(W3)  Writing L(beta) = a - c*beta and G(beta) = L(beta)^t g((d beta - b)/L(beta)),
      one gets   1/w_i = mu * G(beta_i) / L(beta_i).   Hence the n x (t+3) linear system
          (l_0 + l_1 beta_i) = w_i * sum_k h_k beta_i^k
      has a one-dimensional solution space, (l,h) ~ (L, mu G).

(W4)  Every (b,d) with det = l_0 d + l_1 b != 0 then yields a *valid* key
          alpha'_i = (d beta_i + b)/(l_1 beta_i + l_0),
          ghat(x) * det^t = sum_k h_k (l_0 x + b)^k (l_1 x + d)^{t-k},
      normalised to monic.  Indeed ghat is g(s x) for the unknown scale s of (l,h), and
      alpha' = alpha/s; scaling the support and the polynomial together leaves the GRS
      code -- hence the published code -- unchanged.  So the PGL_2 ambiguity that the
      Segre geometry leaves open is *irrelevant*: any choice is accepted.

The fiber-extraction stage returns beta only for the coordinates that were not shortened away;
the missing ones are recovered from (W1): the multiplier vector y' and alpha'(*)y' both
lie in the 110-dimensional space W and are already known on 972 coordinates, which
determines them on the remaining ones.
"""
from __future__ import annotations
import random as rnd

from sage.all import GF, Matrix, PolynomialRing, vector


# ------------------------------------------------------------------ conic handling
def _polarB(cq, Fq):
    """Polar (bilinear) form of the quadratic form with coefficients cq (i<=j)."""
    B = Matrix(Fq, 3, 3)
    c00, c01, c02, c11, c12, c22 = list(cq)
    for (i, j, c) in ((0, 1, c01), (0, 2, c02), (1, 2, c12)):
        B[i, j] = c
        B[j, i] = c
    return B


def _qval(cq, x):
    c00, c01, c02, c11, c12, c22 = list(cq)
    return (c00 * x[0] ** 2 + c01 * x[0] * x[1] + c02 * x[0] * x[2]
            + c11 * x[1] ** 2 + c12 * x[1] * x[2] + c22 * x[2] ** 2)


def conic_params(lams, Fq, rng=None, verbose=True):
    """Parametrise the conic carrying the first Segre factor.

    lams: list of P^2 points (one per coordinate, None allowed).
    Returns list of homogeneous parameters (X:Z), or None on failure.  The parameter is
    determined up to PGL_2 = Aut(conic), which (W4) does not matter.
    """
    rng = rng or rnd.Random(20261)
    ok = [(j, l) for j, l in enumerate(lams) if l is not None]
    if len(ok) < 8:
        return None
    Ev = Matrix(Fq, len(ok), 6)
    for r, (_, l) in enumerate(ok):
        Ev[r] = [l[0] ** 2, l[0] * l[1], l[0] * l[2], l[1] ** 2, l[1] * l[2], l[2] ** 2]
    K = Ev.right_kernel()
    if K.dimension() != 1:
        if verbose:
            print(f"   conic fit: kernel dim {K.dimension()} (!= 1)")
        return None
    cq = K.basis()[0]
    B = _polarB(cq, Fq)
    pts = [l for _, l in ok]

    def prop(a, b):
        for i in range(3):
            if a[i]:
                return b == (b[i] / a[i]) * a
        return False

    # a conic point outside the list, to be used as the centre of the pencil
    P0 = None
    P1 = pts[0]
    for _ in range(400):
        v = vector(Fq, [Fq.random_element() for _ in range(3)])
        qv = _qval(cq, v)
        if qv == 0:
            continue
        s = (P1 * B * v) / qv
        cand = P1 + s * v
        if _qval(cq, cand) != 0:
            continue
        if not any(prop(cand, p) for p in pts):
            P0 = cand
            break
    if P0 is None:
        return None
    # complement of P0: two vectors independent of it
    U = Matrix(Fq, 1, 3, list(P0)).right_kernel().basis_matrix()
    if U.nrows() != 2:
        return None
    T = Matrix(Fq, 3, 3, [list(P0), list(U.row(0)), list(U.row(1))])
    # rows of T are the basis vectors, so coordinates need the inverse TRANSPOSE
    try:
        Tinv = T.transpose().inverse()
    except ZeroDivisionError:
        return None
    out = [None] * len(lams)
    for j, l in ok:
        c = Tinv * l
        if c[1] == 0 and c[2] == 0:
            return None                      # would mean l == P0, excluded above
        out[j] = (c[1], c[2])
    # sanity: distinct coordinates must give distinct parameters
    seen = set()
    for p in out:
        if p is None:
            continue
        X, Z = p
        nrm = (Fq.one(), X / Z) if Z else (Fq.zero(), Fq.one())
        if nrm in seen:
            if verbose:
                print("   duplicate conic parameters")
            return None
        seen.add(nrm)
    return out


def first_factors(FF, fib, Fq):
    """First-factor P^2 point of every fibre (= nu_2(beta_j) in the fibre frame)."""
    out = []
    for F in fib:
        if F is None:
            out.append(None)
            continue
        lam = None
        for r in range(F.nrows()):
            row = vector(Fq, F.row(r))
            if row.is_zero():
                continue
            lam = FF.lambdas(row)
            if lam is not None:
                break
        out.append(lam)
    return out


# --------------------------------------------- (W2'): the squared multiplier by a flag
def multiplier_square(Hbin, params, keep, Fq, t, nsub=240, rng=None, verbose=True):
    """Isolate  y^2 / L^{2t-1}  (up to scale) on a subset of coordinates.

    Y := { w : w (*) beta^u in W  for all u < t }  contains  (y^2/L^{2t-1}) (*) F_{t+1}(beta),
    because W contains (y^2/L^{2t-1}) (*) F_{2t}(beta)  (the "squared Goppa" checks of the
    binary code).  Multiplication by beta therefore induces on Y the flag

        Y  >  (y^2/L^{2t-1}) (*) F_t  >  ...  >  (y^2/L^{2t-1}) (*) F_1  >  span{y^2/L^{2t-1}},

    i.e.  V_{k+1} = { v in V_k : v (*) beta in V_k }  isolates the generator after t steps.
    (Measured: dim Y = t+2, the flag drops 13 -> 11 -> 10 -> ... -> 1 for t = 11.)
    """
    rng = rng or rnd.Random(97531)
    fin = [i for i, p in enumerate(params) if p is not None and p[1] != 0]
    if len(fin) < nsub:
        nsub = len(fin)
    S = sorted(rng.sample(fin, nsub))
    s = len(S)
    orig = [keep[i] for i in S]
    beta = [params[i][0] / params[i][1] for i in S]
    Hs = Hbin.change_ring(Fq).matrix_from_columns(orig)
    A = Hs.right_kernel().basis_matrix()          # (s - mt) x s  = (W|_S)^perp
    if A.nrows() == 0:
        return None
    rows = []
    for u in range(t):
        du = [b ** u for b in beta]
        for r in range(A.nrows()):
            Ar = A.row(r)
            rows.append([Ar[k] * du[k] for k in range(s)])
    V = Matrix(Fq, rows).right_kernel().basis_matrix()
    if verbose:
        print(f"[recover] multiplier space dim {V.nrows()}", flush=True)
    if V.nrows() < t:
        return None
    for _ in range(t + 2):
        if V.nrows() <= 1:
            break
        B = Matrix(Fq, [[V.row(r)[i] * beta[i] for i in range(s)] for r in range(V.nrows())])
        K = B.transpose().augment(-V.transpose()).right_kernel().basis_matrix()
        if K.nrows() == 0:
            return None
        d = V.nrows()
        C = K.matrix_from_columns(range(d))
        N = C * V
        if N.nrows() == d:                        # flag stalled
            return None
        V = N
    if V.nrows() != 1:
        return None
    return S, beta, V.row(0)


# ------------------------------------- (W3'): from the squared multiplier to (L, G)
def goppa_from_square(beta, w2, Fq, t, verbose=True):
    """Solve  L(beta)/w2 = G(beta)^2  with L linear and deg G = t.

    Indeed 1/(y^2/L^{2t-1}) = L^{2t-1} g(alpha)^2 = L^{2t-1} G(beta)^2 / L^{2t}
                            = G(beta)^2 / L(beta) ,
    so the system  (a + b beta_i)/w2_i = sum_r R_r beta_i^r  (deg R <= 2t) is linear in
    (a, b, R) and has a one-dimensional solution space;  R = G^2 has even exponents only,
    which is both a proof of correctness and the way to take the square root.
    """
    s = len(beta)
    inv = [Fq.one() / w2[i] for i in range(s)]
    rows = [[inv[i], beta[i] * inv[i]] + [-beta[i] ** r for r in range(2 * t + 1)]
            for i in range(s)]
    K = Matrix(Fq, rows).right_kernel()
    if K.dimension() < 1:
        if verbose:
            print("   (L,G^2) system has no solution")
        return None
    if K.dimension() > 1:
        # the true solution has R = G^2 with even exponents only: intersect with that
        # subspace (a linear condition) to isolate it.
        cond = Matrix(Fq, t, 2 * t + 3)
        for k in range(t):
            cond[k, 2 + 2 * k + 1] = Fq.one()
        # v = c * K (rows of K span the kernel); the condition cond*v = 0 is linear in c
        Kc = (cond * K.transpose()).right_kernel().basis_matrix() * K
        if Kc.nrows() >= 1:
            K = Kc
    if K.dimension() != 1:
        if verbose:
            print(f"   (L,G^2) space dim {K.dimension()} (!= 1)")
        return None
    v = list(K.basis()[0])
    l0, l1, R = v[0], v[1], v[2:]
    if not all(R[r] == 0 for r in range(1, 2 * t + 1, 2)):
        if verbose:
            print("   G^2 has odd-degree terms: not a square")
        return None
    G = [R[2 * k].sqrt() for k in range(t + 1)]
    if all(c == 0 for c in G) or G[t] == 0:
        return None
    return (l0, l1), G


# ------------------------------------------------- (W4): candidate key from (l,h)
def key_from_lh(l, h, params, Fq, t, verbose=True):
    """All (b,d) candidates: yields (ghat, alpha_list, b, d) for each valid choice."""
    Rx = PolynomialRing(Fq, 'x')
    x = Rx.gen()
    l0, l1 = l
    cands = [(Fq.zero(), Fq.one()), (Fq.one(), Fq.zero())]
    # a couple more, in case the first two give a degenerate degree
    for extra in (Fq.from_integer(2), Fq.from_integer(3)):
        cands.append((extra, Fq.one()))
        cands.append((Fq.one(), extra))
    out = []
    for (b, d) in cands:
        det = l0 * d + l1 * b
        if det == 0:
            continue
        ghat = Rx.zero()
        for k in range(t + 1):
            if h[k]:
                ghat += h[k] * (l0 * x + b) ** k * (l1 * x + d) ** (t - k)
        if ghat.degree() != t:
            continue
        ghat = ghat * (Fq.one() / ghat[t])
        alpha = []
        bad = False
        for p in params:
            if p is None:
                alpha.append(None)
                continue
            X, Z = p
            den = l1 * X + l0 * Z
            if den == 0:
                bad = True
                break
            alpha.append((d * X + b * Z) / den)
        if bad:
            continue
        out.append((ghat, alpha, b, d))
        if len(out) >= 3:
            break
    return out


# ------------------------------- complete the coordinates removed by the shortening
def complete_shortened(Hbin, ghat, alpha, keep, Fq, verbose=True):
    """(W1): y' and alpha'(*)y' lie in W and are known on `keep`; extend to all n.

    alpha is indexed by position in `keep` (None where fiber extraction gave no label).  The
    restriction of W to the labelled coordinates is computed once; its pivot columns give
    an invertible mt x mt block, which determines y' and alpha'(*)y' inside W -- and hence
    their values at the coordinates the shortening removed.
    """
    n = Hbin.ncols()
    W = Hbin.change_ring(Fq)
    mt = W.nrows()
    known = [p for p in range(len(keep)) if alpha[p] is not None]
    if len(known) < mt:
        if verbose:
            print("   completion: too few labelled coordinates")
        return None
    kidx = [keep[p] for p in known]
    Wk = W.matrix_from_columns(kidx)
    piv = Wk.echelon_form().pivots()
    if len(piv) < mt:
        if verbose:
            print(f"   completion: rank of W on labelled coordinates is {len(piv)} < {mt}")
        return None
    Msub = Wk.matrix_from_columns(list(piv)[:mt])
    try:
        Minv = Msub.inverse()
    except ZeroDivisionError:
        return None
    vecs = {}
    for name in ('y', 'a'):
        rhs = vector(Fq, [(Fq.one() / ghat(alpha[known[c]]))
                          * (Fq.one() if name == 'y' else alpha[known[c]])
                          for c in list(piv)[:mt]])
        vecs[name] = (rhs * Minv) * W                    # c = rhs * Msub^-1, element c*W
    yfull, afull = vecs['y'], vecs['a']
    # self-check on the labelled coordinates that were not used
    chk = list(piv)[:mt]
    rest = [c for c in range(len(known)) if c not in set(chk)]
    for c in rest[::37]:
        j = kidx[c]
        if yfull[j] != Fq.one() / ghat(alpha[known[c]]):
            if verbose:
                print("   completion: inconsistent on a labelled coordinate")
            return None
    out = [None] * n
    for p, j in enumerate(keep):
        if alpha[p] is not None:
            out[j] = alpha[p]
    for j in range(n):
        if out[j] is None:
            if yfull[j] == 0:
                return None
            out[j] = afull[j] / yfull[j]
    return out


# ------------------------------------------------------------------- verification
def build_dual_Fq(ghat, alpha, Fq, m, t):
    """Fq-row space of the binary dual of Gamma(alpha, ghat)  (== H_rec row space)."""
    n = len(alpha)
    ginv = [Fq.one() / ghat(a) for a in alpha]
    rows = []
    for u in range(t):
        base = ginv if u == 0 else [ginv[i] * alpha[i] ** u for i in range(n)]
        bits = [[] for _ in range(m)]
        for v in base:
            iv = v._integer_representation()
            for j in range(m):
                bits[j].append((iv >> j) & 1)
        for j in range(m):
            rows.append([Fq.from_integer(b) for b in bits[j]])
    return Matrix(Fq, rows)


def verify_key(Hbin, ghat, alpha, Fq, m, t, verbose=True):
    """Exact acceptance test: span_Fq(dual of Gamma) == span_Fq(published dual)."""
    if ghat.degree() != t or not ghat.is_monic() or not ghat.is_irreducible():
        return False
    if len(set(alpha)) != len(alpha):
        return False
    if any(ghat(a) == 0 for a in alpha):
        return False
    Hc = build_dual_Fq(ghat, alpha, Fq, m, t)
    W = Hbin.change_ring(Fq)
    if Hc.rank() != W.nrows():
        return False
    return Matrix(Fq, list(W.rows()) + list(Hc.rows())).rank() == W.nrows()


# ------------------------------------------------------------------------ driver
def recover_from_fibers(FF, fib, Hbin, keep, Fq, m, t, n, rng=None, verbose=True):
    """Full Goppa polynomial & support recovery from extracted fibers."""
    import time
    t0 = time.time()
    lams = first_factors(FF, fib, Fq)
    if verbose:
        print(f"[recover] first factors: {sum(1 for l in lams if l is not None)}/{len(lams)}"
              f"   {time.time()-t0:.1f}s", flush=True)
    params = conic_params(lams, Fq, rng=rng, verbose=verbose)
    if params is None:
        return None
    if verbose:
        print(f"[recover] conic parametrised   {time.time()-t0:.1f}s", flush=True)
    ms = None
    for k in range(3):
        ms = multiplier_square(Hbin, params, keep, Fq, t,
                               rng=rnd.Random(97531 + 7919 * k), verbose=verbose)
        if ms is not None:
            break
    if ms is None:
        return None
    S, beta, w2 = ms
    if verbose:
        print(f"[recover] squared multiplier on {len(S)} coordinates   "
              f"{time.time()-t0:.1f}s", flush=True)
    fit = goppa_from_square(beta, w2, Fq, t, verbose=verbose)
    if fit is None:
        return None
    l, h = fit
    cands = key_from_lh(l, h, params, Fq, t, verbose=verbose)
    if not cands:
        return None
    for (ghat, alpha972, b, d) in cands:
        alpha = complete_shortened(Hbin, ghat, alpha972, keep, Fq, verbose=verbose)
        if alpha is None:
            continue
        if not verify_key(Hbin, ghat, alpha, Fq, m, t, verbose=verbose):
            continue
        if verbose:
            print(f"[recover] key verified: g monic irreducible of degree {ghat.degree()}, "
                  f"support {len(set(alpha))} distinct   {time.time()-t0:.1f}s", flush=True)
        coeffs = [int(ghat[i]._integer_representation()) for i in range(t + 1)]
        sup = [int(a._integer_representation()) for a in alpha]
        return coeffs, sup
    return None
