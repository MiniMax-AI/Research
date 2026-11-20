"""Higher-Order Vanishing (HOV) machinery for binary Goppa key recovery.

Implements Definition 4/6 of Hemmert-Wiemers (ePrint 2025/1661) specialised to
square-free forms (Lemma 2) and to the derivative orders that survive in
characteristic 2 (Remark 1 / Remark 4 / Lemma 5), plus the parameter rules (3) of
ePrint 2026/1778.

Only p = 3 is supported (the degree used for (m,t) = (10,11)).
"""
from __future__ import annotations
from itertools import combinations
from sage.all import GF, Matrix, vector, binomial


def hov_params(m, t, p):
    """e, s*, kappa, a, b  -- equation (3) of ePrint 2026/1778."""
    # e = 1 + floor(log2((t-1)/(p-1)))
    num, den = t - 1, p - 1
    e = 1
    v = den
    while 2 * v <= num:
        v *= 2
        e += 1
    sstar = (e + 2) * t - (p - 1) * 2 ** e - p
    kappa = m * t - sstar
    a = m - e - 2
    b = (m - 1) * t - sstar
    return e, sstar, kappa, a, b


def shorten(H, cols):
    """Shorten the code spanned by the rows of H on the coordinates `cols`."""
    A = H.matrix_from_columns(list(cols))
    K = A.left_kernel().basis_matrix()
    Hs = K * H
    keep = [c for c in range(H.ncols()) if c not in set(cols)]
    return Hs.matrix_from_columns(keep)


def hov_points(Hs):
    """Definition 6: RREF of Hs, pivot columns removed.  Returns (P, pivots)."""
    R = Hs.echelon_form()
    piv = []
    for i in range(R.nrows()):
        nz = R.nonzero_positions_in_row(i)
        piv.append(nz[0] if nz else None)
    pivset = set(p for p in piv if p is not None)
    keep = [c for c in range(R.ncols()) if c not in pivset]
    return R.matrix_from_columns(keep), piv


def pair_matrix(P):
    """W[a,(i,j)] = P[i,a]*P[j,a] for i<j.  Rows = points, cols = pairs."""
    F = P.base_ring()
    kappa, npts = P.nrows(), P.ncols()
    pairs = list(combinations(range(kappa), 2))
    pidx = {pr: k for k, pr in enumerate(pairs)}
    import numpy as np
    Pnp = np.array(P.rows(), dtype=np.uint8)
    Wnp = np.empty((npts, len(pairs)), dtype=np.uint8)
    for k, (i, j) in enumerate(pairs):
        Wnp[:, k] = Pnp[i] & Pnp[j]
    W = Matrix(F, npts, len(pairs), Wnp.flatten().tolist())
    return W, pairs, pidx


def hov_operator_p3(P):
    """Binary HOV operator for p=3: rows (point, r), columns = square-free cubics.

    T[(a,r),{i,j,l}] = [r in {i,j,l}] * P[i',a]*P[j',a]  with {i',j'} = {i,j,l}\\{r}.

    Construction note.  The naive build materialises kappa blocks by selecting all
    C(kappa,3) columns of W and then stacks them one at a time; that is 4.8e6 column
    copies plus a quadratic re-copy of the whole matrix and costs ~50 s at (m, t) = (10, 11).
    Here the operator is allocated once (m4ri allocation is free) and the
    kappa*C(kappa-1,2) = 194472 non-zero columns are written as ~67k maximal runs of
    consecutive columns with `set_block`, which costs ~6 s.  The matrix produced is
    bit-identical.
    """
    F = P.base_ring()
    kappa, npts = P.nrows(), P.ncols()
    W, pairs, pidx = pair_matrix(P)
    triples = list(combinations(range(kappa), 3))
    ntri = len(triples)
    zc = len(pairs)
    Waug = W.augment(Matrix(F, npts, 1))          # a zero column to scatter in
    colmat = [Waug.matrix_from_columns([k]) for k in range(zc + 1)]

    T = Matrix(F, kappa * npts, ntri)
    for r in range(kappa):
        pos, pr = [], []
        for c, (i, j, l) in enumerate(triples):
            if r == i:
                pos.append(c)
                pr.append(pidx[(j, l)])
            elif r == j:
                pos.append(c)
                pr.append(pidx[(i, l)])
            elif r == l:
                pos.append(c)
                pr.append(pidx[(i, j)])
        ro = r * npts
        a = 0
        while a < len(pos):
            b = a
            while b + 1 < len(pos) and pos[b + 1] == pos[b] + 1:
                b += 1
            if b == a:
                T.set_block(ro, pos[a], colmat[pr[a]])
            else:
                T.set_block(ro, pos[a], Waug.matrix_from_columns(pr[a:b + 1]))
            a = b + 1
    return T, triples


def hov_kernel(Hs, p=3):
    """Binary HOV kernel of the shortened dual Hs.  Returns (basis rows, triples)."""
    P, piv = hov_points(Hs)
    if p != 3:
        raise NotImplementedError("only p=3")
    T, triples = hov_operator_p3(P)
    K = T.right_kernel().basis_matrix()
    return K, triples, P


def kernel_to_forms(K, triples, kappa, F):
    """Rows of K are coefficient vectors over the triple list -> cubic forms.

    Returns a list of dicts {(i,j,l): coeff} over the field F.
    """
    out = []
    for row in K:
        d = {}
        for c, tr in enumerate(triples):
            if row[c]:
                d[tr] = F(row[c])
        out.append(d)
    return out
