#!/usr/bin/env python3
"""Key recovery for a published binary Goppa code by shortened interpolation.

Pipeline (deterministic; no secret input, no unbounded search):

  1.  ``H`` is the 40 x 253 parity check of the published [253, 213] code and
      ``A`` a set of ``n'`` positions.  ``ker(H|_A)`` is the code *shortened*
      to ``A``; it is again a Goppa code with the same ``g`` and minimum
      weight ``2t+1``, so every minimum-weight codeword of it is a
      minimum-weight codeword of the published code.

  2.  For a minimum-weight codeword with support ``S`` the locator
      ``Q_S = prod_{j in S} (X - z_j)`` has degree ``2t+1`` and, in
      characteristic two, ``Q'_S`` is the square of a degree-``t``
      polynomial.  A residue computation on ``Q'/Q`` gives ``g^2 | Q'_S`` for
      *every* codeword, so at weight ``2t+1`` both sides are monic of degree
      ``2t`` and ``Q'_S = g^2`` exactly.  Hence the even elementary symmetric
      functions of the support are constants:  ``e_2(z_S) = g_{t-1}^2``, ...

  3.  ``e_2`` is a quadratic function of the codeword, so its polar form is a
      symmetric bilinear form on the ``k'-1``-dimensional affine chart of odd
      message vectors.  Interpolating the whole space of quadratic relations
      held by the sampled minimum words gives a quotient of dimension ``m+1``
      or ``2m+1``; inside it the Frobenius bit slices of the ``e_2`` polar
      form are exactly the unusually low-rank members.  One further relation
      is known *before* any label: ``sum_{i<j} c_i c_j = C(2t+1,2) = 1``, whose
      polar form ``G0 (J + I) G0^T`` is public and is the unique high-rank
      member of the locus.  Dividing it out leaves a low-rank locus of
      dimension exactly ``m`` with no exceptional direction, which is read off
      an exhaustive, parameter-free census of the ``2^{r-1}``-element quotient.

  4.  Lifting the subspace recovers the ``z^2``-diagonal, hence the labels on
      ``A`` and ``g``.  Both halves of the lift are *linear* solves, not
      searches:

        * the non-pivot weights ``z`` come from one kernel computation.  The
          rank-``m`` lifted form makes any ``k`` of its rows dependent, and
          ``sum_r alpha_r L0[r,j] = 0`` reads
          ``sum_r alpha_r B[r,j] + sum_h (alpha P[R,h]) z_h P[j,h] = 0``.
          Linearising ``y = (alpha P[R,:]) . z`` turns "find ``alpha``" into
          the kernel of a single ``|J| x (k + 41)`` matrix, which yields the
          whole ``k-m`` dimensional dependency space at once, so one block
          determines all 41 bits with probability ``1 - 41 2^{-(k-m)}``.

        * the diagonal comes from the same rank-``m`` condition read on the
          off-diagonal part alone: the left kernel of ``off[R, J]`` for a
          9-set ``R`` *is* the row dependency, and then
          ``L0[j,j] = sum_{r != j} alpha_r off[r,j]``.  This replaces a search
          over ``C(n, m)`` index sets times ``2^m`` diagonal patterns.

      The ``m`` bit slices share one coefficient matrix, so their weight
      vectors come from a *single* elimination with ``m`` right-hand sides
      rather than ``m`` eliminations.

  5.  The abstract algebra carried by the lifted forms is identified with the
      public key's field by one characteristic polynomial: an element whose
      characteristic polynomial on the ``m``-dimensional regular
      representation is irreducible of degree ``m`` generates the algebra, and
      sending it to a root of the same polynomial in the public key's field is a
      field isomorphism.  This replaces the enumeration of all ``2^m - 1``
      nonzero elements of the algebra.

  6.  The remaining labels follow from the published parity check: the key
      fixes one particular basis ``T`` of the published row space through
      ``H_cand = T H`` on ``A``, and then column ``j`` of ``T H`` is a genuine
      Goppa column whose label is ``(z/g)/(1/g) = z``.
"""

import itertools
import os
import random
import struct
import subprocess
import sys
import tempfile
import time

import numpy as np

from sage.all import GF, PolynomialRing, diagonal_matrix, matrix, prod, vector


F2 = GF(2)


# ---------------------------------------------------------------------------
# bitmask linear algebra over F_2 (small systems, python ints)
# ---------------------------------------------------------------------------


def mask_basis(values):
    """Reduced row-echelon basis (as python ints) of the span of ``values``."""

    pivots = {}
    for x in values:
        y = x
        while y:
            p = y.bit_length() - 1
            if p in pivots:
                y ^= pivots[p]
            else:
                pivots[p] = y
                for q in list(pivots):
                    if q != p and (pivots[q] >> p) & 1:
                        pivots[q] ^= y
                break
    return list(pivots.values())


def mask_span(basis):
    """Every element of the F_2 span of a list of ints (as a set)."""

    span = {0}
    for value in basis:
        span |= {x ^ value for x in tuple(span)}
    return span


def left_kernel(rows):
    """Coefficient vectors (ints over the row index space) summing to zero."""

    pivots = {}
    kernel = []
    for index, row in enumerate(rows):
        combo = 1 << index
        x = row
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p][0]
                combo ^= pivots[p][1]
            else:
                pivots[p] = (x, combo)
                combo = 0
                break
        if combo:
            kernel.append(combo)
    return kernel


def as_array(M):
    """A GF(2) matrix as a C-contiguous ``uint8`` numpy array."""

    return np.frombuffer(bytes(M.list()), dtype=np.uint8).reshape(
        M.nrows(), M.ncols())


def rows_as_ints(M):
    """The rows of a GF(2) matrix as python ints, bit ``j`` = column ``j``."""

    packed = np.packbits(as_array(M), axis=1, bitorder="little")
    width = packed.shape[1]
    return [int.from_bytes(packed[i].tobytes(), "little")
            for i in range(packed.shape[0])]


# ---------------------------------------------------------------------------
# parity checks
# ---------------------------------------------------------------------------


def parity_matrix(F, support, poly, m=8, t=5):
    """The m*t x n parity check built from ``1/g(z_j), z_j/g(z_j), ...``."""

    H = matrix(F2, m * t, len(support))
    for j, x in enumerate(support):
        y = F.one() / poly(x)
        xp = F.one()
        for i in range(t):
            v = (y * xp).to_integer()
            for b in range(m):
                H[m * i + b, j] = (v >> b) & 1
            xp *= x
    return H


def _column_values(H, F, j, m=8, t=5):
    """Column ``j`` of a *raw* parity check as ``z^i/g(z)``, i < t."""

    return [F.from_integer(sum(int(H[m * i + b, j]) << b for b in range(m)))
            for i in range(t)]


# ---------------------------------------------------------------------------
# minimum-weight codewords (Lee-Brickell, p = 3)
# ---------------------------------------------------------------------------


def _write_parity(H, path):
    rows, cols = H.nrows(), H.ncols()
    with open(path, "w") as stream:
        for i in range(rows):
            value = 0
            for j in range(cols):
                if H[i, j]:
                    value |= 1 << j
            stream.write("%x\n" % value)


# ---------------------------------------------------------------------------
# the quadratic relation quotient
# ---------------------------------------------------------------------------


def affine_polars(G, words, k=None):
    """Interpolate the quadratic relations of the minimum-weight codewords.

    ``G`` is a basis of the shortened code.  A minimum-weight codeword has odd
    weight, so in systematic message coordinates it lies on the affine
    hyperplane ``ell . m = 1``.  Eliminating one coordinate before
    interpolation removes the whole ``k``-dimensional ideal generated by that
    equation, leaving a quotient of dimension ``m + 1`` or ``2m + 1``.

    Optional ``k`` is verified against ``G.nrows()`` when provided.
    """

    if k is not None and k != G.nrows():
        raise RuntimeError("affine_polars: inconsistent shortened dimension")
    k, n = G.nrows(), G.ncols()
    aff = k - 1
    feats = 1 + aff + aff * (aff - 1) // 2

    def pair_index(i, j):
        if i > j:
            i, j = j, i
        return 1 + aff + i * (2 * aff - i - 1) // 2 + (j - i - 1)

    ell = vector(F2, [sum(G[r, c] for c in range(n)) for r in range(k)])
    eliminated = next(i for i in range(k) if ell[i])
    free = [i for i in range(k) if i != eliminated]

    # m = e_eliminated + E y parameterises ell . m = 1; the column space of E
    # is the tangent hyperplane ell . m = 0, which is also the basis change
    # used by every polar form below.
    E = matrix(F2, k, aff)
    for column, original in enumerate(free):
        E[original, column] = 1
        E[eliminated, column] = ell[original]

    pivots = G.pivots()
    if len(pivots) != k:
        raise RuntimeError("bad shortened-code rank")
    chart_of = [free.index(i) if i != eliminated else -1 for i in range(k)]

    M = matrix(F2, len(words), feats)
    for row, support in enumerate(words):
        support_set = set(support)
        active_index = [i for i in range(k) if pivots[i] in support_set]
        if sum(int(ell[i]) for i in active_index) % 2 != 1:
            raise RuntimeError("minimum word is outside the odd affine hyperplane")
        active = [chart_of[i] for i in active_index]
        M[row, 0] = 1
        for i in active:
            M[row, 1 + i] = 1
        for i, j in itertools.combinations(active, 2):
            M[row, pair_index(i, j)] = 1

    relations = M.right_kernel().basis_matrix()
    return E, relations, aff, pair_index, feats


def polar_forms(relations, aff, pair_index):
    """Symmetric ``aff x aff`` polar form of every relation-quotient member."""

    pairs = list(itertools.combinations(range(aff), 2))
    columns = np.fromiter((pair_index(i, j) for i, j in pairs), dtype=np.intp,
                          count=len(pairs))
    rows_idx = np.fromiter((i for i, _ in pairs), dtype=np.intp, count=len(pairs))
    cols_idx = np.fromiter((j for _, j in pairs), dtype=np.intp, count=len(pairs))
    mats = []
    for relation in relations.rows():
        block = np.frombuffer(bytes(relation.list()), dtype=np.uint8)[columns]
        dense = np.zeros((aff, aff), dtype=np.uint8)
        dense[rows_idx, cols_idx] = block
        dense |= dense.T
        mats.append(matrix(F2, aff, aff, dense.reshape(-1).tolist()))
    return mats


def pack_matrices(mats, aff):
    """Bit-packed row-major dump of ``mats`` for the native rank census."""

    nw = (aff + 63) // 64
    chunks = []
    for B in mats:
        packed = np.packbits(as_array(B), axis=1, bitorder="little")
        pad = nw * 8 - packed.shape[1]
        if pad:
            packed = np.hstack([packed, np.zeros((aff, pad), dtype=np.uint8)])
        chunks.append(np.ascontiguousarray(packed).tobytes())
    return b"".join(chunks)


def rank_census(mats, aff, finder="./rankmask"):
    """``table[mask]`` = rank of the sum of ``mats`` selected by ``mask``."""

    nb = len(mats)
    handle, path = tempfile.mkstemp(prefix="goppa-polars-", suffix=".bin")
    os.close(handle)
    try:
        with open(path, "wb") as stream:
            stream.write(pack_matrices(mats, aff))
        done = subprocess.run([finder, path, str(aff), str(nb)],
                              capture_output=True, text=True)
    finally:
        os.unlink(path)
    if done.returncode != 0:
        raise RuntimeError(f"rank census failed: {done.stderr.strip()}")
    table = [0] * (1 << nb)
    for line in done.stdout.split("\n"):
        if not line:
            continue
        mask, rank = line.split()
        table[int(mask)] = int(rank)
    return table


def punctured_cut(table, nb, maxdim, gap=4):
    """Largest rank cut whose low set is (a punctured) subspace of dim <= maxdim.

    The low-rank locus is a linear subspace with at most one exceptional
    high-rank direction, so ``low(r)`` is a punctured subspace exactly while
    ``r`` stays inside the locus and stops being one as soon as the first
    generic vector enters.  Taking the *largest* such cut needs no rank
    threshold; the ``gap`` requirement (every remaining vector is at least
    ``gap`` ranks higher) rejects the accidental small subspaces that a
    handful of generic vectors can form.  Returns ``(cut, low, basis)`` or
    ``None``.
    """

    by_rank = {}
    for mask in range(1, 1 << nb):
        by_rank.setdefault(table[mask], []).append(mask)
    ranks = sorted(by_rank)
    low = []
    best = None
    for position, rank in enumerate(ranks):
        low.extend(by_rank[rank])
        if len(low) > (1 << maxdim) + 2:
            break
        d = len(mask_basis(low))
        if d > maxdim:
            break
        if d >= 1 and (1 << d) - 2 <= len(low) <= (1 << d) - 1:
            if position + 1 < len(ranks) and ranks[position + 1] >= rank + gap:
                best = (rank, list(low))
    if best is None:
        return None
    cut, low = best
    basis = mask_basis(low)
    if not all(x in mask_span(basis) for x in low):
        return None
    return cut, low, basis


def low_rank_span(mats, aff, m=8, rng=None, finder="./rankmask", probe=8,
                  max_draws=40, patience=12, target=None,
                  log=lambda s: None):
    """The low-rank subspace, found with ``O(r^2 / log r)`` rank computations.

    The locus is a linear subspace ``L`` of dimension ``target`` (``m+1``, or
    ``m`` once the exceptional direction has been divided out) inside the
    ``r``-dimensional quotient.  A uniformly random ``s``-dimensional subspace
    ``W`` meets it in expected dimension ``s dim(L) / r``, so accumulating
    intersections until they span ``L`` costs about ``dim(L) r / s`` draws of
    ``2^s`` ranks each.  With ``s ~ log2 r`` that is a few thousand rank
    computations against ``2^r`` for the exhaustive census, which stays
    available as a fallback.
    """

    if target is None:
        target = m + 1
    nb = len(mats)
    if nb <= target:
        return None
    if rng is None:
        rng = random.Random(0x5EED)
    if nb > 24:
        probe = min(13, nb - target + 1)
        max_draws = max(max_draws, 80)
    probe = min(probe, nb - 1)
    collected = []
    stale = 0
    for _ in range(max_draws):
        current = len(mask_basis(collected))
        if current >= target:
            break
        if current >= target - 1 and stale >= patience:
            break
        basis = []
        while len(basis) < probe:
            grown = mask_basis(basis + [rng.getrandbits(nb)])
            if len(grown) > len(basis):
                basis = grown
        sub = [sum((mats[t] for t in range(nb) if value >> t & 1),
                   matrix(F2, aff, aff)) for value in basis]
        found = punctured_cut(rank_census(sub, aff, finder=finder), probe,
                              min(probe, target))
        gained = False
        if found is not None:
            for mask in found[1]:
                value = _expand(mask, basis)
                if len(mask_basis(collected + [value])) > current:
                    collected.append(value)
                    current += 1
                    gained = True
        stale = 0 if gained else stale + 1
    basis = mask_basis(collected)
    log(f"restricted census: {len(basis)}-dimensional low-rank span")
    if len(basis) < target - 1 or len(basis) > target:
        return None
    return basis


def _expand(mask, basis):
    value = 0
    for i in range(len(basis)):
        if mask >> i & 1:
            value ^= basis[i]
    return value


def _components_from_span(mats, aff, nb, basis, low, cut, m, log, target=None):
    if target is None:
        target = m + 1
    d = len(basis)
    span = mask_span(basis)
    missing = sorted(span - set(low) - {0})
    log(f"low-rank census: {len(low)} vectors at rank <= {cut}, span {d}, "
        f"{len(missing)} exceptional")

    def lift(masks):
        return [sum((mats[t] for t in range(nb) if mask >> t & 1),
                    matrix(F2, aff, aff)) for mask in masks]

    if target == m:
        if d != m or missing:
            raise RuntimeError(f"unexpected low-rank geometry ({d}, "
                               f"{len(missing)})")
        return [lift(basis)], len(low), cut
    if d == m:
        return [lift(basis)], len(low), cut
    if d != m + 1:
        raise RuntimeError(f"unexpected low-rank dimension {d}")

    components = []
    if len(missing) == 1:
        exceptional = missing[0]
        chosen = [exceptional]
        for value in basis:
            if len(mask_basis(chosen + [value])) > len(chosen):
                chosen.append(value)
        if len(chosen) == m + 1:
            components.append(lift(chosen[1:] + [exceptional]))
    elif not missing:
        components.append(lift(basis))
    # The other 510-point geometry is two m-spaces meeting only at zero.
    low_set = set(low)
    seed = low[0]
    first = {0, seed} | {x for x in low if (x ^ seed) in low_set}
    if len(first) == (1 << m):
        other = next(x for x in low if x not in first)
        second = {0, other} | {x for x in low if (x ^ other) in low_set}
        if len(second) == (1 << m) and (first & second) == {0}:
            components.append(lift(mask_basis(first)))
            components.append(lift(mask_basis(second)))
    if not components:
        raise RuntimeError("unexpected low-rank geometry")
    return components, len(low), cut


def low_rank_components(mats, aff, nprime=None, m=8, table=None,
                        finder="./rankmask", restricted=False, target=None,
                        log=lambda s: None):
    """Components of the low-rank locus as lists of ``aff x aff`` matrices."""

    if target is None:
        target = m + 1
    nb = len(mats)
    if table is None and restricted:
        span = low_rank_span(mats, aff, m=m, finder=finder, target=target,
                             log=log)
        if span is not None:
            sub = [sum((mats[t] for t in range(nb) if value >> t & 1),
                       matrix(F2, aff, aff)) for value in span]
            found = punctured_cut(rank_census(sub, aff, finder=finder),
                                  len(span), target)
            # An (target-1)-dimensional span that still hides a high-rank
            # member is a subspace of the true locus, not the locus.
            complete = (found is not None
                        and len(mask_basis(found[1])) == len(span)
                        and (len(span) == target
                             or len(found[1]) == (1 << (target - 1)) - 1))
            if complete:
                cut, low_local, _ = found
                low = [_expand(mask, span) for mask in low_local]
                return _components_from_span(mats, aff, nb, mask_basis(low),
                                             low, cut, m, log, target=target)
            log("restricted census did not close; falling back")
    if table is None:
        table = rank_census(mats, aff, finder=finder)
    found = punctured_cut(table, nb, target)
    if found is None:
        raise RuntimeError("no low-rank subspace found in the census")
    cut, low, basis = found
    return _components_from_span(mats, aff, nb, basis, low, cut, m, log,
                                 target=target)


# ---------------------------------------------------------------------------
# the exceptional direction is public
# ---------------------------------------------------------------------------


def chart_weight_form(G0):
    """Polar form of ``c -> sum_{i<j} c_i c_j`` in the chart basis, as numpy.

    ``sum_{i<j} c_i c_j = C(wt(c), 2)`` is ``1`` on every weight-``2t+1``
    codeword, so this is a relation of the minimum-weight words that is known
    *before* any label is.  Its polar form on the code is ``J + I`` (all ones
    off the diagonal), and pulling that back along ``G0`` gives
    ``G0 (J + I) G0^T``.  It is the unique high-rank member of the low-rank
    locus, so dividing it out leaves a locus of dimension exactly ``m``.
    """

    Gnp = as_array(G0)
    sums = Gnp.sum(axis=1) % 2                       # (G0 1)_a
    T = (np.outer(sums, sums) ^ (Gnp @ Gnp.T)) & 1   # G0 J G0^T + G0 I G0^T
    np.fill_diagonal(T, 0)
    return T


def _pack_row(array):
    flat = np.ascontiguousarray(array).reshape(-1)
    return int.from_bytes(np.packbits(flat, bitorder="little").tobytes(),
                          "little")


def express_in_basis(mats, target):
    """Coefficient mask of ``target`` in the F_2-basis ``mats``, or ``None``.

    The relation quotient has ``2m+1`` members at ``t = 5``, so this is a
    handful of eliminations on integers of ``aff^2`` bits.
    """

    rows = [_pack_row(as_array(B)) for B in mats]
    pivots = {}
    for index, row in enumerate(rows):
        combo = 1 << index
        x = row
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p][0]
                combo ^= pivots[p][1]
            else:
                pivots[p] = (x, combo)
                break
    combo = 0
    x = _pack_row(target)
    while x:
        p = x.bit_length() - 1
        if p not in pivots:
            return None
        x ^= pivots[p][0]
        combo ^= pivots[p][1]
    return combo


# ---------------------------------------------------------------------------
# lifting the low-rank subspace to the hidden labels
# ---------------------------------------------------------------------------


def weight_candidates(B, P, m=8, blocks=(30, 44, 60, 90)):
    """Yield candidate non-pivot diagonal vectors ``z``, cheapest block first.

    ``L0 = B + P diag(z) P^T + diag(w)`` has rank ``m``, so any ``k`` of its
    rows have a ``(k-m)``-dimensional space of dependencies ``alpha``.  For a
    dependency and a column ``j`` outside the block,
    ``sum_r alpha_r L0[r,j] = 0`` reads
    ``sum_r alpha_r B[r,j] + sum_h (alpha P[R,h]) z_h P[j,h] = 0``; with the
    linearisation ``y = (alpha P[R,:]) . z`` this is *linear* in
    ``(alpha, y)``, so the whole dependency space is the kernel of one
    ``|J| x (k + nh)`` matrix and every ``z_h`` with ``(alpha P[R,:])_h = 1``
    for some kernel vector is read off as ``y_h``.
    """

    aff = B.nrows()
    nh = P.ncols()
    Bnp = as_array(B)
    Pnp = as_array(P)
    Prow = rows_as_ints(P)

    def attempt(R):
        k = len(R)
        inside = set(R)
        rest = [j for j in range(aff) if j not in inside]
        # system[i, r]   = B[R[r], rest[i]]
        # system[i, k+h] = P[rest[i], h]
        block = np.empty((len(rest), k + nh), dtype=np.uint8)
        block[:, :k] = Bnp[np.ix_(np.array(R, dtype=np.intp),
                                  np.array(rest, dtype=np.intp))].T
        block[:, k:] = Pnp[np.array(rest, dtype=np.intp), :]
        system = matrix(F2, len(rest), k + nh, block.reshape(-1).tolist())
        known = {}
        for entry in system.right_kernel().basis_matrix().rows():
            q = 0
            for r in range(k):
                if entry[r]:
                    q ^= Prow[R[r]]
            while q:
                low = q & -q
                h = low.bit_length() - 1
                value = int(entry[k + h])
                if known.get(h, value) != value:
                    return None
                known[h] = value
                q ^= low
        if len(known) < nh:
            return None
        return vector(F2, [known[h] for h in range(nh)])

    # Support-restricted dependency solver (works even when aff < nh!)
    width = m + 1
    Brows = rows_as_ints(B)
    Prows = rows_as_ints(P)
    votes = [[0, 0] for _ in range(nh)]
    for start in range(0, aff, max(1, width // 2)):
        if all(v[0] + v[1] >= 2 for v in votes):
            break
        J = [(start + i) % aff for i in range(width)]
        jset = set(J)
        K = [i for i in range(aff) if i not in jset]
        PJ = [Prows[r] for r in J]
        PK_cols = [0] * nh
        BK_cols = [0] * width
        for idx_k, r in enumerate(K):
            br = Brows[r]
            for i, jj in enumerate(J):
                if (br >> jj) & 1:
                    BK_cols[i] |= 1 << idx_k
            tmp = Prows[r]
            while tmp:
                low = tmp & -tmp
                PK_cols[low.bit_length() - 1] |= 1 << idx_k
                tmp ^= low
        for alpha in range(1, 1 << width):
            beta = 0
            rhs = 0
            for i in range(width):
                if (alpha >> i) & 1:
                    beta ^= PJ[i]
                    rhs ^= BK_cols[i]
            active = []
            tmp = beta
            while tmp:
                low = tmp & -tmp
                active.append(low.bit_length() - 1)
                tmp ^= low
            if len(active) == 0 or len(active) > len(K) - 5:
                continue
            pivs = {}
            consistent = True
            for idx_a, h in enumerate(active):
                vec = PK_cols[h]
                comb_m = 1 << idx_a
                while vec:
                    p = vec.bit_length() - 1
                    if p in pivs:
                        vec ^= pivs[p][0]
                        comb_m ^= pivs[p][1]
                    else:
                        pivs[p] = (vec, comb_m)
                        break
                else:
                    consistent = False
                    break
            if not consistent:
                continue
            rem = rhs
            sol_mask = 0
            while rem:
                p = rem.bit_length() - 1
                if p not in pivs:
                    break
                rem ^= pivs[p][0]
                sol_mask ^= pivs[p][1]
            if rem != 0:
                continue
            for idx_a, h in enumerate(active):
                votes[h][(sol_mask >> idx_a) & 1] += 1
    if all(v[0] + v[1] > 0 for v in votes):
        yield vector(F2, [0 if v[0] >= v[1] else 1 for v in votes])
    for size in (m + 2, m + 4, 2 * m, 3 * m, 4 * m, 6 * m, 90):
        if size >= aff or size <= m:
            continue
        for start in range(0, aff, max(1, size // 2)):
            R = sorted({(start + i) % aff for i in range(size)})
            if len(R) <= m:
                continue
            candidate = attempt(R)
            if candidate is not None:
                yield candidate


def complete_diagonal(off, m=8):
    """The diagonal ``d`` with ``off + diag(d)`` of rank ``m``, deterministically.

    The left kernel of ``off[R, J]`` for a 9-set ``R`` is the row dependency of
    ``L0 = off + diag(d)`` restricted to the columns outside ``R``, hence it
    *is* that dependency, and ``d_j = sum_{r != j} alpha_r off[r,j]`` for every
    ``j`` in ``R`` with ``alpha_j = 1``.  Sliding the eight companion rows
    covers every index; no diagonal pattern is ever guessed.
    """

    aff = off.nrows()
    rows = rows_as_ints(off)
    diagonal = [None] * aff
    for j in range(aff):
        pool = [i for i in range(aff) if i != j]
        for start in range(0, len(pool) - m + 1):
            R = [j] + pool[start:start + m]
            mask = 0
            for r in R:
                mask |= 1 << r
            kernel = left_kernel([rows[r] & ~mask for r in R])
            if len(kernel) != 1 or not kernel[0] & 1:
                continue
            alpha = kernel[0]
            value = 0
            for i in range(1, m + 1):
                if alpha >> i & 1:
                    value ^= (rows[R[i]] >> j) & 1
            diagonal[j] = value
            break
        if diagonal[j] is None:
            raise RuntimeError("diagonal is not identifiable")
    return vector(F2, diagonal)


def gf2_solve(A, ncols, R):
    """Solve ``A w = R`` over F_2 for every column of ``R``.

    ``A`` is bit-packed ``uint8`` (bit ``c`` of column byte ``c >> 3``), ``R``
    holds one 0/1 column per right-hand side.
    """

    A = np.ascontiguousarray(A.copy())
    R = np.ascontiguousarray(R.copy())
    nrows = A.shape[0]
    pivots = []
    row = 0
    for c in range(ncols):
        byte, bit = c >> 3, c & 7
        hits = np.flatnonzero((A[row:, byte] >> bit) & 1)
        if hits.size == 0:
            continue
        p = row + int(hits[0])
        if p != row:
            A[[row, p]] = A[[p, row]]
            R[[row, p]] = R[[p, row]]
        sel = np.flatnonzero((A[:, byte] >> bit) & 1)
        sel = sel[sel != row]
        if sel.size:
            A[sel] ^= A[row]
            R[sel] ^= R[row]
        pivots.append(c)
        row += 1
        if row == nrows:
            break
    if len(pivots) != ncols:
        raise RuntimeError("weight system is rank deficient")
    if np.any(R[row:]):
        raise RuntimeError("weight system is inconsistent")
    out = np.zeros((R.shape[1], ncols), dtype=np.uint8)
    for i, c in enumerate(pivots):
        out[:, c] = R[i, :]
    return out


def lift_all(Bs, Gs, nonpivots, pivots, m=8, log=lambda s: None):
    """Lift every polar basis form into one common ``m``-dimensional space.

    Returns ``(W, blocks, I)``: the weight matrix whose columns are the
    ``z_j^2`` bit vectors, the ``m`` lifted forms restricted to the ``m``
    positions ``I`` whose labels form a basis, and ``I`` itself.
    """

    aff = Gs.nrows()
    nprime = Gs.ncols()
    P = Gs.matrix_from_columns(nonpivots)
    Pnp = as_array(P)
    Gnp = as_array(Gs)

    last = None
    for z in weight_candidates(Bs[0], P, m=m):
        off = as_array(Bs[0]) ^ ((Pnp * np.frombuffer(bytes(z.list()),
                                                      dtype=np.uint8))
                                 @ Pnp.T & 1)
        np.fill_diagonal(off, 0)
        try:
            diagonal = complete_diagonal(matrix(F2, aff, aff,
                                               off.reshape(-1).tolist()), m=m)
        except RuntimeError as exc:
            last = exc
            continue
        dense = off.copy()
        dense[np.arange(aff), np.arange(aff)] = np.frombuffer(
            bytes(diagonal.list()), dtype=np.uint8)
        L0 = matrix(F2, aff, aff, dense.reshape(-1).tolist())
        if L0.rank() != m:
            last = RuntimeError(f"lifted form has rank {L0.rank()}")
            continue

        w0 = np.zeros(nprime, dtype=np.uint8)
        w0[np.array(nonpivots, dtype=np.intp)] = np.frombuffer(
            bytes(z.list()), dtype=np.uint8)
        w0[np.array(pivots, dtype=np.intp)] = (
            dense.diagonal() ^ ((Pnp * w0[np.array(nonpivots, dtype=np.intp)])
                                @ Pnp.T & 1).diagonal())
        if not np.array_equal(dense, as_array(Bs[0]) ^ ((Gnp * w0) @ Gnp.T & 1)):
            last = RuntimeError("invalid first diagonal completion")
            continue
        break
    else:
        raise RuntimeError(f"no weight vector lifts: {last!r}")

    kernel = L0.right_kernel().basis_matrix()
    log(f"common rowspace {L0.rank()} kernel {kernel.dimensions()}")

    # Lifted forms must have their rows inside the common space:
    #   (Bs[t] + Gs diag(w) Gs^T) kernel^T = 0.
    Knp = as_array(kernel)
    KG = (Knp @ Gnp) & 1
    system = np.ascontiguousarray(
        (Gnp[:, None, :] & KG[None, :, :]).reshape(aff * kernel.nrows(), nprime))
    Spacked = np.ascontiguousarray(
        np.packbits(system, axis=1, bitorder="little"))

    weights = [w0]
    dense_forms = [dense]
    # Every remaining slice solves the *same* overdetermined system
    # ``(Gs diag(w) Gs^T) K^T = B K^T`` with a different right-hand side, so the
    # elimination is done once for all of them instead of once per slice.
    if len(Bs) > 1:
        rhs = np.empty((aff * kernel.nrows(), len(Bs) - 1), dtype=np.uint8)
        for position, B in enumerate(Bs[1:]):
            rhs[:, position] = ((as_array(B) @ Knp.T) & 1).reshape(-1)
        solved = gf2_solve(Spacked, nprime, np.ascontiguousarray(rhs))
        for position, B in enumerate(Bs[1:]):
            weight = solved[position]
            weights.append(weight)
            dense_forms.append(as_array(B) ^ ((Gnp * weight) @ Gnp.T & 1))

    I = basis_indices(dense_forms, m=m)
    blocks = [[[int(dense_forms[t][I[i], I[j]]) for j in range(m)]
               for i in range(m)] for t in range(len(dense_forms))]
    W = matrix(F2, len(weights), nprime,
               np.stack(weights).reshape(-1).tolist())
    return W, blocks, I


def basis_indices(forms, m=8):
    """Greedy ``m``-set whose labels are an F_2 basis of the abstract field.

    ``forms[t][a,b]`` is bit ``t`` of the product of the labels at ``a`` and
    ``b``, so the columns ``(forms[t][u,c])_t`` for ``c`` in a candidate set
    are the bit vectors of ``zeta_u zeta_c``.  Multiplication by the nonzero
    ``zeta_u`` is an F_2-linear bijection, so those columns are independent
    exactly when the labels themselves are.  Testing that directly (rather
    than testing an ``m x m`` principal block of a single bit slice, whose
    restriction of the form can be degenerate) makes the greedy step exact.
    """

    aff = forms[0].shape[0]
    chosen = []
    for j in range(aff):
        if len(chosen) == m:
            break
        trial = chosen + [j]
        reference = trial[0]
        block = matrix(F2, m, len(trial),
                       [int(forms[t][reference, c]) for t in range(m)
                        for c in trial])
        if block.rank() == len(trial):
            chosen = trial
    if len(chosen) != m:
        raise RuntimeError("the lifted labels do not span the abstract field")
    return chosen


def make_field(W, blocks, I, F, m=8):
    """Read the field structure off the lifted forms and label every column.

    ``blocks[t][i][j]`` is bit ``t`` of the product of the labels at ``I[i]``
    and ``I[j]``, i.e. the multiplication table of an abstract
    ``m``-dimensional algebra.  An element whose characteristic polynomial on
    the regular representation is irreducible of degree ``m`` generates that
    algebra, so the algebra is the field of ``2^m`` elements and sending the
    generator to a root of the same polynomial in the public key's field is a
    field isomorphism.  This replaces enumerating all ``2^m - 1`` elements.
    """

    zero = matrix(F2, m, m)
    slices = [matrix(F2, m, m, [int(blocks[t][i][j]) for i in range(m)
                                for j in range(m)]) for t in range(m)]
    # R[t][j] = bit t of the label at I[j]  (for the unit position)
    R = None
    for unit in range(m):
        candidate = matrix(F2, m, m, [int(blocks[t][unit][j])
                                      for t in range(m) for j in range(m)])
        if candidate.rank() == m:
            R = candidate
            break
    if R is None:
        raise RuntimeError("bilinear algebra has no visible unit isotope")
    Ri = R.inverse()
    # S_i[t][j] = bit t of (label at I[i]) * (label at I[j]) = blocks[t][i][j]
    S = [matrix(F2, m, m, [int(blocks[t][i][j]) for t in range(m)
                           for j in range(m)]) for i in range(m)]
    unit_vector = vector(F2, [1 if i == unit else 0 for i in range(m)])
    basis = [vector(F2, [1 if i == j else 0 for i in range(m)])
             for j in range(m)]

    def Q(x):
        out = zero
        for i in range(m):
            if x[i]:
                out = out + S[i]
        return out

    def product(a, b):
        return (Ri * Q(a)) * b

    if any(product(unit_vector, e) != e or product(e, unit_vector) != e
           for e in basis):
        raise RuntimeError("candidate algebra has no identity")

    def regular(x):
        return Ri * Q(x)

    for i in range(m):
        for j in range(m):
            if regular(basis[i]) * regular(basis[j]) != regular(
                    product(basis[i], basis[j])):
                raise RuntimeError("candidate algebra is not associative")

    modulus = None
    gamma = None
    for value in range(2, 1 << m):
        x = vector(F2, [(value >> i) & 1 for i in range(m)])
        polynomial = regular(x).charpoly()
        if polynomial.degree() == m and polynomial.is_irreducible():
            modulus, gamma = polynomial, x
            break
    if gamma is None:
        raise RuntimeError("candidate algebra is not a field")

    powers = [unit_vector]
    for _ in range(m - 1):
        powers.append(product(powers[-1], gamma))
    change = matrix(F2, m, m, [int(powers[j][i]) for i in range(m)
                               for j in range(m)])
    if change.rank() != m:
        raise RuntimeError("field generator is degenerate")
    to_powers = change.inverse()

    PR = PolynomialRing(F, "X")
    roots = PR([F(int(c)) for c in modulus.list()]).roots()
    if not roots:
        raise RuntimeError("field modulus has no root in the public field")
    root = roots[0][0]
    root_powers = [F.one()]
    for _ in range(m - 1):
        root_powers.append(root_powers[-1] * root)

    if W.nrows() < m:
        raise RuntimeError("weight matrix is too short")
    labels = []
    for j in range(W.ncols()):
        column = vector(F2, [int(W[t, j]) for t in range(m)])
        coeff = to_powers * (Ri * column)
        labels.append(sum((root_powers[i] for i in range(m) if coeff[i]),
                          F.zero()))
    return labels


def goppa_from_word(F, labels, word, t=5):
    """``g`` from one minimum-weight codeword, using ``Q'_S = g^2``."""

    PR = PolynomialRing(F, "X")
    X = PR.gen()
    sigma = prod((X - labels[j] for j in word), PR.one())
    derivative = sigma.derivative()
    # Q' = sum_{i even} e_i X^{2t-i} = g^2, so g[X^{t-j}] = (Q'[X^{2(t-j)}])^{1/2}.
    half = 2 ** (F.degree() - 1)
    return PR([derivative[2 * i] ** half for i in range(t + 1)])


def extend_support(H, F, A, labels_A, poly, m=8, t=5):
    """Label the positions outside ``A`` from the published parity check.

    ``H`` is any basis of the published row space, while the key defines one
    particular basis: the ``t`` blocks of ``z^i/g(z)`` bit slices.  The two
    differ by an invertible ``mt x mt`` matrix ``T``.  The recovered labels on
    ``A`` pin ``T`` down (``H_cand = T H`` on ``mt`` independent columns of
    ``A``), and then every remaining column of ``T H`` is a genuine Goppa
    column whose label is simply ``(z/g)/(1/g) = z``.

    Any representative of the key works here: ``T`` absorbs the Frobenius class
    and the leading coefficient that the lift leaves free, because all of them
    generate the same row space.
    """

    n = H.ncols()
    Hcand = parity_matrix(F, labels_A, poly, m, t)
    known = H.matrix_from_columns(A)
    independent = known.pivots()
    if len(independent) != m * t:
        raise RuntimeError("degenerate shortening")
    T = (Hcand.matrix_from_columns(independent)
         * known.matrix_from_columns(independent).inverse())
    if T * known != Hcand:
        raise RuntimeError("recovered key does not reproduce the published row "
                           "space on A")
    Hraw = T * H

    labels = [None] * n
    for index, j in enumerate(A):
        labels[j] = labels_A[index]
    for j in range(n):
        if labels[j] is not None:
            continue
        values = _column_values(Hraw, F, j, m, t)
        label = values[1] / values[0]
        if (any(values[i] / values[0] != label ** i for i in range(1, t))
                or poly(label) * values[0] != 1):
            raise RuntimeError("cannot extend the support")
        labels[j] = label
    return labels
