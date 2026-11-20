"""Small bit-packed GF(2) helpers (numpy), used by the HOVER selector stage.

Sage's GF(2) matrices are M4RI-backed and excellent for large *square-ish* problems,
but constructing one from a Python list costs ~0.3 us per entry, which dominates when
the matrices are small and numerous (the selector builds a 58500 x 100 system and a
2701 x 740 one).  These helpers keep everything in uint64 bitsets and never leave
numpy, which removes that conversion entirely.
"""
from __future__ import annotations

import numpy as np


def pack(A8):
    """(r, c) uint8 -> (r, ceil(c/64)) uint64, bit j of a row at word j>>6."""
    r, c = A8.shape
    W = (c + 63) // 64
    out = np.zeros((r, W), dtype=np.uint64)
    for j in range(c):
        col = A8[:, j]
        if col.any():
            out[:, j >> 6] |= col.astype(np.uint64) << np.uint64(j & 63)
    return out


def rref(A):
    """In-place reduced row echelon form of a packed (r, W) uint64 matrix.

    Returns the list of pivot column indices.  Rows are permuted so that the pivot of
    pivots[k] sits in row k.
    """
    r, W = A.shape
    ncols = 64 * W
    piv = []
    prow = 0
    for j in range(ncols):
        if prow >= r:
            break
        w = j >> 6
        b = np.uint64(j & 63)
        nz = np.nonzero((A[prow:, w] >> b) & np.uint64(1))[0]
        if nz.size == 0:
            continue
        i = prow + int(nz[0])
        if i != prow:
            A[[prow, i]] = A[[i, prow]]
        mask = ((A[:, w] >> b) & np.uint64(1)).astype(bool)
        mask[prow] = False
        if mask.any():
            A[mask] ^= A[prow]
        piv.append(j)
        prow += 1
    return piv


def row_basis_indices(A8):
    """Indices of a maximal independent set of rows of an (r, c) uint8 matrix.

    Equivalent to `Matrix(GF(2), A8).transpose().pivots()`, without the conversion.
    """
    A = pack(np.ascontiguousarray(A8))
    r = A.shape[0]
    piv = []
    used = np.zeros(r, dtype=bool)
    prow = 0
    ncols = A8.shape[1]
    for j in range(ncols):
        if prow >= r:
            break
        w = j >> 6
        b = np.uint64(j & 63)
        nz = np.nonzero(((A[:, w] >> b) & np.uint64(1)).astype(bool) & ~used)[0]
        if nz.size == 0:
            continue
        i = int(nz[0])
        used[i] = True
        mask = ((A[:, w] >> b) & np.uint64(1)).astype(bool)
        mask[i] = False
        if mask.any():
            A[mask] ^= A[i]
        piv.append(i)
    return piv


def right_kernel(A8):
    """Basis of {x : A8 x = 0} for an (r, c) uint8 matrix; returns (k, c) uint8."""
    A8 = np.ascontiguousarray(A8)
    r, c = A8.shape
    A = pack(A8)
    piv = rref(A)
    pivset = set(piv)
    free = [j for j in range(c) if j not in pivset]
    out = np.zeros((len(free), c), dtype=np.uint8)
    for s, f in enumerate(free):
        out[s, f] = 1
        col = (A[:len(piv), f >> 6] >> np.uint64(f & 63)) & np.uint64(1)
        idx = np.nonzero(col)[0]
        if idx.size:
            out[s, [piv[k] for k in idx]] = 1
    return out
