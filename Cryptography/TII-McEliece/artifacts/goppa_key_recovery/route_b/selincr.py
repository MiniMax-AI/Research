"""Incremental computation of the compatible-endomorphism algebra E_beta.

The selector's `algebra()` imposes, for EVERY z in ker T_e, the rho = m p^2 equations

    sum_{r,a,b} T_e[s,(r,a)] B[a,b] z[r*d+b] = 0 ,   s = 1..rho,          (eq. 11)

i.e. it builds an (nz * rho) x d^2 matrix over F2 and takes its right kernel, where
nz = dim ker T_e = kappa * d - rho = 650 for (m, t) = (10, 11).

But E_beta is *defined* as the kernel of the full stack, and dim E_beta = m, so the
stack has rank exactly d^2 - m = 90:  rho = 90 independent equations already pin the
kernel down.  Since each z contributes rho = 90 rows, a single z generically gives the
whole rank.  We therefore grow the stack one z at a time and stop as soon as the kernel
has dimension m.  Because E_beta is contained in the kernel of every sub-stack, a
sub-stack kernel of dimension exactly m *equals* E_beta -- the stopping test is exact,
not heuristic.  The result is then re-checked against a few further z's (cheap: d^2
unknowns) so that a wrong early stop cannot pass silently.

Cost: (nz * rho) x d^2  ->  (k * rho) x d^2 with k = 1..4 in practice, i.e. the
nz = kappa * d - rho factor drops out of the dominant step.
"""
from __future__ import annotations
import numpy as np


def algebra_incremental(sel, m, max_blocks=64, verbose=False):
    """E_beta of `sel` using as few kernel vectors as possible.  Returns the basis list."""
    from sage.all import Matrix
    from bitsolve import row_basis_indices
    F2 = sel.F2
    d, kap = sel.d, sel.kappa
    T = sel.T
    rho = T.nrows()
    Tnp = np.array(T.rows(), dtype=np.uint8).reshape(rho, kap, d)
    U = np.ascontiguousarray(Tnp.transpose(0, 2, 1)).reshape(rho * d, kap)
    Z = T.right_kernel().basis_matrix()          # nz x (kap*d)
    nz = Z.nrows()

    def zblock(k):
        """(rho, d*d) coefficient block of the k-th kernel vector.  Rows of Z are
        converted one at a time: only the handful that is actually used is needed,
        so the nz = C(kappa,2) - rho factor never enters the conversion cost."""
        zk = np.array(list(Z.row(k)), dtype=np.uint8).reshape(kap, d)
        coef = (np.einsum('xr,rb->xb', U, zk, dtype=np.int64) & 1)
        return np.ascontiguousarray(coef.reshape(rho, d * d)).astype(np.uint8)

    K = None
    used = 0
    stack = None
    for k in range(1, min(max_blocks, nz) + 1):
        block = zblock(k - 1)
        stack = block if stack is None else np.vstack([stack, block])
        used += 1
        idx = row_basis_indices(stack)
        red = stack[np.array(idx, dtype=np.int64)]
        E = Matrix(F2, len(idx), d * d, red.flatten().tolist())
        K = E.right_kernel().basis_matrix()
        if verbose:
            print(f"    z-blocks {used}: rank {len(idx)}, dim E {K.nrows()}", flush=True)
        if K.nrows() <= m:
            break
    if K is None or K.nrows() != m:
        return None, used

    # exactness re-check: the algebra found must annihilate further, unused z's
    Bcoef = np.array([list(r) for r in K.rows()], dtype=np.uint8).reshape(m, d, d)
    check = 0
    for k in range(used, min(used + 8, nz)):
        b3 = zblock(k).reshape(rho, d, d)
        # sum_{a,b} block[s,(a,b)] * B[a,b] must vanish for every s and every basis B
        if int(((np.einsum('sab,mab->sm', b3, Bcoef, dtype=np.int64)) & 1).sum()) != 0:
            return None, used
        check += 1
    sel.Edim = K.nrows()
    sel.Ebasis = [Matrix(F2, d, d, list(row)) for row in K.rows()]
    sel.Echeck = check
    return sel.Ebasis, used
