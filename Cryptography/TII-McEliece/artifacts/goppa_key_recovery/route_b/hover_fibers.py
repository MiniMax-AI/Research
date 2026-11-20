"""Stage 4 of Route B (Frobenius-reduced HOVER): truncated Segre-Veronese fiber extraction."""
from __future__ import annotations
import time
from labels import FibreFrame, tangent, fibre_of

NFIB = 450  # fibers extracted (hover_recover requires >= 240 + 110)


def extract_fibers(cols, Fq, nfib=NFIB, verbose=True):
    """Segre-Veronese fibers of one hidden direction, for the first `nfib` coordinates."""
    tot = cols.ncols()
    use = min(nfib, tot)
    collist = [cols.column(j) for j in range(use)]
    FF = None
    for off in (0, 1, 2, 3, 5, 8):
        if off + 4 > use:
            break
        try:
            FF = FibreFrame(collist, Fq, seeds=(off, off + 1, off + 2, off + 3))
            break
        except Exception as exc:
            if verbose:
                print(f"[fibers] FibreFrame seeds at {off} failed ({exc}); shifting", flush=True)
    if FF is None:
        return None, None
    t0 = time.time()
    fib = [None] * tot
    for j, c in enumerate(collist):
        T = tangent(FF.Ib, FF.idx, c, 9, Fq)
        if T.dimension() != 4:
            continue
        F = fibre_of(T.basis_matrix(), FF.Ib, FF.idx, 9, Fq)
        if F is not None and F.nrows() == 3:
            fib[j] = F
    good = sum(1 for f in fib if f is not None)
    if verbose:
        print(
            f"[fibers] {good}/{use} fibers of dim 3 (of {tot} coordinates), tangent spaces "
            f"dim 4, quadric ideal dim {FF.Ib.nrows()}   {time.time()-t0:.1f}s",
            flush=True,
        )
    if good < 300:
        return None, None
    return FF, fib
