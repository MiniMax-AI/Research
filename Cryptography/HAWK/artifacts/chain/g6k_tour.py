#!/usr/bin/env python3
"""Stage 5: single-thread g6k progressive pump-and-jump BKZ sieve tour.

    python3 chain/g6k_tour.py <basis.txt> <s> <lam1sq> [threads] [schedule]

e.g.  python3 chain/g6k_tour.py Bs.txt 21 512 1 40,50,60,70,80

Runs g6k's pump_n_jump_bkz_tour across the specified sieve dimensions on the
Cholesky-embedded basis Bs.txt (avoiding the fixed-precision Babai nearest-plane
stalls observed in enumeration-based fpylll BKZ at blocksize >= 5), writes the
reduced basis to Bred.sobj in the working directory, and prints RESULT: SUCCESS
once row 0 reaches squared norm <= 1.05 * 2^(2s) * lam1sq.
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from common import WORK, cpu, ensure_work, read_mat

T0 = time.time()


def L(m):
    print(f"[{time.time()-T0:7.1f}s cpu {cpu():7.1f}s] {m}", flush=True)


def main(argv):
    from sage.all import ZZ, matrix
    from fpylll import IntegerMatrix, FPLLL
    from g6k import Siever
    from g6k.siever_params import SieverParams
    from g6k.algorithms.bkz import pump_n_jump_bkz_tour
    from g6k.utils.stats import dummy_tracer

    bfile = argv[0]
    s = int(argv[1])
    target = int(argv[2])
    threads = int(argv[3]) if len(argv) > 3 else 1
    sched = [int(x) for x in argv[4].split(",")] if len(argv) > 4 else [40, 50, 60, 70, 80]

    FPLLL.set_threads(threads)
    B = read_mat(bfile)
    n = B.nrows()
    tgt = (2 ** (2 * s)) * target
    A = IntegerMatrix(n, n)
    for i in range(n):
        for j in range(n):
            A[i, j] = int(B[i, j])

    g = Siever(A, SieverParams(threads=threads))
    g.lll(0, n)
    g.M.update_gso()
    L(f"after LLL ratio={g.M.get_r(0, 0) / tgt:.6g}")

    out = os.path.join(ensure_work(), "Bred.sobj")
    for bs in sched:
        t = time.time()
        try:
            pump_n_jump_bkz_tour(g, dummy_tracer, bs, jump=1)
        except Exception as e:
            L(f"tour bs={bs} FAIL {type(e).__name__} {str(e)[:80]}")
            continue
        g.M.update_gso()
        L(f"tour bs={bs} {time.time()-t:.1f}s ratio={g.M.get_r(0, 0) / tgt:.6g}")
        if g.M.get_r(0, 0) <= tgt * 1.05:
            L("TARGET REACHED")
            matrix(ZZ, n, n, [[A[i][j] for j in range(n)] for i in range(n)]).save(out)
            print(f"wrote {out}", flush=True)
            print("RESULT: SUCCESS", flush=True)
            return 0
    print("RESULT: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
