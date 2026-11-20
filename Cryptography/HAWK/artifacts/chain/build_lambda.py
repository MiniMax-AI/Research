#!/usr/bin/env python3
"""Stages 1-4 of the HAWK-256 key recovery, from the public key alone.

    decode(pk) -> Q
    tau-cocycle constraint matrix cBr (4n x 2n) -> saturated Z-kernel = Z-basis of Lambda_tau
    Lambda' = 2*Lambda_tau + Z*w0 ,  w0 = I (mod 2)          [index asserted == 2^(n-1)]
    exact-Cholesky embedding + flatter + exact transform recovery -> a flat integer Gram

and it checks the two exact certificates that make the construction correct:

    log2 det(Gram_tau)            == HawkParams.logdet_gram_tau   (Prop 4.5, ePrint 2026/1593)
    index(Lambda' in Lambda_tau)  == 2^(n-1)

Runs from the artifacts/ directory (auto-fetching hawk_recover/ on first run if absent).

    python3 chain/build_lambda.py            # reuse the precomputed Lambda_tau basis (~1 min)
    python3 chain/build_lambda.py --full     # rebuild the integer kernel from scratch
"""
from __future__ import annotations

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np
from sage.all import ZZ, matrix, vector, GF

from common import DATA, WORK, cpu, ensure_work, log, read_mat, resolve, write_mat
from hawk_recover import decode, params as P
from hawk_recover.cocycle import _build_cBr, _Q_to_array, _embed_LB, _gram_tau


def build_kernel(pk, p):
    """Saturated Z-kernel of the tau-cocycle constraint system, via HNF + flatter.

    Replaces Sage's LLL (~6 min in ePrint 2026/1593) with flatter reduction
    on the HNF kernel (~143 CPU-s total in Stage 2).
    """
    from scipy.linalg import qr as scipy_qr
    ensure_work()
    cBr = _build_cBr(_Q_to_array(pk), p.dL, p.dK)
    piv = scipy_qr(cBr.T.astype(np.float64), pivoting=True, mode="economic")[-1]
    rows = sorted(int(i) for i in piv[:p.dL])
    MZ = matrix(ZZ, cBr[rows, :].tolist())
    log(f"cBr {cBr.shape}, max|cBr| = {int(np.abs(cBr).max())}; MZ {MZ.nrows()}x{MZ.ncols()}")
    H, U = MZ.transpose().hermite_form(transformation=True)
    nz = sum(1 for i in range(H.nrows()) if not H.row(i).is_zero())
    K = U[nz:, :]
    assert (MZ * K.transpose()).is_zero(), "kernel check failed"
    log(f"HNF kernel: {K.nrows()} rows, height {K.height().nbits()} bits")
    write_mat(os.path.join(WORK, "Khnf.txt"), K)
    subprocess.run(["flatter", os.path.join(WORK, "Khnf.txt"),
                    os.path.join(WORK, "Kflat.txt")], check=True)
    Kf = read_mat(os.path.join(WORK, "Kflat.txt"))
    log(f"flatter-reduced kernel: height {Kf.height().nbits()} bits")
    return Kf


def lam_prime(LB, n):
    """Basis (in LB coordinates) of Lambda' = 2*Lambda_tau + Z*w0, w0 = I mod 2.

    One F_2 solve on the public basis; the determinant assertion below is the
    exact index certificate 2^(n-1).
    """
    Ivec = [0] * (4 * n)
    Ivec[0] = 1
    Ivec[3 * n] = 1
    A2 = matrix(GF(2), [[int(LB[i, j]) % 2 for j in range(4 * n)] for i in range(n)])
    b2 = vector(GF(2), [Ivec[j] % 2 for j in range(4 * n)])
    u = [int(x) for x in A2.solve_left(b2)]
    j0 = u.index(1)
    H = matrix(ZZ, n, n)
    for i in range(n):
        if i == j0:
            for jj in range(n):
                H[i, jj] = u[jj]
        else:
            H[i, i] = 2
    assert H.det().abs() == 2 ** (n - 1), f"index {H.det().abs()} != 2^{n-1}"
    return H


def chol_flatten(G, margin=45):
    """Cholesky-embed G, LLL it with flatter, and return the exact reduced Gram.

    G = L L^T with L lower triangular over a 1600-bit real field; B = round(2^s L)
    is an integer matrix with B B^T = 2^(2s) G + O(2^s).  Reducing B with flatter
    and recovering the unimodular transform U by exact backward substitution
    (B is lower triangular) gives the reduced Gram G1 = U G U^T exactly.

    Preconditions the rank-256 cocycle sublattice whose raw Gram matrix has ~40-bit
    skew, avoiding Gram-mode LLL precision stalls.
    """
    from sage.all import RealField, Integer
    ensure_work()
    n = G.nrows()
    RF = RealField(1600)
    L = matrix(RF, G).cholesky()
    s = int((-min(abs(L[i, i]) for i in range(n)).log(2)).ceil()) + margin
    B = matrix(ZZ, n, n, [[Integer((L[i, j] * RF(2) ** s).round()) for j in range(n)]
                          for i in range(n)])
    log(f"Cholesky embedding: s = {s}, entries {B.height().nbits()} bits")
    write_mat(os.path.join(WORK, "Bc.txt"), B)
    subprocess.run(["flatter", os.path.join(WORK, "Bc.txt"),
                    os.path.join(WORK, "Bcr.txt")], check=True)
    Br = read_mat(os.path.join(WORK, "Bcr.txt"))
    Bl = [[int(B[i, j]) for j in range(n)] for i in range(n)]
    U = matrix(ZZ, n, n)
    for r in range(n):
        b = [int(Br[r, j]) for j in range(n)]
        u = [0] * n
        for i in range(n - 1, -1, -1):
            acc = b[i] - sum(u[j] * Bl[j][i] for j in range(i + 1, n))
            q, rem = divmod(acc, Bl[i][i])
            assert rem == 0, "non-exact division recovering the transform"
            u[i] = q
        for j in range(n):
            U[r, j] = u[j]
    assert U * B == Br, "U*B != Br"
    G1 = U * G * U.transpose()
    log(f"flatter-reduced Gram: height {G1.height().nbits()} bits, diagonal "
        f"{min(int(G1[i,i]) for i in range(n))} .. {max(int(G1[i,i]) for i in range(n))}")
    write_mat(os.path.join(WORK, "G1.txt"), G1)
    L1 = matrix(RF, G1).cholesky()
    Bs = matrix(ZZ, n, n, [[Integer((L1[i, j] * RF(2) ** 21).round()) for j in range(n)]
                           for i in range(n)])
    write_mat(os.path.join(WORK, "Bs.txt"), Bs)
    log(f"s = 21 sieve embedding Bs.txt written ({Bs.height().nbits()} bits)")
    return G1, s


def main(argv):
    full = "--full" in argv
    p = P.HAWK256
    n = p.n
    pk = decode.decode(open(os.path.join(DATA, "pk.bin"), "rb").read(), p)
    log(f"decoded HAWK-256 public key ({len(pk.raw)} bytes), n = {n}")

    if full:
        Kr = build_kernel(pk, p)
    else:
        Kr = read_mat(resolve("Kflat.txt"))
        write_mat(os.path.join(ensure_work(), "Kflat.txt"), Kr)
        log(f"using the vendored Lambda_tau basis ({Kr.nrows()}x{Kr.ncols()}, "
            f"height {Kr.height().nbits()} bits)")

    LB = _embed_LB(Kr, p.dL, p.dK)
    Gt = _gram_tau(LB, p.dL)
    dl2 = float(Gt.det().abs().log(2))
    log(f"log2 det(Gram_tau) = {dl2:.2f}   (HawkParams.logdet_gram_tau = "
        f"{p.logdet_gram_tau:.2f})")
    assert abs(dl2 - p.logdet_gram_tau) < 0.01, "Gram_tau determinant mismatch"

    H = lam_prime(LB, n)
    log(f"Lambda' = 2*Lambda_tau + Z*w0 built; index = 2^{n-1} (asserted)")
    Gp = H * Gt * H.transpose()
    write_mat(os.path.join(ensure_work(), "Gp.txt"), Gp)
    log(f"Gram(Lambda'): {Gp.nrows()}x{Gp.ncols()}, height {Gp.height().nbits()} bits, "
        f"log2 det = {float(Gp.det().abs().log(2)):.2f}")

    G1, s = chol_flatten(Gp)
    d = [int(G1[i, i]) for i in range(n)]
    log(f"lambda_1^2 target = 2n = {p.lam1sq}; best row after LLL is "
        f"{min(d)} = {min(d)/p.lam1sq:.0f}x the target")

    print()
    print("CERTIFICATES")
    print(f"  log2 det(Gram_tau)           = {dl2:.2f}  (expected {p.logdet_gram_tau:.2f})")
    print(f"  index(Lambda' in Lambda_tau) = 2^{n-1}")
    print(f"  Gram(Lambda') log2 det       = {float(Gp.det().abs().log(2)):.2f}  "
          f"(= {p.logdet_gram_tau:.0f} + 2*{n-1})")
    print(f"  total CPU-core-seconds       = {cpu():.1f}")
    print("RESULT: SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
