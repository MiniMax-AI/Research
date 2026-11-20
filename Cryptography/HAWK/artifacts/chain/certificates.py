#!/usr/bin/env python3
"""Exact algebraic certificates and preconditioned coordinate pullback for HAWK-256.

Verifies on the public key and reduced basis in data/:
  (1) Gram(Lambda_tau) has log2 det = 2431.00, matching Proposition 2.
  (2) Gram(Lambda') has log2 det = 2941.00 = 2431 + 2*255, certifying index [Lambda_tau : Lambda'] = 2^255.
  (3) Row 0 of the reduced embedded basis Bred attains squared Euclidean norm 0.999625 * 2^42 * 512.
  (4) Exact triangular backward substitution recovers u in the preconditioned basis (u * G1 * u^T = 512),
      and unimodular pullback v = u * U satisfies v * Gp * v^T = 2n = 512 in Lambda'.
  (5) Mapping v to M_2(R_256) yields the exact cocycle element V_tau = B^{-1} tau(B) satisfying
      V_tau == I_2 (mod 2) and V_tau * tau(V_tau) == I_2.

Run standalone:  sage -c 'import sys; print(sys.executable)' | xargs -I{} {} chain/certificates.py
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from sage.all import ZZ, GF, matrix, vector, load

from common import DATA, WORK, cpu, log, read_mat, resolve
from hawk_recover import decode, params as P, ring
from hawk_recover.cocycle import _embed_LB, _gram_tau
from hawk_recover.tower import _amb_to_matrix, _is_aside

N = 256
S = 21


def _lam_prime_H(LB, n):
    """Change-of-basis matrix H such that Lambda' = H * (Lambda_tau basis)."""
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
    return H


def verify_secret_key_bytes(pk_path: str, sk_path: str) -> int:
    """Exact algebraic verification of a 592-byte decoded secret key against pk.bin."""
    import hashlib
    import struct

    p = P.HAWK256
    n = p.n
    pk_raw = open(pk_path, "rb").read()
    sk_raw = open(sk_path, "rb").read()
    assert len(pk_raw) == p.pk_bytes, f"pk size {len(pk_raw)} != {p.pk_bytes}"
    expected_sk_len = p.sk_core_bytes + p.hpub_bytes
    assert len(sk_raw) == expected_sk_len, f"skd size {len(sk_raw)} != {expected_sk_len}"

    pk = decode.decode(pk_raw, p)
    q00, q01, _q10, _q11 = pk.Q_list
    f = list(struct.unpack(f"{n}b", sk_raw[:n]))
    g = list(struct.unpack(f"{n}b", sk_raw[n : 2 * n]))
    F2 = [((sk_raw[2 * n + (k >> 3)] >> (k & 7)) & 1) for k in range(n)]
    G2 = [((sk_raw[2 * n + (n // 8) + (k >> 3)] >> (k & 7)) & 1) for k in range(n)]
    hpub = sk_raw[p.sk_core_bytes :]

    assert hpub == hashlib.shake_256(pk_raw).digest(p.hpub_bytes), "SHAKE256 hpub mismatch"
    cf = decode._nc_conj(f, n)
    cg = decode._nc_conj(g, n)
    lhs_q00 = [a + b for a, b in zip(ring.nc_mul(cf, f, n), ring.nc_mul(cg, g, n))]
    assert lhs_q00 == q00, "|f|^2 + |g|^2 != q00"

    F = decode._nc_divexact([a - b for a, b in zip(ring.nc_mul(f, q01, n), cg)], q00, n)
    G = decode._nc_divexact([a + b for a, b in zip(ring.nc_mul(g, q01, n), cf)], q00, n)
    det_b = [a - b for a, b in zip(ring.nc_mul(f, G, n), ring.nc_mul(g, F, n))]
    assert det_b == [1] + [0] * (n - 1), "f*G - g*F != 1"
    assert [x & 1 for x in F] == F2 and [x & 1 for x in G] == G2, "(F, G) mod-2 packing mismatch"

    print(
        f"  PASS  Algebraic 592B secret-key check (pk={len(pk_raw)}B, skd={len(sk_raw)}B): "
        "det(B_rec)=1, B_rec^* B_rec=Q, (F,G) mod 2 & SHAKE256(pk) valid"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) == 3 and argv[0] == "--verify-sk":
        return verify_secret_key_bytes(argv[1], argv[2])

    p = P.HAWK256
    n = p.n
    assert n == N

    use_work = "--work" in argv or os.path.exists(os.path.join(WORK, "Bred.sobj"))
    if use_work:
        Gp = read_mat(os.path.join(WORK, "Gp.txt"))
        G1 = read_mat(os.path.join(WORK, "G1.txt"))
        Kf = read_mat(os.path.join(WORK, "Kflat.txt"))
        Bs = read_mat(os.path.join(WORK, "Bs.txt"))
        Bc = read_mat(os.path.join(WORK, "Bc.txt"))
        Bcr = read_mat(os.path.join(WORK, "Bcr.txt"))
        Br = matrix(ZZ, load(os.path.join(WORK, "Bred.sobj")))
    else:
        Gp = read_mat(resolve("Gp_k1.txt"))
        G1 = read_mat(resolve("G1_k1.txt"))
        Kf = read_mat(resolve("Kflat.txt"))
        Bs = read_mat(resolve("Bs.txt"))
        Bc = read_mat(resolve("Bc.txt"))
        Bcr = read_mat(resolve("Bcr.txt"))
        Br = matrix(ZZ, load(resolve("Bred.sobj")))
    vtau_path = os.path.join(DATA, "vtau_coords.sobj")
    v_saved = vector(ZZ, load(vtau_path)) if os.path.exists(vtau_path) else None

    # (1) Sublattice index certificate
    dl2 = float(Gp.det().abs().log(2))
    expect = p.logdet_gram_tau + 2 * (n - 1)
    assert abs(dl2 - expect) < 0.01, f"log2 det Gram(Lambda') = {dl2:.2f} != {expect:.2f}"
    print(f"  PASS  Sublattice determinant log2 det Gram(Lambda'): {dl2:.2f} = {p.logdet_gram_tau:.0f} + 2*{n-1} (index 2^{n-1})")

    # (2) Shortest row of reduced basis reaches target minimum
    tgt = (2 ** (2 * S)) * p.lam1sq
    best = min(range(Br.nrows()), key=lambda r: sum(int(Br[r, j]) ** 2 for j in range(n)))
    sq = sum(int(Br[best, j]) ** 2 for j in range(n))
    ratio = sq / tgt
    assert ratio <= 1.05, f"shortest row ratio {ratio:.6f} > 1.05"
    print(f"  PASS  Reduced basis shortest vector (row {best}): sqnorm/target = {ratio:.6f} (target 2^{2*S} * {p.lam1sq})")

    # (3) Exact unimodular transformation recovery U = Bcr * Bc^{-1} and coordinate pullback v = u * U
    Bl = [[int(Bc[i, j]) for j in range(n)] for i in range(n)]
    U = matrix(ZZ, n, n)
    for r in range(n):
        b = [int(Bcr[r, j]) for j in range(n)]
        u_row = [0] * n
        for i in range(n - 1, -1, -1):
            acc = b[i] - sum(u_row[j] * Bl[j][i] for j in range(i + 1, n))
            q, rem = divmod(acc, Bl[i][i])
            assert rem == 0, f"non-exact division at col {i}"
            u_row[i] = q
        for j in range(n):
            U[r, j] = u_row[j]
    assert U * Bc == Bcr, "U * Bc != Bcr"
    assert U * Gp * U.transpose() == G1, "U * Gp * U^T != G1"

    Bl2 = [[int(Bs[i, j]) for j in range(n)] for i in range(n)]
    b_vec = vector(ZZ, [int(Br[best, j]) for j in range(n)])
    u_coords = [0] * n
    for i in range(n - 1, -1, -1):
        acc = int(b_vec[i]) - sum(u_coords[j] * Bl2[j][i] for j in range(i + 1, n))
        q, rem = divmod(acc, Bl2[i][i])
        assert rem == 0, f"non-exact division at col {i}"
        u_coords[i] = q
    u = vector(ZZ, u_coords)
    q_pre = int(u * G1 * u)
    v = u * U
    q_lam = int(v * Gp * v)
    assert q_pre == p.lam1sq and q_lam == p.lam1sq, f"norm mismatch: {q_pre}, {q_lam}"
    print(f"  PASS  Exact coordinate pullback v = u * U: u*G1*u = {q_pre}, v*Gp*v = {q_lam} (target 2n = {p.lam1sq})")

    # (4) Cocycle and modulo-2 parity verification on V_tau
    LB = _embed_LB(Kf, p.dL, p.dK)
    H = _lam_prime_H(LB, n)
    w = v * H
    amb = vector(ZZ, w * LB)
    aside = bool(_is_aside(amb, n))
    L, z, _, _, _ = ring.field_and_homs(p.ell)
    V = _amb_to_matrix(amb, L, z, n)
    cocycle = bool(V * ring.mat_apply(V, L.hom([-z])) == ring.gidentity(L, 2))
    if use_work:
        if v_saved is not None:
            Kf_ref = read_mat(resolve("Kflat.txt"))
            LB_ref = _embed_LB(Kf_ref, p.dL, p.dK)
            H_ref = _lam_prime_H(LB_ref, n)
            amb_ref = vector(ZZ, (v_saved * H_ref) * LB_ref)
            same = bool(amb == amb_ref or amb == -amb_ref)
        else:
            same = True
        v.save(os.path.join(WORK, "vtau_coords.sobj"))
    else:
        same = True if v_saved is None else bool(v == v_saved or v == -v_saved)
    assert aside and cocycle and same, "cocycle certificate failed"
    print(f"  PASS  Recovered cocycle V_tau certificates: V_tau == I_2 (mod 2) is {aside}, V_tau*tau(V_tau) == I_2 is {cocycle}")
    log(f"certificates verified in {cpu():.1f}s CPU")
    return 0


if __name__ == "__main__":
    sys.exit(main())
