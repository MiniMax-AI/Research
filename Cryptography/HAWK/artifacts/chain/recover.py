#!/usr/bin/env python3
"""Stages 7-8: V_tau -> B_rec via tau-tower descent and 592-byte secret-key encoding.

Loads the coordinate vector v of V_tau in the Lambda' basis (produced by
`certificates.py` in Stage 6, or loaded from data/vtau_coords.sobj), reconstructs
the 2x2 matrix V_tau, runs the cyclotomic tau-tower descent
(`hawk_recover.tower.recover_basis`, fetched via `fetch_third_party.sh` from
github.com/anthropics/cryptography-research-demo @ fa01c398), and encodes B_rec
into the 592-byte reference secret-key format (`recovered_skd.bin`) for C
sign/verify checking via `kat/hawk_shim.c` (invoked automatically by `./run.sh`).

Run:  python3 chain/recover.py [vtau_coords.sobj]
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from sage.all import ZZ, matrix, vector, load

from common import DATA, WORK, cpu, ensure_work, log, read_mat, resolve
from hawk_recover import decode, params as P, ring, tower, encode
from hawk_recover.cocycle import _embed_LB, _gram_tau
from hawk_recover.types import Intertwiner
from hawk_recover.tower import _amb_to_matrix, _is_aside


def lam_prime_H(LB, n):
    from sage.all import GF
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


def main(argv):
    p = P.HAWK256
    n = p.n
    out = ensure_work()
    pk = decode.decode(open(os.path.join(DATA, "pk.bin"), "rb").read(), p)

    use_work = os.path.exists(os.path.join(WORK, "vtau_coords.sobj")) and os.path.exists(
        os.path.join(WORK, "Kflat.txt")
    )
    vfile = argv[0] if argv else (
        os.path.join(WORK, "vtau_coords.sobj")
        if use_work
        else os.path.join(DATA, "vtau_coords.sobj")
    )
    v = vector(ZZ, load(vfile))
    log(f"loaded V_tau coordinates from {os.path.basename(vfile)} "
        f"(len {len(v)}, max|u| {max(abs(x) for x in v)})")

    Kf = read_mat(os.path.join(WORK, "Kflat.txt") if use_work and not argv else resolve("Kflat.txt"))
    LB = _embed_LB(Kf, p.dL, p.dK)
    H = lam_prime_H(LB, n)
    amb = vector(ZZ, (v * H) * LB)

    log(f"is_aside (== I mod 2): {_is_aside(amb, n)}")
    L, z, c, _, _ = ring.field_and_homs(p.ell)
    V = _amb_to_matrix(amb, L, z, n)
    ok = (V * ring.mat_apply(V, L.hom([-z])) == ring.gidentity(L, 2))
    log(f"V_tau * tau(V_tau) == I : {ok}")
    if not ok:
        print("RESULT: FAIL")
        return 1

    itw = Intertwiner(params=p, wt=None, V=V, la_rank=n // 2)
    itw.check(L)
    log("V_tau ACCEPTED as the a-side cocycle")

    t = time.time()
    c0 = cpu()
    rb = tower.recover_basis(itw, pk, verbose=False)
    log(f"tower: wall {time.time()-t:.1f}s cpu {cpu()-c0:.1f}s  "
        f"B*B==Q: {rb.eq_G}   beta_trace={rb.beta_trace}")
    if not rb.eq_G:
        print("RESULT: FAIL")
        return 1

    sk = encode.encode(rb, pk)
    dst = os.path.join(out, "recovered_skd.bin")
    with open(dst, "wb") as f:
        f.write(sk.core + sk.hpub)
    log(f"wrote {dst} ({len(sk.core) + len(sk.hpub)} bytes)")

    print()
    print(f"  total CPU-core-seconds = {cpu():.1f}")
    print("RESULT: SUCCESS")
    print()
    print("Now verify with the reference implementation's own signer/verifier:")
    print(f"  bash ./kat/build_shim.sh && ./kat/hawk_shim signverify data/pk.bin {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
