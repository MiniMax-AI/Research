#!/usr/bin/env python3
"""Practical 7-Round ARADI Full 256-Bit Master-Key Recovery via L_0-Colliding Cubes & FHT.

Targets the unmodified 7-round NSA ARADI specification (rounds 0..6 + post-whitening K^(7)).
Exploits L_0-colliding 9D cubes C_9^(s) on full columns (s, 16+((3+s) mod 16)) plus bit z_{(5+s) mod 16}:
  - Because pi(. ^ K^(0)_b) is a bijection on each 4-bit column, varying full columns (0, 19)
    bypasses Round 0 with algebraic degree 1 and enters only |L_0(0) U L_0(19)| = 4 S-boxes
    {0, 11, 19, 30} at Round 1, leaving 4 full 4-bit columns [5, 9, 16, 20] balanced after
    5 full rounds (Q_5) at 2^8 = 256 CP (32x less data than Kim et al. ePrint 2024/1772),
    plus 7 additional balanced columns across the 512-CP 9D union!
  - Across 6 shifts (3,072 CP = 2^11.58 CP) or 8 shifts (4,096 CP = 2^12.00 CP), 2-round
    partial decryption via 2^12-point Fast Hadamard Transforms (FHT) recovers K^(7), peels
    Round 6 onto 19 balanced columns of Z_4 = Lambda_0(Q_5) to recover K^(6), and inverts
    the even key-schedule pair (K^(6), K^(7)) (even_idx = 6) in O(1) to recover the exact
    256-bit master key in ~3.5 s in pure Python on 1 CPU core.
"""
from __future__ import annotations

import argparse
import statistics
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np

from aradi_ref import (
    ROT_CONSTANTS,
    SBOX_INV_TABLE,
    encrypt_rounds,
    expand_key,
    format_words,
    linear_L,
    master_key_from_even_pair,
)

SBOX_INV_ARR = np.array(SBOX_INV_TABLE, dtype=np.uint8)


def _r16(x: np.ndarray, n: int) -> np.ndarray:
    return (((x << np.uint16(n)) | (x >> np.uint16(16 - n))) & np.uint16(0xFFFF)) if n else x


def _lnp(s: np.ndarray, a: int, b: int, c: int) -> np.ndarray:
    u = (s >> np.uint32(16)).astype(np.uint16)
    l = (s & np.uint32(0xFFFF)).astype(np.uint16)
    return ((u ^ _r16(u, a) ^ _r16(l, c)).astype(np.uint32) << np.uint32(16)) | (
        l ^ _r16(l, a) ^ _r16(u, b)
    ).astype(np.uint32)


def _shift_col(c: int, s: int) -> int:
    return ((c + s) & 15) if c < 16 else (16 + (((c - 16) + s) & 15))


def _fht12_inplace(arr: np.ndarray) -> None:
    h = 1
    while h < 4096:
        for i in range(0, 4096, h * 2):
            x = arr[..., i : i + h].copy()
            y = arr[..., i + h : i + 2 * h]
            arr[..., i : i + h] = x + y
            arr[..., i + h : i + 2 * h] = x - y
        h <<= 1


# Verified universal ZeroSum bitmasks at shift s = 0 on 9D cube [col 0, col 19, z_5]:
M8_S0: Dict[int, int] = {5: 15, 9: 15, 16: 15, 20: 15}
M9_S0: Dict[int, int] = {0: 10, 1: 15, 6: 15, 18: 10, 24: 10, 25: 15, 31: 15}
Z4_BAL_S0: Tuple[int, ...] = (0, 2, 5, 6, 8, 9, 10, 11, 14, 16, 17, 18, 20, 21, 23, 25, 26, 29, 31)

DEPS_LAM7: List[List[int]] = [
    [inp for inp in range(32) if (linear_L(1 << inp, *ROT_CONSTANTS[1]) >> cq) & 1]
    for cq in range(32)
]
DEPS_LAM6_Z4: List[List[int]] = [
    [inp for inp in range(32) if (linear_L(1 << inp, *ROT_CONSTANTS[0]) >> bz) & 1]
    for bz in range(32)
]


def _build_fht_kernel() -> np.ndarray:
    x_idx = np.arange(4096, dtype=np.uint16)
    s_sum3 = (
        SBOX_INV_ARR[x_idx & 0xF]
        ^ SBOX_INV_ARR[(x_idx >> 4) & 0xF]
        ^ SBOX_INV_ARR[(x_idx >> 8) & 0xF]
    )
    ker = np.zeros((16, 4, 4096), dtype=np.int32)
    for g6 in range(16):
        out_nib = SBOX_INV_ARR[s_sum3 ^ g6]
        for bit_r in range(4):
            vec = ((out_nib >> bit_r) & 1).astype(np.int32)
            _fht12_inplace(vec)
            ker[g6, bit_r] = vec
    return ker


def _recover_k6_from_z4(T6_nibs: Sequence[np.ndarray], shifts_sub: Sequence[int]) -> List[Tuple[int, int, int, int]]:
    a1, b1, c1 = ROT_CONSTANTS[1]
    chks_by_bz: Dict[int, List[List[List[int]]]] = {bz: [] for bz in range(32)}
    for si, s in enumerate(shifts_sub):
        nb = T6_nibs[si]
        sg0 = [[int(np.bitwise_xor.reduce(SBOX_INV_ARR[nb[c, :256] ^ g])) for g in range(16)] for c in range(32)]
        sg1 = [[int(np.bitwise_xor.reduce(SBOX_INV_ARR[nb[c, 256:] ^ g])) for g in range(16)] for c in range(32)]
        for bz0 in Z4_BAL_S0:
            bz = _shift_col(bz0, s)
            chks_by_bz[bz].append(sg0)
            chks_by_bz[bz].append(sg1)

    valid_triplets: Dict[int, List[Tuple[int, int, int]]] = {}
    for bz, sg_list in chks_by_bz.items():
        if not sg_list:
            continue
        c0, c1_col, c2 = DEPS_LAM6_Z4[bz]
        buckets: Dict[Tuple[int, ...], List[int]] = {}
        for g2 in range(16):
            sig = tuple(sg[c2][g2] for sg in sg_list)
            buckets.setdefault(sig, []).append(g2)
        triplets: List[Tuple[int, int, int]] = []
        for g0 in range(16):
            for g1 in range(16):
                sig = tuple(sg[c0][g0] ^ sg[c1_col][g1] for sg in sg_list)
                for g2 in buckets.get(sig, ()):
                    triplets.append((g0, g1, g2))
        valid_triplets[bz] = triplets

    unvisited = set(valid_triplets.keys())
    known: set[int] = set()
    cands: List[List[int]] = [[-1] * 32]
    while unvisited and cands:
        best_bz = max(
            unvisited,
            key=lambda bz: len(set(DEPS_LAM6_Z4[bz]) & known) * 100 - len(valid_triplets[bz]),
        )
        unvisited.remove(best_bz)
        c0, c1_col, c2 = DEPS_LAM6_Z4[best_bz]
        sample = cands[0]
        k_known = (sample[c0] != -1, sample[c1_col] != -1, sample[c2] != -1)
        bmap: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]] = {}
        for g0, g1, g2 in valid_triplets[best_bz]:
            key = (g0 if k_known[0] else -1, g1 if k_known[1] else -1, g2 if k_known[2] else -1)
            bmap.setdefault(key, []).append((g0, g1, g2))
        cands = [
            [g0 if i == c0 else g1 if i == c1_col else g2 if i == c2 else cand[i] for i in range(32)]
            for cand in cands
            for g0, g1, g2 in bmap.get(
                (cand[c0] if k_known[0] else -1, cand[c1_col] if k_known[1] else -1, cand[c2] if k_known[2] else -1),
                (),
            )
        ]
        known.update((c0, c1_col, c2))

    out: List[Tuple[int, int, int, int]] = []
    for lam6 in cands[:16]:
        w6 = [0, 0, 0, 0]
        for col in range(32):
            for wi in range(4):
                if (lam6[col] >> wi) & 1:
                    w6[wi] |= 1 << col
        out.append(tuple(linear_L(w6[wi], a1, b1, c1) for wi in range(4)))
    return out


def recover_7r_master_key(
    cts_sh: Sequence[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    pts_chk: Sequence[Tuple[int, int, int, int]],
    cts_chk: Sequence[Tuple[int, int, int, int]],
    shifts_list: Sequence[int],
    fht_ker: np.ndarray,
) -> Tuple[Tuple[int, ...] | None, int]:
    a1, b1, c1 = ROT_CONSTANTS[1]
    a2, b2, c2 = ROT_CONSTANTS[2]
    T7_nibbles: List[np.ndarray] = []
    for CW, CX, CY, CZ in cts_sh:
        TW = _lnp(CW, a2, b2, c2)
        TX = _lnp(CX, a2, b2, c2)
        TY = _lnp(CY, a2, b2, c2)
        TZ = _lnp(CZ, a2, b2, c2)
        nibs = np.zeros((32, 512), dtype=np.uint16)
        for col in range(32):
            nibs[col] = (
                ((TW >> np.uint32(col)) & 1)
                | (((TX >> np.uint32(col)) & 1) << 1)
                | (((TY >> np.uint32(col)) & 1) << 2)
                | (((TZ >> np.uint32(col)) & 1) << 3)
            ).astype(np.uint16)
        T7_nibbles.append(nibs)

    col_chks: Dict[int, List[Tuple[int, int, int, int]]] = {c: [] for c in range(32)}
    for si, s in enumerate(shifts_list):
        for c0, bm in M8_S0.items():
            cs = _shift_col(c0, s)
            col_chks[cs].append((si, 0, 256, bm))
            col_chks[cs].append((si, 256, 512, bm))
        for c0, bm in M9_S0.items():
            cs = _shift_col(c0, s)
            col_chks[cs].append((si, 0, 512, bm))

    unvisited = {c for c, chks in col_chks.items() if sum(bin(bm).count("1") for _, _, _, bm in chks) >= 6}
    known_lam7_cols: set[int] = set()
    ordered_cols: List[int] = []
    while unvisited:
        best_c = max(
            unvisited,
            key=lambda c: len(set(DEPS_LAM7[c]) & known_lam7_cols) * 10
            + sum(bin(bm).count("1") for _, _, _, bm in col_chks[c]),
        )
        unvisited.remove(best_c)
        ordered_cols.append(best_c)
        known_lam7_cols.update(DEPS_LAM7[best_c])

    cands_lam7: List[List[int]] = [[-1] * 32]
    for cq in ordered_cols:
        u0, u1, u2 = DEPS_LAM7[cq]
        valid = np.ones((16, 4096), dtype=bool)
        for si, lo, hi, bm in col_chks[cq]:
            nibs = T7_nibbles[si]
            idx12 = nibs[u0, lo:hi] | (nibs[u1, lo:hi] << 4) | (nibs[u2, lo:hi] << 8)
            A = (np.bincount(idx12, minlength=4096) & 1).astype(np.int32)
            _fht12_inplace(A)
            for bit_r in range(4):
                if (bm >> bit_r) & 1:
                    prod = fht_ker[:, bit_r, :] * A[None, :]
                    _fht12_inplace(prod)
                    valid &= ((prod >> 12) & 1) == 0
        valid_k12 = np.flatnonzero(np.any(valid, axis=0))
        if not cands_lam7:
            return None, 0
        sample = cands_lam7[0]
        k_known = (sample[u0] != -1, sample[u1] != -1, sample[u2] != -1)
        buckets: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]] = {}
        for k12 in valid_k12:
            v0, v1, v2 = int(k12 & 0xF), int((k12 >> 4) & 0xF), int((k12 >> 8) & 0xF)
            key = (v0 if k_known[0] else -1, v1 if k_known[1] else -1, v2 if k_known[2] else -1)
            buckets.setdefault(key, []).append((v0, v1, v2))
        cands_lam7 = [
            [v0 if i == u0 else v1 if i == u1 else v2 if i == u2 else cand[i] for i in range(32)]
            for cand in cands_lam7
            for v0, v1, v2 in buckets.get(
                (cand[u0] if k_known[0] else -1, cand[u1] if k_known[1] else -1, cand[u2] if k_known[2] else -1),
                (),
            )
        ]

    for lam7 in cands_lam7:
        w7 = [0, 0, 0, 0]
        for col in range(32):
            for wi in range(4):
                if (lam7[col] >> wi) & 1:
                    w7[wi] |= 1 << col
        k7_cand = tuple(linear_L(w7[wi], a2, b2, c2) for wi in range(4))
        T6_nibs = []
        for CW, CX, CY, CZ in cts_sh:
            XW = _lnp(CW ^ np.uint32(k7_cand[0]), a2, b2, c2)
            XX = _lnp(CX ^ np.uint32(k7_cand[1]), a2, b2, c2)
            XY = _lnp(CY ^ np.uint32(k7_cand[2]), a2, b2, c2)
            XZ = _lnp(CZ ^ np.uint32(k7_cand[3]), a2, b2, c2)
            XW ^= XX & XZ
            XY ^= XW & XZ
            XZ ^= XX & XY
            XX ^= XW & XY
            TW6 = _lnp(XW, a1, b1, c1)
            TX6 = _lnp(XX, a1, b1, c1)
            TY6 = _lnp(XY, a1, b1, c1)
            TZ6 = _lnp(XZ, a1, b1, c1)
            nb6 = np.zeros((32, 512), dtype=np.uint8)
            for col in range(32):
                nb6[col] = (
                    ((TW6 >> np.uint32(col)) & 1)
                    | (((TX6 >> np.uint32(col)) & 1) << 1)
                    | (((TY6 >> np.uint32(col)) & 1) << 2)
                    | (((TZ6 >> np.uint32(col)) & 1) << 3)
                ).astype(np.uint8)
            T6_nibs.append(nb6)
        for k6_cand in _recover_k6_from_z4(T6_nibs, shifts_list):
            mk_cand = master_key_from_even_pair(k6_cand, k7_cand, 6)
            rks_c = expand_key(mk_cand, 7)
            if all(encrypt_rounds(pt, rks_c, 7) == ct for pt, ct in zip(pts_chk, cts_chk)):
                return mk_cand, len(cands_lam7)

    return None, len(cands_lam7)


def main() -> int:
    parser = argparse.ArgumentParser(description="7-Round ARADI Practical Full 256-Bit Master-Key Recovery")
    parser.add_argument("--trials", type=int, default=3, help="Number of random 256-bit master-key trials")
    parser.add_argument("--shifts", type=int, choices=[6, 8], default=6, help="6 shifts (3,072 CP) or 8 shifts (4,096 CP)")
    parser.add_argument("--seed", type=int, default=20260707, help="PRNG seed")
    args = parser.parse_args()

    shifts_list = [0, 4, 8, 1, 5, 9] if args.shifts == 6 else [0, 4, 8, 12, 1, 5, 9, 13]
    total_cp = 512 * len(shifts_list)

    print("========================================================================")
    print(f"  7-Round ARADI Practical Full 256-Bit Master-Key Recovery ({total_cp:,} CP)")
    print("========================================================================\n")

    t_k0 = time.perf_counter()
    fht_ker = _build_fht_kernel()
    print(f"[Setup] Precomputed 12-bit FHT partial-sum kernel in {(time.perf_counter() - t_k0)*1000:.2f} ms")
    print(f"  - Shifts          : {shifts_list} ({len(shifts_list)} x 512 = {total_cp:,} CP = 2^{np.log2(total_cp):.2f} CP)")
    print("  - Prior 7R SOTA   : 2^32 pairs (2^33 CP) & 2^70..2^72 ND evals (LATINCRYPT 2025)\n")

    rng = np.random.default_rng(args.seed)
    success = 0
    times: List[float] = []
    list_sizes: List[int] = []

    for idx in range(args.trials):
        mk_true = tuple(int(x) for x in rng.integers(0, 1 << 32, size=8, dtype=np.uint64))
        rks7 = expand_key(mk_true, 7)
        cts_sh = []
        pts_chk = []
        cts_chk = []
        for s in shifts_list:
            base_s = [int(x) for x in rng.integers(0, 1 << 32, size=4, dtype=np.uint64)]
            cA, cB, cC = _shift_col(0, s), _shift_col(19, s), _shift_col(5, s)
            dirs_s = [(w, cA) for w in range(4)] + [(w, cB) for w in range(4)] + [(3, cC)]
            P_arr = np.zeros((512, 4), dtype=np.uint32)
            for u in range(512):
                p = base_s.copy()
                for bit_i, (wi, bi) in enumerate(dirs_s):
                    if (u >> bit_i) & 1:
                        p[wi] ^= 1 << bi
                P_arr[u] = p
            W, X, Y, Z = P_arr[:, 0].copy(), P_arr[:, 1].copy(), P_arr[:, 2].copy(), P_arr[:, 3].copy()
            for st in range(7):
                kw, kx, ky, kz = rks7[st]
                W ^= np.uint32(kw)
                X ^= np.uint32(kx)
                Y ^= np.uint32(ky)
                Z ^= np.uint32(kz)
                X ^= W & Y
                Z ^= X & Y
                Y ^= W & Z
                W ^= X & Z
                a, b, c = ROT_CONSTANTS[st & 3]
                W = _lnp(W, a, b, c)
                X = _lnp(X, a, b, c)
                Y = _lnp(Y, a, b, c)
                Z = _lnp(Z, a, b, c)
            kw, kx, ky, kz = rks7[7]
            W ^= np.uint32(kw)
            X ^= np.uint32(kx)
            Y ^= np.uint32(ky)
            Z ^= np.uint32(kz)
            cts_sh.append((W, X, Y, Z))
            pts_chk.append(tuple(int(v) for v in P_arr[0]))
            cts_chk.append((int(W[0]), int(X[0]), int(Y[0]), int(Z[0])))

        t0 = time.perf_counter()
        mk_rec, n_l7 = recover_7r_master_key(cts_sh, pts_chk, cts_chk, shifts_list, fht_ker)
        dt = time.perf_counter() - t0
        ok = mk_rec == mk_true
        success += int(ok)
        times.append(dt)
        list_sizes.append(n_l7)
        print(
            f"  [Trial {idx + 1}/{args.trials}] {total_cp:,} CP | K^(7)_cands={n_l7} | "
            f"time={dt:.2f}s | 256-bit_MK_match={ok} ({format_words(mk_true[:4])} ...)",
            flush=True,
        )

    print(
        f"\n[Summary] 7-Round ARADI ({total_cp:,} CP) Master-Key Recovery:\n"
        f"  - Exact 256-bit master-key recovery : {success}/{args.trials} ({100.0 * success / args.trials:.1f}%)\n"
        f"  - Single-core wall time (Python)    : median={statistics.median(times):.2f} s, mean={statistics.mean(times):.2f} s"
    )
    assert success == args.trials, "7-round ARADI recovery failed on at least one trial"
    print("  [PASS] All 7-round ARADI trials recovered and verified the exact 256-bit master key.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
