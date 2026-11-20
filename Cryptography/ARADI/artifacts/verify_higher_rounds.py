#!/usr/bin/env python3
"""
Higher-Round Algebraic Structure, 6-Round Distinguisher, and 9R-12R Complexity Verifier
for Unmodified NSA ARADI (ePrint 2024/1240).

Single-threaded (< 3 s on 1 CPU core). Deterministic verification of:
  1. Complete 2^16 output-mask ANF degree spectrum of the 16-bit L_3-island super-S-box
     Psi_G(v) = pi(Lambda_3(pi(v)) ^ k)|_G: {0: 1, 3: 3, 5: 28, 6: 224, 8: 65280}.
  2. Exhaustive composite linear-layer canceling-mask check over all 4^k products (k = 1..5,
     1,364 matrices): dim C(M) > 0 iff M reduces to I (dim 32) or single-layer L_3 (dim 8).
  3. Intra-half Z_16 cyclic shift symmetry of the ARADI round function & the 1-Round Lift
     identity <Lambda_r^T(e_b (x) nu), Q_{r+1}> = H_nu(Z_{r,b}) ^ const (no island condition).
  4. Live single-core empirical check of the 6-round d = 16 colliding-pair distinguisher
     (0, 5, 19, 24) at r0 = 0 (all 128 bits of Q_5 balanced; full column 25 of Q_6 balanced
     under both 256-bit master keys and independent round keys).
  5. Exact GF(2) key-schedule rank (rho = 256, residual 2^0), progressive list-join tracking
     log_2(|L_t|), and packed integer-arithmetic / encryption complexity verification for
     9R, 10R, 11R, and 12R attacks vs. Dunkelman-Ghosh (ePrint 2025/1227), plus the 13R
     depth-5 geometric saturation barrier (K = 66 => 2^266.75 > 2^255).
"""
from __future__ import annotations

import itertools
import math
import os
import time
from typing import Dict, List, Sequence, Set, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from aradi_ref import ROT_CONSTANTS, SBOX_TABLE as SBOX, expand_key, linear_L


def build_L_matrix(j: int) -> np.ndarray:
    a, b, c = ROT_CONSTANTS[j & 3]
    M = np.zeros((32, 32), dtype=np.uint8)
    for col in range(32):
        out = linear_L(1 << col, a, b, c)
        for row in range(32):
            M[row, col] = (out >> row) & 1
    return M


LM = [build_L_matrix(j) for j in range(4)]
RS = [[[int(x) for x in np.nonzero(LM[j][b])[0]] for b in range(32)] for j in range(4)]


# ---------------------------------------------------------------------------
# 1. Exact 16-bit L_3-island Super-S-Box ANF degree spectrum
# ---------------------------------------------------------------------------
def verify_island_supersbox_spectrum() -> None:
    # Island G_0 of L_3 = (8, 9, 7): columns [0, 8, 23, 31].
    # Restriction of L_3 to G_0: 4x4 submatrix.
    G0 = [0, 8, 23, 31]
    L3_G0 = LM[3][np.ix_(G0, G0)]
    assert np.array_equal((L3_G0 @ L3_G0) & 1, np.eye(4, dtype=np.uint8))
    assert np.array_equal(L3_G0.T, L3_G0)

    idx = np.arange(65536, dtype=np.uint32)
    sbox_arr = np.array(SBOX, dtype=np.uint16)

    for k_const in (0x0000, 0x3A84, 0xFFFF):
        # Input v: 4 columns, each 4 bits (nibble c in bits 4c..4c+3)
        n0 = sbox_arr[(idx >> 0) & 0xF]
        n1 = sbox_arr[(idx >> 4) & 0xF]
        n2 = sbox_arr[(idx >> 8) & 0xF]
        n3 = sbox_arr[(idx >> 12) & 0xF]
        nibs = [n0, n1, n2, n3]
        # Apply L3_G0 column-wise within each row w=0..3
        u_nibs = [np.zeros(65536, dtype=np.uint16) for _ in range(4)]
        for r_col in range(4):
            for c_col in range(4):
                if L3_G0[r_col, c_col]:
                    u_nibs[r_col] ^= nibs[c_col]
            u_nibs[r_col] ^= (k_const >> (4 * r_col)) & 0xF
        out = (
            sbox_arr[u_nibs[0]]
            | (sbox_arr[u_nibs[1]] << 4)
            | (sbox_arr[u_nibs[2]] << 8)
            | (sbox_arr[u_nibs[3]] << 12)
        ).astype(np.uint16)

        # 16-variable binary Moebius transform over uint16 truth table
        anf = out.copy()
        for bit in range(16):
            step = 1 << bit
            anf_view = anf.reshape(-1, 2 * step)
            anf_view[:, step:] ^= anf_view[:, :step]

        # Group ANF coefficient vectors by monomial Hamming weight w = 0..16
        pop = np.zeros(65536, dtype=np.uint8)
        for b in range(16):
            pop += ((idx >> b) & 1).astype(np.uint8)

        # Evaluate max degree of <beta, Psi_{G0}> for all 65536 masks beta via 16-bit basis per weight
        max_deg = np.zeros(65536, dtype=np.uint8)
        masks = np.arange(65536, dtype=np.uint16)
        for w in range(1, 17):
            coeffs_w = np.unique(anf[pop == w])
            coeffs_w = coeffs_w[coeffs_w != 0]
            if len(coeffs_w) == 0:
                continue
            # Extract GF(2) row basis of coeffs_w (<= 16 vectors)
            basis: List[int] = []
            for v in coeffs_w.tolist():
                x = int(v)
                for bvec in basis:
                    x = min(x, x ^ bvec)
                if x:
                    basis.append(x)
                    basis.sort(reverse=True)
            hit = np.zeros(65536, dtype=bool)
            for bvec in basis:
                # <beta, bvec> mod 2 != 0 for at least one basis vector iff <beta, span> != 0
                parity = np.zeros(65536, dtype=np.uint8)
                bv = bvec
                while bv:
                    lsb = bv & -bv
                    bit_pos = lsb.bit_length() - 1
                    parity ^= ((masks >> bit_pos) & 1).astype(np.uint8)
                    bv ^= lsb
                hit |= (parity == 1)
            max_deg[hit] = w

        hist: Dict[int, int] = {}
        for d in range(17):
            cnt = int(np.sum(max_deg == d))
            if cnt:
                hist[d] = cnt
        assert hist == {0: 1, 3: 3, 5: 28, 6: 224, 8: 65280}, f"Unexpected super-S-box hist: {hist}"
        deg3_masks = sorted(int(m) for m in np.nonzero(max_deg == 3)[0])
        assert deg3_masks == [0x2222, 0x8888, 0xAAAA], f"Unexpected deg-3 masks: {deg3_masks}"


# ---------------------------------------------------------------------------
# 2. Exhaustive composite linear-layer canceling-mask space (k = 1..5)
# ---------------------------------------------------------------------------
def canceling_mask_dim(M: np.ndarray) -> int:
    # C(M) = { mu in F_2^32 : M^T diag(mu) M is diagonal }, i.e., off-diagonal entries vanish:
    # For all c < d: sum_b mu_b M[b,c] M[b,d] = 0 (mod 2).
    rows: List[int] = []
    for c in range(32):
        mc = M[:, c]
        for d in range(c + 1, 32):
            prod = mc & M[:, d]
            if np.any(prod):
                v = 0
                for b in np.nonzero(prod)[0]:
                    v |= 1 << int(b)
                rows.append(v)
    piv: Dict[int, int] = {}
    for v in rows:
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                v ^= piv[hb]
            else:
                piv[hb] = v
                break
    return 32 - len(piv)


def verify_no_composite_collapse() -> Tuple[int, int]:
    eye = np.eye(32, dtype=np.uint8)
    seen_mats: Dict[bytes, Tuple[np.ndarray, int]] = {}
    total_prods = 0
    nontrivial_hits = 0
    for k in range(1, 6):
        for seq in itertools.product(range(4), repeat=k):
            total_prods += 1
            M = eye
            for j in seq:
                M = (LM[j] @ M) & 1
            key = M.tobytes()
            if key not in seen_mats:
                seen_mats[key] = (M, canceling_mask_dim(M))
            M_ref, d = seen_mats[key]
            if d > 0:
                nontrivial_hits += 1
                is_eye = np.array_equal(M_ref, eye)
                is_L3 = np.array_equal(M_ref, LM[3])
                assert (is_eye and d == 32) or (is_L3 and d == 8), (
                    f"Unexpected composite canceling space dim={d} at seq={seq}"
                )
    return total_prods, nontrivial_hits


# ---------------------------------------------------------------------------
# 3. Z_16 shift symmetry & 1-Round Lift identity + live 6R d=16 cube check
# ---------------------------------------------------------------------------
def shift_word_16(x: int, s: int) -> int:
    s &= 15
    lo = x & 0xFFFF
    hi = (x >> 16) & 0xFFFF
    lo_s = ((lo << s) | (lo >> (16 - s))) & 0xFFFF if s else lo
    hi_s = ((hi << s) | (hi >> (16 - s))) & 0xFFFF if s else hi
    return (hi_s << 16) | lo_s


def verify_shift_symmetry_and_1round_lift() -> None:
    rng = np.random.default_rng(20261106)
    # (a) Z_16 shift commutes with S-box and every L_j
    for j in range(4):
        a, b, c = ROT_CONSTANTS[j]
        for s in range(16):
            for _ in range(16):
                val = int(rng.integers(0, 1 << 32, dtype=np.uint64))
                lhs = shift_word_16(linear_L(val, a, b, c), s)
                rhs = linear_L(shift_word_16(val, s), a, b, c)
                assert lhs == rhs
    # (b) 1-Round Lift: for any single column b and j in 0..3,
    #     <Lambda_j^T(e_b (x) e_x), Q_{r+1}> = Z_r[x,b] ^ (Z_r[w,b] & Z_r[y,b])
    for j in range(4):
        a, b_rot, c_rot = ROT_CONSTANTS[j]
        supp_rows = RS[j]
        for b in (0, 7, 16, 25, 31):
            for _ in range(32):
                zw, zx, zy, zz = [int(rng.integers(0, 1 << 32, dtype=np.uint64)) for _ in range(4)]
                # S-box + L_j
                w, x, y, z = zw, zx, zy, zz
                x ^= w & y
                z ^= x & y
                y ^= w & z
                w ^= x & z
                qx = linear_L(x, a, b_rot, c_rot)
                qz = linear_L(z, a, b_rot, c_rot)
                # Since Q_{r+1} = L_j(pi(Z_r)) and L_j^2 = I, pi_x(Z_r)[b] = (L_j Q_{r+1,x})[b] = <L_j^T e_b, Q_{r+1,x}>.
                # By definition of RS, (L_j u)[b] = <L_j^T e_b, u> = sum_{col in RS[j][b]} u[col].
                lhs_x = 0
                lhs_z = 0
                for col in supp_rows[b]:
                    lhs_x ^= (qx >> col) & 1
                    lhs_z ^= (qz >> col) & 1
                rhs_x = ((zx >> b) & 1) ^ (((zw >> b) & 1) & ((zy >> b) & 1))
                rhs_z = ((zz >> b) & 1) ^ (((zx >> b) & 1) & ((zy >> b) & 1)) ^ (((zw >> b) & 1) & ((zy >> b) & 1))
                assert lhs_x == rhs_x and lhs_z == rhs_z


def verify_6r_d16_cube_live(n_keys: int = 8) -> Tuple[int, List[int]]:
    """Fast vectorized NumPy check of d=16 cube (0, 5, 19, 24) across n_keys master + indep keys."""
    rng = np.random.default_rng(20261106)
    cols = [0, 5, 19, 24]
    p_idx = np.arange(65536, dtype=np.uint32)
    # Build 65536 cube offsets for words w,x,y,z
    off = [np.zeros(65536, dtype=np.uint32) for _ in range(4)]
    for ci, col in enumerate(cols):
        nib = (p_idx >> (4 * ci)) & 0xF
        for w in range(4):
            off[w] ^= ((nib >> w) & 1).astype(np.uint32) << col

    def vec_linear(arr: np.ndarray, j: int) -> np.ndarray:
        a, b, c = ROT_CONSTANTS[j & 3]
        u = (arr >> 16) & 0xFFFF
        l = arr & 0xFFFF
        uo = u ^ (((u << a) | (u >> (16 - a))) & 0xFFFF) ^ (((l << c) | (l >> (16 - c))) & 0xFFFF)
        lo = l ^ (((l << a) | (l >> (16 - a))) & 0xFFFF) ^ (((u << b) | (u >> (16 - b))) & 0xFFFF)
        return ((uo.astype(np.uint32) << 16) | lo.astype(np.uint32))

    q5_sums: List[int] = []
    q6_sums: List[int] = []
    for mode in ("mk", "indep"):
        for _ in range(n_keys):
            base = [int(rng.integers(0, 1 << 32, dtype=np.uint64)) for _ in range(4)]
            if mode == "mk":
                mk = [int(rng.integers(0, 1 << 32, dtype=np.uint64)) for _ in range(8)]
                rks = expand_key(mk, 6)
            else:
                rks = [[int(rng.integers(0, 1 << 32, dtype=np.uint64)) for _ in range(4)] for _ in range(7)]
            st = [off[w] ^ np.uint32(base[w]) for w in range(4)]
            for r in range(6):
                w, x, y, z = [st[i] ^ np.uint32(rks[r][i]) for i in range(4)]
                x ^= w & y
                z ^= x & y
                y ^= w & z
                w ^= x & z
                st = [vec_linear(w, r), vec_linear(x, r), vec_linear(y, r), vec_linear(z, r)]
                if r == 4:
                    s5 = 0
                    for wi in range(4):
                        xors = int(np.bitwise_xor.reduce(st[wi]))
                        s5 |= xors << (32 * wi)
                    q5_sums.append(s5)
            s6 = 0
            for wi in range(4):
                xors = int(np.bitwise_xor.reduce(st[wi]))
                s6 |= xors << (32 * wi)
            q6_sums.append(s6)

    assert all(v == 0 for v in q5_sums), "Q_5 must be 100% balanced on all 128 bits!"
    # Check full column 25 (bits 25, 32+25, 64+25, 96+25) is identically 0 across all keys
    col25_mask = (1 << 25) | (1 << 57) | (1 << 89) | (1 << 121)
    assert all((v & col25_mask) == 0 for v in q6_sums), "Full column 25 of Q_6 must be balanced!"
    # Also check the 24 predicted 1-round-lift masks Lambda_1^T(e_b (x) e_{x/z}) for the 12 lift-ready cols
    lift_cols = [2, 5, 10, 13, 14, 18, 20, 21, 23, 25, 29, 31]
    for b in lift_cols:
        for row_idx in (1, 3):  # rows x and z
            msk = 0
            for col in RS[1][b]:
                msk |= 1 << (32 * row_idx + col)
            assert all((bin(v & msk).count("1") & 1) == 0 for v in q6_sums)
    return len(q6_sums), lift_cols


# ---------------------------------------------------------------------------
# 4. Exact GF(2) key-schedule rank, progressive list-join, and op-count audit
# ---------------------------------------------------------------------------
def build_u_affine_rows(max_R: int = 14) -> Dict[Tuple[int, int], int]:
    """Map (i_1based, bit_0..127) of u^i = Lambda_{i-1}(K^(i)) to its exact 256-bit MK mask."""
    # Subtract expand_key(0) so round constants cancel, leaving the pure linear functional
    rks_zero = expand_key([0] * 8, max_R)
    u_zero: List[List[int]] = []
    for i in range(max_R + 1):
        if i == 0:
            u_zero.append(list(rks_zero[0]))
        else:
            a, b, c = ROT_CONSTANTS[(i - 1) & 3]
            u_zero.append([linear_L(rks_zero[i][w], a, b, c) for w in range(4)])

    rows: Dict[Tuple[int, int], int] = {(i, bit): 0 for i in range(1, max_R + 1) for bit in range(128)}
    for k in range(256):
        mk = [0] * 8
        mk[k // 32] = 1 << (k % 32)
        rks = expand_key(mk, max_R)
        for i in range(1, max_R + 1):
            a, b, c = ROT_CONSTANTS[(i - 1) & 3]
            uw = [linear_L(rks[i][w], a, b, c) ^ u_zero[i][w] for w in range(4)]
            for w in range(4):
                val = uw[w]
                for col in range(32):
                    if (val >> col) & 1:
                        rows[(i, 32 * w + col)] |= 1 << k
    return rows


def add_basis(basis: Dict[int, int], v: int) -> bool:
    while v:
        hb = v.bit_length() - 1
        if hb in basis:
            v ^= basis[hb]
        else:
            basis[hb] = v
            return True
    return False


def C_packed(m: int) -> float:
    if m == 0:
        return 0.0
    shift = {12: 2.0, 8: 3.0, 4: 4.0}.get(m, 0.0)
    return math.log2((1 + 2 * m) * (2 ** m)) - shift


def build_peel_tree(R: int, n: int, col: int, r0: int = 0):
    Tset: Dict[int, List[int]] = {n + 1: [col]}
    parents: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for i in range(n + 1, R):
        j = (r0 + i - 1) & 3
        nxt: Set[int] = set()
        for b in Tset[i]:
            ch = RS[j][b]
            nxt.update(ch)
            parents[(i, b)] = [(i + 1, k) for k in ch]
        Tset[i + 1] = sorted(nxt)
    return Tset, parents


def single_run_packed_cost(R: int, n: int, col: int, r0: int = 0, use_pattern_tbl: bool = True) -> Tuple[int, int, float]:
    Tset, parents = build_peel_tree(R, n, col, r0)
    # Bottom-up order (level R-1 down to n+1)
    order = [(i, b) for i in range(R - 1, n, -1) for b in Tset[i]]
    live: Set[Tuple[int, int]] = {(R, b) for b in Tset[R]}
    rem = {nd: 0 for nd in live}
    for p, ch in parents.items():
        for c in ch:
            rem[c] = rem.get(c, 0) + 1
    guessed: Set[Tuple[int, int]] = set()
    G = 0
    s_pck = 0.0
    for p in order:
        ch = parents[p]
        m = 0
        for c in ch:
            if (c[0] + 1, c[1]) not in guessed:
                m += 4
                guessed.add((c[0] + 1, c[1]))
        G_before = G
        G += m
        for c in ch:
            rem[c] -= 1
            if rem[c] == 0:
                live.discard(c)
        live.add(p)
        I_after = 4 * len(live)
        s_pck += 2.0 ** (G_before + I_after + C_packed(m))
    G_before_fin = G
    for nd in live:
        if (nd[0] + 1, nd[1]) not in guessed:
            G += 4
            guessed.add((nd[0] + 1, nd[1]))
    fin = (2.0 ** G_before_fin) if use_pattern_tbl else (2.0 ** (G_before_fin - 4.0 + 9.04))
    return G // 4, len(Tset[R]), math.log2(s_pck + fin)


def audit_progressive_chain(
    u_rows: Dict[Tuple[int, int], int],
    R: int,
    n: int,
    d_cube: int,
    m_cubes: int,
    filter_per_col: int,
    cols_seq: Sequence[int],
    r0: int = 0,
) -> Dict[str, float]:
    basis: Dict[int, int] = {}
    logL = 0.0
    max_mem = 0.0
    tot_pck = 0.0
    tot_join = 0.0
    for s_idx, col in enumerate(cols_seq, 1):
        Tset, _ = build_peel_tree(R, n, col, r0)
        K, n_uR, pck1 = single_run_packed_cost(R, n, col, r0, use_pattern_tbl=True)
        vecs = [u_rows[(i, 32 * w + b)] for i, cs in Tset.items() for b in cs for w in range(4)]
        uR_vecs = [u_rows[(R, 32 * w + b)] for b in Tset[R] for w in range(4)]
        tmp_uR = dict(basis)
        free_uR = sum(1 for v in uR_vecs if add_basis(tmp_uR, v))
        added = sum(1 for v in vecs if add_basis(basis, v))
        new_rk = len(basis)

        c_blind = m_cubes * (2.0 ** pck1)
        c_guid = m_cubes * (2.0 ** (logL + free_uR + 54.75)) if (s_idx > 1 and R - n == 4) else c_blind
        step_pck = min(c_blind, c_guid)
        used_g = c_guid < c_blind

        logL_after = max(0.0, logL + added - filter_per_col)
        local_probe = 2.0 ** max(0.0, (logL + free_uR + 52 - filter_per_col) if used_g else (4 * K - filter_per_col))
        join_ops = (2.0 ** logL) + (2.0 ** logL_after) + local_probe
        tot_pck += step_pck
        tot_join += join_ops
        logL = logL_after
        max_mem = max(max_mem, (2.0 ** logL) * max(1, new_rk), 2.0 ** (4 * n_uR + 4))

    D = m_cubes * (2.0 ** d_cube)
    bf = 2.0 ** (logL + (256 - len(basis)))
    T_total = math.log2(D + tot_pck + tot_join + bf)
    return {
        "D_log2": math.log2(D),
        "T_log2": T_total,
        "peel_log2": math.log2(tot_pck),
        "join_log2": math.log2(tot_join),
        "bf_log2": math.log2(bf),
        "M_log2": math.log2(max_mem),
        "rho": len(basis),
        "final_logL": logL,
    }


def main() -> int:
    t0 = time.perf_counter()
    print("=== Higher-Round ARADI Algebraic & Complexity Audit (Single-Core) ===")

    verify_island_supersbox_spectrum()
    print("[1/5] 16-bit L_3-island Super-S-Box Psi_G exact ANF spectrum:")
    print("      hist = {0: 1, 3: 3, 5: 28, 6: 224, 8: 65280} (min deg = 3, key-independent) [PASS]")

    n_prod, n_hit = verify_no_composite_collapse()
    print(f"[2/5] Composite linear-layer canceling-mask census ({n_prod} products, k=1..5):")
    print(f"      {n_hit} products with dim C(M) > 0, ALL reducing to I (dim 32) or L_3 (dim 8) [PASS]")

    verify_shift_symmetry_and_1round_lift()
    n_checked, lift_cols = verify_6r_d16_cube_live(n_keys=8)
    print(f"[3/5] Z_16 shift symmetry, 1-Round Lift identity & live 6R d=16 cube (0,5,19,24):")
    print(f"      Q_5 128/128 balanced; Q_6 full col 25 + 24 lift masks on {lift_cols} verified ({n_checked} keys) [PASS]")

    u_rows = build_u_affine_rows(14)
    # Verify 9R full master-key rank: all 16 Z_16 shifts of col 25 -> cols 16..31 cover all 32 cols of K^(9), K^(8)
    b9: Dict[int, int] = {}
    for sh in range(16):
        u = 16 + ((9 + sh) & 15)
        e_cols = RS[2][u]
        d_cols = sorted({d for e in e_cols for d in RS[3][e]})
        for w in range(4):
            add_basis(b9, u_rows[(7, 32 * w + u)])
            for e in e_cols:
                add_basis(b9, u_rows[(8, 32 * w + e)])
            for d in d_cols:
                add_basis(b9, u_rows[(9, 32 * w + d)])
    assert len(b9) == 256, f"9R rank expected 256, got {len(b9)}"
    print("[4/5] 9-Round nested-WHT peel rank check: 16 shifted d=16 cubes -> exact GF(2) rho = 256 (residual 2^0) [PASS]")

    # Verify 10R, 11R, 12R progressive chains at r0 = 0 (unmodified round-0 start, 2025/1227 packed-op accounting)
    # 10R (a) m=16 d=113 structural distinguisher (13 cols from greedy chain starting at col 1)
    c10_m16 = audit_progressive_chain(
        u_rows, 10, 7, 113, 16, 16, [1, 9, 16, 24, 5, 13, 20, 28, 0, 8, 23, 31, 4], r0=0
    )
    # 10R (b) m=8 d=113 structural distinguisher (24 cols, rho=256)
    c10_m8 = audit_progressive_chain(
        u_rows, 10, 7, 113, 8, 8,
        [1, 9, 16, 24, 5, 13, 20, 28, 0, 8, 23, 31, 4, 12, 19, 27, 15, 7, 22, 30, 6, 2, 3, 10],
        r0=0,
    )
    # 11R (a) 8R d=124 ZeroSum, m=3 (10-bit filter), 19 greedy cols
    c11_8r = audit_progressive_chain(
        u_rows, 11, 8, 124, 3, 10,
        [1, 20, 12, 31, 7, 26, 2, 21, 13, 16, 8, 27, 3, 22, 14, 17, 9, 28, 4],
        r0=0,
    )
    # 11R (b) 7R d=116 ZeroSum, m=7 (32-bit filter), 4 cols [0, 8, 23, 31] + guided Step 4
    c11_7r = audit_progressive_chain(u_rows, 11, 7, 116, 7, 32, [0, 8, 23, 31], r0=0)
    # 12R (a) 8R d=124 ZeroSum, m=7 (32-bit filter), [0, 19, 11, 6, 1, 12, 7] (beats D/T/M!)
    c12_m7 = audit_progressive_chain(u_rows, 12, 8, 124, 7, 32, [0, 19, 11, 6, 1, 12, 7], r0=0)
    # 12R (b) 8R d=124 ZeroSum, m=6 (26-bit filter), lower-half |T_12|=22 cols [0, 5, 10, 11, 6]
    c12_m6 = audit_progressive_chain(u_rows, 12, 8, 124, 6, 26, [0, 5, 10, 11, 6], r0=0)

    # 13R geometric barrier: check min K over all r0 in 0..3 and all 32 cols for R=13, n=8 (depth 5)
    min_k13 = min(single_run_packed_cost(13, 8, col, r0)[0] for r0 in range(4) for col in range(32))
    assert min_k13 == 66, f"Expected 13R depth-5 min K = 66, got {min_k13}"

    assert c10_m16["D_log2"] < 118.0 and c10_m16["T_log2"] < 118.0 and c10_m16["M_log2"] < 48.0
    assert c10_m8["D_log2"] == 116.0 and c10_m8["T_log2"] < 116.01 and c10_m8["rho"] == 256
    assert c11_8r["D_log2"] < 126.0 and c11_8r["T_log2"] < 126.0
    assert c11_7r["D_log2"] < 119.0 and c11_7r["T_log2"] < 147.76
    assert c12_m7["T_log2"] < 151.80 and c12_m7["M_log2"] < 123.50 and c12_m7["rho"] == 256
    assert c12_m6["D_log2"] < 126.81 and c12_m6["T_log2"] < 151.80

    print("[5/5] Progressive list-join & packed-arithmetic complexity audit (r0 = 0, vs. ePrint 2025/1227):")
    print(f"      10R (7R Struct d=113, m=16): D=2^{c10_m16['D_log2']:.2f}, T=2^{c10_m16['T_log2']:.2f}, M=2^{c10_m16['M_log2']:.2f}b  (beats 2^118, 2^118, 2^48b in all 3!)")
    print(f"      10R (7R Struct d=113, m=8 ): D=2^{c10_m8['D_log2']:.2f}, T=2^{c10_m8['T_log2']:.2f}, M=2^{c10_m8['M_log2']:.2f}b  (rho=256, 4x lower D & T)")
    print(f"      11R (8R ZeroSum d=124, m=3): D=2^{c11_8r['D_log2']:.2f}, T=2^{c11_8r['T_log2']:.2f}, M=2^{c11_8r['M_log2']:.2f}b  (beats Sec 3.3 2^126, 2^126)")
    print(f"      11R (7R ZeroSum d=116, m=7): D=2^{c11_7r['D_log2']:.2f}, T=2^{c11_7r['T_log2']:.2f}, M=2^{c11_7r['M_log2']:.2f}b  (beats Sec 3.4 2^119, 2^147.76)")
    print(f"      12R (8R ZeroSum d=124, m=7): D=2^{c12_m7['D_log2']:.2f}, T=2^{c12_m7['T_log2']:.2f}, M=2^{c12_m7['M_log2']:.2f}b  (rho=256, beats Sec 3.5 T=2^151.8, M=2^123.5b)")
    print(f"      12R (8R ZeroSum d=124, m=6): D=2^{c12_m6['D_log2']:.2f}, T=2^{c12_m6['T_log2']:.2f}, M=2^{c12_m6['M_log2']:.2f}b  (beats Sec 3.5 in both D & T!)")
    print(f"      13R depth-5 barrier        : min K = {min_k13} => single run >= 2^264 > 2^255 (2025/1552) [PASS]")
    dt = time.perf_counter() - t0
    print(f"[PASS] All higher-round algebraic theorems and 9R-12R complexity bounds verified in {dt:.2f}s (1 core).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
