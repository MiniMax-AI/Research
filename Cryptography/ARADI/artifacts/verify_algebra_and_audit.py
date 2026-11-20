#!/usr/bin/env python3
"""Algebraic proofs, LATINCRYPT 2025 blind-column audit, and L_3 mitigation verification.

Verifies all structural and algebraic claims in the paper using pure Python 3
(plus optional NumPy acceleration when available for large empirical rank checks):
  1. Official NSA ARADI test vectors (all 21 vectors + 17 round keys + 16 round states).
  2. Involutive linear maps Lambda_j^2 = I_128 (j = 0..3), uniqueness of L_3 = (8,9,7)
     as a symmetric orthogonal GF(2) matrix (L_3^T L_3 = I_32) with 8 disjoint 4-column
     islands GROUPS, and 0 non-zero cross-term-canceling masks for L_0, L_1, L_2.
  3. Audit of Bellini et al. (LATINCRYPT 2025, ePrint 2025/1443 Section 6.2):
     across all 32 single-column input differences at r_0 = 0, 3-round dead columns
     cover at most 24/32 columns of lambda = Lambda_0(K^(5)), leaving 8 or 13 blind
     columns (32 or 52 unconstrained key bits).
  4. Orthogonal Island Collapse Theorem: analytical (< 30 ms, zero-sampling) derivation
     of all 198 universal 5-round zero-sum equations (99 on row x + 99 on row z) on
     the 48-point prefix of [(2,15),(2,17),(2,22),(2,25),(2,30),(3,30)], 0 violations
     across random 256-bit master keys and independent 640-bit round keys, offset
     sensitivity (r_0 = 0 vs r_0 = 1, 2, 3), and verification that replacing L_3 = (8,9,7)
     with odd-shift non-orthogonal involution L_3' = (7,3,11) drops 48-CP annihilators
     from 198 > 128 to underdetermined 96 < 128.
"""
from __future__ import annotations

from collections import Counter
import math
import random
import time
from typing import Dict, List, Sequence, Tuple

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]

from aradi_ref import (
    GROUPS,
    GROUP_OF,
    MASK32,
    ROT_CONSTANTS,
    SBOX_INV_TABLE,
    SBOX_TABLE,
    apply_linear_state,
    encrypt_rounds,
    expand_key,
    linear_L,
    sbox,
    verify_official_vectors,
)

CUBE_6D_DIRS: List[Tuple[int, int]] = [
    (2, 15),
    (2, 17),
    (2, 22),
    (2, 25),
    (2, 30),
    (3, 30),
]

L3_PRIME_CANDIDATE: Tuple[int, int, int] = (7, 3, 11)
L3_ALT_CANDIDATE: Tuple[int, int, int] = (3, 1, 5)


# ============================================================================
# GF(2) Matrix & Boolean Polynomial Helpers
# ============================================================================

def word_linear_matrix(a: int, b: int, c: int) -> List[List[int]]:
    """Return 32x32 GF(2) matrix M where M[row][col] is bit `row` of L(1 << col)."""
    mat = [[0] * 32 for _ in range(32)]
    for col in range(32):
        out = linear_L(1 << col, a, b, c)
        for row in range(32):
            mat[row][col] = (out >> row) & 1
    return mat


def mat_mul_gf2(A: List[List[int]], B: List[List[int]]) -> List[List[int]]:
    n = len(A)
    return [
        [sum(A[r][k] & B[k][c] for k in range(n)) & 1 for c in range(n)]
        for r in range(n)
    ]


def mat_transpose(A: List[List[int]]) -> List[List[int]]:
    n = len(A)
    return [[A[r][c] for r in range(n)] for c in range(n)]


def is_identity(A: List[List[int]]) -> bool:
    n = len(A)
    return all(A[r][c] == (1 if r == c else 0) for r in range(n) for c in range(n))


def gf2_rref_int(rows: Sequence[int]) -> Dict[int, int]:
    """Return reduced row echelon pivot dictionary {msb_index: row_int}."""
    pivs: Dict[int, int] = {}
    for row in rows:
        v = row
        while v:
            hb = v.bit_length() - 1
            if hb in pivs:
                v ^= pivs[hb]
            else:
                pivs[hb] = v
                v = 0
    for p in sorted(pivs):
        pv = pivs[p]
        for q in sorted(pivs):
            if q != p and ((pivs[q] >> p) & 1):
                pivs[q] ^= pv
    return pivs


def canceling_mask_subspace_dim(M: List[List[int]]) -> int:
    """Dimension of {mu in F_2^32 : M^T diag(mu) M is diagonal over F_2}."""
    rows = []
    for p in range(32):
        for q in range(32):
            if p == q:
                continue
            mask = 0
            for r in range(32):
                if M[r][p] & M[r][q]:
                    mask |= 1 << r
            if mask:
                rows.append(mask)
    return 32 - len(gf2_rref_int(rows))


def anf_4bit(values: Sequence[int]) -> str:
    """Compute exact 4-variable Algebraic Normal Form over (w, x, y, z)."""
    coeffs = list(values)
    for bit in range(4):
        step = 1 << bit
        for m in range(16):
            if m & step:
                coeffs[m] ^= coeffs[m ^ step]
    var_names = ("w", "x", "y", "z")
    terms = []
    for m in range(16):
        if coeffs[m]:
            term = "".join(var_names[b] for b in range(4) if (m >> b) & 1) or "1"
            terms.append(term)
    return " ^ ".join(terms) if terms else "0"


# ============================================================================
# Analytical 198-Annihilator Certificate via Orthogonal Island Collapse
# ============================================================================

_OR_SHIFT_64: List[List[int]] = [[1 << (i | j) for j in range(64)] for i in range(64)]


def _mul_supp64(s1: int, s2: int) -> int:
    if s1 == 0 or s2 == 0:
        return 0
    if s1 == 1:
        return s2
    if s2 == 1:
        return s1
    res = 0
    m1 = s1
    while m1:
        low1 = m1 & -m1
        row = _OR_SHIFT_64[low1.bit_length() - 1]
        m1 ^= low1
        m2 = s2
        while m2:
            low2 = m2 & -m2
            res |= row[low2.bit_length() - 1]
            m2 ^= low2
    return res


def _sbox_supp64(sw: int, sx: int, sy: int, sz: int) -> Tuple[int, int, int, int]:
    uw, ux, uy, uz = sw | 1, sx | 1, sy | 1, sz | 1
    wy = _mul_supp64(uw, uy)
    xy = _mul_supp64(ux, uy)
    xz = _mul_supp64(ux, uz)
    wz = _mul_supp64(uw, uz)
    pw = uw | wy | xy | xz | _mul_supp64(wy, uz)
    px = ux | wy
    py = uy | wy | wz | _mul_supp64(ux, wy)
    pz = uz | xy | wy
    return pw, px, py, pz


def _lin_supp64(word_supp: Sequence[int], M: List[List[int]]) -> List[int]:
    out = [0] * 32
    for row in range(32):
        acc = 0
        mrow = M[row]
        for col in range(32):
            if mrow[col]:
                acc |= word_supp[col]
        out[row] = acc
    return out


def derive_198_certificate(
    cube_dirs: Sequence[Tuple[int, int]] = CUBE_6D_DIRS,
    num_points: int = 48,
) -> Tuple[List[Tuple[int, int, int]], List[int]]:
    """Derive all universal 5-round zero-sum equations on Q_5 ^ K^(5) analytically.

    Returns `(equations, island_nullities)` where each equation is a tuple
    `(word_row_idx, col_mask_32, point_mask_48)` asserting that
        XOR_{p in point_mask_48} parity(u(P_p)[word_row_idx] & col_mask_32) == 0
    universally for all base plaintexts and all independent round keys.
    """
    mats = [word_linear_matrix(*rc) for rc in ROT_CONSTANTS]
    st = [[1] * 32 for _ in range(4)]
    for t, (w_idx, b_idx) in enumerate(cube_dirs):
        st[w_idx][b_idx] |= 1 << (1 << t)

    # Propagate exact monomial supports through Rounds 0, 1, 2 to Q_3.
    for r in range(3):
        nw, nx, ny, nz = [0] * 32, [0] * 32, [0] * 32, [0] * 32
        for b in range(32):
            if (st[0][b] | st[1][b] | st[2][b] | st[3][b]) == 1:
                nw[b] = nx[b] = ny[b] = nz[b] = 1
            else:
                nw[b], nx[b], ny[b], nz[b] = _sbox_supp64(
                    st[0][b], st[1][b], st[2][b], st[3][b]
                )
        M = mats[r]
        st = [
            _lin_supp64(nw, M),
            _lin_supp64(nx, M),
            _lin_supp64(ny, M),
            _lin_supp64(nz, M),
        ]

    # Restriction vectors of the 64 monomials onto the 48-point prefix 0..47.
    v_mono = [0] * 64
    for m in range(64):
        vec = 0
        for p in range(num_points):
            if (p & m) == m:
                vec |= 1 << p
        v_mono[m] = vec

    M0 = mats[0]
    equations: List[Tuple[int, int, int]] = []
    island_nullities: List[int] = []

    for g in GROUPS:
        # Pulling sum_{b in g} Z_{4, row, b} through Q_5 = Lambda_0(Z_4) (L_0^{-1} = L_0).
        col_mask = 0
        for c in range(32):
            if sum(M0[b][c] for b in g) & 1:
                col_mask |= 1 << c

        # Orthogonal Island Collapse: S_3 -> L_3 -> S_4 on rows x and z has degree <= 3 in Q_3.
        supp_g = 1
        for b in g:
            s1 = st[0][b] | st[1][b] | st[2][b] | st[3][b] | 1
            s2 = _mul_supp64(s1, s1)
            s3 = _mul_supp64(s2, s1)
            supp_g |= s3

        rows = [v_mono[m] for m in range(64) if (supp_g >> m) & 1]
        pivs = gf2_rref_int(rows)
        null_basis: List[int] = []
        for f in range(num_points):
            if f not in pivs:
                vec = 1 << f
                for p, pv in pivs.items():
                    if (pv >> f) & 1:
                        vec |= 1 << p
                null_basis.append(vec)

        island_nullities.append(len(null_basis))
        for pm in null_basis:
            equations.append((1, col_mask, pm))  # row x
            equations.append((3, col_mask, pm))  # row z

    return equations, island_nullities


# ============================================================================
# Empirical Per-Word Annihilator Measurement (48-Point Prefix)
# ============================================================================

def measure_word_annihilators_48cp(
    r0: int = 0,
    rot_constants: Sequence[Tuple[int, int, int]] = ROT_CONSTANTS,
    nsamples: int = 1750,
    seed: int = 20260409,
) -> List[int]:
    """Measure exact per-word 5-round zero-sum nullity across 48-point prefix."""
    n_cols = 32 * 47  # 1504 difference coordinates per word

    if np is not None:
        rng = np.random.default_rng(seed + r0)
        base = rng.integers(0, 1 << 32, size=(4, nsamples), dtype=np.uint32)
        rks = rng.integers(0, 1 << 32, size=(5, 4, nsamples), dtype=np.uint32)
        pts = np.empty((48, 4, nsamples), dtype=np.uint32)
        for m in range(48):
            pts[m] = base
            for t, (w_idx, b_idx) in enumerate(CUBE_6D_DIRS):
                if (m >> t) & 1:
                    pts[m, w_idx] ^= np.uint32(1 << b_idx)

        w = pts[:, 0].copy()
        x = pts[:, 1].copy()
        y = pts[:, 2].copy()
        z = pts[:, 3].copy()
        for step in range(5):
            w ^= rks[step, 0]
            x ^= rks[step, 1]
            y ^= rks[step, 2]
            z ^= rks[step, 3]
            x ^= w & y
            z ^= x & y
            y ^= w & z
            w ^= x & z
            a, b, c = rot_constants[(r0 + step) & 3]
            for arr in (w, x, y, z):
                u = (arr >> 16) & np.uint32(0xFFFF)
                l = arr & np.uint32(0xFFFF)
                ua = ((u << a) | (u >> (16 - a))) & np.uint32(0xFFFF) if a else u
                la = ((l << a) | (l >> (16 - a))) & np.uint32(0xFFFF) if a else l
                lc = ((l << c) | (l >> (16 - c))) & np.uint32(0xFFFF) if c else l
                ub = ((u << b) | (u >> (16 - b))) & np.uint32(0xFFFF) if b else u
                arr[:] = ((u ^ ua ^ lc) << 16) | (l ^ la ^ ub)

        nullities = []
        for word_arr in (w, x, y, z):
            diffs = word_arr[1:] ^ word_arr[0:1]  # (47, nsamples)
            mat = np.zeros((nsamples, n_cols), dtype=np.uint8)
            for bit in range(32):
                mat[:, bit * 47 : (bit + 1) * 47] = ((diffs >> bit) & 1).T
            packed = np.packbits(mat, axis=1)
            int_rows = [int.from_bytes(packed[s].tobytes(), "big") for s in range(nsamples)]
            rk = len(gf2_rref_int(int_rows))
            nullities.append(n_cols - rk)
        return nullities

    # Pure-Python stdlib fallback if NumPy is absent.
    py_rng = random.Random(seed + r0)
    word_rows: List[List[int]] = [[] for _ in range(4)]
    for _ in range(nsamples):
        base_pt = tuple(py_rng.getrandbits(32) for _ in range(4))
        rk_list = [tuple(py_rng.getrandbits(32) for _ in range(4)) for _ in range(5)] + [
            (0, 0, 0, 0)
        ]
        outs = []
        for m in range(48):
            p = list(base_pt)
            for t, (wi, bi) in enumerate(CUBE_6D_DIRS):
                if (m >> t) & 1:
                    p[wi] ^= 1 << bi
            outs.append(
                encrypt_rounds((p[0], p[1], p[2], p[3]), rk_list, 5, r0=r0, rot_constants=rot_constants)
            )
        for wi in range(4):
            v = 0
            u0 = outs[0][wi]
            for p_idx in range(1, 48):
                df = outs[p_idx][wi] ^ u0
                for bit in range(32):
                    if (df >> bit) & 1:
                        v |= 1 << (bit * 47 + (p_idx - 1))
            word_rows[wi].append(v)
    return [n_cols - len(gf2_rref_int(word_rows[wi])) for wi in range(4)]


# ============================================================================
# Main Verification Driver
# ============================================================================

def main() -> int:
    print("========================================================================")
    print("  NSA ARADI Structural Algebra, LATINCRYPT 2025 Audit & Mitigation Check")
    print("========================================================================\n")

    # 1. Official test vectors
    print("[1/4] Official NSA Specification & Test Vectors (ePrint 2024/1240)")
    verify_official_vectors(verbose=True)
    print()

    # 2. Linear-layer involutivity, unique orthogonality of L_3, and super-S-box ANF
    print("[2/4] Linear Layer Algebra & Orthogonal Island Collapse Proof")
    mats = [word_linear_matrix(*rc) for rc in ROT_CONSTANTS]
    for j, (a, b, c) in enumerate(ROT_CONSTANTS):
        M = mats[j]
        assert (2 * a) % 16 == (b + c) % 16, f"L_{j} shift condition failed"
        assert is_identity(mat_mul_gf2(M, M)), f"L_{j}^2 != I_32"
        sym = M == mat_transpose(M)
        orth = is_identity(mat_mul_gf2(mat_transpose(M), M))
        cancel_dim = canceling_mask_subspace_dim(M)
        cancel_masks = (1 << cancel_dim) - 1
        print(
            f"  - L_{j} = {(a, b, c)}: involution=True, symmetric={str(sym):5s}, "
            f"orthogonal={str(orth):5s}, canceling_masks={cancel_masks}"
        )
        if j < 3:
            assert not sym and not orth and cancel_dim == 0
        else:
            assert sym and orth and cancel_dim == 8

    # Check 8-island block decomposition of L_3
    M3 = mats[3]
    assert len(GROUPS) == 8 and sorted(c for g in GROUPS for c in g) == list(range(32))
    for r in range(32):
        supp = [c for c in range(32) if M3[r][c]]
        assert len(supp) == 3, f"Row {r} of L_3 does not have weight 3"
        assert len({GROUP_OF[c] for c in supp}) == 1, "L_3 leaks outside 4-column island"
        # Complement of row r inside its 4-column island is a single missing column
        grp = set(GROUPS[GROUP_OF[r]])
        assert set(supp) < grp and len(grp - set(supp)) == 1

    # Verify exact S-box, inverse S-box, and super-S-box ANFs
    anf_w = anf_4bit([(SBOX_TABLE[u] >> 0) & 1 for u in range(16)])
    anf_x = anf_4bit([(SBOX_TABLE[u] >> 1) & 1 for u in range(16)])
    anf_y = anf_4bit([(SBOX_TABLE[u] >> 2) & 1 for u in range(16)])
    anf_z = anf_4bit([(SBOX_TABLE[u] >> 3) & 1 for u in range(16)])
    anf_super_x = anf_4bit(
        [((SBOX_TABLE[u] >> 1) & 1) ^ (((SBOX_TABLE[u] >> 0) & 1) & ((SBOX_TABLE[u] >> 2) & 1)) for u in range(16)]
    )
    anf_super_z = anf_4bit(
        [
            ((SBOX_TABLE[u] >> 3) & 1)
            ^ ((((SBOX_TABLE[u] >> 1) & 1) ^ ((SBOX_TABLE[u] >> 0) & 1)) & ((SBOX_TABLE[u] >> 2) & 1))
            for u in range(16)
        ]
    )
    assert "wxyz" not in anf_super_x and "wxyz" not in anf_super_z
    assert anf_super_x == "x ^ wy ^ xy ^ wz ^ wxz ^ xyz"
    print(f"  - Verified 8 disjoint 4-column islands of L_3: {GROUPS[0]} .. {GROUPS[-1]}")
    print(f"  - S-box ANF: pi_w = {anf_w}; pi_x = {anf_x}; pi_y = {anf_y}; pi_z = {anf_z}")
    print(f"  - Super-S-box collapse (row x): pi_x ^ pi_w*pi_y = {anf_super_x} (deg = 3)")
    print(f"  - Super-S-box collapse (row z): pi_z ^ (pi_x^pi_w)*pi_y = {anf_super_z} (deg = 3)\n")

    # 3. Audit of Bellini et al. (LATINCRYPT 2025, ePrint 2025/1443 Section 6.2)
    print("[3/4] Audit of Bellini et al. (LATINCRYPT 2025, ePrint 2025/1443 Section 6.2)")
    dead_hist: Counter[int] = Counter()
    cov_hist: Counter[int] = Counter()
    blind_hist: Counter[int] = Counter()
    for c0 in range(32):
        active = {c0}
        for r in (0, 1, 2):
            Mr = mats[r]
            active = {row for row in range(32) if any(Mr[row][col] for col in active)}
        dead = sorted(set(range(32)) - active)
        covered = {col for k in dead for col in range(32) if M3[k][col]}
        blind = 32 - len(covered)
        dead_hist[len(dead)] += 1
        cov_hist[len(covered)] += 1
        blind_hist[blind] += 1

    assert dead_hist == Counter({9: 16, 10: 16})
    assert cov_hist == Counter({19: 16, 24: 16})
    assert blind_hist == Counter({13: 16, 8: 16})
    print(f"  - 3-round dead columns across 32 single-column diffs: {dict(sorted(dead_hist.items()))}")
    print(f"  - Covered columns of lambda = Lambda_0(K^(5)) / 32 : {dict(sorted(cov_hist.items()))}")
    print(
        f"  - Blind columns of lambda (unconstrained nibbles)  : {dict(sorted(blind_hist.items()))} "
        "-> 32 or 52 unconstrained key bits!"
    )
    print("  - [PASS] Confirmed: max coverage is 24/32 columns; 1 single-diff 20-CP batch cannot recover K^(5).\n")

    # 4. Analytical 198-annihilator certificate, random-key tests, r0 offsets, and L_3' fix
    print("[4/4] 6-Round 48-CP Analytical Certificate, Offset Check & Zero-Cost L_3' Fix")
    t0 = time.perf_counter()
    cert_eqs, island_nulls = derive_198_certificate()
    cert_ms = (time.perf_counter() - t0) * 1000.0
    row_x_eqs = sum(1 for r, _, _ in cert_eqs if r == 1)
    row_z_eqs = sum(1 for r, _, _ in cert_eqs if r == 3)
    assert island_nulls == [16, 12, 10, 12, 18, 0, 16, 15]
    assert row_x_eqs == 99 and row_z_eqs == 99 and len(cert_eqs) == 198
    print(
        f"  - Analytical certificate derived in {cert_ms:.2f} ms (0 key samples): "
        f"{len(cert_eqs)} universal equations ({row_x_eqs} row x + {row_z_eqs} row z)"
    )
    print(f"  - Per-island zero-sum dimensions: {island_nulls} (sum = {sum(island_nulls)} per quadratic row)")

    # Verify 0 violations across 50 random 256-bit master keys + 50 independent 640-bit round-key tuples
    rng_v = random.Random(20261105)
    violations_mk = 0
    violations_indep = 0
    for trial in range(50):
        base_pt = tuple(rng_v.getrandbits(32) for _ in range(4))
        pts48 = []
        for m in range(48):
            p = list(base_pt)
            for t, (wi, bi) in enumerate(CUBE_6D_DIRS):
                if (m >> t) & 1:
                    p[wi] ^= 1 << bi
            pts48.append((p[0], p[1], p[2], p[3]))

        mk = [rng_v.getrandbits(32) for _ in range(8)]
        rks5_mk = expand_key(mk, 5)
        u_mk = [encrypt_rounds(pt, rks5_mk, 5, r0=0) for pt in pts48]

        rks5_ind = [tuple(rng_v.getrandbits(32) for _ in range(4)) for _ in range(6)]
        u_ind = [encrypt_rounds(pt, rks5_ind, 5, r0=0) for pt in pts48]

        for row_idx, cm, pm in cert_eqs:
            acc_mk = 0
            acc_ind = 0
            for p_idx in range(48):
                if (pm >> p_idx) & 1:
                    acc_mk ^= (u_mk[p_idx][row_idx] & cm).bit_count() & 1
                    acc_ind ^= (u_ind[p_idx][row_idx] & cm).bit_count() & 1
            violations_mk += acc_mk
            violations_indep += acc_ind

    assert violations_mk == 0 and violations_indep == 0
    print(
        f"  - Empirical verification: {violations_mk} violations / {50 * len(cert_eqs)} checks (256-bit master keys), "
        f"{violations_indep} violations / {50 * len(cert_eqs)} checks (independent 640-bit round keys)"
    )

    # Offset sensitivity r_0 = 0 vs r_0 = 1, 2, 3
    nulls_by_r0 = {r0: measure_word_annihilators_48cp(r0=r0) for r0 in range(4)}
    assert nulls_by_r0[0] == [0, 99, 0, 99] and sum(nulls_by_r0[0]) == 198
    assert sum(nulls_by_r0[1]) == 16 and sum(nulls_by_r0[2]) == 0 and sum(nulls_by_r0[3]) == 0
    print(
        "  - Round-offset sensitivity (per-word [w, x, y, z] annihilators on 48-CP cube): "
        + ", ".join(f"r0={r0}: {nulls_by_r0[r0]} (total={sum(nulls_by_r0[r0])})" for r0 in range(4))
    )

    # Zero-cost mitigation L_3' = (7, 3, 11) (and alternative (3, 1, 5))
    for label, cand in (("L_3' = (7, 3, 11)", L3_PRIME_CANDIDATE), ("L_3'' = (3, 1, 5)", L3_ALT_CANDIDATE)):
        a3, b3, c3 = cand
        assert math.gcd(a3, 16) == 1 and (2 * a3) % 16 == (b3 + c3) % 16 and (2 * a3) % 16 != 0
        M3p = word_linear_matrix(a3, b3, c3)
        assert is_identity(mat_mul_gf2(M3p, M3p)), f"{label} is not an involution"
        assert all(sum(M3p[r]) == 3 for r in range(32)), f"{label} row weight != 3"
        assert canceling_mask_subspace_dim(M3p) == 0, f"{label} still has canceling masks"
        rot_mod = [ROT_CONSTANTS[0], ROT_CONSTANTS[1], ROT_CONSTANTS[2], cand]
        nulls_mod = measure_word_annihilators_48cp(r0=0, rot_constants=rot_mod)
        assert nulls_mod == [24, 24, 24, 24] and sum(nulls_mod) == 96
        print(
            f"  - Mitigation {label}: involution=True, XORs/bit=2, canceling_masks=0, "
            f"48-CP annihilators={nulls_mod} (total={sum(nulls_mod)} < 128 -> attack neutralized)"
        )

    print("\n[PASS] All algebraic theorems, SOTA audit claims, and mitigation checks verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
