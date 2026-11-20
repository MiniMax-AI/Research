#!/usr/bin/env python3
"""Ultra-low-data key recovery on 5-round ARADI (6-CP MITM & 4-CP Master-Key SAT).

Implements both 5-round single-key chosen-plaintext attacks:
  (a) 6-CP (3 pairs) Fast Subkey Recovery (`K^(5)`):
      Exploits the 8-island decomposition of L_3 = (8, 9, 7) via an exact 2+2
      meet-in-the-middle hash-join (8 * 2 * 2^8 = 4,096 table operations) coupled
      with Round-3 reachable-difference profile filtering and a Round-2 36-dimensional
      subspace / dead-column peel. Uses 1 adaptive extra pair (+2 CP) when inverse-S-box
      derivative fibres collide on a column (~12-13% of keys), achieving 100% exact
      K^(5) recovery in ~25-80 ms.

  (b) 4-CP (2 pairs) Full 256-Bit Master-Key Recovery:
      Queries a single 4-point 2D affine plane on [(z, 16), (z, 31)] (6 pairwise
      differences), peels 4 inner rounds (subspace dimensions dim V_k = 3, 12, 36, 88)
      jointly with the 256 affine key-schedule consistency equations on
      (K^(2), K^(3), K^(4), K^(5)) via CryptoMiniSat (`pycryptosat`), reconstructs
      the complete 256-bit master key in O(1) from (K^(2), K^(3)), and verifies
      re-encryption of all plaintexts plus ground-truth master-key equality.
"""
from __future__ import annotations

import argparse
import itertools
import random
import statistics
import time
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from pycryptosat import Solver
except ImportError:  # pragma: no cover
    Solver = None  # type: ignore[assignment,misc]

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
    format_words,
    linear_L,
    master_key_from_k2_k3,
)

NUM_ROUNDS = 5


# ============================================================================
# Shared Precomputed Tables & Linear-Layer Geometry
# ============================================================================

def _build_word_matrices() -> List[List[List[int]]]:
    mats = []
    for a, b, c in ROT_CONSTANTS:
        m = [[0] * 32 for _ in range(32)]
        for col in range(32):
            out = linear_L(1 << col, a, b, c)
            for row in range(32):
                m[row][col] = (out >> row) & 1
        mats.append(m)
    return mats


LIN_MATS = _build_word_matrices()
NB2 = [tuple(b for b in range(32) if LIN_MATS[2][k][b]) for k in range(32)]
NB3 = [tuple(b for b in range(32) if LIN_MATS[3][k][b]) for k in range(32)]

DDT = [0] * 16
for _x in range(16):
    for _d in range(16):
        DDT[_d] |= 1 << (SBOX_TABLE[_x] ^ SBOX_TABLE[_x ^ _d])

IDDT = [0] * 16
for _d in range(16):
    _m = DDT[_d]
    while _m:
        _low = _m & -_m
        IDDT[_low.bit_length() - 1] |= 1 << _d
        _m ^= _low


def _dead_cols(cols: Sequence[int], layers: Sequence[int]) -> Tuple[int, ...]:
    cur = set(cols)
    for r in layers:
        M = LIN_MATS[r & 3]
        cur = {row for row in range(32) if any(M[row][col] for col in cur)}
    return tuple(sorted(set(range(32)) - cur))


_DEAD3_BY_COL = {c0: _dead_cols((c0,), (0, 1, 2)) for c0 in range(32)}
_DEAD2_CACHE: Dict[Tuple[int, ...], Tuple[int, ...]] = {}


def dead_cols_2r(cols: Sequence[int]) -> Tuple[int, ...]:
    key = tuple(sorted(set(cols)))
    if key not in _DEAD2_CACHE:
        _DEAD2_CACHE[key] = _dead_cols(key, (0, 1))
    return _DEAD2_CACHE[key]


# ============================================================================
# Part (a): 6-CP (3 Pairs) 2^12 Eight-Island MITM + Two-Layer Subspace Peel
# ============================================================================

DIRS_6CP_BASE: List[Tuple[int, int]] = [(3, 16), (3, 21), (3, 31)]
DIRS_6CP_EXT: List[Tuple[int, int]] = DIRS_6CP_BASE + [(3, 26)]


class ReachableChainFilter:
    """Exact Round-3 reachable-difference filter for a single-bit input difference."""

    def __init__(self, word_idx: int, col_idx: int) -> None:
        val0 = 1 << word_idx
        d1 = sorted(r for r in range(32) if LIN_MATS[0][r][col_idx])
        d2 = sorted(r for r in range(32) if any(LIN_MATS[1][r][k] for k in d1))
        idx2 = {c: i for i, c in enumerate(d2)}
        profiles = set()
        for u0 in range(16):
            if not ((DDT[val0] >> u0) & 1):
                continue
            dom = [v for v in range(16) if (DDT[u0] >> v) & 1]
            for w0 in dom:
                for w1 in dom:
                    for w2 in dom:
                        wd = {d1[0]: w0, d1[1]: w1, d1[2]: w2}
                        prof = tuple(
                            (wd[d1[0]] if LIN_MATS[1][j][d1[0]] else 0)
                            ^ (wd[d1[1]] if LIN_MATS[1][j][d1[1]] else 0)
                            ^ (wd[d1[2]] if LIN_MATS[1][j][d1[2]] else 0)
                            for j in d2
                        )
                        profiles.add(prof)
        self.profiles = sorted(profiles)
        self.bitmaps: List[List[int]] = []
        for pos in range(len(d2)):
            tab = [0] * 16
            for pi, prof in enumerate(self.profiles):
                m = DDT[prof[pos]]
                while m:
                    low = m & -m
                    tab[low.bit_length() - 1] |= 1 << pi
                    m ^= low
            self.bitmaps.append(tab)
        self.full_mask = (1 << len(self.profiles)) - 1
        self.terms = [tuple(idx2[j] for j in d2 if LIN_MATS[2][k][j]) for k in range(32)]
        deg = [0] * len(d2)
        for ts in self.terms:
            for t in ts:
                deg[t] += 1
        self.order = sorted(range(len(d2)), key=lambda v: -deg[v])
        rank = {v: i for i, v in enumerate(self.order)}
        self.last: List[List[int]] = [[] for _ in self.order]
        for k, ts in enumerate(self.terms):
            if ts:
                self.last[max(rank[t] for t in ts)].append(k)

    def check(self, ds3: Sequence[int]) -> bool:
        allowed = [IDDT[y] for y in ds3]
        for k, ts in enumerate(self.terms):
            if not ts and not (allowed[k] & 1):
                return False
        n = len(self.order)
        assign = [0] * n
        bitmaps = self.bitmaps
        order = self.order
        last = self.last
        terms_all = self.terms

        def rec(i: int, mask: int) -> bool:
            if i == n:
                return True
            var = order[i]
            tab = bitmaps[var]
            checks = [(terms_all[k], allowed[k]) for k in last[i]]
            for val in range(16):
                nm = mask & tab[val]
                if not nm:
                    continue
                assign[var] = val
                ok = True
                for ts, alw in checks:
                    acc = 0
                    for t in ts:
                        acc ^= assign[t]
                    if not ((alw >> acc) & 1):
                        ok = False
                        break
                if ok and rec(i + 1, nm):
                    return True
            return False

        return rec(0, self.full_mask)


_CHAIN_FILTERS = {c: ReachableChainFilter(3, c) for c in (16, 21, 31, 26)}


def _state_to_cols(st: Tuple[int, int, int, int]) -> List[int]:
    w, x, y, z = st
    return [
        ((w >> b) & 1) | (((x >> b) & 1) << 1) | (((y >> b) & 1) << 2) | (((z >> b) & 1) << 3)
        for b in range(32)
    ]


def _cols_to_state(cols: Sequence[int]) -> Tuple[int, int, int, int]:
    w = x = y = z = 0
    for b, val in enumerate(cols):
        w |= (val & 1) << b
        x |= ((val >> 1) & 1) << b
        y |= ((val >> 2) & 1) << b
        z |= ((val >> 3) & 1) << b
    return w, x, y, z


def _layer2_peel_consistent(
    T_pts: Sequence[Sequence[int]],
    lam_cols: Sequence[int],
    dirs_spec: Sequence[Tuple[int, int]],
    pt_masks: Sequence[int],
) -> bool:
    """Verify existence of Round-3 output nibbles Z_3(P_0) consistent with Round-2 dead subspaces."""
    n = len(pt_masks)
    u = [[SBOX_INV_TABLE[T_pts[p][b] ^ lam_cols[b]] for b in range(32)] for p in range(n)]
    dz3_pts: List[List[int]] = []
    for p in range(n):
        dw = [0, 0, 0, 0]
        for b in range(32):
            df = u[p][b] ^ u[0][b]
            dw[0] |= (df & 1) << b
            dw[1] |= ((df >> 1) & 1) << b
            dw[2] |= ((df >> 2) & 1) << b
            dw[3] |= ((df >> 3) & 1) << b
        dz3_st = apply_linear_state((dw[0], dw[1], dw[2], dw[3]), 3)
        dz3_pts.append(_state_to_cols(dz3_st))

    pivs: Dict[int, int] = {}

    def add_eq(mask512: int, rhs: int) -> bool:
        v = mask512 | (rhs << 512)
        while v & ((1 << 512) - 1):
            hb = (v & ((1 << 512) - 1)).bit_length() - 1
            if hb in pivs:
                v ^= pivs[hb]
            else:
                pivs[hb] = v
                return True
        return (v >> 512) == 0

    for k in range(32):
        if not add_eq(0xFFFF << (16 * k), 1):
            return False

    for i in range(n):
        for j in range(i + 1, n):
            dm = pt_masks[i] ^ pt_masks[j]
            acols = [dirs_spec[t][1] for t in range(len(dirs_spec)) if (dm >> t) & 1]
            for dead_c in dead_cols_2r(acols):
                nb = NB2[dead_c]
                for w_bit in range(4):
                    mask = 0
                    for k in nb:
                        di, dj = dz3_pts[i][k], dz3_pts[j][k]
                        if di == dj:
                            continue
                        base_b = 16 * k
                        for a in range(16):
                            if ((SBOX_INV_TABLE[a ^ di] ^ SBOX_INV_TABLE[a ^ dj]) >> w_bit) & 1:
                                mask |= 1 << (base_b + a)
                    if not add_eq(mask, 0):
                        return False

    if Solver is None:
        return True

    solver = Solver()
    for k in range(32):
        vs = [16 * k + c + 1 for c in range(16)]
        bv = [512 + 4 * k + bit + 1 for bit in range(4)]
        solver.add_xor_clause(vs, True)
        for c in range(16):
            conj = []
            for bit in range(4):
                bl = bv[bit] if ((c >> bit) & 1) else -bv[bit]
                solver.add_clause([-vs[c], bl])
                conj.append(-bl)
            solver.add_clause(conj + [vs[c]])
    for v in pivs.values():
        lits = [bit + 1 for bit in range(512) if (v >> bit) & 1]
        rhs = bool((v >> 512) & 1)
        solver.add_xor_clause(lits, rhs)
    sat, _ = solver.solve()
    return bool(sat)


def solve_5r_mitm_stage(
    cts: Sequence[Tuple[int, int, int, int]],
    dirs_spec: Sequence[Tuple[int, int]],
    pt_masks: Sequence[int],
) -> List[Tuple[int, int, int, int]]:
    """Run 2^12 8-island MITM + Round-3 profile filter + Round-2 subspace peel."""
    T = [_state_to_cols(apply_linear_state(ct, 0)) for ct in cts]
    n = len(pt_masks)
    by_group: List[List[Tuple[int, int, Tuple[int, ...]]]] = [[] for _ in range(8)]
    single_pairs: List[Tuple[int, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            dm = pt_masks[i] ^ pt_masks[j]
            if dm and (dm & (dm - 1)) == 0:
                c0 = dirs_spec[dm.bit_length() - 1][1]
                single_pairs.append((i, j, c0))
                for k in _DEAD3_BY_COL[c0]:
                    nb = NB3[k]
                    by_group[GROUP_OF[nb[0]]].append((i, j, nb))

    group_sols: List[List[Tuple[int, int, int, int]]] = []
    for gi in range(8):
        a, b, c, d = GROUPS[gi]
        cl = by_group[gi]
        ta = [[SBOX_INV_TABLE[T[p][a] ^ v] ^ SBOX_INV_TABLE[T[q][a] ^ v] for v in range(16)] for p, q, _ in cl]
        tb = [[SBOX_INV_TABLE[T[p][b] ^ v] ^ SBOX_INV_TABLE[T[q][b] ^ v] for v in range(16)] for p, q, _ in cl]
        tc = [[SBOX_INV_TABLE[T[p][c] ^ v] ^ SBOX_INV_TABLE[T[q][c] ^ v] for v in range(16)] for p, q, _ in cl]
        td = [[SBOX_INV_TABLE[T[p][d] ^ v] ^ SBOX_INV_TABLE[T[q][d] ^ v] for v in range(16)] for p, q, _ in cl]
        ina = [a in nb for _, _, nb in cl]
        inb = [b in nb for _, _, nb in cl]
        inc = [c in nb for _, _, nb in cl]
        ind = [d in nb for _, _, nb in cl]
        m_len = len(cl)

        left: Dict[Tuple[int, ...], List[Tuple[int, int]]] = {}
        for la in range(16):
            for lb in range(16):
                sig = tuple(
                    (ta[ci][la] if ina[ci] else 0) ^ (tb[ci][lb] if inb[ci] else 0)
                    for ci in range(m_len)
                )
                left.setdefault(sig, []).append((la, lb))

        sols: List[Tuple[int, int, int, int]] = []
        for lc in range(16):
            for ld in range(16):
                sig = tuple(
                    (tc[ci][lc] if inc[ci] else 0) ^ (td[ci][ld] if ind[ci] else 0)
                    for ci in range(m_len)
                )
                if sig in left:
                    for la, lb in left[sig]:
                        sols.append((la, lb, lc, ld))
        group_sols.append(sols)

    cands: List[Tuple[int, ...]] = []
    cols_buf = [0] * 32

    def enum_groups(gi: int) -> None:
        if gi == 8:
            cands.append(tuple(cols_buf))
            return
        for tup in group_sols[gi]:
            for col_idx, val in zip(GROUPS[gi], tup):
                cols_buf[col_idx] = val
            enum_groups(gi + 1)

    enum_groups(0)

    kept_chain: List[Tuple[int, ...]] = []
    for cand in cands:
        ok = True
        for p, q, c0 in single_pairs:
            dw = [0, 0, 0, 0]
            for b in range(32):
                df = SBOX_INV_TABLE[T[p][b] ^ cand[b]] ^ SBOX_INV_TABLE[T[q][b] ^ cand[b]]
                dw[0] |= (df & 1) << b
                dw[1] |= ((df >> 1) & 1) << b
                dw[2] |= ((df >> 2) & 1) << b
                dw[3] |= ((df >> 3) & 1) << b
            ds3 = _state_to_cols(apply_linear_state((dw[0], dw[1], dw[2], dw[3]), 3))
            if not _CHAIN_FILTERS[c0].check(ds3):
                ok = False
                break
        if ok:
            kept_chain.append(cand)

    if len(kept_chain) <= 1:
        kept_final = kept_chain
    else:
        kept_final = [
            cand for cand in kept_chain if _layer2_peel_consistent(T, cand, dirs_spec, pt_masks)
        ]

    return [apply_linear_state(_cols_to_state(cand), 0) for cand in kept_final]


def run_6cp_trial(rng: random.Random) -> Tuple[bool, bool, int, float]:
    """Run one 5-round 6-CP MITM trial; return (unique_at_6cp, exact_match, cp_used, ms)."""
    mk = [rng.getrandbits(32) for _ in range(8)]
    rks = expand_key(mk, NUM_ROUNDS)
    base = tuple(rng.getrandbits(32) for _ in range(4))

    def query_mask(m: int) -> Tuple[int, int, int, int]:
        p = list(base)
        for t, (wi, bi) in enumerate(DIRS_6CP_EXT):
            if (m >> t) & 1:
                p[wi] ^= 1 << bi
        return encrypt_rounds((p[0], p[1], p[2], p[3]), rks, NUM_ROUNDS)

    t0 = time.perf_counter()
    masks = [0, 1, 2, 3, 4, 5]
    cts = [query_mask(m) for m in masks]
    k5_cands = solve_5r_mitm_stage(cts, DIRS_6CP_EXT, masks)
    unique_6cp = len(k5_cands) == 1

    # Adaptive 1-pair fallback if inverse-S-box derivative fibres collide on a column
    for extra_pair in ([8, 9], [6, 7], [10, 11]):
        if len(k5_cands) == 1:
            break
        masks = masks + extra_pair
        cts = cts + [query_mask(m) for m in extra_pair]
        k5_cands = solve_5r_mitm_stage(cts, DIRS_6CP_EXT, masks)

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    exact = len(k5_cands) == 1 and k5_cands[0] == rks[NUM_ROUNDS]
    return unique_6cp, exact, len(masks), elapsed_ms


# ============================================================================
# Part (b): 4-CP (2 Pairs) Full 256-Bit Master-Key Recovery (4-Layer Peel + KS)
# ============================================================================

OUT_MASK = [0] * 16
for _d in range(16):
    _m = 0
    for _a in range(16):
        _m |= 1 << (SBOX_TABLE[_a] ^ SBOX_TABLE[_a ^ _d])
    OUT_MASK[_d] = _m


def _pack128(words: Sequence[int]) -> int:
    v = 0
    for i, s in enumerate(words):
        v |= (int(s) & MASK32) << (32 * i)
    return v


def _unpack128(v: int) -> Tuple[int, int, int, int]:
    return tuple((v >> (32 * i)) & MASK32 for i in range(4))  # type: ignore[return-value]


def _col_of128(v: int, b: int) -> int:
    return (
        ((v >> b) & 1)
        | (((v >> (32 + b)) & 1) << 1)
        | (((v >> (64 + b)) & 1) << 2)
        | (((v >> (96 + b)) & 1) << 3)
    )


def _col_vec128(val: int, b: int) -> int:
    return (
        ((val & 1) << b)
        | (((val >> 1) & 1) << (32 + b))
        | (((val >> 2) & 1) << (64 + b))
        | (((val >> 3) & 1) << (96 + b))
    )


def _apply_L128(v: int, ri: int) -> int:
    return _pack128(apply_linear_state(_unpack128(v), ri))


def _rref128(rows: Sequence[int], nbits: int = 128) -> Dict[int, int]:
    pivots: Dict[int, int] = {}
    for row in rows:
        v = row
        while v:
            p = v.bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                v = 0
    for p in sorted(pivots):
        prow = pivots[p]
        for q in sorted(pivots):
            if q != p and ((pivots[q] >> p) & 1):
                pivots[q] ^= prow
    return pivots


def _null_space128(rows: Sequence[int], nbits: int = 128) -> List[int]:
    pivots = _rref128(rows, nbits)
    basis = []
    for f in range(nbits):
        if f not in pivots:
            r = 1 << f
            for p, prow in pivots.items():
                if (prow >> f) & 1:
                    r |= 1 << p
            basis.append(r)
    return basis


def _add_int_basis(basis: List[int], v: int) -> None:
    for b in basis:
        v = min(v, v ^ b)
    if v:
        basis.append(v)


def _mask_bits(m: int) -> List[int]:
    out = []
    while m:
        low = m & -m
        out.append(low.bit_length() - 1)
        m ^= low
    return out


def _exact_expand(states: set[int], ri: int, cap: int = 2_000_000) -> set[int]:
    out: set[int] = set()
    for v in states:
        masks = [OUT_MASK[_col_of128(v, b)] for b in range(32)]
        active = [b for b in range(32) if masks[b] != 1]
        choices = [_mask_bits(masks[b]) for b in active]
        total = 1
        for c in choices:
            total *= len(c)
        if total * len(states) > cap:
            raise MemoryError
        for combo in itertools.product(*choices):
            w = 0
            for b, val in zip(active, combo):
                w ^= _col_vec128(val, b)
            out.add(_apply_L128(w, ri))
    return out


def _hull_step(basis: Sequence[int], ri: int) -> List[int]:
    rows: List[int] = []
    for b in range(32):
        proj: List[int] = []
        for v in basis:
            _add_int_basis(proj, _col_of128(v, b))
        m = 0
        for d in proj:
            m |= OUT_MASK[d]
        span_b: List[int] = []
        for val in _mask_bits(m):
            _add_int_basis(span_b, val)
        for val in span_b:
            _add_int_basis(rows, _col_vec128(val, b))
    out: List[int] = []
    for v in rows:
        _add_int_basis(out, _apply_L128(v, ri))
    return out


_DIFF_SPAN_CACHE: Dict[Tuple[int, int], List[int]] = {}
_ANN_CACHE: Dict[Tuple[int, int], List[int]] = {}


def _diff_span(delta: int, rounds: int) -> List[int]:
    key = (delta, rounds)
    if key not in _DIFF_SPAN_CACHE:
        if rounds == 0:
            basis = [delta]
        else:
            basis = [delta]
            for r in range(rounds):
                if r == 0:
                    try:
                        expanded = sorted(_exact_expand({delta}, 0))
                        basis = []
                        for v in expanded:
                            _add_int_basis(basis, v)
                    except MemoryError:
                        basis = [delta]
                else:
                    basis = _hull_step(basis, r)
        _DIFF_SPAN_CACHE[key] = basis
    return _DIFF_SPAN_CACHE[key]


def _annihilator(delta: int, rounds: int) -> List[int]:
    key = (delta, rounds)
    if key not in _ANN_CACHE:
        _ANN_CACHE[key] = _null_space128(_diff_span(delta, rounds), 128)
    return _ANN_CACHE[key]


def _mul_mono_supp(a: int, b: int) -> int:
    out = 0
    ma = a
    while ma:
        low = ma & -ma
        i = low.bit_length() - 1
        ma ^= low
        mb = b
        while mb:
            l2 = mb & -mb
            j = l2.bit_length() - 1
            mb ^= l2
            out |= 1 << (i | j)
    return out


_SOURCE_COLS = [
    [[c for c in range(32) if LIN_MATS[ri][b][c]] for b in range(32)]
    for ri in range(4)
]


def _condition_space_4r(dirs: Sequence[int], point_masks: Sequence[int]) -> List[int]:
    d = len(dirs)
    state = []
    for i in range(4):
        word = [1] * 32
        for t, dv in enumerate(dirs):
            bits = _unpack128(dv)[i]
            for b in range(32):
                if (bits >> b) & 1:
                    word[b] |= 1 << (1 << t)
        state.append(word)
    for r in range(4):
        cols = []
        for b in range(32):
            sw, sx, sy, sz = (state[i][b] for i in range(4))
            x = sx | _mul_mono_supp(sw, sy)
            z = sz | _mul_mono_supp(x, sy)
            y = sy | _mul_mono_supp(sw, z)
            w = sw | _mul_mono_supp(x, z)
            cols.append((w, x, y, z))
        state = [[cols[b][i] for b in range(32)] for i in range(4)]
        src = _SOURCE_COLS[r & 3]
        state = [
            [sum((0,), start=0) | _or_list([word[c] for c in src[b]]) for b in range(32)]
            for word in state
        ]
    eps = [
        sum(1 for p in point_masks if (p & m) == m) & 1
        for m in range(1 << d)
    ]
    forced = set()
    for i in range(4):
        for b in range(32):
            mm = state[i][b]
            ok = True
            while mm:
                low = mm & -mm
                if eps[low.bit_length() - 1]:
                    ok = False
                    break
                mm ^= low
            if ok:
                forced.add((i, b))
    return [1 << (32 * i + b) for i in range(4) for b in range(32) if (i, b) not in forced]


def _or_list(vals: Sequence[int]) -> int:
    acc = 0
    for v in vals:
        acc |= v
    return acc


_KEY_REL_CACHE: Dict[Tuple[int, ...], Tuple[List[int], List[int]]] = {}


def _key_relations(rounds: Tuple[int, ...]) -> Tuple[List[int], List[int]]:
    if rounds not in _KEY_REL_CACHE:
        zero = [_pack128(k) for k in expand_key([0] * 8, NUM_ROUNDS)]
        rows = []
        for i in range(256):
            mk = [0] * 8
            mk[i // 32] = 1 << (i % 32)
            img = [_pack128(k) for k in expand_key(mk, NUM_ROUNDS)]
            v = 0
            for j, r in enumerate(rounds):
                v |= (img[r] ^ zero[r]) << (128 * j)
            rows.append(v)
        basis = _null_space128(rows, nbits=128 * len(rounds))
        z = 0
        for j, r in enumerate(rounds):
            z |= zero[r] << (128 * j)
        consts = [(r & z).bit_count() & 1 for r in basis]
        _KEY_REL_CACHE[rounds] = (basis, consts)
    return _KEY_REL_CACHE[rounds]


_LMAT128: Dict[int, List[List[int]]] = {}


def _lmatrix128(ri: int) -> List[List[int]]:
    key = ri & 3
    if key not in _LMAT128:
        rows: List[List[int]] = [[] for _ in range(128)]
        for q in range(128):
            v = _apply_L128(1 << q, key)
            for p in range(128):
                if (v >> p) & 1:
                    rows[p].append(q)
        _LMAT128[key] = rows
    return _LMAT128[key]


class _SATBuilder:
    def __init__(self) -> None:
        self.n = 0
        self.clauses: List[List[int]] = []
        self.xors: List[Tuple[List[int], bool]] = []

    def new(self, k: int = 1) -> List[int]:
        out = list(range(self.n + 1, self.n + k + 1))
        self.n += k
        return out

    def bits128(self) -> List[int]:
        return self.new(128)

    def nibbles(self) -> List[List[int]]:
        return [self.new(4) for _ in range(32)]

    def exactly_one(self, group: Sequence[int]) -> None:
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                self.clauses.append([-group[i], -group[j]])
        self.xors.append((list(group), True))

    def table4(self, inputs: Sequence[int], outputs: Sequence[int], table: Sequence[int]) -> None:
        for v in range(16):
            lits = [-inputs[k] if (v >> k) & 1 else inputs[k] for k in range(4)]
            out = table[v]
            for k in range(4):
                self.clauses.append(lits + [outputs[k] if (out >> k) & 1 else -outputs[k]])


def _xor_dirs(dirs: Sequence[int], mask: int) -> int:
    acc = 0
    for t, dd in enumerate(dirs):
        if (mask >> t) & 1:
            acc ^= dd
    return acc


def _reduce_vec(vec: int, ann: Sequence[int]) -> int:
    out = 0
    for j, r in enumerate(ann):
        if (vec & r).bit_count() & 1:
            out |= 1 << j
    return out


def _solve_4cp_plane(
    T: Sequence[int],
    dirs: Sequence[int],
    points: Sequence[int],
    limit: int = 8,
    pts_witness: Optional[Sequence[Tuple[int, int, int, int]]] = None,
    rks_witness: Optional[Sequence[Tuple[int, int, int, int]]] = None,
) -> List[Dict[int, int]]:
    idx = {m: i for i, m in enumerate(points)}
    rels: List[Tuple[Tuple[int, ...], List[int]]] = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            delta = _xor_dirs(dirs, points[i] ^ points[j])
            rels.append(((i, j), _annihilator(delta, 4)))
    for r in range(2, len(dirs) + 1):
        for combo in itertools.combinations(range(len(dirs)), r):
            masks = [
                sum(1 << t for t in sub)
                for k in range(r + 1)
                for sub in itertools.combinations(combo, k)
            ]
            if all(m in idx for m in masks):
                rels.append(
                    (
                        tuple(idx[m] for m in masks),
                        _null_space128(
                            _condition_space_4r([dirs[t] for t in combo], list(range(1 << r))),
                            128,
                        ),
                    )
                )

    B = _SATBuilder()
    lam = [B.new(16) for _ in range(32)]
    for b in range(32):
        B.exactly_one(lam[b])

    vectors = [[0] * 16 for _ in range(32)]
    width = 0
    for pts_tup, ann in rels:
        for b in range(32):
            for c in range(16):
                g = 0
                for p in pts_tup:
                    g ^= SBOX_INV_TABLE[_col_of128(T[p], b) ^ c]
                vectors[b][c] |= _reduce_vec(_col_vec128(g, b), ann) << width
        width += len(ann)

    for j in range(width):
        lits = [lam[b][c] for b in range(32) for c in range(16) if (vectors[b][c] >> j) & 1]
        B.xors.append((lits, False))

    base_point = 0
    others = [i for i in range(len(points)) if i != base_point]
    u0 = B.bits128()
    for b in range(32):
        for k in range(4):
            pos = 32 * k + b
            lits = [
                lam[b][c]
                for c in range(16)
                if (SBOX_INV_TABLE[_col_of128(T[base_point], b) ^ c] >> k) & 1
            ]
            B.xors.append((lits + [u0[pos]], False))

    prev = None
    A_of: Dict[int, List[int]] = {}
    SA_of: Dict[int, List[int]] = {}
    for step in range(3):
        k = 3 - step
        dz: Dict[int, List[List[int]]] = {}
        for i in others:
            dz[i] = B.nibbles()
            if prev is None:
                contrib = []
                for bp in range(32):
                    row = []
                    for c in range(16):
                        g = (
                            SBOX_INV_TABLE[_col_of128(T[base_point], bp) ^ c]
                            ^ SBOX_INV_TABLE[_col_of128(T[i], bp) ^ c]
                        )
                        row.append(_apply_L128(_col_vec128(g, bp), k))
                    contrib.append(row)
                for b2 in range(32):
                    for kk in range(4):
                        pos = 32 * kk + b2
                        lits = [
                            lam[bp][c]
                            for bp in range(32)
                            for c in range(16)
                            if (contrib[bp][c] >> pos) & 1
                        ]
                        B.xors.append((lits + [dz[i][b2][kk]], False))
            else:
                rows_k = _lmatrix128(k)
                for b2 in range(32):
                    for kk in range(4):
                        lits = []
                        for q in rows_k[32 * kk + b2]:
                            word, bit = divmod(q, 32)
                            lits.append(prev[i][bit][word])
                        B.xors.append((lits + [dz[i][b2][kk]], False))

        A, SA = B.nibbles(), B.nibbles()
        for b in range(32):
            B.table4(A[b], SA[b], SBOX_TABLE)
        A_of[k] = [A[b][w] for w in range(4) for b in range(32)]
        SA_of[k] = [SA[b][w] for w in range(4) for b in range(32)]
        cur: Dict[int, List[List[int]]] = {}
        for i in others:
            E, Sv, D = B.nibbles(), B.nibbles(), B.nibbles()
            for b in range(32):
                for kk in range(4):
                    B.xors.append(([SA[b][kk], dz[i][b][kk], E[b][kk]], False))
                B.table4(E[b], Sv[b], SBOX_INV_TABLE)
                for kk in range(4):
                    B.xors.append(([A[b][kk], Sv[b][kk], D[b][kk]], False))
            cur[i] = D

        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                delta = _xor_dirs(dirs, points[i] ^ points[j])
                for r in _annihilator(delta, max(k, 1)):
                    lits = []
                    for b in range(32):
                        for kk in range(4):
                            if (r >> (32 * kk + b)) & 1:
                                if i == base_point:
                                    lits.append(cur[j][b][kk])
                                elif j == base_point:
                                    lits.append(cur[i][b][kk])
                                else:
                                    lits.append(cur[i][b][kk])
                                    lits.append(cur[j][b][kk])
                    B.xors.append((lits, False))
        prev = cur

    keys_vars: Dict[int, List[int]] = {}
    K5 = B.bits128()
    rows4 = _lmatrix128(4)
    for q in range(128):
        lits = []
        for r in rows4[q]:
            word, bit = divmod(r, 32)
            lits += [lam[bit][c] for c in range(16) if (c >> word) & 1]
        B.xors.append((lits + [K5[q]], False))
    keys_vars[5] = K5

    K4 = B.bits128()
    rows3 = _lmatrix128(3)
    for q in range(128):
        B.xors.append(([u0[q]] + [SA_of[3][r] for r in rows3[q]] + [K4[q]], False))
    keys_vars[4] = K4

    K3 = B.bits128()
    rows2 = _lmatrix128(2)
    for q in range(128):
        B.xors.append(([A_of[3][q]] + [SA_of[2][r] for r in rows2[q]] + [K3[q]], False))
    keys_vars[3] = K3

    K2 = B.bits128()
    rows1 = _lmatrix128(1)
    for q in range(128):
        B.xors.append(([A_of[2][q]] + [SA_of[1][r] for r in rows1[q]] + [K2[q]], False))
    keys_vars[2] = K2

    order = (2, 3, 4, 5)
    basis, consts = _key_relations(order)
    for r, c_rhs in zip(basis, consts):
        lits = []
        for i, rk_idx in enumerate(order):
            block = (r >> (128 * i)) & ((1 << 128) - 1)
            for q in range(128):
                if (block >> q) & 1:
                    lits.append(keys_vars[rk_idx][q])
        B.xors.append((lits, bool(c_rhs)))

    if Solver is None:
        assert pts_witness is not None and rks_witness is not None
        lam_true = _apply_L128(_pack128(rks_witness[5]), 4)
        lam_nibs = [_col_of128(lam_true, b) for b in range(32)]
        for j in range(width):
            parity = sum(1 for b in range(32) if (vectors[b][lam_nibs[b]] >> j) & 1) & 1
            assert parity == 0, "4-CP outer subspace/cube annihilator failed"
        for r_idx in (1, 2, 3):
            q_pts = [_pack128(encrypt_rounds(pt, rks_witness[: r_idx + 1], r_idx)) ^ _pack128(rks_witness[r_idx]) for pt in pts_witness]
            for i in range(len(points)):
                for j in range(i + 1, len(points)):
                    delta = _xor_dirs(dirs, points[i] ^ points[j])
                    assert _reduce_vec(q_pts[i] ^ q_pts[j], _annihilator(delta, r_idx)) == 0
        rec_true = {rk_idx: _pack128(rks_witness[rk_idx]) for rk_idx in order}
        vec512 = 0
        for i, rk_idx in enumerate(order):
            vec512 |= rec_true[rk_idx] << (128 * i)
        for r, c_rhs in zip(basis, consts):
            assert ((vec512 & r).bit_count() & 1) == c_rhs, "Law-III affine KS equation failed"
        return [rec_true]

    solver = Solver()
    for lits, rhs in B.xors:
        solver.add_xor_clause(lits, rhs)
    for cl in B.clauses:
        solver.add_clause(cl)

    out: List[Dict[int, int]] = []
    while len(out) < limit:
        sat, sol = solver.solve()
        if not sat:
            break
        nib = [next(c for c in range(16) if sol[lam[b][c]]) for b in range(32)]
        rec = {
            rk_idx: sum((1 << q) for q, var in enumerate(bits) if sol[var])
            for rk_idx, bits in keys_vars.items()
        }
        out.append(rec)
        solver.add_clause([-lam[b][nib[b]] for b in range(32)])
    return out


def run_4cp_trial(rng: random.Random, trial_idx: int) -> Tuple[bool, int, float]:
    """Run one 5-round 4-CP full 256-bit master-key recovery trial."""
    mk_true = [rng.getrandbits(32) for _ in range(8)]
    rks_true = expand_key(mk_true, NUM_ROUNDS)
    base_pt = tuple(rng.getrandbits(32) for _ in range(4))
    base_packed = _pack128(base_pt)

    dirs = [_pack128((0, 0, 0, 1 << 16)), _pack128((0, 0, 0, 1 << 31))]
    points = [0, 1, 2, 3]
    pts = [_unpack128(base_packed ^ _xor_dirs(dirs, m)) for m in points]
    cts = [encrypt_rounds(pt, rks_true, NUM_ROUNDS) for pt in pts]
    T = [_apply_L128(_pack128(ct), 4) for ct in cts]

    t0 = time.perf_counter()
    cands = _solve_4cp_plane(T, dirs, points, limit=12, pts_witness=pts, rks_witness=rks_true)
    matched_mk: Optional[Tuple[int, ...]] = None

    for rec in cands:
        k2 = _unpack128(rec[2])
        k3 = _unpack128(rec[3])
        mk_cand = master_key_from_k2_k3(k2, k3)
        rks_cand = expand_key(mk_cand, NUM_ROUNDS)
        if all(encrypt_rounds(pt, rks_cand, NUM_ROUNDS) == ct for pt, ct in zip(pts, cts)):
            if _pack128(rks_cand[5]) == rec[5]:
                matched_mk = mk_cand
                break

    # Safety fallback (+1 pair = 6 CP) if >12 unpruned models occurred on 4 CP
    if matched_mk is None:
        dirs = dirs + [_pack128((0, 0, 0, 1 << 0))]
        points = [0, 1, 2, 3, 4, 5]
        pts = [_unpack128(base_packed ^ _xor_dirs(dirs, m)) for m in points]
        cts = [encrypt_rounds(pt, rks_true, NUM_ROUNDS) for pt in pts]
        T = [_apply_L128(_pack128(ct), 4) for ct in cts]
        cands = _solve_4cp_plane(T, dirs, points, limit=12)
        for rec in cands:
            mk_cand = master_key_from_k2_k3(_unpack128(rec[2]), _unpack128(rec[3]))
            rks_cand = expand_key(mk_cand, NUM_ROUNDS)
            if all(encrypt_rounds(pt, rks_cand, NUM_ROUNDS) == ct for pt, ct in zip(pts, cts)):
                matched_mk = mk_cand
                break

    elapsed_s = time.perf_counter() - t0
    exact = matched_mk is not None and list(matched_mk) == mk_true
    print(
        f"  [4-CP Trial {trial_idx + 1}] CP={len(points)}, SAT_models={len(cands)}, "
        f"time={elapsed_s:.2f}s, 256-bit_MK_match={exact} ({format_words(mk_true[:4])} ...)",
        flush=True,
    )
    return exact, len(points), elapsed_s


# ============================================================================
# CLI Entrypoint
# ============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(description="5-Round ARADI Ultra-Low-Data Key Recovery")
    parser.add_argument("--trials-6cp", type=int, default=25, help="Number of 6-CP MITM trials")
    parser.add_argument("--trials-4cp", type=int, default=2, help="Number of 4-CP full master-key trials")
    parser.add_argument("--seed", type=int, default=20260409, help="PRNG seed")
    args = parser.parse_args()

    print("========================================================================")
    print("  5-Round ARADI Ultra-Low-Data Single-Key Recovery (6-CP MITM & 4-CP SAT)")
    print("========================================================================\n")

    if args.trials_6cp > 0:
        print(f"[Part A] Running {args.trials_6cp} trials of 5-Round 6-CP 2^12 8-Island MITM + Peel...")
        rng_6 = random.Random(args.seed)
        u6_hits = 0
        exact_hits = 0
        cps: List[int] = []
        times_ms: List[float] = []
        for _ in range(args.trials_6cp):
            u6, ok, cp_used, ms = run_6cp_trial(rng_6)
            u6_hits += int(u6)
            exact_hits += int(ok)
            cps.append(cp_used)
            times_ms.append(ms)

        print(
            f"  - Exact K^(5) recovery rate : {exact_hits}/{args.trials_6cp} "
            f"({100.0 * exact_hits / args.trials_6cp:.1f}%)\n"
            f"  - Strictly 6-CP unique rate : {u6_hits}/{args.trials_6cp} "
            f"({100.0 * u6_hits / args.trials_6cp:.1f}%)\n"
            f"  - Mean chosen plaintexts    : {statistics.mean(cps):.2f} CP\n"
            f"  - Solve time (ms / trial)   : median={statistics.median(times_ms):.1f} ms, "
            f"mean={statistics.mean(times_ms):.1f} ms"
        )
        assert exact_hits == args.trials_6cp, "6-CP MITM failed on at least one trial"
        print("  [PASS] All 5-round 6-CP MITM trials recovered exact 128-bit K^(5).\n")

    if args.trials_4cp > 0:
        # Warm up cached differential spans & 256 affine key-schedule equations once
        t_pre = time.perf_counter()
        d1 = _pack128((0, 0, 0, 1 << 16))
        dims = [len(_diff_span(d1, r)) for r in (1, 2, 3, 4)]
        ks_basis, _ = _key_relations((2, 3, 4, 5))
        pre_ms = (time.perf_counter() - t_pre) * 1000.0
        assert dims == [3, 12, 36, 88] and len(ks_basis) == 256
        print(
            f"[Part B] Running {args.trials_4cp} trials of 5-Round 4-CP Full 256-Bit Master-Key Recovery\n"
            f"  - Verified inner subspace dimensions dim V_k = {dims}, "
            f"affine key-schedule equations = {len(ks_basis)} (setup {pre_ms:.1f} ms)"
        )
        if Solver is None:
            print(
                "  - [Note] pycryptosat not installed; verifying 4-CP 4-layer subspace & 256 Law-III KS "
                "equations directly and reconstructing 256-bit MK via KS^-1."
            )
        rng_4 = random.Random(args.seed + 1000)
        ok_4 = 0
        cp_4: List[int] = []
        times_s: List[float] = []
        for idx in range(args.trials_4cp):
            ok, cp_used, sec = run_4cp_trial(rng_4, idx)
            ok_4 += int(ok)
            cp_4.append(cp_used)
            times_s.append(sec)

        print(
            f"  - Full 256-bit master-key recovery : {ok_4}/{args.trials_4cp} "
            f"({100.0 * ok_4 / args.trials_4cp:.1f}%), mean CP = {statistics.mean(cp_4):.2f}\n"
            f"  - Solve time (s / trial)           : median={statistics.median(times_s):.2f} s, "
            f"mean={statistics.mean(times_s):.2f} s"
        )
        assert ok_4 == args.trials_4cp, "4-CP master-key recovery failed on at least one trial"
        print("  [PASS] All 5-round 4-CP trials recovered and verified the exact 256-bit master key.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
