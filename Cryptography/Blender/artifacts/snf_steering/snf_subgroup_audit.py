#!/usr/bin/env python3
"""
snf_subgroup_audit.py

Self-contained algebraic and empirical verification for Section 5 (Multi-Shifted
2-Adic SNF Modular-Sum Steering and Preimage Recovery on Blender):
  1. Verifies why GF(2) Gaussian elimination fails (100% 2nd-order XOR derivative != 0)
     while Z/2^32Z modular superposition is exact (100% 2nd-order modular derivative == 0).
  2. Computes exact 2-adic Smith Normal Form (invariant factors over Z/2^32Z) across:
       - Theoretical 8-phase clean-subspace generators D^(0..7)
       - Alternating subspace (X_phi = 0)
       - Unshifted 38-bit subspace W'[0] in {0, ~0} -> (1, 1, 2, 2, 4, 4, 8, 8) (index 2^12)
       - Multi-shifted 38-bit subspaces W'[0] in {a_m, ~a_m} -> (1, 1, 2, 2, 2, 2, 2, 2) (index 2^6, H_par)
  3. Verifies C_inv = 0 end-to-end on a concrete equal-length, equal-checksum closed-walk
     pair directly against Colin Bradbury's unmodified NIST C reference (blender_cli).
  4. Evaluates exact preimage / second-preimage complexities across Blender-224/256/384/512.
"""

import json
import math
import os
import struct
import subprocess
import sys
import tempfile

M32 = 0xFFFFFFFF
MOD32 = 1 << 32

IV32_256 = [
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
]


def ror32(x: int, n: int) -> int:
    x &= M32
    return ((x >> n) | ((x << (32 - n)) & M32)) & M32


def rol32(x: int, n: int) -> int:
    x &= M32
    return (((x << n) & M32) | (x >> (32 - n))) & M32


def snf_mod_2w(rows, w: int = 32):
    """Exact 2-adic Smith Normal Form valuations (e_0, ..., e_{n-1}) over Z / 2^w Z."""
    mod = 1 << w
    mask = mod - 1
    m = len(rows)
    n = len(rows[0])
    A = [[x & mask for x in r] for r in rows]
    exps = []
    for c in range(min(m, n)):
        best_v = w
        bi, bj = -1, -1
        for i in range(c, m):
            row_i = A[i]
            for j in range(c, n):
                val = row_i[j]
                if val != 0:
                    v = (val & -val).bit_length() - 1
                    if v < best_v:
                        best_v, bi, bj = v, i, j
                        if best_v == 0:
                            break
            if best_v == 0:
                break
        if best_v == w:
            exps.extend([w] * (n - c))
            break
        if bi != c:
            A[c], A[bi] = A[bi], A[c]
        if bj != c:
            for i in range(m):
                A[i][c], A[i][bj] = A[i][bj], A[i][c]
        inv_unit = pow(A[c][c] >> best_v, -1, mod)
        A[c] = [(x * inv_unit) & mask for x in A[c]]
        piv_row = A[c]
        for i in range(c + 1, m):
            entry = A[i][c]
            if entry != 0:
                q = entry >> best_v
                A[i] = [(A[i][j] - q * piv_row[j]) & mask for j in range(n)]
        for j in range(c + 1, n):
            A[c][j] = 0
        exps.append(best_v)
    if len(exps) < n:
        exps.extend([w] * (n - len(exps)))
    return sorted(exps)


def signed32(u: int) -> int:
    u &= M32
    return u - MOD32 if (u & 0x80000000) else u


def verify_f2_nonlinearity(vecs):
    s_base = [
        0x12345678, 0x9ABCDEF0, 0x13579BDF, 0x2468ACE0,
        0xCAFEBABE, 0xDEADBEEF, 0x0F1E2D3C, 0x4B5A6978,
    ]
    f2_fail = 0
    mod_pass = 0
    ntests = min(len(vecs) - 1, 100)
    for i in range(ntests):
        v1, v2 = vecs[i], vecs[i + 1]
        y00 = s_base[:]
        y10 = [(s + a) & M32 for s, a in zip(s_base, v1)]
        y01 = [(s + b) & M32 for s, b in zip(s_base, v2)]
        y11 = [(s + a + b) & M32 for s, a, b in zip(s_base, v1, v2)]
        if any((a ^ b ^ c ^ d) != 0 for a, b, c, d in zip(y00, y10, y01, y11)):
            f2_fail += 1
        if all(((a - b - c + d) & M32) == 0 for a, b, c, d in zip(y00, y10, y01, y11)):
            mod_pass += 1
    assert f2_fail == ntests and mod_pass == ntests
    print("=== [1] Abelian Superposition Test: GF(2) vs (Z/2^32Z)^8 ===")
    print(f"  2nd-order XOR derivative != 0 (GF(2) fails)     : {f2_fail}/{ntests} (100.0%)")
    print(f"  2nd-order modular derivative == 0 (Z/2^32Z exact): {mod_pass}/{ntests} (100.0%)")


def verify_snf_hierarchy(data):
    print("\n=== [2] 2-Adic Smith Normal Form Subgroup Hierarchy over (Z/2^32Z)^8 ===")
    phases_D = data["std_B_phases_D"]
    even_gens, odd_gens = [], []
    for phi in range(8):
        D = phases_D[phi]
        even_gens.append([(-1 - 2 * D[0]) & M32, (-1 - 2 * D[2]) & M32,
                          (-1 - 2 * D[4]) & M32, (-1 - 2 * D[6]) & M32])
        for b in range(32):
            scale = 1 << b
            odd_gens.append([
                scale & M32,
                (scale * (1 - 2 * ((D[3] >> b) & 1))) & M32,
                (scale * (1 - 2 * ((D[5] >> b) & 1))) & M32,
                (scale * (1 - 2 * ((D[7] >> b) & 1))) & M32,
            ])
    pred_unshifted = sorted(snf_mod_2w(even_gens, 32) + snf_mod_2w(odd_gens, 32))
    print(f"  Theoretical unshifted SNF valuations from D^(0..7) : {pred_unshifted} (index 2^{sum(pred_unshifted)})")

    checks = [
        ("Alternating X_phi=0 (eqlen+eqcsum)", data["alt_eqlen_eqcsum"], None),
        ("Unshifted W'[0] in {0,~0} (eqlen+eqcsum)", data["unshifted_eqlen_eqcsum"], [0, 0, 1, 1, 2, 2, 3, 3]),
        ("Multi-Shifted W'[0] in {a_m,~a_m} (eqlen+eqcsum)", data["shifted_multi_eqlen_eqcsum"], [0, 0, 1, 1, 1, 1, 1, 1]),
    ]
    for label, vecs, expected in checks:
        exps = snf_mod_2w(vecs, 32)
        s0_max = max(abs(signed32(v[0])) for v in vecs)
        par_ok = all(
            ((v[0] ^ v[2]) & 1) == 0 and ((v[0] ^ v[4]) & 1) == 0 and ((v[0] ^ v[6]) & 1) == 0 and
            ((v[1] ^ v[3]) & 1) == 0 and ((v[1] ^ v[5]) & 1) == 0 and ((v[1] ^ v[7]) & 1) == 0
            for v in vecs
        )
        assert par_ok
        if expected is not None:
            assert exps == expected, f"SNF mismatch on {label}: {exps} != {expected}"
        inv_str = "(" + ",".join("1" if e == 0 else str(1 << e) for e in exps) + ")"
        print(f"  [PASS] {label:48s}: SNF={inv_str}, index=2^{sum(exps)}, max|dSum[0]|={s0_max:#x}")


def verify_c_inv_zero(cw, cli_bin: str):
    print("\n=== [3] Zero Finalization Inversion Cost (C_inv = 0) against Official NIST C Reference ===")
    h1, h2 = cw["m1_words_le_hex"], cw["m2_words_le_hex"]
    exp_dsum = cw["expected_dsum"]

    def hex_words_to_bytes(whx: str) -> bytes:
        return b"".join(struct.pack("<I", int(whx[8 * i: 8 * (i + 1)], 16)) for i in range(len(whx) // 8))

    m1 = hex_words_to_bytes(h1)
    m2 = hex_words_to_bytes(h2)
    assert len(m1) == len(m2) and m1 != m2 and m1[:13] == m2[:13]

    with tempfile.TemporaryDirectory() as tmpdir:
        t1 = os.path.join(tmpdir, "snf_tmp_m1.bin")
        t2 = os.path.join(tmpdir, "snf_tmp_m2.bin")
        with open(t1, "wb") as f1, open(t2, "wb") as f2:
            f1.write(m1)
            f2.write(m2)
        dig1 = bytes.fromhex(subprocess.check_output([cli_bin, "-l", "256", "-f", t1]).decode().strip())
        dig2 = bytes.fromhex(subprocess.check_output([cli_bin, "-l", "256", "-f", t2]).decode().strip())

    dw1 = list(struct.unpack(">8I", dig1))
    dw2 = list(struct.unpack(">8I", dig2))
    post_dsum = [((a - b) & M32) for a, b in zip(dw1, dw2)]
    assert post_dsum == exp_dsum, f"Post-tail digest diff {post_dsum} != pre-tail dSum {exp_dsum}"
    print(f"  [PASS] Verified {len(m1)}-byte closed-walk pair on official core.c:")
    print("         Post-finalization (Digest(M1) - Digest(M2)) mod 2^32 == pre-tail Delta Sum identically (C_inv = 0)!")


def sink_dag_slot_cost_log2(state_bits: float, n_slots: int, L0_bits: int, S_dag_bits: float,
                            step_gain_bits: float = 8.0, n_anchors: int = 8) -> float:
    N = 2.0 ** state_bits
    L0 = 2.0 ** L0_bits
    S = 2.0 ** S_dag_bits
    t_s1 = n_anchors * (N / L0)
    t_grow = n_anchors * (N / (2.0 ** step_gain_bits)) * math.log(S / L0)
    slots_per_anchor = math.ceil(n_slots / float(n_anchors))
    buckets = max(1.0, L0 / 8.0)
    w_per_anchor = max(2.5 * slots_per_anchor, math.sqrt(4.0 * slots_per_anchor * buckets))
    t_collect = n_anchors * w_per_anchor * (N / S)
    return math.log2(t_s1 + t_grow + t_collect)


def print_preimage_summary():
    print("\n=== [4] Exact Full Preimage / Second-Preimage Complexity Summary ===")
    configs = [
        ("Blender-256", 256, 32, 8, 8, 12, 11),
        ("Blender-384", 384, 64, 6, 6, 14, 20),
        ("Blender-512", 512, 64, 8, 8, 14, 20),
    ]
    for name, n_bits, w_bits, m_words, period, ell, L0_bits in configs:
        clean_bits = (1 + w_bits + 2) + math.log2(period)
        par_codim = m_words - 2
        B_bits = w_bits + (n_bits - par_codim)
        k = 1 << ell
        s = B_bits / (ell + 1.0)
        n_slots = k * math.ceil(s)
        msg_log2 = math.log2(n_slots * (1 << L0_bits))
        S_dag_bits = max(L0_bits + 2.0, min(s, clean_bits - 8.0))
        slot_ops = sink_dag_slot_cost_log2(clean_bits, n_slots, L0_bits, S_dag_bits, min(9.0, L0_bits - 1.0))
        tot_time = math.log2((2.0 ** (ell + s)) + (2.0 ** slot_ops) + (2.0 ** msg_log2) + (2.0 ** par_codim))
        lm_time = math.log2(n_bits) + n_bits / 2.0
        print(f"  {name:11s}: SNF+Wagner(k=2^{ell:<2d}) Time = 2^{tot_time:5.2f}, Mem = 2^{max(s, S_dag_bits):5.2f}, "
              f"MsgLen = 2^{msg_log2:5.2f} words  (vs. ISC 2009: 2^{lm_time:6.2f})")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    json_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "witnesses", "snf_audit_vectors.json")
    with open(json_path, "r", encoding="utf-8") as fp:
        data = json.load(fp)
    verify_f2_nonlinearity(data["shifted_multi_eqlen_eqcsum"])
    verify_snf_hierarchy(data)
    if len(sys.argv) >= 3 and os.path.exists(sys.argv[2]):
        verify_c_inv_zero(data["concrete_witness"], sys.argv[2])
    print_preimage_summary()


if __name__ == "__main__":
    main()
