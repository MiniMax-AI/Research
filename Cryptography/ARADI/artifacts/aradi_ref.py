#!/usr/bin/env python3
"""Reference implementation and official test vectors for NSA ARADI.

Specification:
  Patricia Greene, Mark Motley, and Bryan Weeks,
  "ARADI and LLAMA: Low-Latency Cryptography for Memory Encryption",
  IACR Cryptology ePrint Archive, Report 2024/1240 (revised 2024-09-05),
  Sections 3.1-3.4 and 4.2 (Figures 3.4, 3.5, 4.3, 4.4).

Conventions matching ePrint 2024/1240 verbatim:
  - 128-bit state is a 4-tuple of 32-bit words (w, x, y, z), where w is the
    most significant word and z is the least significant word.
  - Column b in {0, ..., 31} consists of bit b ((word >> b) & 1) across
    (w, x, y, z), packed as w_b | (x_b << 1) | (y_b << 2) | (z_b << 3).
  - Round i in {0, ..., R - 1} computes Q_{i+1} = Lambda_{(r0+i) mod 4}(pi(Q_i ^ K^(i))),
    followed after round R - 1 by post-whitening C = Q_R ^ K^(R).
  - 256-bit master key is an 8-element sequence K[0..7] of 32-bit words.

Pure Python 3 standard library only; zero external dependencies.
"""
from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

MASK32 = 0xFFFFFFFF
MASK16 = 0xFFFF

# Shift parameters (a_j, b_j, c_j) for L_j, j = 0, 1, 2, 3 (2024/1240, p. 5 & Fig. 3.5).
ROT_CONSTANTS: List[Tuple[int, int, int]] = [
    (11, 8, 14),
    (10, 9, 11),
    (9, 4, 14),
    (8, 9, 7),
]

# Eight disjoint 4-column islands induced by L_3 = (8, 9, 7) since gcd(8, 16) = 8.
GROUPS: List[List[int]] = [
    [0, 8, 23, 31],
    [1, 9, 16, 24],
    [2, 10, 17, 25],
    [3, 11, 18, 26],
    [4, 12, 19, 27],
    [5, 13, 20, 28],
    [6, 14, 21, 29],
    [7, 15, 22, 30],
]
GROUP_OF = {col: gi for gi, grp in enumerate(GROUPS) for col in grp}


def rol32(x: int, n: int) -> int:
    """32-bit left circular shift S_32^n(x)."""
    n &= 31
    x &= MASK32
    return ((x << n) | (x >> (32 - n))) & MASK32 if n else x


def ror32(x: int, n: int) -> int:
    """32-bit right circular shift S_32^{-n}(x)."""
    return rol32(x, (32 - (n & 31)) & 31)


def rol16(x: int, n: int) -> int:
    """16-bit left circular shift S_16^n(x)."""
    n &= 15
    x &= MASK16
    return ((x << n) | (x >> (16 - n))) & MASK16 if n else x


def sbox(w: int, x: int, y: int, z: int) -> Tuple[int, int, int, int]:
    """Bit-sliced 4-bit ARADI S-box pi (4 sequential Toffoli gates)."""
    w &= MASK32
    x &= MASK32
    y &= MASK32
    z &= MASK32
    x ^= w & y
    z ^= x & y
    y ^= w & z
    w ^= x & z
    return w & MASK32, x & MASK32, y & MASK32, z & MASK32


def sbox_inv(w: int, x: int, y: int, z: int) -> Tuple[int, int, int, int]:
    """Bit-sliced inverse S-box pi^{-1} (reverse Toffoli order)."""
    w &= MASK32
    x &= MASK32
    y &= MASK32
    z &= MASK32
    w ^= x & z
    y ^= w & z
    z ^= x & y
    x ^= w & y
    return w & MASK32, x & MASK32, y & MASK32, z & MASK32


SBOX_TABLE: List[int] = [0] * 16
SBOX_INV_TABLE: List[int] = [0] * 16
for _u in range(16):
    _ow, _ox, _oy, _oz = sbox(_u & 1, (_u >> 1) & 1, (_u >> 2) & 1, (_u >> 3) & 1)
    _v = _ow | (_ox << 1) | (_oy << 2) | (_oz << 3)
    SBOX_TABLE[_u] = _v
    SBOX_INV_TABLE[_v] = _u


def linear_L(word: int, a: int, b: int, c: int) -> int:
    """32-bit word linear map L(a, b, c, u || l) (involutive when 2a = b + c mod 16)."""
    u = (word >> 16) & MASK16
    l = word & MASK16
    nu = u ^ rol16(u, a) ^ rol16(l, c)
    nl = l ^ rol16(l, a) ^ rol16(u, b)
    return ((nu << 16) | nl) & MASK32


def apply_linear_state(
    state: Tuple[int, int, int, int],
    layer_idx: int,
    rot_constants: Sequence[Tuple[int, int, int]] = ROT_CONSTANTS,
) -> Tuple[int, int, int, int]:
    """Apply Lambda_{layer_idx mod 4} word-wise to (w, x, y, z)."""
    a, b, c = rot_constants[layer_idx & 3]
    return (
        linear_L(state[0], a, b, c),
        linear_L(state[1], a, b, c),
        linear_L(state[2], a, b, c),
        linear_L(state[3], a, b, c),
    )


def key_M(ri: int, rj: int, X: int, Y: int) -> Tuple[int, int]:
    """Key-schedule pair mixing M(i, j, X, Y) from Fig. 3.5."""
    rY = rol32(Y, ri)
    return (rY ^ rol32(X, rj) ^ X) & MASK32, (rY ^ X) & MASK32


def key_M_inv(ri: int, rj: int, A: int, B: int) -> Tuple[int, int]:
    """Exact inverse of key_M(ri, rj, X, Y) -> (A, B)."""
    X = ror32((A ^ B) & MASK32, rj)
    Y = ror32((B ^ X) & MASK32, ri)
    return X, Y


def expand_key(master_key: Sequence[int], num_rounds: int = 16) -> List[Tuple[int, int, int, int]]:
    """Expand 8-word master key K[0..7] into num_rounds + 1 128-bit round keys.

    Returns [K^(0), K^(1), ..., K^(num_rounds)], where K^(num_rounds) is the
    post-whitening key drawn from half j = num_rounds mod 2.
    """
    if len(master_key) != 8:
        raise ValueError(f"Expected 8 32-bit key words, got {len(master_key)}")
    K = [int(k) & MASK32 for k in master_key]
    round_keys: List[Tuple[int, int, int, int]] = []
    for i in range(num_rounds):
        j = i & 1
        base = 4 * j
        round_keys.append((K[base], K[base + 1], K[base + 2], K[base + 3]))
        K[1], K[0] = key_M(1, 3, K[1], K[0])
        K[3], K[2] = key_M(9, 28, K[3], K[2])
        K[5], K[4] = key_M(1, 3, K[5], K[4])
        K[7], K[6] = key_M(9, 28, K[7], K[6])
        K[7] = (K[7] ^ i) & MASK32
        if j == 0:
            K[1], K[2] = K[2], K[1]
            K[5], K[6] = K[6], K[5]
        else:
            K[1], K[4] = K[4], K[1]
            K[3], K[6] = K[6], K[3]
    j = num_rounds & 1
    base = 4 * j
    round_keys.append((K[base], K[base + 1], K[base + 2], K[base + 3]))
    return round_keys


def inv_key_update(state_after: Sequence[int], step_idx: int) -> List[int]:
    """Invert step `step_idx` of the ARADI key-schedule register update."""
    K = [int(w) & MASK32 for w in state_after]
    j = step_idx & 1
    if j == 0:
        K[1], K[2] = K[2], K[1]
        K[5], K[6] = K[6], K[5]
    else:
        K[1], K[4] = K[4], K[1]
        K[3], K[6] = K[6], K[3]
    K[7] = (K[7] ^ step_idx) & MASK32
    K[1], K[0] = key_M_inv(1, 3, K[1], K[0])
    K[3], K[2] = key_M_inv(9, 28, K[3], K[2])
    K[5], K[4] = key_M_inv(1, 3, K[5], K[4])
    K[7], K[6] = key_M_inv(9, 28, K[7], K[6])
    return K


def master_key_from_even_pair(
    k_even: Sequence[int],
    k_odd: Sequence[int],
    even_idx: int,
) -> Tuple[int, int, int, int, int, int, int, int]:
    """Reconstruct the 256-bit master key K[0..7] in O(1) from (K^(2m), K^(2m+1)).

    At even step m = even_idx (j = 0):
      - K^(m) = state_m[0..3]
      - Step m updates the upper half (s4, s5, s6, s7) independently via M_0, M_1,
        XORs m into word 7, and swaps words (5, 6), so K^(m+1) = state_{m+1}[4..7]
        uniquely determines state_m[4..7].
      - Inverting steps m - 1, ..., 0 recovers K[0..7].
    """
    if even_idx < 0 or (even_idx & 1) != 0:
        raise ValueError(f"even_idx must be a non-negative even integer, got {even_idx}")
    u4, u6, u5, u7 = (int(w) & MASK32 for w in k_odd)
    s5, s4 = key_M_inv(1, 3, u5, u4)
    s7, s6 = key_M_inv(9, 28, (u7 ^ even_idx) & MASK32, u6)
    state = [int(w) & MASK32 for w in k_even] + [s4, s5, s6, s7]
    for step in range(even_idx - 1, -1, -1):
        state = inv_key_update(state, step)
    return tuple(state)  # type: ignore[return-value]


def master_key_from_k2_k3(
    k2: Sequence[int],
    k3: Sequence[int],
) -> Tuple[int, int, int, int, int, int, int, int]:
    """Closed-form O(1) master-key recovery from (K^(2), K^(3)) (even_idx = 2)."""
    return master_key_from_even_pair(k2, k3, 2)


def master_key_from_k4_k5(
    k4: Sequence[int],
    k5: Sequence[int],
) -> Tuple[int, int, int, int, int, int, int, int]:
    """Closed-form O(1) master-key recovery from (K^(4), K^(5)) (even_idx = 4)."""
    return master_key_from_even_pair(k4, k5, 4)


def encrypt_rounds(
    pt: Tuple[int, int, int, int],
    round_keys: Sequence[Tuple[int, int, int, int]],
    num_rounds: int = 16,
    r0: int = 0,
    rot_constants: Sequence[Tuple[int, int, int]] = ROT_CONSTANTS,
) -> Tuple[int, int, int, int]:
    """Encrypt a 128-bit block across `num_rounds` rounds starting at offset `r0`.

    `round_keys` must contain at least `num_rounds + 1` 128-bit keys:
    `round_keys[0 .. num_rounds - 1]` are applied before the S-box in steps
    `0 .. num_rounds - 1` (with linear layer `(r0 + step) mod 4`), and
    `round_keys[num_rounds]` is XORed as post-whitening.
    """
    if len(round_keys) < num_rounds + 1:
        raise ValueError(
            f"Need at least {num_rounds + 1} round keys, got {len(round_keys)}"
        )
    w, x, y, z = (int(v) & MASK32 for v in pt)
    for step in range(num_rounds):
        kw, kx, ky, kz = round_keys[step]
        w, x, y, z = sbox(w ^ kw, x ^ kx, y ^ ky, z ^ kz)
        a, b, c = rot_constants[(r0 + step) & 3]
        w = linear_L(w, a, b, c)
        x = linear_L(x, a, b, c)
        y = linear_L(y, a, b, c)
        z = linear_L(z, a, b, c)
    kw, kx, ky, kz = round_keys[num_rounds]
    return (w ^ kw) & MASK32, (x ^ kx) & MASK32, (y ^ ky) & MASK32, (z ^ kz) & MASK32


def decrypt_rounds(
    ct: Tuple[int, int, int, int],
    round_keys: Sequence[Tuple[int, int, int, int]],
    num_rounds: int = 16,
    r0: int = 0,
    rot_constants: Sequence[Tuple[int, int, int]] = ROT_CONSTANTS,
) -> Tuple[int, int, int, int]:
    """Exact inverse of `encrypt_rounds` (uses L_j^{-1} = L_j when involutive)."""
    if len(round_keys) < num_rounds + 1:
        raise ValueError(
            f"Need at least {num_rounds + 1} round keys, got {len(round_keys)}"
        )
    kw, kx, ky, kz = round_keys[num_rounds]
    w, x, y, z = (
        (ct[0] ^ kw) & MASK32,
        (ct[1] ^ kx) & MASK32,
        (ct[2] ^ ky) & MASK32,
        (ct[3] ^ kz) & MASK32,
    )
    for step in range(num_rounds - 1, -1, -1):
        a, b, c = rot_constants[(r0 + step) & 3]
        w = linear_L(w, a, b, c)
        x = linear_L(x, a, b, c)
        y = linear_L(y, a, b, c)
        z = linear_L(z, a, b, c)
        w, x, y, z = sbox_inv(w, x, y, z)
        kw, kx, ky, kz = round_keys[step]
        w = (w ^ kw) & MASK32
        x = (x ^ kx) & MASK32
        y = (y ^ ky) & MASK32
        z = (z ^ kz) & MASK32
    return w, x, y, z


def bytes_to_words(data: bytes) -> Tuple[int, int, int, int]:
    """Unpack 16 big-endian bytes into (w, x, y, z)."""
    if len(data) != 16:
        raise ValueError(f"Expected 16 bytes, got {len(data)}")
    return tuple(int.from_bytes(data[4 * i : 4 * i + 4], "big") for i in range(4))  # type: ignore[return-value]


def words_to_bytes(words: Iterable[int]) -> bytes:
    """Pack 32-bit words into big-endian bytes."""
    return b"".join((int(w) & MASK32).to_bytes(4, "big") for w in words)


def format_words(words: Iterable[int]) -> str:
    """Format 32-bit words in lowercase 8-hex-digit blocks matching 2024/1240."""
    return " ".join(f"{int(w) & MASK32:08x}" for w in words)


# ============================================================================
# Official NSA Test Vectors from ePrint 2024/1240
# ============================================================================

# Figure 3.4 (p. 7): K[0..7], P = (0, 0, 0, 0), all 17 round keys, all 16
# intermediate round states (after Lambda_{i mod 4}), and final ciphertext.
FIG_3_4_KEY: List[int] = [
    0x03020100,
    0x07060504,
    0x0B0A0908,
    0x0F0E0D0C,
    0x13121110,
    0x17161514,
    0x1B1A1918,
    0x1F1E1D1C,
]
FIG_3_4_PT: Tuple[int, int, int, int] = (0x00000000, 0x00000000, 0x00000000, 0x00000000)
FIG_3_4_CT: Tuple[int, int, int, int] = (0x3F09ABF4, 0x00E3BD74, 0x03260DEF, 0xB7C53912)

FIG_3_4_ROUND_KEYS: List[Tuple[int, int, int, int]] = [
    (0x03020100, 0x07060504, 0x0B0A0908, 0x0F0E0D0C),
    (0x31323734, 0x2B2C2D2A, 0x89829F94, 0xEADDCCFB),
    (0x19181312, 0x49484342, 0xBFB2B5B8, 0xEFE2E5E8),
    (0x93D8DD96, 0x49BBF102, 0x12918D0E, 0x2CAF0292),
    (0x7C795E5B, 0x6E0A4A2F, 0x708952AB, 0x0FB51EB7),
    (0x73BE37F3, 0xB12DE15C, 0x6D10261A, 0x63FA1FB1),
    (0x30E1A565, 0x56518EBA, 0x38A4DC70, 0x43B62B6B),
    (0x6FF94BF4, 0xA1525D49, 0x960D690A, 0xF40AC5E6),
    (0x652B43FA, 0x7EA0CAA1, 0x8356ECA6, 0xEED8D0CA),
    (0x1E8816B8, 0xEAF40402, 0xBF1911DB, 0xD2ED83C3),
    (0x2AED0767, 0xD7E42972, 0x0DDCAC43, 0xE0CE34BD),
    (0xE587DB6F, 0xD93A728E, 0xE7A79043, 0x54E47C4C),
    (0x5DEAFDDF, 0x1235C451, 0xB9420597, 0x1BC4FB83),
    (0xF95881FC, 0xA9CBAE8E, 0x266A00C2, 0x64230546),
    (0xCC0FAB2E, 0x5B7AAD77, 0x32495539, 0xB022810A),
    (0x71C5C046, 0x8AB9AA02, 0xD8FB0856, 0xB7DFA119),
    (0xA443053B, 0x69322A8E, 0xE8ABFB4F, 0x41CF0CA8),
]

FIG_3_4_ROUND_STATES: List[Tuple[int, int, int, int]] = [
    (0x367F232B, 0x25252020, 0x4A4A4040, 0x7C35636B),
    (0xEE64F20F, 0xF9BCE360, 0x418D0976, 0x1042F571),
    (0xE65DA996, 0x564E30AA, 0x8EBFFAD6, 0x2CFEA43D),
    (0x19FB2D3E, 0x2EA0FF0A, 0x5F80E087, 0xEAB056A4),
    (0x48FF5CCA, 0x5747215C, 0x587A96C3, 0x5C895983),
    (0x99F5DB6E, 0xB376D237, 0x35C04785, 0x11C1FBE7),
    (0x05DD6B05, 0xBA589C3F, 0x9705656E, 0xC46926D9),
    (0xA4A55EF4, 0x9A71C3E1, 0x239F293B, 0xC7AB0EBA),
    (0xF2CA9329, 0xAC68354A, 0xCBA990DC, 0xEFEC06A6),
    (0x2B4F661F, 0x1F94AECD, 0x8572FAE6, 0x79CCB74A),
    (0xEC3A6302, 0x9CA4753C, 0x91C92F12, 0xA0FF38F9),
    (0x4205949D, 0x0E2828C7, 0xBBA29CDE, 0x7BF46C7F),
    (0x7EA3E1A5, 0x4F7FAFE6, 0x6673F583, 0xE469266B),
    (0x27E6107F, 0x1A3E9F60, 0xE6F1261C, 0xAD5374A4),
    (0xC621BE33, 0xD8AA33DC, 0xCF025FB6, 0x93C87CDA),
    (0x9B4AAECF, 0x69D197FA, 0xEB8DF6A0, 0xF60A35BA),
]

# Figure 4.3 (p. 15): 10 LMAC counter-block encryptions P = (0, n, 0, 0), n = 0..9 under Key A.
FIG_4_3_KEY_A: List[int] = list(FIG_3_4_KEY)
FIG_4_3_VECTORS: List[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]] = [
    ((0, 0, 0, 0), (0x3F09ABF4, 0x00E3BD74, 0x03260DEF, 0xB7C53912)),
    ((0, 1, 0, 0), (0x27BFA00D, 0x83349615, 0x3B7095FD, 0xE633EFB8)),
    ((0, 2, 0, 0), (0x4CC1C23C, 0x04A204A4, 0xC514133E, 0x329286F8)),
    ((0, 3, 0, 0), (0x5C858E64, 0x45D745AF, 0xDAC0923E, 0x7B2D3272)),
    ((0, 4, 0, 0), (0x66E19501, 0xB781A2BE, 0x44C06864, 0x1BFCCBE0)),
    ((0, 5, 0, 0), (0xD7000BB1, 0xB3690D84, 0xCB846C76, 0xB4A1D8B3)),
    ((0, 6, 0, 0), (0x8B633523, 0xD66CC10A, 0xA6F166CD, 0xE45D4704)),
    ((0, 7, 0, 0), (0x92D7052D, 0xA5E13D30, 0x26DB3C87, 0x0F53B391)),
    ((0, 8, 0, 0), (0xFEBD86A9, 0x60E03CE7, 0x5199E241, 0x84747297)),
    ((0, 9, 0, 0), (0xAD670931, 0xA473A30E, 0x9D272FAA, 0x37959E5D)),
]

# Figure 4.4 (p. 16): 10 LMAC counter-block encryptions P = (0, n, 0, 0), n = 0..9 under Key B.
FIG_4_4_KEY_B: List[int] = [
    0x2C2D2E2F,
    0x28292A2B,
    0x24252627,
    0x20212223,
    0x1C1D1E1F,
    0x18191A1B,
    0x14151617,
    0x10111213,
]
FIG_4_4_VECTORS: List[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]] = [
    ((0, 0, 0, 0), (0xC974550C, 0x11722BDE, 0x40E7D38A, 0xC2A94E16)),
    ((0, 1, 0, 0), (0xF3E34230, 0x647EA2AF, 0x6D3640AC, 0xCE145AA9)),
    ((0, 2, 0, 0), (0x3B72E799, 0xD991F1FD, 0x07B15CD9, 0xBC5AC5B6)),
    ((0, 3, 0, 0), (0x286CC3B6, 0x4C814DC0, 0x741E9CEB, 0x20517185)),
    ((0, 4, 0, 0), (0x14331840, 0x913663C2, 0xC2AD0EF1, 0xC017F4B0)),
    ((0, 5, 0, 0), (0xF5CC3C8B, 0x9861CFDD, 0x16903048, 0xA37BB4E3)),
    ((0, 6, 0, 0), (0xE0E21C2E, 0x4559F34C, 0x8097C345, 0x21547D6B)),
    ((0, 7, 0, 0), (0x5F3109AA, 0xDFA8FEBD, 0xA10FE784, 0xF9F41E73)),
    ((0, 8, 0, 0), (0x0A695D1E, 0x4A8BFDF4, 0x89F7F326, 0x7CB82365)),
    ((0, 9, 0, 0), (0xBA163957, 0x14959F3E, 0xA6550B33, 0xC910158A)),
]


def verify_official_vectors(verbose: bool = True) -> bool:
    """Verify all 21 official NSA vectors, 17 round keys, 16 round states, and key inversion."""
    rks_a = expand_key(FIG_3_4_KEY, 16)
    assert rks_a == FIG_3_4_ROUND_KEYS, "Fig. 3.4 round-key schedule mismatch"

    st = FIG_3_4_PT
    for i in range(16):
        kw, kx, ky, kz = rks_a[i]
        st = sbox(st[0] ^ kw, st[1] ^ kx, st[2] ^ ky, st[3] ^ kz)
        st = apply_linear_state(st, i)
        assert st == FIG_3_4_ROUND_STATES[i], f"Fig. 3.4 intermediate state mismatch at round {i}"

    ct_3_4 = encrypt_rounds(FIG_3_4_PT, rks_a, 16, r0=0)
    assert ct_3_4 == FIG_3_4_CT, "Fig. 3.4 final ciphertext mismatch"
    assert decrypt_rounds(ct_3_4, rks_a, 16, r0=0) == FIG_3_4_PT, "Fig. 3.4 decrypt mismatch"

    for idx, (pt, expected_ct) in enumerate(FIG_4_3_VECTORS):
        got = encrypt_rounds(pt, rks_a, 16, r0=0)
        assert got == expected_ct, f"Fig. 4.3 vector {idx} mismatch: {got} != {expected_ct}"
        assert decrypt_rounds(got, rks_a, 16, r0=0) == pt

    rks_b = expand_key(FIG_4_4_KEY_B, 16)
    for idx, (pt, expected_ct) in enumerate(FIG_4_4_VECTORS):
        got = encrypt_rounds(pt, rks_b, 16, r0=0)
        assert got == expected_ct, f"Fig. 4.4 vector {idx} mismatch: {got} != {expected_ct}"
        assert decrypt_rounds(got, rks_b, 16, r0=0) == pt

    # Verify O(1) inverse key schedule across all even pairs (2m, 2m+1) for both keys.
    for mk, rks in ((FIG_3_4_KEY, rks_a), (FIG_4_4_KEY_B, rks_b)):
        for even_idx in range(0, 16, 2):
            rec = master_key_from_even_pair(rks[even_idx], rks[even_idx + 1], even_idx)
            assert list(rec) == mk, f"master_key_from_even_pair failed at even_idx={even_idx}"
        assert list(master_key_from_k2_k3(rks[2], rks[3])) == mk
        assert list(master_key_from_k4_k5(rks[4], rks[5])) == mk

    if verbose:
        print(
            "[PASS] Official NSA ARADI vectors verified (ePrint 2024/1240):\n"
            f"  - Fig. 3.4: 17/17 round keys, 16/16 intermediate states, CT = {format_words(ct_3_4)}\n"
            "  - Fig. 4.3: 10/10 Key-A counter blocks (n = 0..9) + decryption round-trip\n"
            "  - Fig. 4.4: 10/10 Key-B counter blocks (n = 0..9) + decryption round-trip\n"
            "  - Inverse key schedule: verified exact O(1) recovery for even_idx in {0,2,4,6,8,10,12,14}"
        )
    return True


if __name__ == "__main__":
    verify_official_vectors(verbose=True)
