#!/usr/bin/env python3
"""Derive a Phi_n-component solver from kemsolve.cpp.

The canonical LOCKER/ROLLO ring is R = F[Z]/(Z^n-1).  Its factor Z+1 carries a
trivial partial Frobenius (so every MinRank row coincides there and all the
(r+1)x(r+1) minors vanish identically); the attack lives on the cyclotomic
component R' = F[Z]/Phi_n.  For n prime with ord_n(2) = n-1 the polynomial Phi_n
is irreducible over GF(2) of degree n-1, so R' = F (x)_{F2} GF(2^{n-1}) is a
field and rho = id_F (x) sigma fixes F pointwise exactly as in the
irreducible-modulus LOCKER-I build.  The projection pi : R -> R' is the ring
homomorphism Z |-> Z, i.e. f_i |-> f_i + f_{n-1} on the n F-coefficients.

Reduction mod Phi_n is cheap: Z^n - 1 = (Z-1)Phi_n gives Z^n = 1, so a product
of degree <= 2n-4 needs ONE fold of period n plus one Z^{n-1} = 1+...+Z^{n-2}
step -- cheaper than the sparse modulus of the irreducible-modulus LOCKER-I build.
"""
import sys
from pathlib import Path

NIN = int(sys.argv[1])          # degree of the input ring Z^NIN - 1
N = NIN - 1                     # degree of Phi_NIN
out_path = sys.argv[2]
src_path = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).resolve().parent.parent / "schur_locker1/kemsolve.cpp"
src = src_path.read_text()

def mask(k):
    return f"(((u128)1 << {k}) - 1)"

old_arith = """static constexpr u128 MASK83 = (u128(1) << 83) - 1;
static inline Fe reduce(u128 lo, u128 hi) {  // value = lo + hi*2^128, reduce mod y^83+..
    const u128 low = lo & MASK83;
    const u128 high = (lo >> 83) | (hi << 45);
    u128 z = low ^ high ^ (high << 2) ^ (high << 4) ^ (high << 7);
    const u128 top = z >> 83;
    z = (z ^ top ^ (top << 2) ^ (top << 4) ^ (top << 7)) & MASK83;
    return {(uint64_t)z, (uint64_t)(z >> 64)};
}"""
new_arith = f"""static constexpr u128 MASKN = {mask(NIN)};
static constexpr u128 MASKP = {mask(N)};
// Phi_{NIN} = 1 + y + ... + y^{N}; y^{NIN} = 1 and y^{N} = 1 + ... + y^{N-1}.
static inline Fe reduce(u128 lo, u128 hi) {{  // value = lo + hi*2^128, degree <= {2*NIN-4}
    u128 z = (lo & MASKN) ^ ((lo >> {NIN}) | (hi << {128-NIN}));
    const u128 top = (z >> {N}) & 1;
    z &= MASKP;
    if (top) z ^= MASKP;
    return {{(uint64_t)z, (uint64_t)(z >> 64)}};
}}"""
assert old_arith in src
src = src.replace(old_arith, new_arith)

src = src.replace("""    Fe result{1, 0}, power = a;
    for (unsigned bit = 0; bit < 83; ++bit) {""",
f"""    Fe result{{1, 0}}, power = a;
    for (unsigned bit = 0; bit < {N}; ++bit) {{""")

src = src.replace("""static inline void mul_x(Fe &a) {  // multiply by y
    const uint64_t top = (a.hi >> 18) & 1;
    a = {a.lo << 1, ((a.hi << 1) | (a.lo >> 63)) & ((uint64_t(1) << 19) - 1)};
    if (top) a.lo ^= 0x95;  // y^83 = y^7+y^4+y^2+1
}""",
f"""static inline void mul_x(Fe &a) {{  // multiply by y
    const uint64_t top = (a.hi >> {N-1-64}) & 1;
    a = {{a.lo << 1, ((a.hi << 1) | (a.lo >> 63)) & ((uint64_t(1) << {N-64}) - 1)}};
    if (top) {{ a.lo ^= ~uint64_t(0); a.hi ^= (uint64_t(1) << {N-64}) - 1; }}
}}""")

src = src.replace("constexpr int N = 83;   // [L:F2] = ring index count",
f"constexpr int N = {N};   // [L:F2] = deg Phi_{NIN} (component carrying the attack)\n"
f"constexpr int NIN = {NIN}; // degree of the input ring Z^{NIN}-1 (projected away by pi)")

src = src.replace("""static Qe transpose_ring(const std::vector<unsigned char> &raw) {
    if (raw.size() < (size_t)N * 9) throw std::runtime_error("short ring element");
    Qe out = qzero();
    for (int j = 0; j < N; ++j) {
        uint64_t lo = 0;
        for (int b = 0; b < 8; ++b) lo |= (uint64_t)raw[9 * j + b] << (8 * b);
        const uint8_t hi = raw[9 * j + 8];""",
f"""// pi : F[Z]/(Z^{NIN}-1) -> F[Z]/Phi_{NIN},  f_i |-> f_i + f_{NIN-1}  (i < {N}).
static Qe transpose_ring(const std::vector<unsigned char> &raw) {{
    if (raw.size() < (size_t)NIN * 9) throw std::runtime_error("short ring element");
    uint64_t plo[NIN]; uint8_t phi[NIN];
    for (int j = 0; j < NIN; ++j) {{
        uint64_t lo = 0;
        for (int b = 0; b < 8; ++b) lo |= (uint64_t)raw[9 * j + b] << (8 * b);
        plo[j] = lo; phi[j] = raw[9 * j + 8];
    }}
    for (int j = 0; j < N; ++j) {{ plo[j] ^= plo[{N}]; phi[j] ^= phi[{N}]; }}
    Qe out = qzero();
    for (int j = 0; j < N; ++j) {{
        const uint64_t lo = plo[j];
        const uint8_t hi = phi[j];""")

src = src.replace("""// Itoh-Tsujii: a^(2^5893 - 2)
static Qe qinv(const Qe &a) {
    if (qiszero(a)) throw std::runtime_error("inverse of zero in Q");
    const int m = 5893 - 1;  // compute a^(2^m - 1) then square""",
f"""// Itoh-Tsujii: a^(2^(M*N) - 2),  M*N = 71*{N} = {71*N}
static Qe qinv(const Qe &a) {{
    if (qiszero(a)) throw std::runtime_error("inverse of zero in Q");
    const int m = M * N - 1;  // compute a^(2^m - 1) then square""")

for a, b in [("if (pk.size() != (size_t)N * 9)", "if (pk.size() != (size_t)NIN * 9)"),
             ("if (ct.size() != (size_t)N * 9 + 128)", "if (ct.size() != (size_t)NIN * 9 + 128)"),
             ("const std::vector<unsigned char> ctu(ct.begin(), ct.begin() + N * 9);",
              "const std::vector<unsigned char> ctu(ct.begin(), ct.begin() + NIN * 9);"),
             ("for (int i = 0; i < 64; ++i) m[i] = ct[N * 9 + i] ^ mask[i];",
              "for (int i = 0; i < 64; ++i) m[i] = ct[NIN * 9 + i] ^ mask[i];"),
             ("if (dig[i] != ct[N * 9 + 64 + i]) digest_ok = false;",
              "if (dig[i] != ct[NIN * 9 + 64 + i]) digest_ok = false;"),
             ("mc.insert(mc.end(), ct.begin() + N * 9, ct.begin() + N * 9 + 64);",
              "mc.insert(mc.end(), ct.begin() + NIN * 9, ct.begin() + NIN * 9 + 64);")]:
    assert a in src, a
    src = src.replace(a, b)

open(out_path, 'w').write(src)
print(f"wrote {out_path}: NIN={NIN} N={N} M*N={71*N}")
