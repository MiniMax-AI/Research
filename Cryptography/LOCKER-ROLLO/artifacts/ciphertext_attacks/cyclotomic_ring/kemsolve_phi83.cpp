// kemsolve_phi83.cpp -- recover the real KEM shared secret on the canonical
// Z^83 - 1 cyclotomic-ring variant of LOCKER-I (projected onto Phi_83).
//
// Construction (canonical-ring LOCKER-I variant compatible with reference/native/cli):
//   F  = GF(2^71) = F2[t]/(t^71+t^6+1)          (NTL BuildSparseIrred_GF2X(71))
//   R  = F[Z]/(Z^83 - 1),  Z^83 - 1 = (Z + 1) * Phi_83(Z)
// Projecting R onto R' = F[Z]/(Phi_83(Z)) via f_i |-> f_i + f_82 yields a field
// R' = F (x)_{F2} L with L = GF(2^82) = F2[y]/(Phi_83(y)), since Phi_83 is
// irreducible over GF(2) and gcd(82, 71) = 1.  Writing an element of R' as a
// 71-tuple over L ("transpose", tau) is an F2-algebra isomorphism, and the
// partial Frobenius rho = id_F (x) sigma_L fixes F pointwise.
//
// Ciphertext: u = e1 + e2*h with e1 = sum_k a_k eps_k, e2 = sum_k a_k eta_k,
// a_1..a_5 an F2-basis of the support E (dim 5 in F), eps_k,eta_k in L.  Applying
// rho^q (which fixes every a_k in F) gives, for every q,
//       rho^q(u) = sum_k a_k ( rho^q(eps_k) + rho^q(h) rho^q(eta_k) ).
// Hence the 9x6 matrix with row q = [ rho^q(u) | Z_q^(1) .. Z_q^(5) ], where
// Z_q^(k) = rho^q(eps_k) + H_q rho^q(eta_k) and H_q = rho^q(h), has rank <= 5.
// The MaxMinors linearisation solves for the Plucker coordinates of the 5x18
// matrix [v0|v1|rho v0|rho v1|...]; the Frobenius relations p(S+1,T+1)=p(S,T)^2
// isolate the true vector.  From it we read off pi(e2), hence E, then
// m = v xor SHA512(E) and ss = SHA512(m || u || v).
//
// Build: g++ -O3 -march=native -fopenmp -std=c++17 -o kemsolve_phi83 kemsolve_phi83.cpp

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#include <immintrin.h>
#endif
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using u128 = unsigned __int128;

// ---------------------------------------------------------------- GF(2^82) --
struct alignas(16) Fe {
    uint64_t lo, hi;
};
static inline bool nonzero(const Fe &a) { return (a.lo | a.hi) != 0; }
static inline Fe add(const Fe &a, const Fe &b) { return {a.lo ^ b.lo, a.hi ^ b.hi}; }

// Carry-less multiply: use PCLMULQDQ when the compiler targets it, otherwise a
// portable bit loop (selected automatically via __PCLMUL__ or NO_PCLMUL).
// Scalar fallback when compiled without -mavx2 / -mpclmul.
#if defined(NO_PCLMUL) || !defined(__PCLMUL__)
#define KEM_NO_PCLMUL 1
#endif

// 64x64 -> 128 carry-less product (pclmul, with a portable fallback)
static inline void clmul64(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
#ifdef KEM_NO_PCLMUL
    uint64_t r0 = 0, r1 = 0;
    for (int i = 0; i < 64; ++i)
        if ((b >> i) & 1) {
            r0 ^= a << i;
            if (i) r1 ^= a >> (64 - i);
        }
    lo = r0;
    hi = r1;
#else
    const __m128i p = _mm_clmulepi64_si128(_mm_set_epi64x(0, (long long)a),
                                           _mm_set_epi64x(0, (long long)b), 0x00);
    lo = (uint64_t)_mm_cvtsi128_si64(p);
    hi = (uint64_t)_mm_extract_epi64(p, 1);
#endif
}

// full 83x83 carry-less product, returned as (low 128 bits, bits 128..)
static inline void clmul(const Fe &a, const Fe &b, u128 &pl, u128 &ph) {
#ifdef KEM_NO_PCLMUL
    uint64_t l00, h00, l01, h01, l10, h10, l11, h11;
    clmul64(a.lo, b.lo, l00, h00);
    clmul64(a.lo, b.hi, l01, h01);
    clmul64(a.hi, b.lo, l10, h10);
    clmul64(a.hi, b.hi, l11, h11);
    const u128 mid = ((u128)l01 | ((u128)h01 << 64)) ^ ((u128)l10 | ((u128)h10 << 64));
    pl = ((u128)l00 | ((u128)h00 << 64)) ^ (mid << 64);
    ph = (mid >> 64) ^ ((u128)l11 | ((u128)h11 << 64));
#else
    const __m128i av = _mm_set_epi64x((long long)a.hi, (long long)a.lo);
    const __m128i bv = _mm_set_epi64x((long long)b.hi, (long long)b.lo);
    const __m128i p00 = _mm_clmulepi64_si128(av, bv, 0x00);
    const __m128i cr = _mm_xor_si128(_mm_clmulepi64_si128(av, bv, 0x01),
                                     _mm_clmulepi64_si128(av, bv, 0x10));
    const __m128i p11 = _mm_clmulepi64_si128(av, bv, 0x11);
    const u128 lo = (u128)(uint64_t)_mm_cvtsi128_si64(p00) |
                    ((u128)(uint64_t)_mm_extract_epi64(p00, 1) << 64);
    const u128 mid = (u128)(uint64_t)_mm_cvtsi128_si64(cr) |
                     ((u128)(uint64_t)_mm_extract_epi64(cr, 1) << 64);
    const u128 hi = (u128)(uint64_t)_mm_cvtsi128_si64(p11) |
                    ((u128)(uint64_t)_mm_extract_epi64(p11, 1) << 64);
    pl = lo ^ (mid << 64);
    ph = (mid >> 64) ^ hi;
#endif
}

static constexpr u128 MASKN = (((u128)1 << 83) - 1);
static constexpr u128 MASKP = (((u128)1 << 82) - 1);
// Phi_83 = 1 + y + ... + y^82; y^83 = 1 and y^82 = 1 + ... + y^81.
static inline Fe reduce(u128 lo, u128 hi) {  // value = lo + hi*2^128, degree <= 162
    u128 z = (lo & MASKN) ^ ((lo >> 83) | (hi << 45));
    const u128 top = (z >> 82) & 1;
    z &= MASKP;
    if (top) z ^= MASKP;
    return {(uint64_t)z, (uint64_t)(z >> 64)};
}
static inline Fe mul(const Fe &a, const Fe &b) {
    u128 pl, ph;
    clmul(a, b, pl, ph);
    return reduce(pl, ph);
}
static inline Fe square(const Fe &a) { return mul(a, a); }

static Fe inverse(const Fe &a) {
    if (!nonzero(a)) throw std::runtime_error("inverse of zero in GF(2^83)");
    Fe result{1, 0}, power = a;
    for (unsigned bit = 0; bit < 82; ++bit) {
        if (bit) result = mul(result, power);
        power = square(power);
    }
    return result;
}
static inline void mul_x(Fe &a) {  // multiply by y
    const uint64_t top = (a.hi >> 17) & 1;
    a = {a.lo << 1, ((a.hi << 1) | (a.lo >> 63)) & ((uint64_t(1) << 18) - 1)};
    if (top) { a.lo ^= ~uint64_t(0); a.hi ^= (uint64_t(1) << 18) - 1; }
}

// ------------------------------------------------- Q = L[T]/(T^71+T^6+1) -----
constexpr int M = 71;   // [F:F2]
constexpr int N = 82;   // [L:F2] = deg Phi_83 (component carrying the attack)
constexpr int NIN = 83; // degree of the input ring Z^83-1 (projected away by pi)
constexpr int R = 5;    // support dimension
constexpr int PAIRS = 9;

struct Qe {
    Fe c[M];
};

static Qe qzero() {
    Qe z;
    for (int i = 0; i < M; ++i) z.c[i] = {0, 0};
    return z;
}
static Qe qone() {
    Qe z = qzero();
    z.c[0] = {1, 0};
    return z;
}
static bool qiszero(const Qe &a) {
    for (int i = 0; i < M; ++i)
        if (nonzero(a.c[i])) return false;
    return true;
}
static Qe qadd(const Qe &a, const Qe &b) {
    Qe out;
    for (int i = 0; i < M; ++i) out.c[i] = add(a.c[i], b.c[i]);
    return out;
}
static Qe qmul(const Qe &a, const Qe &b) {
    Fe p[2 * M - 1];
    for (int i = 0; i < 2 * M - 1; ++i) p[i] = {0, 0};
    for (int i = 0; i < M; ++i) {
        if (!nonzero(a.c[i])) continue;
        for (int j = 0; j < M; ++j) {
            if (!nonzero(b.c[j])) continue;
            p[i + j] = add(p[i + j], mul(a.c[i], b.c[j]));
        }
    }
    for (int d = 2 * M - 2; d >= M; --d) {  // T^71 = T^6 + 1
        if (!nonzero(p[d])) continue;
        p[d - M] = add(p[d - M], p[d]);
        p[d - M + 6] = add(p[d - M + 6], p[d]);
    }
    Qe out;
    for (int i = 0; i < M; ++i) out.c[i] = p[i];
    return out;
}
// partial Frobenius: squares the L-coefficients, fixes T  (this is the Galois
// automorphism of R over F used by the attack -- NOT the field square)
static Qe rho(const Qe &a) {
    Qe out;
    for (int i = 0; i < M; ++i) out.c[i] = square(a.c[i]);
    return out;
}
// genuine field square in Q
static Qe qsquare(const Qe &a) {
    Fe p[2 * M - 1];
    for (int i = 0; i < 2 * M - 1; ++i) p[i] = {0, 0};
    for (int i = 0; i < M; ++i)
        if (nonzero(a.c[i])) p[2 * i] = square(a.c[i]);
    for (int d = 2 * M - 2; d >= M; --d) {
        if (!nonzero(p[d])) continue;
        p[d - M] = add(p[d - M], p[d]);
        p[d - M + 6] = add(p[d - M + 6], p[d]);
    }
    Qe out;
    for (int i = 0; i < M; ++i) out.c[i] = p[i];
    return out;
}

static Qe qpow2k(Qe a, int k) {
    for (int i = 0; i < k; ++i) a = qsquare(a);
    return a;
}
// Itoh-Tsujii: a^(2^(M*N) - 2),  M*N = 71*82 = 5822
static Qe qinv(const Qe &a) {
    if (qiszero(a)) throw std::runtime_error("inverse of zero in Q");
    const int m = M * N - 1;  // compute a^(2^m - 1) then square
    int top = 31;
    while (!((m >> top) & 1)) --top;
    Qe result = a;
    int k = 1;
    for (int i = top - 1; i >= 0; --i) {
        result = qmul(qpow2k(result, k), result);
        k *= 2;
        if ((m >> i) & 1) {
            result = qsquare(result);
            result = qmul(result, a);
            k += 1;
        }
    }
    return qsquare(result);
}

// ------------------------------------------------------------- utilities ----
static double elapsed(const std::chrono::steady_clock::time_point &t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}
static int threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
static int thread_id() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

static std::string read_file(const char *path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error(std::string("cannot open ") + path);
    return std::string(std::istreambuf_iterator<char>(in), {});
}
static std::string json_string(const std::string &text, const std::string &key) {
    size_t at = text.find('"' + key + '"');
    if (at == std::string::npos) throw std::runtime_error("missing JSON key " + key);
    at = text.find(':', at + key.size() + 2);
    at = text.find('"', at);
    const size_t end = text.find('"', at + 1);
    std::string v = text.substr(at + 1, end - at - 1);
    if (v.rfind("0x", 0) == 0 || v.rfind("0X", 0) == 0) v.erase(0, 2);
    return v;
}
static unsigned hex_digit(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    throw std::runtime_error("invalid hex digit");
}
static std::vector<unsigned char> unhex(std::string s) {
    if (s.rfind("0x", 0) == 0 || s.rfind("0X", 0) == 0) s.erase(0, 2);
    if (s.size() & 1) throw std::runtime_error("odd hex string");
    std::vector<unsigned char> out(s.size() / 2);
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (unsigned char)((hex_digit(s[2 * i]) << 4) | hex_digit(s[2 * i + 1]));
    return out;
}

// ---------------------------------------------------------------- SHA-512 ---
namespace sha512 {
static const uint64_t K[80] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL, 0xe9b5dba58189dbbcULL,
    0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL, 0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL,
    0xd807aa98a3030242ULL, 0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL, 0xc19bf174cf692694ULL,
    0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL, 0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL,
    0x2de92c6f592b0275ULL, 0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL, 0xbf597fc7beef0ee4ULL,
    0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL, 0x06ca6351e003826fULL, 0x142929670a0e6e70ULL,
    0x27b70a8546d22ffcULL, 0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
    0x650a73548baf63deULL, 0x766a0abb3c77b2a8ULL, 0x81c2c92e47edaee6ULL, 0x92722c851482353bULL,
    0xa2bfe8a14cf10364ULL, 0xa81a664bbc423001ULL, 0xc24b8b70d0f89791ULL, 0xc76c51a30654be30ULL,
    0xd192e819d6ef5218ULL, 0xd69906245565a910ULL, 0xf40e35855771202aULL, 0x106aa07032bbd1b8ULL,
    0x19a4c116b8d2d0c8ULL, 0x1e376c085141ab53ULL, 0x2748774cdf8eeb99ULL, 0x34b0bcb5e19b48a8ULL,
    0x391c0cb3c5c95a63ULL, 0x4ed8aa4ae3418acbULL, 0x5b9cca4f7763e373ULL, 0x682e6ff3d6b2b8a3ULL,
    0x748f82ee5defb2fcULL, 0x78a5636f43172f60ULL, 0x84c87814a1f0ab72ULL, 0x8cc702081a6439ecULL,
    0x90befffa23631e28ULL, 0xa4506cebde82bde9ULL, 0xbef9a3f7b2c67915ULL, 0xc67178f2e372532bULL,
    0xca273eceea26619cULL, 0xd186b8c721c0c207ULL, 0xeada7dd6cde0eb1eULL, 0xf57d4f7fee6ed178ULL,
    0x06f067aa72176fbaULL, 0x0a637dc5a2c898a6ULL, 0x113f9804bef90daeULL, 0x1b710b35131c471bULL,
    0x28db77f523047d84ULL, 0x32caab7b40c72493ULL, 0x3c9ebe0a15c9bebcULL, 0x431d67c49c100d4cULL,
    0x4cc5d4becb3e42b6ULL, 0x597f299cfc657e2aULL, 0x5fcb6fab3ad6faecULL, 0x6c44198c4a475817ULL};
static inline uint64_t ror(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
static void digest(const unsigned char *msg, size_t len, unsigned char out[64]) {
    uint64_t h[8] = {0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL,
                     0xa54ff53a5f1d36f1ULL, 0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
                     0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL};
    const size_t total = ((len + 16) / 128 + 1) * 128;
    std::vector<unsigned char> buf(total, 0);
    if (len) memcpy(buf.data(), msg, len);
    buf[len] = 0x80;
    const uint64_t bits = (uint64_t)len * 8;
    for (int i = 0; i < 8; ++i) buf[total - 1 - i] = (unsigned char)(bits >> (8 * i));
    for (size_t off = 0; off < total; off += 128) {
        uint64_t w[80];
        for (int i = 0; i < 16; ++i) {
            uint64_t v = 0;
            for (int j = 0; j < 8; ++j) v = (v << 8) | buf[off + 8 * i + j];
            w[i] = v;
        }
        for (int i = 16; i < 80; ++i) {
            const uint64_t s0 = ror(w[i - 15], 1) ^ ror(w[i - 15], 8) ^ (w[i - 15] >> 7);
            const uint64_t s1 = ror(w[i - 2], 19) ^ ror(w[i - 2], 61) ^ (w[i - 2] >> 6);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint64_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 80; ++i) {
            const uint64_t S1 = ror(e, 14) ^ ror(e, 18) ^ ror(e, 41);
            const uint64_t ch = (e & f) ^ (~e & g);
            const uint64_t t1 = hh + S1 + ch + K[i] + w[i];
            const uint64_t S0 = ror(a, 28) ^ ror(a, 34) ^ ror(a, 39);
            const uint64_t mj = (a & b) ^ (a & c) ^ (b & c);
            const uint64_t t2 = S0 + mj;
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) out[8 * i + j] = (unsigned char)(h[i] >> (56 - 8 * j));
}
}  // namespace sha512

// ------------------------------------------------------- transposed views ---
// pi : F[Z]/(Z^83-1) -> F[Z]/Phi_83,  f_i |-> f_i + f_82  (i < 82).
static Qe transpose_ring(const std::vector<unsigned char> &raw) {
    if (raw.size() < (size_t)NIN * 9) throw std::runtime_error("short ring element");
    uint64_t plo[NIN]; uint8_t phi[NIN];
    for (int j = 0; j < NIN; ++j) {
        uint64_t lo = 0;
        for (int b = 0; b < 8; ++b) lo |= (uint64_t)raw[9 * j + b] << (8 * b);
        plo[j] = lo; phi[j] = raw[9 * j + 8];
    }
    for (int j = 0; j < N; ++j) { plo[j] ^= plo[82]; phi[j] ^= phi[82]; }
    Qe out = qzero();
    for (int j = 0; j < N; ++j) {
        const uint64_t lo = plo[j];
        const uint8_t hi = phi[j];
        for (int bit = 0; bit < M; ++bit) {
            const bool set = bit < 64 ? ((lo >> bit) & 1) : ((hi >> (bit - 64)) & 1);
            if (!set) continue;
            if (j < 64) out.c[bit].lo |= uint64_t(1) << j;
            else out.c[bit].hi |= uint64_t(1) << (j - 64);
        }
    }
    return out;
}
// coefficients live in F = GF(2^71): they need 71 bits, hence u128 storage
static std::vector<u128> transpose_back(const Qe &v) {
    std::vector<u128> out(N, 0);
    for (int bit = 0; bit < M; ++bit) {
        const uint64_t lo = v.c[bit].lo, hi = v.c[bit].hi;
        if (!(lo | hi)) continue;
        const u128 mbit = (u128)1 << bit;
        for (int j = 0; j < 64; ++j)
            if ((lo >> j) & 1) out[j] |= mbit;
        for (int j = 0; j < N - 64; ++j)
            if ((hi >> j) & 1) out[64 + j] |= mbit;
    }
    return out;
}

// ---------------------------------------------------------- combinations ----
using Combo5 = std::array<int, 5>;
using Combo6 = std::array<int, 6>;
static void combos5(int next, int depth, Combo5 &cur, std::vector<Combo5> &out) {
    if (depth == 5) { out.push_back(cur); return; }
    for (int x = next; x <= PAIRS - (5 - depth); ++x) {
        cur[depth] = x;
        combos5(x + 1, depth + 1, cur, out);
    }
}
static void combos6(int next, int depth, Combo6 &cur, std::vector<Combo6> &out) {
    if (depth == 6) { out.push_back(cur); return; }
    for (int x = next; x <= PAIRS - (6 - depth); ++x) {
        cur[depth] = x;
        combos6(x + 1, depth + 1, cur, out);
    }
}
static uint32_t code5(const Combo5 &c) {
    uint32_t z = 0;
    for (int i = 0; i < 5; ++i) z |= (uint32_t)c[i] << (4 * i);
    return z;
}

// ------------------------------------------------ stage 1: MaxMinors matrix -
static void build_matrix(const std::vector<Combo6> &six, const std::vector<Combo5> &five,
                         const std::map<uint32_t, uint32_t> &five_index,
                         const std::array<Qe, PAIRS> &H, const std::array<Qe, PAIRS> &U,
                         std::vector<Fe> &matrix, uint32_t &nrows, uint32_t &ncols) {
    std::vector<Qe> products(1 << PAIRS);
    products[0] = qone();
    for (unsigned mask = 1; mask < products.size(); ++mask) {
        const unsigned bit = __builtin_ctz(mask);
        products[mask] = qmul(products[mask ^ (1u << bit)], H[bit]);
    }
    ncols = (uint32_t)five.size() * 32;
    nrows = (uint32_t)six.size() * M;
    matrix.assign((size_t)nrows * ncols, Fe{0, 0});
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t eq = 0; eq < (int64_t)six.size(); ++eq) {
        const Combo6 &sel = six[eq];
        for (int omit = 0; omit < 6; ++omit) {
            Combo5 pairs{};
            for (int src = 0, dst = 0; src < 6; ++src)
                if (src != omit) pairs[dst++] = sel[src];
            const uint32_t block = five_index.at(code5(pairs));
            for (unsigned side = 0; side < 32; ++side) {
                unsigned hmask = 0;
                for (int k = 0; k < 5; ++k)
                    if ((side >> k) & 1) hmask |= 1u << pairs[k];
                const Qe coeff = qmul(U[sel[omit]], products[hmask]);
                const uint32_t col = block * 32 + side;
                for (int deg = 0; deg < M; ++deg)
                    matrix[((size_t)eq * M + deg) * ncols + col] = coeff.c[deg];
            }
        }
    }
}

// ------------------------------------------- stage 1: reduced MaxMinors ----
// Every 6-subset used contains pair 0, so the 5-subset obtained by dropping
// pair 0 (the "leaf" block) occurs in exactly one equation, while a 5-subset
// containing 0 occurs in PAIRS-5 of them.  Eliminating the 2^R leaf columns
// inside their own equation block is therefore an exact Schur complement: it
// removes 2^R columns and (M - 2^R) of the M rows of that block and leaves a
// system with exactly the same kernel on the surviving Plucker coordinates.
//   full   : C(PAIRS-1,5)*M      x C(PAIRS,5)*2^R     = 3976 x 4032
//   reduced: C(PAIRS-1,5)*(M-2^R) x C(PAIRS-1,4)*2^R   = 2184 x 2240
// Both the flop count and (decisively, since the elimination is bandwidth
// bound) the streamed footprint drop by ~3.5x.
struct Reduced {
    std::vector<Fe> mat;
    uint32_t nrows = 0, ncols = 0, nleaf = 0, nfull = 0;
    std::vector<uint32_t> gcol;    // full Plucker column -> reduced column
    std::vector<uint32_t> roff;    // per equation: first consistency row
    std::vector<uint32_t> eoff;    // per equation: offset into loc2c/pivcol
    std::vector<uint32_t> poff;    // per equation: offset into pivrow (in Fe)
    std::vector<uint32_t> loc2c;   // 6*32 local columns -> reduced column
    std::vector<uint32_t> pivcol;  // determined leaf columns (full indexing)
    std::vector<Fe> pivrow;        // pivot rows, 6*32 entries each
};

// rank (and the pivot pattern) of the 71 x 32 leaf block of one equation
static void leaf_pivots(const Combo5 &leaf, const Qe &u0, const std::vector<Qe> &products,
                        bool piv[32], int &rank) {
    std::vector<Fe> blk((size_t)M * 32, Fe{0, 0});
    for (int side = 0; side < 32; ++side) {
        unsigned hmask = 0;
        for (int k = 0; k < 5; ++k)
            if ((side >> k) & 1) hmask |= 1u << leaf[k];
        const Qe c = qmul(u0, products[hmask]);   // exactly pass 1's leaf column
        for (int deg = 0; deg < M; ++deg) blk[(size_t)deg * 32 + side] = c.c[deg];
    }
    rank = 0;
    for (int c = 0; c < 32; ++c) piv[c] = false;
    for (int c = 0; c < 32; ++c) {
        int r = rank;
        while (r < M && !nonzero(blk[(size_t)r * 32 + c])) ++r;
        if (r == M) continue;
        if (r != rank)
            for (int j = 0; j < 32; ++j) std::swap(blk[(size_t)r * 32 + j], blk[(size_t)rank * 32 + j]);
        const Fe inv = inverse(blk[(size_t)rank * 32 + c]);
        for (int j = 0; j < 32; ++j) blk[(size_t)rank * 32 + j] = mul(blk[(size_t)rank * 32 + j], inv);
        for (int r2 = 0; r2 < M; ++r2) {
            if (r2 == rank) continue;
            const Fe f = blk[(size_t)r2 * 32 + c];
            if (!nonzero(f)) continue;
            for (int j = 0; j < 32; ++j)
                blk[(size_t)r2 * 32 + j] = add(blk[(size_t)r2 * 32 + j], mul(f, blk[(size_t)rank * 32 + j]));
        }
        piv[c] = true;
        ++rank;
    }
}

static void build_reduced(const std::vector<Combo6> &six, const std::vector<Combo5> &five,
                          const std::map<uint32_t, uint32_t> &five_index,
                          const std::array<Qe, PAIRS> &H, const std::array<Qe, PAIRS> &U,
                          Reduced &R) {
    std::vector<Qe> products(1 << PAIRS);
    products[0] = qone();
    for (unsigned mask = 1; mask < products.size(); ++mask) {
        const unsigned bit = __builtin_ctz(mask);
        products[mask] = qmul(products[mask ^ (1u << bit)], H[bit]);
    }
    const uint32_t neq = (uint32_t)six.size();
    const uint32_t nblocks = (uint32_t)five.size();
    R.nleaf = 0;
    for (const Combo5 &x : five) if (x[0] != 0) ++R.nleaf;   // leaves come first
    R.nfull = nblocks * 32;
    R.gcol.assign(R.nfull, UINT32_MAX);
    for (uint32_t b = R.nleaf; b < nblocks; ++b)
        for (int s = 0; s < 32; ++s) R.gcol[b * 32 + s] = (b - R.nleaf) * 32 + s;

    // ---- pass 0: leaf ranks, pivot pattern, column/row budget
    std::vector<std::array<bool, 32>> piv(neq);
    std::vector<int> lrank(neq);
    std::vector<uint32_t> leafblk(neq);
    R.roff.assign(neq + 1, 0);
    R.eoff.assign(neq + 1, 0);
    R.poff.assign(neq + 1, 0);
    uint32_t next = (nblocks - R.nleaf) * 32;
    for (uint32_t eq = 0; eq < neq; ++eq) {
        Combo5 leaf{};
        for (int src = 1, dst = 0; src < 6; ++src) leaf[dst++] = six[eq][src];
        leafblk[eq] = five_index.at(code5(leaf));
        leaf_pivots(leaf, U[six[eq][0]], products, piv[eq].data(), lrank[eq]);
        for (int s = 0; s < 32; ++s)
            if (!piv[eq][s]) R.gcol[leafblk[eq] * 32 + s] = next++;
        R.roff[eq + 1] = R.roff[eq] + (uint32_t)(M - lrank[eq]);
        R.eoff[eq + 1] = R.eoff[eq] + 192;
        R.poff[eq + 1] = R.poff[eq] + (uint32_t)lrank[eq];
    }
    R.nrows = R.roff[neq];
    R.ncols = next;
    R.mat.assign((size_t)R.nrows * R.ncols, Fe{0, 0});
    R.loc2c.assign((size_t)neq * 192, 0);
    R.pivcol.assign(R.poff[neq], 0);
    R.pivrow.assign((size_t)R.poff[neq] * 192, Fe{0, 0});

    // ---- pass 1: build each 71 x 192 block, eliminate its leaf columns
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t eq = 0; eq < (int64_t)neq; ++eq) {
        const Combo6 &sel = six[eq];
        std::vector<Fe> blk((size_t)M * 192, Fe{0, 0});
        uint32_t *l2c = &R.loc2c[(size_t)eq * 192];
        for (int omit = 0; omit < 6; ++omit) {
            Combo5 pairs{};
            for (int src = 0, dst = 0; src < 6; ++src)
                if (src != omit) pairs[dst++] = sel[src];
            const uint32_t block = five_index.at(code5(pairs));
            for (int side = 0; side < 32; ++side) {
                unsigned hmask = 0;
                for (int k = 0; k < 5; ++k)
                    if ((side >> k) & 1) hmask |= 1u << pairs[k];
                const Qe coeff = qmul(U[sel[omit]], products[hmask]);
                const int lc = omit * 32 + side;
                l2c[lc] = R.gcol[block * 32 + side];
                for (int deg = 0; deg < M; ++deg) blk[(size_t)deg * 192 + lc] = coeff.c[deg];
            }
        }
        int row = 0;
        for (int c = 0; c < 32; ++c) {
            if (!piv[eq][c]) continue;
            int r = row;
            while (r < M && !nonzero(blk[(size_t)r * 192 + c])) ++r;
            if (r == M) continue;  // cannot happen: pass 0 saw a pivot here
            if (r != row)
                for (int j = 0; j < 192; ++j)
                    std::swap(blk[(size_t)r * 192 + j], blk[(size_t)row * 192 + j]);
            Fe *p = &blk[(size_t)row * 192];
            const Fe inv = inverse(p[c]);
            // a leaf column that was not pivoted may sit before c, so the whole
            // row has to be scaled / eliminated, not just the tail
            p[c] = Fe{1, 0};
            for (int j = 0; j < 192; ++j)
                if (j != c) p[j] = mul(p[j], inv);
            for (int r2 = 0; r2 < M; ++r2) {
                if (r2 == row) continue;
                Fe *q = &blk[(size_t)r2 * 192];
                const Fe f = q[c];
                if (!nonzero(f)) continue;
                q[c] = Fe{0, 0};
                for (int j = 0; j < 192; ++j)
                    if (j != c) q[j] = add(q[j], mul(f, p[j]));
            }
            R.pivcol[R.poff[eq] + (uint32_t)row] = leafblk[eq] * 32 + (uint32_t)c;
            ++row;
        }
        // store the pivot rows only now: after the whole loop each of them has a
        // single 1 among the leaf columns, so back substitution is a direct read
        for (int r = 0; r < row; ++r)
            std::copy(&blk[(size_t)r * 192], &blk[(size_t)r * 192 + 192],
                      &R.pivrow[(size_t)(R.poff[eq] + (uint32_t)r) * 192]);
        for (int r = row; r < M; ++r) {
            if (R.roff[eq] + (uint32_t)(r - row) >= R.roff[eq + 1]) break;  // defensive
            Fe *dst = &R.mat[(size_t)(R.roff[eq] + (uint32_t)(r - row)) * R.ncols];
            const Fe *q = &blk[(size_t)r * 192];
            for (int j = 0; j < 192; ++j)
                if (l2c[j] != UINT32_MAX) dst[l2c[j]] = q[j];
        }
    }
}

// expand a reduced kernel basis back to the full Plucker coordinate vector
static void expand_basis(const Reduced &R, uint32_t ncols_red,
                         const std::vector<Fe> &red, size_t dim, uint32_t neq,
                         std::vector<Fe> &full, uint32_t nfull) {
    full.assign((size_t)dim * nfull, Fe{0, 0});
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t f = 0; f < (int64_t)dim; ++f) {
        Fe *out = &full[(size_t)f * nfull];
        const Fe *x = &red[(size_t)f * ncols_red];
        for (uint32_t j = 0; j < nfull; ++j)
            if (R.gcol[j] != UINT32_MAX) out[j] = x[R.gcol[j]];
        std::vector<Fe> v(192);
        for (uint32_t eq = 0; eq < neq; ++eq) {
            const uint32_t *l2c = &R.loc2c[(size_t)eq * 192];
            for (int j = 0; j < 192; ++j) v[j] = (l2c[j] == UINT32_MAX) ? Fe{0, 0} : x[l2c[j]];
            for (uint32_t pi = R.poff[eq]; pi < R.poff[eq + 1]; ++pi) {
                const Fe *pr = &R.pivrow[(size_t)pi * 192];
                Fe acc = Fe{0, 0};
                for (int j = 0; j < 192; ++j)
                    if (nonzero(v[j]) && nonzero(pr[j])) acc = add(acc, mul(pr[j], v[j]));
                out[R.pivcol[pi]] = acc;
            }
        }
    }
}

// ------------------------------- stage 2: kernel over GF(2^83) (blocked) -----
struct KernelResult {
    std::vector<Fe> basis;
    std::vector<uint32_t> free_cols;
};

static KernelResult lkernel(std::vector<Fe> &mat, uint32_t nrows, uint32_t ncols, int B) {
    std::vector<uint32_t> pivots;
    pivots.reserve(std::min(nrows, ncols));
    const int nth = threads();
    std::vector<Fe> coef((size_t)nth * B);
    uint32_t rank = 0, col = 0;
    std::vector<uint32_t> bcols, brows;
    while (rank < nrows && col < ncols) {
        const uint32_t base = rank;
        bcols.clear();
        brows.clear();
        uint32_t c = col;
        // ---- panel: pick up to B pivots, Gauss-Jordan inside the block columns
        while ((int)bcols.size() < B && rank < nrows && c < ncols) {
            bool accepted = false;
            uint32_t cand = rank;
            while (cand < nrows) {
                uint32_t found = cand;
                while (found < nrows && !nonzero(mat[(size_t)found * ncols + c])) ++found;
                if (found == nrows) break;
                if (found != rank)
                    for (uint32_t j = 0; j < ncols; ++j)
                        std::swap(mat[(size_t)found * ncols + j], mat[(size_t)rank * ncols + j]);
                Fe *p = mat.data() + (size_t)rank * ncols;
                for (size_t j = 0; j < bcols.size(); ++j) {  // reduce vs earlier pivots
                    const Fe *pj = mat.data() + (size_t)brows[j] * ncols;
                    const Fe f = p[bcols[j]];
                    if (!nonzero(f)) continue;
                    p[bcols[j]] = Fe{0, 0};
                    for (uint32_t jj = bcols[j] + 1; jj < ncols; ++jj)
                        p[jj] = add(p[jj], mul(f, pj[jj]));
                }
                if (nonzero(p[c])) { accepted = true; break; }
                cand = rank + 1;  // reduction killed column c: keep row, try the next
            }
            if (!accepted) { ++c; continue; }
            Fe *p = mat.data() + (size_t)rank * ncols;
            const Fe scale = inverse(p[c]);
            p[c] = Fe{1, 0};
            for (uint32_t j = c + 1; j < ncols; ++j) p[j] = mul(p[j], scale);
            for (uint32_t k = base; k < rank; ++k) {  // clear c from earlier pivots
                Fe *pk = mat.data() + (size_t)k * ncols;
                const Fe f = pk[c];
                if (!nonzero(f)) continue;
                pk[c] = Fe{0, 0};
                for (uint32_t j = c + 1; j < ncols; ++j) pk[j] = add(pk[j], mul(f, p[j]));
            }
            bcols.push_back(c);
            brows.push_back(rank);
            pivots.push_back(c);
            ++rank;
            ++c;
        }
        const int nb = (int)bcols.size();
        if (!nb) break;
        // ---- apply: one streaming pass per remaining row for the whole block
        const uint32_t lo = bcols.front();
        const uint32_t rr0 = rank;
#pragma omp parallel for schedule(static)
        for (int64_t rr = rr0; rr < (int64_t)nrows; ++rr) {
            Fe *row = mat.data() + (size_t)rr * ncols;
            Fe *cf = coef.data() + (size_t)thread_id() * B;
            bool any = false;
            for (int i = 0; i < nb; ++i) {
                cf[i] = row[bcols[i]];
                if (nonzero(cf[i])) any = true;
            }
            if (!any) continue;
            uint32_t j = lo;
#if defined(__AVX512F__) && defined(__VPCLMULQDQ__) && defined(__AVX512BW__)
            // four coefficients at a time: accumulate unreduced 166-bit products
            const __m512i zero = _mm512_setzero_si512();
            u128 lo4[4], hi4[4];
            for (; j + 4 <= ncols; j += 4) {
                __m512i accA = _mm512_loadu_si512((const void *)(row + j));
                __m512i accB = zero;
                for (int i = 0; i < nb; ++i) {
                    if (!nonzero(cf[i])) continue;
                    const __m512i a = _mm512_broadcast_i32x4(
                        _mm_set_epi64x((long long)cf[i].hi, (long long)cf[i].lo));
                    const __m512i b =
                        _mm512_loadu_si512((const void *)(mat.data() + (size_t)brows[i] * ncols + j));
                    const __m512i p00 = _mm512_clmulepi64_epi128(a, b, 0x00);
                    const __m512i cr = _mm512_xor_si512(_mm512_clmulepi64_epi128(a, b, 0x01),
                                                        _mm512_clmulepi64_epi128(a, b, 0x10));
                    const __m512i p11 = _mm512_clmulepi64_epi128(a, b, 0x11);
                    accA = _mm512_xor_si512(accA,
                        _mm512_xor_si512(p00, _mm512_bslli_epi128(cr, 8)));
                    accB = _mm512_xor_si512(accB,
                        _mm512_xor_si512(_mm512_bsrli_epi128(cr, 8), p11));
                }
                _mm512_storeu_si512((void *)lo4, accA);
                _mm512_storeu_si512((void *)hi4, accB);
                for (int t = 0; t < 4; ++t) row[j + t] = reduce(lo4[t], hi4[t]);
            }
#endif
            for (; j < ncols; ++j) {
                u128 al = (u128)row[j].lo | ((u128)row[j].hi << 64), ah = 0;
                for (int i = 0; i < nb; ++i) {
                    if (!nonzero(cf[i])) continue;
                    u128 pl, ph;
                    clmul(cf[i], mat[(size_t)brows[i] * ncols + j], pl, ph);
                    al ^= pl;
                    ah ^= ph;
                }
                row[j] = reduce(al, ah);
            }
        }
        if (getenv("KEMINV")) {
            long bad1 = 0, bad2 = 0;
            for (uint32_t r = rank; r < nrows; ++r)
                for (size_t j = 0; j < pivots.size(); ++j)
                    if (nonzero(mat[(size_t)r * ncols + pivots[j]])) ++bad1;
            for (size_t i = 0; i < pivots.size(); ++i)
                for (size_t j = 0; j <= i; ++j) {
                    const Fe v = mat[i * ncols + pivots[j]];
                    const bool want = (i == j);
                    if (nonzero(v) != want) ++bad2;
                }
            fprintf(stderr, "block end: pivots=%zu pool_nonzero=%ld pivot_not_id=%ld\n",
                    pivots.size(), bad1, bad2);
            if (bad1 || bad2) exit(3);
        }
        col = c;
    }
    // ---- free columns and kernel basis by back substitution
    std::vector<uint32_t> is_pivot(ncols, 0);
    for (uint32_t p : pivots) is_pivot[p] = 1;
    KernelResult res;
    for (uint32_t j = 0; j < ncols; ++j)
        if (!is_pivot[j]) res.free_cols.push_back(j);
    const size_t dim = res.free_cols.size();
    res.basis.assign(dim * ncols, Fe{0, 0});
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t fi = 0; fi < (int64_t)dim; ++fi) {
        const uint32_t f = res.free_cols[fi];
        res.basis[(size_t)fi * ncols + f] = Fe{1, 0};
        for (int64_t i = (int64_t)pivots.size() - 1; i >= 0; --i) {
            const Fe *p = mat.data() + (size_t)i * ncols;
            // a pivot row may be nonzero at a column that was skipped inside its
            // own block, so the free column has to be accounted for separately
            Fe acc = (f < pivots[i]) ? p[f] : Fe{0, 0};
            for (uint32_t j = pivots[i] + 1; j < ncols; ++j) {
                const Fe v = res.basis[(size_t)fi * ncols + j];
                if (nonzero(v) && nonzero(p[j])) acc = add(acc, mul(p[j], v));
            }
            res.basis[(size_t)fi * ncols + pivots[i]] = acc;
        }
    }
    return res;
}

// -------------------- stage 3: Frobenius relations (binary linear algebra) --
static inline void xor_field_at(uint64_t *row, size_t offset, const Fe &a) {
    const size_t blk = offset >> 6;
    const unsigned sh = offset & 63;
    if (sh == 0) {
        row[blk] ^= a.lo;
        row[blk + 1] ^= a.hi;
    } else {
        row[blk] ^= a.lo << sh;
        row[blk + 1] ^= (a.lo >> (64 - sh)) ^ (a.hi << sh);
        row[blk + 2] ^= a.hi >> (64 - sh);
    }
}

// kernel over GF(2) of the matrix whose rows are `rows` (nrows x ncols bits)
static std::vector<std::vector<uint64_t>> binary_kernel(std::vector<uint64_t> &rows,
                                                       size_t nrows, size_t ncols,
                                                       size_t words) {
    std::vector<uint32_t> pivots;
    size_t rank = 0;
    for (size_t col = 0; col < ncols && rank < nrows; ++col) {
        size_t found = rank;
        while (found < nrows &&
               !((rows[found * words + (col >> 6)] >> (col & 63)) & 1))
            ++found;
        if (found == nrows) continue;
        if (found != rank)
            for (size_t w = 0; w < words; ++w)
                std::swap(rows[found * words + w], rows[rank * words + w]);
        const uint64_t *p = rows.data() + rank * words;
        const size_t w0 = col >> 6;
#pragma omp parallel for schedule(static)
        for (int64_t ii = (int64_t)rank + 1; ii < (int64_t)nrows; ++ii) {
            uint64_t *r = rows.data() + (size_t)ii * words;
            if (!((r[w0] >> (col & 63)) & 1)) continue;
            for (size_t w = w0; w < words; ++w) r[w] ^= p[w];
        }
        pivots.push_back((uint32_t)col);
        ++rank;
    }
    std::vector<uint32_t> is_pivot(ncols, 0);
    for (uint32_t p : pivots) is_pivot[p] = 1;
    std::vector<std::vector<uint64_t>> basis;
    for (size_t f = 0; f < ncols; ++f) {
        if (is_pivot[f]) continue;
        std::vector<uint64_t> x(words, 0);
        x[f >> 6] |= uint64_t(1) << (f & 63);
        for (int64_t i = (int64_t)pivots.size() - 1; i >= 0; --i) {
            const uint64_t *p = rows.data() + (size_t)i * words;
            uint64_t acc = (p[f >> 6] >> (f & 63)) & 1;
            for (size_t gj = (size_t)i + 1; gj < pivots.size(); ++gj) {
                const uint32_t pc = pivots[gj];
                if ((x[pc >> 6] >> (pc & 63)) & 1)
                    acc ^= (p[pc >> 6] >> (pc & 63)) & 1;
            }
            if (acc) x[pivots[i] >> 6] ^= uint64_t(1) << (pivots[i] & 63);
        }
        basis.push_back(std::move(x));
    }
    return basis;
}

// --------------------------------------------- stage 4: rebuild generator ---
struct Plucker {
    const std::vector<Combo5> *five;
    const std::map<uint32_t, uint32_t> *index;
    std::vector<Fe> values;
    size_t ncols;
    Fe at(const Combo5 &s, unsigned mask) const {
        return values[(size_t)index->at(code5(s)) * 32 + mask];
    }
};

static Fe det5(Fe m[5][5]) {
    // Gaussian elimination over GF(2^83)
    Fe det{1, 0};
    for (int c = 0; c < 5; ++c) {
        int piv = -1;
        for (int r = c; r < 5; ++r)
            if (nonzero(m[r][c])) { piv = r; break; }
        if (piv < 0) return Fe{0, 0};
        if (piv != c) {
            for (int j = 0; j < 5; ++j) std::swap(m[c][j], m[piv][j]);
        }
        const Fe inv = inverse(m[c][c]);
        det = mul(det, m[c][c]);
        for (int j = c; j < 5; ++j) m[c][j] = mul(m[c][j], inv);
        for (int r = 0; r < 5; ++r) {
            if (r == c) continue;
            const Fe f = m[r][c];
            if (!nonzero(f)) continue;
            for (int j = c; j < 5; ++j) m[r][j] = add(m[r][j], mul(f, m[c][j]));
        }
    }
    return det;
}

// index of a column selection {(pair,side)} in the Plucker coordinate vector
static size_t coord_index(const Plucker &P, const std::vector<std::pair<int, int>> &sel) {
    Combo5 s{};
    unsigned mask = 0;
    for (int i = 0; i < 5; ++i) {
        s[i] = sel[i].first;
        if (sel[i].second) mask |= 1u << i;
    }
    return (size_t)P.index->at(code5(s)) * 32 + mask;
}

int main(int argc, char **argv) {
    const auto t0 = std::chrono::steady_clock::now();
    std::string pk_hex, ct_hex;
    int max_rel = getenv("MAXREL") ? atoi(getenv("MAXREL")) : 0;  // 0 = all
    if (argc >= 3) {
        pk_hex = argv[1];
        ct_hex = argv[2];
        if (argc >= 4) max_rel = atoi(argv[3]);
    } else if (argc == 2) {
        const std::string text = read_file(argv[1]);
        pk_hex = json_string(text, "public_key");
        ct_hex = json_string(text, "encapsulation");
    } else {
        std::string line;
        while (std::getline(std::cin, line)) {
            const size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            std::string key = line.substr(0, colon), val = line.substr(colon + 1);
            while (!key.empty() && isspace((unsigned char)key.back())) key.pop_back();
            while (!val.empty() && isspace((unsigned char)val.front())) val.erase(0, 1);
            if (key == "public_key" || key == "public-key") pk_hex = val;
            if (key == "encapsulation") ct_hex = val;
        }
    }
    if (pk_hex.empty() || ct_hex.empty()) {
        std::cerr << "usage: kemsolve PUBLIC_KEY_HEX CIPHERTEXT_HEX [MAX_RELATIONS]\n";
        return 2;
    }
    const std::vector<unsigned char> pk = unhex(pk_hex);
    const std::vector<unsigned char> ct = unhex(ct_hex);
    if (pk.size() != (size_t)NIN * 9) throw std::runtime_error("unexpected public key size");
    if (ct.size() != (size_t)NIN * 9 + 128) throw std::runtime_error("unexpected ciphertext size");
    const std::vector<unsigned char> ctu(ct.begin(), ct.begin() + NIN * 9);

    std::array<Qe, PAIRS> H, U;
    H[0] = transpose_ring(pk);
    U[0] = transpose_ring(ctu);
    for (int q = 1; q < PAIRS; ++q) {
        H[q] = rho(H[q - 1]);
        U[q] = rho(U[q - 1]);
    }

    std::vector<Combo5> five_all, five;
    std::vector<Combo6> six_all, six;
    Combo5 c5{};
    Combo6 c6{};
    combos5(0, 0, c5, five_all);
    combos6(0, 0, c6, six_all);
    for (const Combo5 &x : five_all) if (x[0] != 0) five.push_back(x);  // leaves first
    for (const Combo5 &x : five_all) if (x[0] == 0) five.push_back(x);
    for (const Combo6 &x : six_all)
        if (x[0] == 0) six.push_back(x);
    std::map<uint32_t, uint32_t> five_index;
    for (uint32_t i = 0; i < five.size(); ++i) five_index[code5(five[i])] = i;

    Reduced red;
    build_reduced(six, five, five_index, H, U, red);
    const uint32_t ncols = red.nfull;
    std::cerr << "reduced matrix " << red.nrows << "x" << red.ncols
              << " (full " << six.size() * M << "x" << ncols << ") t=" << elapsed(t0) << "\n";

    KernelResult ker0 = lkernel(red.mat, red.nrows, red.ncols,
                                getenv("BLK") ? atoi(getenv("BLK")) : 16);
    red.mat.clear();
    red.mat.shrink_to_fit();
    KernelResult ker;
    ker.free_cols = ker0.free_cols;
    expand_basis(red, red.ncols, ker0.basis, ker0.free_cols.size(), (uint32_t)six.size(),
                 ker.basis, ncols);
    ker0.basis.clear();
    ker0.basis.shrink_to_fit();
    if (const char *kf = getenv("KERNFILE")) {  // debug: load an external kernel basis
        std::ifstream in(kf, std::ios::binary);
        uint32_t nc, dm;
        in.read((char *)&nc, 4);
        in.read((char *)&dm, 4);
        ker.free_cols.assign(dm, 0);
        ker.basis.assign((size_t)dm * nc, Fe{0, 0});
        in.read((char *)ker.basis.data(), (std::streamsize)dm * nc * sizeof(Fe));
        std::cerr << "loaded external kernel " << dm << "x" << nc << "\n";
    }
    const size_t dim = ker.free_cols.size();
    std::cerr << "L-kernel dim=" << dim << " t=" << elapsed(t0) << "\n";
    if (const char *kf = getenv("KEMDUMPKERN")) {
        FILE *fp = fopen(kf, "wb");
        fwrite(&ncols, 4, 1, fp);
        uint32_t dm = (uint32_t)dim;
        fwrite(&dm, 4, 1, fp);
        fwrite(ker.basis.data(), sizeof(Fe), ker.basis.size(), fp);
        fclose(fp);
    }
    if (dim == 0 || dim > 512) throw std::runtime_error("unexpected kernel dimension");

    // ---- Frobenius relations: v[second] = v[first]^2 for every shiftable pair
    std::vector<std::pair<uint32_t, uint32_t>> relations;
    for (uint32_t ci = 0; ci < five.size(); ++ci) {
        if (five[ci].back() == PAIRS - 1) continue;
        Combo5 sh = five[ci];
        for (int &x : sh) ++x;
        const uint32_t si = five_index.at(code5(sh));
        for (uint32_t mask = 0; mask < 32; ++mask)
            relations.emplace_back(ci * 32 + mask, si * 32 + mask);
    }
    size_t nrel = relations.size();
    if (max_rel > 0 && (size_t)max_rel < nrel) nrel = max_rel;
    else if (max_rel == 0 && getenv("MAXREL") == nullptr) nrel = 256;
    const size_t unknowns = dim * N;
    const size_t words = (unknowns + 63) / 64 + 2;
    const size_t brows = nrel * N;
    std::vector<uint64_t> bin(brows * words, 0);
    // A[(r,b)][(f,j)] = bit_b( y^j * left_r + y^(2j) * right_r )
    // with left_r = basis_f[col(S+1,T+1)], right_r = square(basis_f[col(S,T)]).
    // Parallel over relations: each thread owns a disjoint block of rows.
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t r = 0; r < (int64_t)nrel; ++r) {
        std::vector<Fe> v(N);
        for (size_t f = 0; f < dim; ++f) {
            Fe l = ker.basis[f * ncols + relations[r].second];
            Fe rt = square(ker.basis[f * ncols + relations[r].first]);
            for (int j = 0; j < N; ++j) {
                v[j] = add(l, rt);
                mul_x(l);
                mul_x(rt);
                mul_x(rt);
            }
            for (int b = 0; b < N; ++b) {
                uint64_t lo = 0, hi = 0;
                for (int j = 0; j < N; ++j) {
                    const bool set = b < 64 ? ((v[j].lo >> b) & 1) : ((v[j].hi >> (b - 64)) & 1);
                    if (!set) continue;
                    if (j < 64) lo |= uint64_t(1) << j;
                    else hi |= uint64_t(1) << (j - 64);
                }
                xor_field_at(bin.data() + ((size_t)r * N + b) * words, (size_t)f * N, Fe{lo, hi});
            }
        }
    }
    std::cerr << "binary system " << brows << "x" << unknowns << " t=" << elapsed(t0) << "\n";
    if (const char *df = getenv("KEMDUMP")) {
        FILE *fp = fopen(df, "wb");
        const size_t dw = (unknowns + 63) / 64;
        std::vector<uint64_t> buf(dw);
        for (size_t r = 0; r < brows; ++r) {
            for (size_t w = 0; w < dw; ++w) buf[w] = bin[r * words + w];
            fwrite(buf.data(), 8, dw, fp);
        }
        fclose(fp);
    }
    std::vector<std::vector<uint64_t>> kbin = binary_kernel(bin, brows, unknowns, words);
    std::cerr << "binary kernel dim=" << kbin.size() << " t=" << elapsed(t0) << "\n";
    if (kbin.size() != 1 && nrel < relations.size()) {
        // not enough relations to isolate the solution: use all of them
        nrel = relations.size();
        const size_t brows2 = nrel * N;
        bin.assign(brows2 * words, 0);
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t r = 0; r < (int64_t)nrel; ++r) {
            std::vector<Fe> v(N);
            for (size_t f = 0; f < dim; ++f) {
                Fe l = ker.basis[f * ncols + relations[r].second];
                Fe rt = square(ker.basis[f * ncols + relations[r].first]);
                for (int j = 0; j < N; ++j) {
                    v[j] = add(l, rt);
                    mul_x(l);
                    mul_x(rt);
                    mul_x(rt);
                }
                for (int b = 0; b < N; ++b) {
                    uint64_t lo = 0, hi = 0;
                    for (int j = 0; j < N; ++j) {
                        const bool set = b < 64 ? ((v[j].lo >> b) & 1) : ((v[j].hi >> (b - 64)) & 1);
                        if (!set) continue;
                        if (j < 64) lo |= uint64_t(1) << j;
                        else hi |= uint64_t(1) << (j - 64);
                    }
                    xor_field_at(bin.data() + ((size_t)r * N + b) * words, (size_t)f * N, Fe{lo, hi});
                }
            }
        }
        kbin = binary_kernel(bin, brows2, unknowns, words);
        std::cerr << "binary kernel (all relations) dim=" << kbin.size()
                  << " t=" << elapsed(t0) << "\n";
    }
    bin.clear();
    bin.shrink_to_fit();
    if (kbin.empty()) throw std::runtime_error("Frobenius system has trivial kernel");

    // ---- combine: lambda_f from the kernel vector, then the Plucker vector
    Plucker P;
    P.five = &five;
    P.index = &five_index;
    P.ncols = ncols;
    P.values.assign(ncols, Fe{0, 0});
    for (size_t k = 0; k < kbin.size(); ++k) {
        P.values.assign(ncols, Fe{0, 0});
        std::vector<Fe> lambda(dim, Fe{0, 0});
        for (size_t f = 0; f < dim; ++f) {
            uint64_t lo = 0, hi = 0;
            for (int b = 0; b < N; ++b)
                if ((kbin[k][(f * N + b) >> 6] >> ((f * N + b) & 63)) & 1) {
                    if (b < 64) lo |= uint64_t(1) << b;
                    else hi |= uint64_t(1) << (b - 64);
                }
            lambda[f] = Fe{lo, hi};
        }
#pragma omp parallel for schedule(static)
        for (int64_t j = 0; j < (int64_t)ncols; ++j) {
            Fe v{0, 0};
            for (size_t f = 0; f < dim; ++f) {
                const Fe b = ker.basis[f * ncols + j];
                if (nonzero(lambda[f]) && nonzero(b)) v = add(v, mul(lambda[f], b));
            }
            P.values[j] = v;
        }

        // ---- reconstruct the 5x18 generator (chart: pairs 0..4, side 1)
        std::vector<std::pair<int, int>> chart;
        for (int i = 0; i < 5; ++i) chart.push_back({i, 1});
        const Fe chart_minor = P.values[coord_index(P, chart)];
        if (!nonzero(chart_minor)) {
            std::cerr << "chart minor is zero\n";
            continue;
        }
        Fe Gm[5][2 * PAIRS];
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 2 * PAIRS; ++j) Gm[i][j] = Fe{0, 0};
        for (int i = 0; i < 5; ++i) Gm[i][2 * i + 1] = Fe{1, 0};
        const Fe inv_chart = inverse(chart_minor);
        for (int pair = 5; pair < PAIRS; ++pair)
            for (int side = 0; side < 2; ++side)
                for (int row = 0; row < 5; ++row) {
                    std::vector<std::pair<int, int>> sel;
                    for (int i = 0; i < 5; ++i)
                        if (i != row) sel.push_back({i, 1});
                    sel.push_back({pair, side});
                    std::sort(sel.begin(), sel.end());
                    Gm[row][2 * pair + side] = mul(P.values[coord_index(P, sel)], inv_chart);
                }
        for (int i = 0; i < 5; ++i) {
            std::vector<std::pair<int, int>> sel;
            for (int j = 0; j < 5; ++j) sel.push_back({j, j == i ? 0 : 1});
            std::sort(sel.begin(), sel.end());
            Gm[i][2 * i] = mul(P.values[coord_index(P, sel)], inv_chart);
        }
        bool ok = true;
        for (int i = 0; i < 5 && ok; ++i) {
            int pair = -1, side = -1;
            for (int p = 5; p < PAIRS && pair < 0; ++p)
                for (int s = 0; s < 2; ++s)
                    if (nonzero(Gm[i][2 * p + s])) { pair = p; side = s; break; }
            if (pair < 0) { ok = false; break; }
            for (int row = 0; row < 5; ++row) {
                if (row == i) continue;
                std::vector<std::pair<int, int>> sel;
                for (int j = 0; j < 5; ++j)
                    if (j != i && j != row) sel.push_back({j, 1});
                sel.push_back({i, 0});
                sel.push_back({pair, side});
                std::sort(sel.begin(), sel.end());
                const Fe d = mul(P.values[coord_index(P, sel)], inv_chart);
                const Fe num = add(d, mul(Gm[i][2 * i], Gm[row][2 * pair + side]));
                Gm[row][2 * i] = mul(num, inverse(Gm[i][2 * pair + side]));
            }
        }
        if (!ok) { std::cerr << "no nonsingular alternate chart\n"; continue; }
        // consistency: every maximal minor must match
        bool good = true;
        for (size_t ci = 0; ci < five.size() && good; ++ci)
            for (unsigned mask = 0; mask < 32 && good; ++mask) {
                Fe m[5][5];
                for (int i = 0; i < 5; ++i) {
                    const int col = 2 * five[ci][i] + ((mask >> i) & 1);
                    for (int r = 0; r < 5; ++r) m[r][i] = Gm[r][col];
                }
                const Fe d = mul(det5(m), chart_minor);
                if (add(d, P.values[ci * 32 + mask]).lo || add(d, P.values[ci * 32 + mask]).hi)
                    good = false;
            }
        std::cerr << "generator check=" << good << " t=" << elapsed(t0) << "\n";
        if (!good) continue;

        // ---- solve for the a_k:  U_q = sum_k c_k ( G[k,2q] + H_q G[k,2q+1] )
        Qe A[PAIRS][5];
        for (int q = 0; q < PAIRS; ++q)
            for (int k = 0; k < 5; ++k) {
                Qe g0 = qzero(), g1 = qzero();
                g0.c[0] = Gm[k][2 * q];
                g1.c[0] = Gm[k][2 * q + 1];
                A[q][k] = qadd(g0, qmul(H[q], g1));
            }
        // Gaussian elimination on the 9x5 system
        Qe rhs[PAIRS];
        for (int q = 0; q < PAIRS; ++q) rhs[q] = U[q];
        int prow[5];
        int rank = 0;
        for (int c = 0; c < 5 && rank < PAIRS; ++c) {
            int piv = -1;
            for (int r = rank; r < PAIRS; ++r)
                if (!qiszero(A[r][c])) { piv = r; break; }
            if (piv < 0) continue;
            for (int j = 0; j < 5; ++j) std::swap(A[rank][j], A[piv][j]);
            std::swap(rhs[rank], rhs[piv]);
            const Qe inv = qinv(A[rank][c]);
            if (getenv("KEMINV2")) {
                const Qe chk = qmul(A[rank][c], inv);
                std::cerr << "qinv check: c[0]=" << std::hex << chk.c[0].lo << " c[1]=" << chk.c[1].lo
                          << " others_nonzero=" << std::dec;
                int cnt = 0;
                for (int i = 1; i < M; ++i) if (nonzero(chk.c[i])) ++cnt;
                std::cerr << cnt << "\n";
            }
            for (int j = 0; j < 5; ++j) A[rank][j] = qmul(A[rank][j], inv);
            rhs[rank] = qmul(rhs[rank], inv);
            for (int r = 0; r < PAIRS; ++r) {
                if (r == rank) continue;
                if (qiszero(A[r][c])) continue;
                const Qe f = A[r][c];
                for (int j = 0; j < 5; ++j) A[r][j] = qadd(A[r][j], qmul(f, A[rank][j]));
                rhs[r] = qadd(rhs[r], qmul(f, rhs[rank]));
            }
            prow[c] = rank;
            ++rank;
        }
        if (rank != 5) { std::cerr << "singular 9x5 system\n"; continue; }
        Qe coef[5];
        for (int c = 0; c < 5; ++c) coef[c] = rhs[prow[c]];
        // residual check
        bool resid = true;
        for (int q = 0; q < PAIRS && resid; ++q) {
            Qe acc = qzero();
            for (int k = 0; k < 5; ++k) {
                Qe g0 = qzero(), g1 = qzero();
                g0.c[0] = Gm[k][2 * q];
                g1.c[0] = Gm[k][2 * q + 1];
                acc = qadd(acc, qmul(coef[k], qadd(g0, qmul(H[q], g1))));
            }
            if (!qiszero(qadd(acc, U[q]))) resid = false;
        }
        std::cerr << "9x5 solve residual ok=" << resid << " t=" << elapsed(t0) << "\n";
        if (!resid) continue;

        // ---- e2 = sum_k c_k * v1_k ; then the support
        Qe B = qzero();
        for (int k = 0; k < 5; ++k) {
            Qe g1 = qzero();
            g1.c[0] = Gm[k][1];
            B = qadd(B, qmul(coef[k], g1));
        }
        // Frobenius self-check
        bool frob = true;
        Qe cur = B;
        for (int q = 0; q < PAIRS && frob; ++q) {
            Qe want = qzero();
            for (int k = 0; k < 5; ++k) {
                Qe g1 = qzero();
                g1.c[0] = Gm[k][2 * q + 1];
                want = qadd(want, qmul(coef[k], g1));
            }
            if (!qiszero(qadd(cur, want))) frob = false;
            cur = rho(cur);
        }
        std::cerr << "frobenius check=" << frob << "\n";
        const std::vector<u128> e2 = transpose_back(B);
        if (getenv("KEMB")) {
            fprintf(stderr, "B:");
            for (int i = 0; i < M; ++i)
                fprintf(stderr, " %llx:%llx", (unsigned long long)B.c[i].hi, (unsigned long long)B.c[i].lo);
            fprintf(stderr, "\n");
        }

        // reduced row echelon over GF(2) of the N coordinates
        std::vector<u128> rows(e2);
        std::vector<u128> supp;
        {
            std::vector<u128> r(rows);
            int rk = 0;
            for (int c = 0; c < M && rk < (int)r.size(); ++c) {
                int piv = -1;
                for (int i = rk; i < (int)r.size(); ++i)
                    if ((r[i] >> c) & 1) { piv = i; break; }
                if (piv < 0) continue;
                std::swap(r[rk], r[piv]);
                for (int i = 0; i < (int)r.size(); ++i)
                    if (i != rk && ((r[i] >> c) & 1)) r[i] ^= r[rk];
                ++rk;
            }
            supp.assign(r.begin(), r.begin() + rk);
        }
        std::cerr << "support dim=" << supp.size() << " t=" << elapsed(t0) << "\n";
        if (supp.size() != R) continue;

        // ---- m = v xor SHA512(support), ss = SHA512(m || u || v)
        unsigned char supp_bytes[R * 9];
        for (int i = 0; i < R; ++i)
            for (int b = 0; b < 9; ++b) supp_bytes[9 * i + b] = (unsigned char)(supp[i] >> (8 * b));
        unsigned char mask[64], dig[64];
        sha512::digest(supp_bytes, sizeof(supp_bytes), mask);
        unsigned char m[64];
        for (int i = 0; i < 64; ++i) m[i] = ct[NIN * 9 + i] ^ mask[i];
        sha512::digest(m, 64, dig);
        bool digest_ok = true;
        for (int i = 0; i < 64; ++i)
            if (dig[i] != ct[NIN * 9 + 64 + i]) digest_ok = false;
        std::cerr << "digest ok=" << digest_ok << "\n";
        if (!digest_ok) continue;
        std::vector<unsigned char> mc;
        mc.insert(mc.end(), m, m + 64);
        mc.insert(mc.end(), ctu.begin(), ctu.end());
        mc.insert(mc.end(), ct.begin() + NIN * 9, ct.begin() + NIN * 9 + 64);
        unsigned char ss[64];
        sha512::digest(mc.data(), mc.size(), ss);
        for (int i = 0; i < 64; ++i) printf("%02x", ss[i]);
        printf("\n");
        std::cerr << "total t=" << elapsed(t0) << "\n";
        return 0;
    }
    std::cerr << "no candidate verified\n";
    return 1;
}
