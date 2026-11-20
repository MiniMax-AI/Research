// gf2n_pclmul.hpp -- GF(2^n) arithmetic for 1 <= n <= 512 using PCLMULQDQ.
//
// Element layout: W = ceil(n/64) little-endian 64-bit words (flat, stride W).
// Multiplication : carry-less schoolbook over words with _mm_clmulepi64_si128
//                  (Karatsuba for W == 2), then table-driven reduction.
//
// Reduction tables (built once per field):
//   BAS[t]      = x^(n+t) mod f(x),  t = 0 .. 128W-1-n            (W words each)
//   TH[p][c]    = XOR_{b: bit b of c} BAS[64W + 8p + b - n]      (high words)
//   TL[p][c]    = XOR_{b: bit b of c} BAS[8*(n/8) + 8p + b - n]  (top partial word)
// Both tables are indexed by a byte of the accumulator, so a full reduction
// costs (8W + NL) table lookups of W words each, no data-dependent branching.
#pragma once
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <immintrin.h>

#define GF2N_WMAX 8

struct Gf2n {
    int n = 0, W = 0;
    int NL = 0;                       // number of low (top-partial-word) table bytes
    uint64_t mod[GF2N_WMAX] = {0};    // low W words of f(x); bit n is implicit
    uint64_t mask_top = ~0ULL;        // mask for the top word
    std::vector<uint64_t> TH;         // (8W) * 256 * W
    std::vector<uint64_t> TL;         // NL  * 256 * W
    std::vector<uint64_t> BAS;        // (128W-n) * W   (also used for inversion)

    // ---- construction -----------------------------------------------------
    // fbits: low W words of the modulus polynomial f (degree n, monic).
    void init(int n_, const uint64_t *fbits) {
        n = n_; W = (n + 63) / 64;
        for (int i = 0; i < GF2N_WMAX; ++i) mod[i] = 0;
        for (int i = 0; i < W; ++i) mod[i] = fbits[i];
        int topbits = n - 64 * (W - 1);
        mask_top = (topbits == 64) ? ~0ULL : ((1ULL << topbits) - 1);

        const int maxT = 128 * W - 1 - n;          // t ranges 0..maxT
        BAS.assign((size_t)(maxT + 1) * W, 0);
        // BAS[0] = x^n mod f = mod
        for (int i = 0; i < W; ++i) BAS[i] = mod[i];
        const int topshift = (n - 1) % 64;          // bit n-1 lives in word W-1
        for (int t = 1; t <= maxT; ++t) {
            uint64_t *prev = &BAS[(size_t)(t - 1) * W];
            uint64_t *cur  = &BAS[(size_t)t * W];
            uint64_t top = (prev[W - 1] >> topshift) & 1ULL;   // bit n-1 of prev
            for (int i = W - 1; i > 0; --i) cur[i] = (prev[i] << 1) | (prev[i - 1] >> 63);
            cur[0] = prev[0] << 1;
            if (top) for (int i = 0; i < W; ++i) cur[i] ^= mod[i];   // x^n == mod
            cur[W - 1] &= mask_top;
        }
        // high table: byte positions 64W + 8p, p = 0..8W-1
        TH.assign((size_t)(8 * W) * 256 * W + 2, 0);
        for (int p = 0; p < 8 * W; ++p) {
            int t0 = 64 * W + 8 * p - n;
            for (int c = 1; c < 256; ++c) {
                uint64_t *dst = &TH[((size_t)p * 256 + c) * W];
                int b = __builtin_ctz(c);
                const uint64_t *src = &BAS[(size_t)(t0 + b) * W];
                for (int i = 0; i < W; ++i) dst[i] = src[i];
                int rest = c & (c - 1);
                if (rest) {
                    const uint64_t *s2 = &BAS[(size_t)(t0 + __builtin_ctz(rest)) * W];
                    for (int i = 0; i < W; ++i) dst[i] ^= s2[i];
                    // generic loop for >2 bits
                    int r2 = rest & (rest - 1);
                    while (r2) {
                        const uint64_t *s3 = &BAS[(size_t)(t0 + __builtin_ctz(r2)) * W];
                        for (int i = 0; i < W; ++i) dst[i] ^= s3[i];
                        r2 &= r2 - 1;
                    }
                }
            }
        }
        // low table: byte positions 8*(n/8) + 8p, p = 0..NL-1, covering up to 64W-1
        int b0 = 8 * (n / 8);
        NL = (64 * W - b0) / 8;
        TL.assign((size_t)NL * 256 * W + 2, 0);
        for (int p = 0; p < NL; ++p) {
            int t0 = b0 + 8 * p - n;               // may be negative for p == 0
            int skip = (t0 < 0) ? -t0 : 0;        // bits below n in the first byte
            for (int c = 1; c < 256; ++c) {
                uint64_t *dst = &TL[((size_t)p * 256 + c) * W];
                int cc = c & ~((1 << skip) - 1);
                if (!cc) continue;
                int r2 = cc;
                while (r2) {
                    int b = __builtin_ctz(r2);
                    const uint64_t *s3 = &BAS[(size_t)(t0 + b) * W];
                    for (int i = 0; i < W; ++i) dst[i] ^= s3[i];
                    r2 &= r2 - 1;
                }
            }
        }
    }

    // ---- primitives -------------------------------------------------------
    inline void zero(uint64_t *a) const { for (int i = 0; i < W; ++i) a[i] = 0; }
    inline void copy(uint64_t *d, const uint64_t *s) const { for (int i = 0; i < W; ++i) d[i] = s[i]; }
    inline void add(uint64_t *d, const uint64_t *s) const { for (int i = 0; i < W; ++i) d[i] ^= s[i]; }
    inline bool iszero(const uint64_t *a) const {
        uint64_t z = 0; for (int i = 0; i < W; ++i) z |= a[i]; return z == 0;
    }
    inline bool eq(const uint64_t *a, const uint64_t *b) const {
        for (int i = 0; i < W; ++i) if (a[i] != b[i]) return false; return true;
    }

    // ---- reduction --------------------------------------------------------
    // 4-way interleaved accumulators break the load->xor dependency chain that
    // otherwise serialises the (8W + NL) table lookups.
    template <int W_>
    inline void reduce_t(uint64_t *acc) const {
        constexpr int NV = (W_ + 1) / 2;
        __m128i t[4][NV];
        for (int u = 0; u < 4; ++u)
            for (int v = 0; v < NV; ++v) t[u][v] = _mm_setzero_si128();
        int k = 0;
        const uint64_t *th = TH.data();
        for (int j = 2 * W_ - 1; j >= W_; --j) {
            uint64_t c = acc[j];
            acc[j] = 0;
            const uint64_t *base = th + (size_t)(8 * (j - W_)) * 256 * W_;
            for (int b = 0; b < 8; ++b) {
                unsigned byte = (unsigned)((c >> (8 * b)) & 0xff);
                const uint64_t *s = base + ((size_t)b * 256 + byte) * W_;
                __m128i *d = t[k & 3];
                for (int v = 0; v < NV; ++v)
                    d[v] = _mm_xor_si128(d[v], _mm_loadu_si128((const __m128i *)(s + 2 * v)));
                ++k;
            }
        }
        const int b0 = 8 * (n / 8);
        const unsigned lm0 = (unsigned)(0xffULL << (n % 8)) & 0xffu;
        const uint64_t *tl = TL.data();
        for (int p = 0; p < NL; ++p) {
            int word = (b0 + 8 * p) / 64, sh = (b0 + 8 * p) % 64;
            unsigned byte = (unsigned)((acc[word] >> sh) & 0xff);
            if (p == 0) byte &= lm0;
            acc[word] ^= ((uint64_t)byte) << sh;
            const uint64_t *s = tl + ((size_t)p * 256 + byte) * W_;
            __m128i *d = t[k & 3];
            for (int v = 0; v < NV; ++v)
                d[v] = _mm_xor_si128(d[v], _mm_loadu_si128((const __m128i *)(s + 2 * v)));
            ++k;
        }
        uint64_t tmp[2 * NV + 2];
        for (int v = 0; v < NV; ++v) {
            __m128i r = _mm_xor_si128(_mm_xor_si128(t[0][v], t[1][v]), _mm_xor_si128(t[2][v], t[3][v]));
            _mm_storeu_si128((__m128i *)(tmp + 2 * v), r);
        }
        for (int i = 0; i < W_; ++i) acc[i] ^= tmp[i];
        acc[W_ - 1] &= mask_top;
    }

    inline void reduce(uint64_t *acc) const {          // acc has 2W words
        switch (W) {
            case 1: reduce_t<1>(acc); break;
            case 2: reduce_t<2>(acc); break;
            case 3: reduce_t<3>(acc); break;
            case 4: reduce_t<4>(acc); break;
            case 5: reduce_t<5>(acc); break;
            case 6: reduce_t<6>(acc); break;
            case 7: reduce_t<7>(acc); break;
            default: reduce_t<8>(acc); break;
        }
    }

    inline void mul(uint64_t *out, const uint64_t *a, const uint64_t *b) const {
        uint64_t acc[2 * GF2N_WMAX];
        for (int i = 0; i < 2 * W; ++i) acc[i] = 0;
        if (W == 2) {                                   // Karatsuba, 3 clmuls
            __m128i A = _mm_set_epi64x((long long)a[1], (long long)a[0]);
            __m128i B = _mm_set_epi64x((long long)b[1], (long long)b[0]);
            __m128i p0 = _mm_clmulepi64_si128(A, B, 0x00);
            __m128i p1 = _mm_clmulepi64_si128(A, B, 0x11);
            __m128i ax = _mm_xor_si128(_mm_srli_si128(A, 8), A);
            __m128i bx = _mm_xor_si128(_mm_srli_si128(B, 8), B);
            __m128i p2 = _mm_clmulepi64_si128(ax, bx, 0x00);
            __m128i mid = _mm_xor_si128(_mm_xor_si128(p2, p0), p1);
            acc[0] = (uint64_t)_mm_cvtsi128_si64(p0);
            acc[1] = (uint64_t)_mm_cvtsi128_si64(_mm_xor_si128(_mm_srli_si128(p0, 8), mid));
            acc[2] = (uint64_t)_mm_cvtsi128_si64(_mm_xor_si128(_mm_srli_si128(mid, 8), p1));
            acc[3] = (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p1, 8));
        } else {
            __m128i A[GF2N_WMAX], B[GF2N_WMAX], ACC[2 * GF2N_WMAX];
            for (int i = 0; i < W; ++i) { A[i] = _mm_set_epi64x(0, (long long)a[i]); B[i] = _mm_set_epi64x(0, (long long)b[i]); }
            for (int i = 0; i < 2 * W - 1; ++i) ACC[i] = _mm_setzero_si128();
            for (int i = 0; i < W; ++i)
                for (int j = 0; j < W; ++j)
                    ACC[i + j] = _mm_xor_si128(ACC[i + j], _mm_clmulepi64_si128(A[i], B[j], 0x00));
            for (int i = 0; i < 2 * W - 1; ++i) {
                acc[i]     ^= (uint64_t)_mm_cvtsi128_si64(ACC[i]);
                acc[i + 1] ^= (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(ACC[i], 8));
            }
        }
        reduce(acc);
        for (int i = 0; i < W; ++i) out[i] = acc[i];
    }

    // out += a*b   (fused multiply-add, the SpMV inner kernel)
    inline void muladd(uint64_t *out, const uint64_t *a, const uint64_t *b) const {
        uint64_t acc[2 * GF2N_WMAX];
        for (int i = 0; i < 2 * W; ++i) acc[i] = 0;
        if (W == 2) {
            __m128i A = _mm_set_epi64x((long long)a[1], (long long)a[0]);
            __m128i B = _mm_set_epi64x((long long)b[1], (long long)b[0]);
            __m128i p0 = _mm_clmulepi64_si128(A, B, 0x00);
            __m128i p1 = _mm_clmulepi64_si128(A, B, 0x11);
            __m128i ax = _mm_xor_si128(_mm_srli_si128(A, 8), A);
            __m128i bx = _mm_xor_si128(_mm_srli_si128(B, 8), B);
            __m128i p2 = _mm_clmulepi64_si128(ax, bx, 0x00);
            __m128i mid = _mm_xor_si128(_mm_xor_si128(p2, p0), p1);
            acc[0] = (uint64_t)_mm_cvtsi128_si64(p0);
            acc[1] = (uint64_t)_mm_cvtsi128_si64(_mm_xor_si128(_mm_srli_si128(p0, 8), mid));
            acc[2] = (uint64_t)_mm_cvtsi128_si64(_mm_xor_si128(_mm_srli_si128(mid, 8), p1));
            acc[3] = (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p1, 8));
        } else {
            __m128i A[GF2N_WMAX], B[GF2N_WMAX], ACC[2 * GF2N_WMAX];
            for (int i = 0; i < W; ++i) { A[i] = _mm_set_epi64x(0, (long long)a[i]); B[i] = _mm_set_epi64x(0, (long long)b[i]); }
            for (int i = 0; i < 2 * W - 1; ++i) ACC[i] = _mm_setzero_si128();
            for (int i = 0; i < W; ++i)
                for (int j = 0; j < W; ++j)
                    ACC[i + j] = _mm_xor_si128(ACC[i + j], _mm_clmulepi64_si128(A[i], B[j], 0x00));
            for (int i = 0; i < 2 * W - 1; ++i) {
                acc[i]     ^= (uint64_t)_mm_cvtsi128_si64(ACC[i]);
                acc[i + 1] ^= (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(ACC[i], 8));
            }
        }
        reduce(acc);
        for (int i = 0; i < W; ++i) out[i] ^= acc[i];
    }

    inline void sqr(uint64_t *out, const uint64_t *a) const { mul(out, a, a); }

    inline void mulx(uint64_t *a) const {                 // a <- x*a mod f
        uint64_t top = (a[W - 1] >> ((n - 1) % 64)) & 1ULL;
        for (int i = W - 1; i > 0; --i) a[i] = (a[i] << 1) | (a[i - 1] >> 63);
        a[0] <<= 1;
        if (top) for (int i = 0; i < W; ++i) a[i] ^= mod[i];
        a[W - 1] &= mask_top;
    }

    // inverse via Fermat + left-to-right addition chain:
    //   a^(2^n - 2) = ( a^(2^(n-1) - 1) )^2 ,  and  a^(2^k - 1) = (a^(2^(k-1)-1))^2 * a
    void inv(uint64_t *out, const uint64_t *a) const {
        uint64_t r[GF2N_WMAX], t[GF2N_WMAX];
        copy(r, a);                              // r = a^(2^1 - 1)
        for (int i = 1; i <= n - 2; ++i) {
            sqr(t, r);
            mul(r, t, a);                        // r = a^(2^(i+1) - 1)
        }
        sqr(out, r);                             // a^(2^n - 2)
    }
};
