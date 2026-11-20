// attack_param.cpp -- parameterised version of the Frobenius/MaxMinors
// support-recovery attack, used to measure how the attack cost scales with
//   m  = [F:F2]        (support field degree; F = GF(2^m))
//   n  = [L:F2]        (ring degree;      L = GF(2^n), R = F (x) L)
//   r  = dim_{F2} E    (support dimension)
//   p  = number of partial-Frobenius indices used by the MaxMinors model
//
// Instance model (all arithmetic done in the transposed world Q = L[T]/(f)):
//   u = sum_k a_k ( v0_k + h v1_k ),  a_k an F2-basis of E (subset of F),
//   v0_k, v1_k random elements of L, h random in R.
// The attack must recover span{a_k} from (h, u) alone.
//
// Build:
//   g++ -O3 -march=native -fopenmp -std=c++17 -DPM=71 -DPN=83 -DPR=5 -DPP=9 \
//       -DFMOD=0x41ULL -DRMOD=0x95ULL -o attack_param attack_param.cpp

#include <immintrin.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using u128 = unsigned __int128;

#ifndef PM
#define PM 71
#endif
#ifndef PN
#define PN 83
#endif
#ifndef PR
#define PR 5
#endif
#ifndef PP
#define PP 9
#endif
#ifndef FMOD
#define FMOD 0x41ULL /* t^71 + t^6 + 1 */
#endif
#ifndef RMOD
#define RMOD 0x95ULL /* Z^83 + Z^7 + Z^4 + Z^2 + 1 */
#endif

static constexpr int M = PM, N = PN, R = PR, PAIRS = PP;
static constexpr u128 MASKN = ((u128)1 << N) - 1;
static constexpr uint64_t RMODL = RMOD;  // low part of the L-modulus
static constexpr uint64_t FMODL = FMOD;  // low part of the T-modulus

// ---------------------------------------------------------------- L = GF(2^n)
struct alignas(16) Fe {
    uint64_t lo, hi;
};
static inline bool nonzero(const Fe &a) { return (a.lo | a.hi) != 0; }
static inline Fe add(const Fe &a, const Fe &b) { return {a.lo ^ b.lo, a.hi ^ b.hi}; }

static inline void clmul64(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
    const __m128i p = _mm_clmulepi64_si128(_mm_set_epi64x(0, (long long)a),
                                           _mm_set_epi64x(0, (long long)b), 0x00);
    lo = (uint64_t)_mm_cvtsi128_si64(p);
    hi = (uint64_t)_mm_extract_epi64(p, 1);
}
// full n x n carry-less product as (low 128 bits, high 128 bits)
static inline void clmul(const Fe &a, const Fe &b, u128 &pl, u128 &ph) {
    uint64_t l00, h00, l01, h01, l10, h10, l11, h11;
    clmul64(a.lo, b.lo, l00, h00);
    clmul64(a.lo, b.hi, l01, h01);
    clmul64(a.hi, b.lo, l10, h10);
    clmul64(a.hi, b.hi, l11, h11);
    const u128 mid = ((u128)l01 | ((u128)h01 << 64)) ^ ((u128)l10 | ((u128)h10 << 64));
    pl = ((u128)l00 | ((u128)h00 << 64)) ^ (mid << 64);
    ph = (mid >> 64) ^ ((u128)l11 | ((u128)h11 << 64));
}
// fold the bits above n back down using the sparse modulus
static inline u128 fold(u128 z) {
    // z = low + top*2^n with top the overflow; z == low + top*sum_k x^k (mod P).
    // The top*2^n part must be dropped before folding, otherwise it is folded twice.
    for (int pass = 0; pass < 4; ++pass) {
        const u128 top = z >> N;
        if (!top) break;
        z &= MASKN;
        for (uint64_t mm = RMODL; mm; mm &= mm - 1) z ^= top << __builtin_ctzll(mm);
    }
    return z & MASKN;
}
static inline Fe reduce(u128 lo, u128 hi) {
    const u128 low = lo & MASKN;
    const u128 high = (lo >> N) | (N < 128 ? (hi << (128 - N)) : 0);
    u128 z = low;
    for (uint64_t mm = RMODL; mm; mm &= mm - 1) z ^= high << __builtin_ctzll(mm);
    z = fold(z);
    return Fe{(uint64_t)z, (uint64_t)(z >> 64)};
}
static inline Fe mul(const Fe &a, const Fe &b) {
    u128 pl, ph;
    clmul(a, b, pl, ph);
    return reduce(pl, ph);
}
static inline Fe square(const Fe &a) { return mul(a, a); }
static Fe inverse(const Fe &a) {
    if (!nonzero(a)) throw std::runtime_error("inverse of zero");
    Fe r{1, 0}, p = a;
    for (int b = 0; b < N; ++b) {
        if (b) r = mul(r, p);
        p = mul(p, p);
    }
    return r;
}
static inline void mul_x(Fe &a) {  // multiply by the generator y of L
    Fe t = a;
    a.lo = t.lo << 1;
    a.hi = ((t.hi << 1) | (t.lo >> 63));
    if (N < 128) a.hi &= (N - 64 >= 64) ? ~0ULL : (((uint64_t)1 << (N - 64)) - 1);
    if (N >= 65 && ((t.hi >> (N - 65)) & 1)) a.lo ^= RMODL;  // y^n = sum of low terms
    if (N < 65 && ((t.lo >> (N - 1)) & 1)) a.lo ^= RMODL;
}

// ------------------------------------------------- Q = L[T]/(T^m + FMOD) -----
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
    Qe o;
    for (int i = 0; i < M; ++i) o.c[i] = add(a.c[i], b.c[i]);
    return o;
}
static Qe qmul(const Qe &a, const Qe &b) {
    std::vector<Fe> p(2 * M - 1, Fe{0, 0});
    for (int i = 0; i < M; ++i) {
        if (!nonzero(a.c[i])) continue;
        for (int j = 0; j < M; ++j) {
            if (!nonzero(b.c[j])) continue;
            p[i + j] = add(p[i + j], mul(a.c[i], b.c[j]));
        }
    }
    for (int d = 2 * M - 2; d >= M; --d) {
        if (!nonzero(p[d])) continue;
        for (uint64_t mm = FMODL; mm; mm &= mm - 1)
            p[d - M + __builtin_ctzll(mm)] = add(p[d - M + __builtin_ctzll(mm)], p[d]);
    }
    Qe o;
    for (int i = 0; i < M; ++i) o.c[i] = p[i];
    return o;
}
static Qe rho(const Qe &a) {  // partial Frobenius id_F (x) sigma_L
    Qe o;
    for (int i = 0; i < M; ++i) o.c[i] = square(a.c[i]);
    return o;
}
static Qe qsquare(const Qe &a) {
    std::vector<Fe> p(2 * M - 1, Fe{0, 0});
    for (int i = 0; i < M; ++i)
        if (nonzero(a.c[i])) p[2 * i] = square(a.c[i]);
    for (int d = 2 * M - 2; d >= M; --d) {
        if (!nonzero(p[d])) continue;
        for (uint64_t mm = FMODL; mm; mm &= mm - 1)
            p[d - M + __builtin_ctzll(mm)] = add(p[d - M + __builtin_ctzll(mm)], p[d]);
    }
    Qe o;
    for (int i = 0; i < M; ++i) o.c[i] = p[i];
    return o;
}
static Qe qpow2k(Qe a, int k) {
    for (int i = 0; i < k; ++i) a = qsquare(a);
    return a;
}
static Qe qinv(const Qe &a) {
    const int n = M * N;
    const int m = n - 1;
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

static double elapsed(const std::chrono::steady_clock::time_point &t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}
static int thread_id() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}
static int threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// ------------------------------------------------------------------ combos --
static void combos_rec(int next, int need, int total, std::vector<int> &cur,
                       std::vector<std::vector<int>> &out) {
    if (!need) { out.push_back(cur); return; }
    for (int x = next; x <= total - need; ++x) {
        cur.push_back(x);
        combos_rec(x + 1, need - 1, total, cur, out);
        cur.pop_back();
    }
}
static uint32_t code_of(const std::vector<int> &c) {
    uint32_t z = 0;
    for (size_t i = 0; i < c.size(); ++i) z |= (uint32_t)c[i] << (4 * i);
    return z;
}

int main(int argc, char **argv) {
    const auto t0 = std::chrono::steady_clock::now();
    unsigned long seed = argc > 1 ? strtoul(argv[1], nullptr, 10) : 12345;
    int BLK = argc > 2 ? atoi(argv[2]) : 16;
    size_t maxrel = argc > 3 ? (size_t)atoll(argv[3]) : 256;
    std::mt19937_64 rng(seed);
    auto randN = [&]() {
        Fe a{rng(), rng()};
        if (N < 64) a.lo &= ((uint64_t)1 << N) - 1, a.hi = 0;
        else if (N < 128) a.hi &= ((uint64_t)1 << (N - 64)) - 1;
        return a;
    };
    auto randF = [&]() { return (u128)(rng() & (((u128)1 << M) - 1)); };

    // ---- instance: h random, u = sum_k a_k (v0_k + h v1_k)
    Qe h = qzero();
    for (int i = 0; i < M; ++i) h.c[i] = randN();
    std::vector<u128> a(R);
    std::vector<Fe> v0(R), v1(R);
    for (int k = 0; k < R; ++k) {
        a[k] = randF();
        v0[k] = randN();
        v1[k] = randN();
    }
    // a_k are F2-independent (retry until the rank is r)
    auto rank_of = [&](const std::vector<u128> &v) {
        std::vector<u128> r(v);
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
        return rk;
    };
    while (rank_of(a) != R)
        for (int k = 0; k < R; ++k) a[k] = randF();

    Qe u = qzero();
    for (int k = 0; k < R; ++k) {
        Qe ak = qzero();
        for (int i = 0; i < M; ++i)
            if ((a[k] >> i) & 1) ak.c[i] = Fe{1, 0};
        Qe g0 = qzero(), g1 = qzero();
        g0.c[0] = v0[k];
        g1.c[0] = v1[k];
        u = qadd(u, qmul(ak, qadd(g0, qmul(h, g1))));
    }

    std::vector<Qe> H(PAIRS), U(PAIRS);
    H[0] = h;
    U[0] = u;
    for (int q = 1; q < PAIRS; ++q) {
        H[q] = rho(H[q - 1]);
        U[q] = rho(U[q - 1]);
    }

    std::vector<std::vector<int>> five, six_all, six;
    std::vector<int> cur;
    combos_rec(0, R, PAIRS, cur, five);
    combos_rec(0, R + 1, PAIRS, cur, six_all);
    for (auto &c : six_all)
        if (c[0] == 0) six.push_back(c);
    std::map<uint32_t, uint32_t> five_index;
    for (uint32_t i = 0; i < five.size(); ++i) five_index[code_of(five[i])] = i;
    const int NMASK = 1 << R;

    const uint32_t ncols = (uint32_t)five.size() * NMASK;
    const uint32_t nrows = (uint32_t)six.size() * M;

    // ---- MaxMinors matrix
    std::vector<Qe> products(1 << PAIRS);
    products[0] = qone();
    for (unsigned mask = 1; mask < products.size(); ++mask) {
        const unsigned bit = __builtin_ctz(mask);
        products[mask] = qmul(products[mask ^ (1u << bit)], H[bit]);
    }
    std::vector<Fe> mat((size_t)nrows * ncols, Fe{0, 0});
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t eq = 0; eq < (int64_t)six.size(); ++eq) {
        const std::vector<int> &sel = six[eq];
        for (int omit = 0; omit < R + 1; ++omit) {
            std::vector<int> pairs;
            for (int s = 0; s < R + 1; ++s)
                if (s != omit) pairs.push_back(sel[s]);
            const uint32_t block = five_index.at(code_of(pairs));
            for (int side = 0; side < NMASK; ++side) {
                unsigned hmask = 0;
                for (int k = 0; k < R; ++k)
                    if ((side >> k) & 1) hmask |= 1u << pairs[k];
                const Qe coeff = qmul(U[sel[omit]], products[hmask]);
                const uint32_t col = block * NMASK + side;
                for (int deg = 0; deg < M; ++deg)
                    mat[((size_t)eq * M + deg) * ncols + col] = coeff.c[deg];
            }
        }
    }
    const double t_build = elapsed(t0);
    if (getenv("SELFTEST")) {
        // field sanity: commutativity, associativity, inverse, mul_x == mul by y
        std::mt19937_64 g(7);
        auto rnd = [&]() { Fe a{g(), g()}; a.hi &= ((uint64_t)1 << (N - 64)) - 1; return a; };
        Fe y{2, 0};
        int bad = 0, b_comm = 0, b_assoc = 0, b_inv = 0, b_mulx = 0;
        for (int t = 0; t < 200; ++t) {
            Fe a = rnd(), b = rnd(), c = rnd();
            if (add(mul(a, b), mul(b, a)).lo || add(mul(a, b), mul(b, a)).hi) ++b_comm;
            if (add(mul(a, mul(b, c)), mul(mul(a, b), c)).lo ||
                add(mul(a, mul(b, c)), mul(mul(a, b), c)).hi) ++b_assoc;
            Fe ai = inverse(a);
            Fe one = mul(a, ai);
            if (one.lo != 1 || one.hi != 0) ++b_inv;
            Fe s1 = a; mul_x(s1);
            Fe s2 = mul(a, y);
            if (s1.lo != s2.lo || s1.hi != s2.hi) ++b_mulx;
            ++bad;
        }
        fprintf(stderr, "selftest comm=%d assoc=%d inv=%d mulx=%d\n", b_comm, b_assoc, b_inv, b_mulx);
    }
    std::vector<Fe> orig;

    orig = mat;
    // ---- kernel over L (blocked)
    std::vector<uint32_t> pivots, bcols, brows;
    const int nth = threads();
    std::vector<Fe> coef((size_t)nth * BLK);
    uint32_t rank = 0, col = 0;
    while (rank < nrows && col < ncols) {
        const uint32_t base = rank;
        bcols.clear();
        brows.clear();
        uint32_t c = col;
        while ((int)bcols.size() < BLK && rank < nrows && c < ncols) {
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
                for (size_t j = 0; j < bcols.size(); ++j) {
                    const Fe *pj = mat.data() + (size_t)brows[j] * ncols;
                    const Fe f = p[bcols[j]];
                    if (!nonzero(f)) continue;
                    p[bcols[j]] = Fe{0, 0};
                    for (uint32_t jj = bcols[j] + 1; jj < ncols; ++jj)
                        p[jj] = add(p[jj], mul(f, pj[jj]));
                }
                if (nonzero(p[c])) { accepted = true; break; }
                cand = rank + 1;
            }
            if (!accepted) { ++c; continue; }
            Fe *p = mat.data() + (size_t)rank * ncols;
            const Fe scale = inverse(p[c]);
            p[c] = Fe{1, 0};
            for (uint32_t j = c + 1; j < ncols; ++j) p[j] = mul(p[j], scale);
            for (uint32_t k = base; k < rank; ++k) {
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
        const uint32_t lo = bcols.front(), rr0 = rank;
#pragma omp parallel for schedule(static)
        for (int64_t rr = rr0; rr < (int64_t)nrows; ++rr) {
            Fe *row = mat.data() + (size_t)rr * ncols;
            Fe *cf = coef.data() + (size_t)thread_id() * BLK;
            bool any = false;
            for (int i = 0; i < nb; ++i) {
                cf[i] = row[bcols[i]];
                if (nonzero(cf[i])) any = true;
            }
            if (!any) continue;
            uint32_t j = lo;
#if defined(__AVX512F__) && defined(__VPCLMULQDQ__) && defined(__AVX512BW__)
            {
                u128 lo4[4], hi4[4];
                for (; j + 4 <= ncols; j += 4) {
                    __m512i accA = _mm512_loadu_si512((const void *)(row + j));
                    __m512i accB = _mm512_setzero_si512();
                    for (int i = 0; i < nb; ++i) {
                        if (!nonzero(cf[i])) continue;
                        const __m512i a = _mm512_broadcast_i32x4(
                            _mm_set_epi64x((long long)cf[i].hi, (long long)cf[i].lo));
                        const __m512i b = _mm512_loadu_si512(
                            (const void *)(mat.data() + (size_t)brows[i] * ncols + j));
                        const __m512i p00 = _mm512_clmulepi64_epi128(a, b, 0x00);
                        const __m512i cr = _mm512_xor_si512(
                            _mm512_clmulepi64_epi128(a, b, 0x01),
                            _mm512_clmulepi64_epi128(a, b, 0x10));
                        const __m512i p11 = _mm512_clmulepi64_epi128(a, b, 0x11);
                        accA = _mm512_xor_si512(
                            accA, _mm512_xor_si512(p00, _mm512_bslli_epi128(cr, 8)));
                        accB = _mm512_xor_si512(
                            accB, _mm512_xor_si512(_mm512_bsrli_epi128(cr, 8), p11));
                    }
                    _mm512_storeu_si512((void *)lo4, accA);
                    _mm512_storeu_si512((void *)hi4, accB);
                    for (int t = 0; t < 4; ++t) row[j + t] = reduce(lo4[t], hi4[t]);
                }
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
        col = c;
    }
    const double t_kernel = elapsed(t0);
    std::vector<uint32_t> is_pivot(ncols, 0);
    for (uint32_t p : pivots) is_pivot[p] = 1;
    std::vector<uint32_t> free_cols;
    for (uint32_t j = 0; j < ncols; ++j)
        if (!is_pivot[j]) free_cols.push_back(j);
    const size_t dim = free_cols.size();
    std::vector<Fe> basis(dim * ncols, Fe{0, 0});
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t fi = 0; fi < (int64_t)dim; ++fi) {
        const uint32_t f = free_cols[fi];
        basis[(size_t)fi * ncols + f] = Fe{1, 0};
        for (int64_t i = (int64_t)pivots.size() - 1; i >= 0; --i) {
            const Fe *p = mat.data() + (size_t)i * ncols;
            Fe acc = (f < pivots[i]) ? p[f] : Fe{0, 0};
            for (uint32_t j = pivots[i] + 1; j < ncols; ++j) {
                const Fe v = basis[(size_t)fi * ncols + j];
                if (nonzero(v) && nonzero(p[j])) acc = add(acc, mul(p[j], v));
            }
            basis[(size_t)fi * ncols + pivots[i]] = acc;
        }
    }
    mat.clear();
    mat.shrink_to_fit();

    // ---- Frobenius relations (binary)
    std::vector<std::pair<uint32_t, uint32_t>> relations;
    for (uint32_t ci = 0; ci < five.size(); ++ci) {
        if (five[ci].back() == PAIRS - 1) continue;
        std::vector<int> sh = five[ci];
        for (int &x : sh) ++x;
        const uint32_t si = five_index.at(code_of(sh));
        for (int mask = 0; mask < NMASK; ++mask)
            relations.emplace_back(ci * NMASK + mask, si * NMASK + mask);
    }
    size_t nrel = std::min(maxrel, relations.size());
    const size_t unknowns = dim * N;
    const size_t words = (unknowns + 63) / 64 + 2;
    auto build_bin = [&](size_t nr) {
        const size_t brows2 = nr * N;
        std::vector<uint64_t> bin(brows2 * words, 0);
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t r = 0; r < (int64_t)nr; ++r) {
            std::vector<Fe> v(N);
            for (size_t f = 0; f < dim; ++f) {
                Fe l = basis[f * ncols + relations[r].second];
                Fe rt = square(basis[f * ncols + relations[r].first]);
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
                    uint64_t *row = bin.data() + ((size_t)r * N + b) * words;
                    const size_t off = f * N, blk = off >> 6;
                    const unsigned sh = off & 63;
                    if (sh == 0) {
                        row[blk] ^= lo;
                        row[blk + 1] ^= hi;
                    } else {
                        row[blk] ^= lo << sh;
                        row[blk + 1] ^= (lo >> (64 - sh)) ^ (hi << sh);
                        row[blk + 2] ^= hi >> (64 - sh);
                    }
                }
            }
        }
        return bin;
    };
    auto bin_kernel = [&](std::vector<uint64_t> &rows, size_t nrows2) {
        std::vector<uint32_t> pv;
        size_t rk = 0;
        for (size_t c = 0; c < unknowns && rk < nrows2; ++c) {
            size_t found = rk;
            while (found < nrows2 && !((rows[found * words + (c >> 6)] >> (c & 63)) & 1)) ++found;
            if (found == nrows2) continue;
            if (found != rk)
                for (size_t w = 0; w < words; ++w)
                    std::swap(rows[found * words + w], rows[rk * words + w]);
            const uint64_t *p = rows.data() + rk * words;
            const size_t w0 = c >> 6;
#pragma omp parallel for schedule(static)
            for (int64_t ii = (int64_t)rk + 1; ii < (int64_t)nrows2; ++ii) {
                uint64_t *r = rows.data() + (size_t)ii * words;
                if (!((r[w0] >> (c & 63)) & 1)) continue;
                for (size_t w = w0; w < words; ++w) r[w] ^= p[w];
            }
            pv.push_back((uint32_t)c);
            ++rk;
        }
        std::vector<uint32_t> isp(unknowns, 0);
        for (uint32_t x : pv) isp[x] = 1;
        std::vector<std::vector<uint64_t>> out;
        for (size_t f = 0; f < unknowns; ++f) {
            if (isp[f]) continue;
            std::vector<uint64_t> x(words, 0);
            x[f >> 6] |= uint64_t(1) << (f & 63);
            for (int64_t i = (int64_t)pv.size() - 1; i >= 0; --i) {
                const uint64_t *p = rows.data() + (size_t)i * words;
                uint64_t acc = (p[f >> 6] >> (f & 63)) & 1;
                for (size_t gj = (size_t)i + 1; gj < pv.size(); ++gj) {
                    const uint32_t pc = pv[gj];
                    if ((x[pc >> 6] >> (pc & 63)) & 1) acc ^= (p[pc >> 6] >> (pc & 63)) & 1;
                }
                if (acc) x[pv[i] >> 6] ^= uint64_t(1) << (pv[i] & 63);
            }
            out.push_back(std::move(x));
        }
        return out;
    };
    std::vector<uint64_t> bin = build_bin(nrel);
    std::vector<std::vector<uint64_t>> kbin = bin_kernel(bin, nrel * N);
    if (kbin.size() != 1 && nrel < relations.size()) {
        nrel = relations.size();
        bin = build_bin(nrel);
        kbin = bin_kernel(bin, nrel * N);
    }
    const double t_bin = elapsed(t0);
    if (getenv("KCHECK")) {
        long bad = 0;
#pragma omp parallel for reduction(+ : bad) schedule(dynamic, 1)
        for (int64_t f = 0; f < (int64_t)dim; ++f)
            for (uint32_t i = 0; i < nrows; ++i) {
                Fe v{0, 0};
                const Fe *row = &orig[(size_t)i * ncols];
                for (uint32_t j = 0; j < ncols; ++j)
                    if (nonzero(row[j]) && nonzero(basis[f * ncols + j]))
                        v = add(v, mul(row[j], basis[f * ncols + j]));
                bad += nonzero(v);
            }
        fprintf(stderr, "kernel check bad=%ld of %zu\n", bad, dim * nrows);
    }
    orig.clear();
    orig.shrink_to_fit();

    printf("{\"m\":%d,\"n\":%d,\"r\":%d,\"p\":%d,\"matrix\":[%u,%u],\"kernel_dim\":%zu,"
           "\"binary\":[%zu,%zu],\"bin_kernel_dim\":%zu",
           M, N, R, PAIRS, nrows, ncols, dim, nrel * N, unknowns, kbin.size());
    if (kbin.size() != 1) {
        printf(",\"status\":\"kernel_dim!=1\",\"t_build\":%.2f,\"t_kernel\":%.2f,\"t_bin\":%.2f}\n",
               t_build, t_kernel, t_bin);
        return 1;
    }

    // ---- Plucker vector and generator reconstruction
    std::vector<Fe> pv(ncols, Fe{0, 0});
    std::vector<Fe> lambda(dim, Fe{0, 0});
    for (size_t f = 0; f < dim; ++f) {
        uint64_t lo = 0, hi = 0;
        for (int b = 0; b < N; ++b)
            if ((kbin[0][(f * N + b) >> 6] >> ((f * N + b) & 63)) & 1) {
                if (b < 64) lo |= uint64_t(1) << b;
                else hi |= uint64_t(1) << (b - 64);
            }
        lambda[f] = Fe{lo, hi};
    }
#pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < (int64_t)ncols; ++j) {
        Fe v{0, 0};
        for (size_t f = 0; f < dim; ++f) {
            const Fe b = basis[f * ncols + j];
            if (nonzero(lambda[f]) && nonzero(b)) v = add(v, mul(lambda[f], b));
        }
        pv[j] = v;
    }
    auto coord = [&](std::vector<std::pair<int, int>> sel) {
        std::sort(sel.begin(), sel.end());
        std::vector<int> s;
        int mask = 0;
        for (int i = 0; i < R; ++i) {
            s.push_back(sel[i].first);
            if (sel[i].second) mask |= 1 << i;
        }
        return (size_t)five_index.at(code_of(s)) * NMASK + mask;
    };
    std::vector<std::pair<int, int>> chart;
    for (int i = 0; i < R; ++i) chart.push_back({i, 1});
    const Fe chart_minor = pv[coord(chart)];
    if (!nonzero(chart_minor)) { printf(",\"status\":\"singular chart\"}\n"); return 1; }
    const Fe inv_chart = inverse(chart_minor);
    std::vector<std::vector<Fe>> G(R, std::vector<Fe>(2 * PAIRS, Fe{0, 0}));
    for (int i = 0; i < R; ++i) G[i][2 * i + 1] = Fe{1, 0};
    for (int pair = R; pair < PAIRS; ++pair)
        for (int side = 0; side < 2; ++side)
            for (int row = 0; row < R; ++row) {
                std::vector<std::pair<int, int>> sel;
                for (int i = 0; i < R; ++i)
                    if (i != row) sel.push_back({i, 1});
                sel.push_back({pair, side});
                G[row][2 * pair + side] = mul(pv[coord(sel)], inv_chart);
            }
    for (int i = 0; i < R; ++i) {
        std::vector<std::pair<int, int>> sel;
        for (int j = 0; j < R; ++j) sel.push_back({j, j == i ? 0 : 1});
        G[i][2 * i] = mul(pv[coord(sel)], inv_chart);
    }
    bool ok = true;
    for (int i = 0; i < R && ok; ++i) {
        int pair = -1, side = -1;
        for (int p = R; p < PAIRS && pair < 0; ++p)
            for (int s = 0; s < 2; ++s)
                if (nonzero(G[i][2 * p + s])) { pair = p; side = s; break; }
        if (pair < 0) { ok = false; break; }
        for (int row = 0; row < R; ++row) {
            if (row == i) continue;
            std::vector<std::pair<int, int>> sel;
            for (int j = 0; j < R; ++j)
                if (j != i && j != row) sel.push_back({j, 1});
            sel.push_back({i, 0});
            sel.push_back({pair, side});
            const Fe d = mul(pv[coord(sel)], inv_chart);
            const Fe num = add(d, mul(G[i][2 * i], G[row][2 * pair + side]));
            G[row][2 * i] = mul(num, inverse(G[i][2 * pair + side]));
        }
    }
    if (!ok) { printf(",\"status\":\"no alternate chart\"}\n"); return 1; }

    // ---- solve for the a_k
    std::vector<std::vector<Qe>> A(PAIRS, std::vector<Qe>(R));
    for (int q = 0; q < PAIRS; ++q)
        for (int k = 0; k < R; ++k) {
            Qe g0 = qzero(), g1 = qzero();
            g0.c[0] = G[k][2 * q];
            g1.c[0] = G[k][2 * q + 1];
            A[q][k] = qadd(g0, qmul(H[q], g1));
        }
    std::vector<Qe> rhs(U.begin(), U.end());
    std::vector<int> prow(R, -1);
    int rk = 0;
    for (int c = 0; c < R && rk < PAIRS; ++c) {
        int piv = -1;
        for (int r = rk; r < PAIRS; ++r)
            if (!qiszero(A[r][c])) { piv = r; break; }
        if (piv < 0) continue;
        std::swap(A[rk], A[piv]);
        std::swap(rhs[rk], rhs[piv]);
        const Qe inv = qinv(A[rk][c]);
        for (int j = 0; j < R; ++j) A[rk][j] = qmul(A[rk][j], inv);
        rhs[rk] = qmul(rhs[rk], inv);
        for (int r = 0; r < PAIRS; ++r) {
            if (r == rk || qiszero(A[r][c])) continue;
            const Qe f = A[r][c];
            for (int j = 0; j < R; ++j) A[r][j] = qadd(A[r][j], qmul(f, A[rk][j]));
            rhs[r] = qadd(rhs[r], qmul(f, rhs[rk]));
        }
        prow[c] = rk;
        ++rk;
    }
    if (rk != R) { printf(",\"status\":\"singular system\"}\n"); return 1; }
    Qe B = qzero();
    for (int k = 0; k < R; ++k) {
        Qe g1 = qzero();
        g1.c[0] = G[k][1];
        B = qadd(B, qmul(rhs[prow[k]], g1));
    }
    // ---- support of the recovered e2
    std::vector<u128> e2(N, 0);
    for (int bit = 0; bit < M; ++bit) {
        const uint64_t lo = B.c[bit].lo, hi = B.c[bit].hi;
        for (int j = 0; j < 64; ++j)
            if ((lo >> j) & 1) e2[j] |= (u128)1 << bit;
        for (int j = 0; j < N - 64; ++j)
            if ((hi >> j) & 1) e2[64 + j] |= (u128)1 << bit;
    }
    auto rref = [&](std::vector<u128> v) {
        int r2 = 0;
        for (int c = 0; c < M && r2 < (int)v.size(); ++c) {
            int piv = -1;
            for (int i = r2; i < (int)v.size(); ++i)
                if ((v[i] >> c) & 1) { piv = i; break; }
            if (piv < 0) continue;
            std::swap(v[r2], v[piv]);
            for (int i = 0; i < (int)v.size(); ++i)
                if (i != r2 && ((v[i] >> c) & 1)) v[i] ^= v[r2];
            ++r2;
        }
        v.resize(r2);
        return v;
    };
    const std::vector<u128> got = rref(e2), want = rref(a);
    const bool success = (got.size() == (size_t)R) && (got == want);
    const double t_total = elapsed(t0);
    printf(",\"status\":\"%s\",\"t_build\":%.2f,\"t_kernel\":%.2f,\"t_bin\":%.2f,\"t_total\":%.2f,"
           "\"threads\":%d}\n",
           success ? "RECOVERED" : "FAILED", t_build, t_kernel, t_bin, t_total, threads());
    return success ? 0 : 1;
}
