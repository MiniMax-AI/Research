// tsm.cpp -- Transposed Support-Minor (TSM) attack on ideal-LRPC KEMs
// (ROLLO-I / LOCKER-I family), runtime-parameterised over (m, n, r, P).
//
// Algebraic model:
//   R  = GF(2^m)[X]/(f(X)), deg f = n,  h = x^-1 y,  u = e1 + e2 h,
//   supp(e1) = supp(e2) = E,  dim_{GF(2)} E = r.
//   Transposition Phi: R -> Q = GF(2^n)[T]/(g(T)), g the GF(2^m) modulus; valid
//   because gcd(m,n)=1 makes Q the field GF(2^(mn)).  Then
//   Phi(e1) = sum_k a_k E~_k, Phi(e2) = sum_k b_k E~_k with a_k,b_k in GF(2^n)
//   and E~_k in GF(2^m) subset Q, fixed by sigma (Frobenius of GF(2^n)).
//   Applying sigma^q for q = 0..P-1 gives P equations
//       U_q = sum_k E~_k ( A_{q,k} + H_q B_{q,k} ),
//   whose support-minor linearisation (maximal minors of the r x 2P matrix
//   G = [A_0 B_0 ... A_{P-1} B_{P-1}]) has C(P,r)*2^r unknowns over GF(2^n)
//   and C(P-1,r)*m equations; the residual kernel is killed by the Frobenius
//   relations p_{Q+1} = sigma(p_Q), which are GF(2)-linear in the kernel
//   coordinates.
//
// usage:
//   tsm --m 71 --n 83 --r 5 --pairs 9 --ring 83,7,4,2,0 --field 71,6,0 \
//       --pk HEX --ct HEX --out PLUCKER.bin [--stats]
//
// build: g++ -O3 -march=native -fopenmp -std=c++17 tsm.cpp \
//            $(pkg-config --cflags m4ri) -o tsm $(pkg-config --libs m4ri)

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#include <immintrin.h>
#endif
#include <m4ri/m4ri.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
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

static int N_BITS = 83;      // degree of the base field GF(2^n)
static int N_WORDS = 2;
static const int MAXW = 4;   // supports n <= 256
struct Fe { uint64_t w[MAXW]; };

static std::vector<int> LOWPOLY;              // exponents < n of x^n + sum x^low
static std::vector<Fe> XRED;                 // XRED[e-n] = x^e mod modulus, n<=e<=2n-2
static int RED_CHUNKS = 0;
static std::vector<Fe> REDTAB;

static inline Fe zero_fe() { Fe a; for (int i = 0; i < MAXW; ++i) a.w[i] = 0; return a; }
static inline bool nonzero(const Fe &a) { uint64_t z = 0; for (int i = 0; i < N_WORDS; ++i) z |= a.w[i]; return z != 0; }
static inline Fe add(const Fe &a, const Fe &b) { Fe c; for (int i = 0; i < N_WORDS; ++i) c.w[i] = a.w[i] ^ b.w[i]; return c; }
static inline void xor_into(Fe &a, const Fe &b) { for (int i = 0; i < N_WORDS; ++i) a.w[i] ^= b.w[i]; }
static inline void mask_top(Fe &a) { if (N_BITS & 63) a.w[N_WORDS - 1] &= (uint64_t(1) << (N_BITS & 63)) - 1; }

static void build_tables() {
    const int top = 2 * N_BITS - 2;
    XRED.assign(top - N_BITS + 1, zero_fe());
    for (int e = N_BITS; e <= top; ++e) {
        Fe v = zero_fe();
        v.w[e >> 6] = uint64_t(1) << (e & 63);           // start from x^e
        for (int b = e; b >= N_BITS; --b) {
            if (!((v.w[b >> 6] >> (b & 63)) & 1)) continue;
            v.w[b >> 6] &= ~(uint64_t(1) << (b & 63));
            for (int j : LOWPOLY) {
                int pos = b - N_BITS + j;
                v.w[pos >> 6] ^= uint64_t(1) << (pos & 63);
            }
        }
        XRED[e - N_BITS] = v;
    }
    RED_CHUNKS = (N_BITS - 1 + 7) / 8;
    REDTAB.assign(size_t(RED_CHUNKS) * 256, zero_fe());
    for (int c = 0; c < RED_CHUNKS; ++c)
        for (int v = 1; v < 256; ++v) {
            Fe acc = zero_fe();
            for (int k = 0; k < 8; ++k) {
                if (!((v >> k) & 1)) continue;
                int e = N_BITS + 8 * c + k;
                if (e > top) continue;
                xor_into(acc, XRED[e - N_BITS]);
            }
            REDTAB[size_t(c) * 256 + v] = acc;
        }
}

static inline void clmul64(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
#if defined(__PCLMUL__)
    const __m128i p = _mm_clmulepi64_si128(_mm_set_epi64x(0, (long long)a),
                                           _mm_set_epi64x(0, (long long)b), 0x00);
    lo = (uint64_t)_mm_cvtsi128_si64(p);
    hi = (uint64_t)_mm_extract_epi64(p, 1);
#else
    uint64_t r0 = 0, r1 = 0;
    for (int i = 0; i < 64; ++i)
        if ((b >> i) & 1) {
            r0 ^= a << i;
            if (i) r1 ^= a >> (64 - i);
        }
    lo = r0;
    hi = r1;
#endif
}

__attribute__((always_inline)) static inline Fe mul(const Fe &a, const Fe &b) {
    uint64_t p[2 * MAXW];
    for (int i = 0; i < 2 * MAXW; ++i) p[i] = 0;
    for (int i = 0; i < N_WORDS; ++i) {
        if (!a.w[i]) continue;
        for (int j = 0; j < N_WORDS; ++j) {
            if (!b.w[j]) continue;
            uint64_t lo, hi;
            clmul64(a.w[i], b.w[j], lo, hi);
            p[i + j]     ^= lo;
            p[i + j + 1] ^= hi;
        }
    }
    Fe out = zero_fe();
    for (int i = 0; i < N_WORDS; ++i) out.w[i] = p[i];
    mask_top(out);
    for (int c = RED_CHUNKS - 1; c >= 0; --c) {
        const int base = N_BITS + 8 * c;
        const int widx = base >> 6, off = base & 63;
        uint64_t lo = p[widx] >> off;
        if (off) lo |= p[widx + 1] << (64 - off);
        const uint32_t byte = uint32_t(lo) & 0xFFu;
        if (!byte) continue;
        xor_into(out, REDTAB[size_t(c) * 256 + byte]);
    }
    return out;
}

static inline Fe square(const Fe &a) { return mul(a, a); }

static inline Fe mul_x(const Fe &a) {
    Fe z = zero_fe();
    uint64_t carry = 0;
    for (int i = 0; i < N_WORDS; ++i) {
        z.w[i] = (a.w[i] << 1) | carry;
        carry = a.w[i] >> 63;
    }
    if ((N_BITS & 63) && (z.w[N_WORDS - 1] >> (N_BITS & 63)) & 1) {
        // the bit that left the field is at position N_BITS
        for (int j : LOWPOLY) z.w[j >> 6] ^= uint64_t(1) << (j & 63);
    } else if (!(N_BITS & 63) && carry) {
        for (int j : LOWPOLY) z.w[j >> 6] ^= uint64_t(1) << (j & 63);
    }
    mask_top(z);
    return z;
}

static Fe inverse(const Fe &a) {
    if (!nonzero(a)) throw std::runtime_error("inverse of zero");
    Fe result = zero_fe(); result.w[0] = 1;
    Fe power = a;
    for (int bit = 0; bit < N_BITS; ++bit) {
        if (bit) result = mul(result, power);
        power = square(power);
    }
    return result;
}

// ------------------------------------------------------------ GF(2^n)[T]/(g) --
static int M_BITS = 71;
static std::vector<int> GLOW;
struct Qe { std::vector<Fe> c; };

static Qe qzero() { Qe q; q.c.assign(M_BITS, zero_fe()); return q; }
static Qe qone() { Qe q = qzero(); q.c[0].w[0] = 1; return q; }
static Qe qmul(const Qe &a, const Qe &b) {
    std::vector<Fe> p(size_t(2) * M_BITS - 1, zero_fe());
    for (int i = 0; i < M_BITS; ++i) {
        if (!nonzero(a.c[i])) continue;
        for (int j = 0; j < M_BITS; ++j) {
            if (!nonzero(b.c[j])) continue;
            xor_into(p[i + j], mul(a.c[i], b.c[j]));
        }
    }
    for (int d = 2 * M_BITS - 2; d >= M_BITS; --d) {
        if (!nonzero(p[d])) continue;
        for (int j : GLOW) xor_into(p[d - M_BITS + j], p[d]);
    }
    Qe out; out.c.assign(p.begin(), p.begin() + M_BITS); return out;
}
static Qe qrho(const Qe &a) { Qe o; o.c.resize(M_BITS); for (int i = 0; i < M_BITS; ++i) o.c[i] = square(a.c[i]); return o; }

// ------------------------------------------------------------------ helpers --
static double elapsed(const std::chrono::steady_clock::time_point &s) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - s).count();
}
static std::vector<unsigned char> unhex(std::string s) {
    if (s.rfind("0x", 0) == 0 || s.rfind("0X", 0) == 0) s.erase(0, 2);
    std::vector<unsigned char> o(s.size() / 2);
    auto d = [](char c) { return c <= '9' ? c - '0' : (c | 32) - 'a' + 10; };
    for (size_t i = 0; i < o.size(); ++i) o[i] = (unsigned char)((d(s[2*i]) << 4) | d(s[2*i+1]));
    return o;
}
static std::vector<int> parse_list(const std::string &s) {
    std::vector<int> v; size_t i = 0;
    while (i <= s.size()) { size_t j = s.find(',', i); if (j == std::string::npos) j = s.size();
        if (j > i) v.push_back(atoi(s.substr(i, j - i).c_str())); if (j == s.size()) break; i = j + 1; }
    return v;
}
static void combos_rec(int next, int need, int total, std::vector<int> &cur, std::vector<std::vector<int>> &out) {
    if (!need) { out.push_back(cur); return; }
    for (int x = next; x <= total - need; ++x) { cur.push_back(x); combos_rec(x + 1, need - 1, total, cur, out); cur.pop_back(); }
}
static uint32_t combo_code(const std::vector<int> &c) { uint32_t z = 0; for (size_t i = 0; i < c.size(); ++i) z |= uint32_t(c[i]) << (5 * i); return z; }

static Qe transpose_ring(const std::vector<unsigned char> &raw, int nelem, int mbytes) {
    Qe out = qzero();
    for (int j = 0; j < nelem; ++j)
        for (int bit = 0; bit < M_BITS; ++bit)
            if ((raw[mbytes * j + (bit >> 3)] >> (bit & 7)) & 1)
                out.c[bit].w[j >> 6] ^= uint64_t(1) << (j & 63);
    return out;
}

// xor an n-bit field element into a bit-row at bit offset 'offset'
static inline void xor_field_at(word *row, size_t offset, const Fe &a) {
    size_t block = offset >> 6;
    const unsigned shift = unsigned(offset & 63);
    if (shift == 0) { for (int i = 0; i < N_WORDS; ++i) row[block + i] ^= a.w[i]; return; }
    for (int i = 0; i < N_WORDS; ++i) {
        row[block + i]     ^= a.w[i] << shift;
        row[block + i + 1] ^= a.w[i] >> (64 - shift);
    }
}

int main(int argc, char **argv) {
    const auto started = std::chrono::steady_clock::now();
    int R = 5, PAIRS = 9, NELEM = 0, MBYTES = 0;
    std::string pk_hex, ct_hex, out_path = "plucker.bin";
    std::vector<int> ring, field;
    bool stats = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--m") M_BITS = atoi(nx().c_str());
        else if (a == "--n") N_BITS = atoi(nx().c_str());
        else if (a == "--r") R = atoi(nx().c_str());
        else if (a == "--pairs") PAIRS = atoi(nx().c_str());
        else if (a == "--ring") ring = parse_list(nx());
        else if (a == "--field") field = parse_list(nx());
        else if (a == "--pk") pk_hex = nx();
        else if (a == "--ct") ct_hex = nx();
        else if (a == "--out") out_path = nx();
        else if (a == "--stats") stats = true;
        else { std::cerr << "unknown option " << a << "\n"; return 2; }
    }
    NELEM = ring.empty() ? N_BITS : ring[0];
    N_WORDS = (N_BITS + 63) / 64;
    if (N_WORDS > MAXW) { std::cerr << "n too large\n"; return 2; }
    MBYTES = (M_BITS + 7) / 8;
    for (size_t i = 1; i < ring.size(); ++i) LOWPOLY.push_back(ring[i]);
    for (size_t i = 1; i < field.size(); ++i) GLOW.push_back(field[i]);
    build_tables();

    const auto pk = unhex(pk_hex);
    const auto ct = unhex(ct_hex);
    if (int(pk.size()) < NELEM * MBYTES || int(ct.size()) < NELEM * MBYTES)
        throw std::runtime_error("input too short for the declared parameters");

    std::vector<Qe> H(PAIRS), U(PAIRS);
    H[0] = transpose_ring(pk, NELEM, MBYTES);
    U[0] = transpose_ring(ct, NELEM, MBYTES);
    for (int q = 1; q < PAIRS; ++q) { H[q] = qrho(H[q-1]); U[q] = qrho(U[q-1]); }
    std::vector<Qe> products(size_t(1) << PAIRS, qzero());
    products[0] = qone();
    for (unsigned mask = 1; mask < (1u << PAIRS); ++mask) {
        const unsigned bit = __builtin_ctz(mask);
        products[mask] = qmul(products[mask ^ (1u << bit)], H[bit]);
    }
    std::cerr << "products elapsed=" << elapsed(started) << "\n";

    std::vector<std::vector<int>> fives, sixes_all, sixes;
    { std::vector<int> cur; combos_rec(0, R, PAIRS, cur, fives); }
    { std::vector<int> cur; combos_rec(0, R + 1, PAIRS, cur, sixes_all); }
    for (auto &c : sixes_all) if (c[0] == 0) sixes.push_back(c);
    std::map<uint32_t, uint32_t> five_index;
    for (uint32_t i = 0; i < fives.size(); ++i) five_index[combo_code(fives[i])] = i;

    const uint32_t ncols = uint32_t(fives.size()) << R;
    const uint32_t nrows = uint32_t(sixes.size()) * M_BITS;
    std::vector<Fe> matrix(size_t(nrows) * ncols, zero_fe());

#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t eq = 0; eq < int64_t(sixes.size()); ++eq) {
        const std::vector<int> &sel = sixes[eq];
        for (int omit = 0; omit < R + 1; ++omit) {
            std::vector<int> pairs;
            for (int s = 0; s < R + 1; ++s) if (s != omit) pairs.push_back(sel[s]);
            const uint32_t block = five_index.at(combo_code(pairs));
            for (unsigned sm = 0; sm < (1u << R); ++sm) {
                unsigned hmask = 0;
                for (int k = 0; k < R; ++k) if ((sm >> k) & 1) hmask |= 1u << pairs[k];
                const Qe coef = qmul(U[sel[omit]], products[hmask]);
                const uint32_t col = block * (1u << R) + sm;
                for (int d = 0; d < M_BITS; ++d)
                    matrix[(size_t(eq) * M_BITS + d) * ncols + col] = coef.c[d];
            }
        }
    }
    std::cerr << "matrix " << nrows << "x" << ncols << " elapsed=" << elapsed(started) << "\n";

    // ---------------------------------------------- echelon over GF(2^n) --
    std::vector<uint32_t> order(nrows); for (uint32_t i = 0; i < nrows; ++i) order[i] = i;
    std::vector<uint32_t> pivots; pivots.reserve(std::min(nrows, ncols));
    uint32_t rank = 0;
    for (uint32_t col = 0; col < ncols && rank < nrows; ++col) {
        uint32_t found = rank;
        while (found < nrows && !nonzero(matrix[size_t(order[found]) * ncols + col])) ++found;
        if (found == nrows) continue;
        std::swap(order[rank], order[found]);
        Fe *pivot = matrix.data() + size_t(order[rank]) * ncols;
        const Fe sc = inverse(pivot[col]);
        pivot[col] = zero_fe(); pivot[col].w[0] = 1;
        for (uint32_t j = col + 1; j < ncols; ++j) pivot[j] = mul(pivot[j], sc);
#pragma omp parallel for schedule(static)
        for (int64_t ii = int64_t(rank) + 1; ii < int64_t(nrows); ++ii) {
            Fe *row = matrix.data() + size_t(order[ii]) * ncols;
            const Fe factor = row[col];
            if (!nonzero(factor)) continue;
            row[col] = zero_fe();
            for (uint32_t j = col + 1; j < ncols; ++j)
                if (nonzero(pivot[j])) xor_into(row[j], mul(factor, pivot[j]));
        }
        pivots.push_back(col); ++rank;
    }
    const uint32_t nullity = ncols - rank;
    std::cerr << "echelon rank=" << rank << " nullity=" << nullity << " elapsed=" << elapsed(started) << "\n";
    if (!nullity) { std::cerr << "trivial kernel\n"; return 1; }

    std::vector<unsigned char> is_pivot(ncols, 0);
    for (uint32_t c : pivots) is_pivot[c] = 1;
    std::vector<uint32_t> free_cols;
    for (uint32_t c = 0; c < ncols; ++c) if (!is_pivot[c]) free_cols.push_back(c);
    std::vector<Fe> basis(size_t(nullity) * ncols, zero_fe());
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t f = 0; f < int64_t(nullity); ++f) {
        Fe *ker = basis.data() + size_t(f) * ncols;
        ker[free_cols[f]].w[0] = 1;
        for (int64_t ii = int64_t(rank) - 1; ii >= 0; --ii) {
            const uint32_t col = pivots[ii];
            const Fe *row = matrix.data() + size_t(order[ii]) * ncols;
            Fe value = zero_fe();
            for (uint32_t j = col + 1; j < ncols; ++j)
                if (nonzero(ker[j]) && nonzero(row[j])) xor_into(value, mul(row[j], ker[j]));
            ker[col] = value;
        }
    }
    std::cerr << "kernel basis elapsed=" << elapsed(started) << "\n";
    matrix.clear(); matrix.shrink_to_fit();

    // ------------------------------- Frobenius relations p_{Q+1} = sigma(p_Q) --
    std::vector<std::pair<uint32_t, uint32_t>> rel;
    for (uint32_t ci = 0; ci < fives.size(); ++ci) {
        if (fives[ci].back() == PAIRS - 1) continue;
        std::vector<int> shifted = fives[ci];
        for (int &x : shifted) ++x;
        const uint32_t si = five_index.at(combo_code(shifted));
        for (unsigned sm = 0; sm < (1u << R); ++sm)
            rel.emplace_back(ci * (1u << R) + sm, si * (1u << R) + sm);
    }
    const size_t nrel = rel.size();
    const size_t UBits = size_t(nullity) * N_BITS;
    const size_t EqBits = nrel * N_BITS;
    std::cerr << "relations=" << nrel << " frobenius system " << EqBits << "x" << UBits
              << " elapsed=" << elapsed(started) << "\n";

    mzd_t *trans = mzd_init(UBits, EqBits);
#pragma omp parallel for schedule(static)
    for (int64_t f = 0; f < int64_t(nullity); ++f) {
        std::vector<Fe> left(nrel), right(nrel);
        for (size_t r = 0; r < nrel; ++r) {
            left[r]  = basis[size_t(f) * ncols + rel[r].second];
            right[r] = square(basis[size_t(f) * ncols + rel[r].first]);
        }
        for (int bit = 0; bit < N_BITS; ++bit) {
            word *row = mzd_row(trans, rci_t(size_t(f) * N_BITS + bit));
            for (size_t r = 0; r < nrel; ++r)
                xor_field_at(row, r * N_BITS, add(left[r], right[r]));
            for (size_t r = 0; r < nrel; ++r) { left[r] = mul_x(left[r]); right[r] = mul_x(mul_x(right[r])); }
        }
    }
    std::cerr << "frobenius transpose built elapsed=" << elapsed(started) << "\n";
    mzd_t *system = mzd_transpose(nullptr, trans);
    mzd_free(trans);
    mzd_t *kernel = mzd_kernel_left_pluq(system, 0);
    mzd_free(system);
    if (!kernel) { std::cerr << "frobenius system has trivial kernel\n"; return 3; }
    std::cerr << "frobenius kernel " << kernel->nrows << "x" << kernel->ncols
              << " elapsed=" << elapsed(started) << "\n";
    if (kernel->ncols != 1) { std::cerr << "frobenius kernel dimension != 1\n"; mzd_free(kernel); return 3; }

    std::vector<Fe> pluck(ncols, zero_fe());
    {
        std::vector<Fe> params(nullity, zero_fe());
        for (uint32_t f = 0; f < nullity; ++f)
            for (int bit = 0; bit < N_BITS; ++bit)
                if (mzd_read_bit(kernel, rci_t(size_t(f) * N_BITS + bit), 0))
                    params[f].w[bit >> 6] ^= uint64_t(1) << (bit & 63);
        mzd_free(kernel);
#pragma omp parallel for schedule(static)
        for (int64_t j = 0; j < int64_t(ncols); ++j) {
            Fe value = zero_fe();
            for (uint32_t f = 0; f < nullity; ++f)
                if (nonzero(params[f]) && nonzero(basis[size_t(f) * ncols + j]))
                    xor_into(value, mul(params[f], basis[size_t(f) * ncols + j]));
            pluck[j] = value;
        }
    }

    std::ofstream out(out_path, std::ios::binary);
    const uint32_t hdr[4] = {ncols, uint32_t(R), uint32_t(PAIRS), uint32_t(N_WORDS)};
    out.write((const char *)hdr, sizeof hdr);
    out.write((const char *)pluck.data(), std::streamsize(pluck.size() * sizeof(Fe)));
    if (!out) throw std::runtime_error("failed writing Pluecker vector");
    if (stats)
        std::cerr << "{\"pairs\":" << PAIRS << ",\"nrows\":" << nrows << ",\"ncols\":" << ncols
                  << ",\"rank\":" << rank << ",\"nullity\":" << nullity
                  << ",\"relations\":" << nrel << ",\"seconds\":" << elapsed(started) << "}\n";
    std::cerr << "total elapsed=" << elapsed(started) << "\n";
    return 0;
}
