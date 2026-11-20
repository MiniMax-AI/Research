// Rank census over a small space of symmetric matrices.
//
// usage: rankmask <file> <N> <NB> [threads]
//
// The file holds NB symmetric N x N matrices over F_2, row-major, each row
// packed into ceil(N/64) little-endian 64-bit words (bit j of the row lives in
// word j/64 at bit position j%64).
//
// Prints "<mask> <rank>" for every nonzero mask in [1, 2^NB).  The hidden
// low-rank locus of the polar quotient is an 8-dimensional subspace whose
// members have rank about m + n'/2, far below the rank of a generic member of
// the quotient, so an exhaustive census plus a gap search identifies it without
// any structural guesswork.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: rankmask <file> <N> <NB> [threads]\n");
        return 2;
    }
    const int N = std::atoi(argv[2]);
    const int NB = std::atoi(argv[3]);
    const int NW = (N + 63) / 64;
    if (N <= 0 || NB <= 0 || NB > 24) return 2;

    std::vector<uint64_t> mats((size_t)NB * N * NW);
    {
        std::ifstream in(argv[1], std::ios::binary);
        in.read(reinterpret_cast<char*>(mats.data()), (std::streamsize)(mats.size() * 8));
        if (!in || (size_t)in.gcount() != mats.size() * 8) {
            std::fprintf(stderr, "rankmask: short read\n");
            return 3;
        }
    }
#ifdef _OPENMP
    if (argc > 4) omp_set_num_threads(std::atoi(argv[4]));
#endif

    const unsigned total = 1u << NB;
    std::vector<int> ranks(total, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<uint64_t> a((size_t)N * NW);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
        for (int mask = 1; mask < (int)total; ++mask) {
            std::fill(a.begin(), a.end(), 0ULL);
            for (int t = 0; t < NB; ++t) {
                if (!((mask >> t) & 1)) continue;
                const uint64_t* src = &mats[(size_t)t * N * NW];
                for (size_t i = 0; i < (size_t)N * NW; ++i) a[i] ^= src[i];
            }
            int rank = 0;
            for (int col = 0; col < N; ++col) {
                const int w = col >> 6, b = col & 63;
                int p = rank;
                while (p < N && ((a[(size_t)p * NW + w] >> b) & 1ULL) == 0) ++p;
                if (p == N) continue;
                if (p != rank)
                    for (int q = 0; q < NW; ++q)
                        std::swap(a[(size_t)p * NW + q], a[(size_t)rank * NW + q]);
                for (int r = rank + 1; r < N; ++r) {
                    if ((a[(size_t)r * NW + w] >> b) & 1ULL) {
                        for (int q = w; q < NW; ++q)
                            a[(size_t)r * NW + q] ^= a[(size_t)rank * NW + q];
                    }
                }
                ++rank;
            }
            ranks[mask] = rank;
        }
    }

    std::string out;
    out.reserve((size_t)total * 10);
    char buf[32];
    for (unsigned m = 1; m < total; ++m) {
        int len = std::snprintf(buf, sizeof buf, "%u %d\n", m, ranks[m]);
        out.append(buf, len);
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
