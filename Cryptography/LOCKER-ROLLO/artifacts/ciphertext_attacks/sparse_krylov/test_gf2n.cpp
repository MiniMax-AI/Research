// test_gf2n.cpp -- correctness + micro-benchmark of the GF(2^n) PCLMUL kernel.
#include "gf2n_pclmul.hpp"
#include <chrono>
#include <random>
#include <cassert>

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// reference schoolbook GF(2^n) multiply on raw bit-polynomials (big-endian-free,
// uses std::vector<uint64_t>), used only to validate the fast kernel.
struct Ref {
    int n; std::vector<uint64_t> f;
    Ref(int n_, std::vector<uint64_t> f_) : n(n_), f(f_) {}
    std::vector<uint64_t> mul(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) const {
        int W = (n + 63) / 64;
        std::vector<uint64_t> acc(2 * W + 2, 0);
        for (int i = 0; i < n; ++i) {
            if (!((a[i / 64] >> (i % 64)) & 1ULL)) continue;
            for (int j = 0; j < n; ++j) {
                if ((b[j / 64] >> (j % 64)) & 1ULL) acc[(i + j) / 64] ^= 1ULL << ((i + j) % 64);
            }
        }
        for (int q = 2 * W * 64 - 1; q >= n; --q) {
            if ((acc[q / 64] >> (q % 64)) & 1ULL) {
                acc[q / 64] ^= 1ULL << (q % 64);
                for (int i = 0; i < n; ++i)
                    if ((f[i / 64] >> (i % 64)) & 1ULL) acc[(q - n + i) / 64] ^= 1ULL << ((q - n + i) % 64);
            }
        }
        acc.resize(W);
        return acc;
    }
};

int main(int argc, char **argv) {
    // (n, modulus exponents) for every field used by the target parameter sets.
    struct Cfg { int n; std::vector<int> exps; const char *tag; };
    std::vector<Cfg> cfgs = {
        {83,  {83,7,4,2,0},        "LOCKER-I/II  n=83"},
        {89,  {89,6,5,3,0},        "LOCKER-III   n=89"},
        {101, {101,7,6,1,0},       "LOCKER-IV    n=101"},
        {103, {103,9,0},           "LOCKER-V/VI  n=103"},
        {107, {107,9,7,4,0},       "g(2^107)     n=107"},
        {109, {109,5,4,2,0},       "g(2^109)     n=109"},
        {113, {113,9,0},           "ROLLO-I-256  n=113"},
        {149, {149,10,9,7,0},      "LOCKER-VII   n=149"},
        {157, {157,6,5,2,0},       "LOCKER-IX    n=157"},
        {189, {189,6,5,2,0},       "ROLLO-II-128 n=189"},
        {193, {193,15,0},          "ROLLO-II-192 n=193"},
        {211, {211,11,10,8,0},     "ROLLO-II-256 n=211"},
    };
    std::mt19937_64 rng(12345);
    printf("%-24s %3s %8s %14s %12s %10s\n", "field", "W", "verified", "ns/mul", "ns/muladd", "cyc@3GHz");
    for (auto &c : cfgs) {
        int W = (c.n + 63) / 64;
        std::vector<uint64_t> fb(W, 0);
        for (int e : c.exps) if (e < c.n) fb[e / 64] |= 1ULL << (e % 64);
        Gf2n F; F.init(c.n, fb.data());
        Ref R(c.n, fb);

        // --- correctness: 200 random products + inverse checks
        int ok = 0;
        for (int t = 0; t < 200; ++t) {
            std::vector<uint64_t> a(W), b(W), got(W);
            for (int i = 0; i < W; ++i) { a[i] = rng(); b[i] = rng(); }
            a[W - 1] &= F.mask_top; b[W - 1] &= F.mask_top;
            F.mul(got.data(), a.data(), b.data());
            auto ref = R.mul(a, b);
            bool good = true;
            for (int i = 0; i < W; ++i) if (got[i] != ref[i]) good = false;
            if (!good) { printf("MUL MISMATCH n=%d trial=%d\n", c.n, t); return 1; }
            // muladd consistency
            std::vector<uint64_t> acc(W, 0), acc2(W, 0);
            for (int i = 0; i < W; ++i) { acc[i] = rng(); acc2[i] = acc[i]; }
            acc[W-1] &= F.mask_top; acc2[W-1] &= F.mask_top;
            F.muladd(acc2.data(), a.data(), b.data());
            for (int i = 0; i < W; ++i) acc[i] ^= got[i];
            for (int i = 0; i < W; ++i) if (acc[i] != acc2[i]) { printf("MULADD MISMATCH n=%d\n", c.n); return 1; }
            ok++;
        }
        // inverse
        for (int t = 0; t < 20; ++t) {
            std::vector<uint64_t> a(W), ai(W), one(W, 0);
            for (int i = 0; i < W; ++i) a[i] = rng();
            a[W - 1] &= F.mask_top;
            if (F.iszero(a.data())) continue;
            F.inv(ai.data(), a.data());
            std::vector<uint64_t> pr(W);
            F.mul(pr.data(), a.data(), ai.data());
            one[0] = 1;
            for (int i = 0; i < W; ++i) if (pr[i] != one[i]) { printf("INV MISMATCH n=%d\n", c.n); return 1; }
        }

        // --- benchmark
        const int NB = 4096;
        std::vector<uint64_t> A((size_t)NB * W), B((size_t)NB * W), Cv(W, 0);
        for (size_t i = 0; i < A.size(); ++i) { A[i] = rng(); B[i] = rng(); }
        for (int i = 0; i < NB; ++i) { A[(size_t)i*W+W-1] &= F.mask_top; B[(size_t)i*W+W-1] &= F.mask_top; }
        const long iters = 2000000;
        const int MASK = NB - 1;
        int ia = 0, ib = 1;
        double t0 = now();
        for (long it = 0; it < iters; ++it) {
            F.mul(Cv.data(), &A[(size_t)ia * W], &B[(size_t)ib * W]);
            ia = (ia + 1) & MASK; ib = (ib + 3) & MASK;
        }
        double t1 = now();
        double ns_mul = (t1 - t0) * 1e9 / iters;
        ia = 0; ib = 1;
        t0 = now();
        for (long it = 0; it < iters; ++it) {
            F.muladd(Cv.data(), &A[(size_t)ia * W], &B[(size_t)ib * W]);
            ia = (ia + 1) & MASK; ib = (ib + 3) & MASK;
        }
        t1 = now();
        double ns_ma = (t1 - t0) * 1e9 / iters;
        printf("%-24s %3d %8d %14.2f %12.2f %10.1f\n", c.tag, W, ok, ns_mul, ns_ma, ns_ma * 3.0);
        fflush(stdout);
    }
    printf("ALL GF(2^n) KERNEL TESTS PASSED\n");
    return 0;
}
