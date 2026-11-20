// sparse_bw.hpp -- sparse matrices over GF(2^n) + (block) SpMV + Wiedemann kernel solver.
#pragma once
#include "gf2n_pclmul.hpp"
#include <vector>
#include <random>
#include <cstdio>
#include <omp.h>
#include <chrono>

static inline double wtime() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------- sparse matrix
struct SparseMat {
    size_t nrows = 0, ncols = 0, nnz = 0;
    int W = 0;
    std::vector<size_t> rowptr;      // nrows + 1
    std::vector<uint32_t> col;       // nnz
    std::vector<uint64_t> val;       // nnz * W
    size_t bytes() const {
        return rowptr.size() * 8 + col.size() * 4 + val.size() * 8;
    }
};

// y = A x   (single RHS).  x, y are ncols / nrows vectors of W words.
inline void spmv(const Gf2n &F, const SparseMat &A, const uint64_t *x, uint64_t *y) {
    const int W = F.W;
    const size_t M = A.nrows;
#pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)M; ++i) {
        uint64_t acc[GF2N_WMAX];
        F.zero(acc);
        for (size_t p = A.rowptr[i]; p < A.rowptr[i + 1]; ++p)
            F.muladd(acc, &A.val[p * W], &x[(size_t)A.col[p] * W]);
        F.copy(&y[(size_t)i * W], acc);
    }
}

// Y = A X with beta right-hand sides.  Layout: X[j*beta + k], Y[i*beta + k].
inline void spmv_blk(const Gf2n &F, const SparseMat &A, const uint64_t *X, uint64_t *Y, int beta) {
    const int W = F.W;
    const size_t M = A.nrows;
#pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)M; ++i) {
        uint64_t acc[64 * GF2N_WMAX];
        for (int k = 0; k < beta * W; ++k) acc[k] = 0;
        for (size_t p = A.rowptr[i]; p < A.rowptr[i + 1]; ++p) {
            const uint64_t *a = &A.val[p * W];
            const uint64_t *xv = &X[(size_t)A.col[p] * beta * W];
            for (int k = 0; k < beta; ++k) F.muladd(acc + k * W, a, xv + k * W);
        }
        for (int k = 0; k < beta * W; ++k) Y[(size_t)i * beta * W + k] = acc[k];
    }
}

// y = A^T x  (single RHS) via per-thread accumulation buffers.
inline void spmvT(const Gf2n &F, const SparseMat &A, const uint64_t *x, uint64_t *y) {
    const int W = F.W;
    const size_t N = A.ncols, M = A.nrows;
    int nt = omp_get_max_threads();
    std::vector<uint64_t> buf((size_t)nt * N * W, 0);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint64_t *b = &buf[(size_t)tid * N * W];
#pragma omp for schedule(static)
        for (long i = 0; i < (long)M; ++i) {
            const uint64_t *xi = &x[(size_t)i * W];
            for (size_t p = A.rowptr[i]; p < A.rowptr[i + 1]; ++p)
                F.muladd(&b[(size_t)A.col[p] * W], &A.val[p * W], xi);
        }
    }
#pragma omp parallel for schedule(static)
    for (long j = 0; j < (long)N; ++j) {
        uint64_t acc[GF2N_WMAX];
        F.zero(acc);
        for (int t = 0; t < nt; ++t) F.add(acc, &buf[((size_t)t * N + j) * W]);
        F.copy(&y[(size_t)j * W], acc);
    }
}

// Y = A^T X, beta right-hand sides, same layout as spmv_blk.
inline void spmvT_blk(const Gf2n &F, const SparseMat &A, const uint64_t *X, uint64_t *Y, int beta) {
    const int W = F.W;
    const size_t N = A.ncols, M = A.nrows;
    int nt = omp_get_max_threads();
    std::vector<uint64_t> buf((size_t)nt * N * beta * W, 0);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint64_t *b = &buf[(size_t)tid * N * beta * W];
#pragma omp for schedule(static)
        for (long i = 0; i < (long)M; ++i) {
            const uint64_t *xi = &X[(size_t)i * beta * W];
            for (size_t p = A.rowptr[i]; p < A.rowptr[i + 1]; ++p) {
                uint64_t *dst = &b[(size_t)A.col[p] * beta * W];
                const uint64_t *a = &A.val[p * W];
                for (int k = 0; k < beta; ++k) F.muladd(dst + k * W, a, xi + k * W);
            }
        }
    }
#pragma omp parallel for schedule(static)
    for (long j = 0; j < (long)N; ++j) {
        for (int k = 0; k < beta * W; ++k) {
            uint64_t s = 0;
            for (int t = 0; t < nt; ++t) s ^= buf[((size_t)t * N + j) * beta * W + k];
            Y[(size_t)j * beta * W + k] = s;
        }
    }
}

// ------------------------------------------------------------- Berlekamp-Massey
// Returns the MINIMAL ANNIHILATING POLYNOMIAL of s[0..L-1] in the form
//   mu(t) = sum_{i=0..d} mu_i t^i ,  mu_d = 1 ,  sum_i mu_i s[k+i] = 0 for all k.
// (Massey's connection polynomial C satisfies sum_i C_i s[n-i] = 0 with C_0 = 1,
//  so mu is exactly the reversal of C.)
inline std::vector<uint64_t> berlekamp_massey(const Gf2n &F, const std::vector<uint64_t> &s, size_t L) {
    const int W = F.W;
    std::vector<uint64_t> C((L + 2) * W, 0), B((L + 2) * W, 0), T((L + 2) * W, 0);
    C[0] = 1; B[0] = 1;                       // C = B = 1, Lc = 0, m = 1
    size_t Lc = 0, m = 1;
    uint64_t b[GF2N_WMAX]; F.zero(b); b[0] = 1;
    uint64_t d[GF2N_WMAX], tmp[GF2N_WMAX], inv[GF2N_WMAX];
    for (size_t N_ = 0; N_ < L; ++N_) {
        // d = s[N_] + sum_{i=1..Lc} C[i]*s[N_-i]
        F.copy(d, &s[N_ * W]);
        for (size_t i = 1; i <= Lc && i <= N_; ++i)
            if (!F.iszero(&C[i * W])) { F.mul(tmp, &C[i * W], &s[(N_ - i) * W]); F.add(d, tmp); }
        if (F.iszero(d)) { ++m; }
        else {
            T = C;
            F.inv(inv, b);
            F.mul(tmp, d, inv);
            // C = C - (d/b) x^m B   (char 2: XOR)
            for (size_t i = 0; i + m <= L + 1; ++i)
                if (!F.iszero(&B[i * W])) { uint64_t u[GF2N_WMAX]; F.mul(u, tmp, &B[i * W]); F.add(&C[(i + m) * W], u); }
            if (2 * Lc <= N_) { Lc = N_ + 1 - Lc; B = T; F.copy(b, d); m = 1; }
            else ++m;
        }
    }
    std::vector<uint64_t> out((Lc + 1) * W, 0);
    for (size_t i = 0; i <= Lc; ++i) F.copy(&out[i * W], &C[(Lc - i) * W]);   // reverse
    return out;
}

// ------------------------------------------------------------------- utilities
inline void rand_vec(const Gf2n &F, uint64_t *v, size_t n, std::mt19937_64 &rng) {
    for (size_t i = 0; i < n * F.W; ++i) v[i] = rng();
    for (size_t i = 0; i < n; ++i) v[i * F.W + F.W - 1] &= F.mask_top;
}
inline void zero_vec(const Gf2n &F, uint64_t *v, size_t n) {
    for (size_t i = 0; i < n * F.W; ++i) v[i] = 0;
}
inline void dot(const Gf2n &F, uint64_t *out, const uint64_t *a, const uint64_t *b, size_t n) {
    uint64_t acc[GF2N_WMAX]; F.zero(acc);
    uint64_t t[GF2N_WMAX];
    for (size_t i = 0; i < n; ++i) { F.mul(t, &a[i * F.W], &b[i * F.W]); F.add(acc, t); }
    F.copy(out, acc);
}
