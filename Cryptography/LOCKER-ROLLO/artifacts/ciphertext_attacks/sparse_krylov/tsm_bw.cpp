// tsm_bw.cpp -- build the Degree-1 TSM ("MaxMinors") sparse system over GF(2^n) for a
// real LOCKER/ROLLO instance, then solve its kernel with (block) Wiedemann / Lanczos
// style Krylov methods and benchmark the sparse kernels.
//
//   ./tsm_bw build <inst.bin> <out.mat> [--stats]
//   ./tsm_bw solve <mat> [maxfrac]
//   ./tsm_bw bench <mat> <beta> [reps]
//
// Matrix file: 'TSMM', u64 nrows,ncols,nnz,W, u64 rowptr[nrows+1], u32 col[nnz], u64 val[nnz*W]
#include "sparse_bw.hpp"
#include <algorithm>
#include <cstring>
#include <string>

// ------------------------------------------------------------------ Q = L[T]/(g)
struct QRing {
    const Gf2n *F; int m; std::vector<uint64_t> g;   // g: m bits (0/1 coefficients)
    void init(const Gf2n *F_, int m_, uint64_t gbits) {
        F = F_; m = m_; g.assign(m, 0);
        for (int j = 0; j < m; ++j) g[j] = (gbits >> j) & 1ULL;
    }
    // c = a*b in L[T]/(g); a,b,c are m elements of L (W words each)
    void mul(uint64_t *c, const uint64_t *a, const uint64_t *b, uint64_t *scratch) const {
        const int W = F->W, m_ = m;
        for (int i = 0; i < (2 * m_ - 1) * W; ++i) scratch[i] = 0;
        for (int i = 0; i < m_; ++i) {
            if (F->iszero(&a[i * W])) continue;
            for (int j = 0; j < m_; ++j) {
                if (F->iszero(&b[j * W])) continue;
                F->muladd(&scratch[(i + j) * W], &a[i * W], &b[j * W]);
            }
        }
        for (int k = 2 * m_ - 2; k >= m_; --k) {
            if (F->iszero(&scratch[k * W])) continue;
            for (int j = 0; j < m_; ++j)
                if (g[j]) F->add(&scratch[(k - m_ + j) * W], &scratch[k * W]);
            F->zero(&scratch[k * W]);
        }
        for (int i = 0; i < m_ * W; ++i) c[i] = scratch[i];
    }
};

// ---------------------------------------------------------------- combinations
static void combos(int start, int need, int total, std::vector<int> &cur,
                   std::vector<std::vector<int>> &out) {
    if ((int)cur.size() == need) { out.push_back(cur); return; }
    for (int v = start; v < total; ++v) {
        cur.push_back(v);
        combos(v + 1, need, total, cur, out);
        cur.pop_back();
    }
}

struct Inst {
    int m, n, r, P, W;
    std::vector<uint64_t> H, U;     // P * m * W each
};

static Inst read_inst(const char *path, uint64_t &gbits_out) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    char magic[4]; fread(magic, 1, 4, f);
    if (memcmp(magic, "TSMI", 4)) { fprintf(stderr, "bad magic\n"); exit(1); }
    int hdr[5]; fread(hdr, 4, 5, f);
    Inst I; I.m = hdr[0]; I.n = hdr[1]; I.r = hdr[2]; I.P = hdr[3]; I.W = hdr[4];
    size_t sz = (size_t)I.P * I.m * I.W;
    I.H.resize(sz); I.U.resize(sz);
    fread(I.H.data(), 8, sz, f); fread(I.U.data(), 8, sz, f);
    fclose(f);
    gbits_out = 0;   // filled by caller from params
    return I;
}

static void write_mat(const char *path, const SparseMat &A) {
    FILE *f = fopen(path, "wb");
    uint64_t h[4] = {A.nrows, A.ncols, A.nnz, (uint64_t)A.W};
    fwrite("TSMM", 1, 4, f); fwrite(h, 8, 4, f);
    fwrite(A.rowptr.data(), 8, A.nrows + 1, f);
    fwrite(A.col.data(), 4, A.nnz, f);
    fwrite(A.val.data(), 8, A.nnz * A.W, f);
    fclose(f);
}
static SparseMat read_mat(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    char magic[4]; fread(magic, 1, 4, f);
    if (memcmp(magic, "TSMM", 4)) { fprintf(stderr, "bad matrix magic\n"); exit(1); }
    uint64_t h[4]; fread(h, 8, 4, f);
    SparseMat A; A.nrows = h[0]; A.ncols = h[1]; A.nnz = h[2]; A.W = (int)h[3];
    A.rowptr.resize(A.nrows + 1); A.col.resize(A.nnz); A.val.resize(A.nnz * A.W);
    fread(A.rowptr.data(), 8, A.nrows + 1, f);
    fread(A.col.data(), 4, A.nnz, f);
    fread(A.val.data(), 8, A.nnz * A.W, f);
    fclose(f);
    return A;
}

// ------------------------------------------------------------------- build
static SparseMat build_maxminors(const Inst &I, const Gf2n &F, uint64_t gbits, bool stats) {
    const int m = I.m, r = I.r, P = I.P, W = F.W;
    QRing Q; Q.init(&F, m, gbits);

    // products[mask] = prod_{q in mask} H_q
    const size_t NM = (size_t)1 << P;
    std::vector<uint64_t> prod(NM * (size_t)m * W, 0);
    prod[0] = 1;                                     // multiplicative unit
    std::vector<uint64_t> scratch((size_t)(2 * m) * W, 0);
    for (size_t mask = 1; mask < NM; ++mask) {
        int b = __builtin_ctzll(mask);
        Q.mul(&prod[mask * m * W], &prod[(mask ^ ((size_t)1 << b)) * m * W], &I.H[(size_t)b * m * W], scratch.data());
    }
    if (stats) fprintf(stderr, "products built (2^%d = %zu Q-elements)\n", P, NM);

    // r-subsets of {0..P-1}  -> column blocks
    std::vector<std::vector<int>> subs; { std::vector<int> cur; combos(0, r, P, cur, subs); }
    const size_t nblk = subs.size();
    // (r+1)-subsets containing 0  ==  r-subsets of {1..P-1}  -> equations
    std::vector<std::vector<int>> eqs; { std::vector<int> cur; combos(1, r, P, cur, eqs); }
    const size_t neq = eqs.size();

    // hash map from subset code -> block index
    std::vector<uint64_t> code(nblk);
    for (size_t i = 0; i < nblk; ++i) {
        uint64_t c = 0; for (int v : subs[i]) c |= 1ULL << v; code[i] = c;
    }
    std::sort(code.begin(), code.end());
    auto block_of = [&](uint64_t c) -> size_t {
        return std::lower_bound(code.begin(), code.end(), c) - code.begin();
    };

    const size_t SIDES = (size_t)1 << r;
    const size_t ncols = nblk * SIDES;
    const size_t nrows = neq * m;
    const size_t w = (size_t)(r + 1) * SIDES;         // row weight

    SparseMat A; A.nrows = nrows; A.ncols = ncols; A.W = W; A.nnz = nrows * w;
    A.rowptr.resize(nrows + 1);
    for (size_t i = 0; i <= nrows; ++i) A.rowptr[i] = i * w;
    A.col.resize(A.nnz); A.val.assign(A.nnz * W, 0);

    double t0 = wtime();
#pragma omp parallel
    {
        std::vector<uint64_t> scr((size_t)(2 * m) * W, 0), coef((size_t)m * W, 0);
#pragma omp for schedule(dynamic, 1)
        for (long e = 0; e < (long)neq; ++e) {
            int S[16]; S[0] = 0;
            for (int k = 0; k < r; ++k) S[k + 1] = eqs[e][k];
            for (int omit = 0; omit <= r; ++omit) {
                int T[16]; int ti = 0;
                for (int k = 0; k <= r; ++k) if (k != omit) T[ti++] = S[k];
                uint64_t tc = 0; for (int k = 0; k < r; ++k) tc |= 1ULL << T[k];
                size_t blk = block_of(tc);
                const uint64_t *Uq = &I.U[(size_t)S[omit] * m * W];
                for (size_t side = 0; side < SIDES; ++side) {
                    size_t hmask = 0;
                    for (int k = 0; k < r; ++k) if ((side >> k) & 1) hmask |= (size_t)1 << T[k];
                    Q.mul(coef.data(), Uq, &prod[hmask * m * W], scr.data());
                    size_t col = blk * SIDES + side;
                    for (int dd = 0; dd < m; ++dd) {
                        size_t row = (size_t)e * m + dd;
                        size_t pos = A.rowptr[row] + (size_t)omit * SIDES + side;
                        A.col[pos] = (uint32_t)col;
                        F.copy(&A.val[pos * W], &coef[(size_t)dd * W]);
                    }
                }
            }
        }
    }
    if (stats) fprintf(stderr, "matrix %zu x %zu, nnz=%zu (row weight %zu) built in %.1fs\n",
                       nrows, ncols, A.nnz, w, wtime() - t0);
    return A;
}

// ------------------------------------------------------------------- solve
struct SolveStats { size_t spmv = 0; double t = 0; size_t deg = 0; size_t e = 0; bool ok = false; };

static SolveStats wiedemann_kernel(const Gf2n &F, const SparseMat &A, uint64_t seed,
                                   double maxfrac, std::vector<uint64_t> &kernel_out) {
    const int W = F.W;
    const size_t N = A.ncols, M = A.nrows;
    std::mt19937_64 rng(seed);
    SolveStats st;
    size_t Lmax = (size_t)(2.0 * N * maxfrac) + 8;
    std::vector<uint64_t> u(N * W), x0(N * W), tmp(M * W), tmp2(N * W);
    rand_vec(F, u.data(), N, rng);      // u0 = y
    rand_vec(F, x0.data(), N, rng);     // x0 = x  (FIXED probe vector: s_i = x^T B^i y)
    std::vector<uint64_t> s(Lmax * W, 0);
    double t0 = wtime();
    size_t L = 0;
    for (size_t i = 0; i < Lmax; ++i) {
        dot(F, &s[i * W], x0.data(), u.data(), N);      // s_i = x^T B^i y
        // u <- B u = A^T (A u)
        spmv(F, A, u.data(), tmp.data());   ++st.spmv;
        spmvT(F, A, tmp.data(), tmp2.data()); ++st.spmv;
        std::swap(u, tmp2);
        L = i + 1;
        if ((i & 255) == 255) {
            // early exit: try BM on what we have
            auto mu = berlekamp_massey(F, s, L);
            size_t d = mu.size() / W - 1;
            if (2 * d <= L - 1) break;
        }
    }
    st.t = wtime() - t0;
    auto mu = berlekamp_massey(F, s, L);
    size_t d = mu.size() / W - 1;
    st.deg = d;
    // strip t^e
    size_t e = 0;
    while (e <= d && F.iszero(&mu[e * W])) ++e;
    if (e == 0 || e > d) { st.ok = false; return st; }
    // w = nu(B) u0   via Horner (pass 2)
    std::vector<uint64_t> w(N * W, 0), acc(N * W, 0), t1v(M * W), t2v(N * W);
    // regenerate u0 deterministically
    { std::mt19937_64 r2(seed); rand_vec(F, u.data(), N, r2); }
    // w <- mu_d * u0   (Horner seed for nu(B) u0, nu = mu / t^e)
#pragma omp parallel for schedule(static)
    for (long j = 0; j < (long)N; ++j)
        F.mul(&w[(size_t)j * W], &mu[d * W], &u[(size_t)j * W]);
    for (size_t i = d; i-- > e; ) {
        // w <- B w
        spmv(F, A, w.data(), t1v.data());   ++st.spmv;
        spmvT(F, A, t1v.data(), t2v.data()); ++st.spmv;
        std::swap(w, t2v);
        if (!F.iszero(&mu[i * W])) {
            // w += mu_i * u0
#pragma omp parallel for schedule(static)
            for (long j = 0; j < (long)N; ++j) {
                uint64_t t[GF2N_WMAX];
                F.mul(t, &mu[i * W], &u[(size_t)j * W]);
                F.add(&w[(size_t)j * W], t);
            }
        }
    }
    // v = B^{e-1} w
    for (size_t i = 0; i + 1 < e; ++i) {
        spmv(F, A, w.data(), t1v.data());   ++st.spmv;
        spmvT(F, A, t1v.data(), t2v.data()); ++st.spmv;
        std::swap(w, t2v);
    }
    st.t = wtime() - t0;
    // verify A v == 0
    std::vector<uint64_t> chk(M * W, 0);
    spmv(F, A, w.data(), chk.data()); ++st.spmv;
    size_t nzc = 0;
    for (size_t i = 0; i < M; ++i) if (!F.iszero(&chk[i * W])) ++nzc;
    fprintf(stderr, "[dbg] L=%zu d=%zu e=%zu  A*v nonzero rows=%zu/%zu\n", L, d, e, nzc, M);
    st.ok = (nzc == 0);
    size_t nzw = 0;
    for (size_t i = 0; i < N; ++i) if (!F.iszero(&w[i * W])) ++nzw;
    fprintf(stderr, "[dbg] v nonzero coords=%zu/%zu\n", nzw, N);
    if (nzw == 0) st.ok = false;
    st.e = e;
    kernel_out = w;
    return st;
}

// -------------------------------------------------------------------- main
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: tsm_bw build|solve|bench ...\n"); return 1; }
    std::string cmd = argv[1];
    if (cmd == "build") {
        bool stats = true;
        uint64_t dummy;
        Inst I = read_inst(argv[2], dummy);
        uint64_t gbits = 0;
        { std::string g = argv[4]; size_t pos = 0;
          while (pos < g.size()) { size_t c = g.find(',', pos); if (c == std::string::npos) c = g.size();
            int e = atoi(g.substr(pos, c - pos).c_str()); if (e < I.m) gbits |= 1ULL << e; pos = c + 1; } }
        std::vector<uint64_t> fb(I.W, 0);
        std::string ex = argv[5];
        size_t pos = 0;
        while (pos < ex.size()) {
            size_t c = ex.find(',', pos);
            if (c == std::string::npos) c = ex.size();
            int e = atoi(ex.substr(pos, c - pos).c_str());
            if (e < I.n) fb[e / 64] |= 1ULL << (e % 64);
            pos = c + 1;
        }
        Gf2n F; F.init(I.n, fb.data());
        SparseMat A = build_maxminors(I, F, gbits, stats);
        write_mat(argv[3], A);
        // exact row-weight audit
        size_t zeros = 0;
        for (size_t p = 0; p < A.nnz; ++p) if (F.iszero(&A.val[p * A.W])) ++zeros;
        printf("{\"nrows\":%zu,\"ncols\":%zu,\"nnz\":%zu,\"rowweight\":%zu,\"zero_entries\":%zu,\"bytes\":%zu}\n",
               A.nrows, A.ncols, A.nnz, A.nnz / A.nrows, zeros, A.bytes());
        return 0;
    }
    // solve / bench need the field parameters: argv pattern  <mat> <nexps> ...
    if (cmd == "solve") {
        SparseMat A = read_mat(argv[2]);
        std::vector<uint64_t> fb(A.W, 0);
        std::string ex = argv[3]; size_t pos = 0;
        while (pos < ex.size()) {
            size_t c = ex.find(',', pos); if (c == std::string::npos) c = ex.size();
            int e = atoi(ex.substr(pos, c - pos).c_str());
            if (e < A.W * 64) fb[e / 64] |= 1ULL << (e % 64);
            pos = c + 1;
        }
        Gf2n F; F.init(atoi(argv[4]), fb.data());
        double maxfrac = (argc > 5) ? atof(argv[5]) : 1.0;
        std::vector<uint64_t> kv;
        SolveStats st = wiedemann_kernel(F, A, 12345, maxfrac, kv);
        printf("{\"nrows\":%zu,\"ncols\":%zu,\"nnz\":%zu,\"mindeg\":%zu,\"t_ord\":%zu,\"spmvs\":%zu,\"seconds\":%.2f,\"kernel_verified\":%s}\n",
               A.nrows, A.ncols, A.nnz, st.deg, st.e, st.spmv, st.t, st.ok ? "true" : "false");
        return st.ok ? 0 : 2;
    }
    if (cmd == "bench") {
        SparseMat A = read_mat(argv[2]);
        std::vector<uint64_t> fb(A.W, 0);
        std::string ex = argv[3]; size_t pos = 0;
        while (pos < ex.size()) {
            size_t c = ex.find(',', pos); if (c == std::string::npos) c = ex.size();
            int e = atoi(ex.substr(pos, c - pos).c_str());
            if (e < A.W * 64) fb[e / 64] |= 1ULL << (e % 64);
            pos = c + 1;
        }
        Gf2n F; F.init(atoi(argv[4]), fb.data());
        int beta = atoi(argv[5]);
        int reps = (argc > 6) ? atoi(argv[6]) : 5;
        std::mt19937_64 rng(999);
        std::vector<uint64_t> X(A.ncols * beta * A.W), Y(A.nrows * beta * A.W), Z(A.ncols * beta * A.W);
        rand_vec(F, X.data(), A.ncols * beta, rng);
        // warm-up
        spmv_blk(F, A, X.data(), Y.data(), beta);
        double t0 = wtime();
        for (int i = 0; i < reps; ++i) spmv_blk(F, A, X.data(), Y.data(), beta);
        double t1 = wtime();
        for (int i = 0; i < reps; ++i) spmvT_blk(F, A, Y.data(), Z.data(), beta);
        double t2 = wtime();
        double per = (t1 - t0) / reps, perT = (t2 - t1) / reps;
        double mad = (double)A.nnz * beta;
        printf("{\"beta\":%d,\"spmv_s\":%.4f,\"spmvT_s\":%.4f,\"ns_per_muladd\":%.3f,\"ns_per_muladd_T\":%.3f,\"muladds_per_s\":%.3e}\n",
               beta, per, perT, per * 1e9 / mad, perT * 1e9 / mad, mad / per);
        return 0;
    }
    fprintf(stderr, "unknown command\n");
    return 1;
}
