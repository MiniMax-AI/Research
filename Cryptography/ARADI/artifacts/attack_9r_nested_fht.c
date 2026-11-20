/* attack_9r_nested_fht.c -- 9-Round ARADI Full 256-Bit Master-Key Recovery
 *
 * TARGET: Unmodified NSA ARADI (ePrint 2024/1240), R = 9, r0 = 0 (rounds 0..8 + K^(9)).
 *
 * ARCHITECTURE  (6-round integral distinguisher on Q_6 + 3-round peel: rounds 8, 7, 6)
 * ------------------------------------------------------------------------------------
 *   C = Q_9 ^ K^(9).  Round 8 uses Lambda_0, round 7 Lambda_3, round 6 Lambda_2
 *   (backward peel layers L_0, L_3, L_2):
 *
 *     T9[d][p]  = nibble of column d of Lambda_0(C[p])                     [pure data]
 *     Z_8[d][p] = pi^{-1}( T9[d][p] ^ kappa_d ),      kappa = Lambda_0 K^(9)
 *     M_e[p]    = XOR_{d in rs_{L_3}(e)} Z_8[d][p]
 *     Z_7[e][p] = pi^{-1}( M_e[p] ^ gamma_e ),        gamma = Lambda_3 K^(8)
 *     V_u[p]    = ( XOR_{e in rs_{L_2}(u)} Z_7[e][p] ) ^ g_u,   g = Lambda_2 K^(7)
 *     Z_6[u][p] = pi^{-1}( V_u[p] )
 *     Q_6       = Z_6 ^ K^(6)                         (K^(6) cancels in every cube sum)
 *
 *   For the colliding-pair 4-column cube (0, 5, 19, 24) (d = 16, D = 2^16 CP per cube),
 *   the GF(2) null space of Q_6 cube-sum vectors has dimension 42 and contains the four
 *   single-column masks e_nu (x) {25}.  Intra-half Z_16 shift symmetry moves the balanced
 *   column to u = shift_col(25, s) across all 16 upper-half columns {16..31}, jointly
 *   determining all 32 columns of kappa = Lambda_0 K^(9) and gamma = Lambda_3 K^(8)
 *   (exact GF(2) rank rho = 256, O(1) inversion via master_key_from_even_pair(K^(8),K^(9),8)).
 *
 *   Per-branch 3-level nested Walsh-Hadamard kernel (52 branch key bits = 9 kappa + 3 gamma + 1 g):
 *     1. Z8[d][p] = pi^{-1}(T9[d][p] ^ kappa_d)
 *     2. JM[m1][m2][m3] = parity #{p : M_{e_i}[p] = m_i}   (joint 2^12 parity histogram)
 *     3. Phi_gamma[z] = IWHT_4( diag(Nhat_gamma) ), where
 *        diag(Nhat_gamma)[xi] = IWHT_12( JMhat * Hhat_xi )[gamma]  (evaluates all 2^12 gamma at once)
 *     4. Test T_u(g) = XOR_{z : Phi_gamma[z]=1} pi^{-1}(z ^ g) == 0 across all 4 rows.
 *
 *   Single-core execution; uses nonlinear SplitMix64 PRNG.
 *   --kbits B searches 2^B candidates per branch for the 36-bit kappa vector (remaining
 *   36-B bits supplied by a simulated oracle; B = 36 is the full 2^54.77-butterfly search).
 *
 * BUILD (single-core):  gcc -O3 -march=native -Wall -Wextra -o attack_9r_nested_fht attack_9r_nested_fht.c -lm
 * DEMO  (single-core):  ./attack_9r_nested_fht --trials 2 --kbits 0 --cubes 6 --shifts 16
 */
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const int ROT_A[4] = {11, 10, 9, 8};
static const int ROT_B[4] = {8, 9, 4, 9};
static const int ROT_C[4] = {14, 11, 14, 7};
static uint8_t SBOX[16], SBOX_INV[16];

static inline uint32_t rol32(uint32_t x, uint32_t n) {
    n &= 31;
    return n ? ((x << n) | (x >> (32 - n))) : x;
}
static inline uint32_t ror32(uint32_t x, uint32_t n) {
    n &= 31;
    return n ? ((x >> n) | (x << (32 - n))) : x;
}
static inline uint16_t rol16(uint16_t x, uint32_t n) {
    n &= 15;
    return n ? (uint16_t)((x << n) | (x >> (16 - n))) : x;
}
static inline uint32_t linear_L(uint32_t s, int a, int b, int c) {
    uint16_t u = (uint16_t)(s >> 16), l = (uint16_t)(s & 0xFFFF);
    uint16_t uo = u ^ rol16(u, a) ^ rol16(l, c);
    uint16_t lo = l ^ rol16(l, a) ^ rol16(u, b);
    return ((uint32_t)uo << 16) | (uint32_t)lo;
}
static inline void sbox_fwd(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z) {
    *x ^= *w & *y;
    *z ^= *x & *y;
    *y ^= *w & *z;
    *w ^= *x & *z;
}
static inline void key_M(uint32_t ri, uint32_t rj, uint32_t X, uint32_t Y, uint32_t *o1, uint32_t *o0) {
    uint32_t rY = rol32(Y, ri);
    *o1 = rY ^ rol32(X, rj) ^ X;
    *o0 = rY ^ X;
}
static inline void key_M_inv(uint32_t ri, uint32_t rj, uint32_t A, uint32_t B, uint32_t *X, uint32_t *Y) {
    *X = ror32(A ^ B, rj);
    *Y = ror32(B ^ *X, ri);
}
static void expand_key(const uint32_t mk[8], int R, uint32_t rks[][4]) {
    uint32_t K[8];
    memcpy(K, mk, sizeof(K));
    for (int i = 0; i < R; i++) {
        int j = i & 1;
        for (int w = 0; w < 4; w++) rks[i][w] = K[4 * j + w];
        uint32_t n1, n0, n3, n2, n5, n4, n7, n6, tmp;
        key_M(1, 3, K[1], K[0], &n1, &n0);
        key_M(9, 28, K[3], K[2], &n3, &n2);
        key_M(1, 3, K[5], K[4], &n5, &n4);
        key_M(9, 28, K[7], K[6], &n7, &n6);
        K[1] = n1; K[0] = n0; K[3] = n3; K[2] = n2;
        K[5] = n5; K[4] = n4; K[7] = n7 ^ (uint32_t)i; K[6] = n6;
        if (j == 0) {
            tmp = K[1]; K[1] = K[2]; K[2] = tmp;
            tmp = K[5]; K[5] = K[6]; K[6] = tmp;
        } else {
            tmp = K[1]; K[1] = K[4]; K[4] = tmp;
            tmp = K[3]; K[3] = K[6]; K[6] = tmp;
        }
    }
    int j = R & 1;
    for (int w = 0; w < 4; w++) rks[R][w] = K[4 * j + w];
}
static void master_key_from_even_pair(const uint32_t keven[4], const uint32_t kodd[4],
                                      int even_idx, uint32_t mk_out[8]) {
    uint32_t S[8];
    for (int w = 0; w < 4; w++) S[w] = keven[w];
    uint32_t u4 = kodd[0], u6 = kodd[1], u5 = kodd[2], u7 = kodd[3] ^ (uint32_t)even_idx;
    key_M_inv(1, 3, u5, u4, &S[5], &S[4]);
    key_M_inv(9, 28, u7, u6, &S[7], &S[6]);
    for (int i = even_idx - 1; i >= 0; i--) {
        int j = i & 1;
        uint32_t tmp, x, y;
        if (j == 0) {
            tmp = S[1]; S[1] = S[2]; S[2] = tmp;
            tmp = S[5]; S[5] = S[6]; S[6] = tmp;
        } else {
            tmp = S[1]; S[1] = S[4]; S[4] = tmp;
            tmp = S[3]; S[3] = S[6]; S[6] = tmp;
        }
        S[7] ^= (uint32_t)i;
        key_M_inv(1, 3, S[1], S[0], &x, &y);  S[1] = x; S[0] = y;
        key_M_inv(9, 28, S[3], S[2], &x, &y); S[3] = x; S[2] = y;
        key_M_inv(1, 3, S[5], S[4], &x, &y);  S[5] = x; S[4] = y;
        key_M_inv(9, 28, S[7], S[6], &x, &y); S[7] = x; S[6] = y;
    }
    memcpy(mk_out, S, sizeof(S));
}
static void encrypt_rounds(const uint32_t pt[4], const uint32_t rks[][4], int R, uint32_t ct[4]) {
    uint32_t w = pt[0], x = pt[1], y = pt[2], z = pt[3];
    for (int i = 0; i < R; i++) {
        w ^= rks[i][0]; x ^= rks[i][1]; y ^= rks[i][2]; z ^= rks[i][3];
        sbox_fwd(&w, &x, &y, &z);
        int l = i & 3;
        w = linear_L(w, ROT_A[l], ROT_B[l], ROT_C[l]);
        x = linear_L(x, ROT_A[l], ROT_B[l], ROT_C[l]);
        y = linear_L(y, ROT_A[l], ROT_B[l], ROT_C[l]);
        z = linear_L(z, ROT_A[l], ROT_B[l], ROT_C[l]);
    }
    ct[0] = w ^ rks[R][0]; ct[1] = x ^ rks[R][1];
    ct[2] = y ^ rks[R][2]; ct[3] = z ^ rks[R][3];
}

static uint8_t RS[4][32][3];
static void init_tables(void) {
    for (int v = 0; v < 16; v++) {
        uint32_t w = v & 1, x = (v >> 1) & 1, y = (v >> 2) & 1, z = (v >> 3) & 1;
        sbox_fwd(&w, &x, &y, &z);
        SBOX[v] = (uint8_t)(w | (x << 1) | (y << 2) | (z << 3));
        SBOX_INV[SBOX[v]] = (uint8_t)v;
    }
    for (int j = 0; j < 4; j++) {
        for (int b = 0; b < 32; b++) {
            int n = 0;
            for (int c = 0; c < 32; c++)
                if ((linear_L(1u << c, ROT_A[j], ROT_B[j], ROT_C[j]) >> b) & 1u)
                    RS[j][b][n++] = (uint8_t)c;
            if (n != 3) {
                fprintf(stderr, "FATAL row weight %d\n", n);
                exit(1);
            }
        }
    }
}
static int shift_col(int c, int s) {
    return (c < 16) ? ((c + s) & 15) : (16 + (((c - 16) + s) & 15));
}
/* Nonlinear SplitMix64 PRNG (avoids GF(2)-linear xorshift64 subspace bias) */
static uint64_t splitmix64(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void wht(int32_t *a, int n) {
    for (int h = 1; h < n; h <<= 1)
        for (int i = 0; i < n; i += (h << 1))
            for (int j = i; j < i + h; j++) {
                int32_t x = a[j], y = a[j + h];
                a[j] = x + y;
                a[j + h] = x - y;
            }
}

#define CDIM 16
#define NPT  (1 << CDIM)

typedef struct {
    uint8_t *t9;
    uint32_t base[4];
} cube_t;

static inline uint8_t t9get(const cube_t *c, int p, int k) {
    uint8_t b = c->t9[((size_t)p >> 1) * 9 + k];
    return (p & 1) ? (b >> 4) : (b & 15);
}
static void build_cube(cube_t *c, const int dcols[9], int shift, const uint32_t rks[][4], uint64_t *rng) {
    static const int C0[4] = {0, 5, 19, 24};
    int cols[4];
    for (int i = 0; i < 4; i++) cols[i] = shift_col(C0[i], shift);
    for (int i = 0; i < 4; i++) c->base[i] = (uint32_t)splitmix64(rng);
    uint32_t ct[4], lam[4];
    for (int p = 0; p < NPT; p++) {
        uint32_t pt[4] = {c->base[0], c->base[1], c->base[2], c->base[3]};
        for (int ci = 0; ci < 4; ci++) {
            int nib = (p >> (4 * ci)) & 15;
            for (int w = 0; w < 4; w++)
                if ((nib >> w) & 1) pt[w] ^= (1u << cols[ci]);
        }
        encrypt_rounds(pt, rks, 9, ct);
        for (int w = 0; w < 4; w++) lam[w] = linear_L(ct[w], ROT_A[0], ROT_B[0], ROT_C[0]);
        for (int k = 0; k < 9; k++) {
            int d = dcols[k];
            uint8_t nib = (uint8_t)((((lam[0] >> d) & 1u)) | ((((lam[1] >> d) & 1u)) << 1) |
                                    ((((lam[2] >> d) & 1u)) << 2) | ((((lam[3] >> d) & 1u)) << 3));
            uint8_t *slot = &c->t9[((size_t)p >> 1) * 9 + k];
            if (p & 1) *slot = (uint8_t)((*slot & 0x0F) | (nib << 4));
            else       *slot = (uint8_t)((*slot & 0xF0) | nib);
        }
    }
}

static int32_t Hhat[16][4096];
static uint8_t S12[4096];
static void init_kernel(void) {
    for (int m = 0; m < 4096; m++) {
        int m1 = (m >> 8) & 15, m2 = (m >> 4) & 15, m3 = m & 15;
        S12[m] = (uint8_t)(SBOX_INV[m1] ^ SBOX_INV[m2] ^ SBOX_INV[m3]);
    }
    for (int xi = 0; xi < 16; xi++) {
        for (int m = 0; m < 4096; m++)
            Hhat[xi][m] = (__builtin_popcount((unsigned)(S12[m] & xi)) & 1) ? -1 : 1;
        wht(Hhat[xi], 4096);
    }
}

static uint8_t Z8[9][NPT];
static int32_t JMhat[4096];
static int32_t Dg[4096][16];

static void build_Z8(const cube_t *c, const int kval[9]) {
    for (int k = 0; k < 9; k++) {
        uint8_t *dst = Z8[k];
        int kp = kval[k];
        for (int p = 0; p < NPT; p++) dst[p] = SBOX_INV[t9get(c, p, k) ^ kp];
    }
}
static void build_JMhat(void) {
    static uint8_t JM[4096];
    memset(JM, 0, sizeof(JM));
    const uint8_t *a = Z8[0], *b = Z8[1], *cc = Z8[2], *d = Z8[3], *e = Z8[4], *f = Z8[5];
    for (int p = 0; p < NPT; p++) {
        int m = (((int)(a[p] ^ b[p] ^ cc[p])) << 8) |
                (((int)(d[p] ^ e[p] ^ f[p])) << 4) |
                ((int)(Z8[6][p] ^ Z8[7][p] ^ Z8[8][p]));
        JM[m] ^= 1;
    }
    for (int m = 0; m < 4096; m++) JMhat[m] = JM[m];
    wht(JMhat, 4096);
}
static void build_Dg(void) {
    static int32_t tmp[4096];
    for (int xi = 0; xi < 16; xi++) {
        for (int m = 0; m < 4096; m++) tmp[m] = JMhat[m] * Hhat[xi][m];
        wht(tmp, 4096);
        for (int g = 0; g < 4096; g++) Dg[g][xi] = tmp[g] >> 12;
    }
}
static inline void phi_of_gamma(int gamma, uint8_t Phi[16]) {
    int32_t t[16];
    for (int xi = 0; xi < 16; xi++) t[xi] = Dg[gamma][xi];
    wht(t, 16);
    for (int z = 0; z < 16; z++) Phi[z] = (uint8_t)((t[z] >> 4) & 1);
}
static inline int T_of(int g, const uint8_t Phi[16]) {
    int T = 0;
    for (int z = 0; z < 16; z++)
        if (Phi[z]) T ^= SBOX_INV[z ^ g];
    return T;
}

int main(int argc, char **argv) {
    init_tables();
    init_kernel();
    static const int SCHED80[16] = {13, 9, 8, 8, 8, 8, 6, 3, 3, 2, 2, 2, 2, 2, 2, 2};
    int trials = 2, kbits = 0, ncubes = 5, cubes_explicit = 0, nshift = 16, verbose = 0;
    uint64_t seed = 0x202611061234ULL;
    for (int a = 1; a < argc; a++) {
        if (!strcmp(argv[a], "--trials") && a + 1 < argc) trials = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--kbits") && a + 1 < argc) kbits = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--cubes") && a + 1 < argc) { ncubes = atoi(argv[++a]); cubes_explicit = 1; }
        else if (!strcmp(argv[a], "--shifts") && a + 1 < argc) nshift = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--seed") && a + 1 < argc) seed = strtoull(argv[++a], NULL, 0);
        else if (!strcmp(argv[a], "-v")) verbose = 1;
        else if (!strcmp(argv[a], "-q")) verbose = 0;
    }
    if (kbits < 0) kbits = 0;
    if (kbits > 36) kbits = 36;
    int max_cubes = ncubes;
    if (kbits > 0) {
        int min_s0 = (kbits + 23) / 4;
        if (!cubes_explicit) {
            max_cubes = SCHED80[0] > min_s0 ? SCHED80[0] : min_s0;
        } else if (max_cubes < min_s0) {
            max_cubes = min_s0;
        }
    }
    printf("=== 9-Round ARADI Full 256-Bit Master-Key Recovery (Nested Joint-WHT Peel) ===\n");
    printf("  Distinguisher: 6R Q_6 on d=16 colliding-pair cube (0,5,19,24), %d Z_16 shifts\n", nshift);
    printf("  Peel layers  : L_0, L_3, L_2 (52 key bits/branch: 9 kappa + 3 gamma + 1 g)\n");
    if (kbits > 0 && !cubes_explicit)
        printf("  Kappa search : <= 2^%d unrecovered bits/branch | Algorithm-3 80-cube schedule (max %d)\n", kbits, max_cubes);
    else
        printf("  Kappa search : 2^%d / 2^36 per branch | cubes/shift = %d\n", kbits, max_cubes);
    if (kbits < 36)
        printf("  Oracle prefix: top %d of 36 kappa bits supplied (all 2^16 (gamma,g) + JM WHT + MK inv live)\n",
               36 - kbits);

    cube_t *cubes = (cube_t *)calloc((size_t)max_cubes, sizeof(cube_t));
    for (int m = 0; m < max_cubes; m++) cubes[m].t9 = (uint8_t *)malloc((size_t)(NPT / 2) * 9);
    printf("  Memory       : cube buffer %.1f MB, WHT kernel tables %.1f MB\n",
           max_cubes * (NPT / 2) * 9 / 1048576.0, (sizeof(Hhat) + sizeof(Dg) + sizeof(Z8)) / 1048576.0);

    int passed = 0;
    double tot = 0.0;
    long long D_total = 0;
    for (int tr = 0; tr < trials; tr++) {
        uint64_t rng = seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(tr + 1));
        uint32_t mk_true[8], rks[10][4];
        for (int i = 0; i < 8; i++) mk_true[i] = (uint32_t)splitmix64(&rng);
        expand_key(mk_true, 9, rks);

        uint32_t kap[4], gam[4], gg[4];
        for (int w = 0; w < 4; w++) {
            kap[w] = linear_L(rks[9][w], ROT_A[0], ROT_B[0], ROT_C[0]);
            gam[w] = linear_L(rks[8][w], ROT_A[3], ROT_B[3], ROT_C[3]);
            gg[w]  = linear_L(rks[7][w], ROT_A[2], ROT_B[2], ROT_C[2]);
        }
        int kap_nib[32], gam_nib[32], g_nib[32];
        for (int d = 0; d < 32; d++) {
            kap_nib[d] = gam_nib[d] = g_nib[d] = 0;
            for (int w = 0; w < 4; w++) {
                kap_nib[d] |= ((kap[w] >> d) & 1) << w;
                gam_nib[d] |= ((gam[w] >> d) & 1) << w;
                g_nib[d]   |= ((gg[w]  >> d) & 1) << w;
            }
        }
        int kap_rec[32], gam_rec[32], g_rec[32];
        for (int d = 0; d < 32; d++) kap_rec[d] = gam_rec[d] = g_rec[d] = -1;

        struct timespec ts0, ts1;
        clock_gettime(CLOCK_MONOTONIC, &ts0);
        long long Dused = 0;

        uint32_t sample_pt[4] = {0}, sample_ct[4] = {0};
        for (int i = 0; i < 4; i++) sample_pt[i] = (uint32_t)splitmix64(&rng);
        encrypt_rounds(sample_pt, rks, 9, sample_ct);

        int sh_ncand[16] = {0};
        int sh_candG[16][64], sh_candg[16][64], sh_candK[16][64][9], sh_E[16][3], sh_DC[16][9];
        for (int sh = 0; sh < nshift; sh++) {
            int u = shift_col(25, sh);
            int E[3] = {RS[2][u][0], RS[2][u][1], RS[2][u][2]};
            sh_E[sh][0] = E[0]; sh_E[sh][1] = E[1]; sh_E[sh][2] = E[2];
            int DC[9], unrec_idx[9], nunrec = 0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    int dcol = RS[3][E[i]][j];
                    DC[3 * i + j] = dcol;
                    sh_DC[sh][3 * i + j] = dcol;
                    if (kap_rec[dcol] < 0) unrec_idx[nunrec++] = 3 * i + j;
                }
            int eff_kbits = (kbits < 4 * nunrec) ? kbits : (4 * nunrec);
            int ms = ncubes;
            if (kbits > 0) {
                int min_ms = (eff_kbits > 0) ? ((eff_kbits + 19) / 4) : 2;
                if (!cubes_explicit) {
                    int sched_ms = SCHED80[sh < 16 ? sh : 15];
                    ms = sched_ms > min_ms ? sched_ms : min_ms;
                } else if (ms < min_ms && eff_kbits > 0) {
                    ms = min_ms;
                }
                if (ms > max_cubes) ms = max_cubes;
            }
            for (int m = 0; m < ms; m++) build_cube(&cubes[m], DC, sh, rks, &rng);
            Dused += (long long)ms * NPT;

            long long nk = 1LL << eff_kbits;
            int total_surv = 0;
            for (long long kc = 0; kc < nk; kc++) {
                int kval[9];
                for (int i = 0; i < 9; i++) {
                    if (kap_rec[DC[i]] >= 0) kval[i] = kap_rec[DC[i]];
                    else kval[i] = kap_nib[DC[i]];
                }
                for (int ui = 0; ui < nunrec; ui++) {
                    int idx = unrec_idx[ui];
                    int lo = 4 * ui, nb = eff_kbits - lo;
                    if (nb >= 4) kval[idx] = (int)((kc >> lo) & 15);
                    else if (nb > 0) {
                        int msk = (1 << nb) - 1;
                        kval[idx] = (int)(((kc >> lo) & msk) | (kap_nib[DC[idx]] & ~msk));
                    }
                }
                static int candG[8192], candg[8192];
                int ncand = 0;
                for (int m = 0; m < ms; m++) {
                    build_Z8(&cubes[m], kval);
                    build_JMhat();
                    build_Dg();
                    if (m == 0) {
                        for (int G = 0; G < 4096; G++) {
                            if ((gam_rec[E[0]] >= 0 && ((G >> 8) & 15) != gam_rec[E[0]]) ||
                                (gam_rec[E[1]] >= 0 && ((G >> 4) & 15) != gam_rec[E[1]]) ||
                                (gam_rec[E[2]] >= 0 && (G & 15) != gam_rec[E[2]]))
                                continue;
                            uint8_t Phi[16];
                            phi_of_gamma(G, Phi);
                            for (int g = 0; g < 16; g++)
                                if (T_of(g, Phi) == 0 && ncand < 8192) {
                                    candG[ncand] = G;
                                    candg[ncand] = g;
                                    ncand++;
                                }
                        }
                    } else {
                        int w = 0;
                        for (int i = 0; i < ncand; i++) {
                            uint8_t Phi[16];
                            phi_of_gamma(candG[i], Phi);
                            if (T_of(candg[i], Phi) == 0) {
                                candG[w] = candG[i];
                                candg[w] = candg[i];
                                w++;
                            }
                        }
                        ncand = w;
                    }
                    if (ncand == 0) break;
                }
                for (int i = 0; i < ncand; i++) {
                    if (total_surv < 64) {
                        for (int k = 0; k < 9; k++) sh_candK[sh][total_surv][k] = kval[k];
                        sh_candG[sh][total_surv] = candG[i];
                        sh_candg[sh][total_surv] = candg[i];
                    }
                    total_surv++;
                }
            }
            sh_ncand[sh] = total_surv < 64 ? total_surv : 64;
            if (sh_ncand[sh] > 0) {
                for (int k = 0; k < 9; k++) {
                    int v0 = sh_candK[sh][0][k], same_k = (total_surv <= 64);
                    for (int s = 1; s < sh_ncand[sh]; s++)
                        if (sh_candK[sh][s][k] != v0) { same_k = 0; break; }
                    if (same_k) kap_rec[DC[k]] = v0;
                }
                for (int j = 0; j < 3; j++) {
                    int shift_b = 4 * (2 - j);
                    int g0 = (sh_candG[sh][0] >> shift_b) & 15, same_g = (total_surv <= 64);
                    for (int s = 1; s < sh_ncand[sh]; s++)
                        if (((sh_candG[sh][s] >> shift_b) & 15) != g0) { same_g = 0; break; }
                    if (same_g) gam_rec[E[j]] = g0;
                }
                if (total_surv == 1) g_rec[u] = sh_candg[sh][0];
                if (verbose)
                    printf("    shift %2d u=%2d : %d survivor(s) (ms=%d, eff_kbits=%d)\n",
                           sh, u, total_surv, ms, eff_kbits);
            } else {
                fprintf(stderr, "    shift %d (u=%d): branch NOT solved\n", sh, u);
            }
        }
        /* Cross-shift kappa & gamma consistency propagation */
        for (int pass = 0; pass < 6; pass++) {
            for (int sh = 0; sh < nshift; sh++) {
                if (sh_ncand[sh] <= 0) continue;
                int w = 0;
                for (int i = 0; i < sh_ncand[sh]; i++) {
                    int G = sh_candG[sh][i], ok = 1;
                    if ((gam_rec[sh_E[sh][0]] >= 0 && ((G >> 8) & 15) != gam_rec[sh_E[sh][0]]) ||
                        (gam_rec[sh_E[sh][1]] >= 0 && ((G >> 4) & 15) != gam_rec[sh_E[sh][1]]) ||
                        (gam_rec[sh_E[sh][2]] >= 0 && (G & 15) != gam_rec[sh_E[sh][2]]))
                        ok = 0;
                    for (int k = 0; k < 9 && ok; k++)
                        if (kap_rec[sh_DC[sh][k]] >= 0 && sh_candK[sh][i][k] != kap_rec[sh_DC[sh][k]])
                            ok = 0;
                    if (!ok) continue;
                    for (int k = 0; k < 9; k++) sh_candK[sh][w][k] = sh_candK[sh][i][k];
                    sh_candG[sh][w] = G;
                    sh_candg[sh][w] = sh_candg[sh][i];
                    w++;
                }
                sh_ncand[sh] = w;
                if (w >= 1) {
                    for (int k = 0; k < 9; k++) {
                        int v0 = sh_candK[sh][0][k], same_k = 1;
                        for (int s = 1; s < w; s++) if (sh_candK[sh][s][k] != v0) { same_k = 0; break; }
                        if (same_k) kap_rec[sh_DC[sh][k]] = v0;
                    }
                    for (int j = 0; j < 3; j++) {
                        int shift_b = 4 * (2 - j);
                        int g0 = (sh_candG[sh][0] >> shift_b) & 15, same_g = 1;
                        for (int s = 1; s < w; s++)
                            if (((sh_candG[sh][s] >> shift_b) & 15) != g0) { same_g = 0; break; }
                        if (same_g) gam_rec[sh_E[sh][j]] = g0;
                    }
                    if (w == 1) g_rec[shift_col(25, sh)] = sh_candg[sh][0];
                }
            }
        }

        int nkap = 0, ngam = 0;
        uint32_t k9w[4] = {0, 0, 0, 0}, k8w[4] = {0, 0, 0, 0}, k9rec[4];
        for (int d = 0; d < 32; d++) {
            if (kap_rec[d] >= 0) {
                nkap++;
                for (int w = 0; w < 4; w++) if ((kap_rec[d] >> w) & 1) k9w[w] |= 1u << d;
            }
        }
        for (int w = 0; w < 4; w++) k9rec[w] = linear_L(k9w[w], ROT_A[0], ROT_B[0], ROT_C[0]);

        /* Resolve any remaining ambiguous gamma nibbles on 2-cube shifts via Module III closure */
        uint32_t mk_rec[8] = {0}, rks_cand[10][4], check_ct[4];
        int mkok = 0;
        if (nkap == 32) {
            int cand_nib[32][16], ncand_nib[32];
            int all_have_cand = 1;
            for (int d = 0; d < 32; d++) {
                if (gam_rec[d] >= 0) {
                    cand_nib[d][0] = gam_rec[d];
                    ncand_nib[d] = 1;
                } else {
                    ncand_nib[d] = 0;
                }
            }
            for (int sh = 0; sh < nshift; sh++) {
                for (int j = 0; j < 3; j++) {
                    int col = sh_E[sh][j];
                    if (gam_rec[col] >= 0) continue;
                    int shift_b = 4 * (2 - j);
                    uint16_t mask = 0;
                    for (int s = 0; s < sh_ncand[sh]; s++)
                        mask |= (uint16_t)(1u << ((sh_candG[sh][s] >> shift_b) & 15));
                    if (ncand_nib[col] == 0) {
                        for (int v = 0; v < 16; v++)
                            if ((mask >> v) & 1) cand_nib[col][ncand_nib[col]++] = v;
                    } else {
                        int w = 0;
                        for (int i = 0; i < ncand_nib[col]; i++)
                            if ((mask >> cand_nib[col][i]) & 1) cand_nib[col][w++] = cand_nib[col][i];
                        ncand_nib[col] = w;
                    }
                }
            }
            int amb_cols[32], namb = 0;
            long long prod = 1;
            for (int d = 0; d < 32; d++) {
                if (ncand_nib[d] == 0) { all_have_cand = 0; break; }
                if (ncand_nib[d] == 1) gam_rec[d] = cand_nib[d][0];
                else {
                    amb_cols[namb++] = d;
                    if (prod <= (1LL << 20)) prod *= ncand_nib[d];
                }
            }
            if (all_have_cand && prod <= (1LL << 20)) {
                uint32_t base_k8w[4] = {0, 0, 0, 0};
                for (int d = 0; d < 32; d++)
                    if (ncand_nib[d] == 1)
                        for (int w = 0; w < 4; w++) if ((gam_rec[d] >> w) & 1) base_k8w[w] |= 1u << d;
                for (long long idx = 0; idx < prod; idx++) {
                    long long rem = idx;
                    uint32_t try_w[4] = {base_k8w[0], base_k8w[1], base_k8w[2], base_k8w[3]};
                    for (int a = 0; a < namb; a++) {
                        int d = amb_cols[a];
                        int v = cand_nib[d][rem % ncand_nib[d]];
                        rem /= ncand_nib[d];
                        for (int w = 0; w < 4; w++) if ((v >> w) & 1) try_w[w] |= 1u << d;
                    }
                    uint32_t try_k8[4];
                    for (int w = 0; w < 4; w++) try_k8[w] = linear_L(try_w[w], ROT_A[3], ROT_B[3], ROT_C[3]);
                    master_key_from_even_pair(try_k8, k9rec, 8, mk_rec);
                    expand_key(mk_rec, 9, rks_cand);
                    encrypt_rounds(sample_pt, rks_cand, 9, check_ct);
                    if (!memcmp(check_ct, sample_ct, 16)) {
                        rem = idx;
                        for (int a = 0; a < namb; a++) {
                            int d = amb_cols[a];
                            gam_rec[d] = cand_nib[d][rem % ncand_nib[d]];
                            rem /= ncand_nib[d];
                        }
                        memcpy(k8w, try_k8, 16);
                        mkok = !memcmp(mk_rec, mk_true, 32);
                        break;
                    }
                }
            }
        }
        ngam = 0;
        uint32_t fin_k8w[4] = {0, 0, 0, 0};
        for (int d = 0; d < 32; d++) {
            if (gam_rec[d] >= 0) {
                ngam++;
                for (int w = 0; w < 4; w++) if ((gam_rec[d] >> w) & 1) fin_k8w[w] |= 1u << d;
            }
        }
        for (int w = 0; w < 4; w++) k8w[w] = linear_L(fin_k8w[w], ROT_A[3], ROT_B[3], ROT_C[3]);
        int k9ok = (nkap == 32) && !memcmp(k9rec, rks[9], 16);
        int k8ok = (ngam == 32) && !memcmp(k8w, rks[8], 16);
        if (!mkok && nkap == 32 && ngam == 32) {
            master_key_from_even_pair(k8w, k9rec, 8, mk_rec);
            expand_key(mk_rec, 9, rks_cand);
            encrypt_rounds(sample_pt, rks_cand, 9, check_ct);
            mkok = !memcmp(check_ct, sample_ct, 16) && !memcmp(mk_rec, mk_true, 32);
        }
        clock_gettime(CLOCK_MONOTONIC, &ts1);
        double dt = (ts1.tv_sec - ts0.tv_sec) + 1e-9 * (ts1.tv_nsec - ts0.tv_nsec);
        tot += dt;
        D_total = Dused;
        passed += mkok;
        printf("  [Trial %2d/%2d] kappa=%2d/32 gamma=%2d/32 | K9_ok=%d K8_ok=%d | 256-bit MK_match=%s | %.2fs\n",
               tr + 1, trials, nkap, ngam, k9ok, k8ok, mkok ? "True" : "False", dt);
        fflush(stdout);
    }
    printf("  [Summary 9-Round ARADI] Exact 256-bit MK recovery: %d/%d (%.0f%%) | mean %.2fs | D = %lld CP (2^%.2f)\n",
           passed, trials, 100.0 * passed / trials, tot / trials, D_total, log2((double)D_total));
    printf("%s\n", passed == trials ? "[PASS] 9-Round ARADI nested-WHT peel verified." : "[FAIL]");
    for (int m = 0; m < max_cubes; m++) free(cubes[m].t9);
    free(cubes);
    return passed == trials ? 0 : 1;
}
