#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static void *alloc64(size_t sz) {
    void *p = NULL;
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L
    if (posix_memalign(&p, 64, sz) == 0) return p;
#endif
    return malloc(sz);
}

static const uint32_t ROT_A[4] = {11, 10, 9, 8};
static const uint32_t ROT_B[4] = {8, 9, 4, 9};
static const uint32_t ROT_C[4] = {14, 11, 14, 7};
static uint8_t SBOX_INV[16];

static inline uint32_t rol32(uint32_t x, uint32_t n) {
    return n ? ((x << n) | (x >> (32 - n))) : x;
}
static inline uint32_t ror32(uint32_t x, uint32_t n) {
    return n ? ((x >> n) | (x << (32 - n))) : x;
}
static inline uint16_t rol16(uint16_t x, uint32_t n) {
    return n ? (uint16_t)((x << n) | (x >> (16 - n))) : x;
}
static inline uint32_t linear_L(uint32_t s, uint32_t a, uint32_t b, uint32_t c) {
    uint16_t u = (uint16_t)(s >> 16);
    uint16_t l = (uint16_t)(s & 0xFFFF);
    uint16_t uo = u ^ rol16(u, a) ^ rol16(l, c);
    uint16_t lo = l ^ rol16(l, a) ^ rol16(u, b);
    return ((uint32_t)uo << 16) | (uint32_t)lo;
}
static inline void sbox_fwd(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z) {
    *x ^= (*w & *y);
    *z ^= (*x & *y);
    *y ^= (*w & *z);
    *w ^= (*x & *z);
}
static inline void sbox_bwd(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z) {
    *w ^= (*x & *z);
    *y ^= (*w & *z);
    *z ^= (*x & *y);
    *x ^= (*w & *y);
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
        for (int w = 0; w < 4; w++) rks[i][w] = K[4*j + w];
        uint32_t n1, n0, n3, n2, n5, n4, n7, n6, tmp;
        key_M(1, 3, K[1], K[0], &n1, &n0); K[1] = n1; K[0] = n0;
        key_M(9, 28, K[3], K[2], &n3, &n2); K[3] = n3; K[2] = n2;
        key_M(1, 3, K[5], K[4], &n5, &n4); K[5] = n5; K[4] = n4;
        key_M(9, 28, K[7], K[6], &n7, &n6); K[7] = n7 ^ (uint32_t)i; K[6] = n6;
        if (j == 0) {
            tmp = K[1]; K[1] = K[2]; K[2] = tmp;
            tmp = K[5]; K[5] = K[6]; K[6] = tmp;
        } else {
            tmp = K[1]; K[1] = K[4]; K[4] = tmp;
            tmp = K[3]; K[3] = K[6]; K[6] = tmp;
        }
    }
    int j = R & 1;
    for (int w = 0; w < 4; w++) rks[R][w] = K[4*j + w];
}

static void master_key_from_even_pair(const uint32_t keven[4], const uint32_t kodd[4], int even_idx, uint32_t mk_out[8]) {
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
        key_M_inv(1, 3, S[1], S[0], &x, &y); S[1] = x; S[0] = y;
        key_M_inv(9, 28, S[3], S[2], &x, &y); S[3] = x; S[2] = y;
        key_M_inv(1, 3, S[5], S[4], &x, &y); S[5] = x; S[4] = y;
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
    ct[0] = w ^ rks[R][0]; ct[1] = x ^ rks[R][1]; ct[2] = y ^ rks[R][2]; ct[3] = z ^ rks[R][3];
}

static void fht12(int32_t *a) {
    for (int h = 1; h < 4096; h <<= 1) {
        for (int i = 0; i < 4096; i += (h << 1)) {
            for (int j = i; j < i + h; j++) {
                int32_t x = a[j], y = a[j + h];
                a[j] = x + y;
                a[j + h] = x - y;
            }
        }
    }
}

static int32_t F_KERNEL[16][4][4096];
static void init_fht_kernel(void) {
    for (int v = 0; v < 16; v++) {
        uint32_t w = v & 1, x = (v >> 1) & 1, y = (v >> 2) & 1, z = (v >> 3) & 1;
        sbox_bwd(&w, &x, &y, &z);
        SBOX_INV[v] = (uint8_t)(w | (x << 1) | (y << 2) | (z << 3));
    }
    for (int g = 0; g < 16; g++) {
        for (int r = 0; r < 4; r++) {
            for (int x = 0; x < 4096; x++) {
                uint8_t s = SBOX_INV[x & 15] ^ SBOX_INV[(x >> 4) & 15] ^ SBOX_INV[(x >> 8) & 15];
                uint8_t out = SBOX_INV[s ^ g];
                F_KERNEL[g][r][x] = (out >> r) & 1;
            }
            fht12(F_KERNEL[g][r]);
        }
    }
}

/* Compute 4*M-bit signature table F_out[65536] for branch u with lam8 deps (v0, v1, v2) across M 4096-pt cubes! */
static void build_branch_sigs(uint8_t (*T8_nibs)[4096][32], int M, int v0, int v1, int v2, uint64_t F_out[65536]) {
    memset(F_out, 0, 65536 * sizeof(uint64_t));
    int32_t A[4096], tmp[4096];
    for (int m = 0; m < M; m++) {
        memset(A, 0, sizeof(A));
        for (int p = 0; p < 4096; p++) {
            int idx = T8_nibs[m][p][v0] | (T8_nibs[m][p][v1] << 4) | (T8_nibs[m][p][v2] << 8);
            A[idx] ^= 1;
        }
        fht12(A);
        for (int g7 = 0; g7 < 16; g7++) {
            for (int r = 0; r < 4; r++) {
                for (int x = 0; x < 4096; x++) tmp[x] = A[x] * F_KERNEL[g7][r][x];
                fht12(tmp);
                uint64_t bit_mask = 1ULL << (4*m + r);
                int base_idx = g7 << 12;
                for (int x = 0; x < 4096; x++) {
                    if ((tmp[x] >> 12) & 1) {
                        F_out[base_idx | x] |= bit_mask;
                    }
                }
            }
        }
    }
}

static uint64_t xorshift64(uint64_t *st) {
    uint64_t x = *st;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return (*st = x);
}

static int shift_col(int c, int s) {
    return (c < 16) ? ((c + s) & 15) : (16 + (((c - 16) + s) & 15));
}

static void get_deps(int lin_idx, int target_col, int out3[3]) {
    int k = 0;
    for (int inp = 0; inp < 32; inp++) {
        uint32_t y = linear_L(1u << inp, ROT_A[lin_idx], ROT_B[lin_idx], ROT_C[lin_idx]);
        if ((y >> target_col) & 1u) out3[k++] = inp;
    }
}

/* Hash table for 65536 entries of F_3 */
#define HT_BITS 18
#define HT_SIZE (1 << HT_BITS)
#define HT_MASK (HT_SIZE - 1)
typedef struct { uint64_t sig; int32_t k3; int32_t next; } ht_node_t;

int main(int argc, char **argv) {
    init_fht_kernel();
    const int M0 = 11; /* 11 * 4096 = 45,056 CP at shift s=0 */
    const int M8 = 11; /* 11 * 4096 = 45,056 CP at shift s=8 -> TOTAL = 22 * 4096 = 90,112 CP (2^16.46 CP)! */
    uint8_t (*T8_s0)[4096][32] = alloc64(M0 * sizeof(*T8_s0));
    uint8_t (*T8_s8)[4096][32] = alloc64(M8 * sizeof(*T8_s8));
    uint32_t (*CT_s0)[4096][4] = alloc64(M0 * sizeof(*CT_s0));
    uint32_t (*CT_s8)[4096][4] = alloc64(M8 * sizeof(*CT_s8));
    uint64_t *bloom = alloc64((1 << 18) * sizeof(uint64_t)); /* 2^24-bit = 2 MB L3-resident bitset! */
    int32_t *ht_head = malloc(HT_SIZE * sizeof(int32_t));
    ht_node_t *ht_nodes = malloc(65536 * sizeof(ht_node_t));

    uint64_t rng = 0x202611068888ULL;
    int trials = 10, passed = 0;
    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--trials") == 0 && a + 1 < argc) trials = atoi(argv[++a]);
        else if (strcmp(argv[a], "--seed") == 0 && a + 1 < argc) rng = strtoull(argv[++a], NULL, 0) ^ 0x9E3779B97F4A7C15ULL;
    }
    double total_solve_sec = 0.0;

    printf("========================================================================\n");
    printf("  8-Round ARADI Practical Full 256-Bit Master-Key Recovery (90,112 CP)\n");
    printf("========================================================================\n\n");

    for (int tr = 0; tr < trials; tr++) {
        uint32_t mk_true[8], rks[9][4];
        for (int i = 0; i < 8; i++) mk_true[i] = (uint32_t)xorshift64(&rng);
        expand_key(mk_true, 8, rks);

        uint32_t sample_pt[4], sample_ct[4];
        /* Collect M0 cubes at s=0 [cols (0, 5, 19)] and M8 cubes at s=8 [cols (8, 13, 27)] */
        for (int phase = 0; phase < 2; phase++) {
            int s = (phase == 0) ? 0 : 8;
            int M = (phase == 0) ? M0 : M8;
            int c0 = shift_col(0, s), c1 = shift_col(5, s), c2 = shift_col(19, s);
            int cols[3] = {c0, c1, c2};
            uint8_t (*T8)[4096][32] = (phase == 0) ? T8_s0 : T8_s8;
            uint32_t (*CT)[4096][4] = (phase == 0) ? CT_s0 : CT_s8;
            for (int m = 0; m < M; m++) {
                uint32_t base[4] = {(uint32_t)xorshift64(&rng), (uint32_t)xorshift64(&rng),
                                    (uint32_t)xorshift64(&rng), (uint32_t)xorshift64(&rng)};
                for (int u = 0; u < 4096; u++) {
                    uint32_t pt[4] = {base[0], base[1], base[2], base[3]};
                    for (int ci = 0; ci < 3; ci++) {
                        int nib = (u >> (4*ci)) & 15;
                        for (int w = 0; w < 4; w++) {
                            if ((nib >> w) & 1) pt[w] ^= (1u << cols[ci]);
                        }
                    }
                    uint32_t ct[4];
                    encrypt_rounds(pt, rks, 8, ct);
                    if (phase == 0 && m == 0 && u == 0) {
                        memcpy(sample_pt, pt, 16);
                        memcpy(sample_ct, ct, 16);
                    }
                    memcpy(CT[m][u], ct, 16);
                    /* Apply Lambda_3 (since round 7 uses L_{7 mod 4} = L_3) */
                    uint32_t tw[4];
                    for (int w = 0; w < 4; w++) tw[w] = linear_L(ct[w], ROT_A[3], ROT_B[3], ROT_C[3]);
                    for (int b = 0; b < 32; b++) {
                        T8[m][u][b] = (uint8_t)(((tw[0] >> b) & 1u) | (((tw[1] >> b) & 1u) << 1) |
                                                (((tw[2] >> b) & 1u) << 2) | (((tw[3] >> b) & 1u) << 3));
                    }
                }
            }
        }

        struct timespec ts0, ts1;
        clock_gettime(CLOCK_MONOTONIC, &ts0);

        /* === ONLINE 8-ROUND KEY RECOVERY === */
        /* Step 1 (Shift s=0): Solve b=25 -> filter with b=31 -> filter with b=18 -> filter with b=10! */
        int b_s0[4] = {25, 31, 18, 10};
        static uint64_t F_s0[4][3][65536];
        int u_s0[4][3], v_s0[4][3][3];
        for (int bi = 0; bi < 4; bi++) {
            get_deps(1, b_s0[bi], u_s0[bi]);
            for (int br = 0; br < 3; br++) {
                get_deps(2, u_s0[bi][br], v_s0[bi][br]);
                build_branch_sigs(T8_s0, M0, v_s0[bi][br][0], v_s0[bi][br][1], v_s0[bi][br][2], F_s0[bi][br]);
            }
        }

        /* Index F_s0[0][2] (b=25, branch 2: u=31, v=[1,22,31]) in 2MB 24-bit Bloom + hash table */
        memset(bloom, 0, (1 << 18) * sizeof(uint64_t));
        memset(ht_head, -1, HT_SIZE * sizeof(int32_t));
        for (int k3 = 0; k3 < 65536; k3++) {
            uint64_t sig = F_s0[0][2][k3];
            uint32_t lo24 = (uint32_t)(sig & 0xFFFFFFULL);
            bloom[lo24 >> 6] |= (1ULL << (lo24 & 63));
            uint32_t h = (uint32_t)((sig * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
            ht_nodes[k3].sig = sig;
            ht_nodes[k3].k3 = k3;
            ht_nodes[k3].next = ht_head[h];
            ht_head[h] = k3;
        }

        /* Also index F_s0[1][1] (b=31, branch 1: u=21, v=[7,21,28]) in a second hash table for instant O(256) verification! */
        int32_t *ht2_head = malloc(HT_SIZE * sizeof(int32_t));
        ht_node_t *ht2_nodes = malloc(65536 * sizeof(ht_node_t));
        memset(ht2_head, -1, HT_SIZE * sizeof(int32_t));
        for (int k = 0; k < 65536; k++) {
            uint64_t sig = F_s0[1][1][k];
            uint32_t h = (uint32_t)((sig * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
            ht2_nodes[k].sig = sig; ht2_nodes[k].k3 = k; ht2_nodes[k].next = ht2_head[h]; ht2_head[h] = k;
        }

        int lam8[32];
        for (int i = 0; i < 32; i++) lam8[i] = -1;

        /* 2^32 L3-resident Bloom scan on b=25: F_14[k1] ^ F_25[k2] == F_31[k3] */
        const uint64_t *F1 = F_s0[0][0], *F2 = F_s0[0][1];
        int found_s0 = 0;
        for (int k1 = 0; k1 < 65536 && !found_s0; k1++) {
            uint64_t s1 = F1[k1];
            for (int k2 = 0; k2 < 65536; k2++) {
                uint64_t target = s1 ^ F2[k2];
                uint32_t lo24 = (uint32_t)(target & 0xFFFFFFULL);
                if (!((bloom[lo24 >> 6] >> (lo24 & 63)) & 1ULL)) continue;
                uint32_t h = (uint32_t)((target * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
                for (int idx = ht_head[h]; idx != -1; idx = ht_nodes[idx].next) {
                    if (ht_nodes[idx].sig != target) continue;
                    int k3 = ht_nodes[idx].k3;
                    /* Cross-filter b=25 candidate against b=31 (u_s0[1]=[4,21,31],
                       v_s0[1]=[[4,11,16], [7,21,28], [1,22,31]]): branch 2 reuses k3
                       and branch 0 shares lam8[11], lam8[16] from k2. */
                    int l11 = k2 & 15, l16 = (k2 >> 4) & 15;
                    uint64_t s_b31_br2 = F_s0[1][2][k3]; /* shared u=31, v=[1,22,31] */
                    for (int g7_4 = 0; g7_4 < 16; g7_4++) {
                        for (int l4 = 0; l4 < 16; l4++) {
                            int idx_b31_0 = (g7_4 << 12) | l4 | (l11 << 4) | (l16 << 8);
                            uint64_t need21 = F_s0[1][0][idx_b31_0] ^ s_b31_br2;
                            uint32_t h2 = (uint32_t)((need21 * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
                            for (int idx2 = ht2_head[h2]; idx2 != -1; idx2 = ht2_nodes[idx2].next) {
                                if (ht2_nodes[idx2].sig != need21) continue;
                                int k21 = ht2_nodes[idx2].k3;
                                /* Cascade to b=18 (v=[[7,14,19], [4,18,25], [10,24,31]])
                                   using known nibbles l7, l14, l4, l25, l31. */
                                int l7 = k21 & 15, l14 = (k1 >> 4) & 15, l25 = (k2 >> 8) & 15, l31 = (k3 >> 8) & 15;
                                /* Match remaining unknowns across the three branches of b=18 */
                                for (int g18_0 = 0; g18_0 < 16 && !found_s0; g18_0++) {
                                    for (int l19 = 0; l19 < 16 && !found_s0; l19++) {
                                        uint64_t sb0 = F_s0[2][0][(g18_0 << 12) | l7 | (l14 << 4) | (l19 << 8)];
                                        for (int g18_1 = 0; g18_1 < 16 && !found_s0; g18_1++) {
                                            for (int l18 = 0; l18 < 16 && !found_s0; l18++) {
                                                uint64_t need_br2 = sb0 ^ F_s0[2][1][(g18_1 << 12) | l4 | (l18 << 4) | (l25 << 8)];
                                                for (int g18_2 = 0; g18_2 < 16; g18_2++) {
                                                    for (int l10 = 0; l10 < 16; l10++) {
                                                        for (int l24 = 0; l24 < 16; l24++) {
                                                            if (F_s0[2][2][(g18_2 << 12) | l10 | (l24 << 4) | (l31 << 8)] == need_br2) {
                                                                /* Cascade to b=10 (v=[[0,7,28], [1,10,22], [3,17,24]])
                                                                   using known nibbles l7, l28, l1, l10, l22, l24. */
                                                                int l28 = (k21 >> 8) & 15, l1 = k3 & 15, l22 = (k3 >> 4) & 15;
                                                                for (int g10_1 = 0; g10_1 < 16 && !found_s0; g10_1++) {
                                                                    uint64_t s10_1 = F_s0[3][1][(g10_1 << 12) | l1 | (l10 << 4) | (l22 << 8)];
                                                                    for (int g10_0 = 0; g10_0 < 16 && !found_s0; g10_0++) {
                                                                        for (int l0 = 0; l0 < 16 && !found_s0; l0++) {
                                                                            uint64_t need10_2 = s10_1 ^ F_s0[3][0][(g10_0 << 12) | l0 | (l7 << 4) | (l28 << 8)];
                                                                            for (int g10_2 = 0; g10_2 < 16; g10_2++) {
                                                                                for (int l3 = 0; l3 < 16; l3++) {
                                                                                    for (int l17 = 0; l17 < 16; l17++) {
                                                                                        if (F_s0[3][2][(g10_2 << 12) | l3 | (l17 << 4) | (l24 << 8)] == need10_2) {
                                                                                            /* Record all 20 recovered columns at shift s=0 */
                                                                                            lam8[5] = k1 & 15; lam8[14] = l14; lam8[26] = (k1 >> 8) & 15;
                                                                                            lam8[11] = l11; lam8[16] = l16; lam8[25] = l25;
                                                                                            lam8[1] = l1; lam8[22] = l22; lam8[31] = l31;
                                                                                            lam8[4] = l4; lam8[7] = l7; lam8[21] = (k21 >> 4) & 15; lam8[28] = l28;
                                                                                            lam8[19] = l19; lam8[18] = l18; lam8[10] = l10; lam8[24] = l24;
                                                                                            lam8[0] = l0; lam8[3] = l3; lam8[17] = l17;
                                                                                            found_s0 = 1;
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        free(ht2_head); free(ht2_nodes);

        /* Step 2 (Shift s=8): Jointly solve remaining 12 columns of lam8 across b in {17, 23, 26, 2} (176-bit filter on 96 bits -> 2^-80 FP)! */
        int b_s8[4] = {17, 23, 26, 2};
        int u8_all[4][3], v8_all[4][3][3];
        static uint64_t F8_all[4][3][65536];
        for (int bi = 0; bi < 4; bi++) {
            get_deps(1, b_s8[bi], u8_all[bi]);
            for (int br = 0; br < 3; br++) {
                get_deps(2, u8_all[bi][br], v8_all[bi][br]);
                build_branch_sigs(T8_s8, M8, v8_all[bi][br][0], v8_all[bi][br][1], v8_all[bi][br][2], F8_all[bi][br]);
            }
        }
        /* Enumerate all b=17 solutions (~1.0002 expected) and require consistency across b=23, 26, 2 */
        int base_lam8[32];
        memcpy(base_lam8, lam8, sizeof(base_lam8));
        int v0_cnt[3] = {0, 0, 0}, *v0_list[3];
        for (int br = 0; br < 3; br++) {
            v0_list[br] = malloc(65536 * sizeof(int));
            for (int k = 0; k < 65536; k++) {
                int n0 = k & 15, n1 = (k >> 4) & 15, n2 = (k >> 8) & 15;
                if (base_lam8[v8_all[0][br][0]] != -1 && base_lam8[v8_all[0][br][0]] != n0) continue;
                if (base_lam8[v8_all[0][br][1]] != -1 && base_lam8[v8_all[0][br][1]] != n1) continue;
                if (base_lam8[v8_all[0][br][2]] != -1 && base_lam8[v8_all[0][br][2]] != n2) continue;
                v0_list[br][v0_cnt[br]++] = k;
            }
        }
        memset(ht_head, -1, HT_SIZE * sizeof(int32_t));
        for (int i = 0; i < v0_cnt[2]; i++) {
            int k = v0_list[2][i];
            uint64_t sig = F8_all[0][2][k];
            uint32_t h = (uint32_t)((sig * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
            ht_nodes[i].sig = sig; ht_nodes[i].k3 = k; ht_nodes[i].next = ht_head[h]; ht_head[h] = i;
        }
        int found_s8 = 0;
        for (int i0 = 0; i0 < v0_cnt[0] && !found_s8; i0++) {
            int k0 = v0_list[0][i0];
            uint64_t s0_v = F8_all[0][0][k0];
            for (int i1 = 0; i1 < v0_cnt[1] && !found_s8; i1++) {
                int k1 = v0_list[1][i1];
                uint64_t need0 = s0_v ^ F8_all[0][1][k1];
                uint32_t h0 = (uint32_t)((need0 * 0x9E3779B97F4A7C15ULL) >> (64 - HT_BITS));
                for (int idx0 = ht_head[h0]; idx0 != -1 && !found_s8; idx0 = ht_nodes[idx0].next) {
                    if (ht_nodes[idx0].sig != need0) continue;
                    int cand_lam8[32];
                    memcpy(cand_lam8, base_lam8, sizeof(cand_lam8));
                    int ks0[3] = {k0, k1, ht_nodes[idx0].k3};
                    for (int br = 0; br < 3; br++) {
                        cand_lam8[v8_all[0][br][0]] = ks0[br] & 15;
                        cand_lam8[v8_all[0][br][1]] = (ks0[br] >> 4) & 15;
                        cand_lam8[v8_all[0][br][2]] = (ks0[br] >> 8) & 15;
                    }
                    int ok_chain = 1;
                    for (int bi = 1; bi < 4 && ok_chain; bi++) {
                        int vc[3] = {0, 0, 0}, *vl[3];
                        for (int br = 0; br < 3; br++) {
                            vl[br] = malloc(65536 * sizeof(int));
                            for (int k = 0; k < 65536; k++) {
                                int n0 = k & 15, n1 = (k >> 4) & 15, n2 = (k >> 8) & 15;
                                if (cand_lam8[v8_all[bi][br][0]] != -1 && cand_lam8[v8_all[bi][br][0]] != n0) continue;
                                if (cand_lam8[v8_all[bi][br][1]] != -1 && cand_lam8[v8_all[bi][br][1]] != n1) continue;
                                if (cand_lam8[v8_all[bi][br][2]] != -1 && cand_lam8[v8_all[bi][br][2]] != n2) continue;
                                vl[br][vc[br]++] = k;
                            }
                        }
                        int found_bi = 0;
                        for (int a0 = 0; a0 < vc[0] && !found_bi; a0++) {
                            uint64_t sa0 = F8_all[bi][0][vl[0][a0]];
                            for (int a1 = 0; a1 < vc[1] && !found_bi; a1++) {
                                uint64_t need_bi = sa0 ^ F8_all[bi][1][vl[1][a1]];
                                for (int a2 = 0; a2 < vc[2]; a2++) {
                                    if (F8_all[bi][2][vl[2][a2]] == need_bi) {
                                        int ks[3] = {vl[0][a0], vl[1][a1], vl[2][a2]};
                                        for (int br = 0; br < 3; br++) {
                                            cand_lam8[v8_all[bi][br][0]] = ks[br] & 15;
                                            cand_lam8[v8_all[bi][br][1]] = (ks[br] >> 4) & 15;
                                            cand_lam8[v8_all[bi][br][2]] = (ks[br] >> 8) & 15;
                                        }
                                        found_bi = 1; break;
                                    }
                                }
                            }
                        }
                        for (int br = 0; br < 3; br++) free(vl[br]);
                        if (!found_bi) ok_chain = 0;
                    }
                    if (ok_chain) {
                        memcpy(lam8, cand_lam8, sizeof(lam8));
                        found_s8 = 1;
                    }
                }
            }
        }
        for (int br = 0; br < 3; br++) free(v0_list[br]);

        /* Reconstruct K^(8) = Lambda_3(lam8), peel Round 7 on 2 cubes of s=0 and s=8 to get K^(7) & K^(6), then KS^{-1}(K^(6), K^(7), 6)! */
        uint32_t w8[4] = {0,0,0,0}, k8_rec[4];
        for (int b = 0; b < 32; b++)
            for (int w = 0; w < 4; w++)
                if ((lam8[b] >> w) & 1) w8[w] |= (1u << b);
        for (int w = 0; w < 4; w++) k8_rec[w] = linear_L(w8[w], ROT_A[3], ROT_B[3], ROT_C[3]);

        /* Peel Round 7 across all M0 = 11 cubes at shift s = 0 to form T7_peeled = Lambda_2(C_7). */
        uint8_t (*T7_peeled)[4096][32] = malloc(M0 * sizeof(*T7_peeled));
        uint32_t (*C7_peeled)[4096][4] = malloc(M0 * sizeof(*C7_peeled));
        for (int m = 0; m < M0; m++) {
            for (int u = 0; u < 4096; u++) {
                uint32_t xw[4];
                for (int w = 0; w < 4; w++) xw[w] = linear_L(CT_s0[m][u][w] ^ k8_rec[w], ROT_A[3], ROT_B[3], ROT_C[3]);
                sbox_bwd(&xw[0], &xw[1], &xw[2], &xw[3]);
                memcpy(C7_peeled[m][u], xw, 16);
                uint32_t tw7[4];
                for (int w = 0; w < 4; w++) tw7[w] = linear_L(xw[w], ROT_A[2], ROT_B[2], ROT_C[2]);
                for (int b = 0; b < 32; b++)
                    T7_peeled[m][u][b] = (uint8_t)(((tw7[0]>>b)&1u)|(((tw7[1]>>b)&1u)<<1)|(((tw7[2]>>b)&1u)<<2)|(((tw7[3]>>b)&1u)<<3));
            }
        }
        int lam7[32], lam6[32];
        for (int i = 0; i < 32; i++) { lam7[i] = -1; lam6[i] = 0; }
        /* Pass 1: Recover all 32 nibbles of lam7 across the 26 full-4bit balanced Q_5 columns at s=0 and s=8.
           Check uniqueness of the 12-bit lam7 tuple (k & 0xFFF) so multiset-even columns {10,18,25,31} (where all 16 g6
           guesses vanish for the true 12-bit tuple) and multi-g6 hits on affine S-box outputs are retained cleanly! */
        static const int Q5_BAL[26] = {0,1,2,4,5,6,8,9,10,11,13,14,15,16,17,18,19,20,21,23,24,25,28,29,30,31};
        uint64_t (*F7_all)[65536] = malloc(32 * sizeof(*F7_all));
        uint64_t full44_mask = (1ULL << (4 * M0)) - 1ULL;
        for (int qi = 0; qi < 26; qi++) {
            int cq = Q5_BAL[qi];
            int d7[3];
            get_deps(1, cq, d7);
            build_branch_sigs(T7_peeled, M0, d7[0], d7[1], d7[2], F7_all[cq]);
            int uniq_low12 = -1, conflict = 0;
            for (int k = 0; k < 65536; k++) {
                if ((F7_all[cq][k] & full44_mask) == 0) {
                    int low12 = k & 0xFFF;
                    if (uniq_low12 == -1) uniq_low12 = low12;
                    else if (uniq_low12 != low12) { conflict = 1; break; }
                }
            }
            if (uniq_low12 != -1 && !conflict) {
                lam7[d7[0]] = uniq_low12 & 15;
                lam7[d7[1]] = (uniq_low12 >> 4) & 15;
                lam7[d7[2]] = (uniq_low12 >> 8) & 15;
            }
        }
        /* Reconstruct K^(7) = Lambda_2(lam7), then peel Rounds 7 and 6 on the 256-pt 8D subcubes at s = 0 and s = 8 to recover K^(6). */
        uint32_t w7[4]={0,0,0,0}, w6[4]={0,0,0,0}, k7_rec[4], k6_rec[4], mk_rec[8];
        for (int b = 0; b < 32; b++)
            for (int w = 0; w < 4; w++)
                if ((lam7[b] >> w) & 1) w7[w] |= (1u << b);
        for (int w = 0; w < 4; w++) k7_rec[w] = linear_L(w7[w], ROT_A[2], ROT_B[2], ROT_C[2]);

        /* Peel both s=0 and s=8 on the 256-pt 8D subcube ((u>>4)&15)==0 (cols (0,19) at s=0 and cols (8,27) at s=8!) to C_6 -> T_6 = L_1(C_6)!
           On that 256-pt subcube, Z_4 = L_0(U_5) has 19 full-4bit balanced cols at s=0 AND 19 at s=8, giving a 44-bit filter on just 12 bits (4,096 checks per col, < 0.2 ms total)! */
        static const int Z4_BAL[19] = {0,2,5,6,8,9,10,11,14,16,17,18,20,21,23,25,26,29,31};
        for (int phase6 = 0; phase6 < 2; phase6++) {
            int s = (phase6 == 0) ? 0 : 8;
            uint32_t (*CT8_src)[4096][4] = (phase6 == 0) ? CT_s0 : CT_s8;
            int M6 = (phase6 == 0) ? M0 : M8;
            /* Build S_sig[col][g6] (44-bit XOR sum over 256-pt subcube across 11 cubes!) */
            uint64_t S6_sig[32][16];
            memset(S6_sig, 0, sizeof(S6_sig));
            for (int m = 0; m < M6; m++) { /* 11 * 4 = 44-bit filter on 12-bit key! */
                uint8_t xor_nib[32][16];
                memset(xor_nib, 0, sizeof(xor_nib));
                for (int u0 = 0; u0 < 16; u0++) {
                    for (int u2 = 0; u2 < 16; u2++) {
                        int u = u0 | (u2 << 8); /* col 5 fixed to 0 -> exact 256-pt 8D cube (s, 16+((3+s)&15))! */
                        uint32_t xw[4];
                        for (int w = 0; w < 4; w++) xw[w] = linear_L(CT8_src[m][u][w] ^ k8_rec[w], ROT_A[3], ROT_B[3], ROT_C[3]);
                        sbox_bwd(&xw[0], &xw[1], &xw[2], &xw[3]); /* C_7 */
                        for (int w = 0; w < 4; w++) xw[w] = linear_L(xw[w] ^ k7_rec[w], ROT_A[2], ROT_B[2], ROT_C[2]);
                        sbox_bwd(&xw[0], &xw[1], &xw[2], &xw[3]); /* C_6 */
                        uint32_t tw6[4];
                        for (int w = 0; w < 4; w++) tw6[w] = linear_L(xw[w], ROT_A[1], ROT_B[1], ROT_C[1]); /* T_6 */
                        for (int b = 0; b < 32; b++) {
                            uint8_t tn = (uint8_t)(((tw6[0]>>b)&1u)|(((tw6[1]>>b)&1u)<<1)|(((tw6[2]>>b)&1u)<<2)|(((tw6[3]>>b)&1u)<<3));
                            for (int g = 0; g < 16; g++) xor_nib[b][g] ^= SBOX_INV[tn ^ g];
                        }
                    }
                }
                for (int b = 0; b < 32; b++)
                    for (int g = 0; g < 16; g++)
                        S6_sig[b][g] |= ((uint64_t)xor_nib[b][g]) << (4*m);
            }
            for (int zi = 0; zi < 19; zi++) {
                int bz = shift_col(Z4_BAL[zi], s);
                int d6[3];
                get_deps(0, bz, d6);
                int nsol6 = 0, bg0 = 0, bg1 = 0, bg2 = 0;
                for (int g0 = 0; g0 < 16; g0++) {
                    uint64_t v0 = S6_sig[d6[0]][g0];
                    for (int g1 = 0; g1 < 16; g1++) {
                        uint64_t v01 = v0 ^ S6_sig[d6[1]][g1];
                        for (int g2 = 0; g2 < 16; g2++) {
                            if ((v01 ^ S6_sig[d6[2]][g2]) == 0ULL) {
                                nsol6++;
                                bg0 = g0; bg1 = g1; bg2 = g2;
                            }
                        }
                    }
                }
                if (nsol6 == 1) {
                    lam6[d6[0]] = bg0; lam6[d6[1]] = bg1; lam6[d6[2]] = bg2;
                }
            }
        }
        free(T7_peeled); free(C7_peeled); free(F7_all);
        for (int b = 0; b < 32; b++)
            for (int w = 0; w < 4; w++)
                if ((lam6[b] >> w) & 1) w6[w] |= (1u << b);
        for (int w = 0; w < 4; w++) k6_rec[w] = linear_L(w6[w], ROT_A[1], ROT_B[1], ROT_C[1]);
        master_key_from_even_pair(k6_rec, k7_rec, 6, mk_rec);

        clock_gettime(CLOCK_MONOTONIC, &ts1);
        double dt = (ts1.tv_sec - ts0.tv_sec) + 1e-9 * (ts1.tv_nsec - ts0.tv_nsec);
        total_solve_sec += dt;

        int match = (memcmp(mk_rec, mk_true, 32) == 0);
        uint32_t chk_ct[4], chk_rks[9][4];
        expand_key(mk_rec, 8, chk_rks);
        encrypt_rounds(sample_pt, chk_rks, 8, chk_ct);
        int ct_ok = (memcmp(chk_ct, sample_ct, 16) == 0);
        passed += (match && ct_ok);

        printf("  [Trial %d/%d] 90,112 CP | Online solve = %.3f s | K8_ok=%d K7_ok=%d K6_ok=%d | 256-bit_MK_match=%s (%08x %08x ...)\n",
               tr + 1, trials, dt,
               memcmp(k8_rec, rks[8], 16) == 0,
               memcmp(k7_rec, rks[7], 16) == 0,
               memcmp(k6_rec, rks[6], 16) == 0,
               (match && ct_ok) ? "TRUE" : "FALSE",
               mk_rec[0], mk_rec[1]);
        fflush(stdout);
    }

    printf("\n[SUMMARY] 8-Round ARADI Full 256-Bit Master-Key Recovery: %d/%d (%.0f%%), mean online time = %.3f s (1 core)!\n",
           passed, trials, 100.0 * passed / trials, total_solve_sec / trials);
    free(T8_s0); free(T8_s8); free(CT_s0); free(CT_s8); free(bloom); free(ht_head); free(ht_nodes);
    return (passed == trials) ? 0 : 1;
}
