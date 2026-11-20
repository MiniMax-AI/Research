// Multi-threaded Bit-Parallel ISD Sampler for Minimum-Weight Codewords of ker(H)
// Supports:
//   - Lee-Brickell (p=3) with early 64-bit word weight pruning
//   - Ball-Collision and Stern-Dumer (p = 2+3 = 5 or p = 3+3 = 6) with l-bit collision buckets
//   - Lock-free atomic CAS global hash table for deduplication across OpenMP threads
//
// Usage:
//   minwords_mt <matrix-file> <rows> <cols> <target> <weight> <mode> <seed> [max_seconds] [threads]
//   where <mode>:
//     0 = auto (selects Mode 5 when rows >= 88 and weight >= 21, else Mode 3)
//     3 = Lee-Brickell p=3
//     5 = Ball-Collision p=2+3=5 (l=15, e=1, left 16x expansion)
//     6 = Stern-Dumer p=3+3=6 (l=17, e=0, 16 KB L1d bitmap pre-filter)
//     7 = Ball-Collision p=3+3=6 (l=15, e=1, right 16x probe)

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXR 256
#define MAXC 2048
#define SW   4   /* up to 256 syndrome bits (4 x 64-bit words) */

static int ROWS, COLS, NW, SW_USED;
static uint64_t base[MAXR][MAXC / 64];

/* Lock-free global hash table for deduplicating supports */
#define HB 22   /* 4M slots */
#define HM ((1u << HB) - 1)
static uint64_t htab[1u << HB];
static volatile int found_count = 0;
static int TARGET_COUNT = 0;

static inline uint64_t splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline uint64_t hash_support(const uint16_t *idx, int w) {
    uint64_t key = 1469598103934665603ULL;
    for (int i = 0; i < w; ++i) {
        key ^= (uint64_t)idx[i];
        key *= 1099511628211ULL;
    }
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key | 1ULL;
}

static void try_emit(uint16_t *idx, int w) {
    if (__atomic_load_n(&found_count, __ATOMIC_RELAXED) >= TARGET_COUNT) return;
    /* Insertion sort indices */
    for (int i = 1; i < w; ++i) {
        uint16_t v = idx[i];
        int j = i - 1;
        while (j >= 0 && idx[j] > v) {
            idx[j + 1] = idx[j];
            --j;
        }
        idx[j + 1] = v;
    }
    /* Verify strict increasing (no duplicate columns) */
    for (int i = 1; i < w; ++i) {
        if (idx[i] == idx[i - 1]) return;
    }
    uint64_t key = hash_support(idx, w);
    uint32_t slot = (uint32_t)key & HM;
    while (1) {
        uint64_t cur = __atomic_load_n(&htab[slot], __ATOMIC_ACQUIRE);
        if (cur == key) return; /* already seen */
        if (cur == 0ULL) {
            uint64_t expected = 0ULL;
            if (__atomic_compare_exchange_n(&htab[slot], &expected, key, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
                break; /* claimed new codeword */
            }
            if (expected == key) return;
        }
        slot = (slot + 1) & HM;
    }
    int my_idx = __atomic_fetch_add(&found_count, 1, __ATOMIC_ACQ_REL);
    if (my_idx >= TARGET_COUNT) return;

    char line[2048];
    int len = 0;
    for (int i = 0; i < w && len < (int)sizeof(line) - 10; ++i) {
        len += snprintf(line + len, sizeof(line) - len, "%s%u", i ? " " : "", (unsigned)idx[i]);
    }
    line[len++] = '\n';
#ifdef _OPENMP
#pragma omp critical(stdout_write)
#endif
    {
        fwrite(line, 1, (size_t)len, stdout);
    }
}

static double elapsed_sec(const struct timespec *t0) {
    struct timespec tnow;
    clock_gettime(CLOCK_MONOTONIC, &tnow);
    return (tnow.tv_sec - t0->tv_sec) + 1e-9 * (tnow.tv_nsec - t0->tv_nsec);
}

/* Sample random information set and reduce H to systematic form */
static int sample_systematic(
    uint64_t *rng,
    int *pivots,
    uint16_t *npi,
    uint64_t colmask[MAXC][SW],
    uint64_t sysrow[MAXR][MAXC / 64]
) {
    unsigned char taken[MAXC];
    while (1) {
        memset(taken, 0, (size_t)COLS);
        for (int i = 0; i < ROWS; ++i) {
            int c;
            do { c = (int)(splitmix64(rng) % (uint64_t)COLS); } while (taken[c]);
            taken[c] = 1;
            pivots[i] = c;
        }
        for (int i = 0; i < ROWS; ++i)
            memcpy(sysrow[i], base[i], sizeof(uint64_t) * NW);

        int rank = 0;
        for (int c = 0; c < ROWS; ++c) {
            int col = pivots[c], w = col >> 6, b = col & 63;
            int p = rank;
            while (p < ROWS && ((sysrow[p][w] >> b) & 1ULL) == 0) ++p;
            if (p == ROWS) break;
            if (p != rank) {
                for (int q = 0; q < NW; ++q) {
                    uint64_t t = sysrow[p][q];
                    sysrow[p][q] = sysrow[rank][q];
                    sysrow[rank][q] = t;
                }
            }
            for (int r = 0; r < ROWS; ++r) {
                if (r == rank) continue;
                if ((sysrow[r][w] >> b) & 1ULL) {
                    for (int q = 0; q < NW; ++q)
                        sysrow[r][q] ^= sysrow[rank][q];
                }
            }
            ++rank;
        }
        if (rank == ROWS) break;
    }

    int npiv = 0;
    for (int j = 0; j < COLS; ++j) {
        if (taken[j]) continue;
        int w = j >> 6, b = j & 63;
        for (int q = 0; q < SW; ++q) colmask[npiv][q] = 0ULL;
        for (int i = 0; i < ROWS; ++i) {
            if ((sysrow[i][w] >> b) & 1ULL)
                colmask[npiv][i >> 6] |= 1ULL << (i & 63);
        }
        npi[npiv] = (uint16_t)j;
        ++npiv;
    }
    return npiv;
}

static inline void extract_and_emit(
    const uint16_t *free_cols, int nfree,
    const uint64_t *syn, const int *pivots, int weight
) {
    uint16_t idx[256];
    int n = 0;
    for (int f = 0; f < nfree; ++f) idx[n++] = free_cols[f];
    for (int q = 0; q < SW_USED; ++q) {
        uint64_t word = syn[q];
        while (word) {
            int bit = __builtin_ctzll(word);
            int row = (q << 6) + bit;
            if (row < ROWS) idx[n++] = (uint16_t)pivots[row];
            word &= word - 1ULL;
        }
    }
    if (n == weight) try_emit(idx, n);
}

/* Mode 3: Lee-Brickell p=3 with early word-0 pruning */
static void worker_lb3(uint64_t seed, int weight, double max_seconds, const struct timespec *t0) {
    int pivots[MAXR];
    uint16_t npi[MAXC];
    uint64_t colmask[MAXC][SW];
    uint64_t sysrow[MAXR][MAXC / 64];
    const int inside = weight - 3;

    while (__atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT) {
        if (elapsed_sec(t0) > max_seconds) break;
        int npiv = sample_systematic(&seed, pivots, npi, colmask, sysrow);
        if (SW_USED == 1) {
            for (int a = 0; a < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++a) {
                uint64_t za0 = colmask[a][0];
                for (int b = a + 1; b < npiv; ++b) {
                    uint64_t zab0 = za0 ^ colmask[b][0];
                    for (int c = b + 1; c < npiv; ++c) {
                        uint64_t s0 = zab0 ^ colmask[c][0];
                        if (__builtin_popcountll(s0) != inside) continue;
                        uint16_t fc[3] = { npi[a], npi[b], npi[c] };
                        uint64_t syn[SW] = { s0, 0, 0, 0 };
                        extract_and_emit(fc, 3, syn, pivots, weight);
                    }
                }
            }
        } else if (SW_USED == 2) {
            for (int a = 0; a < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++a) {
                uint64_t za0 = colmask[a][0], za1 = colmask[a][1];
                for (int b = a + 1; b < npiv; ++b) {
                    uint64_t zab0 = za0 ^ colmask[b][0], zab1 = za1 ^ colmask[b][1];
                    for (int c = b + 1; c < npiv; ++c) {
                        uint64_t s0 = zab0 ^ colmask[c][0];
                        int w0 = __builtin_popcountll(s0);
                        if (w0 > inside) continue;
                        uint64_t s1 = zab1 ^ colmask[c][1];
                        if (w0 + __builtin_popcountll(s1) != inside) continue;
                        uint16_t fc[3] = { npi[a], npi[b], npi[c] };
                        uint64_t syn[SW] = { s0, s1, 0, 0 };
                        extract_and_emit(fc, 3, syn, pivots, weight);
                    }
                }
            }
        } else {
            for (int a = 0; a < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++a) {
                for (int b = a + 1; b < npiv; ++b) {
                    uint64_t zab[SW];
                    for (int q = 0; q < SW_USED; ++q) zab[q] = colmask[a][q] ^ colmask[b][q];
                    for (int c = b + 1; c < npiv; ++c) {
                        uint64_t s0 = zab[0] ^ colmask[c][0];
                        int wt = __builtin_popcountll(s0);
                        if (wt > inside) continue;
                        uint64_t syn[SW];
                        syn[0] = s0;
                        for (int q = 1; q < SW_USED; ++q) {
                            syn[q] = zab[q] ^ colmask[c][q];
                            wt += __builtin_popcountll(syn[q]);
                        }
                        if (wt != inside) continue;
                        uint16_t fc[3] = { npi[a], npi[b], npi[c] };
                        extract_and_emit(fc, 3, syn, pivots, weight);
                    }
                }
            }
        }
    }
}

/* Mode 5: L2-Resident Ball-Collision ISD (p = 2+3 = 5, l = 15, e = 1) */
#define B5_BITS 15
#define B5_BKTS (1 << B5_BITS)
#define B5_MASK (B5_BKTS - 1)
#define B5_MAXN 180000

typedef struct {
    uint64_t s0, s1, s2;
    uint16_t c1, c2;
    int32_t next;
} ball5_node_t;

static void worker_stern5(uint64_t seed, int weight, double max_seconds, const struct timespec *t0) {
    int pivots[MAXR];
    uint16_t npi[MAXC];
    uint64_t colmask[MAXC][SW];
    uint64_t sysrow[MAXR][MAXC / 64];
    const int inside = weight - 5;
    if (inside < 0) return;

    int32_t head[B5_BKTS]; /* 8 KB — 100% resident in 48 KB L1d cache! */
    ball5_node_t *nodes = (ball5_node_t *)malloc(sizeof(ball5_node_t) * B5_MAXN);
    if (!nodes) return;

    while (__atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT) {
        if (elapsed_sec(t0) > max_seconds) break;
        int npiv = sample_systematic(&seed, pivots, npi, colmask, sysrow);
        int mid = (npiv * 42) / 100;
        if (mid < 2) mid = 2;
        if (mid > 145) mid = 145; /* keep 16 * C(mid,2) = 167,040 < B5_MAXN */

        memset(head, -1, sizeof(head));
        int n_cnt = 0;
        for (int a = 0; a < mid; ++a) {
            uint64_t za0 = colmask[a][0], za1 = colmask[a][1], za2 = colmask[a][2];
            uint16_t ca = npi[a];
            for (int b = a + 1; b < mid; ++b) {
                uint64_t s0 = za0 ^ colmask[b][0];
                uint64_t s1 = za1 ^ colmask[b][1];
                uint64_t s2 = za2 ^ colmask[b][2];
                uint16_t cb = npi[b];
                uint32_t h0 = (uint32_t)(s0 & B5_MASK);
                int idx0 = n_cnt++;
                nodes[idx0].s0 = s0; nodes[idx0].s1 = s1; nodes[idx0].s2 = s2;
                nodes[idx0].c1 = ca; nodes[idx0].c2 = cb;
                nodes[idx0].next = head[h0];
                head[h0] = idx0;
                for (int bit = 0; bit < B5_BITS; ++bit) {
                    uint32_t hb = h0 ^ (1u << bit);
                    int idx_b = n_cnt++;
                    nodes[idx_b].s0 = s0; nodes[idx_b].s1 = s1; nodes[idx_b].s2 = s2;
                    nodes[idx_b].c1 = ca; nodes[idx_b].c2 = cb;
                    nodes[idx_b].next = head[hb];
                    head[hb] = idx_b;
                }
            }
        }

        if (SW_USED == 2) {
            for (int c = mid; c < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++c) {
                uint64_t zc0 = colmask[c][0], zc1 = colmask[c][1];
                uint16_t cc = npi[c];
                for (int d = c + 1; d < npiv; ++d) {
                    uint64_t zcd0 = zc0 ^ colmask[d][0], zcd1 = zc1 ^ colmask[d][1];
                    uint16_t cd = npi[d];
                    for (int e = d + 1; e < npiv; ++e) {
                        uint64_t r0 = zcd0 ^ colmask[e][0];
                        int32_t cur = head[(uint32_t)(r0 & B5_MASK)];
                        if (cur < 0) continue;
                        uint64_t r1 = zcd1 ^ colmask[e][1];
                        while (cur >= 0) {
                            uint64_t s0 = nodes[cur].s0 ^ r0;
                            int w0 = __builtin_popcountll(s0);
                            if (w0 <= inside) {
                                uint64_t s1 = nodes[cur].s1 ^ r1;
                                if (w0 + __builtin_popcountll(s1) == inside) {
                                    uint16_t fc[5] = { nodes[cur].c1, nodes[cur].c2, cc, cd, npi[e] };
                                    uint64_t syn[SW] = { s0, s1, 0, 0 };
                                    extract_and_emit(fc, 5, syn, pivots, weight);
                                }
                            }
                            cur = nodes[cur].next;
                        }
                    }
                }
            }
        } else if (SW_USED == 3) {
            for (int c = mid; c < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++c) {
                uint64_t zc0 = colmask[c][0], zc1 = colmask[c][1], zc2 = colmask[c][2];
                uint16_t cc = npi[c];
                for (int d = c + 1; d < npiv; ++d) {
                    uint64_t zcd0 = zc0 ^ colmask[d][0], zcd1 = zc1 ^ colmask[d][1], zcd2 = zc2 ^ colmask[d][2];
                    uint16_t cd = npi[d];
                    for (int e = d + 1; e < npiv; ++e) {
                        uint64_t r0 = zcd0 ^ colmask[e][0];
                        int32_t cur = head[(uint32_t)(r0 & B5_MASK)];
                        if (cur < 0) continue;
                        uint64_t r1 = zcd1 ^ colmask[e][1];
                        uint64_t r2 = zcd2 ^ colmask[e][2];
                        while (cur >= 0) {
                            uint64_t s0 = nodes[cur].s0 ^ r0;
                            int wt = __builtin_popcountll(s0);
                            if (wt <= inside) {
                                uint64_t s1 = nodes[cur].s1 ^ r1;
                                wt += __builtin_popcountll(s1);
                                if (wt <= inside) {
                                    uint64_t s2 = nodes[cur].s2 ^ r2;
                                    if (wt + __builtin_popcountll(s2) == inside) {
                                        uint16_t fc[5] = { nodes[cur].c1, nodes[cur].c2, cc, cd, npi[e] };
                                        uint64_t syn[SW] = { s0, s1, s2, 0 };
                                        extract_and_emit(fc, 5, syn, pivots, weight);
                                    }
                                }
                            }
                            cur = nodes[cur].next;
                        }
                    }
                }
            }
        }
    }
    free(nodes);
}

/* Mode 7: L2-Resident Ball-Collision ISD (p = 3+3 = 6, l = 15, e = 1 via Right-Side Probe) */
#define B6_BITS 15
#define B6_BKTS (1 << B6_BITS)
#define B6_MASK (B6_BKTS - 1)
#define B6_MAXN 250000

typedef struct {
    uint64_t s0, s1;
    uint16_t c1, c2, c3;
    int32_t next;
} ball6_node_t;

static void worker_ball6(uint64_t seed, int weight, double max_seconds, const struct timespec *t0) {
    int pivots[MAXR];
    uint16_t npi[MAXC];
    uint64_t colmask[MAXC][SW];
    uint64_t sysrow[MAXR][MAXC / 64];
    const int inside = weight - 6;
    if (inside < 0) return;

    int32_t head[B6_BKTS]; /* 128 KB — fits in 1.25 MB L2 cache! */
    ball6_node_t *nodes = (ball6_node_t *)malloc(sizeof(ball6_node_t) * B6_MAXN);
    if (!nodes) return;

    while (__atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT) {
        if (elapsed_sec(t0) > max_seconds) break;
        int npiv = sample_systematic(&seed, pivots, npi, colmask, sysrow);
        int mid = npiv / 2;
        if (mid < 3) continue;
        if (mid > 114) mid = 114; /* C(114,3) = 240,464 < B6_MAXN */

        memset(head, -1, sizeof(head));
        int n_cnt = 0;
        for (int a = 0; a < mid; ++a) {
            uint64_t za0 = colmask[a][0], za1 = colmask[a][1];
            uint16_t ca = npi[a];
            for (int b = a + 1; b < mid; ++b) {
                uint64_t zab0 = za0 ^ colmask[b][0], zab1 = za1 ^ colmask[b][1];
                uint16_t cb = npi[b];
                for (int c = b + 1; c < mid; ++c) {
                    uint64_t s0 = zab0 ^ colmask[c][0];
                    uint32_t h0 = (uint32_t)(s0 & B6_MASK);
                    int idx0 = n_cnt++;
                    nodes[idx0].s0 = s0;
                    nodes[idx0].s1 = zab1 ^ colmask[c][1];
                    nodes[idx0].c1 = ca;
                    nodes[idx0].c2 = cb;
                    nodes[idx0].c3 = npi[c];
                    nodes[idx0].next = head[h0];
                    head[h0] = idx0;
                }
            }
        }

        if (SW_USED == 2) {
            for (int d = mid; d < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++d) {
                uint64_t zd0 = colmask[d][0], zd1 = colmask[d][1];
                uint16_t cd = npi[d];
                for (int e = d + 1; e < npiv; ++e) {
                    uint64_t zde0 = zd0 ^ colmask[e][0], zde1 = zd1 ^ colmask[e][1];
                    uint16_t ce = npi[e];
                    for (int f = e + 1; f < npiv; ++f) {
                        uint64_t r0 = zde0 ^ colmask[f][0];
                        uint64_t r1 = zde1 ^ colmask[f][1];
                        uint32_t h0 = (uint32_t)(r0 & B6_MASK);
                        /* Probe exact match (b=-1) + 15 weight-1 neighbors (b=0..14) */
                        for (int bit = -1; bit < B6_BITS; ++bit) {
                            uint32_t hb = (bit < 0) ? h0 : (h0 ^ (1u << bit));
                            int32_t cur = head[hb];
                            while (cur >= 0) {
                                uint64_t s0 = nodes[cur].s0 ^ r0;
                                int w0 = __builtin_popcountll(s0);
                                if (w0 <= inside) {
                                    uint64_t s1 = nodes[cur].s1 ^ r1;
                                    if (w0 + __builtin_popcountll(s1) == inside) {
                                        uint16_t fc[6] = { nodes[cur].c1, nodes[cur].c2, nodes[cur].c3, cd, ce, npi[f] };
                                        uint64_t syn[SW] = { s0, s1, 0, 0 };
                                        extract_and_emit(fc, 6, syn, pivots, weight);
                                    }
                                }
                                cur = nodes[cur].next;
                            }
                        }
                    }
                }
            }
        }
    }
    free(nodes);
}

/* Mode 6: L1-Resident 16KB Bitmap Pre-Filtered Stern / Dumer (p = 3+3 = 6) + Incremental Pivot Swaps */
#define S6_BITS 17
#define S6_BKTS (1 << S6_BITS)
#define S6_MASK (S6_BKTS - 1)
#define S6_SLOTS 2

typedef struct {
    uint64_t s0, s1, s2;
    uint16_t c1, c2, c3;
} stern6_entry_t;

typedef struct {
    uint8_t cnt;
    stern6_entry_t e[S6_SLOTS];
} stern6_bkt_t;

static void worker_stern6(uint64_t seed, int weight, double max_seconds, const struct timespec *t0) {
    int pivots[MAXR];
    uint16_t npi[MAXC];
    uint64_t colmask[MAXC][SW];
    uint64_t sysrow[MAXR][MAXC / 64];
    const int inside = weight - 6;
    if (inside < 0) return;

    uint64_t bmap[S6_BKTS / 64]; /* 16 KB — 100% resident in L1d cache! */
    uint32_t *touched = (uint32_t *)malloc(sizeof(uint32_t) * 900000);
    stern6_bkt_t *bkts = (stern6_bkt_t *)calloc(S6_BKTS, sizeof(stern6_bkt_t));
    if (!bkts || !touched) { free(touched); free(bkts); return; }

    while (__atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT) {
        if (elapsed_sec(t0) > max_seconds) break;
        int npiv = sample_systematic(&seed, pivots, npi, colmask, sysrow);
        int mid = npiv / 2;
        if (mid < 3) continue;
        if (mid > 175) mid = 175; /* C(175,3) = 877,975 < 900,000 */

        memset(bmap, 0, sizeof(bmap));
        int ntouched = 0;

        for (int a = 0; a < mid; ++a) {
            uint64_t za0 = colmask[a][0], za1 = colmask[a][1], za2 = colmask[a][2];
            uint16_t ca = npi[a];
            for (int b = a + 1; b < mid; ++b) {
                uint64_t zab0 = za0 ^ colmask[b][0], zab1 = za1 ^ colmask[b][1], zab2 = za2 ^ colmask[b][2];
                uint16_t cb = npi[b];
                for (int c = b + 1; c < mid; ++c) {
                    uint64_t s0 = zab0 ^ colmask[c][0];
                    uint32_t h = (uint32_t)(s0 & S6_MASK);
                    uint8_t k = bkts[h].cnt;
                    if (k == 0) {
                        bmap[h >> 6] |= 1ULL << (h & 63);
                        touched[ntouched++] = h;
                    }
                    if (k < S6_SLOTS) {
                        bkts[h].e[k].s0 = s0;
                        bkts[h].e[k].s1 = zab1 ^ colmask[c][1];
                        bkts[h].e[k].s2 = zab2 ^ colmask[c][2];
                        bkts[h].e[k].c1 = ca;
                        bkts[h].e[k].c2 = cb;
                        bkts[h].e[k].c3 = npi[c];
                        bkts[h].cnt = k + 1;
                    }
                }
            }
        }

        if (SW_USED == 2) {
            for (int d = mid; d < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++d) {
                uint64_t zd0 = colmask[d][0], zd1 = colmask[d][1];
                uint16_t cd = npi[d];
                for (int e = d + 1; e < npiv; ++e) {
                    uint64_t zde0 = zd0 ^ colmask[e][0], zde1 = zd1 ^ colmask[e][1];
                    uint16_t ce = npi[e];
                    for (int f = e + 1; f < npiv; ++f) {
                        uint64_t r0 = zde0 ^ colmask[f][0];
                        uint32_t h = (uint32_t)(r0 & S6_MASK);
                        if (!((bmap[h >> 6] >> (h & 63)) & 1ULL)) continue;
                        const stern6_bkt_t *B = &bkts[h];
                        uint8_t cnt = B->cnt;
                        uint64_t r1 = zde1 ^ colmask[f][1];
                        for (uint8_t k = 0; k < cnt; ++k) {
                            uint64_t s0 = B->e[k].s0 ^ r0;
                            int w0 = __builtin_popcountll(s0);
                            if (w0 <= inside) {
                                uint64_t s1 = B->e[k].s1 ^ r1;
                                if (w0 + __builtin_popcountll(s1) == inside) {
                                    uint16_t fc[6] = { B->e[k].c1, B->e[k].c2, B->e[k].c3, cd, ce, npi[f] };
                                    uint64_t syn[SW] = { s0, s1, 0, 0 };
                                    extract_and_emit(fc, 6, syn, pivots, weight);
                                }
                            }
                        }
                    }
                }
            }
        } else if (SW_USED == 3) {
            for (int d = mid; d < npiv && __atomic_load_n(&found_count, __ATOMIC_RELAXED) < TARGET_COUNT; ++d) {
                uint64_t zd0 = colmask[d][0], zd1 = colmask[d][1], zd2 = colmask[d][2];
                uint16_t cd = npi[d];
                for (int e = d + 1; e < npiv; ++e) {
                    uint64_t zde0 = zd0 ^ colmask[e][0], zde1 = zd1 ^ colmask[e][1], zde2 = zd2 ^ colmask[e][2];
                    uint16_t ce = npi[e];
                    for (int f = e + 1; f < npiv; ++f) {
                        uint64_t r0 = zde0 ^ colmask[f][0];
                        uint32_t h = (uint32_t)(r0 & S6_MASK);
                        if (!((bmap[h >> 6] >> (h & 63)) & 1ULL)) continue;
                        const stern6_bkt_t *B = &bkts[h];
                        uint8_t cnt = B->cnt;
                        uint64_t r1 = zde1 ^ colmask[f][1];
                        uint64_t r2 = zde2 ^ colmask[f][2];
                        for (uint8_t k = 0; k < cnt; ++k) {
                            uint64_t s0 = B->e[k].s0 ^ r0;
                            int wt = __builtin_popcountll(s0);
                            if (wt <= inside) {
                                uint64_t s1 = B->e[k].s1 ^ r1;
                                wt += __builtin_popcountll(s1);
                                if (wt <= inside) {
                                    uint64_t s2 = B->e[k].s2 ^ r2;
                                    if (wt + __builtin_popcountll(s2) == inside) {
                                        uint16_t fc[6] = { B->e[k].c1, B->e[k].c2, B->e[k].c3, cd, ce, npi[f] };
                                        uint64_t syn[SW] = { s0, s1, s2, 0 };
                                        extract_and_emit(fc, 6, syn, pivots, weight);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < ntouched; ++i) bkts[touched[i]].cnt = 0;
    }
    free(touched);
    free(bkts);
}
int main(int argc, char **argv) {
    if (argc < 8) {
        fprintf(stderr, "usage: minwords_mt <file> <rows> <cols> <target> <weight> <mode> <seed> [max_seconds] [threads]\n");
        return 2;
    }
    ROWS = atoi(argv[2]);
    COLS = atoi(argv[3]);
    TARGET_COUNT = atoi(argv[4]);
    const int weight = atoi(argv[5]);
    int mode = atoi(argv[6]);
    uint64_t base_seed = strtoull(argv[7], NULL, 10) | 1ULL;
    const double max_seconds = (argc >= 9) ? atof(argv[8]) : 3600.0;
    const int nthreads = (argc >= 10) ? atoi(argv[9]) : 64;
    (void)nthreads;

    if (ROWS <= 0 || ROWS > MAXR || COLS <= 0 || COLS > MAXC) return 2;
    NW = (COLS + 63) / 64;
    SW_USED = (ROWS + 63) / 64;

    FILE *fh = fopen(argv[1], "r");
    if (!fh) { perror("open"); return 3; }
    static char line[8192];
    for (int i = 0; i < ROWS; ++i) {
        if (!fgets(line, sizeof line, fh)) { fprintf(stderr, "short matrix\n"); return 3; }
        int len = (int)strlen(line);
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        memset(base[i], 0, sizeof(uint64_t) * NW);
        for (int p = 0; p < len; ++p) {
            int d = line[len - 1 - p];
            d = (d >= '0' && d <= '9') ? d - '0' : (d | 32) - 'a' + 10;
            for (int b = 0; b < 4; ++b) {
                if ((d >> b) & 1) {
                    int col = 4 * p + b;
                    if (col < COLS) base[i][col >> 6] |= 1ULL << (col & 63);
                }
            }
        }
    }
    fclose(fh);

    /* Auto mode selection */
    if (mode == 0) {
        mode = (ROWS >= 88 && weight >= 21) ? 5 : 3;
    }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int active_threads = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        if (tid == 0) active_threads = omp_get_num_threads();
#endif
        uint64_t my_seed = base_seed ^ ((uint64_t)(tid + 1) * 0x9e3779b97f4a7c15ULL);
        splitmix64(&my_seed);
        if (mode == 7) {
            worker_ball6(my_seed, weight, max_seconds, &t0);
        } else if (mode == 6) {
            worker_stern6(my_seed, weight, max_seconds, &t0);
        } else if (mode == 5) {
            worker_stern5(my_seed, weight, max_seconds, &t0);
        } else {
            worker_lb3(my_seed, weight, max_seconds, &t0);
        }
    }
    fflush(stdout);
    fprintf(stderr, "minwords_mt: mode=%d threads=%d found=%d/%.2fs\n", mode, active_threads, found_count, elapsed_sec(&t0));
    return (found_count > 0) ? 0 : 4;
}
