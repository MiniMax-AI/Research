/*
 * cycle_swap_solver.c
 *
 * Single-Stage Overlapping Cycle-Swap Collision Solver for Full Blender-256
 * (Shifted 38-Bit Clean Subspace W'[0] in {x, ~x}).
 *
 * Mechanism:
 *   1. From the standard SHA-256 IV (Word = IV_256, Sum = 0, c1 = c2 = 0),
 *      an 8-word clean entry reaches the 38-bit invariant subspace where
 *      W'[0] in {x, ~x} every round, making the 6-word stride-2 XOR chain
 *      message-independent with period 8.
 *   2. Multiple threads run a DFS sharing a deterministic trunk prefix over
 *      the first K levels. Because any vertex along a closed walk serves as
 *      a valid base point by cyclic rotation, ANY two ancestor path collisions
 *      [a1, b1) and [a2, b2) whose intervals overlap (max(a1, a2) < min(b1, b2))
 *      immediately yield two distinct closed walks U != V at the shared anchor
 *      s_m (m = max(a1, a2)) with |U| == |V| == 0 (mod 8), eliminating Stage 2.
 *   3. Swapping M1 = prefix || U || V and M2 = prefix || V || U preserves
 *      message byte length, the 13-byte self-repeating prefix, ending state s_m,
 *      modular accumulator Sum[0..7], and both message checksums (C1, C2).
 *
 * Usage:
 *   ./cycle_swap_solver [nthreads=8] [seed=0] [out_m1.bin] [out_m2.bin] [x_hex=0x0]
 */

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

typedef uint32_t u32;
typedef uint64_t u64;

static inline u32 ror32(u32 x, int n) { return (x >> n) | (x << (32 - n)); }
static inline u32 rol32(u32 x, int n) { return (x << n) | (x >> (32 - n)); }

static double wall_sec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static const u32 HIV[8] = {
    0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
    0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
};

typedef struct {
    u32 W[8];
    int c1, c2;
} St;

static inline void step32(const u32 *W, int c1, int c2, u32 f, u32 *nw, int *nc1, int *nc2) {
    u32 T1 = ror32(f, 8) ^ W[0];
    u32 tmp = ror32(W[2], 8) ^ W[4];
    u64 s = (u64)T1 + tmp + (u32)c2;
    u32 T2 = (u32)s;
    int k2 = (int)(s >> 32);

    s = (u64)(f ^ W[5]) + (rol32(W[3], 8) ^ W[1]) + (u32)c1;
    T1 = (u32)s;
    int k1 = (int)(s >> 32);

    if ((k1 ^ k2) != 0)      { T1 = rol32(T1, 7); T2 = ror32(T2, 7); }
    else if (k1 == 0)        { T1 = rol32(T1, 8); T2 = ror32(T2, 8); }
    else                     { T1 = rol32(T1, 6); T2 = ror32(T2, 6); }

    u32 t = ror32(W[0], 7);
    nw[0] = W[1] ^ T2; nw[1] = W[2] ^ T1; nw[2] = W[3] ^ T2; nw[3] = W[4] ^ T1;
    nw[4] = W[5] ^ T2; nw[5] = W[6] ^ T1; nw[6] = W[7] ^ T2; nw[7] = t ^ T1;
    *nc1 = k1;
    *nc2 = k2;
}

static u32 X_CONST = 0u;

static inline int children(const St *s, St *ns, u32 *fr) {
    int n = 0;
    const u32 tgt[2] = { X_CONST, ~X_CONST };
    for (int b = 0; b < 2; b++) {
        u32 goal = s->W[1] ^ tgt[b];
        u32 K2 = ror32(s->W[2], 8) ^ s->W[4];
        const int rs[3] = {6, 7, 8};
        for (int i = 0; i < 3; i++) {
            u32 f = rol32((rol32(goal, rs[i]) - K2 - (u32)s->c2) ^ s->W[0], 8);
            u32 nw[8];
            int k1, k2;
            step32(s->W, s->c1, s->c2, f, nw, &k1, &k2);
            if (nw[0] == tgt[b]) {
                memcpy(ns[n].W, nw, 32);
                ns[n].c1 = k1;
                ns[n].c2 = k2;
                fr[n] = f;
                n++;
                break;
            }
        }
    }
    return n;
}

#define NENTRY 8
#define LIMITW 4088

static St BASE;
static u32 ENTRY[NENTRY];

static void build_entry(void) {
    St stk_s[40];
    u32 stk_f[40];
    St kids[40][2];
    u32 kfr[40][2];
    int nch[40], idx[40];
    memcpy(stk_s[0].W, HIV, 32);
    stk_s[0].c1 = 0;
    stk_s[0].c2 = 0;
    nch[0] = -1;
    int sp = 0;
    while (sp >= 0) {
        if (sp == NENTRY + 16) {
            BASE = stk_s[NENTRY];
            memcpy(ENTRY, stk_f, NENTRY * sizeof(u32));
            return;
        }
        if (nch[sp] < 0) {
            nch[sp] = children(&stk_s[sp], kids[sp], kfr[sp]);
            idx[sp] = 0;
        }
        if (idx[sp] >= nch[sp]) { sp--; continue; }
        int p = idx[sp]++;
        stk_f[sp] = kfr[sp][p];
        stk_s[sp + 1] = kids[sp][p];
        nch[sp + 1] = -1;
        sp++;
    }
    fprintf(stderr, "entry DFS failed\n");
    exit(1);
}

static inline int seq(const St *a, const St *b) {
    return a->c1 == b->c1 && a->c2 == b->c2 && memcmp(a->W, b->W, 32) == 0;
}

static inline u64 hkey(const St *s) {
    u64 h = 1469598103934665603ULL;
    for (int i = 0; i < 8; i++) { h ^= s->W[i]; h *= 1099511628211ULL; }
    h ^= (u64)(s->c1 * 2 + s->c2);
    h *= 1099511628211ULL;
    return h;
}

typedef struct {
    St *s;
    int *val;
    unsigned char *used;
    u32 mask;
    int last;
} Set;

static void set_init(Set *S, u32 bits) {
    S->mask = (1u << bits) - 1u;
    S->s = (St *)malloc(sizeof(St) * (size_t)(S->mask + 1u));
    S->val = (int *)malloc(sizeof(int) * (size_t)(S->mask + 1u));
    S->used = (unsigned char *)calloc(S->mask + 1u, 1);
}

static void set_free(Set *S) {
    free(S->s);
    free(S->val);
    free(S->used);
}

static int set_put(Set *S, const St *s, int val, int *old) {
    u32 q = (u32)(hkey(s) & S->mask);
    for (u32 p = 0; p <= S->mask; p++) {
        u32 j = (q + p) & S->mask;
        if (!S->used[j]) {
            S->used[j] = 1;
            S->s[j] = *s;
            S->val[j] = val;
            S->last = (int)j;
            return 0;
        }
        if (seq(&S->s[j], s)) {
            *old = S->val[j];
            S->val[j] = val;
            return 1;
        }
    }
    return 0;
}

typedef struct {
    St s;
    int next, nch;
    St kid[2];
    u32 kfr[2];
} Frame;

#define MAXC 512
typedef struct { int a, b; u32 *w; } Col;
static Col COL[MAXC];
static int NCOL = 0;
static pthread_mutex_t cmu = PTHREAD_MUTEX_INITIALIZER;
static volatile int DONE = 0;
static int R_M, R_LU, R_LV;
static u32 *R_PRE = NULL, *R_U = NULL, *R_V = NULL;
static long long S1NODES = 0;

static int NT = 8, D1 = 2040, TRUNKK = 1976;
static long long NLIM = 80000000LL;
static u64 SEEDOFF = 0;

static inline u64 trunkbit(int d, u64 seed) {
    u64 z = (u64)d * 0x9E3779B97F4A7C15ULL + 0xDEADBEEF12345678ULL + seed * 0xBF58476D1CE4E5B9ULL;
    z ^= z >> 30; z *= 0xBF58476D1CE4E5B9ULL;
    z ^= z >> 27; z *= 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return z & 1ULL;
}

static void offer(int a, int b, const u32 *w) {
    pthread_mutex_lock(&cmu);
    for (int i = 0; i < NCOL && !DONE; i++) {
        int a2 = COL[i].a, b2 = COL[i].b;
        int m = (a > a2) ? a : a2;
        int hi = (b < b2) ? b : b2;
        if (m >= hi) continue;
        if (NENTRY + m + (b - a) + (b2 - a2) > LIMITW) continue;
        if (memcmp(w, COL[i].w, (size_t)m * sizeof(u32)) != 0) continue;

        /* Ensure U and V are genuinely distinct word sequences when |U| == |V| */
        int lu = b - a, lv = b2 - a2;
        u32 *u_tmp = (u32 *)malloc(sizeof(u32) * (size_t)lu);
        u32 *v_tmp = (u32 *)malloc(sizeof(u32) * (size_t)lv);
        int t = 0;
        for (int k = m; k < b; k++) u_tmp[t++] = w[k];
        for (int k = a; k < m; k++) u_tmp[t++] = w[k];
        t = 0;
        for (int k = m; k < b2; k++) v_tmp[t++] = COL[i].w[k];
        for (int k = a2; k < m; k++) v_tmp[t++] = COL[i].w[k];

        /* Check U|V != V|U */
        int same = 1;
        for (int idx = 0; idx < lu + lv; idx++) {
            u32 wu = (idx < lu) ? u_tmp[idx] : v_tmp[idx - lu];
            u32 wv = (idx < lv) ? v_tmp[idx] : u_tmp[idx - lv];
            if (wu != wv) { same = 0; break; }
        }
        if (same) { free(u_tmp); free(v_tmp); continue; }

        R_M = m; R_LU = lu; R_LV = lv;
        R_PRE = (u32 *)malloc(sizeof(u32) * (size_t)m);
        memcpy(R_PRE, w, sizeof(u32) * (size_t)m);
        R_U = u_tmp;
        R_V = v_tmp;
        DONE = 1;
        break;
    }
    if (!DONE && NCOL < MAXC) {
        COL[NCOL].a = a;
        COL[NCOL].b = b;
        COL[NCOL].w = (u32 *)malloc(sizeof(u32) * (size_t)b);
        memcpy(COL[NCOL].w, w, sizeof(u32) * (size_t)b);
        NCOL++;
    }
    pthread_mutex_unlock(&cmu);
}

static void *worker(void *arg) {
    long id = (long)arg;
    u64 rng = 0x9E3779B97F4A7C15ULL * (u64)(id + 1) + 12345ULL + SEEDOFF;
    for (int i = 0; i < 32; i++) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; }

    Set H;
    set_init(&H, 14);
    Frame *stk = (Frame *)malloc(sizeof(Frame) * (size_t)(D1 + 2));
    u32 *fr = (u32 *)malloc(sizeof(u32) * (size_t)(D1 + 2));
    int *isnew = (int *)malloc(sizeof(int) * (size_t)(D1 + 2));
    int *slot = (int *)malloc(sizeof(int) * (size_t)(D1 + 2));
    long long nodes = 0;

    int sp = 0;
    stk[0].s = BASE;
    stk[0].nch = -1;
    stk[0].next = 0;
    while (sp >= 0 && nodes < NLIM && !DONE) {
        Frame *f = &stk[sp];
        if (f->nch < 0) {
            nodes++;
            int old = -1;
            if (set_put(&H, &f->s, sp, &old)) {
                if (sp > old) offer(old, sp, fr);
                sp--;
                continue;
            } else {
                isnew[sp] = 1;
                slot[sp] = H.last;
            }
            if (sp >= D1 || nodes >= NLIM || DONE) {
                if (isnew[sp]) H.used[slot[sp]] = 0;
                sp--;
                continue;
            }
            St ns[2]; u32 f2[2];
            int n = children(&f->s, ns, f2);
            f->nch = n;
            f->next = 0;
            for (int i = 0; i < n; i++) { f->kid[i] = ns[i]; f->kfr[i] = f2[i]; }
            if (n == 2) {
                rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                int sw = (sp < TRUNKK) ? (int)trunkbit(sp, SEEDOFF) : (int)(rng & 1ULL);
                if (sw) {
                    St ts = f->kid[0]; f->kid[0] = f->kid[1]; f->kid[1] = ts;
                    u32 tf = f->kfr[0]; f->kfr[0] = f->kfr[1]; f->kfr[1] = tf;
                }
            }
            if (n == 0) {
                if (isnew[sp]) H.used[slot[sp]] = 0;
                sp--;
                continue;
            }
        }
        if (f->next >= f->nch) {
            if (isnew[sp]) H.used[slot[sp]] = 0;
            sp--;
            continue;
        }
        int pick = f->next++;
        fr[sp] = f->kfr[pick];
        stk[sp + 1].s = f->kid[pick];
        stk[sp + 1].nch = -1;
        stk[sp + 1].next = 0;
        sp++;
    }
    __sync_fetch_and_add(&S1NODES, nodes);
    set_free(&H);
    free(stk); free(fr); free(isnew); free(slot);
    return NULL;
}

static void write_words_le(const char *path, const u32 *w1, int n1, const u32 *w2, int n2, const u32 *w3, int n3, const u32 *w4, int n4) {
    FILE *fp = fopen(path, "wb");
    if (!fp) { perror("fopen"); exit(1); }
    const u32 *arrs[4] = {w1, w2, w3, w4};
    int lens[4] = {n1, n2, n3, n4};
    for (int a = 0; a < 4; a++) {
        for (int i = 0; i < lens[a]; i++) {
            u32 v = arrs[a][i];
            unsigned char le[4] = {(unsigned char)v, (unsigned char)(v >> 8), (unsigned char)(v >> 16), (unsigned char)(v >> 24)};
            fwrite(le, 1, 4, fp);
        }
    }
    fclose(fp);
}

int main(int argc, char **argv) {
    NT = (argc > 1) ? atoi(argv[1]) : 8;
    if (NT < 1) NT = 1;
    if (NT > 64) NT = 64;
    u64 base_seed = (argc > 2) ? strtoull(argv[2], NULL, 0) : 0ULL;
    const char *out_m1 = (argc > 3) ? argv[3] : NULL;
    const char *out_m2 = (argc > 4) ? argv[4] : NULL;
    if (argc > 5) X_CONST = (u32)strtoul(argv[5], NULL, 0);

    double t0 = wall_sec();
    build_entry();

    for (int attempt = 0; attempt < 8 && !DONE; attempt++) {
        SEEDOFF = base_seed + (u64)attempt * 7919ULL;
        pthread_t th[64];
        for (long i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, (void *)i);
        for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    }
    double dt = wall_sec() - t0;

    if (!DONE) {
        fprintf(stderr, "Search failed to pair overlapping cycles\n");
        return 1;
    }

    int total_words = NENTRY + R_M + R_LU + R_LV;
    printf("[Subspace Cycle-Swap] Solved in %.3f s (%d threads, %lld DFS nodes, x=0x%08x)\n",
           dt, NT, S1NODES, X_CONST);
    printf("  Anchor prefix : %d words (%d entry + %d trunk)\n", NENTRY + R_M, NENTRY, R_M);
    printf("  Closed walks  : |U| = %d words, |V| = %d words (both == 0 mod 8: %s)\n",
           R_LU, R_LV, ((R_LU % 8 == 0) && (R_LV % 8 == 0)) ? "YES" : "NO");
    printf("  Message size  : %d words (%d bytes)\n", total_words, total_words * 4);

    if (out_m1 && out_m2) {
        write_words_le(out_m1, ENTRY, NENTRY, R_PRE, R_M, R_U, R_LU, R_V, R_LV);
        write_words_le(out_m2, ENTRY, NENTRY, R_PRE, R_M, R_V, R_LV, R_U, R_LU);
    }

    /* Also save U, V specification if requested via env var */
    const char *uv_path = getenv("SAVE_UV_SPEC");
    if (uv_path) {
        FILE *fp = fopen(uv_path, "w");
        if (fp) {
            fprintf(fp, "entry");
            for (int i = 0; i < NENTRY; i++) fprintf(fp, " %08x", ENTRY[i]);
            fprintf(fp, "\nprefix %d", R_M);
            for (int i = 0; i < R_M; i++) fprintf(fp, " %08x", R_PRE[i]);
            fprintf(fp, "\nU %d", R_LU);
            for (int i = 0; i < R_LU; i++) fprintf(fp, " %08x", R_U[i]);
            fprintf(fp, "\nV %d", R_LV);
            for (int i = 0; i < R_LV; i++) fprintf(fp, " %08x", R_V[i]);
            fprintf(fp, "\n");
            fclose(fp);
        }
    }

    for (int i = 0; i < NCOL; i++) free(COL[i].w);
    free(R_PRE); free(R_U); free(R_V);
    return 0;
}
