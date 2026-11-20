/*
 * snf_herding_solver.c
 *
 * Pure-C Proof-of-Concept: Fixed-IV Partial-Digest Preimage (Herding Sum) on Full Blender-256
 * via Multi-Shifted 2-Adic SNF Steering + 38-Bit Clean Subspace (W'[0] in {a, ~a}) + Growing-Sink-DAG + Wagner k-Tree.
 *
 * Supports:
 *   Mode 32: Hits any arbitrary 32-bit digest prefix Digest[0..31] == T0* (solving the 64-bit
 *            joint modular constraint (Delta C1, Delta Sum[0]) == (0, delta_0) mod 2^32 via a
 *            2-level 4-list Wagner tree in ~1.5 seconds).
 *   Mode 64: Hits any arbitrary 64-bit digest prefix Digest[0..63] == (T0*, T1*) (TWO full 32-bit
 *            digest words simultaneously, solving the 96-bit joint modular constraint
 *            (Delta C1, Delta Sum[0], Delta Sum[1]) == (0, delta_0, delta_1) mod 2^32 via a
 *            streaming 3-level 8-list Wagner tree in ~5 seconds).
 *
 * Usage:
 *   ./snf_herding_solver 32 0xdeadbeef out_32.bin [nthreads]
 *   ./snf_herding_solver 64 0xdeadbeefcafebabe out_64.bin [nthreads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef uint16_t u16;
typedef uint8_t  u8;

static double wall_sec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static inline u32 ror32(u32 x, int n) { return (x >> n) | (x << (32 - n)); }
static inline u32 rol32(u32 x, int n) { return (x << n) | (x >> (32 - n)); }

static const u32 HIV[8] = {
    0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
    0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
};

/* ==================== Exact Full Blender-256 Reference Engine ==================== */
typedef struct {
    u32 Word[8];
    u32 Sum[8];
    u32 C1, C2;
    u8 carries;
    int roundsToGo;
    u64 srcLen;
    u64 srcLen2[3];
    u8 residue[8];
    int residueLength;
    u8 firstBlock[13];
    int firstByteCount;
} Blender256State;

static void b256_init(Blender256State *ss) {
    memset(ss, 0, sizeof(*ss));
    memcpy(ss->Word, HIV, 32);
    ss->roundsToGo = 32;
}

static inline void step32(const u32 *W, int c1, int c2, u32 f,
                          u32 *nw, int *nc1, int *nc2) {
    u32 T1 = ror32(f, 8) ^ W[0];
    u32 tmp = ror32(W[2], 8) ^ W[4];
    u64 s = (u64)T1 + tmp + (u32)c2;
    u32 T2 = (u32)s;
    int k2 = (int)(s >> 32);

    s = (u64)(f ^ W[5]) + (rol32(W[3], 8) ^ W[1]) + (u32)c1;
    T1 = (u32)s;
    int k1 = (int)(s >> 32);

    if ((k1 ^ k2) != 0) { T1 = rol32(T1, 7); T2 = ror32(T2, 7); }
    else if (k1 == 0)   { T1 = rol32(T1, 8); T2 = ror32(T2, 8); }
    else                { T1 = rol32(T1, 6); T2 = ror32(T2, 6); }

    u32 t = ror32(W[0], 7);
    nw[0] = W[1] ^ T2; nw[1] = W[2] ^ T1; nw[2] = W[3] ^ T2; nw[3] = W[4] ^ T1;
    nw[4] = W[5] ^ T2; nw[5] = W[6] ^ T1; nw[6] = W[7] ^ T2; nw[7] = t ^ T1;
    *nc1 = k1; *nc2 = k2;
}

static void b256_do_hash(Blender256State *ss, const u8 *dta, size_t bytecount) {
    if (bytecount == 0) return;
    int c1 = ss->carries & 1;
    int c2 = (ss->carries >> 1) & 1;
    size_t pos = 0;
    while (pos < bytecount) {
        u32 frag = (u32)dta[pos] | ((u32)dta[pos + 1] << 8) |
                   ((u32)dta[pos + 2] << 16) | ((u32)dta[pos + 3] << 24);
        pos += 4;
        ss->C1 += frag;
        ss->C2 += ~frag;
        u32 nw[8]; int nc1, nc2;
        step32(ss->Word, c1, c2, frag, nw, &nc1, &nc2);
        for (int i = 0; i < 8; i++) {
            ss->Word[i] = nw[i];
            ss->Sum[i] += nw[i];
        }
        c1 = nc1; c2 = nc2;
        if (--ss->roundsToGo == 0) ss->roundsToGo = 16;
    }
    ss->carries = (u8)((c2 << 1) | (c1 & 1));
}

static void b256_do_update(Blender256State *ss, const u8 *dta, size_t bytecount) {
    if (ss->firstByteCount < 13) {
        size_t need = (size_t)(13 - ss->firstByteCount);
        if (bytecount < need) need = bytecount;
        memcpy(ss->firstBlock + ss->firstByteCount, dta, need);
        ss->firstByteCount += (int)need;
    }
    b256_do_hash(ss, dta, bytecount);
}

static void b256_update_words(Blender256State *ss, const u32 *words, size_t nwords) {
    /* Word-aligned message update (exact match with byte-level little-endian update) */
    u64 bitlen = (u64)nwords * 32ULL;
    u64 old = ss->srcLen;
    ss->srcLen = old + bitlen;
    if (ss->srcLen < (old | bitlen)) ss->srcLen2[0]++;

    u8 chunk[4096];
    size_t idx = 0;
    while (idx < nwords) {
        size_t take = nwords - idx;
        if (take > 1024) take = 1024;
        for (size_t i = 0; i < take; i++) {
            u32 w = words[idx + i];
            chunk[4 * i + 0] = (u8)(w);
            chunk[4 * i + 1] = (u8)(w >> 8);
            chunk[4 * i + 2] = (u8)(w >> 16);
            chunk[4 * i + 3] = (u8)(w >> 24);
        }
        b256_do_update(ss, chunk, take * 4);
        idx += take;
    }
}

static void b256_final(Blender256State *ss, u32 out_words_be[8], u8 out_bytes[32]) {
    u8 buffer[64];
    memset(buffer, 0, sizeof(buffer));
    int lenlen = 0;
    u64 datalen;
    if (ss->srcLen2[2] != 0)      { lenlen = 24; datalen = ss->srcLen2[2]; }
    else if (ss->srcLen2[1] != 0) { lenlen = 16; datalen = ss->srcLen2[1]; }
    else if (ss->srcLen2[0] != 0) { lenlen = 8;  datalen = ss->srcLen2[0]; }
    else                          { datalen = ss->srcLen; }
    do { lenlen++; datalen >>= 8; } while (datalen != 0);

    int spaceleft = (ss->roundsToGo - 2) * 4;
    int spaceneeded = lenlen + 2 + ss->residueLength;
    while (spaceleft < spaceneeded) spaceleft += 64;
    int fillsize = spaceleft + 1 - spaceneeded;

    int i = ss->residueLength;
    if (i > 0) memcpy(buffer, ss->residue, (size_t)i);
    int j = ss->firstByteCount;
    if (j == 0) j = 13;
    int filllen = j;
    if (fillsize > j && (j << 1) <= 13) {
        int k = j;
        while (j + k <= 13) {
            memcpy(ss->firstBlock + k, ss->firstBlock, (size_t)j);
            k += j;
        }
        filllen = k;
    }
    int hwm = ((64 - filllen) / 4) * 4;
    while (fillsize > 0) {
        if (filllen > fillsize) filllen = fillsize;
        memcpy(buffer + i, ss->firstBlock, (size_t)filllen);
        fillsize -= filllen;
        i += filllen;
        if (i > hwm) {
            b256_do_update(ss, buffer, (size_t)hwm);
            i -= hwm;
            memmove(buffer, buffer + hwm, (size_t)i);
        }
    }
    datalen = ss->srcLen;
    j = lenlen;
    int k = 0;
    while (1) {
        buffer[i++] = (u8)(datalen & 0xFFu);
        if (--j == 0) break;
        datalen >>= 8;
        if ((k & 7) == 7) datalen = ss->srcLen2[k >> 3];
        k++;
        if (i > hwm && (i + j + 1) > 64) {
            b256_do_update(ss, buffer, (size_t)hwm);
            i -= hwm;
            memmove(buffer, buffer + hwm, (size_t)i);
        }
    }
    buffer[i++] = (u8)lenlen;
    b256_do_update(ss, buffer, (size_t)i);

    u32 c1_bar = ~ss->C1;
    u32 c2_val = ss->C2;
    for (int b = 0; b < 4; b++) buffer[b]     = (u8)(c1_bar >> (8 * b));
    for (int b = 0; b < 4; b++) buffer[4 + b] = (u8)(c2_val >> (8 * b));
    b256_do_update(ss, buffer, 8);

    for (int w = 0; w < 8; w++) {
        out_words_be[w] = ss->Sum[w];
        out_bytes[4 * w + 0] = (u8)(ss->Sum[w] >> 24);
        out_bytes[4 * w + 1] = (u8)(ss->Sum[w] >> 16);
        out_bytes[4 * w + 2] = (u8)(ss->Sum[w] >> 8);
        out_bytes[4 * w + 3] = (u8)(ss->Sum[w]);
    }
}

/* ==================== Shifted 38-Bit Clean Subspace + Sink-DAG ==================== */
typedef struct { u32 W[8]; int c1, c2; } St;

static inline int st_eq(const St *a, const St *b) {
    return a->c1 == b->c1 && a->c2 == b->c2 && memcmp(a->W, b->W, 32) == 0;
}
static inline u64 st_hash(const St *s) {
    u64 h = 1469598103934665603ULL;
    for (int i = 0; i < 8; i++) { h ^= s->W[i]; h *= 1099511628211ULL; }
    h ^= (u64)(s->c1 * 2 + s->c2);
    h *= 1099511628211ULL;
    return h;
}
static inline int children_a(const St *s, u32 a_const, St *ns, u32 *fr) {
    int n = 0;
    const u32 tgt[2] = { a_const, ~a_const };
    for (int b = 0; b < 2; b++) {
        u32 goal = s->W[1] ^ tgt[b];
        u32 K2 = ror32(s->W[2], 8) ^ s->W[4];
        const int rs[3] = {6, 7, 8};
        for (int i = 0; i < 3; i++) {
            u32 f = rol32((rol32(goal, rs[i]) - K2 - (u32)s->c2) ^ s->W[0], 8);
            u32 nw[8]; int k1, k2;
            step32(s->W, s->c1, s->c2, f, nw, &k1, &k2);
            if (nw[0] == tgt[b]) {
                memcpy(ns[n].W, nw, 32);
                ns[n].c1 = k1; ns[n].c2 = k2;
                fr[n] = f; n++; break;
            }
        }
    }
    return n;
}

static int build_shifted_entry_try(const St *start, u32 a_const, St *out_state, u32 *entry_frags, int nsteps) {
    const int lookahead = nsteps + 24;
    St stk_s[64]; u32 stk_f[64];
    St kids[64][2]; u32 kfr[64][2]; int nch[64], idx[64];
    stk_s[0] = *start; nch[0] = -1;
    int sp = 0;
    while (sp >= 0) {
        if (sp == lookahead) {
            *out_state = stk_s[nsteps];
            memcpy(entry_frags, stk_f, sizeof(u32) * (size_t)nsteps);
            return 1;
        }
        if (nch[sp] < 0) {
            nch[sp] = children_a(&stk_s[sp], a_const, kids[sp], kfr[sp]);
            idx[sp] = 0;
        }
        if (idx[sp] >= nch[sp]) { sp--; continue; }
        int p = idx[sp]++;
        stk_f[sp] = kfr[sp][p];
        stk_s[sp + 1] = kids[sp][p];
        nch[sp + 1] = -1;
        sp++;
    }
    return 0;
}

typedef struct {
    int len;
    u32 c1, sum[8];
    u32 *frags;
} WalkInfo;

static void eval_walk(const St *B, const u32 *frags, int len, WalkInfo *out) {
    St s = *B;
    u32 c1 = 0, sum[8] = {0};
    for (int i = 0; i < len; i++) {
        u32 nw[8]; int k1, k2;
        c1 += frags[i];
        step32(s.W, s.c1, s.c2, frags[i], nw, &k1, &k2);
        for (int j = 0; j < 8; j++) sum[j] += nw[j];
        memcpy(s.W, nw, 32); s.c1 = k1; s.c2 = k2;
    }
    if (!st_eq(&s, B)) { fprintf(stderr, "Walk failed to close!\n"); exit(1); }
    out->len = len; out->c1 = c1; memcpy(out->sum, sum, 32);
    out->frags = malloc(sizeof(u32) * (size_t)len);
    memcpy(out->frags, frags, sizeof(u32) * (size_t)len);
}

typedef struct { St s; int dist_to_B, next_idx; u32 edge_frag; } DagNode;
typedef struct { u32 *slot_idx; u32 mask; DagNode *nodes; int n_nodes, cap_nodes; } SinkDag;

static void dag_init(SinkDag *D, int hbits, int cap) {
    D->mask = (1u << hbits) - 1u;
    D->slot_idx = calloc((size_t)(D->mask + 1u), sizeof(u32));
    D->nodes = malloc(sizeof(DagNode) * (size_t)cap);
    D->n_nodes = 0; D->cap_nodes = cap;
}
static void dag_free(SinkDag *D) { free(D->slot_idx); free(D->nodes); }
static int dag_lookup(const SinkDag *D, const St *s) {
    u32 q = (u32)(st_hash(s) & D->mask);
    for (u32 p = 0; p <= D->mask; p++) {
        u32 j = (q + p) & D->mask;
        u32 v = __atomic_load_n(&D->slot_idx[j], __ATOMIC_ACQUIRE);
        if (v == 0) return -1;
        if (st_eq(&D->nodes[v - 1].s, s)) return (int)(v - 1);
    }
    return -1;
}
static int dag_insert(SinkDag *D, const St *s, int dist, int nxt, u32 frag) {
    if (D->n_nodes >= D->cap_nodes) return -1;
    u32 q = (u32)(st_hash(s) & D->mask);
    for (u32 p = 0; p <= D->mask; p++) {
        u32 j = (q + p) & D->mask;
        u32 v = __atomic_load_n(&D->slot_idx[j], __ATOMIC_ACQUIRE);
        if (v == 0) {
            int idx = D->n_nodes++;
            D->nodes[idx].s = *s; D->nodes[idx].dist_to_B = dist;
            D->nodes[idx].next_idx = nxt; D->nodes[idx].edge_frag = frag;
            __atomic_store_n(&D->slot_idx[j], (u32)(idx + 1), __ATOMIC_RELEASE);
            return idx;
        }
        if (st_eq(&D->nodes[v - 1].s, s)) return (int)(v - 1);
    }
    return -1;
}

typedef struct { St s; int next, nch; St kid[2]; u32 kfr[2]; } DfsFrame;

typedef struct {
    St start; u32 a_const; int d_cap;
    volatile int *done; pthread_mutex_t *mu;
    int *out_d1, *out_L; u32 **out_frags; u64 seed;
} S1Arg;

static void *s1_worker(void *p) {
    S1Arg *a = (S1Arg *)p;
    u64 rng = a->seed;
    #define RND() (rng ^= rng << 13, rng ^= rng >> 7, rng ^= rng << 17, rng)
    for (int i = 0; i < 32; i++) RND();

    const u32 hmask = (1u << 16) - 1u;
    int *head = malloc(sizeof(int) * (size_t)(hmask + 1u));
    for (u32 i = 0; i <= hmask; i++) head[i] = -1;
    int D = a->d_cap;
    DfsFrame *stk = malloc(sizeof(DfsFrame) * (size_t)(D + 2));
    u32 *fr = malloc(sizeof(u32) * (size_t)(D + 2));
    int *prev = malloc(sizeof(int) * (size_t)(D + 2));
    u32 *hb = malloc(sizeof(u32) * (size_t)(D + 2));
    int *pushed = calloc((size_t)(D + 2), sizeof(int));

    stk[0].s = a->start; stk[0].nch = -1;
    int sp = 0;
    while (sp >= 0 && !*(a->done)) {
        DfsFrame *f = &stk[sp];
        if (f->nch < 0) {
            u32 h = (u32)(st_hash(&f->s) & hmask);
            int hit = -1;
            for (int cur = head[h]; cur >= 0; cur = prev[cur])
                if (st_eq(&stk[cur].s, &f->s)) { hit = cur; break; }
            if (hit >= 0) {
                pushed[sp] = 0;
                int L = sp - hit;
                if (L > 0 && (L & 7) == 0 && L <= 2400) {
                    pthread_mutex_lock(a->mu);
                    if (!*(a->done)) {
                        *(a->done) = 1; *(a->out_d1) = hit; *(a->out_L) = L;
                        *(a->out_frags) = malloc(sizeof(u32) * (size_t)sp);
                        memcpy(*(a->out_frags), fr, sizeof(u32) * (size_t)sp);
                    }
                    pthread_mutex_unlock(a->mu);
                }
                sp--; continue;
            }
            hb[sp] = h; prev[sp] = head[h]; head[h] = sp; pushed[sp] = 1;
            if (sp >= D || *(a->done)) { head[hb[sp]] = prev[sp]; pushed[sp] = 0; sp--; continue; }
            f->nch = children_a(&f->s, a->a_const, f->kid, f->kfr);
            f->next = 0;
            if (f->nch == 2 && (RND() & 1)) {
                St ts = f->kid[0]; f->kid[0] = f->kid[1]; f->kid[1] = ts;
                u32 tf = f->kfr[0]; f->kfr[0] = f->kfr[1]; f->kfr[1] = tf;
            }
            if (f->nch == 0) { head[hb[sp]] = prev[sp]; pushed[sp] = 0; sp--; continue; }
        }
        if (f->next >= f->nch) {
            if (pushed[sp]) { head[hb[sp]] = prev[sp]; pushed[sp] = 0; }
            sp--; continue;
        }
        int pick = f->next++;
        fr[sp] = f->kfr[pick];
        stk[sp + 1].s = f->kid[pick]; stk[sp + 1].nch = -1;
        sp++;
    }
    #undef RND
    free(head); free(stk); free(fr); free(prev); free(hb); free(pushed);
    return NULL;
}

typedef struct {
    St B; u32 a_const;
    SinkDag *dag; pthread_mutex_t *dag_mu;
    WalkInfo *walks; volatile int *n_walks; int target_walks;
    pthread_mutex_t *w_mu; u64 seed;
} S2Arg;

static void *s2_worker(void *p) {
    S2Arg *a = (S2Arg *)p;
    u64 rng = a->seed;
    #define RND() (rng ^= rng << 13, rng ^= rng >> 7, rng ^= rng << 17, rng)
    for (int i = 0; i < 32; i++) RND();

    const int D2 = 512, MIN_SP = 24, MAX_DIST = 2800;
    DfsFrame *stk = malloc(sizeof(DfsFrame) * (size_t)(D2 + 4));
    u32 *fr = malloc(sizeof(u32) * (size_t)(D2 + 4));
    u32 *tw = malloc(sizeof(u32) * (size_t)(D2 + MAX_DIST + 16));

    while (__atomic_load_n(a->n_walks, __ATOMIC_RELAXED) < a->target_walks) {
        stk[0].s = a->B; stk[0].nch = -1;
        int sp = 0, local = 0;
        while (sp >= 0 && local < 65536 &&
               __atomic_load_n(a->n_walks, __ATOMIC_RELAXED) < a->target_walks) {
            DfsFrame *f = &stk[sp];
            if (f->nch < 0) {
                local++;
                if (sp >= MIN_SP) {
                    int hit = dag_lookup(a->dag, &f->s);
                    int tlen = (hit >= 0) ? a->dag->nodes[hit].dist_to_B : -1;
                    if (hit >= 0 && ((sp + tlen) & 7) == 0) {
                        int tot = sp + tlen;
                        memcpy(tw, fr, sizeof(u32) * (size_t)sp);
                        int cur = hit, pos = sp;
                        while (cur > 0) { tw[pos++] = a->dag->nodes[cur].edge_frag; cur = a->dag->nodes[cur].next_idx; }
                        pthread_mutex_lock(a->dag_mu);
                        int cidx = hit, dc = tlen;
                        for (int t = sp - 1; t >= MIN_SP; t--) {
                            if (++dc > MAX_DIST) break;
                            int ni = dag_insert(a->dag, &stk[t].s, dc, cidx, fr[t]);
                            if (ni < 0) break;
                            cidx = ni;
                        }
                        pthread_mutex_unlock(a->dag_mu);

                        if (pos == tot) {
                            WalkInfo wi;
                            eval_walk(&a->B, tw, tot, &wi);
                            pthread_mutex_lock(a->w_mu);
                            int nw = *(a->n_walks), dup = 0;
                            for (int k = 0; k < nw; k++) {
                                if (a->walks[k].len == wi.len && a->walks[k].c1 == wi.c1 &&
                                    a->walks[k].sum[0] == wi.sum[0]) { dup = 1; break; }
                            }
                            if (!dup && nw < a->target_walks) {
                                a->walks[nw] = wi;
                                __atomic_store_n(a->n_walks, nw + 1, __ATOMIC_RELEASE);
                            } else {
                                free(wi.frags);
                            }
                            pthread_mutex_unlock(a->w_mu);
                        }
                        sp--; continue;
                    }
                }
                if (sp >= D2) { sp--; continue; }
                f->nch = children_a(&f->s, a->a_const, f->kid, f->kfr);
                f->next = 0;
                if (f->nch == 2 && (RND() & 1)) {
                    St ts = f->kid[0]; f->kid[0] = f->kid[1]; f->kid[1] = ts;
                    u32 tf = f->kfr[0]; f->kfr[0] = f->kfr[1]; f->kfr[1] = tf;
                }
                if (f->nch == 0) { sp--; continue; }
            }
            if (f->next >= f->nch) { sp--; continue; }
            int pick = f->next++;
            fr[sp] = f->kfr[pick];
            stk[sp + 1].s = f->kid[pick]; stk[sp + 1].nch = -1;
            sp++;
        }
    }
    #undef RND
    free(stk); free(fr); free(tw);
    return NULL;
}

/* ==================== Binary Sub-Slot Representation ==================== */
typedef struct {
    int anchor;
    int u_idx, v_idx;
    u32 dc1, ds0, ds1;
} BinSlot;

/* ==================== Mode 32: 4-List 2-Level Wagner Solver (64 bits: C1, Sum[0]) ==================== */
typedef struct { u32 c1, s0, mask22; } L32Item;
typedef struct { u32 s0; u16 c1_hi11; u16 m_hi; u32 m_lo; } M32Item; /* 44-bit mask in m_lo(32) + m_hi(12) */

static int solve_wagner_32bit(const BinSlot *slots, u32 off_s0, u8 out_bits[88]) {
    /* 4 lists of 22 binary slots each: L_0, L_1, L_2, L_3 (size 2^22 = 4,194,304) */
    const u32 N = 1u << 22;
    const u32 K1_BITS = 21;
    const u32 K1_BUCKETS = 1u << K1_BITS;
    const u32 K1_MASK = K1_BUCKETS - 1u;
    const u32 CAP_M = 9000000u;

    /* Helper to precompute 2^11 + 2^11 tables for list r (slots[22*r .. 22*r+21]) */
    u32 tlo_c1[4][2048], tlo_s0[4][2048], thi_c1[4][2048], thi_s0[4][2048];
    for (int r = 0; r < 4; r++) {
        for (u32 m = 0; m < 2048u; m++) {
            u32 c1a = 0, s0a = (r == 3) ? off_s0 : 0u;
            u32 c1b = 0, s0b = 0;
            for (int b = 0; b < 11; b++) {
                if ((m >> b) & 1u) {
                    c1a += slots[22 * r + b].dc1;      s0a += slots[22 * r + b].ds0;
                    c1b += slots[22 * r + 11 + b].dc1; s0b += slots[22 * r + 11 + b].ds0;
                }
            }
            tlo_c1[r][m] = c1a; tlo_s0[r][m] = s0a;
            thi_c1[r][m] = c1b; thi_s0[r][m] = s0b;
        }
    }

    /* Merge L_0 + L_1 -> M_0 */
    u32 *cnt = calloc((size_t)(K1_BUCKETS + 1u), sizeof(u32));
    L32Item *A = malloc(sizeof(L32Item) * (size_t)N);
    for (u32 i = 0; i < N; i++) {
        u32 c1 = tlo_c1[0][i & 2047u] + thi_c1[0][i >> 11];
        cnt[(c1 & K1_MASK) + 1u]++;
    }
    for (u32 b = 1; b <= K1_BUCKETS; b++) cnt[b] += cnt[b - 1];
    u32 *cur = malloc(sizeof(u32) * (size_t)K1_BUCKETS);
    memcpy(cur, cnt, sizeof(u32) * (size_t)K1_BUCKETS);
    for (u32 i = 0; i < N; i++) {
        u32 c1 = tlo_c1[0][i & 2047u] + thi_c1[0][i >> 11];
        u32 s0 = tlo_s0[0][i & 2047u] + thi_s0[0][i >> 11];
        u32 pos = cur[c1 & K1_MASK]++;
        A[pos] = (L32Item){ c1, s0, i };
    }
    free(cur);

    M32Item *M0 = malloc(sizeof(M32Item) * (size_t)CAP_M);
    u32 nM0 = 0;
    for (u32 j = 0; j < N && nM0 < CAP_M; j++) {
        u32 yc1 = tlo_c1[1][j & 2047u] + thi_c1[1][j >> 11];
        u32 ys0 = tlo_s0[1][j & 2047u] + thi_s0[1][j >> 11];
        u32 want = (-yc1) & K1_MASK;
        for (u32 p = cnt[want]; p < cnt[want + 1u] && nM0 < CAP_M; p++) {
            u32 sc1 = A[p].c1 + yc1;
            u32 ss0 = A[p].s0 + ys0;
            u64 m44 = ((u64)j << 22) | (u64)A[p].mask22;
            M0[nM0++] = (M32Item){ ss0, (u16)(sc1 >> 21), (u16)(m44 >> 32), (u32)m44 };
        }
    }
    free(A); free(cnt);

    /* Index M_0 by 24-bit hash of (s0 & 0xFFFFFF) */
    const u32 H2_BUCKETS = 1u << 24;
    const u32 H2_MASK = H2_BUCKETS - 1u;
    u32 *cnt2 = calloc((size_t)(H2_BUCKETS + 1u), sizeof(u32));
    for (u32 i = 0; i < nM0; i++) cnt2[(M0[i].s0 & H2_MASK) + 1u]++;
    for (u32 b = 1; b <= H2_BUCKETS; b++) cnt2[b] += cnt2[b - 1];
    u32 *cur2 = malloc(sizeof(u32) * (size_t)H2_BUCKETS);
    memcpy(cur2, cnt2, sizeof(u32) * (size_t)H2_BUCKETS);
    M32Item *M0s = malloc(sizeof(M32Item) * (size_t)nM0);
    for (u32 i = 0; i < nM0; i++) M0s[cur2[M0[i].s0 & H2_MASK]++] = M0[i];
    free(cur2); free(M0);

    /* Sort L_2 and stream L_3 to probe M0s directly! */
    cnt = calloc((size_t)(K1_BUCKETS + 1u), sizeof(u32));
    A = malloc(sizeof(L32Item) * (size_t)N);
    for (u32 i = 0; i < N; i++) {
        u32 c1 = tlo_c1[2][i & 2047u] + thi_c1[2][i >> 11];
        cnt[(c1 & K1_MASK) + 1u]++;
    }
    for (u32 b = 1; b <= K1_BUCKETS; b++) cnt[b] += cnt[b - 1];
    cur = malloc(sizeof(u32) * (size_t)K1_BUCKETS);
    memcpy(cur, cnt, sizeof(u32) * (size_t)K1_BUCKETS);
    for (u32 i = 0; i < N; i++) {
        u32 c1 = tlo_c1[2][i & 2047u] + thi_c1[2][i >> 11];
        u32 s0 = tlo_s0[2][i & 2047u] + thi_s0[2][i >> 11];
        A[cur[c1 & K1_MASK]++] = (L32Item){ c1, s0, i };
    }
    free(cur);

    int solved = 0;
    for (u32 j = 0; j < N && !solved; j++) {
        u32 yc1 = tlo_c1[3][j & 2047u] + thi_c1[3][j >> 11];
        u32 ys0 = tlo_s0[3][j & 2047u] + thi_s0[3][j >> 11];
        u32 want1 = (-yc1) & K1_MASK;
        for (u32 p = cnt[want1]; p < cnt[want1 + 1u] && !solved; p++) {
            u32 rc1_hi11 = (A[p].c1 + yc1) >> 21;
            u32 rs0 = A[p].s0 + ys0;
            u32 want_s0 = -rs0;
            u16 want_c1hi = (u16)((0u - rc1_hi11) & 0x7FFu);
            u32 bk = want_s0 & H2_MASK;
            for (u32 q = cnt2[bk]; q < cnt2[bk + 1u]; q++) {
                if (M0s[q].s0 == want_s0 && M0s[q].c1_hi11 == want_c1hi) {
                    u64 left44  = ((u64)M0s[q].m_hi << 32) | (u64)M0s[q].m_lo;
                    u64 right44 = ((u64)j << 22) | (u64)A[p].mask22;
                    for (int b = 0; b < 44; b++) out_bits[b]      = (u8)((left44 >> b) & 1u);
                    for (int b = 0; b < 44; b++) out_bits[44 + b] = (u8)((right44 >> b) & 1u);
                    solved = 1;
                    break;
                }
            }
        }
    }
    free(A); free(cnt); free(M0s); free(cnt2);
    return solved;
}

/* ==================== Mode 64: Streaming 8-List 3-Level Wagner Solver (96 bits) ==================== */
typedef struct { u32 c1, s0, s1, mask24; } L64Item;
typedef struct { u32 s0, s1; u64 pack; } M64Item; /* pack: mask48 in [0..47], c1_hi9 in [48..56] */
/* Level-2 entry (20 bytes): 32-bit state word s1, 96-bit index mask (m0, m1, m2), and upper 16 bits of s0. */
typedef struct {
    u32 s1;
    u32 m0, m1, m2; /* 96 bits: m0 = bits 0..31, m1 = bits 32..63, m2 = bits 64..95 */
    u16 s0_hi16;
    u16 pad;
} N64Full;

typedef struct {
    int tid, nthreads;
    const u32 (*tlo_c1)[4096], (*tlo_s0)[4096], (*tlo_s1)[4096];
    const u32 (*thi_c1)[4096], (*thi_s0)[4096], (*thi_s1)[4096];
    int listB_id;
    const L64Item *Asorted;
    const u32 *Acnt;
    M64Item *outM;
    volatile u32 *n_outM;
    u32 cap_outM;
} StreamL1Arg;

static void *stream_l1_worker(void *p) {
    StreamL1Arg *a = (StreamL1Arg *)p;
    const u32 N = 1u << 24;
    const u32 K1_MASK = (1u << 23) - 1u;
    int r = a->listB_id;
    u32 chunk = (N + (u32)a->nthreads - 1u) / (u32)a->nthreads;
    u32 beg = (u32)a->tid * chunk;
    u32 end = beg + chunk; if (end > N) end = N;

    M64Item local_buf[256];
    int nb = 0;
    for (u32 j = beg; j < end; j++) {
        u32 lo = j & 4095u, hi = j >> 12;
        u32 yc1 = a->tlo_c1[r][lo] + a->thi_c1[r][hi];
        u32 ys0 = a->tlo_s0[r][lo] + a->thi_s0[r][hi];
        u32 ys1 = a->tlo_s1[r][lo] + a->thi_s1[r][hi];
        u32 want = (-yc1) & K1_MASK;
        for (u32 q = a->Acnt[want]; q < a->Acnt[want + 1u]; q++) {
            u32 sc1 = a->Asorted[q].c1 + yc1;
            u32 ss0 = a->Asorted[q].s0 + ys0;
            u32 ss1 = a->Asorted[q].s1 + ys1;
            u64 m48 = ((u64)j << 24) | (u64)a->Asorted[q].mask24;
            u64 pack = m48 | ((u64)((sc1 >> 23) & 0x1FFu) << 48);
            local_buf[nb++] = (M64Item){ ss0, ss1, pack };
            if (nb == 256) {
                u32 dst = __atomic_fetch_add(a->n_outM, 256u, __ATOMIC_RELAXED);
                if (dst + 256u <= a->cap_outM) memcpy(a->outM + dst, local_buf, sizeof(local_buf));
                nb = 0;
            }
        }
    }
    if (nb > 0) {
        u32 dst = __atomic_fetch_add(a->n_outM, (u32)nb, __ATOMIC_RELAXED);
        if (dst + (u32)nb <= a->cap_outM) memcpy(a->outM + dst, local_buf, sizeof(M64Item) * (size_t)nb);
    }
    return NULL;
}

/* Stream L_3 against L_2 and immediately probe sorted M_0 to emit N_0 */
typedef struct {
    int tid, nthreads;
    const u32 (*tlo_c1)[4096], (*tlo_s0)[4096], (*tlo_s1)[4096];
    const u32 (*thi_c1)[4096], (*thi_s0)[4096], (*thi_s1)[4096];
    const L64Item *L2s; const u32 *L2cnt;
    const M64Item *M0s; const u32 *M0cnt;
    N64Full *N0; volatile u32 *nN0; u32 capN0;
} StreamN0Arg;

static void *stream_n0_worker(void *p) {
    StreamN0Arg *a = (StreamN0Arg *)p;
    const u32 N = 1u << 24;
    const u32 K1_MASK = (1u << 23) - 1u;
    u32 chunk = (N + (u32)a->nthreads - 1u) / (u32)a->nthreads;
    u32 beg = (u32)a->tid * chunk, end = beg + chunk; if (end > N) end = N;

    N64Full lbuf[128]; int nb = 0;
    for (u32 j = beg; j < end; j++) {
        u32 lo = j & 4095u, hi = j >> 12;
        u32 yc1 = a->tlo_c1[3][lo] + a->thi_c1[3][hi];
        u32 ys0 = a->tlo_s0[3][lo] + a->thi_s0[3][hi];
        u32 ys1 = a->tlo_s1[3][lo] + a->thi_s1[3][hi];
        u32 w1 = (-yc1) & K1_MASK;
        for (u32 q = a->L2cnt[w1]; q < a->L2cnt[w1 + 1u]; q++) {
            u32 rc1_hi9 = ((a->L2s[q].c1 + yc1) >> 23) & 0x1FFu;
            u32 rs0 = a->L2s[q].s0 + ys0;
            u32 rs1 = a->L2s[q].s1 + ys1;
            u64 m1_48 = ((u64)j << 24) | (u64)a->L2s[q].mask24;
            u32 w2 = (((0u - rc1_hi9) & 0x1FFu) << 16) | ((-rs0) & 0xFFFFu);
            for (u32 u = a->M0cnt[w2]; u < a->M0cnt[w2 + 1u]; u++) {
                u32 tot_s0 = a->M0s[u].s0 + rs0;
                u32 tot_s1 = a->M0s[u].s1 + rs1;
                u64 m0_48 = a->M0s[u].pack & 0xFFFFFFFFFFFFULL;
                N64Full item;
                item.s1 = tot_s1;
                item.s0_hi16 = (u16)(tot_s0 >> 16);
                item.pad = 0;
                item.m0 = (u32)(m0_48 & 0xFFFFFFFFULL);
                item.m1 = (u32)((m0_48 >> 32) | ((m1_48 & 0xFFFFULL) << 16));
                item.m2 = (u32)(m1_48 >> 16);
                lbuf[nb++] = item;
                if (nb == 128) {
                    u32 dst = __atomic_fetch_add(a->nN0, 128u, __ATOMIC_RELAXED);
                    if (dst + 128u <= a->capN0) memcpy(a->N0 + dst, lbuf, sizeof(lbuf));
                    nb = 0;
                }
            }
        }
    }
    if (nb > 0) {
        u32 dst = __atomic_fetch_add(a->nN0, (u32)nb, __ATOMIC_RELAXED);
        if (dst + (u32)nb <= a->capN0) memcpy(a->N0 + dst, lbuf, sizeof(N64Full) * (size_t)nb);
    }
    return NULL;
}

/* Final Step: Stream L_7 against L_6 -> M_3 -> probe M_2 -> probe N_0 for exact 96-bit zero! */
typedef struct {
    int tid, nthreads;
    const u32 (*tlo_c1)[4096], (*tlo_s0)[4096], (*tlo_s1)[4096];
    const u32 (*thi_c1)[4096], (*thi_s0)[4096], (*thi_s1)[4096];
    const L64Item *L6s; const u32 *L6cnt;
    const M64Item *M2s; const u32 *M2cnt;
    const N64Full *N0s; const u32 *N0cnt;
    volatile int *solved;
    pthread_mutex_t *sol_mu;
    u32 *shared_out_w6; /* 192 bits */
} StreamFinalArg;

static void *stream_final_worker(void *p) {
    StreamFinalArg *a = (StreamFinalArg *)p;
    const u32 N = 1u << 24;
    const u32 K1_MASK = (1u << 23) - 1u;
    const u32 H3_MASK = (1u << 25) - 1u;
    u32 chunk = (N + (u32)a->nthreads - 1u) / (u32)a->nthreads;
    u32 beg = (u32)a->tid * chunk, end = beg + chunk; if (end > N) end = N;

    for (u32 j = beg; j < end && !__atomic_load_n(a->solved, __ATOMIC_RELAXED); j++) {
        u32 lo = j & 4095u, hi = j >> 12;
        u32 yc1 = a->tlo_c1[7][lo] + a->thi_c1[7][hi];
        u32 ys0 = a->tlo_s0[7][lo] + a->thi_s0[7][hi];
        u32 ys1 = a->tlo_s1[7][lo] + a->thi_s1[7][hi];
        u32 w1 = (-yc1) & K1_MASK;
        for (u32 q = a->L6cnt[w1]; q < a->L6cnt[w1 + 1u]; q++) {
            u32 rc1_hi9 = ((a->L6s[q].c1 + yc1) >> 23) & 0x1FFu;
            u32 rs0 = a->L6s[q].s0 + ys0;
            u32 rs1 = a->L6s[q].s1 + ys1;
            u64 m3_48 = ((u64)j << 24) | (u64)a->L6s[q].mask24;
            u32 w2 = (((0u - rc1_hi9) & 0x1FFu) << 16) | ((-rs0) & 0xFFFFu);
            for (u32 u = a->M2cnt[w2]; u < a->M2cnt[w2 + 1u]; u++) {
                u32 zs0 = a->M2s[u].s0 + rs0;
                u32 zs1 = a->M2s[u].s1 + rs1;
                u16 want_s0hi = (u16)(0u - (zs0 >> 16));
                u32 want_s1   = 0u - zs1;
                u32 bk3 = want_s1 & H3_MASK;
                for (u32 v = a->N0cnt[bk3]; v < a->N0cnt[bk3 + 1u]; v++) {
                    if (a->N0s[v].s1 == want_s1 && a->N0s[v].s0_hi16 == want_s0hi) {
                        pthread_mutex_lock(a->sol_mu);
                        if (!*(a->solved)) {
                            u64 m2_48 = a->M2s[u].pack & 0xFFFFFFFFFFFFULL;
                            a->shared_out_w6[0] = a->N0s[v].m0;
                            a->shared_out_w6[1] = a->N0s[v].m1;
                            a->shared_out_w6[2] = a->N0s[v].m2;
                            a->shared_out_w6[3] = (u32)(m2_48 & 0xFFFFFFFFULL);
                            a->shared_out_w6[4] = (u32)((m2_48 >> 32) | ((m3_48 & 0xFFFFULL) << 16));
                            a->shared_out_w6[5] = (u32)(m3_48 >> 16);
                            __atomic_store_n(a->solved, 1, __ATOMIC_RELEASE);
                        }
                        pthread_mutex_unlock(a->sol_mu);
                        return NULL;
                    }
                }
            }
        }
    }
    return NULL;
}

static void build_sorted_leaf64(int r, const u32 (*tlo_c1)[4096], const u32 (*tlo_s0)[4096], const u32 (*tlo_s1)[4096],
                                const u32 (*thi_c1)[4096], const u32 (*thi_s0)[4096], const u32 (*thi_s1)[4096],
                                L64Item **out_s, u32 **out_cnt) {
    const u32 N = 1u << 24;
    const u32 KB = 1u << 23;
    const u32 KM = KB - 1u;
    u32 *cnt = calloc((size_t)(KB + 1u), sizeof(u32));
    for (u32 i = 0; i < N; i++) {
        u32 c1 = tlo_c1[r][i & 4095u] + thi_c1[r][i >> 12];
        cnt[(c1 & KM) + 1u]++;
    }
    for (u32 b = 1; b <= KB; b++) cnt[b] += cnt[b - 1];
    u32 *cur = malloc(sizeof(u32) * (size_t)KB);
    memcpy(cur, cnt, sizeof(u32) * (size_t)KB);
    L64Item *A = malloc(sizeof(L64Item) * (size_t)N);
    for (u32 i = 0; i < N; i++) {
        u32 lo = i & 4095u, hi = i >> 12;
        u32 c1 = tlo_c1[r][lo] + thi_c1[r][hi];
        u32 s0 = tlo_s0[r][lo] + thi_s0[r][hi];
        u32 s1 = tlo_s1[r][lo] + thi_s1[r][hi];
        A[cur[c1 & KM]++] = (L64Item){ c1, s0, s1, i };
    }
    free(cur);
    *out_s = A; *out_cnt = cnt;
}

static void sort_m64_in_place(M64Item *M, u32 nM, M64Item **out_s, u32 **out_cnt) {
    const u32 KB = 1u << 25;
    u32 *cnt = calloc((size_t)(KB + 1u), sizeof(u32));
    for (u32 i = 0; i < nM; i++) {
        u32 k2 = (((u32)((M[i].pack >> 48) & 0x1FFu)) << 16) | (M[i].s0 & 0xFFFFu);
        cnt[k2 + 1u]++;
    }
    for (u32 b = 1; b <= KB; b++) cnt[b] += cnt[b - 1];
    u32 *cur = malloc(sizeof(u32) * (size_t)KB);
    memcpy(cur, cnt, sizeof(u32) * (size_t)KB);
    M64Item *Ms = malloc(sizeof(M64Item) * (size_t)nM);
    for (u32 i = 0; i < nM; i++) {
        u32 k2 = (((u32)((M[i].pack >> 48) & 0x1FFu)) << 16) | (M[i].s0 & 0xFFFFu);
        Ms[cur[k2]++] = M[i];
    }
    free(cur); free(M);
    *out_s = Ms; *out_cnt = cnt;
}

static int solve_wagner_64bit(const BinSlot *slots, u32 off_s0, u32 off_s1,
                              int nthreads, u8 out_bits[192]) {
    u32 (*tlo_c1)[4096] = malloc(sizeof(u32[8][4096]));
    u32 (*tlo_s0)[4096] = malloc(sizeof(u32[8][4096]));
    u32 (*tlo_s1)[4096] = malloc(sizeof(u32[8][4096]));
    u32 (*thi_c1)[4096] = malloc(sizeof(u32[8][4096]));
    u32 (*thi_s0)[4096] = malloc(sizeof(u32[8][4096]));
    u32 (*thi_s1)[4096] = malloc(sizeof(u32[8][4096]));

    for (int r = 0; r < 8; r++) {
        for (u32 m = 0; m < 4096u; m++) {
            u32 c1a = 0, s0a = (r == 7) ? off_s0 : 0u, s1a = (r == 7) ? off_s1 : 0u;
            u32 c1b = 0, s0b = 0, s1b = 0;
            for (int b = 0; b < 12; b++) {
                if ((m >> b) & 1u) {
                    c1a += slots[24 * r + b].dc1;      s0a += slots[24 * r + b].ds0;      s1a += slots[24 * r + b].ds1;
                    c1b += slots[24 * r + 12 + b].dc1; s0b += slots[24 * r + 12 + b].ds0; s1b += slots[24 * r + 12 + b].ds1;
                }
            }
            tlo_c1[r][m] = c1a; tlo_s0[r][m] = s0a; tlo_s1[r][m] = s1a;
            thi_c1[r][m] = c1b; thi_s0[r][m] = s0b; thi_s1[r][m] = s1b;
        }
    }

    const u32 CAP_M = 34500000u;
    const u32 CAP_N = 34500000u;
    pthread_t th[32];

    /* Step A: Build M_0 = L_0 + L_1 and sort by Level-2 25-bit key */
    L64Item *L0s; u32 *L0cnt;
    build_sorted_leaf64(0, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1, &L0s, &L0cnt);
    M64Item *M0raw = malloc(sizeof(M64Item) * (size_t)CAP_M);
    volatile u32 nM0 = 0;
    StreamL1Arg sl1[32];
    for (int t = 0; t < nthreads; t++) {
        sl1[t] = (StreamL1Arg){ t, nthreads, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1,
                                1, L0s, L0cnt, M0raw, &nM0, CAP_M };
        pthread_create(&th[t], NULL, stream_l1_worker, &sl1[t]);
    }
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    free(L0s); free(L0cnt);
    if (nM0 > CAP_M) nM0 = CAP_M;
    M64Item *M0s; u32 *M0cnt;
    sort_m64_in_place(M0raw, nM0, &M0s, &M0cnt);
    fprintf(stderr, "  [Wagner64] |M_0| = %u\n", nM0);

    /* Step B: Sort L_2, stream L_3 -> probe M_0 -> materialize N_0 */
    L64Item *L2s; u32 *L2cnt;
    build_sorted_leaf64(2, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1, &L2s, &L2cnt);
    N64Full *N0raw = malloc(sizeof(N64Full) * (size_t)CAP_N);
    volatile u32 nN0 = 0;
    StreamN0Arg sn0[32];
    for (int t = 0; t < nthreads; t++) {
        sn0[t] = (StreamN0Arg){ t, nthreads, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1,
                                L2s, L2cnt, M0s, M0cnt, N0raw, &nN0, CAP_N };
        pthread_create(&th[t], NULL, stream_n0_worker, &sn0[t]);
    }
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    free(L2s); free(L2cnt); free(M0s); free(M0cnt);
    if (nN0 > CAP_N) nN0 = CAP_N;
    fprintf(stderr, "  [Wagner64] |N_0| = %u\n", nN0);

    /* Step C: Sort N_0 by (s1 & ((1<<25)-1)) */
    const u32 H3B = 1u << 25;
    const u32 H3M = H3B - 1u;
    u32 *N0cnt = calloc((size_t)(H3B + 1u), sizeof(u32));
    for (u32 i = 0; i < nN0; i++) N0cnt[(N0raw[i].s1 & H3M) + 1u]++;
    for (u32 b = 1; b <= H3B; b++) N0cnt[b] += N0cnt[b - 1];
    u32 *cur3 = malloc(sizeof(u32) * (size_t)H3B);
    memcpy(cur3, N0cnt, sizeof(u32) * (size_t)H3B);
    N64Full *N0s = malloc(sizeof(N64Full) * (size_t)nN0);
    for (u32 i = 0; i < nN0; i++) N0s[cur3[N0raw[i].s1 & H3M]++] = N0raw[i];
    free(cur3); free(N0raw);

    /* Step D: Build M_2 = L_4 + L_5 and sort by Level-2 25-bit key */
    L64Item *L4s; u32 *L4cnt;
    build_sorted_leaf64(4, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1, &L4s, &L4cnt);
    M64Item *M2raw = malloc(sizeof(M64Item) * (size_t)CAP_M);
    volatile u32 nM2 = 0;
    for (int t = 0; t < nthreads; t++) {
        sl1[t] = (StreamL1Arg){ t, nthreads, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1,
                                5, L4s, L4cnt, M2raw, &nM2, CAP_M };
        pthread_create(&th[t], NULL, stream_l1_worker, &sl1[t]);
    }
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
    free(L4s); free(L4cnt);
    if (nM2 > CAP_M) nM2 = CAP_M;
    M64Item *M2s; u32 *M2cnt;
    sort_m64_in_place(M2raw, nM2, &M2s, &M2cnt);
    fprintf(stderr, "  [Wagner64] |M_2| = %u\n", nM2);

    /* Step E: Sort L_6, stream L_7 -> probe M_2 -> probe N_0 */
    L64Item *L6s; u32 *L6cnt;
    build_sorted_leaf64(6, tlo_c1, tlo_s0, tlo_s1, thi_c1, thi_s0, thi_s1, &L6s, &L6cnt);
    volatile int solved = 0;
    pthread_mutex_t sol_mu = PTHREAD_MUTEX_INITIALIZER;
    u32 shared_w6[6] = {0};
    StreamFinalArg sf[32];
    for (int t = 0; t < nthreads; t++) {
        sf[t].tid = t; sf[t].nthreads = nthreads;
        sf[t].tlo_c1 = tlo_c1; sf[t].tlo_s0 = tlo_s0; sf[t].tlo_s1 = tlo_s1;
        sf[t].thi_c1 = thi_c1; sf[t].thi_s0 = thi_s0; sf[t].thi_s1 = thi_s1;
        sf[t].L6s = L6s; sf[t].L6cnt = L6cnt;
        sf[t].M2s = M2s; sf[t].M2cnt = M2cnt;
        sf[t].N0s = N0s; sf[t].N0cnt = N0cnt;
        sf[t].solved = &solved; sf[t].sol_mu = &sol_mu;
        sf[t].shared_out_w6 = shared_w6;
        pthread_create(&th[t], NULL, stream_final_worker, &sf[t]);
    }
    for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);

    if (solved) {
        for (int b = 0; b < 192; b++) {
            out_bits[b] = (u8)((shared_w6[b >> 5] >> (b & 31)) & 1u);
        }
    }

    free(L6s); free(L6cnt); free(M2s); free(M2cnt); free(N0s); free(N0cnt);
    free(tlo_c1); free(tlo_s0); free(tlo_s1);
    free(thi_c1); free(thi_s0); free(thi_s1);
    return solved;
}

static size_t serialize_multi_anchor_msg(
    const u32 bridge_frags[8][9],
    const int d1[8],
    u32 *const s1_frags[8],
    WalkInfo *const anchor_walks[8],
    const BinSlot *chosen,
    int nslots,
    const u8 *bits, /* NULL => all 0 (reference message) */
    u32 **out_words)
{
    size_t total_words = 0;
    for (int m = 0; m < 8; m++) {
        total_words += 9 + (size_t)d1[m];
    }
    for (int i = 0; i < nslots; i++) {
        int m = chosen[i].anchor;
        int wi = (bits && bits[i]) ? chosen[i].v_idx : chosen[i].u_idx;
        total_words += (size_t)anchor_walks[m][wi].len;
    }

    u32 *buf = malloc(sizeof(u32) * total_words);
    size_t pos = 0;
    for (int m = 0; m < 8; m++) {
        memcpy(buf + pos, bridge_frags[m], 9 * sizeof(u32));
        pos += 9;
        if (d1[m] > 0) {
            memcpy(buf + pos, s1_frags[m], (size_t)d1[m] * sizeof(u32));
            pos += (size_t)d1[m];
        }
        for (int i = 0; i < nslots; i++) {
            if (chosen[i].anchor == m) {
                int wi = (bits && bits[i]) ? chosen[i].v_idx : chosen[i].u_idx;
                int len = anchor_walks[m][wi].len;
                memcpy(buf + pos, anchor_walks[m][wi].frags, (size_t)len * sizeof(u32));
                pos += (size_t)len;
            }
        }
    }
    *out_words = buf;
    return total_words;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <32|64> <target_hex> <out_msg.bin> [nthreads]\n", argv[0]);
        return 1;
    }
    int mode_bits = atoi(argv[1]);
    u64 target_val = strtoull(argv[2], NULL, 0);
    const char *out_path = argv[3];
    int nthreads = (argc > 4) ? atoi(argv[4]) : 8;
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 32) nthreads = 32;
    if (mode_bits != 32 && mode_bits != 64) {
        fprintf(stderr, "mode_bits must be 32 or 64\n");
        return 1;
    }

    double t_start = wall_sec();
    /*
     * 8 distinct shifted clean-subspace constants W'[0] in {a_m, ~a_m}:
     * Breaks BOTH the unshifted Sum[0] Archimedean bound AND the Single-Anchor
     * Multiplicative-Progression Trap (Delta S_0 == -(2*a_m+1)*Delta K_m mod 2^32)!
     */
    const u32 A_CONSTS[8] = {
        0x9E3779B9u, 0x7F4A7C15u, 0xD1B54A32u, 0x94D049BBu,
        0xBF58476Du, 0x1CE4E5B9u, 0xA24BAED4u, 0x9FB21C65u
    };
    int nslots_needed = (mode_bits == 32) ? 88 : 192;
    int pairs_per_anchor_needed = (nslots_needed / 8) + 4; /* 15 for Mode 32, 28 for Mode 64 */
    int init_walks_per_anchor   = (mode_bits == 32) ? 90 : 160;

    u32 bridge_frags[8][9];
    St ent_st[8], B[8];
    int d1[8] = {0}, L0[8] = {0};
    u32 *s1_frags[8] = {NULL};
    WalkInfo *anchor_walks[8] = {NULL};
    int anchor_nwalks[8] = {0};
    BinSlot *anchor_pairs[8] = {NULL};
    int anchor_npairs[8] = {0};

    double t_s1_tot = 0.0, t_s2_tot = 0.0;
    St cur_st;
    memcpy(cur_st.W, HIV, 32);
    cur_st.c1 = 0; cur_st.c2 = 0;

    for (int m = 0; m < 8; m++) {
        { int ok_b = 0; for (u32 w0 = 0; w0 < 256 && !ok_b; w0++) { St st1; int k1, k2; step32(cur_st.W, cur_st.c1, cur_st.c2, w0, st1.W, &k1, &k2); st1.c1 = k1; st1.c2 = k2; if (build_shifted_entry_try(&st1, A_CONSTS[m], &ent_st[m], bridge_frags[m] + 1, 8)) { bridge_frags[m][0] = w0; ok_b = 1; } } if (!ok_b) { fprintf(stderr, "bridge failed\n"); return 1; } }

        /* Stage 1 at anchor m */
        double ts1_0 = wall_sec();
        volatile int s1_done = 0;
        pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
        pthread_t th[32]; S1Arg s1a[32];
        for (int t = 0; t < nthreads; t++) {
            s1a[t] = (S1Arg){ ent_st[m], A_CONSTS[m], 12000, &s1_done, &mu,
                              &d1[m], &L0[m], &s1_frags[m],
                              0x9E3779B97F4A7C15ULL * (u64)(32 * m + t + 1) + 777ULL };
            pthread_create(&th[t], NULL, s1_worker, &s1a[t]);
        }
        for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);
        if (!s1_done) { fprintf(stderr, "Stage 1 failed at anchor %d!\n", m); return 1; }
        t_s1_tot += wall_sec() - ts1_0;

        St s = ent_st[m];
        for (int i = 0; i < d1[m]; i++) {
            u32 nw[8]; int k1, k2;
            step32(s.W, s.c1, s.c2, s1_frags[m][i], nw, &k1, &k2);
            memcpy(s.W, nw, 32); s.c1 = k1; s.c2 = k2;
        }
        B[m] = s;
        St *CS = malloc(sizeof(St) * (size_t)L0[m]);
        u32 *CF = s1_frags[m] + d1[m];
        for (int i = 0; i < L0[m]; i++) {
            CS[i] = s;
            u32 nw[8]; int k1, k2;
            step32(s.W, s.c1, s.c2, CF[i], nw, &k1, &k2);
            memcpy(s.W, nw, 32); s.c1 = k1; s.c2 = k2;
        }

        /* Stage 2 at anchor m: grow Sink-DAG and collect >= pairs_per_anchor_needed equal-length pairs */
        double ts2_0 = wall_sec();
        SinkDag dag; dag_init(&dag, 19, 320000);
        dag.nodes[0] = (DagNode){ B[m], 0, -1, 0 }; dag.n_nodes = 1;
        int nxt = 0;
        for (int i = L0[m] - 1; i >= 1; i--) nxt = dag_insert(&dag, &CS[i], L0[m] - i, nxt, CF[i]);
        {
            u32 q = (u32)(st_hash(&B[m]) & dag.mask);
            for (u32 p = 0; p <= dag.mask; p++) {
                u32 j = (q + p) & dag.mask;
                if (dag.slot_idx[j] == 0) { dag.slot_idx[j] = 1; break; }
            }
        }
        free(CS);

        int cap_walks = 600;
        anchor_walks[m] = calloc((size_t)cap_walks, sizeof(WalkInfo));
        volatile int n_walks = 0;
        eval_walk(&B[m], CF, L0[m], &anchor_walks[m][0]); n_walks = 1;

        pthread_mutex_t dag_mu = PTHREAD_MUTEX_INITIALIZER, w_mu = PTHREAD_MUTEX_INITIALIZER;
        int goal_walks = init_walks_per_anchor;
        while (1) {
            S2Arg s2a[32];
            for (int t = 0; t < nthreads; t++) {
                s2a[t] = (S2Arg){ B[m], A_CONSTS[m], &dag, &dag_mu, anchor_walks[m], &n_walks,
                                  goal_walks, &w_mu,
                                  0xBF58476D1CE4E5B9ULL * (u64)(64 * m + t + goal_walks) + 12345ULL };
                pthread_create(&th[t], NULL, s2_worker, &s2a[t]);
            }
            for (int t = 0; t < nthreads; t++) pthread_join(th[t], NULL);

            free(anchor_pairs[m]);
            anchor_pairs[m] = malloc(sizeof(BinSlot) * (size_t)n_walks);
            int *used = calloc((size_t)n_walks, sizeof(int));
            int np = 0;
            for (int i = 0; i < n_walks; i++) {
                if (used[i]) continue;
                for (int j = i + 1; j < n_walks; j++) {
                    if (!used[j] && anchor_walks[m][i].len == anchor_walks[m][j].len) {
                        u32 dc1 = anchor_walks[m][j].c1 - anchor_walks[m][i].c1;
                        u32 ds0 = anchor_walks[m][j].sum[0] - anchor_walks[m][i].sum[0];
                        u32 ds1 = anchor_walks[m][j].sum[1] - anchor_walks[m][i].sum[1];
                        if (dc1 != 0 || ds0 != 0 || ds1 != 0) {
                            used[i] = used[j] = 1;
                            anchor_pairs[m][np++] = (BinSlot){ m, i, j, dc1, ds0, ds1 };
                            break;
                        }
                    }
                }
            }
            free(used);
            anchor_npairs[m] = np;
            anchor_nwalks[m] = n_walks;
            if (np >= pairs_per_anchor_needed || goal_walks + 60 >= cap_walks) break;
            goal_walks += 60;
        }
        dag_free(&dag);
        t_s2_tot += wall_sec() - ts2_0;
        if (anchor_npairs[m] < (nslots_needed / 8)) {
            fprintf(stderr, "Not enough pairs at anchor %d (%d)\n", m, anchor_npairs[m]);
            return 1;
        }
        cur_st = B[m];
    }

    int tot_pairs = 0;
    for (int m = 0; m < 8; m++) tot_pairs += anchor_npairs[m];
    fprintf(stderr, "[Stage 1+2] 8 shifted anchors ready: total_pairs=%d (need %d), s1=%.2fs, s2=%.2fs\n",
            tot_pairs, nslots_needed, t_s1_tot, t_s2_tot);

    double t_s2_done = wall_sec();
    u8 win_bits[192] = {0};
    BinSlot chosen[192];
    int solved = 0;
    for (int attempt = 0; attempt < 4 && !solved; attempt++) {
        /* Interleave anchors m = i % 8 so EVERY leaf list spans all 8 multipliers -(2*a_m+1)! */
        int cursor[8];
        for (int m = 0; m < 8; m++) cursor[m] = attempt % 4;
        for (int i = 0; i < nslots_needed; i++) {
            int m = i & 7;
            int idx_in_m = (cursor[m]++) % anchor_npairs[m];
            chosen[i] = anchor_pairs[m][idx_in_m];
        }

        u32 *ref_words = NULL;
        size_t total_words = serialize_multi_anchor_msg(
            bridge_frags, d1, s1_frags, anchor_walks, chosen, nslots_needed, NULL, &ref_words);

        Blender256State bss; b256_init(&bss);
        b256_update_words(&bss, ref_words, total_words);
        u32 dref[8]; u8 dref_bytes[32];
        b256_final(&bss, dref, dref_bytes);
        free(ref_words);

        if (mode_bits == 32) {
            u32 t0 = (u32)(target_val & 0xFFFFFFFFULL);
            u32 off_s0 = dref[0] - t0;
            solved = solve_wagner_32bit(chosen, off_s0, win_bits);
        } else {
            u32 t0 = (u32)(target_val >> 32);
            u32 t1 = (u32)(target_val & 0xFFFFFFFFULL);
            u32 off_s0 = dref[0] - t0;
            u32 off_s1 = dref[1] - t1;
            solved = solve_wagner_64bit(chosen, off_s0, off_s1, nthreads, win_bits);
        }
        if (!solved) fprintf(stderr, "  [Wagner] Attempt %d missed Poisson draw, reshuffling...\n", attempt);
    }
    if (!solved) { fprintf(stderr, "Wagner solver failed!\n"); return 1; }
    double t_wagner = wall_sec();

    /* Assemble final winning message M(win_bits), hash from scratch, verify, and save! */
    u32 *win_words = NULL;
    size_t total_words = serialize_multi_anchor_msg(
        bridge_frags, d1, s1_frags, anchor_walks, chosen, nslots_needed, win_bits, &win_words);

    Blender256State bwin; b256_init(&bwin);
    b256_update_words(&bwin, win_words, total_words);
    u32 dwin[8]; u8 dwin_bytes[32];
    b256_final(&bwin, dwin, dwin_bytes);

    FILE *fp = fopen(out_path, "wb");
    if (!fp) { perror("fopen"); return 1; }
    for (size_t i = 0; i < total_words; i++) {
        u8 le[4] = { (u8)win_words[i], (u8)(win_words[i] >> 8),
                     (u8)(win_words[i] >> 16), (u8)(win_words[i] >> 24) };
        fwrite(le, 1, 4, fp);
    }
    fclose(fp);

    double t_end = wall_sec();
    printf("=== Fixed-IV %d-Bit Partial-Digest Preimage FOUND ===\n", mode_bits);
    printf("  Target %d-bit slice  : %s\n", mode_bits, argv[2]);
    printf("  Message length       : %zu words (%zu bytes = %.2f KB)\n",
           total_words, total_words * 4, (double)(total_words * 4) / 1024.0);
    printf("  Wall time breakdown  : Stage1=%.2fs, Stage2=%.2fs, Wagner=%.2fs, Total=%.2fs\n",
           t_s1_tot, t_s2_tot, t_wagner - t_s2_done, t_end - t_start);
    printf("  Full Blender-256 hex : ");
    for (int i = 0; i < 32; i++) printf("%02x", dwin_bytes[i]);
    printf("\n");

    if (mode_bits == 32) {
        u32 t0 = (u32)(target_val & 0xFFFFFFFFULL);
        if (dwin[0] != t0) { fprintf(stderr, "MISMATCH: %08x != %08x\n", dwin[0], t0); return 2; }
    } else {
        u32 t0 = (u32)(target_val >> 32), t1 = (u32)(target_val & 0xFFFFFFFFULL);
        if (dwin[0] != t0 || dwin[1] != t1) {
            fprintf(stderr, "MISMATCH: %08x_%08x != %08x_%08x\n", dwin[0], dwin[1], t0, t1);
            return 2;
        }
    }
    printf("  Internal C check     : [PASS] Exact %d-bit match from standard SHA-256 IV!\n", mode_bits);

    free(win_words);
    for (int m = 0; m < 8; m++) {
        for (int i = 0; i < anchor_nwalks[m]; i++) free(anchor_walks[m][i].frags);
        free(anchor_walks[m]);
        free(anchor_pairs[m]);
        free(s1_frags[m]);
    }
    return 0;
}
