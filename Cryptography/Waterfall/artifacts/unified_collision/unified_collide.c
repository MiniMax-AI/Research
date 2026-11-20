/*
 * unified_collide.c -- Unified Full-State Collision Solver for Waterfall-256
 *
 * Synthesises:
 *   1. Maximal 13-dimensional GF(2)-affine differential fixed-point coset
 *      T_beta = c + V (|T_beta| = 2^13 = 8192, V = H_24 + (F_2^8 << 24) with
 *      |H_24| = 32 and free top byte).
 *   2. Single-excursion in-window pool parity cancellation and active
 *      block-counter step steering (via closed-form Xbox^{-1} inversion and
 *      coset filtering), reducing the minimal full-state collision length from
 *      44 blocks (2816 B) and 9 blocks (576 B) down to the provably minimal
 *      7 blocks (119 steps = 448 bytes) and 8 blocks (136 steps = 512 bytes).
 *   3. Byte-cascade linear decomposition of the 4-round Xbox permutation,
 *      reducing every 2^13-coset scan from 8192 Xbox evaluations to 32
 *      low-24-bit coset steps with O(1) top-byte table lookup.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)) && defined(__BMI2__)
#include <immintrin.h>
#else
static inline uint32_t _pext_u32(uint32_t val, uint32_t mask) {
    uint32_t res = 0;
    for (uint32_t bb = 1; mask; bb <<= 1) {
        if (val & (mask & (0u - mask))) res |= bb;
        mask &= mask - 1;
    }
    return res;
}
static inline uint32_t _pdep_u32(uint32_t val, uint32_t mask) {
    uint32_t res = 0;
    for (uint32_t bb = 1; mask; bb <<= 1) {
        if (val & bb) res |= mask & (0u - mask);
        mask &= mask - 1;
    }
    return res;
}
#endif
#include "model.h"
#include "xdiff6.h"
#include PATTERN_FILE

#define BETA 0x00f9e70eu
#undef ROT
#define ROT(w) ((((uint32_t)(w)) << 24) | (((uint32_t)(w)) >> 8))
#define ROTL8(w) (((w) >> 24) | ((w) << 8))

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint32_t xbox_inv(uint32_t w) {
    for (int k = 0; k < 4; k++) {
        w = ROTL8(w);
        w ^= XB[SBINV[w & 0xFF]];
    }
    return w;
}

/* ---------------- T, V, sigma, and H_24 decomposition ---------------- */
static uint32_t VB[32];
static int PIV[32];
static int NV;
static uint32_t NPMASK, TAU;
static uint32_t *CO;       /* all 2^13 elements of V */
static uint32_t H24[32];   /* the 32 low-24-bit elements of H_24 (V = H24 + (F_2^8 << 24)) */
static uint32_t CSHIFT;    /* coset representative: T_beta = CSHIFT ^ V */

/* Byte-cascade linear filter tables for O(1) top-byte resolution */
static uint32_t PHI24[256];      /* phi(ds) = ROT(XB[i] ^ XB[i^e]) & 0xFFFFFF */
static int16_t  PHI_BY_SYN[1 << 19]; /* maps synidx(ROT(XB[i]^XB[i^e])) -> ds (or -1) */
static int16_t  ROTXB_BY_SYN[1 << 19]; /* maps synidx(ROT(XB[i])) -> i (or -1) */
static uint8_t  DDT_BYTE[256][256];    /* DDT_BYTE[e][ds] = an index i with sbox[i]^sbox[i^e]==ds, or 0xFF */
static uint8_t  DDT_HAS[256][256];

static inline uint32_t synrep(uint32_t x) {
    uint32_t s = x;
    for (int i = 0; i < NV; i++) s ^= VB[i] & (uint32_t)(-(int32_t)((x >> PIV[i]) & 1));
    return s;
}
static inline uint32_t synidx(uint32_t x) { return _pext_u32(synrep(x), NPMASK); }
static inline int      inV(uint32_t x)    { return synrep(x) == 0; }
static inline uint32_t fromsyn(uint32_t k){ return _pdep_u32(k, NPMASK); }

static void init_cascade_tables(void) {
    memset(PHI_BY_SYN, -1, sizeof(PHI_BY_SYN));
    memset(ROTXB_BY_SYN, -1, sizeof(ROTXB_BY_SYN));
    memset(DDT_HAS, 0, sizeof(DDT_HAS));
    for (int e = 0; e < 256; e++) {
        for (int i = 0; i < 256; i++) {
            int ds = sbox[i] ^ sbox[i ^ e];
            if (!DDT_HAS[e][ds]) {
                DDT_HAS[e][ds] = 1;
                DDT_BYTE[e][ds] = (uint8_t)i;
            }
        }
    }
    for (int i = 0; i < 256; i++) {
        uint32_t rx = ROT(XB[i]);
        ROTXB_BY_SYN[synidx(rx)] = (int16_t)i;
        for (int d = 0; d < 256; d++) {
            int ds = sbox[i] ^ sbox[i ^ d];
            uint32_t rxd = ROT(XB[i] ^ XB[i ^ d]);
            PHI24[ds] = rxd;
            PHI_BY_SYN[synidx(rxd)] = (int16_t)ds;
        }
    }
}

/* ---------------- state accessors ---------------- */
static hstate SA[PT + 1], SB[PT + 1];
static uint32_t X[PT];

static inline uint32_t Pa(const hstate *s, int t)  { return s->stream1[(t + 15) % 16] ^ s->stream1[t % 16]; }
static inline uint32_t Pb(const hstate *s, int t)  { return s->stream2[(t + 6) % 7]   ^ s->stream2[t % 7]; }
static inline uint32_t Pc(const hstate *s, int t)  { return s->stream3[(t + 5) % 6]   ^ s->stream3[t % 6]; }
static inline uint32_t KAf(const hstate *s, int t) { return s->stream1[(t + 1) % 16]  ^ s->stream2[(t + 1) % 7] ^ 1u; }
static inline uint32_t KCf(const hstate *s, int t) { return s->stream1[(t + 1) % 16]  ^ s->stream3[(t + 1) % 6] ^ 2u; }
static inline uint32_t KAPf(const hstate *s, int t){ return s->stream2[(t + 1) % 7]   ^ s->stream3[(t + 1) % 6] ^ 3u; }
static inline int actof(int s) { return (PU[s] ? 1 : 0) | (PV[s] ? 2 : 0) | (PW[s] ? 4 : 0); }
static inline int nact(int s)  { return PU[s] + PV[s] + PW[s]; }
static inline int iscounter(int t) { return t % 17 == 16; }

/* Check whether step t owns an x_t-independent 19-bit condition:
 * Either t is a non-counter step with >=2 active streams, OR t is a counter step with 1 active stream. */

static int ring_preset_target_ctr(int tk, int *stream_idx) {
    const int KS[3] = {16, 7, 6};
    for (int i = 0; i < 3; i++) {
        int tc = tk + KS[i];
        if (tc < PT && iscounter(tc) && actof(tc) == (1 << i) && actof(tc - 1) == (1 << i)) {
            *stream_idx = i;
            return tc;
        }
    }
    return 0;
}

static int is_locked_pre_ctr(int t) {
    return (t + 1 < PT && iscounter(t + 1) && nact(t + 1) == 1 && actof(t) == actof(t + 1));
}

static inline int has_xfree_cond(int t) {
    int s_idx;
    if (ring_preset_target_ctr(t, &s_idx)) return 1;
    if (iscounter(t) && nact(t) == 1 && actof(t - 1) == actof(t)) return 0;
    return iscounter(t) ? (nact(t) >= 1) : (nact(t) >= 2);
}

static void place(int t, uint32_t x) {
    X[t] = x;
    SA[t + 1] = SA[t]; ustep(&SA[t + 1], x);
    SB[t + 1] = SB[t]; ustep(&SB[t + 1], x ^ (PD[t] ? BETA : 0u));
}

static uint32_t shift_of(const hstate *s, int t) {
    switch (actof(t)) {
    case 2: return Pb(s, t) ^ 1u;
    case 4: return Pc(s, t) ^ 2u;
    case 6: return Pb(s, t) ^ 1u;
    default: return Pa(s, t);
    }
}

static int ATTEMPT = 0;
static inline uint32_t filler(int t, uint64_t n) {
    uint64_t h = 0x9E3779B97F4A7C15ull * (uint64_t)(t + 1)
               ^ (0xBF58476D1CE4E5B9ull * (n + 1))
               ^ (0xD6E8FEB86659FD93ull * (uint64_t)(ATTEMPT + 1));
    h ^= h >> 31; h *= 0x94D049BB133111EBull; h ^= h >> 29;
    return (uint32_t)h;
}

static int conds_ok(int t, uint32_t x) {
    const hstate *s = &SA[t];
    uint32_t ya = Pa(s, t) ^ x, yb = Pb(s, t) ^ x ^ 1u, yc = Pc(s, t) ^ x ^ 2u;
    int a = actof(t), ok = 1;
    if (a & 1) ok &= (synidx(ya) == TAU);
    if (a & 2) ok &= (synidx(yb) == TAU);
    if (a & 4) ok &= (synidx(yc) == TAU);
    return ok;
}

static int xfree_ok(int t) {
    if (iscounter(t)) {
        return conds_ok(t, (uint32_t)(t / 17 + 1));
    }
    const hstate *s = &SA[t];
    uint32_t A = Pa(s, t) ^ Pb(s, t) ^ 1u, C = Pa(s, t) ^ Pc(s, t) ^ 2u;
    switch (actof(t)) {
    case 3: return inV(A);
    case 5: return inV(C);
    case 6: return inV(A ^ C);
    case 7: return inV(A) && inV(C);
    default: return 1;
    }
}

/* ============================================================ */
/*  3-stream injection step                                     */
/* ============================================================ */
static int injection(int T0) {
    for (int t = 0; t <= T0 - 5; t++) {
        if (iscounter(t)) place(t, (uint32_t)(t / 17 + 1));
        else place(t, filler(t, 0));
    }

    int r4 = T0 - 4, r3 = T0 - 3, r2 = T0 - 2, r1 = T0 - 1;
    uint32_t G4 = Pa(&SA[r4], r4) ^ Pb(&SA[r4], r4) ^ 1u ^ (Pa(&SA[r4], r4) ^ Pc(&SA[r4], r4) ^ 2u);
    uint32_t A4 = Pa(&SA[r4], r4) ^ Pb(&SA[r4], r4) ^ 1u;
    uint32_t kap4 = KAPf(&SA[r4], r4);
    uint32_t kap3 = KAPf(&SA[r4], r3);

    uint32_t z4 = 0; int done = 0;
    unsigned long long ngood = xd_ngood(kap3);
    for (unsigned long long n = 0; n < ngood && n < 4096; n++) {
        uint32_t Gg = xd_goodG(kap3, n);
        xchar c;
        if (!xd_prep(G4, Gg ^ kap4, &c)) continue;
        uint32_t zeta = xd_word(&c, 0);
        if (PHIW(zeta, G4) != (Gg ^ kap4)) continue;
        z4 = zeta ^ A4;
        done = 1; break;
    }
    if (!done) return 0;
    place(r4, z4 ^ Pa(&SA[r4], r4));

    uint32_t G3 = Pa(&SA[r3], r3) ^ Pb(&SA[r3], r3) ^ 1u ^ (Pa(&SA[r3], r3) ^ Pc(&SA[r3], r3) ^ 2u);
    uint32_t A3 = Pa(&SA[r3], r3) ^ Pb(&SA[r3], r3) ^ 1u;
    uint32_t kap3b = KAPf(&SA[r3], r3);
    xchar c3;
    if (!xd_prep(G3, kap3b, &c3)) return 0;
    uint32_t eta = xd_word(&c3, 0);
    if (PHIW(eta, G3) != kap3b) return 0;
    place(r3, (eta ^ A3) ^ Pa(&SA[r3], r3));
    if (!inV(Pa(&SA[r2], r2) ^ Pb(&SA[r2], r2) ^ 1u ^ (Pa(&SA[r2], r2) ^ Pc(&SA[r2], r2) ^ 2u)))
        return 0;

    uint32_t D    = Pa(&SA[r2], r2) ^ Pb(&SA[r2], r2) ^ 1u;
    uint32_t g1   = KAPf(&SA[r2], r2);
    uint32_t kap1 = KAPf(&SA[r2], r1);
    uint32_t kc1  = KCf(&SA[r2], r1);
    uint32_t kc2  = KCf(&SA[r2], r2);
    uint32_t Z = 0, Cp = 0, z2 = 0; int got = 0;
    for (int k = 0; k < (1 << 13) && !got; k++) {
        xchar c1;
        if (!xd_prep(g1, kap1 ^ CO[k], &c1)) continue;
        Z = xd_word(&c1, 0);
        if (synidx(PHIW(Z, g1) ^ kap1) != 0) continue;
        uint32_t zp = Z ^ g1;
        for (int k2 = 0; k2 < (1 << 13); k2++) {
            uint32_t cand;
            if (!xd_diff(zp, kc1 ^ CO[k2], &cand)) continue;
            if (synidx(PHIW(zp, cand) ^ kc1) != 0) continue;
            xchar c2;
            if (!xd_prep(D, cand ^ kc2, &c2)) continue;
            uint32_t zz = xd_word(&c2, 0);
            if (PHIW(zz, D) != (cand ^ kc2)) continue;
            Cp = cand; z2 = zz; got = 1; break;
        }
    }
    if (!got) return 0;

    uint32_t A1 = Cp ^ g1;
    uint32_t z1 = Z ^ A1;
    place(r2, z2 ^ Pa(&SA[r2], r2));
    place(r1, z1 ^ Pa(&SA[r1], r1));

    if (!xfree_ok(T0)) return 0;
    uint32_t x0 = fromsyn(TAU ^ synidx(shift_of(&SA[T0], T0)));
    place(T0, x0);
    return conds_ok(T0, X[T0]);
}

/* ============================================================ */
/*  Unified walk: 2-stream steps + 1-stream counter steps       */
/*  with 32-element byte-cascade O(1) top-byte projection       */
/* ============================================================ */
static long SWEEP;
static long SWEEPM[5];
static int INJ = -1;

static int place_plain(int r, uint64_t n) {
    if (iscounter(r)) {
        uint32_t ctr = (uint32_t)(r / 17 + 1);
        if (!conds_ok(r, ctr)) return 0;
        place(r, ctr);
        return 1;
    }
    if (!xfree_ok(r)) return 0;
    uint32_t x;
    if (is_locked_pre_ctr(r)) x = shift_of(&SA[r], r) ^ CSHIFT;
    else if (actof(r) == 0)   x = filler(r, n);
    else                      x = fromsyn(TAU ^ synidx(shift_of(&SA[r], r))) ^ CO[n & ((1 << 13) - 1)];
    place(r, x);
    return conds_ok(r, x);
}

static int pick_k(int t, int *mode) {
    int s_idx;
    if (ring_preset_target_ctr(t, &s_idx)) {
        if (nact(t) == 0) { *mode = 1; return 1; }
        for (int k = 1; k <= 3; k++) {
            int s = t - k;
            if (s < 0 || iscounter(s) || s == INJ) continue;
            int ok = 1;
            for (int r = s; r < t; r++) if (has_xfree_cond(r)) { ok = 0; break; }
            if (ok) { *mode = 2; return k; }
        }
        return 0;
    }
    if (t < 1 || iscounter(t - 1)) return 0;
    if (nact(t - 1) == 0) { *mode = 1; return 1; }
    for (int k = 2; k <= 3; k++) {
        int s = t - k;
        if (s < 0 || iscounter(s) || s == INJ) continue;
        int ok = 1;
        for (int r = s; r < t; r++) if (has_xfree_cond(r)) { ok = 0; break; }
        if (ok) { *mode = 2; return k; }
    }
    for (int k = 2; k <= 3; k++) {
        int s = t - k;
        if (s < 0 || iscounter(s) || s == INJ) continue;
        int ok = 1;
        for (int r = s + 1; r < t; r++) if (has_xfree_cond(r)) { ok = 0; break; }
        if (ok) { *mode = 3; return k; }
    }
    if (!has_xfree_cond(t - 1)) { *mode = 4; return 1; }
    return 0;
}

static int solve_z(uint32_t base, uint32_t K, uint32_t *out) {
    for (int j = 0; j < (1 << 13); j++) {
        xchar c;
        if (!xd_prep(base, K ^ CO[j], &c)) continue;
        uint32_t z = xd_word(&c, 0);
        if (PHIW(z, base) != (K ^ CO[j])) continue;
        *out = z; return 1;
    }
    return 0;
}

/* Fast 32-step byte-cascade coset solver for 2-stream step t:
 * Instead of scanning 2^13 = 8192 elements of z in zb ^ V, iterates only over
 * the 32 low-24-bit values z0 = (zb & 0xFFFFFF) ^ H24[m] and resolves the top
 * byte b3 in O(1) via PHI_BY_SYN and DDT_BYTE! */
static int scan_coset_fast(int t, int kind, uint32_t base, uint32_t K) {
    uint32_t A = Pa(&SA[t - 1], t - 1) ^ Pb(&SA[t - 1], t - 1) ^ 1u;
    uint32_t Pa_t1 = Pa(&SA[t - 1], t - 1);
    uint32_t zb = Pa_t1 ^ fromsyn(TAU ^ synidx(shift_of(&SA[t - 1], t - 1)));
    if (kind == 6) zb ^= A;
    uint32_t zb_low = zb & 0xFFFFFFu;
    uint8_t base_b3 = (uint8_t)(base >> 24);

    for (int m = 0; m < 32; m++) {
        uint32_t z0 = zb_low ^ H24[m];
        /* Evaluate rounds 1..3 of xboxf(z0) and xboxf(z0 ^ base) with top byte of z0 = 0 */
        uint32_t w0 = z0, w1 = z0 ^ base;
        for (int r = 0; r < 3; r++) {
            w0 ^= XB[w0 & 0xFF]; w0 = ROT(w0);
            w1 ^= XB[w1 & 0xFF]; w1 = ROT(w1);
        }
        uint8_t c3_0 = (uint8_t)(w0 & 0xFF);
        uint8_t c3_1 = (uint8_t)(w1 & 0xFF);
        uint8_t e4 = c3_0 ^ c3_1;
        uint32_t phi_base = ROT(w0) ^ ROT(w1);
        uint32_t need_syn = synidx(phi_base ^ K);
        int16_t ds4 = PHI_BY_SYN[need_syn];
        if (ds4 >= 0 && DDT_HAS[e4][ds4]) {
            uint8_t i4 = DDT_BYTE[e4][ds4];
            uint8_t b3 = i4 ^ c3_0;
            uint32_t z = z0 | ((uint32_t)b3 << 24);
            if (inV(PHIW(z, base) ^ K)) {
                place(t - 1, ((kind == 6) ? (z ^ A) : z) ^ Pa_t1);
                return 1;
            }
        }
    }

    return 0;
}

/* Steering a 1-stream counter step t (t % 17 == 16, actof(t) in {1, 2, 4}):
 *   - Mode 1 (nact(t-1) == 0): O(1) exact inversion via xbox_inv!
 *   - Mode 2/3/4 (nact(t-1) == 1): 32-step byte-cascade lookup via ROTXB_BY_SYN! */
static int steer_counter_step(int t, int mode) {
    uint32_t ctr = (uint32_t)(t / 17 + 1);
    int act_t = actof(t);
    uint32_t K_ctr, P_t1;
    if (act_t == 1) {
        K_ctr = SA[t - 1].stream1[t % 16] ^ ctr;
        P_t1  = Pa(&SA[t - 1], t - 1);
    } else if (act_t == 2) {
        K_ctr = SA[t - 1].stream2[t % 7] ^ ctr ^ 1u;
        P_t1  = Pb(&SA[t - 1], t - 1) ^ 1u;
    } else {
        K_ctr = SA[t - 1].stream3[t % 6] ^ ctr ^ 2u;
        P_t1  = Pc(&SA[t - 1], t - 1) ^ 2u;
    }

    if (mode == 1) {
        /* x_{t-1} is completely unconstrained: invert Xbox in O(1)! */
        uint32_t Y_t1 = xbox_inv(K_ctr ^ CSHIFT);
        place(t - 1, Y_t1 ^ P_t1);
        return conds_ok(t, ctr);
    }

    /* nact(t-1) == 1: x_{t-1} must lie in xb_base ^ V.
     * In terms of Y_i(t-1) = P_t1 ^ x_{t-1}, Y_i(t-1) lies in yb ^ V where
     * yb = P_t1 ^ fromsyn(TAU ^ synidx(shift_of(&SA[t-1], t-1))).
     * For each of the 32 low-24-bit values y0 = (yb & 0xFFFFFF) ^ H24[m],
     * xboxf(y0 | (b3 << 24)) = ROT(w0) ^ ROT(XB[i4]) with i4 = b3 ^ (w0 & 0xFF),
     * so ROTXB_BY_SYN resolves i4 (and hence b3) in O(1)! */
    uint32_t xb_base = fromsyn(TAU ^ synidx(shift_of(&SA[t - 1], t - 1)));
    uint32_t yb_low  = (P_t1 ^ xb_base) & 0xFFFFFFu;
    for (int m = 0; m < 32; m++) {
        uint32_t y0 = yb_low ^ H24[m];
        uint32_t w0 = y0;
        for (int r = 0; r < 3; r++) {
            w0 ^= XB[w0 & 0xFF]; w0 = ROT(w0);
        }
        uint8_t c3 = (uint8_t)(w0 & 0xFF);
        uint32_t w_base = ROT(w0);
        uint32_t need_syn = TAU ^ synidx(w_base ^ K_ctr);
        int16_t i4 = ROTXB_BY_SYN[need_syn];
        if (i4 >= 0) {
            uint8_t b3 = (uint8_t)i4 ^ c3;
            uint32_t Y_t1 = y0 | ((uint32_t)b3 << 24);
            uint32_t x_t1 = Y_t1 ^ P_t1;
            if (conds_ok(t - 1, x_t1)) {
                place(t - 1, x_t1);
                if (conds_ok(t, ctr)) return 1;
            }
        }
    }
    return 0;
}


static int steer_ring_preset(int tk, int tc, int stream_i) {
    uint32_t ctr = (uint32_t)(tc / 17 + 1);
    uint32_t K_target = xboxf(CSHIFT) ^ ctr ^ (uint32_t)stream_i;
    uint32_t P_tk = (stream_i == 0) ? Pa(&SA[tk], tk)
                  : (stream_i == 1) ? (Pb(&SA[tk], tk) ^ 1u)
                                    : (Pc(&SA[tk], tk) ^ 2u);
    if (nact(tk) == 0) {
        uint32_t Y_tk = xbox_inv(K_target ^ CSHIFT);
        place(tk, Y_tk ^ P_tk);
        return 1;
    }
    uint32_t xb_base = fromsyn(TAU ^ synidx(shift_of(&SA[tk], tk)));
    uint32_t yb_low  = (P_tk ^ xb_base) & 0xFFFFFFu;
    for (int m = 0; m < 32; m++) {
        uint32_t y0 = yb_low ^ H24[m];
        uint32_t w0 = y0;
        for (int r = 0; r < 3; r++) {
            w0 ^= XB[w0 & 0xFF]; w0 = ROT(w0);
        }
        uint8_t c3 = (uint8_t)(w0 & 0xFF);
        uint32_t w_base = ROT(w0);
        uint32_t need_syn = TAU ^ synidx(w_base ^ K_target);
        int16_t i4 = ROTXB_BY_SYN[need_syn];
        if (i4 >= 0) {
            uint8_t b3 = (uint8_t)i4 ^ c3;
            uint32_t Y_tk = y0 | ((uint32_t)b3 << 24);
            uint32_t x_tk = Y_tk ^ P_tk;
            if (conds_ok(tk, x_tk) && synidx(xboxf(Y_tk) ^ K_target) == TAU) {
                place(tk, x_tk);
                return 1;
            }
        }
    }
    return 0;
}

static int steer(int t, int k, int mode) {
    int ring_s_idx;
    int ring_tc = ring_preset_target_ctr(t, &ring_s_idx);
    if (ring_tc) {
        int s = t - k;
        if (mode == 1) return steer_ring_preset(t, ring_tc, ring_s_idx);
        uint64_t omax = (actof(s) == 0) ? (1ull << 16) : (1ull << 14);
        for (uint64_t n = 0; n < omax; n++) {
            SWEEP++;
            uint32_t xs = (actof(s) == 0) ? filler(s, n) : fromsyn(TAU ^ synidx(shift_of(&SA[s], s))) ^ CO[n & ((1 << 13) - 1)];
            place(s, xs);
            if (!conds_ok(s, xs) || !xfree_ok(s)) continue;
            int bad = 0;
            for (int r = s + 1; r < t; r++) if (!place_plain(r, 0)) { bad = 1; break; }
            if (bad) continue;
            if (steer_ring_preset(t, ring_tc, ring_s_idx)) return 1;
        }
        return 0;
    }
    int s = t - k;
    if (s < 0 || iscounter(s) || s == INJ) return 0;
    int is_ctr_target = iscounter(t);
    int kind = actof(t);

    if (mode == 1) {
        if (is_ctr_target) {
            if (!steer_counter_step(t, 1)) return 0;
            SWEEPM[1]++;
            return 1;
        }
        uint32_t A = Pa(&SA[t - 1], t - 1) ^ Pb(&SA[t - 1], t - 1) ^ 1u;
        uint32_t C = Pa(&SA[t - 1], t - 1) ^ Pc(&SA[t - 1], t - 1) ^ 2u;
        uint32_t base = (kind == 3) ? A : (kind == 5) ? C : A ^ C;
        uint32_t K    = (kind == 3) ? KAf(&SA[t - 1], t - 1) : (kind == 5) ? KCf(&SA[t - 1], t - 1)
                                                                            : KAPf(&SA[t - 1], t - 1);
        uint32_t z;
        if (!solve_z(base, K, &z)) return 0;
        SWEEPM[1]++;
        place(t - 1, ((kind == 6) ? (z ^ A) : z) ^ Pa(&SA[t - 1], t - 1));
        return 1;
    }
    if (mode == 4) {
        SWEEP++;
        if (is_ctr_target) {
            if (!steer_counter_step(t, 4)) return 0;
            SWEEPM[4]++;
            return 1;
        }
        uint32_t A = Pa(&SA[t - 1], t - 1) ^ Pb(&SA[t - 1], t - 1) ^ 1u;
        uint32_t C = Pa(&SA[t - 1], t - 1) ^ Pc(&SA[t - 1], t - 1) ^ 2u;
        uint32_t base = (kind == 3) ? A : (kind == 5) ? C : A ^ C;
        uint32_t K    = (kind == 3) ? KAf(&SA[t - 1], t - 1) : (kind == 5) ? KCf(&SA[t - 1], t - 1)
                                                                            : KAPf(&SA[t - 1], t - 1);
        if (!scan_coset_fast(t, kind, base, K)) return 0;
        SWEEPM[4]++;
        return 1;
    }

    uint64_t outer_max = (mode == 3) ? (1ull << 16)
                       : (actof(s) == 0) ? (1ull << 16) : (1ull << 14);
    for (uint64_t n = 0; n < outer_max; n++) {
        SWEEP++;
        uint32_t xs = (actof(s) == 0) ? filler(s, n)
                    : fromsyn(TAU ^ synidx(shift_of(&SA[s], s))) ^ CO[n & ((1 << 13) - 1)];
        place(s, xs);
        if (!conds_ok(s, xs)) continue;
        if (!xfree_ok(s))     continue;
        int r, bad = 0;
        for (r = s + 1; r < t; r++) if (!place_plain(r, 0)) { bad = 1; break; }
        if (bad) continue;
        if (is_ctr_target) {
            if (steer_counter_step(t, mode)) { SWEEPM[mode]++; return 1; }
        } else {
            uint32_t A = Pa(&SA[t - 1], t - 1) ^ Pb(&SA[t - 1], t - 1) ^ 1u;
            uint32_t C = Pa(&SA[t - 1], t - 1) ^ Pc(&SA[t - 1], t - 1) ^ 2u;
            uint32_t base = (kind == 3) ? A : (kind == 5) ? C : A ^ C;
            uint32_t K    = (kind == 3) ? KAf(&SA[t - 1], t - 1) : (kind == 5) ? KCf(&SA[t - 1], t - 1)
                                                                                : KAPf(&SA[t - 1], t - 1);
            if (scan_coset_fast(t, kind, base, K)) { SWEEPM[mode]++; return 1; }
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    int start_attempt = (argc > 1) ? atoi(argv[1]) : 0;
    int hbits = (argc > 2) ? atoi(argv[2]) : 256;
    double t0 = now();
    xd_init();

    xchar ct;
    if (!xd_prep(BETA, BETA, &ct)) return 1;
    long nT = (long)ct.total;
    uint32_t *Tl = malloc(sizeof(uint32_t) * nT);
    for (long n = 0; n < nT; n++) Tl[n] = xd_word(&ct, (unsigned long long)n);

    uint32_t c0 = Tl[0]; NV = 0;
    for (long i = 1; i < nT; i++) {
        uint32_t x = Tl[i] ^ c0;
        for (int j = 0; j < NV; j++) if ((x >> PIV[j]) & 1) x ^= VB[j];
        if (x) {
            int p = __builtin_ctz(x); VB[NV] = x; PIV[NV] = p;
            for (int j = 0; j < NV; j++) if ((VB[j] >> p) & 1) VB[j] ^= x;
            NV++;
        }
    }
    NPMASK = 0;
    for (int b = 0; b < 32; b++) {
        int isp = 0;
        for (int i = 0; i < NV; i++) if (PIV[i] == b) isp = 1;
        if (!isp) NPMASK |= 1u << b;
    }
    TAU = synidx(c0);
    CO = malloc(sizeof(uint32_t) * (1 << 13));
    CO[0] = 0;
    for (int i = 0; i < NV; i++)
        for (int j = 0; j < (1 << i); j++) CO[(1 << i) + j] = CO[j] ^ VB[i];
    CSHIFT = fromsyn(TAU);

    /* Extract the 32 low-24-bit elements H24 of V */
    int nh24 = 0;
    for (int i = 0; i < (1 << 13); i++) {
        if ((CO[i] >> 24) == 0) H24[nh24++] = CO[i];
    }
    if (nh24 != 32) { fprintf(stderr, "unexpected |H24|=%d\n", nh24); return 1; }
    free(Tl);

    init_cascade_tables();

    int T0 = -1;
    for (int t = 0; t < PT; t++) if (PU[t] && PV[t] && PW[t]) T0 = t;
    INJ = T0;

    double t4 = now();
    int nsteer = 0, tries = 0, built = 0;
    for (ATTEMPT = start_attempt; ATTEMPT < start_attempt + 4000 && !built; ATTEMPT++) {
        tries++;
        hinit(&SA[0], hbits); hinit(&SB[0], hbits);
        if (!injection(T0)) continue;
        nsteer = 0; int ok = 1;
        for (int t = T0 + 1; t < PT; t++) {
            int ring_s_idx;
            if (ring_preset_target_ctr(t, &ring_s_idx)) {
                int mode = 0, k = pick_k(t, &mode);
                if (!k || !steer(t, k, mode)) { ok = 0; break; }
                nsteer++;
                continue;
            }
            if (has_xfree_cond(t)) {
                int mode = 0, k = pick_k(t, &mode);
                if (!k || !steer(t, k, mode)) { ok = 0; break; }
                nsteer++;
            }
            if (!xfree_ok(t) || !place_plain(t, 0)) { ok = 0; break; }
        }
        if (ok) built = 1;
    }
    if (!built) { fprintf(stderr, "sweep failed\n"); return 1; }

    int same = !memcmp(SA[PT].stream1, SB[PT].stream1, sizeof SA[PT].stream1)
            && !memcmp(SA[PT].stream2, SB[PT].stream2, sizeof SA[PT].stream2)
            && !memcmp(SA[PT].stream3, SB[PT].stream3, sizeof SA[PT].stream3)
            && !memcmp(SA[PT].pool2,   SB[PT].pool2,   sizeof SA[PT].pool2)
            && !memcmp(SA[PT].pool3,   SB[PT].pool3,   sizeof SA[PT].pool3);
    if (!same) { fprintf(stderr, "state mismatch\n"); return 1; }

    static unsigned char m1[NBLK * 64], m2[NBLK * 64];
    for (int t = 0; t < PT; t++) {
        if (iscounter(t)) continue;
        int off = (t / 17) * 64 + (t % 17) * 4;
        uint32_t a = X[t], b = X[t] ^ (PD[t] ? BETA : 0u);
        m1[off] = a; m1[off + 1] = a >> 8; m1[off + 2] = a >> 16; m1[off + 3] = a >> 24;
        m2[off] = b; m2[off + 1] = b >> 8; m2[off + 2] = b >> 16; m2[off + 3] = b >> 24;
    }
    unsigned char d1[64], d2[64];
    hhash(hbits, m1, sizeof m1, d1); hhash(hbits, m2, sizeof m2, d2);
    if (memcmp(d1, d2, hbits / 8) || !memcmp(m1, m2, sizeof m1)) return 1;

    fprintf(stderr, "[%d blocks = %d B, PT=%d, T0=%d] solved in %.4f s (walk %.4f s, %d attempts, %ld outer trials)\n",
            NBLK, NBLK * 64, PT, T0, now() - t0, now() - t4, tries, SWEEP);
    fprintf(stderr, "Waterfall-%d digest: ", hbits);
    for (int i = 0; i < hbits / 8; i++) fprintf(stderr, "%02x", d1[i]);
    fprintf(stderr, "\n");
    for (size_t i = 0; i < sizeof m1; i++) printf("%02x", m1[i]);
    printf("\n");
    for (size_t i = 0; i < sizeof m2; i++) printf("%02x", m2[i]);
    printf("\n");
    return 0;
}
