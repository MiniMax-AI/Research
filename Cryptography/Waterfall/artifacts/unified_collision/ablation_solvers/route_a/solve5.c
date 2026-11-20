/* Realise the differential excursion twice (544 steps apart) and print two
 * colliding messages as hex.
 *
 * Model.  The message difference at step t is EV[t]*DIFF and every stream
 * carries exactly DIFF while it is active, so with
 *      a_i[t] = s_i[t-1] ^ s_i[t-K_i] ^ m[t] ^ c_i
 * we get  Delta a_i[t] = b_i[t]*DIFF  automatically, and the only requirement is
 *      a_i[t] in F_DIFF(DIFF)   for every active stream i.
 * The top byte of such an element is free, so that is the 24-bit test
 * (a & 0xffffff) in S32.  F = S32 with free top byte is affine, F^F = H, so a
 * step with w active streams costs w-1 conditions "gamma in H" (19 bits each)
 * when the message word can be chosen, and w when it cannot (counter steps).
 *
 * Injection (the only expensive step).  The first message difference activates
 * all three streams at once, so two conditions (38 bits) must hold at the same
 * time on a single 32-bit unknown: with v the X-input of stream 1 one step
 * earlier,  low24(X(v) ^ X(v^D_i) ^ T_i) in HH  for i = 1,2.  Brute force costs
 * 2^32 per instance and ~2^6 instances.  Instead this uses the byte cascade of
 * the X-box (see gen_data3.py):
 *
 *      low24( X(v) ^ X(v^D) ) = c + phi(y3),
 *      c = ( y0^L3(y1)^L2(y2), y1^L3(y2), y2 ),   phi(y) = (L1 y, L2 y, L3 y)
 *
 * where (y0,y1,y2) depends on the low three bytes of v only and y3 on the top
 * byte only.  Membership in HH therefore implies the linear filter
 *      M c = M T,   ker M = U = HH + im(phi)   (dim U = 13, M has 11 rows),
 * and for fixed low 16 bits of v that filter is 11 equations in the 8 bits of
 * y2, resolved by the table SOLVE[].  The top byte then follows in closed form
 * from the projection LP (U = HH + im(phi) is a direct sum).  Cost per instance
 * drops from 2^32 X-box evaluations to 2^16 table steps: the search exponent
 * falls from 38 to 22.
 *
 * The remaining conditions (w = 2 steps later in the excursion) are single
 * 19-bit tests and are solved by re-randomising a minimal window of message
 * words in front of them (>= 24 bits of freedom, see PLAN column 6 in fq_data2.h:
 * wfull is the window materialised once, with redraw start umin = wfull - 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "fq_data2.h"

#define OFF 20
static const int KK[3] = {16, 7, 6};
static const uint32_t CC[3] = {0u, 1u, 2u};

static uint32_t SA[2][3][NST + OFF];
static uint32_t WD[2][NST];
static uint32_t DM[NST];
static uint32_t *bitF, *bitH;
/* DDTP[d][y] packs (up to four) x with SBOXC[x]^SBOXC[x^d] == y, 0xff padded */
static uint32_t DDTP[256][256];

static uint64_t rs = 0x9e3779b97f4a7c15ULL;
static inline uint32_t rnd32(void) {
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (uint32_t)(rs >> 32);
}
static inline uint32_t Xf(uint32_t p) {
    for (int i = 0; i < 4; i++) { p ^= XBOXC[p & 0xFF]; p = (p << 24) | (p >> 8); }
    return p;
}
static inline int inF(uint32_t v) { v &= 0xFFFFFFu; return (bitF[v >> 5] >> (v & 31)) & 1u; }
static inline int inH(uint32_t v) { v &= 0xFFFFFFu; return (bitH[v >> 5] >> (v & 31)) & 1u; }
static inline uint32_t randF(void) { return S32[rnd32() % NS32] | ((rnd32() & 0xFF) << 24); }

static void step(int m, int t) {
    uint32_t in = WD[m][t];
    SA[m][0][t + OFF] = Xf(SA[m][0][t - 1 + OFF] ^ SA[m][0][t - 16 + OFF] ^ in);
    SA[m][1][t + OFF] = Xf(SA[m][1][t - 1 + OFF] ^ SA[m][1][t - 7 + OFF] ^ in ^ 1u);
    SA[m][2][t + OFF] = Xf(SA[m][2][t - 1 + OFF] ^ SA[m][2][t - 6 + OFF] ^ in ^ 2u);
}
static void rerun(int from, int to) {
    if (from < 0) from = 0;
    for (int t = from; t <= to; t++) { step(0, t); step(1, t); }
}
static void setword(int t, uint32_t v) { WD[0][t] = v; WD[1][t] = v ^ DM[t]; }
static inline uint32_t Ai(int i, int t) {
    return SA[0][i][t - 1 + OFF] ^ SA[0][i][t - KK[i] + OFF] ^ CC[i];
}
static inline uint32_t gmm(int i, int j, int t) { return Ai(i, t) ^ Ai(j, t); }
static inline int isctr(int t) { return (t % 17) == 16; }

static const int *BP[3] = {B1, B2, B3};
static int wof(int rel) { return (rel >= 0 && rel < TT) ? BP[0][rel] + BP[1][rel] + BP[2][rel] : 0; }

/* re-randomise steps a..b of an excursion: an inactive step gets a free word, a
 * step with an active stream gets m = A_i ^ f with f in F (always valid, and
 * the extra conditions of a w >= 2 step are relations between the A values that
 * do not involve this word), a counter step is left alone */
static void replay(int start, int a, int b) {
    for (int t = a; t <= b; t++) {
        int rel = t - start;
        if (rel >= 0 && rel < PT && !isctr(t)) {
            int w = wof(rel);
            if (w == 0) {
                setword(t, rnd32());
            } else {
                int i = BP[0][rel] ? 0 : (BP[1][rel] ? 1 : 2);
                setword(t, Ai(i, t) ^ randF());
            }
        }
        rerun(t, t);
    }
}

/* verify the whole excursion: Delta s_i[t] must be b_i[rel]*DIFF */
static void chkall(int start) {
    for (int t = start - 2; t < start + PT + 6 && t < NST; t++) {
        int rel = t - start;
        for (int i = 0; i < 3; i++) {
            int act = (rel >= 0 && rel < TT) ? BP[i][rel] : 0;
            uint32_t want = act ? DIFF : 0u;
            uint32_t got = SA[0][i][t + OFF] ^ SA[1][i][t + OFF];
            if (got != want) {
                fprintf(stderr, "  !! chkall start=%d t=%d rel=%d stream=%d want=%08x got=%08x\n",
                        start, t, rel, i, want, got);
                exit(3);
            }
        }
    }
}

static void chk(int start, int t) {
    int rel = t - start;
    for (int i = 0; i < 3; i++) {
        int act = (rel >= 0 && rel < TT) ? BP[i][rel] : 0;
        uint32_t want = act ? DIFF : 0u;
        uint32_t got = SA[0][i][t + OFF] ^ SA[1][i][t + OFF];
        if (got != want) {
            uint32_t a0 = Ai(i, t) ^ WD[0][t];
            fprintf(stderr, "  !! start=%d t=%d rel=%d stream=%d want=%08x got=%08x inF=%d\n",
                    start, t, rel, i, want, got, inF(a0));
            exit(3);
        }
    }
}

/* ------------------------------------------------------------- injection */
/* Find v with low24(X(v) ^ X(v^D[i]) ^ T[i]) in HH for i = 0,1.
 * 2^16 table steps: the two low bytes of v are enumerated, the filter
 * M c = M T fixes y2 (hence the third byte) and the projection LP fixes y3
 * (hence the top byte).  Returns 1 on success. */
static int inj_scan(const uint32_t D[2], const uint32_t T[2], uint32_t *pv) {
    uint32_t dl[2][4], sig[2];
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 4; k++) dl[i][k] = (D[i] >> (8 * k)) & 0xFFu;
        uint32_t t = T[i] & 0xFFFFFFu;
        sig[i] = SYN0[t & 0xFF] ^ SYN1[(t >> 8) & 0xFF] ^ SYN2[(t >> 16) & 0xFF];
    }
    for (uint32_t a0 = 0; a0 < 256; a0++) {
        uint32_t y0[2], du2[2];
        for (int i = 0; i < 2; i++) {
            y0[i] = SBOXC[a0] ^ SBOXC[a0 ^ dl[i][0]];
            du2[i] = dl[i][1] ^ LM1[y0[i]];
        }
        uint32_t b1 = XB1[a0], b2 = XB2[a0], b3 = XB3[a0];
        for (uint32_t a1 = 0; a1 < 256; a1++) {
            uint32_t u2 = a1 ^ b1, s2 = SBOXC[u2];
            uint32_t cA[2], y1[2], y2s[2];
            int ok = 1;
            for (int i = 0; i < 2; i++) {
                y1[i] = s2 ^ SBOXC[u2 ^ du2[i]];
                uint32_t ca = (y0[i] ^ LM3[y1[i]]) | (y1[i] << 8);
                uint32_t tau = sig[i] ^ SYN0[ca & 0xFF] ^ SYN1[(ca >> 8) & 0xFF]
                              ^ SYN2[(ca >> 16) & 0xFF];
                uint32_t s = SOLVE[tau];
                if (s == NOSOLVE) { ok = 0; break; }
                cA[i] = ca; y2s[i] = s;
            }
            if (!ok) continue;
            uint32_t du3[2];
            for (int i = 0; i < 2; i++) du3[i] = dl[i][2] ^ LM2[y0[i]] ^ LM1[y1[i]];
            uint32_t pre = DDTP[du3[0]][y2s[0]];
            for (int k = 0; k < 4; k++) {
                uint32_t u3 = (pre >> (8 * k)) & 0xFFu;
                if (u3 == 0xFFu) break;
                uint32_t s3 = SBOXC[u3];
                if ((s3 ^ SBOXC[u3 ^ du3[1]]) != y2s[1]) continue;
                uint32_t a2 = u3 ^ b2 ^ XB1[u2], base3 = b3 ^ XB2[u2] ^ XB1[u3];
                for (int i = 0; i < 2; i++) {
                    uint32_t R = ((cA[i] ^ PSI24[y2s[i]]) ^ T[i]) & 0xFFFFFFu;
                    uint32_t y3n = LP0[R & 0xFF] ^ LP1[(R >> 8) & 0xFF] ^ LP2[(R >> 16) & 0xFF];
                    uint32_t du4 = dl[i][3] ^ LM3[y0[i]] ^ LM2[y1[i]] ^ LM1[y2s[i]];
                    uint32_t p4 = DDTP[du4][y3n];
                    for (int m = 0; m < 4; m++) {
                        uint32_t u4 = (p4 >> (8 * m)) & 0xFFu;
                        if (u4 == 0xFFu) break;
                        uint32_t v = a0 | (a1 << 8) | (a2 << 16) | ((u4 ^ base3) << 24);
                        uint32_t x = Xf(v);
                        if (!inH(x ^ Xf(v ^ D[0]) ^ T[0])) continue;
                        if (!inH(x ^ Xf(v ^ D[1]) ^ T[1])) continue;
                        *pv = v;
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

static void do_inject(int start) {
    int t = start;
    uint32_t pre0[3], Pj[3], T[2];
    for (int j = 0; j < 3; j++) {
        pre0[j] = SA[0][j][t - 1 - KK[j] + OFF] ^ CC[j];
        Pj[j]   = SA[0][j][t - 3 + OFF] ^ SA[0][j][t - 2 - KK[j] + OFF];
    }
    T[0] = SA[0][0][t - KK[0] + OFF] ^ SA[0][1][t - KK[1] + OFF] ^ CC[0] ^ CC[1];
    T[1] = SA[0][0][t - KK[0] + OFF] ^ SA[0][2][t - KK[2] + OFF] ^ CC[0] ^ CC[2];

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    uint32_t v = 0, m2 = 0, base0 = 0;
    uint64_t inst = 0;
    for (; inst < (1ULL << 24); inst++) {
        m2 = (uint32_t)(inst * 0x9e3779b1u + 0x7f4a7c15u);
        uint32_t base[3];
        for (int j = 0; j < 3; j++) base[j] = Xf(Pj[j] ^ m2 ^ CC[j]) ^ pre0[j];
        uint32_t D[2] = { base[0] ^ base[1], base[0] ^ base[2] };
        if (inj_scan(D, T, &v)) { base0 = base[0]; break; }
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fprintf(stderr, "  injection at %d: instance %llu, %.3f s\n", start,
            (unsigned long long)inst,
            (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec));
    if (inst >= (1ULL << 24)) { fprintf(stderr, "injection search failed\n"); exit(1); }

    setword(t - 2, m2);
    setword(t - 1, v ^ base0);
    rerun(t - 2, t - 1);
    if (!inH(gmm(0, 1, t)) || !inH(gmm(0, 2, t))) { fprintf(stderr, "inject verify failed\n"); exit(1); }
    setword(t, Ai(0, t) ^ randF());
    rerun(t, t);
    for (int k = 1; k < 3; k++)
        if (!inF(Ai(k, t) ^ WD[0][t])) { fprintf(stderr, "inject a_i not in F\n"); exit(1); }
    /* the words of steps t-2..t changed, so every state after t is stale */
    rerun(t + 1, start + PT + 8 < NST - 1 ? start + PT + 8 : NST - 1);
    chk(start, t);
}

/* ------------------------------------------------------------- conditions */
/* A w = 2 condition at step t with active streams (i,j) is exactly the injection
 * equation again.  With a = Ai(i,t-1) ^ m[t-1] the X-input of stream i at t-1,
 *
 *      gmm(i,j,t) = s_i[t-1] ^ s_j[t-1] ^ s_i[t-K_i] ^ s_j[t-K_j] ^ c_i ^ c_j
 *                 = X(a) ^ X(a ^ delta) ^ Tc,
 *      delta = Ai(i,t-1) ^ Ai(j,t-1),   Tc = s_i[t-K_i] ^ s_j[t-K_j] ^ c_i ^ c_j,
 *
 * so the 19-bit test "gmm in H" is the byte-cascade condition
 * low24(X(a) ^ X(a^delta) ^ Tc) in HH and is solved by the same filter
 * M c = M Tc instead of by re-randomising the window.
 *
 * The only difference from the injection is the domain of a.  Every stream that
 * is active at t-1 constrains a to a coset of F (13 bits), and 13 < 19, so:
 *   - if no stream is active at t-1, a is free in 2^32, 2^13 of which work: the
 *     two low bytes are enumerated and SOLVE[] pins the third, a handful of steps;
 *   - otherwise a must stay in the coset of the first active stream at t-1 (the
 *     other cosets coincide because the condition met there is a relation between
 *     the A values).  The coset holds 2^13 values against a 19-bit test, so the
 *     word at t-2 and everything the window allows before it is redrawn until the
 *     coset meets the condition: 2^6 attempts on average, each costing 32 syndrome
 *     tests on the 32 admissible low-24 patterns rather than 2^13 hash evaluations.
 * Both branches are deterministic scans of an explicit set, so the only
 * randomness left is which instance of the window is drawn. */
static int finish_a(uint32_t a0, uint32_t a1, uint32_t a2, uint32_t u2, uint32_t u3,
                    uint32_t d, uint32_t Tc, uint32_t y0, uint32_t y1, uint32_t y2,
                    uint32_t *pa) {
    uint32_t c = (y0 ^ LM3[y1] ^ LM2[y2]) | ((y1 ^ LM3[y2]) << 8) | (y2 << 16);
    uint32_t Rr = (c ^ Tc) & 0xFFFFFFu;
    uint32_t y3 = LP0[Rr & 0xFF] ^ LP1[(Rr >> 8) & 0xFF] ^ LP2[(Rr >> 16) & 0xFF];
    uint32_t du4 = ((d >> 24) & 0xFFu) ^ LM3[y0] ^ LM2[y1] ^ LM1[y2];
    uint32_t p4 = DDTP[du4][y3];
    uint32_t base3 = XB3[a0] ^ XB2[u2] ^ XB1[u3];
    uint32_t low = a0 | (a1 << 8) | (a2 << 16);
    for (int m = 0; m < 4; m++) {
        uint32_t u4 = (p4 >> (8 * m)) & 0xFFu;
        if (u4 == 0xFFu) break;
        uint32_t a = low | ((u4 ^ base3) << 24);
        if (!inH(Xf(a) ^ Xf(a ^ d) ^ Tc)) continue;
        *pa = a;
        return 1;
    }
    return 0;
}

/* a restricted to the coset base24 ^ F (32 admissible low-24 patterns) */
static int scan_coset(uint32_t d, uint32_t Tc, uint32_t base24, uint32_t *pa) {
    uint32_t sig = SYN0[Tc & 0xFF] ^ SYN1[(Tc >> 8) & 0xFF] ^ SYN2[(Tc >> 16) & 0xFF];
    uint32_t d0 = d & 0xFF, d1 = (d >> 8) & 0xFF, d2 = (d >> 16) & 0xFF;
    for (int q = 0; q < NS32; q++) {
        uint32_t low = (base24 ^ S32[q]) & 0xFFFFFFu;
        uint32_t a0 = low & 0xFF, a1 = (low >> 8) & 0xFF, a2 = (low >> 16) & 0xFF;
        uint32_t y0 = SBOXC[a0] ^ SBOXC[a0 ^ d0];
        uint32_t u2 = a1 ^ XB1[a0];
        uint32_t y1 = SBOXC[u2] ^ SBOXC[u2 ^ (d1 ^ LM1[y0])];
        uint32_t u3 = a2 ^ XB2[a0] ^ XB1[u2];
        uint32_t y2 = SBOXC[u3] ^ SBOXC[u3 ^ (d2 ^ LM2[y0] ^ LM1[y1])];
        uint32_t c = (y0 ^ LM3[y1] ^ LM2[y2]) | ((y1 ^ LM3[y2]) << 8) | (y2 << 16);
        if ((SYN0[c & 0xFF] ^ SYN1[(c >> 8) & 0xFF] ^ SYN2[(c >> 16) & 0xFF]) != sig) continue;
        if (finish_a(a0, a1, a2, u2, u3, d, Tc, y0, y1, y2, pa)) return 1;
    }
    return 0;
}

/* a free in all of 2^32.  Exhaustive, but not always satisfiable: when the low 24
 * bits of d vanish the cascade collapses to y0 = y1 = y2 = 0, c = 0, and the
 * condition degenerates to "Tc in U", which usually fails.  The caller then draws
 * a fresh window. */
static int scan_free(uint32_t d, uint32_t Tc, uint32_t *pa) {
    uint32_t sig = SYN0[Tc & 0xFF] ^ SYN1[(Tc >> 8) & 0xFF] ^ SYN2[(Tc >> 16) & 0xFF];
    uint32_t d0 = d & 0xFF, d1 = (d >> 8) & 0xFF, d2 = (d >> 16) & 0xFF;
    for (uint32_t a0 = 0; a0 < 256; a0++) {
        uint32_t y0 = SBOXC[a0] ^ SBOXC[a0 ^ d0];
        uint32_t du2 = d1 ^ LM1[y0];
        uint32_t b1 = XB1[a0], b2 = XB2[a0];
        for (uint32_t a1 = 0; a1 < 256; a1++) {
            uint32_t u2 = a1 ^ b1;
            uint32_t y1 = SBOXC[u2] ^ SBOXC[u2 ^ du2];
            uint32_t ca = (y0 ^ LM3[y1]) | (y1 << 8);
            uint32_t tau = sig ^ SYN0[ca & 0xFF] ^ SYN1[(ca >> 8) & 0xFF];
            uint32_t y2 = SOLVE[tau];
            if (y2 == NOSOLVE) continue;
            uint32_t du3 = d2 ^ LM2[y0] ^ LM1[y1];
            uint32_t p3 = DDTP[du3][y2];
            for (int k = 0; k < 4; k++) {
                uint32_t u3 = (p3 >> (8 * k)) & 0xFFu;
                if (u3 == 0xFFu) break;
                uint32_t a2 = u3 ^ b2 ^ XB1[u2];
                if (finish_a(a0, a1, a2, u2, u3, d, Tc, y0, y1, y2, pa)) return 1;
            }
        }
    }
    return 0;
}

/* ---------------------------------------------------------------- driver */
static void solve_exc(int start) {
    do_inject(start);
    int wstart = start + PLAN[NPLAN - 1][0];
    for (int j = 0; j < NPLAN; j++) {
        int rel = PLAN[j][0], wt = PLAN[j][1], ctr = PLAN[j][5];
        int act[3] = {PLAN[j][2], PLAN[j][3], PLAN[j][4]};
        int t = start + rel;
        /* wfull: every step since the previous condition, because each of them may
         * carry a single active stream whose word has to be placed inside F at all.
         * umin = wfull-1 is the previous condition step, still redrawable because its
         * own condition is a relation between A values that does not involve m. */
        int wfull = start + PLAN[j][6], umin = start + PLAN[j][6] - 1;
        if (ctr || umin > t - 2) {   /* guaranteed by the asserts in gen_data3.py */
            fprintf(stderr, "unsupported plan entry at rel %d\n", rel); exit(1);
        }
        int i = act[0], jj = act[1];
        int rprev = t - 1 - start;
        int kprev = wof(rprev) ? (BP[0][rprev] ? 0 : (BP[1][rprev] ? 1 : 2)) : -1;
        uint32_t a = 0;
        long long inst = 0;
        int ok = 0;
        struct timespec c0, c1;
        replay(start, wfull, t - 1);
        clock_gettime(CLOCK_MONOTONIC, &c0);
        for (; inst < (1LL << 24); inst++) {
            if (inst) replay(start, umin, t - 2);   /* fresh instance of the window */
            uint32_t Ai1 = Ai(i, t - 1);
            uint32_t d = Ai1 ^ Ai(jj, t - 1);
            uint32_t Tc = SA[0][i][t - KK[i] + OFF] ^ SA[0][jj][t - KK[jj] + OFF]
                        ^ CC[i] ^ CC[jj];
            /* both branches are exhaustive scans of their candidate set, but either
             * can come up empty (a coset is 6 bits short of the condition, and a
             * degenerate d collapses the cascade), so both redraw on failure */
            ok = (kprev < 0) ? scan_free(d, Tc, &a)
                             : scan_coset(d, Tc, (Ai(kprev, t - 1) ^ Ai1) & 0xFFFFFFu, &a);
            if (ok) break;
        }
        clock_gettime(CLOCK_MONOTONIC, &c1);
        fprintf(stderr, "  cond rel %2d: %lld instances (%s), %.4f s\n", rel, inst + 1,
                kprev < 0 ? "free" : "coset",
                (c1.tv_sec - c0.tv_sec) + 1e-9 * (c1.tv_nsec - c0.tv_nsec));
        if (!ok) { fprintf(stderr, "condition search failed at rel %d\n", rel); exit(1); }
        setword(t - 1, Ai(i, t - 1) ^ a);
        rerun(t - 1, t - 1);
        for (int k = 1; k < wt; k++) if (!inH(gmm(act[0], act[k], t))) {
            fprintf(stderr, "cond verify failed at rel %d\n", rel); exit(1);
        }
        if (!ctr) setword(t, Ai(act[0], t) ^ randF());
        rerun(t, t);
        for (int k = 0; k < wt; k++) if (!inF(Ai(act[k], t) ^ WD[0][t])) {
            fprintf(stderr, "cond a_i not in F at rel %d\n", rel); exit(1);
        }
        /* m[t-1] and m[t] changed: refresh the states that depend on them */
        rerun(t + 1, start + PT + 8 < NST - 1 ? start + PT + 8 : NST - 1);
        chk(start, t);
    }
    /* the tail of the excursion still has single-stream steps to satisfy */
    replay(start, wstart, start + PT + 8 < NST - 1 ? start + PT + 8 : NST - 1);
    chkall(start);
}

int main(int argc, char **argv) {
    if (argc > 1) rs = strtoull(argv[1], 0, 0) * 2654435761u + 88172645463325252ULL;

    bitF = calloc(1 << 19, 4);
    bitH = calloc(1 << 19, 4);
    for (int i = 0; i < NS32; i++) bitF[S32[i] >> 5] |= 1u << (S32[i] & 31);
    for (int i = 0; i < NH; i++) bitH[HH[i] >> 5] |= 1u << (HH[i] & 31);
    memset(DDTP, 0xFF, sizeof DDTP);
    for (int d = 0; d < 256; d++)
        for (int x = 0; x < 256; x++) {
            int y = SBOXC[x] ^ SBOXC[x ^ d];
            uint32_t *slot = &DDTP[d][y];
            for (int k = 0; k < 4; k++)
                if (((*slot >> (8 * k)) & 0xFFu) == 0xFFu) {
                    *slot = (*slot & ~(0xFFu << (8 * k))) | ((uint32_t)x << (8 * k));
                    break;
                }
        }

    for (int m = 0; m < 2; m++)
        for (int s = 0; s < 3; s++)
            for (int i = 0; i < OFF; i++) SA[m][s][i] = 0;
    for (int m = 0; m < 2; m++) SA[m][1][-7 + OFF] = 256u;

    int st[2] = {EXC1, EXC2};
    for (int t = 0; t < NST; t++) {
        if (isctr(t)) { DM[t] = 0; WD[0][t] = WD[1][t] = (uint32_t)(t / 17 + 1); continue; }
        uint32_t d = 0;
        for (int ee = 0; ee < 2; ee++) {
            int r = t - st[ee];
            if (r >= 0 && r < PT && EV[r]) d = DIFF;
        }
        DM[t] = d;
        WD[0][t] = rnd32();
        WD[1][t] = WD[0][t] ^ d;
    }
    rerun(0, NST - 1);

    struct timespec s0, s1;
    clock_gettime(CLOCK_MONOTONIC, &s0);
    for (int e = 0; e < 2; e++) {
        fprintf(stderr, "excursion %d at step %d\n", e, st[e]);
        solve_exc(st[e]);
        rerun(0, NST - 1);
    }
    clock_gettime(CLOCK_MONOTONIC, &s1);
    fprintf(stderr, "search time %.3f s\n",
            (s1.tv_sec - s0.tv_sec) + 1e-9 * (s1.tv_nsec - s0.tv_nsec));

    int bad = 0;
    for (int s = 0; s < 3; s++)
        for (int t = NST - 40; t < NST; t++)
            if (SA[0][s][t + OFF] != SA[1][s][t + OFF]) bad++;
    uint32_t p2[32] = {0}, p3[32] = {0};
    for (int t = 0; t < NST; t++) {
        p2[t % 32] ^= SA[0][1][t + OFF] ^ SA[1][1][t + OFF];
        p3[t % 32] ^= SA[0][2][t + OFF] ^ SA[1][2][t + OFF];
    }
    for (int i = 0; i < 32; i++) if (p2[i] || p3[i]) bad++;
    fprintf(stderr, "state+pool mismatches: %d\n", bad);
    if (bad) return 1;

    unsigned char *out = malloc(NB * 64);
    for (int e = 0; e < 2; e++) {
        for (int blk = 0; blk < NB; blk++)
            for (int j = 0; j < 16; j++) {
                uint32_t w = WD[e][blk * 17 + j];
                unsigned char *q = out + blk * 64 + 4 * j;
                q[0] = w; q[1] = w >> 8; q[2] = w >> 16; q[3] = w >> 24;
            }
        for (int i = 0; i < NB * 64; i++) printf("%02x", out[i]);
        printf("\n");
    }
    return 0;
}
