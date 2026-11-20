/*
 * xdiff6.h -- closed-form differential algebra of the 4-round byte SPN `xbox'.
 *
 * xbox(w) is  w ^= XB[w&0xff]; w = ROR8(w)  applied four times.  Writing the
 * state as four bytes b0..b3 (b0 = w&0xff) and i_k for the sbox index of round
 * k, one round is
 *      b0' = b1 ^ XB1[i],  b1' = b2 ^ XB2[i],  b2' = b3 ^ XB3[i],  b3' = sbox[i]
 * with i = b0  (because XB0[i] = sbox[i]^i).  Since
 *      XB[i]^XB[i^d] = ( ds^d, A1(ds), A2(ds), A3(ds) )   with ds = sbox[i]^sbox[i^d]
 * (verified exhaustively), the *difference* trajectory obeys, with e_k the index
 * difference of round k and ds_k its sbox difference,
 *      e_{k+1} = d1_k ^ A1(ds_k),  d1_{k+1} = d2_k ^ A2(ds_k),
 *      d2_{k+1} = d3_k ^ A3(ds_k),  d3_{k+1} = ds_k.
 * Chaining backwards from the target output difference TG determines ds_4, ds_3,
 * ds_2, ds_1 uniquely, and chaining forwards from the input difference G then
 * determines e_1..e_4.  Hence
 *
 *      #{ w : Phi_w(G) = TG }  =  prod_k DDT[e_k][ds_k]            (closed form)
 *
 * and, more importantly, every such w can be *constructed* in O(1) from a
 * choice of round indices i_k in S_k = { i : sbox[i]^sbox[i^e_k] = ds_k }:
 *      b0 = i1,  b1 = i2 ^ XB1[i1],  b2 = i3 ^ XB1[i2] ^ XB2[i1],
 *      b3 = i4 ^ XB1[i3] ^ XB2[i2] ^ XB3[i1].
 * Dually, for a *fixed* word w the difference e with Phi_w(e) = TG is unique and
 * is obtained by the same backward chain (each S_k degenerates to one element
 * because sbox is a bijection).
 *
 * Nothing here searches: T = {Y : Phi_Y(BETA)=BETA}, the injection alignment and
 * every steering step of the walk are all produced by these two constructors.
 */
#ifndef XDIFF6_H
#define XDIFF6_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

#define B0(w) ((unsigned char)((w)     ))
#define B1(w) ((unsigned char)((w) >>  8))
#define B2(w) ((unsigned char)((w) >> 16))
#define B3(w) ((unsigned char)((w) >> 24))
#define XB0(i) ((unsigned char)(XB[i]     ))
#define XB1(i) ((unsigned char)(XB[i] >>  8))
#define XB2(i) ((unsigned char)(XB[i] >> 16))
#define XB3(i) ((unsigned char)(XB[i] >> 24))

static unsigned short DDT[256][256];
static unsigned char SBINV[256];
static unsigned char A1t[256], A2t[256], A3t[256];
static unsigned char Nz[256];        /* Nz[ds] = #{e : DDT[e][ds] > 0}          */
static unsigned char Lz[256][256];   /* Lz[ds][.] = those e                     */
static int xd_ready = 0;

static void xd_init(void){
    if (xd_ready) return;
    xb_init();
    memset(DDT, 0, sizeof DDT);
    for (int i = 0; i < 256; i++)
        for (int d = 0; d < 256; d++) DDT[d][sbox[i] ^ sbox[i ^ d]]++;
    for (int i = 0; i < 256; i++) SBINV[sbox[i]] = (unsigned char)i;
    memset(A1t, 0, sizeof A1t); memset(A2t, 0, sizeof A2t); memset(A3t, 0, sizeof A3t);
    unsigned char seen[256]; memset(seen, 0, sizeof seen);
    for (int i = 0; i < 256; i++) for (int d = 0; d < 256; d++){
        unsigned ds = sbox[i] ^ sbox[i ^ d];
        unsigned char a1 = XB1(i) ^ XB1(i ^ d), a2 = XB2(i) ^ XB2(i ^ d),
                      a3 = XB3(i) ^ XB3(i ^ d);
        if (!seen[ds]) { A1t[ds] = a1; A2t[ds] = a2; A3t[ds] = a3; seen[ds] = 1; }
        else if (A1t[ds] != a1 || A2t[ds] != a2 || A3t[ds] != a3){
            fprintf(stderr, "xd_init: A_j not a function of ds (%u)\n", ds); exit(1);
        }
    }
    for (int ds = 0; ds < 256; ds++){
        int n = 0;
        for (int e = 0; e < 256; e++) if (DDT[e][ds]) Lz[ds][n++] = (unsigned char)e;
        Nz[ds] = (unsigned char)n;
    }
    xd_ready = 1;
}

/* the 4-round characteristic forced by the target output difference */
static void xd_char(uint32_t TG, unsigned char ds[4]){
    unsigned char t0 = B0(TG), t1 = B1(TG), t2 = B2(TG), t3 = B3(TG);
    ds[3] = t3;
    ds[2] = (unsigned char)(t2 ^ A3t[ds[3]]);
    ds[1] = (unsigned char)(t1 ^ A2t[ds[3]] ^ A3t[ds[2]]);
    ds[0] = (unsigned char)(t0 ^ A1t[ds[3]] ^ A2t[ds[2]] ^ A3t[ds[1]]);
}

/* forward chain: round index differences implied by the input difference */
static void xd_ek(uint32_t G, const unsigned char ds[4], unsigned char e[4]){
    e[0] = B0(G);
    e[1] = (unsigned char)(B1(G) ^ A1t[ds[0]]);
    e[2] = (unsigned char)(B2(G) ^ A2t[ds[0]] ^ A1t[ds[1]]);
    e[3] = (unsigned char)(B3(G) ^ A3t[ds[0]] ^ A2t[ds[1]] ^ A1t[ds[2]]);
}

typedef struct {
    uint32_t G, TG;
    unsigned char ds[4], e[4];
    unsigned char idx[4][256];
    int cnt[4];
    unsigned long long total;
} xchar;

static int xd_prep(uint32_t G, uint32_t TG, xchar *c){
    xd_init();
    c->G = G; c->TG = TG;
    xd_char(TG, c->ds);
    xd_ek(G, c->ds, c->e);
    unsigned long long tot = 1;
    for (int k = 0; k < 4; k++){
        int n = 0;
        for (int i = 0; i < 256; i++)
            if ((unsigned)(sbox[i] ^ sbox[i ^ c->e[k]]) == c->ds[k]) c->idx[k][n++] = (unsigned char)i;
        c->cnt[k] = n;
        if (!n){ c->total = 0; return 0; }
        tot *= (unsigned long long)n;
    }
    c->total = tot;
    return 1;
}

/* n-th word w with Phi_w(G) = TG  (0 <= n < total) */
static uint32_t xd_word(const xchar *c, unsigned long long n){
    unsigned i1 = c->idx[0][ n % (unsigned)c->cnt[0] ]; n /= (unsigned)c->cnt[0];
    unsigned i2 = c->idx[1][ n % (unsigned)c->cnt[1] ]; n /= (unsigned)c->cnt[1];
    unsigned i3 = c->idx[2][ n % (unsigned)c->cnt[2] ]; n /= (unsigned)c->cnt[2];
    unsigned i4 = c->idx[3][ n % (unsigned)c->cnt[3] ];
    unsigned b1 = i2 ^ XB1(i1);
    unsigned b2 = i3 ^ XB1(i2) ^ XB2(i1);
    unsigned b3 = i4 ^ XB1(i3) ^ XB2(i2) ^ XB3(i1);
    return (uint32_t)(i1 | (b1 << 8) | (b2 << 16) | (b3 << 24));
}

/* the round index bytes of a concrete word */
static void xd_idx(uint32_t w, unsigned char i[4]){
    for (int k = 0; k < 4; k++){ i[k] = (unsigned char)(w & 0xff); w ^= XB[i[k]]; w = ROT(w); }
}

/* the unique difference e with Phi_w(e) = TG, for a *fixed* word w */
static int xd_diff(uint32_t w, uint32_t TG, uint32_t *out){
    xd_init();
    unsigned char ds[4], iv[4], e[4];
    xd_char(TG, ds);
    xd_idx(w, iv);
    for (int k = 0; k < 4; k++)
        e[k] = (unsigned char)(iv[k] ^ SBINV[(unsigned)(sbox[iv[k]] ^ ds[k])]);
    unsigned b0 = e[0];
    unsigned b1 = e[1] ^ A1t[ds[0]];
    unsigned b2 = e[2] ^ A2t[ds[0]] ^ A1t[ds[1]];
    unsigned b3 = e[3] ^ A3t[ds[0]] ^ A2t[ds[1]] ^ A1t[ds[2]];
    uint32_t E = (uint32_t)(b0 | (b1 << 8) | (b2 << 16) | (b3 << 24));
    if ((xboxf(w) ^ xboxf(w ^ E)) != TG) return 0;
    *out = E;
    return 1;
}

/* n-th input difference G admitting at least one w with Phi_w(G) = TG */
static uint32_t xd_goodG(uint32_t TG, unsigned long long n){
    xd_init();
    unsigned char ds[4];
    xd_char(TG, ds);
    unsigned e1 = Lz[ds[0]][ n % Nz[ds[0]] ]; n /= Nz[ds[0]];
    unsigned e2 = Lz[ds[1]][ n % Nz[ds[1]] ]; n /= Nz[ds[1]];
    unsigned e3 = Lz[ds[2]][ n % Nz[ds[2]] ]; n /= Nz[ds[2]];
    unsigned e4 = Lz[ds[3]][ n % Nz[ds[3]] ];
    unsigned g0 = e1;
    unsigned g1 = e2 ^ A1t[ds[0]];
    unsigned g2 = e3 ^ A2t[ds[0]] ^ A1t[ds[1]];
    unsigned g3 = e4 ^ A3t[ds[0]] ^ A2t[ds[1]] ^ A1t[ds[2]];
    return (uint32_t)(g0 | (g1 << 8) | (g2 << 16) | (g3 << 24));
}
static unsigned long long xd_ngood(uint32_t TG){
    xd_init();
    unsigned char ds[4];
    xd_char(TG, ds);
    return (unsigned long long)Nz[ds[0]] * Nz[ds[1]] * Nz[ds[2]] * Nz[ds[3]];
}

#define PHIW(w,e) (xboxf(w) ^ xboxf((w) ^ (e)))
#endif
