/* =========================================================================
 * joc2024_colm0.c -- Direct, line-by-line C implementation of the 2024
 * Journal of Cryptology specification of COLM_0:
 *   E. Andreeva, A. Bogdanov, N. Datta, A. Luykx, B. Mennink, M. Nandi,
 *   E. Tischhauser, K. Yasuda, "The COLM Authenticated Encryption Scheme",
 *   Journal of Cryptology 37:15 (2024), Section 3.1.1, Figure 2 (COLM_0.enc,
 *   proc_AD, generate_mask_enc_0) and Figure 3 (COLM_0.dec, generate_mask_dec_0).
 *
 * Supports both recommended blockcipher instantiations (Section 3.1.3, p. 10
 * and Section 5.1.2, p. 27):
 *   -DCOLM_KEY_BITS=128  ->  AES-COLM_0-128 (128-bit key, AES-128)
 *   -DCOLM_KEY_BITS=256  ->  AES-COLM_0-256 (256-bit key, AES-256)
 * ========================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "crypto_aead.h"
#include <openssl/aes.h>

#ifndef COLM_KEY_BITS
#define COLM_KEY_BITS 128
#endif

typedef unsigned char u8;
typedef unsigned long long u64;
typedef unsigned char block[16];

/* GF(2^128) modulo f(x) = x^128 + x^7 + x^2 + x + 1 (Section 2, p. 6; Section 3.2, p. 13) */
static void xor_blk(block d, const block s1, const block s2) {
    for (int i = 0; i < 16; i++) d[i] = s1[i] ^ s2[i];
}
static void copy_blk(block out, const block in) {
    memcpy(out, in, 16);
}
static void mult_2(block b, const block s) {
    block d;
    u8 msb = s[0] >> 7;
    for (int i = 0; i < 15; i++) d[i] = (u8)((s[i] << 1) | (s[i + 1] >> 7));
    d[15] = (u8)(s[15] << 1);
    if (msb) d[15] ^= 0x87;
    copy_blk(b, d);
}
static void mult_3(block d, const block s) {
    block t;
    mult_2(t, s);
    xor_blk(d, t, s);
}
static void mult_7(block d, const block s) {
    block t;
    mult_2(t, s);
    mult_3(t, t);
    xor_blk(d, t, s);
}

/* ozp(X): one-zero padding (Section 2, p. 6) */
static void ozp_load(block dst, const u8 *src, u64 len) {
    if (len >= 16) {
        memcpy(dst, src, 16);
    } else {
        memcpy(dst, src, (size_t)len);
        dst[len] = 0x80;
        for (u64 i = len + 1; i < 16; i++) dst[i] = 0x00;
    }
}

/* Linear mixing rho(x, st) = (x + 3*st, x + 2*st) [Section 3.1, p. 9] */
static void rho_mix(block Y, block W_next, const block X, const block W) {
    block w2, w3;
    mult_2(w2, W);
    mult_3(w3, W);
    xor_blk(Y, X, w3);
    xor_blk(W_next, X, w2);
}

/* Inverse linear mixing mu(y, st) = (y + 3*st, y + st) [Section 3.1, p. 9] */
static void mu_mix(block X, block W_next, const block Y, const block W) {
    block w3;
    mult_3(w3, W);
    xor_blk(X, Y, w3);
    xor_blk(W_next, Y, W);
}

/* Algorithm proc_AD(A, npub, param, Delta_1) [Figure 2, p. 11] */
static void proc_AD_joc2024(
    block IV, const u8 *ad, u64 adlen,
    const u8 *npub, const u8 *param, const block L, const AES_KEY *ek)
{
    block Delta1, blk, AA, Z, Wp;
    mult_3(Delta1, L); /* Delta_1[0] <- 3 * L */

    memcpy(blk, npub, 8);
    memcpy(blk + 8, param, 8);
    xor_blk(AA, blk, Delta1);
    AES_encrypt(AA, Wp, ek); /* line 2: W'[0] <- E_K((npub || param) xor Delta_1[0]) */

    if (adlen == 0) {
        copy_blk(IV, Wp);
        return;
    }

    u64 a = (adlen + 15) / 16;
    for (u64 i = 1; i <= a - 1; i++) {
        mult_2(Delta1, Delta1); /* Delta_1[i] <- 2 * Delta_1[i-1] */
        memcpy(blk, ad + 16 * (i - 1), 16);
        xor_blk(AA, blk, Delta1);
        AES_encrypt(AA, Z, ek);
        xor_blk(Wp, Wp, Z);
    }

    u64 rem = adlen - 16 * (a - 1);
    if (rem == 16) {
        mult_2(Delta1, Delta1); /* c_1 = 2 when |A*[a]| == n */
    } else {
        mult_7(Delta1, Delta1); /* c_1 = 7 when |A*[a]| < n */
    }
    ozp_load(blk, ad + 16 * (a - 1), rem);
    xor_blk(AA, blk, Delta1);
    AES_encrypt(AA, Z, ek);
    xor_blk(IV, Wp, Z);
}

/* Algorithm COLM_0.enc(A, M, npub, param) [Figure 2, p. 11] */
int crypto_aead_encrypt(
    unsigned char *c, unsigned long long *clen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *ad, unsigned long long adlen,
    const unsigned char *nsec,
    const unsigned char *npub,
    const unsigned char *k)
{
    (void)nsec;
    u8 param[8] = {0,0,0,0,0,0,0,0};
    AES_KEY ek;
    AES_set_encrypt_key(k, COLM_KEY_BITS, &ek);

    block zero = {0}, L, W, Delta2, Delta3, Mi, Ml, MM, X, Y, CC, Ci;
    AES_encrypt(zero, L, &ek); /* line 1: L <- E_K(0) */

    copy_blk(Delta2, L);       /* Delta_2[0] <- L */
    mult_3(Delta3, L);
    mult_3(Delta3, Delta3);    /* Delta_3[0] <- 3^2 * L */

    proc_AD_joc2024(W, ad, adlen, npub, param, L, &ek); /* lines 3-4: W[0] <- IV */

    u64 ell = (mlen == 0) ? 1 : (mlen + 15) / 16;
    u64 last_len = (mlen == 0) ? 0 : (mlen - 16 * (ell - 1));
    memset(Ml, 0, 16);

    /* lines 6-11: for i = 1 ... (ell - 1) */
    for (u64 i = 1; i <= ell - 1; i++) {
        mult_2(Delta2, Delta2);
        mult_2(Delta3, Delta3);
        memcpy(Mi, m + 16 * (i - 1), 16);
        xor_blk(Ml, Ml, Mi);
        xor_blk(MM, Mi, Delta2);
        AES_encrypt(MM, X, &ek);
        rho_mix(Y, W, X, W);
        AES_encrypt(Y, CC, &ek);
        xor_blk(Ci, CC, Delta3);
        memcpy(c + 16 * (i - 1), Ci, 16);
    }

    /* line 12: M[ell] <- M[1] xor ... xor ozp(M*[ell]) */
    ozp_load(Mi, m + 16 * (ell - 1), last_len);
    xor_blk(Ml, Ml, Mi);

    /* generate_mask_enc_0 lines 11-13: c_2 = (|M*[ell]| == n) ? 7 : 7^2 */
    mult_7(Delta2, Delta2);
    mult_7(Delta3, Delta3);
    if (last_len < 16) {
        mult_7(Delta2, Delta2);
        mult_7(Delta3, Delta3);
    }

    /* lines 13-17: block ell */
    xor_blk(MM, Ml, Delta2);
    AES_encrypt(MM, X, &ek);
    rho_mix(Y, W, X, W);
    AES_encrypt(Y, CC, &ek);
    xor_blk(Ci, CC, Delta3);
    memcpy(c + 16 * (ell - 1), Ci, 16);

    /* lines 18-23: block ell + 1 (M[ell + 1] <- M[ell]) */
    mult_2(Delta2, Delta2);
    mult_2(Delta3, Delta3);
    xor_blk(MM, Ml, Delta2);
    AES_encrypt(MM, X, &ek);
    rho_mix(Y, W, X, W);
    AES_encrypt(Y, CC, &ek);
    xor_blk(Ci, CC, Delta3);

    /* line 24: return |C[1] || ... || C[ell + 1]|_{|M| + n} */
    memcpy(c + 16 * ell, Ci, (size_t)last_len);
    *clen = mlen + 16;
    return 0;
}

/* Algorithm COLM_0.dec(A, C, npub, param) [Figure 3, p. 12] */
int crypto_aead_decrypt(
    unsigned char *m, unsigned long long *mlen,
    unsigned char *nsec,
    const unsigned char *c, unsigned long long clen,
    const unsigned char *ad, unsigned long long adlen,
    const unsigned char *npub,
    const unsigned char *k)
{
    (void)nsec;
    if (clen < 16) return -1;

    u8 param[8] = {0,0,0,0,0,0,0,0};
    AES_KEY ek, dk;
    AES_set_encrypt_key(k, COLM_KEY_BITS, &ek);
    AES_set_decrypt_key(k, COLM_KEY_BITS, &dk);

    block zero = {0}, L, W, Delta2, Delta3, CC, Y, X, MM, Mi, Ml, Mstar_l, Ctag;
    AES_encrypt(zero, L, &ek); /* line 1: L <- E_K(0) */

    copy_blk(Delta2, L);
    mult_3(Delta3, L);
    mult_3(Delta3, Delta3);

    proc_AD_joc2024(W, ad, adlen, npub, param, L, &ek); /* lines 3-4: W[0] <- IV */

    u64 msg_bytes = clen - 16;
    u64 ell = (msg_bytes == 0) ? 1 : (msg_bytes + 15) / 16;
    u64 n_prime = (msg_bytes == 0) ? 0 : (msg_bytes - 16 * (ell - 1));
    memset(Mstar_l, 0, 16);

    /* lines 6-11: for i = 1 ... ell */
    for (u64 i = 1; i <= ell; i++) {
        if (i < ell) {
            mult_2(Delta2, Delta2);
            mult_2(Delta3, Delta3);
        } else {
            mult_7(Delta2, Delta2);
            mult_7(Delta3, Delta3);
            if (n_prime < 16) {
                mult_7(Delta2, Delta2);
                mult_7(Delta3, Delta3);
            }
        }
        memcpy(CC, c + 16 * (i - 1), 16);
        xor_blk(CC, CC, Delta3);
        AES_decrypt(CC, Y, &dk);
        mu_mix(X, W, Y, W);
        AES_decrypt(X, MM, &dk);
        xor_blk(Mi, MM, Delta2);

        if (i < ell) {
            memcpy(m + 16 * (i - 1), Mi, 16);
            xor_blk(Mstar_l, Mstar_l, Mi);
        } else {
            copy_blk(Ml, Mi);               /* M[ell] */
            xor_blk(Mstar_l, Mstar_l, Ml);  /* line 12: M*[ell] <- M[1] xor ... xor M[ell] */
        }
    }

    /* lines 13-18: M[ell + 1] <- M[ell]; recompute C[ell + 1] */
    mult_2(Delta2, Delta2);
    mult_2(Delta3, Delta3);
    xor_blk(MM, Ml, Delta2);
    AES_encrypt(MM, X, &ek);
    rho_mix(Y, W, X, W);
    AES_encrypt(Y, CC, &ek);
    xor_blk(Ctag, CC, Delta3);

    /* lines 19-25: tag & padding verification */
    if (n_prime == 16) {
        if (memcmp(Ctag, c + 16 * ell, 16) != 0) return -1;
        memcpy(m + 16 * (ell - 1), Mstar_l, 16);
    } else {
        if (memcmp(Ctag, c + 16 * ell, (size_t)n_prime) != 0) return -1;
        if (Mstar_l[n_prime] != 0x80) return -1;
        for (u64 j = n_prime + 1; j < 16; j++) {
            if (Mstar_l[j] != 0x00) return -1;
        }
        memcpy(m + 16 * (ell - 1), Mstar_l, (size_t)n_prime);
    }
    *mlen = msg_bytes;
    return 0;
}
