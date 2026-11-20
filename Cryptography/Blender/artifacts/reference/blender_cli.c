/*
 * blender_cli.c
 *
 * Offline command-line verifier linking directly against Colin Bradbury's
 * unmodified NIST SHA-3 Round 1 C reference implementation (core.c).
 *
 * Modes:
 *   -l <bits>        Digest bit length (e.g., 224, 256, 384, 512; default 256)
 *   -p <hex>         Message input as hexadecimal string
 *   -f <path>        Message input as raw binary file
 *   --sfs-iv <hex>   Semi-free-start initial state: loads custom Word[0..w-1],
 *                    Sum[0..w-1], and (carry1, carry2) in big-endian word format
 *                    (2*w*b + 2 bytes) and invokes official Update() + Final()
 *   --spec-384       Optional: initializes Word[5] = 0x8eb44a8768581511 in
 *                    Blender-384 (official core.c Init(384) initializes 5 of 6 words)
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"

static int parse_hex(const char *s, uint8_t **out_buf, size_t *out_len) {
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) s += 2;
    size_t n = strlen(s);
    if (n % 2 != 0) return -1;
    size_t blen = n / 2;
    uint8_t *buf = (uint8_t *)malloc(blen ? blen : 1);
    if (!buf) return -1;
    for (size_t i = 0; i < blen; i++) {
        unsigned int v = 0;
        if (sscanf(s + 2 * i, "%2x", &v) != 1) {
            free(buf);
            return -1;
        }
        buf[i] = (uint8_t)v;
    }
    *out_buf = buf;
    *out_len = blen;
    return 0;
}

static int read_binary_file(const char *path, uint8_t **out_buf, size_t *out_len) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return -1; }
    long sz = ftell(fp);
    if (sz < 0) { fclose(fp); return -1; }
    rewind(fp);
    uint8_t *buf = (uint8_t *)malloc(sz ? (size_t)sz : 1);
    if (!buf) { fclose(fp); return -1; }
    if (sz > 0 && fread(buf, 1, (size_t)sz, fp) != (size_t)sz) {
        free(buf);
        fclose(fp);
        return -1;
    }
    fclose(fp);
    *out_buf = buf;
    *out_len = (size_t)sz;
    return 0;
}

static uint64_t read_be(const uint8_t *p, int b) {
    uint64_t v = 0;
    for (int i = 0; i < b; i++) v = (v << 8) | p[i];
    return v;
}

int main(int argc, char **argv) {
    int hashbitlen = 256;
    int spec_384 = 0;
    const char *msg_hex = NULL;
    const char *file_path = NULL;
    const char *sfs_hex = NULL;

    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--hashbitlen") == 0) && i + 1 < argc) {
            hashbitlen = atoi(argv[++i]);
        } else if ((strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--plaintext") == 0) && i + 1 < argc) {
            msg_hex = argv[++i];
        } else if ((strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--file") == 0) && i + 1 < argc) {
            file_path = argv[++i];
        } else if (strcmp(argv[i], "--sfs-iv") == 0 && i + 1 < argc) {
            sfs_hex = argv[++i];
        } else if (strcmp(argv[i], "--spec-384") == 0) {
            spec_384 = 1;
        }
    }

    if (!msg_hex && !file_path) {
        fprintf(stderr, "Usage: %s -l <bits> (-p <hex> | -f <file>) [--sfs-iv <hex>] [--spec-384]\n", argv[0]);
        return 1;
    }

    uint8_t *msg = NULL;
    size_t msg_len = 0;
    int rc_in = file_path ? read_binary_file(file_path, &msg, &msg_len)
                          : parse_hex(msg_hex, &msg, &msg_len);
    if (rc_in != 0) {
        fprintf(stderr, "Failed to read input message\n");
        return 1;
    }

    uint8_t out[128] = {0};
    HashReturn ret;

    if (!sfs_hex && !spec_384) {
        ret = Hash(hashbitlen, msg, (DataLength)msg_len * 8ULL, out);
    } else {
        hashState st;
        ret = Init(&st, hashbitlen);
        if (ret != SUCCESS) {
            fprintf(stderr, "Init(%d) failed: %d\n", hashbitlen, ret);
            free(msg);
            return 1;
        }
        if (spec_384 && hashbitlen == 384) {
            st.Word[5].ui64 = 0x8eb44a8768581511ULL;
        }
        if (sfs_hex) {
            uint8_t *ivbuf = NULL;
            size_t ivlen = 0;
            int w = st.wordCount;
            int b = st.wordSizeBytes;
            size_t expected = (size_t)(2 * w * b + 2);
            if (parse_hex(sfs_hex, &ivbuf, &ivlen) != 0 || ivlen != expected) {
                fprintf(stderr, "Expected %zu-byte sfs-iv for Blender-%d, got %zu bytes\n",
                        expected, hashbitlen, ivlen);
                free(ivbuf);
                free(msg);
                return 1;
            }
            for (int i = 0; i < w; i++) {
                uint64_t wv = read_be(ivbuf + i * b, b);
                uint64_t sv = read_be(ivbuf + (w + i) * b, b);
                if (b == 8) {
                    st.Word[i].ui64 = wv;
                    st.Sum[i].ui64 = sv;
                } else if (b == 4) {
                    st.Word[i].ui32 = (uint32_t)wv;
                    st.Sum[i].ui32 = (uint32_t)sv;
                } else {
                    st.Word[i].ui16 = (uint16_t)wv;
                    st.Sum[i].ui16 = (uint16_t)sv;
                }
            }
            uint8_t c1 = ivbuf[2 * w * b] & 1;
            uint8_t c2 = ivbuf[2 * w * b + 1] & 1;
            st.carries = (uint8_t)((c2 << 1) | c1);
            free(ivbuf);
        }
        ret = Update(&st, msg, (DataLength)msg_len * 8ULL);
        if (ret == SUCCESS) ret = Final(&st, out);
    }

    free(msg);
    if (ret != SUCCESS) {
        fprintf(stderr, "Hash(%d) failed: %d\n", hashbitlen, ret);
        return 1;
    }

    for (int i = 0; i < hashbitlen / 8; i++) printf("%02x", out[i]);
    printf("\n");
    return 0;
}
