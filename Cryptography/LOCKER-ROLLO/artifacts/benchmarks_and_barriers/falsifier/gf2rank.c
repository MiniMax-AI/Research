/* gf2rank.c -- exact GF(2) rank of a dense bit-matrix via M4RI.
 * input file layout (little-endian):
 *   uint64 nrows, uint64 ncols
 *   nrows * ceil(ncols/64) uint64 words; word j of row i holds bits 64j..64j+63 (LSB = bit 64j)
 * prints: "<rank> <nrows> <ncols>"
 * build: gcc -O3 -o gf2rank gf2rank.c $(pkg-config --cflags --libs m4ri)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <m4ri/m4ri.h>

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: gf2rank file\n"); return 2; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("open"); return 2; }
    uint64_t nr, nc;
    if (fread(&nr, 8, 1, f) != 1 || fread(&nc, 8, 1, f) != 1) return 2;
    uint64_t wpr = (nc + 63) / 64;
    mzd_t *M = mzd_init((rci_t)nr, (rci_t)nc);
    uint64_t *buf = (uint64_t *)malloc(wpr * 8);
    for (uint64_t i = 0; i < nr; ++i) {
        if (fread(buf, 8, wpr, f) != wpr) { fprintf(stderr, "short read\n"); return 2; }
        for (uint64_t j = 0; j < wpr; ++j) {
            uint64_t w = buf[j];
            while (w) {
                uint64_t b = w & (~w + 1);
                int k = __builtin_ctzll(w);
                mzd_write_bit(M, (rci_t)i, (rci_t)(64 * j + k), 1);
                w ^= b;
            }
        }
    }
    fclose(f);
    rci_t rk = mzd_echelonize(M, 1);
    printf("%u %llu %llu\n", (unsigned)rk, (unsigned long long)nr, (unsigned long long)nc);
    mzd_free(M); free(buf);
    return 0;
}
