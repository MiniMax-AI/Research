// Ground-truth dumper: replicates keygen/encaps internals and prints them as hex.
#include "api.h"
#include "parameters.h"
#include "rng.h"
#include "hash.h"
#include "parsing.h"
#include "kemset_types.h"
#include "ffi_field.h"
#include "ffi_vec.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);
int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);
int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);
void randombytes_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength);

static void phex(const char *label, const unsigned char *d, size_t n) {
  printf("%s: 0x", label);
  for (size_t i = 0; i < n; i++) printf("%02x", d[i]);
  printf("\n");
}

int main(int argc, char **argv) {
  unsigned char entropy[48];
  FILE *u = fopen("/dev/urandom", "rb");
  if (fread(entropy, 1, 48, u) != 48) return 1;
  fclose(u);
  randombytes_init(entropy, NULL, 256);

  ffi_field_init();
  ffi_vec_init_mulmod();

  unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
  if (crypto_kem_keypair(pk, sk)) return 1;

  secretKey sk_tmp;
  kemset_secret_key_from_string(sk_tmp, sk);
  phex("F", (const unsigned char *)nullptr, 0);
  unsigned char buf[PARAM_N * GF2MBYTES];
  ffi_vec_to_string(buf, sk_tmp.F, PARAM_D);
  phex("F", buf, PARAM_D * GF2MBYTES);
  ffi_vec_to_string(buf, sk_tmp.x, PARAM_N);
  phex("x", buf, PARAM_N * GF2MBYTES);
  ffi_vec_to_string(buf, sk_tmp.y, PARAM_N);
  phex("y", buf, PARAM_N * GF2MBYTES);
  phex("pk", pk, CRYPTO_PUBLICKEYBYTES);
  phex("sk", sk, CRYPTO_SECRETKEYBYTES);

  // ---- replicate encapsulation with full visibility ----
  publicKey pk_tmp;
  kemset_public_key_from_string(pk_tmp, pk);

  unsigned char m[CRYPTO_BYTES];
  randombytes(m, CRYPTO_BYTES);
  unsigned char seed_G[CRYPTO_BYTES];
  memcpy(seed_G, m, CRYPTO_BYTES);
  AES_XOF_struct *G = (AES_XOF_struct *)malloc(sizeof(AES_XOF_struct));
  seedexpander_init(G, seed_G, seed_G + 32, SEEDEXPANDER_MAX_LENGTH);
  unsigned char theta[SEEDEXPANDER_SEED_BYTES];
  seedexpander(G, theta, SEEDEXPANDER_SEED_BYTES);
  AES_XOF_struct *enc = (AES_XOF_struct *)malloc(sizeof(AES_XOF_struct));
  seedexpander_init(enc, theta, theta + 32, SEEDEXPANDER_MAX_LENGTH);

  ffi_vec E, polyE1, polyE2;
  ffi_vec_set_random_full_rank_using_seedexpander(E, PARAM_R, enc);
  ffi_vec_to_string(buf, E, PARAM_R);
  phex("E_raw", buf, PARAM_R * GF2MBYTES);
  ffi_vec_set_random_from_support_using_seedexpander(polyE1, PARAM_N, E, PARAM_R, enc);
  ffi_vec_set_random_from_support_using_seedexpander(polyE2, PARAM_N, E, PARAM_R, enc);
  ffi_vec_to_string(buf, polyE1, PARAM_N);
  phex("e1", buf, PARAM_N * GF2MBYTES);
  ffi_vec_to_string(buf, polyE2, PARAM_N);
  phex("e2", buf, PARAM_N * GF2MBYTES);

  ciphertext c_tmp;
  // NOTE: ffi_vec_add(o,v1,v2) zeroes o first, so passing o aliases v1 destroys
  // v1 (u becomes e1).  Compute the syndrome without aliasing so that this
  // generator emits the *specified* ciphertext u = e1 + e2*h.
  ffi_vec e2h; ffi_vec_mul(e2h, polyE2, pk_tmp.h, PARAM_N);
  ffi_vec_to_string(buf, e2h, PARAM_N);
  phex("e2h", buf, PARAM_N * GF2MBYTES);
  ffi_vec_add(c_tmp.u, e2h, polyE1, PARAM_N);
  ffi_vec_to_string(buf, c_tmp.u, PARAM_N);
  phex("u_fix", buf, PARAM_N * GF2MBYTES);

  ffi_vec_echelonize(E, PARAM_R);
  ffi_vec_to_string(buf, E, PARAM_R);
  phex("E_echelon", buf, PARAM_R * GF2MBYTES);

  unsigned char support[PARAM_R * GF2MBYTES], hashSupp[CRYPTO_BYTES];
  ffi_vec_to_string(support, E, PARAM_R);
  sha512(hashSupp, support, PARAM_R * GF2MBYTES);
  for (int i = 0; i < CRYPTO_BYTES; i++) c_tmp.v[i] = m[i] ^ hashSupp[i];
  sha512(c_tmp.d, m, CRYPTO_BYTES);
  unsigned char mc[CIPHERTEXT_BYTES], ss[CRYPTO_BYTES];
  memcpy(mc, m, SHA512_BYTES);
  ffi_vec_to_string(mc + SHA512_BYTES, c_tmp.u, PARAM_N);
  memcpy(mc + SHA512_BYTES + PARAM_N * GF2MBYTES, c_tmp.v, SHA512_BYTES);
  sha512(ss, mc, CIPHERTEXT_BYTES);
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
  kemset_ciphertext_to_string(ct, c_tmp);

  phex("m", m, CRYPTO_BYTES);
  phex("ct", ct, CRYPTO_CIPHERTEXTBYTES);
  phex("ss", ss, CRYPTO_BYTES);

  // sanity: real decapsulation
  unsigned char ss2[CRYPTO_BYTES];
  int rc = crypto_kem_dec(ss2, ct, sk);
  printf("dec_rc: %d\n", rc);
  phex("ss2", ss2, CRYPTO_BYTES);
  printf("match: %d\n", memcmp(ss, ss2, CRYPTO_BYTES) == 0);
  return 0;
}
