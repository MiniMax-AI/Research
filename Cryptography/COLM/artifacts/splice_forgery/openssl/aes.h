/* Self-contained FIPS-197 AES-128 and AES-256 implementation providing the
 * OpenSSL symbols used by COLM's reference implementations:
 *     AES_KEY, AES_set_encrypt_key, AES_set_decrypt_key, AES_encrypt, AES_decrypt
 * Validated against FIPS-197 Appendix C.1 (AES-128) and Appendix C.3 (AES-256).
 */
#ifndef TINY_AES_SHIM_H
#define TINY_AES_SHIM_H
#include <stdint.h>
#include <string.h>

typedef struct { uint32_t rk[60]; int rounds; } AES_KEY;

static const uint8_t _sbox[256] = {
0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16};
static uint8_t _rsbox[256];
static int _rsbox_ready = 0;
static const uint8_t _Rcon[11] = {0x8d,0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1b,0x36};

static inline uint8_t _xtime(uint8_t x){ return (uint8_t)((x<<1) ^ (((x>>7)&1)*0x1b)); }
static inline uint8_t _gmul(uint8_t a, uint8_t b){
  uint8_t p=0; for(int i=0;i<8;i++){ if(b&1) p^=a; uint8_t hi=a&0x80; a<<=1; if(hi) a^=0x1b; b>>=1; } return p; }

static void _init_rsbox(void){ if(_rsbox_ready) return; for(int i=0;i<256;i++) _rsbox[_sbox[i]]=(uint8_t)i; _rsbox_ready=1; }

static void _key_expansion(const uint8_t *key, int bits, AES_KEY *k){
  int Nk = (bits == 256) ? 8 : 4;
  int Nr = (bits == 256) ? 14 : 10;
  int Nw = 4 * (Nr + 1);
  k->rounds = Nr;
  uint8_t *rk = (uint8_t*)k->rk;
  memcpy(rk, key, 4 * Nk);
  for(int i = Nk; i < Nw; i++){
    uint8_t t[4]; memcpy(t, rk + 4*(i-1), 4);
    if(i % Nk == 0){
      uint8_t tmp=t[0]; t[0]=_sbox[t[1]]; t[1]=_sbox[t[2]]; t[2]=_sbox[t[3]]; t[3]=_sbox[tmp];
      t[0] ^= _Rcon[i / Nk];
    } else if(Nk == 8 && (i % Nk == 4)){
      t[0]=_sbox[t[0]]; t[1]=_sbox[t[1]]; t[2]=_sbox[t[2]]; t[3]=_sbox[t[3]];
    }
    for(int j=0;j<4;j++) rk[4*i+j] = rk[4*(i-Nk)+j] ^ t[j];
  }
}

static void _add_round_key(uint8_t *s, const uint8_t *rk, int r){ for(int i=0;i<16;i++) s[i]^=rk[16*r+i]; }

static void _sub_bytes(uint8_t *s, const uint8_t *box){ for(int i=0;i<16;i++) s[i]=box[s[i]]; }

static void _shift_rows(uint8_t *s){
  uint8_t t;
  t=s[1]; s[1]=s[5]; s[5]=s[9]; s[9]=s[13]; s[13]=t;
  t=s[2]; s[2]=s[10]; s[10]=t; t=s[6]; s[6]=s[14]; s[14]=t;
  t=s[15]; s[15]=s[11]; s[11]=s[7]; s[7]=s[3]; s[3]=t;
}
static void _inv_shift_rows(uint8_t *s){
  uint8_t t;
  t=s[13]; s[13]=s[9]; s[9]=s[5]; s[5]=s[1]; s[1]=t;
  t=s[2]; s[2]=s[10]; s[10]=t; t=s[6]; s[6]=s[14]; s[14]=t;
  t=s[3]; s[3]=s[7]; s[7]=s[11]; s[11]=s[15]; s[15]=t;
}
static void _mix_columns(uint8_t *s){
  for(int c=0;c<4;c++){
    uint8_t a0=s[4*c],a1=s[4*c+1],a2=s[4*c+2],a3=s[4*c+3];
    uint8_t all=a0^a1^a2^a3;
    s[4*c+0]=a0^all^_xtime(a0^a1);
    s[4*c+1]=a1^all^_xtime(a1^a2);
    s[4*c+2]=a2^all^_xtime(a2^a3);
    s[4*c+3]=a3^all^_xtime(a3^a0);
  }
}
static void _inv_mix_columns(uint8_t *s){
  for(int c=0;c<4;c++){
    uint8_t a0=s[4*c],a1=s[4*c+1],a2=s[4*c+2],a3=s[4*c+3];
    s[4*c+0]=_gmul(a0,14)^_gmul(a1,11)^_gmul(a2,13)^_gmul(a3,9);
    s[4*c+1]=_gmul(a0,9)^_gmul(a1,14)^_gmul(a2,11)^_gmul(a3,13);
    s[4*c+2]=_gmul(a0,13)^_gmul(a1,9)^_gmul(a2,14)^_gmul(a3,11);
    s[4*c+3]=_gmul(a0,11)^_gmul(a1,13)^_gmul(a2,9)^_gmul(a3,14);
  }
}

static void _cipher(const uint8_t *in, uint8_t *out, const AES_KEY *k, int decrypt){
  uint8_t s[16]; memcpy(s,in,16);
  const uint8_t *rk=(const uint8_t*)k->rk;
  int Nr = (k->rounds == 14) ? 14 : 10;
  if(!decrypt){
    _add_round_key(s,rk,0);
    for(int r=1;r<Nr;r++){ _sub_bytes(s,_sbox); _shift_rows(s); _mix_columns(s); _add_round_key(s,rk,r); }
    _sub_bytes(s,_sbox); _shift_rows(s); _add_round_key(s,rk,Nr);
  } else {
    _init_rsbox();
    _add_round_key(s,rk,Nr);
    for(int r=Nr-1;r>=1;r--){ _inv_shift_rows(s); _sub_bytes(s,_rsbox); _add_round_key(s,rk,r); _inv_mix_columns(s); }
    _inv_shift_rows(s); _sub_bytes(s,_rsbox); _add_round_key(s,rk,0);
  }
  memcpy(out,s,16);
}

static inline int AES_set_encrypt_key(const uint8_t *key, int bits, AES_KEY *k){ _key_expansion(key,bits,k); return 0; }
static inline int AES_set_decrypt_key(const uint8_t *key, int bits, AES_KEY *k){ _key_expansion(key,bits,k); return 0; }
static inline void AES_encrypt(const uint8_t *in, uint8_t *out, const AES_KEY *k){ _cipher(in,out,k,0); }
static inline void AES_decrypt(const uint8_t *in, uint8_t *out, const AES_KEY *k){ _cipher(in,out,k,1); }
#endif
