#ifndef MODEL_H
#define MODEL_H
#include <stdint.h>
#include <string.h>

#define STREAM1 16
#define STREAM2 7
#define STREAM3 6
#define POOLFACTOR 2
#define POOL (STREAM1*POOLFACTOR)
#define FINALUPDATES 16
#define FINALSTEPS 4

static const unsigned char sbox[256] = {
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
 0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};

static uint32_t XB[256];
static int xb_ready = 0;
static void xb_init(void){
  if (xb_ready) return;
  for (int i=0;i<256;i++){
    unsigned char w1 = sbox[i];
    unsigned char w2;
    XB[i] = (uint32_t)(w1 ^ i);
    w2 = w1 ^ (w1<<2) ^ (w1>>6) ^ (w1<<3);
    w2 = (w2<<7) | (w2>>1) ^ 0x3b;
    XB[i] |= ((uint32_t)w2)<<8;
    w2 = w1 ^ (w1<<5) ^ (w1>>3) ^ (w1>>7) ^ 0x95;
    XB[i] |= ((uint32_t)w2)<<16;
    w2 = w1 ^ (w1<<7) ^ (w1>>1) ^ (w1>>5) ^ 0x6a;
    XB[i] |= ((uint32_t)w2)<<24;
  }
  xb_ready = 1;
}

#define ROT(w) ((w<<24)|(w>>8))
static inline uint32_t xboxf(uint32_t w){
  w ^= XB[w & 0xFF]; w = ROT(w);
  w ^= XB[w & 0xFF]; w = ROT(w);
  w ^= XB[w & 0xFF]; w = ROT(w);
  w ^= XB[w & 0xFF]; w = ROT(w);
  return w;
}

typedef struct {
  uint32_t hashbitlen;
  unsigned char inputbuffer[64];
  uint32_t inputbits, inputblocks, inputslabs;
  uint32_t istream1, stream1[STREAM1];
  uint32_t istream2, stream2[STREAM2];
  uint32_t istream3, stream3[STREAM3];
  uint32_t ipool, pool2[POOL], pool3[POOL];
  /* instrumentation: step counter */
  uint64_t nsteps;
} hstate;

static int hinit(hstate *s, int hashbitlen){
  xb_init();
  if (hashbitlen<=0 || hashbitlen>32*STREAM1 || hashbitlen%32) return 2;
  memset(s,0,sizeof(*s));
  s->hashbitlen = hashbitlen;
  s->stream2[0] = hashbitlen;
  s->istream1 = STREAM1-1; s->istream2 = STREAM2-1; s->istream3 = STREAM3-1;
  s->ipool = POOL-1;
  return 0;
}

static inline void ustep(hstate *s, uint32_t input){
  uint32_t w1 = s->stream1[s->istream1];
  uint32_t w2 = s->stream2[s->istream2];
  uint32_t w3 = s->stream3[s->istream3];
  if(++s->istream1==STREAM1) s->istream1=0;
  if(++s->istream2==STREAM2) s->istream2=0;
  if(++s->istream3==STREAM3) s->istream3=0;
  w1 ^= s->stream1[s->istream1] ^ input; w1 = xboxf(w1); s->stream1[s->istream1]=w1;
  w2 ^= s->stream2[s->istream2] ^ input ^ 1u; w2 = xboxf(w2); s->stream2[s->istream2]=w2;
  w3 ^= s->stream3[s->istream3] ^ input ^ 2u; w3 = xboxf(w3); s->stream3[s->istream3]=w3;
  if(++s->ipool==POOL) s->ipool=0;
  s->pool2[s->ipool] ^= w2;
  s->pool3[s->ipool] ^= w3;
  s->nsteps++;
}

static void hupdate(hstate *s, const unsigned char *data, uint64_t nbytes){
  int i = s->inputbits/8;
  while (nbytes--) {
    s->inputbits += 8;
    s->inputbuffer[i++] = *data++;
    if (i==64){
      for (int j=0;j<16;j++){
        uint32_t in = s->inputbuffer[4*j] | (s->inputbuffer[4*j+1]<<8)
                    | (s->inputbuffer[4*j+2]<<16) | ((uint32_t)s->inputbuffer[4*j+3]<<24);
        ustep(s,in);
      }
      s->inputblocks++;
      if(!s->inputblocks) s->inputslabs++;
      ustep(s, s->inputblocks);
      i=0; s->inputbits=0;
    }
  }
}

static void finalstep(uint32_t *wi, uint32_t *pp, uint32_t *wo){
  for (int i=0;i<STREAM1;i++)
    wo[i] = (wi[i]&0x14141414u) | (wi[(i+3)%STREAM1]&0x28282828u)
          | (wi[(i+7)%STREAM1]&0x42424242u) | (wi[(i+12)%STREAM1]&0x81818181u);
  for (int i=0;i<STREAM1;i++) wo[i] = xboxf(wo[i]);
  if (pp) for (int i=0;i<STREAM1;i++) wo[i] ^= pp[i];
}

static void hfinal(hstate *s, unsigned char *hashval){
  uint32_t work1[STREAM1], work2[STREAM1], *pin=work1, *pout=work2, *tmp;
  uint32_t *pp;
  uint32_t i, pass;
  if (s->inputbits){
    i = (s->inputbits+7)/8;
    for(; i<64; i++) s->inputbuffer[i]=0;
    s->inputblocks++;
    if(!s->inputblocks) s->inputslabs++;
    for(i=0;i<16;i++){
      uint32_t in = s->inputbuffer[4*i] | (s->inputbuffer[4*i+1]<<8)
                  | (s->inputbuffer[4*i+2]<<16) | ((uint32_t)s->inputbuffer[4*i+3]<<24);
      ustep(s,in);
    }
    ustep(s, s->inputblocks);
  }
  ustep(s, s->inputbits);
  ustep(s, s->inputslabs);
  for(i=0;i<FINALUPDATES;i++) ustep(s,i);
  for(i=0;i<STREAM1;i++) work1[i]=s->stream1[i];
  pp = s->pool2;
  for(pass=0;pass<POOLFACTOR;pass++){ finalstep(pin,pp,pout); tmp=pin;pin=pout;pout=tmp; pp+=STREAM1; }
  pp = s->pool3;
  for(pass=0;pass<POOLFACTOR;pass++){ finalstep(pin,pp,pout); tmp=pin;pin=pout;pout=tmp; pp+=STREAM1; }
  for(pass=0;pass<FINALSTEPS;pass++){ finalstep(pin,NULL,pout); tmp=pin;pin=pout;pout=tmp; }
  for(i=0;i<STREAM1;i++) pin[i]^=s->stream1[i];
  pp = s->pool2;
  for(pass=0;pass<POOLFACTOR;pass++){ for(i=0;i<STREAM1;i++) pin[i]^=pp[i]; pp+=STREAM1; }
  pp = s->pool3;
  for(pass=0;pass<POOLFACTOR;pass++){ for(i=0;i<STREAM1;i++) pin[i]^=pp[i]; pp+=STREAM1; }
  for(i=0;i<(s->hashbitlen/32);i++){
    hashval[4*i]   = (unsigned char)pin[i];
    hashval[4*i+1] = (unsigned char)(pin[i]>>8);
    hashval[4*i+2] = (unsigned char)(pin[i]>>16);
    hashval[4*i+3] = (unsigned char)(pin[i]>>24);
  }
}

static void hhash(int hashbitlen, const unsigned char *data, uint64_t nbytes, unsigned char *out){
  hstate s; hinit(&s,hashbitlen); hupdate(&s,data,nbytes); hfinal(&s,out);
}
#endif
