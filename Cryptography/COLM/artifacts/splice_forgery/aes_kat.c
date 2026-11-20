#include <stdio.h>
#include <stdlib.h>
#include <openssl/aes.h>
int main(int argc, char**argv){
  /* FIPS-197 C.1 */
  unsigned char key[16]; for(int i=0;i<16;i++) key[i]=i;
  unsigned char pt[16]={0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88,0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xff};
  unsigned char ct[16], back[16];
  AES_KEY k1,k2; AES_set_encrypt_key(key,128,&k1); AES_set_decrypt_key(key,128,&k2);
  AES_encrypt(pt,ct,&k1); AES_decrypt(ct,back,&k2);
  printf("FIPS197 ct = "); for(int i=0;i<16;i++) printf("%02x",ct[i]); printf("\n");
  printf("expect     = 69c4e0d86a7b0430d8cdb78070b4c55a\n");
  int ok = 1; for(int i=0;i<16;i++) if(back[i]!=pt[i]) ok=0;
  printf("roundtrip  = %s\n", ok?"OK":"BAD");
  /* random vectors for cross-check against pycryptodome */
  if(argc>1){
    FILE*f=fopen(argv[1],"r"); char line[256];
    while(fgets(line,sizeof line,f)){
      char k[40],p[40],c[40]; 
      if(sscanf(line,"%32s %32s %32s", k,p,c)!=3) continue;
      unsigned char kb[16],pb[16],cb[16];
      for(int i=0;i<16;i++){ sscanf((char*)k+2*i,"%2hhx",&kb[i]); sscanf((char*)p+2*i,"%2hhx",&pb[i]); sscanf((char*)c+2*i,"%2hhx",&cb[i]); }
      unsigned char out[16]; AES_KEY kk; AES_set_encrypt_key(kb,128,&kk); AES_encrypt(pb,out,&kk);
      for(int i=0;i<16;i++) if(out[i]!=cb[i]){ printf("MISMATCH at %s\n",line); return 1; }
      AES_set_decrypt_key(kb,128,&kk); AES_decrypt(cb,out,&kk);
      for(int i=0;i<16;i++) if(out[i]!=pb[i]){ printf("DEC MISMATCH at %s\n",line); return 1; }
    }
    printf("cross-check vs pycryptodome: all vectors OK\n");
  }
  return 0;
}
