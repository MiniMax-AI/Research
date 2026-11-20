/* One-query forgery baseline strategies against COLM_0.
 * Evaluates eight concrete single-query ciphertext manipulation strategies
 * against the unmodified SUPERCOP colm0v1/ref implementation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "crypto_aead.h"
#define ADLEN 32
static void fill(unsigned char*b,long n,unsigned s){long i;for(i=0;i<n;i++){s=s*1103515245u+12345u;b[i]=(unsigned char)(s>>16);} }
int main(int argc,char**argv){
  int trials=(argc>1)?atoi(argv[1]):200, L=(argc>2)?atoi(argv[2]):130;
  long MLEN=16L*L, CLEN=MLEN+16;
  const char *names[8]={"truncate to k blocks, keep the tag","swap two interior blocks",
    "duplicate one interior block","submit under a different nonce","submit under different AD",
    "rotate the ciphertext by one block","drop the tag block","all-zero ciphertext"};
  int acc[8]={0,0,0,0,0,0,0,0};
  for(int t=0;t<trials;t++){
    unsigned char key[16],np[8],np2[8],ad[ADLEN],ad2[ADLEN];
    unsigned char *m=malloc(MLEN),*c=malloc(CLEN),*f=malloc(CLEN),*o=malloc(CLEN);
    unsigned long long cl=0,mo=0;
    fill(key,16,t*7919u+1); fill(np,8,t*104729u+2); fill(np2,8,t*104729u+3);
    fill(ad,ADLEN,t*15485863u+4); fill(ad2,ADLEN,t*15485863u+5); fill(m,MLEN,t*32452843u+6);
    if(crypto_aead_encrypt(c,&cl,m,MLEN,ad,ADLEN,NULL,np,key)!=0){printf("enc fail\n");return 1;}
    if(t==0){ unsigned long long ml=0; int r=crypto_aead_decrypt(o,&ml,NULL,c,cl,ad,ADLEN,np,key);
      printf("KAT honest round trip: %s\n\n",(r==0&&ml==(unsigned long long)MLEN&&!memcmp(o,m,MLEN))?"PASS":"FAIL");
      if(r!=0) return 1; }
    /* 1. truncate to k blocks but keep the original tag block */
    { int k=L/2; memcpy(f,c,CLEN); memmove(f+16L*k, c+16L*L, 16);
      if(crypto_aead_decrypt(o,&mo,NULL,f,16L*k+16,ad,ADLEN,np,key)==0) acc[0]++; }
    /* 2. swap two interior blocks */
    { memcpy(f,c,CLEN); unsigned char tmp[16];
      memcpy(tmp,f+16L*10,16); memcpy(f+16L*10,f+16L*60,16); memcpy(f+16L*60,tmp,16);
      if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad,ADLEN,np,key)==0) acc[1]++; }
    /* 3. duplicate one interior block */
    { memcpy(f,c,CLEN); memcpy(f+16L*10,f+16L*60,16);
      if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad,ADLEN,np,key)==0) acc[2]++; }
    /* 4. same ciphertext, different nonce */
    { memcpy(f,c,CLEN); if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad,ADLEN,np2,key)==0) acc[3]++; }
    /* 5. same ciphertext, different associated data */
    { memcpy(f,c,CLEN); if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad2,ADLEN,np,key)==0) acc[4]++; }
    /* 6. rotate the ciphertext left by one block */
    { memcpy(f,c+16,CLEN-16); memcpy(f+CLEN-16,c,16);
      if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad,ADLEN,np,key)==0) acc[5]++; }
    /* 7. drop the tag block */
    { memcpy(f,c,CLEN); if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN-16,ad,ADLEN,np,key)==0) acc[6]++; }
    /* 8. the all-zero ciphertext */
    { memset(f,0,CLEN); if(crypto_aead_decrypt(o,&mo,NULL,f,CLEN,ad,ADLEN,np,key)==0) acc[7]++; }
    free(m);free(c);free(f);free(o);
  }
  printf("one-query strategies against COLM_0, l = %d blocks, %d trials\n",L,trials);
  int any=0;
  for(int i=0;i<8;i++){ printf("  %-42s accepted %d / %d\n",names[i],acc[i],trials); any+=acc[i]; }
  printf("\ntotal acceptances across all eight strategies: %d\n",any);
  return any;
}
