/* xcheck.c -- cross-implementation consistency check.
 *
 * Links TWO of the authors' own SUPERCOP implementations of COLM_0 into one
 * binary (symbol-renamed at compile time via -D macros in run.sh) and checks that
 * they produce BYTE-IDENTICAL ciphertexts for the same inputs.  This is the
 * known-answer test that says the two builds really are the same function, so
 * an attack that works on one works on the other for a reason and not by luck.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ref_encrypt(unsigned char*,unsigned long long*,const unsigned char*,unsigned long long,const unsigned char*,unsigned long long,const unsigned char*,const unsigned char*,const unsigned char*);
int ref_decrypt(unsigned char*,unsigned long long*,unsigned char*,const unsigned char*,unsigned long long,const unsigned char*,unsigned long long,const unsigned char*,const unsigned char*);
int ni_encrypt(unsigned char*,unsigned long long*,const unsigned char*,unsigned long long,const unsigned char*,unsigned long long,const unsigned char*,const unsigned char*,const unsigned char*);
int ni_decrypt(unsigned char*,unsigned long long*,unsigned char*,const unsigned char*,unsigned long long,const unsigned char*,unsigned long long,const unsigned char*,const unsigned char*);

static void fill(unsigned char*b,long n,unsigned s){long i;for(i=0;i<n;i++){s=s*1103515245u+12345u;b[i]=(unsigned char)(s>>16);} }

int main(int argc,char**argv){
  int trials=(argc>1)?atoi(argv[1]):20;
  unsigned seed=(argc>2)?(unsigned)strtoul(argv[2],0,10):1u;
  int same=0, kat=0, niacc=0, refacc=0;
  /* the splice index set of the fade variant, which both builds must honour */
  int I0[5]={1,2,3,8,129};
  for(int t=0;t<trials;t++){
    unsigned char key[16],np1[8],np2[8],ad[32];
    long L=130, MLEN=16*L, CLEN=MLEN+16;
    unsigned char *m=malloc(MLEN),*c1=malloc(CLEN),*c2=malloc(CLEN);
    unsigned char *d1=malloc(CLEN),*d2=malloc(CLEN),*cm=malloc(CLEN),*o1=malloc(CLEN),*o2=malloc(CLEN);
    unsigned long long x1=0,x2=0,y1=0,y2=0,z1=0,z2=0;
    fill(key,16,seed+t*7919u); fill(np1,8,seed+t*104729u+1); fill(np2,8,seed+t*104729u+2);
    fill(ad,32,seed+t*15485863u+3); fill(m,MLEN,seed+t*32452843u+4);
    ref_encrypt(c1,&x1,m,MLEN,ad,32,NULL,np1,key);
    ni_encrypt (d1,&y1,m,MLEN,ad,32,NULL,np1,key);
    ref_encrypt(c2,&x2,m,MLEN,ad,32,NULL,np2,key);
    ni_encrypt (d2,&y2,m,MLEN,ad,32,NULL,np2,key);
    if(x1==y1 && x2==y2 && !memcmp(c1,d1,x1) && !memcmp(c2,d2,x2)) same++;
    /* honest round trip through the AES-NI build */
    { unsigned long long ml=0; int r=ni_decrypt(o1,&ml,NULL,d1,y1,ad,32,np1,key);
      if(r==0 && ml==(unsigned long long)MLEN && !memcmp(o1,m,MLEN)) kat++; }
    /* the SAME forgery, submitted to BOTH builds */
    memcpy(cm,c1,CLEN);
    for(int i=0;i<5;i++) memcpy(cm+16L*(I0[i]-1),c2+16L*(I0[i]-1),16);
    { unsigned long long ml=0; if(ref_decrypt(o1,&ml,NULL,cm,CLEN,ad,32,np1,key)==0) refacc++; }
    { unsigned long long ml=0; if(ni_decrypt (o2,&ml,NULL,cm,CLEN,ad,32,np1,key)==0) niacc++;  }
    free(m);free(c1);free(c2);free(d1);free(d2);free(cm);free(o1);free(o2);
  }
  printf("  identical ciphertexts from the two builds : %d / %d\n",same,trials);
  printf("  AES-NI build round trips its own output   : %d / %d\n",kat,trials);
  printf("  forgery ACCEPTED by the reference build  : %d / %d\n",refacc,trials);
  printf("  forgery ACCEPTED by the AES-NI build     : %d / %d\n",niacc,trials);
  int ok=(same==trials)&&(kat==trials)&&(refacc==trials)&&(niacc==trials);
  printf(ok?"XCHECK: PASS\n":"XCHECK: FAIL\n");
  return ok?0:1;
}
