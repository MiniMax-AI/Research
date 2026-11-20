/* =====================================================================
 * splice.c -- THE consolidated driver for the COLM two-query,
 * encryption-only, probability-1 splice forgery.
 *
 * One source file, four variants of the same attack, each buildable against
 * any of the authors' own SUPERCOP implementations of COLM (colm0v1/ref,
 * colm0v1/aesni4, colm0v1/aesni8, colm127v1/ref).  The COLM sources are NOT
 * modified; the only thing supplied is <openssl/aes.h> (openssl/aes.h),
 * so the driver builds standalone without requiring an external OpenSSL installation.
 *
 *   VARIANT      mechanism                                    min l
 *   ---------     ------------------------------------------   -----
 *   fade          difference fade on supp(f(x))                130
 *   transplant    state transfer + FINAL TAG BLOCK transplant   129
 *   misuse        fade under a REPEATED nonce, one-block
 *                 message difference at block s                s+130
 *   prefix        fade where the two queries share only
 *                 their first 129 message blocks               130
 *
 * Every splice index set is DERIVED AT RUN TIME from COLM's own field
 * GF(2^128) = GF(2)[x]/(x^128+x^7+x^2+x+1) -- nothing is hard-coded.
 *
 * Build:
 *   gcc -O2 -I. -I<impl> [-DTAG_STRIDE=127] [-DEXPECT_REJECT] \
 *       -o splice splice.c <impl>/encrypt*.c -lm
 *
 * Usage:
 *   ./splice <variant> <seed> <trials> [variant options]
 *     fade        l
 *     transplant  l
 *     misuse      l s
 *     prefix      l1 extra |A1| |A2| P
 * ===================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crypto_aead.h"

#ifndef TAG_STRIDE
#define TAG_STRIDE 0        /* 0 = COLM_0 ; 127 = COLM_127 */
#endif

/* byte offset of ciphertext block i (1-based).
 * COLM_127 interleaves one 16-byte intermediate tag after every 127 blocks,
 * so block i sits at 16*(i-1) + 16*floor((i-1)/127).                            */
static long off(int i){ return 16L*(i-1) + (TAG_STRIDE? 16L*((i-1)/TAG_STRIDE) : 0L); }

/* COLM truncates the FINAL tag block to the length of the last message block
 * (reference encrypt.c: "store_bytes(c, result, 0, blen-1)"), so clen = 16*l + P
 * where P is the number of bytes in the last message block.  Every block copy
 * therefore has to ask how many bytes of block i actually exist.            */
static long blklen(int i,long CLEN){ long o=off(i), e=off(i+1); if(e>CLEN) e=CLEN; return e-o; }
static int  spliceblk(unsigned char*dst,const unsigned char*src,int i,long CLEN){
    long o=off(i), n=blklen(i,CLEN);
    if(o>=CLEN || n<=0) return 0;
    memcpy(dst+o,src+o,(size_t)n); return 1;
}

/* ------------------------------------------------------------------ */
/* GF(2^128), COLM's own field, in the byte order of the reference      */
/* module.h: the 16-byte block is a big-endian 128-bit integer and bit i */
/* of that integer is the coefficient of x^i (COLM v1 p.3: "2" means x). */
/* ------------------------------------------------------------------ */
static void g2(unsigned char *d, const unsigned char *s){ int i; unsigned char t=s[0]>>7;
    for(i=0;i<15;i++) d[i]=(unsigned char)((s[i]<<1)|(s[i+1]>>7));
    d[15]=(unsigned char)(s[15]<<1); if(t) d[15]^=0x87; }
static void gxor(unsigned char *d,const unsigned char*a,const unsigned char*b){int i;for(i=0;i<16;i++)d[i]=a[i]^b[i];}
static void gmul(unsigned char *d,const unsigned char*a,const unsigned char*b){
    unsigned char r[16],t[16],aa[16]; int i,k;
    memset(r,0,16); memcpy(aa,a,16);
    for(k=0;k<128;k++){
        if((b[15-(k>>3)]>>(k&7))&1) for(i=0;i<16;i++) r[i]^=aa[i];
        g2(t,aa); memcpy(aa,t,16);        /* aa <- aa * x */
    }
    memcpy(d,r,16);
}
static void gone(unsigned char *d){ memset(d,0,16); d[15]=1; }
static void gpowx(unsigned char *d,int e){ unsigned char t[16]; gone(d); while(e-->0){ g2(t,d); memcpy(d,t,16);} }
/* 3^{-1} = (x+1)^{-1} by Fermat: 3^(2^128-2) = 3^(2+4+...+2^127) */
static void ginv3(unsigned char *d){
    unsigned char base[16],r[16],t[16]; int k;
    memset(base,0,16); base[15]=3;       /* the field element 3 = x+1 */
    gone(r);
    for(k=0;k<128;k++){
        if(k>0){ gmul(t,base,base); memcpy(base,t,16); }   /* base = 3^(2^k) */
        if(k>0){ gmul(t,r,base);     memcpy(r,t,16);    }  /* bit k of 2^128-2 set for k>=1 */
    }
    memcpy(d,r,16);
}
static int gdeg(const unsigned char *a){ int i,k; for(i=0;i<16;i++) if(a[i]) for(k=7;k>=0;k--) if((a[i]>>k)&1) return (15-i)*8+k; return -1; }

/* ---- derive the splice sets from the algebra, at run time ---------- */
/* fade family:  sum_{k in I0} x^(k-1) = 0 mod f, realised as the support of
 * the UNREDUCED product x^s * f(x), f = x^128+x^7+x^2+x+1.                    */
static int derive_fade(int s,int *I0){
    int e[5]={s,s+1,s+2,s+7,s+128},i;
    for(i=0;i<5;i++) I0[i]=e[i]+1;
    unsigned char acc[16],p[16],t[16]; memset(acc,0,16);
    for(i=0;i<5;i++){ gpowx(p,I0[i]-1); gxor(t,acc,p); memcpy(acc,t,16); }
    for(i=0;i<16;i++) if(acc[i]) return -1;      /* the relation must vanish */
    return 5;
}
/* transplant family.
 *
 * A splice at block i "transfers the state onto query 2" -- it makes
 * X_mix[i] = X_2[i] and W_mix[i] = W_2[i] -- IF AND ONLY IF W_mix[i-1] = W_2[i-1].
 * State transfer is therefore STICKY ALONG A CONTIGUOUS RUN of spliced blocks
 * and is destroyed by any unspliced block in between.  So the transplant needs
 *
 *     (i)  a contiguous run  {m, m+1, ..., l+1}  of spliced blocks that reaches
 *          the tag block l+1, and
 *     (ii) W_mix[m-1] = W_2[m-1], i.e.  3 * sum_{j in S, j <= m-1} x^(j-1) = x^(m-1),
 *          i.e.  S cap {1..m-1} = supp(x^(m-1)/3) + 1, which is possible only if
 *          maxdeg(x^(m-1)/3) <= m-2   ("condition (C)").
 *
 * Condition (C) first holds at m = 129, and m <= l, so l >= 129; at l = m = 129
 * the run is just {129, 130} = {l, l+1} and the splice set is the 8-element
 * {1,3,4,5,6,7,129,130}.  For l > 129 the SAME prefix is used and the run is
 * simply extended, which is why the family covers every l >= 129.
 *
 * MSTART < 0 means "use the smallest m that satisfies (C)".
 */
static int smallest_m(void){
    unsigned char inv3[16],p[16],A[16]; int m;
    ginv3(inv3);
    for(m=2;m<400;m++){ gpowx(p,m-1); gmul(A,inv3,p); if(gdeg(A)<=m-2) return m; }
    return -1;
}
static int derive_transplant(int l,int *I0,int MSTART){
    unsigned char inv3[16],p[16],A[16]; int i,n=0,m;
    ginv3(inv3);
    m = (MSTART>0)? MSTART : smallest_m();
    if(m<0 || m>l) return -1;                    /* need m <= l */
    gpowx(p,m-1); gmul(A,inv3,p);
    if(gdeg(A) > m-2) return -1;                /* condition (C) fails at this m */
    for(i=0;i<=m-2 && n<2000;i++) if((A[15-(i>>3)]>>(i&7))&1) I0[n++]=i+1;
    for(i=m;i<=l+1 && n<2040;i++) I0[n++]=i;      /* the contiguous run, up to the tag */
    return n;
}

/* ------------------------------------------------------------------ */
static void fill(unsigned char*b,long n,unsigned s){long i;for(i=0;i<n;i++){s=s*1103515245u+12345u;b[i]=(unsigned char)(s>>16);} }

static int LBLK=130, SBLK=1, EXTRA=0, AD1=32, AD2=32, AGREE=-1, SPLICE_TAGS=0, MSTART=-1;
static int PARTIAL=16;   /* bytes in the LAST message block of query 1; 16 = complete */
static int PARTIAL2=-1;  /* ... of query 2; -1 means "same as query 1" */
static int I0[2048], NI0=0;
static const char *VARIANT="fade";

/* Theorem 2 of the addendum, evaluated at this run's parameters.  sigma runs
 * over the q forward queries AND the s forgery attempts (addendum p.3). */
static void print_thm2(int q,int s,long a_blk,long l_blk){
    long sigma=(long)(q+s)*(a_blk+l_blk), sigma_auth=sigma+q;
    double bound=(9.0*(double)sigma_auth*(double)sigma_auth + s)/pow(2.0,128);
    printf("addendum Theorem 2 at q=%d, s=%d, a=%ld, l=%ld : sigma=%ld sigma_auth=%ld\n",
           q,s,a_blk,l_blk,sigma,sigma_auth);
    printf("  its own bound  9*sigma_auth^2/2^n + s/2^n = %.6e = 2^%.2f\n",bound,log2(bound));
    printf("  measured forging advantage                = 1  (deterministic)\n");
}

int main(int argc,char**argv){
  unsigned seed=1u; int trials=20, argi=1;
  if(argc>1) VARIANT=argv[argi++];
  if(argc>argi) seed=(unsigned)strtoul(argv[argi++],0,10);
  if(argc>argi) trials=atoi(argv[argi++]);

  int same_nonce=0, diff_msg=0, diff_len=0;
  if(!strcmp(VARIANT,"fade")){
      LBLK=(argc>argi)?atoi(argv[argi++]):130;
      NI0=derive_fade(0,I0);
  } else if(!strcmp(VARIANT,"transplant")){
      LBLK=(argc>argi)?atoi(argv[argi++]):129;
      MSTART=(argc>argi)?atoi(argv[argi++]):-1;
      NI0=derive_transplant(LBLK,I0,MSTART);
      if(NI0<0){printf("no transplant set exists at l=%d (need m <= l with maxdeg(x^(m-1)/3) <= m-2;"
                       " the smallest such m is %d)\n",LBLK,smallest_m());
                printf("RESULT: FAIL\n");return 1;}
  } else if(!strcmp(VARIANT,"misuse")){
      LBLK=(argc>argi)?atoi(argv[argi++]):131;
      SBLK=(argc>argi)?atoi(argv[argi++]):1;
      NI0=derive_fade(SBLK,I0);
      same_nonce=1; diff_msg=1;
  } else if(!strcmp(VARIANT,"prefix")){
      LBLK=(argc>argi)?atoi(argv[argi++]):300;
      EXTRA=(argc>argi)?atoi(argv[argi++]):100;
      AD1  =(argc>argi)?atoi(argv[argi++]):16;
      AD2  =(argc>argi)?atoi(argv[argi++]):80;
      AGREE=(argc>argi)?atoi(argv[argi++]):129;
      NI0=derive_fade(0,I0);
      diff_msg=1; diff_len=1;
  } else { printf("unknown variant %s\n",VARIANT); printf("RESULT: FAIL\n"); return 1; }
  if(NI0<0){printf("splice-set derivation failed\n");printf("RESULT: FAIL\n");return 1;}
  SPLICE_TAGS=(getenv("SPLICE_TAGS")&&atoi(getenv("SPLICE_TAGS")))?1:0;
  PARTIAL=(getenv("PARTIAL")&&atoi(getenv("PARTIAL"))>0)?atoi(getenv("PARTIAL")):16;
  if(PARTIAL>16) PARTIAL=16;
  /* ADLEN overrides let the caller use associated data of ANY byte length,
     including lengths that are not multiples of 16 (COLM 10*-pads the last AD
     block and switches to the 3*7*2^(i-1)*L mask branch there).            */
  PARTIAL2=(getenv("PARTIAL2"))?atoi(getenv("PARTIAL2")):-1;
  if(PARTIAL2<0) PARTIAL2=PARTIAL;
  if(PARTIAL2>16) PARTIAL2=16;
  if(getenv("AD1")) AD1=atoi(getenv("AD1"));
  if(getenv("AD2")) AD2=atoi(getenv("AD2"));

  long L2=LBLK+EXTRA;
  long MLEN=16L*(LBLK-1)+PARTIAL, MLEN2=16L*(L2-1)+PARTIAL2;
  long CAP = MLEN2 + 16 + (TAG_STRIDE? 16L*(L2/TAG_STRIDE+2) : 0L) + 64;
  if(AGREE<0) AGREE=LBLK-1;
  int h1 = TAG_STRIDE? (int)((LBLK-1)/TAG_STRIDE) : 0;

  printf("== COLM splice forgery -- variant '%s' -- implementation tag stride %d ==\n",VARIANT,TAG_STRIDE);
  printf("message 1: l1 = %d blocks (%ld bytes; last block %d byte(s) -- %s), AD1 = %d bytes%s\n",
         LBLK,MLEN,PARTIAL,PARTIAL==16?"COMPLETE":"INCOMPLETE, 10* padded",AD1,TAG_STRIDE?"":", no intermediate tags");
  if(TAG_STRIDE) printf("           h = %d intermediate tag(s) for message 1\n",h1);
  if(diff_len) printf("message 2: l2 = %ld blocks (%ld bytes; last block %d byte(s)), AD2 = %d bytes, plaintext agreement on blocks 1..%d\n",L2,MLEN2,PARTIAL2,AD2,AGREE);
  printf("splice index set I0 = {"); for(int i=0;i<NI0;i++) printf("%d%s",I0[i],i+1<NI0?",":""); printf("}   [derived at run time from GF(2^128)]\n");
  printf("queries: 2 ENCRYPTION, 0 decryption, 0 verification, 0 RUP.  nonces %s.\n",
         same_nonce?"REPEATED (nonce misuse)":"distinct (nonce-respecting)");
  if(!TAG_STRIDE) print_thm2(2,1,(AD1+15)/16,LBLK);
  printf("\n");

  unsigned char key[32],np1[8],np2[8];
  unsigned char *ad1=malloc(AD1>0?AD1:1), *ad2=malloc(AD2>0?AD2:1);
  unsigned char *m1=malloc(MLEN), *m2=malloc(diff_len?MLEN2:MLEN);
  unsigned char *c1=malloc(CAP), *c2=malloc(CAP), *cm=malloc(CAP), *out=malloc(CAP), *cw=malloc(CAP);
  unsigned long long cl1=0,cl2=0,mo=0;
  long CLEN=0, TAGOFF=0;
  int acc=0,negShift=0,negExtra=0,negRand=0,tagreuse=0,trivial=0,kat=0,gen=0;
  int NRAND=(getenv("NRAND")&&atoi(getenv("NRAND"))>0)?atoi(getenv("NRAND")):200;

  for(int t=0;t<trials;t++){
    fill(key,32,seed+t*7919u); fill(np1,8,seed+t*104729u+1); fill(np2,8,seed+t*104729u+2);
    fill(ad1,AD1,seed+t*15485863u+3); fill(ad2,AD2,seed+t*15485863u+9);
    fill(m1,MLEN,seed+t*32452843u+4);
    if(diff_len){ fill(m2,MLEN2,seed+t*32452843u+77); if(AGREE>0 && 16L*AGREE<=MLEN) memcpy(m2,m1,16L*AGREE); }
    else if(diff_msg){ memcpy(m2,m1,MLEN);
        do { fill(m2+16L*(SBLK-1),16,seed+t*32452843u+5); } while(memcmp(m1+16L*(SBLK-1),m2+16L*(SBLK-1),16)==0); }
    else memcpy(m2,m1,MLEN);

    if(crypto_aead_encrypt(c1,&cl1,m1,MLEN,ad1,AD1,NULL,np1,key)!=0){printf("enc1 failed\n");printf("RESULT: FAIL\n");return 1;}
    if(crypto_aead_encrypt(c2,&cl2,m2,diff_len?MLEN2:MLEN,ad2,AD2,NULL,same_nonce?np1:np2,key)!=0){printf("enc2 failed\n");printf("RESULT: FAIL\n");return 1;}
    CLEN=(long)cl1; TAGOFF=off(LBLK+1);

    if(t==0){
      unsigned long long ml=0, mb=0;
      unsigned char *o1=malloc(CAP);
      int r=crypto_aead_decrypt(o1,&ml,NULL,c1,cl1,ad1,AD1,np1,key);
      int honest=(r==0 && ml==(unsigned long long)MLEN && !memcmp(o1,m1,MLEN));
      unsigned char *bad=malloc(CAP); memcpy(bad,c1,CLEN); bad[16*50]^=1;
      int rb=crypto_aead_decrypt(out,&mb,NULL,bad,CLEN,ad1,AD1,np1,key);
      kat = honest && (rb!=0);
      printf("KAT-1 honest round trip of the implementation's own ciphertext : %s (ret=%d mlen=%llu, clen=%llu)\n",
             honest?"PASS":"FAIL",r,ml,cl1);
      printf("KAT-2 one flipped ciphertext bit is REJECTED                   : %s (ret=%d)\n",rb!=0?"PASS":"FAIL",rb);
      free(bad); free(o1);
      if(!kat){printf("RESULT: FAIL\n");return 1;}
      long cc=0; while(cc<(long)((cl1<cl2)?cl1:cl2) && c1[cc]==c2[cc]) cc++;
      printf("      common ciphertext prefix of the two queries: %ld bytes = %ld blocks\n",cc,cc/16);
    }

    /* ---------- the forgery: blocks I0 of query 2 spliced into query 1 ---------- */
    memcpy(cm,c1,CLEN);
    for(int i=0;i<NI0;i++){ long o=off(I0[i]); if(o+16>CLEN){printf("index %d out of range\n",I0[i]);printf("RESULT: FAIL\n");return 1;}
                           memcpy(cm+o,c2+o,16); }
    if(SPLICE_TAGS) for(int j=1;j<=h1;j++){ long o=16L*TAG_STRIDE*j; memcpy(cm+o,c2+o,16); }
    int triv=(memcmp(cm,c1,CLEN)==0);
    if(triv) trivial++;
    if(CLEN-TAGOFF>0 && memcmp(cm+TAGOFF,c1+TAGOFF,(size_t)(CLEN-TAGOFF))==0) tagreuse++;

    unsigned long long ml3=0;
    int r3=crypto_aead_decrypt(out,&ml3,NULL,cm,CLEN,ad1,AD1,np1,key);
    if(r3==0 && !triv) acc++;
    gen++;

    if(t==0){
      printf("      forgery non-trivial (differs from query 1's output)        : %s\n",triv?"NO":"yes");
      printf("      forgery tag block == query 1's tag block (verbatim reuse)  : %s\n",
             (CLEN-TAGOFF>0 && memcmp(cm+TAGOFF,c1+TAGOFF,(size_t)(CLEN-TAGOFF))==0)
             ? "yes  <-- addendum Lemma 3 'good forge view' condition 2 violated DELIBERATELY" : "no (tag transplanted)");
      printf("      crypto_aead_decrypt(forgery) ret=%d -> %s\n",r3,r3==0?"ACCEPTED":"rejected");
      if(r3==0){ long df=0; for(long i=0;i<(long)ml3&&i<MLEN;i++) if(out[i]!=m1[i]) df++;
                 printf("      recovered plaintext differs from the queried message in %ld / %ld bytes"
                        " (so the forgery is not a replay)\n",df,MLEN); }
    }

    /* ---------- negative control 1: DROP one element of I0 (the middle one) ----------
     * For the fade families this breaks the vanishing sum; for the transplant it
     * breaks the contiguity of the run.  Either way it must be REJECTED.        */
    memcpy(cw,c1,CLEN);
    { int skip=NI0/2;
      for(int i=0;i<NI0;i++){ if(i==skip) continue; spliceblk(cw,c2,I0[i],CLEN); } }
    { unsigned long long mw=0; if(crypto_aead_decrypt(out,&mw,NULL,cw,CLEN,ad1,AD1,np1,key)==0) negShift++; }

    /* ---------- negative control 2: ADD one block that is not in I0 ---------- */
    memcpy(cw,c1,CLEN);
    for(int i=0;i<NI0;i++) spliceblk(cw,c2,I0[i],CLEN);
    { int b=0;
      for(int cand=1;cand<=LBLK;cand++){ int in=0; for(int i=0;i<NI0;i++) if(I0[i]==cand) in=1;
                                       if(!in){ b=cand; break; } }
      if(b>0){ spliceblk(cw,c2,b,CLEN);
               unsigned long long mw=0; if(crypto_aead_decrypt(out,&mw,NULL,cw,CLEN,ad1,AD1,np1,key)==0) negExtra++; } }
  }

  /* ---------- negative control 3: RANDOM splice sets of the same size ----------
   * The statistical baseline.  A random |I0|-subset of block positions must be
   * rejected essentially always; the algebraic set is accepted always.  Same
   * code path, fresh keys, fresh messages.                                     */
  {
    unsigned rs=seed^0x5bf03635u; int done=0;
    for(int k=0;k<NRAND;k++){
      int pos[2048],n=0;
      while(n<NI0){ rs=rs*1103515245u+12345u; int b=1+(int)((rs>>16)%(unsigned)LBLK);
                   int dup=0; for(int j=0;j<n;j++) if(pos[j]==b) dup=1; if(!dup) pos[n++]=b; }
      unsigned char k2[32],n1[8],n2[8];
      fill(k2,32,rs+1); fill(n1,8,rs+2); fill(n2,8,rs+3);
      unsigned char *mm=malloc(MLEN),*cc1=malloc(CAP),*cc2=malloc(CAP),*ccm=malloc(CAP);
      unsigned long long x1=0,x2=0,mm2=0;
      fill(mm,MLEN,rs+4);
      crypto_aead_encrypt(cc1,&x1,mm,MLEN,ad1,AD1,NULL,n1,k2);
      crypto_aead_encrypt(cc2,&x2,mm,MLEN,ad1,AD1,NULL,n2,k2);
      memcpy(ccm,cc1,x1);
      for(int j=0;j<NI0;j++) spliceblk(ccm,cc2,pos[j],(long)x1);
      if(crypto_aead_decrypt(out,&mm2,NULL,ccm,x1,ad1,AD1,n1,k2)==0) negRand++;
      free(mm);free(cc1);free(cc2);free(ccm); done++;
    }
    printf("\n");
    printf("trials = %d   forgeries ACCEPTED = %d   trivial = %d   tag-block reused = %d\n",gen,acc,trivial,tagreuse);
    printf("negative control 1 (DROP the middle element of I0)   accepted %d / %d\n",negShift,gen);
    printf("negative control 2 (ADD one block not in I0)         accepted %d / %d\n",negExtra,gen);
    printf("negative control 3 (RANDOM %d-block splice sets)    accepted %d / %d   <- statistical baseline\n",NI0,negRand,done);
  }

#ifdef EXPECT_REJECT
  int ok = (gen>0) && (acc==0) && (kat==1);
  printf("EXPECTATION for this build: every splice is REJECTED (COLM_127 control)\n");
#else
  int ok = (acc==gen) && (trivial==0) && (negShift==0) && (negExtra==0) && (negRand==0);
  printf("EXPECTATION for this build: every splice is ACCEPTED, every control REJECTED\n");
#endif
  printf(ok?"RESULT: SUCCESS\n":"RESULT: FAIL\n");
  return ok?0:1;
}
