/* independent audit driver: links ONE Waterfall variant at a time */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned char BitSequence;
typedef unsigned long long DataLength;
typedef enum { SUCCESS = 0, FAIL = 1, BAD_HASHBITLEN = 2 } HashReturn;
extern HashReturn Hash(int, const BitSequence*, DataLength, BitSequence*);
extern HashReturn Init(void*, int);
int main(int argc, char**argv){
  int hbl = argc>2 ? atoi(argv[2]) : 256;
  FILE*f=fopen(argv[1],"rb"); if(!f){perror("open");return 2;}
  static unsigned char buf[1<<20]; size_t n=fread(buf,1,sizeof buf,f); fclose(f);
  unsigned char out[64];
  HashReturn r = Hash(hbl, buf, (DataLength)n*8ULL, out);
  if(r!=SUCCESS){fprintf(stderr,"Hash ret %d\n",r);return 3;}
  for(int i=0;i<hbl/8;i++) printf("%02x", out[i]);
  printf("\n");
  return 0;
}
