/* forge.c -- Zero-patch extraction, signing, and verification driver for
 *            the official SBC_VOLE / SBC_VOLE^M reference implementation
 *            (IACR artifact asiacrypt/2025/a10, Huth & Joux, ASIACRYPT 2025).
 *
 * SBCvole.c and SBCtools.c are compiled 100% UNMODIFIED (zero patches applied;
 * only -Dmain=sbcvole_artifact_main_unused is passed so this file's main()
 * provides the CLI entry point).
 *
 * Commands:
 *   driver keygen <passphrase>
 *       Derives a 256-bit secret key SecX = SHA256(passphrase).
 *   driver sign <w0> <w1> <w2> <w3> <message> <sigfile> <pkfile>
 *       Calls the unmodified CreateKeyFromSecret() and DoVOLESignature(),
 *       writing ONLY the raw binary signature <sigfile> and public key <pkfile>.
 *   driver extract <message> <sigfile> <pkfile> <outfile>
 *       Reads ONLY the public inputs (<message>, <sigfile>, <pkfile>), calls
 *       the unmodified CheckVOLESignature(messhash, UV, CompSig, workMem),
 *       deserialises B from CompSig using the artifact's own F257FromSignature,
 *       reads the verifier's reconstructed evaluations XOpen = F_x(Delta*),
 *       YOpen = F_y(Delta*), and GlobDelta = Delta* from workMem, and writes
 *       the public transcript (DELTA, B, FX, FY, UV) to <outfile>.
 *   driver verify <message> <sigfile> <pkfile>
 *       Calls the unmodified CheckVOLESignature() and prints FORGE_VERIFY: 0|1.
 */
#include <inttypes.h>
#include "SBCunified.h"

extern void CreateKeyFromSecret(u64 SecX[4], u64 UV[2 * 130 * 5], u64 UVxy[4 * 5]);
extern int DoVOLESignature(uchar *messhash, uchar *randomness, u64 *SecX,
                           u64 *UV, u64 *UVxy, uchar *CompSig, int *sigsize,
                           u64 *workMemory);
extern int CheckVOLESignature(uchar *messhash, u64 UV[2 * 130 * 5],
                              uchar *CompSig, u64 *workMemory);
extern void U64FromSignature(u64 *buff, int size, uchar *CompSig, int *sigsize);
extern void F257FromSignature(u64 *buff, int num, uchar *CompSig, int *sigsize,
                              int *leftoverpos, int *leftovernum);

static u64 *alloc_workmem(void)
{
  /* Verbatim from SBCvole.c:1284-1294 */
  int workMemSize = 0;
  workMemSize += 2 * PLAYER_SHARE_MEM * NbRounds * LogPlayers;
  workMemSize += PLAYER_OFFSET_MEM * NbRounds;
  workMemSize += 4 * ExpRounds;
  workMemSize += 4 * LEAVES_SIZE;
  workMemSize += 2 * EXPAND_RATIO * LEAVES_SIZE;
  workMemSize += 4 * EXPAND_RATIO * NbRounds * LogPlayers;
  workMemSize += 4 * DELTA_LEAVES;
  workMemSize += 10 * DELTA_LEAVES;
  workMemSize += 5 * 5 * Lambda;
  return (u64 *)calloc(workMemSize, sizeof(u64));
}

static void hash_message(const char *msg, uchar messhash[32])
{
  SHA256_CTX hash_ctx;
  SHA256Init(&hash_ctx);
  SHA256Update(&hash_ctx, (uchar *)msg, strlen(msg));
  SHA256Final(messhash, &hash_ctx);
}

static void extract_public_transcript(uchar *CompSig, u64 UV[2 * 130 * 5],
                                      u64 *workMemory, FILE *out)
{
  /* 1. Locate Tvec, XOpen, YOpen in workMemory using the exact pointer layout
   *    of CheckVOLESignature (SBCvole.c:1044-1055). */
  u64 *TmpAllocPtr = workMemory;
  u64 *HyperMemorys = AssignNewU64Ptr(TmpAllocPtr, 2 * PLAYER_SHARE_MEM * NbRounds * LogPlayers);
  u64 *Offsets      = AssignNewU64Ptr(TmpAllocPtr, PLAYER_OFFSET_MEM * NbRounds);
  u64 *PreTree      = AssignNewU64Ptr(TmpAllocPtr, 4 * ExpRounds);
  u64 *FullTree     = AssignNewU64Ptr(TmpAllocPtr, 4 * LEAVES_SIZE);
  u64 *ExpandedTree = AssignNewU64Ptr(TmpAllocPtr, 2 * EXPAND_RATIO * LEAVES_SIZE);
  u64 *FoldedExpTree= AssignNewU64Ptr(TmpAllocPtr, 4 * EXPAND_RATIO * NbRounds * LogPlayers);
  u64 *TTree        = AssignNewU64Ptr(TmpAllocPtr, 4 * DELTA_LEAVES);
  u64 *Tvec         = AssignNewU64Ptr(TmpAllocPtr, 10 * DELTA_LEAVES);
  u64 *XOpen        = AssignNewU64Ptr(TmpAllocPtr, 5 * Lambda);
  u64 *YOpen        = AssignNewU64Ptr(TmpAllocPtr, 5 * Lambda);
  (void)HyperMemorys; (void)Offsets; (void)PreTree; (void)FullTree;
  (void)ExpandedTree; (void)FoldedExpTree; (void)TTree;

  /* 2. Recompute GlobDelta = Delta* from Queries[r] and Tvec (SBCvole.c:1058-1065, 736-757). */
  u64 GlobDelta[5] = {0, 0, 0, 0, 0};
  int r, i, lvl, pos, q;
  for (r = 0; r < NbRounds; r++) {
#ifdef L8
    int qr = CompSig[r];
#else
    int qr = ((((int)CompSig[2 * r]) << 8) + CompSig[2 * r + 1]) & ((1 << LogPlayers) - 1);
#endif
    for (i = 0; i < LogPlayers; i++) {
      int good = ((qr >> i) & 1) ^ 1;
      int counter = r * LogPlayers + i;
      if (good) {
        for (pos = 0; pos < 5; pos++) {
          GlobDelta[pos] ^= Tvec[5 * Lambda + 5 * counter + pos];
        }
      }
    }
  }

  /* 3. Deserialise B from CompSig using the artifact's own U64FromSignature and
   *    F257FromSignature (SBCvole.c:1067-1157). */
  int sigoffset = 32;
  u64 dummy[10];
#ifdef HALFKEYSEED
  U64FromSignature(dummy, 2, CompSig, &sigoffset);
#else
  U64FromSignature(dummy, 4, CompSig, &sigoffset);
#endif
  for (r = 0; r < NbRounds; r++) {
    for (lvl = 0; lvl < LogPlayers; lvl++) {
      U64FromSignature(dummy, 2, CompSig, &sigoffset);
    }
  }
  for (r = 0; r < NbRounds; r++) {
    U64FromSignature(dummy, 2, CompSig, &sigoffset);
#ifndef MINMAXCOMMIT
#ifdef SINGLECOMMITSIZE
    U64FromSignature(dummy, 2, CompSig, &sigoffset);
#else
    U64FromSignature(dummy, 4, CompSig, &sigoffset);
#endif
#endif
  }
  int leftoverpos = 0, leftovernum = 8;
  for (r = 1; r < NbRounds; r++) {
    F257FromSignature(dummy, 1, CompSig, &sigoffset, &leftoverpos, &leftovernum);
  }
  u64 B_and_b[10];
  F257FromSignature(B_and_b, 2, CompSig, &sigoffset, &leftoverpos, &leftovernum);

  /* 4. Write the 5 public transcript lines (DELTA, B, FX, FY, UV). */
  fprintf(out, "DELTA"); for (q = 0; q < 5; q++) fprintf(out, " %016" PRIx64, GlobDelta[q]); fprintf(out, "\n");
  fprintf(out, "B");     for (q = 0; q < 5; q++) fprintf(out, " %016" PRIx64, B_and_b[q]);   fprintf(out, "\n");
  fprintf(out, "FX");    for (q = 0; q < 5 * Lambda; q++) fprintf(out, " %016" PRIx64, XOpen[q]); fprintf(out, "\n");
  fprintf(out, "FY");    for (q = 0; q < 5 * Lambda; q++) fprintf(out, " %016" PRIx64, YOpen[q]); fprintf(out, "\n");
  fprintf(out, "UV");    for (q = 0; q < 2 * 130 * 5; q++) fprintf(out, " %016" PRIx64, UV[q]);   fprintf(out, "\n");
}

int main(int argc, char **argv)
{
  SHA256_CTX hash_ctx;
  u64 SecX[4];
  u64 UV[2 * 130 * 5];
  u64 UVxy[4 * 5];
  uchar messhash[32];
  uchar randomness[32];
  uchar CompSig[10000];
  u64 *workMem;
  int sigsize = 0, abort = 1, ok, nbsign = 0;
  int i;
  FILE *f;

  if (argc < 3) {
    fprintf(stderr,
            "usage: %s keygen <passphrase>\n"
            "       %s sign <w0> <w1> <w2> <w3> <message> <sigfile> <pkfile>\n"
            "       %s extract <message> <sigfile> <pkfile> <outfile>\n"
            "       %s verify <message> <sigfile> <pkfile>\n",
            argv[0], argv[0], argv[0], argv[0]);
    return 2;
  }

  if (!strcmp(argv[1], "keygen")) {
    SHA256Init(&hash_ctx);
    SHA256Update(&hash_ctx, (uchar *)argv[2], strlen(argv[2]));
    SHA256Final((uchar *)SecX, &hash_ctx);
    printf("SECX %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n",
           SecX[0], SecX[1], SecX[2], SecX[3]);
    return 0;
  }

  if (!strcmp(argv[1], "sign")) {
    if (argc < 9) { fprintf(stderr, "sign: bad args\n"); return 2; }
    for (i = 0; i < 4; i++) SecX[i] = strtoull(argv[2 + i], NULL, 16);
    hash_message(argv[6], messhash);

    SHA256Init(&hash_ctx);
    SHA256Update(&hash_ctx, (uchar *)"BlaBlaRandom", 12);
    SHA256Final(randomness, &hash_ctx);

    workMem = alloc_workmem();
    CreateKeyFromSecret(SecX, UV, UVxy);
    while (abort == 1) {
      abort = DoVOLESignature(messhash, randomness, SecX, UV, UVxy, CompSig,
                              &sigsize, workMem);
      nbsign++;
      randomness[0]++;
      if (nbsign > 10000) { fprintf(stderr, "too many aborts\n"); return 3; }
    }
    ok = CheckVOLESignature(messhash, UV, CompSig, workMem);

    f = fopen(argv[7], "wb"); fwrite(CompSig, 1, sigsize, f); fclose(f);
    f = fopen(argv[8], "wb"); fwrite(UV, sizeof(u64), 2 * 130 * 5, f); fclose(f);

    printf("SIGN_D=%d SIGN_TAU=%d SIGSIZE=%d\n", LogPlayers, NbRounds, sigsize);
    printf("SIGN_MSGHASH ");
    for (i = 0; i < 32; i++) printf("%02x", messhash[i]);
    printf("\n");
    printf("SIGN_PK0 %016" PRIx64 "\n", UV[0]);
    printf("SIGN_SELFVERIFY: %d\n", ok);
    free(workMem);
    return 0;
  }

  if (!strcmp(argv[1], "extract")) {
    if (argc < 6) { fprintf(stderr, "extract: bad args\n"); return 2; }
    hash_message(argv[2], messhash);
    memset(CompSig, 0, sizeof(CompSig));
    f = fopen(argv[3], "rb");
    if (!f) { fprintf(stderr, "no sigfile\n"); return 2; }
    sigsize = (int)fread(CompSig, 1, sizeof(CompSig), f);
    fclose(f);
    memset(UV, 0, sizeof(UV));
    f = fopen(argv[4], "rb");
    if (!f) { fprintf(stderr, "no pkfile\n"); return 2; }
    if (fread(UV, sizeof(u64), 2 * 130 * 5, f) != 2 * 130 * 5) {
      fprintf(stderr, "short pkfile\n"); return 2;
    }
    fclose(f);
    workMem = alloc_workmem();
    ok = CheckVOLESignature(messhash, UV, CompSig, workMem);
    f = fopen(argv[5], "w");
    if (!f) { fprintf(stderr, "cannot open outfile\n"); free(workMem); return 2; }
    extract_public_transcript(CompSig, UV, workMem, f);
    fclose(f);
    free(workMem);
    if (!ok) {
      fprintf(stderr, "extract: CheckVOLESignature rejected input signature\n");
      return 1;
    }
    return 0;
  }

  if (!strcmp(argv[1], "verify")) {
    if (argc < 5) { fprintf(stderr, "verify: bad args\n"); return 2; }
    hash_message(argv[2], messhash);
    memset(CompSig, 0, sizeof(CompSig));
    f = fopen(argv[3], "rb");
    if (!f) { fprintf(stderr, "no sigfile\n"); return 2; }
    sigsize = (int)fread(CompSig, 1, sizeof(CompSig), f);
    fclose(f);
    memset(UV, 0, sizeof(UV));
    f = fopen(argv[4], "rb");
    if (!f) { fprintf(stderr, "no pkfile\n"); return 2; }
    if (fread(UV, sizeof(u64), 2 * 130 * 5, f) != 2 * 130 * 5) {
      fprintf(stderr, "short pkfile\n"); return 2;
    }
    fclose(f);
    workMem = alloc_workmem();
    ok = CheckVOLESignature(messhash, UV, CompSig, workMem);
    printf("VERIFY_D=%d VERIFY_TAU=%d SIGSIZE=%d\n", LogPlayers, NbRounds, sigsize);
    printf("VERIFY_MSGHASH ");
    for (i = 0; i < 32; i++) printf("%02x", messhash[i]);
    printf("\n");
    printf("VERIFY_PK0 %016" PRIx64 "\n", UV[0]);
    printf("FORGE_VERIFY: %d\n", ok);
    free(workMem);
    return ok ? 0 : 1;
  }

  fprintf(stderr, "unknown mode %s\n", argv[1]);
  return 2;
}
