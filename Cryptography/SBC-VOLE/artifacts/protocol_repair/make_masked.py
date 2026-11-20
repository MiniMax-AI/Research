#!/usr/bin/env python3
"""
Build the MASKED variant of the authors' ASIACRYPT-2025 artifact.

Source of truth: asiacrypt-2025-a10.zip (sha256 e275d529c750bad077ca7bbe0ff355e9a183a4038df08d7869ef7ca95d3447a7),
fetched on demand by ../fetch_third_party.sh and unpacked read-only into a
temporary staging directory. This script writes a *new* source tree and never
edits the original archive.

What -DMASKB adds, and why each line is there
---------------------------------------------
The paper (Fig. 9 step 17, p.23) transmits the Fiat-Shamir check coefficient B
in the clear, and B is an F_2-affine functional of the 256-bit secret key, so a
single signature yields the key by one 257x256 GF(2) elimination
(Section 4, 2^25.2 bit ops).

The repair is the masking step QuickSilver (ePrint 2021/076 Sec 3.1) and
VOLEitH (ePrint 2023/996 Sec 6.2) both perform and this construction omits: add
one more global cGGM row whose value W is forced equal across all tau trees by
per-tree offsets, transmit B~ = B + W instead of B, and let the verifier cancel
the mask using its own punctured view of the new row.

Algebra (all in F_{2^257}, char 2, so + == XOR)
  signer   : Mtot = XOR_r M0_r ,  B~ = B + W ,  A~ = A + Mtot        (A~ hashed into h_1)
  verifier : Qtot = XOR_r Q_r  ,  A'' = Phi_open + GlobDelta*B~ + Qtot
  identity : Q_r = M0_r + W*RoundDelta_r   =>   Qtot = Mtot + W*GlobDelta
  check    : A'' = (A + GlobDelta*B) + GlobDelta*(B+W) + Qtot
                 = A + GlobDelta*W + Mtot + W*GlobDelta = A + Mtot = A~   OK

-DMASKB_NOH0 omits the mask offsets dw_j from h_0 (Section 6.3 of the paper):
dw_j enters the verifier's check linearly with coefficient Delta*_j and nothing
else in the transcript depends on it, allowing a cheating prover to solve for
it after seeing the challenge.  Used by the soundness-control driver.
"""
import hashlib
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
FETCH_SCRIPT = os.path.normpath(os.path.join(HERE, "..", "fetch_third_party.sh"))
ZIP = os.path.join(HERE, "vendor", "asiacrypt-2025-a10.zip")
ZIP_SHA = "e275d529c750bad077ca7bbe0ff355e9a183a4038df08d7869ef7ca95d3447a7"

EXPECTED_PAT_C_SHA = "60b754925a062076c528c6dadaed191e8a55aa0a7bd6379ee73eb5684e7133de"
EXPECTED_PAT_H_SHA = "daefe26cf96acaea602a425b165539ae067a8c85053387c3b7bc1947daa0794b"


def apply_line_ops(src: str, ops, expected_sha: str, label: str):
    lines = src.splitlines(keepends=True)
    log = []
    # Apply from bottom to top so 0-based line indices into the pristine file remain invariant.
    for i1, i2, payload, tag in sorted(ops, key=lambda x: (x[0], x[1]), reverse=True):
        lines[i1:i2] = payload.splitlines(keepends=True)
        log.append((tag, 1))
    log.reverse()
    out = "".join(lines)
    got = hashlib.sha256(out.encode("utf-8")).hexdigest()
    assert got == expected_sha, f"{label} output sha256 mismatch: {got} != {expected_sha}"
    return out, log


def edits_sbcvole(src: str):
    h0_block = (
        "#if defined(MASKB) && !defined(MASKB_NOH0)\n"
        "  // Bind per-tree mask offsets dw_j into h_0; omitting them leaves dw_j linear\n"
        "  // with coefficient Delta*_j and allows post-challenge forgery (Section 6.3).\n"
        "  for(r=1;r<NbRounds;r++) {\n"
        "    SHA256Update(&hash_ctx, (uchar *)(Offsets+r*PLAYER_OFFSET_MEM+mOFFSET), 5*sizeof(u64));\n"
        "  }\n"
        "#endif\n"
    )
    ops = [
        (592, 592,
         "#ifdef MASKB\n"
         "u64 MaskAtilde[5];   /* A~ = A + Mtot, stashed for the FORGE_NOH0 prover */\n"
         "#endif\n\n",
         "stash-atilde"),
        (627, 627,
         "#ifdef MASKB\n"
         "  u64 M0Round[5], Mtot[5];\n"
         "  for(i=0;i<5;i++) Mtot[i]=0;\n"
         "#endif\n",
         "signer-decl"),
        (656, 656,
         "#ifdef MASKB\n"
         "      for(i=0;i<5;i++) M0Round[i]=0;\n"
         "#endif\n",
         "signer-roundinit"),
        (672, 672,
         "#ifdef MASKB\n"
         "      FFMulAdd(M0Round,HyperMemorys+PLAYER_SHARE_MEM*(2*counter)+mOFFSET,TDeltavec+5*Lambda+5*counter);\n"
         "#endif\n",
         "signer-ffmuladd"),
        (684, 684,
         "#ifdef MASKB\n"
         "    for(i=0;i<5;i++) Mtot[i]^=M0Round[i];   /* signer: accumulate Mtot */\n"
         "#endif\n",
         "signer-mtot"),
        (706, 707,
         "#ifdef MASKB\n"
         "  // W = Offsets[mOFFSET..] is round 0's raw mask value; the force-equal loop\n"
         "  // below has left it untouched and made Offsets[r][m] = raw_r + W for r>=1.\n"
         "  {\n"
         "    u64 Atilde[5], Braw[5];\n"
         "    for(i=0;i<5;i++) Braw[i]=Scalars[25+i];\n"
         "    for(i=0;i<5;i++) {\n"
         "      Atilde[i]=Scalars[20+i]^Mtot[i];      // A~ = A + Mtot   (hashed, not sent)\n"
         "      MaskAtilde[i]=Atilde[i];\n"
         "      Scalars[25+i]^=Offsets[mOFFSET+i];    // B~ = B + W      (sent)\n"
         "    }\n"
         "    SHA256Update(hash_ptr, (uchar *)Atilde,5*sizeof(u64));\n"
         "#ifdef DUMPIDENT\n"
         "    if (getenv(\"DUMPIDENT\")) {\n"
         "      FILE *f=fopen(getenv(\"DUMPIDENT\"),\"a\"); int q;\n"
         "      fprintf(f,\"S_A\");    for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Scalars[20+q]); fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"S_BRAW\"); for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Braw[q]); fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"S_AT\");   for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Atilde[q]);    fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"S_BT\");   for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Scalars[25+q]); fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"S_MTOT\"); for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Mtot[q]);      fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"S_W\");    for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Offsets[mOFFSET+q]); fprintf(f,\"\\n\");\n"
         "      fclose(f);\n"
         "    }\n"
         "#endif\n"
         "  }\n"
         "#else\n",
         "signer-maskAB"),
        (708, 708, "#endif\n", "signer-maskAB-end"),
        (718, 719, "", "verifier-decl-gap"),
        (720, 720,
         "#ifdef MASKB\n"
         "  u64 MOpenRound[5], Qtot[5];\n"
         "  for(i=0;i<5;i++) Qtot[i]=0;\n"
         "#endif\n\n",
         "verifier-decl"),
        (729, 729,
         "#ifdef MASKB\n"
         "    for(int pos=0;pos<5;pos++) MOpenRound[pos]=0;\n"
         "#endif\n",
         "verifier-roundinit"),
        (750, 750,
         "#ifdef MASKB\n"
         "      FFMulAdd(MOpenRound,HyperMemorys+PLAYER_SHARE_MEM*(2*counter+good)+mOFFSET,TDeltavec+5*Lambda+5*counter);\n"
         "#endif\n",
         "verifier-ffmuladd"),
        (764, 764,
         "#ifdef MASKB\n"
         "    for(i=0;i<5;i++) Qtot[i]^=MOpenRound[i];  /* verifier: accumulate Qtot */\n"
         "#endif\n",
         "verifier-qtot"),
        (785, 785,
         "#ifdef DUMPIDENT\n"
         "  if (getenv(\"DUMPIDENT\")) {\n"
         "    FILE *f=fopen(getenv(\"DUMPIDENT\"),\"a\"); int q;\n"
         "    fprintf(f,\"V_PHIOPEN\"); for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Scalars[20+q]); fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"V_GDEL\");    for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)GlobDelta[q]);   fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"V_BT\");      for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Scalars[25+q]);  fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"V_QTOT\");    for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Qtot[q]);        fprintf(f,\"\\n\");\n"
         "    fclose(f);\n"
         "  }\n"
         "#endif\n",
         "verifier-dumpident"),
        (786, 786,
         "#ifdef MASKB\n"
         "  // A'' = Phi_open + GlobDelta*B~ + Qtot  ==  A + Mtot  ==  A~\n"
         "  for(i=0;i<5;i++) Scalars[20+i]^=Qtot[i];\n"
         "#endif\n",
         "verifier-addQ"),
        (789, 789,
         "#ifdef DUMPATTACK\n"
         "  if (getenv(\"DUMPATTACK\")) {\n"
         "    FILE *f=fopen(getenv(\"DUMPATTACK\"),\"a\");\n"
         "    int q;\n"
         "    fprintf(f,\"DELTA\"); for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)GlobDelta[q]); fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"B\");     for(q=0;q<5;q++) fprintf(f,\" %016llx\",(unsigned long long)Scalars[25+q]); fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"FX\");    for(q=0;q<5*Lambda;q++) fprintf(f,\" %016llx\",(unsigned long long)XOpen[q]); fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"FY\");    for(q=0;q<5*Lambda;q++) fprintf(f,\" %016llx\",(unsigned long long)YOpen[q]); fprintf(f,\"\\n\");\n"
         "    fprintf(f,\"UV\");    for(q=0;q<2*130*5;q++) fprintf(f,\" %016llx\",(unsigned long long)UV[q]); fprintf(f,\"\\n\");\n"
         "    fclose(f);\n"
         "  }\n"
         "#endif\n",
         "dumpattack-verifier"),
        (900, 900, "#ifndef MASKB\n", "signer-forcequal-guard"),
        (901, 901, "#endif\n", "signer-forcequal-guard-end"),
        (903, 903,
         "\n#ifdef MASKB\n"
         "  // Force equality of the mask values across trees: dw_r := raw_r + raw_0\n"
         "  for(r=1;r<NbRounds;r++) {\n"
         "      for(i=0;i<5;i++) Offsets[r*PLAYER_OFFSET_MEM+mOFFSET+i]^=Offsets[mOFFSET+i];\n"
         "    }\n"
         "#endif\n",
         "signer-forcequal"),
        (924, 924, h0_block, "hash-dw-into-h0-sign"),
        (954, 954,
         "#if defined(FORGE_NOH0) && defined(MASKB)\n"
         "#include \"forge_noh0.inc\"\n"
         "#endif\n\n",
         "forge-include"),
        (1013, 1013,
         "#ifdef MASKB\n"
         "    for(r=1;r<NbRounds;r++) {\n"
         "      F257ToSignature(Offsets+r*PLAYER_OFFSET_MEM+mOFFSET, 1, CompSig, sigsize, &leftoverpos, &leftovernum);\n"
         "    }\n"
         "#endif\n",
         "signer-serialise"),
        (1154, 1154,
         "#ifdef MASKB\n"
         "    for(r=1;r<NbRounds;r++) {\n"
         "      F257FromSignature(Offsets+r*PLAYER_OFFSET_MEM+mOFFSET, 1,\n"
         "\t\t\tCompSig, &sigoffset, &leftoverpos, &leftovernum);\n"
         "    }\n"
         "#endif\n",
         "verifier-deserialise"),
        (1171, 1171,
         "#ifdef MASKB\n"
         "\tfor(pos=0;pos<5;pos++) {\n"
         "\t  HyperMemorys[PLAYER_SHARE_MEM*(2*(r*LogPlayers+bitpos)+1)+mOFFSET+pos]^=\n"
         "\t    Offsets[r*PLAYER_OFFSET_MEM+mOFFSET+pos];\n"
         "\t}\n"
         "#endif\n",
         "verifier-patch"),
        (1209, 1209, h0_block, "hash-dw-into-h0-verify"),
        (1296, 1296,
         "#ifdef DUMPATTACK\n"
         "    if (getenv(\"DUMPATTACK\")) {\n"
         "      FILE *f=fopen(getenv(\"DUMPATTACK\"),\"a\"); int q;\n"
         "      fprintf(f,\"SECX\"); for(q=0;q<4;q++) fprintf(f,\" %016llx\",(unsigned long long)SecX[q]); fprintf(f,\"\\n\");\n"
         "      fprintf(f,\"UVXY\"); for(q=0;q<20;q++) fprintf(f,\" %016llx\",(unsigned long long)UVxy[q]); fprintf(f,\"\\n\");\n"
         "      fclose(f);\n"
         "    }\n"
         "#endif\n"
         "#if defined(FORGE_NOH0) || defined(FALSEWITNESS)\n"
         "    // FALSE WITNESS.  Keep the victim's public key (UV, UVxy) and throw the\n"
         "    // secret away: (x,y) no longer solves NSBC[u,v], so Phi(Delta) is genuinely\n"
         "    // quadratic and the honest check cannot pass.  Everything the attacker does\n"
         "    // from here uses the public key only.\n"
         "    if (getenv(\"FORGE_NOH0\")) {\n"
         "      SecX[0]^=0x0123456789abcdefUL; SecX[1]^=0xfedcba9876543210UL;\n"
         "      SecX[2]^=0x5555555555555555UL; SecX[3]^=0xaaaaaaaaaaaaaaaaUL;\n"
         "      fprintf(stderr,\"FORGE_NOH0: signing with a FALSE witness against the victim public key\\n\");\n"
         "    }\n"
         "#endif\n",
         "dumpattack-keygen"),
    ]
    return apply_line_ops(src, ops, EXPECTED_PAT_C_SHA, "SBCvole.c")


def edits_header(src: str):
    ops = [
        (25, 25,
         "#ifdef MASKB\n"
         "// one more 257-bit row per leaf: 2*NBOFF128 + 4*NBOFF257 + 1 <= 2*EXPAND_RATIO\n"
         "#ifdef SINGLECOMMITSIZE\n"
         "#define EXPAND_RATIO 8\n"
         "#else\n"
         "#define EXPAND_RATIO 9\n"
         "#endif\n"
         "#else\n",
         "expand_ratio"),
        (29, 29, "#endif\n", "expand_ratio_end"),
        (40, 40,
         "#ifdef MASKB\n"
         "#define NBOFF257 2\n"
         "#else\n",
         "nboff257"),
        (41, 41, "#endif\n", "nboff257_end"),
        (51, 51,
         "#ifdef MASKB\n"
         "#define PLAYER_SHARE_MEM 21\n"
         "#define PLAYER_OFFSET_MEM 21\n"
         "#define mOFFSET (zOFFSET+5)\n"
         "#else\n",
         "share_mem"),
        (53, 53, "#endif\n", "share_mem_end"),
    ]
    out, _ = apply_line_ops(src, ops, EXPECTED_PAT_H_SHA, "SBCunified.h")
    return out


def ensure_vendor_zip():
    if not os.path.isfile(ZIP):
        subprocess.run(["bash", FETCH_SCRIPT, "--ensure"], check=True)


def main(dest):
    ensure_vendor_zip()
    with open(ZIP, "rb") as fh:
        blob = fh.read()
    got = hashlib.sha256(blob).hexdigest()
    assert got == ZIP_SHA, "vendor zip sha256 mismatch: %s" % got

    stage = os.path.join(dest, "_stage")
    if os.path.isdir(stage):
        shutil.rmtree(stage)
    os.makedirs(stage)
    try:
        zipfile.ZipFile(ZIP).extractall(stage)
        srcdir = os.path.join(stage, "ArtifactSBCVOLE")

        with open(os.path.join(srcdir, "SBCvole.c")) as fh:
            c, log = edits_sbcvole(fh.read())
        with open(os.path.join(srcdir, "SBCunified.h")) as fh:
            h = edits_header(fh.read())

        os.makedirs(dest, exist_ok=True)
        for name in ("SBCvole.c", "SBCunified.h", "SBCtools.c", "SBCunified-CPUspecific.h"):
            with open(os.path.join(dest, name), "w") as fh:
                fh.write(c if name == "SBCvole.c" else h if name == "SBCunified.h"
                         else open(os.path.join(srcdir, name)).read())
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    print("patched %s: %d verified edits" % (os.path.basename(os.path.normpath(dest)), len(log)))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "masked")
