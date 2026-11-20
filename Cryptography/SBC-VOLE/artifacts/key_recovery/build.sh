#!/bin/bash
# build.sh -- unpack and build the AUTHORS' OWN reference implementation
#             with ZERO modifications to SBCvole.c and SBCtools.c.
#
#   ./build.sh <workdir> <D> <basic|minmax>
#
# Produces, in <workdir>/work:
#   pristine    the artifact's own main(), built from the 100% unmodified source
#   driver      forge.c linked against the 100% unmodified SBCvole.c (with
#               -Dmain=sbcvole_artifact_main_unused so forge.c provides main())
set -e
W="$1"; D="$2"; MM="$3"
HERE="$(cd "$(dirname "$0")" && pwd)"
[ "$MM" = "minmax" ] && MFLAG=-DMINMAXCOMMIT || MFLAG=

bash "$HERE/../fetch_third_party.sh" --ensure

if [ ! -f "$W/work/SBCvole.c" ]; then
  mkdir -p "$W"
  python3 -c "import zipfile,sys;zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" \
      "$HERE/vendor/asiacrypt-2025-a10.zip" "$W"
  mkdir -p "$W/crypto-algorithms"
  cp "$HERE"/vendor/crypto-algorithms/* "$W/crypto-algorithms/"
  cp -r "$W/ArtifactSBCVOLE" "$W/work"
  cp "$HERE/forge.c" "$W/work/"
  # Verify that SBCvole.c and SBCtools.c are byte-for-byte identical to the zip
  cmp -s "$W/ArtifactSBCVOLE/SBCvole.c" "$W/work/SBCvole.c"
  cmp -s "$W/ArtifactSBCVOLE/SBCtools.c" "$W/work/SBCtools.c"
fi

cd "$W/work"
CA="$W/crypto-algorithms"
CC="${CC:-$(command -v gcc >/dev/null 2>&1 && echo gcc || echo cc)}"
CFLAGS="-O2 -DL$D $MFLAG -DNBREPEATS=1 -DGENERICAES -DGENERICSHA"
# 1. the artifact's own main(), 100% unmodified source
"$CC" $CFLAGS -DSILENT -o pristine \
    SBCvole.c SBCtools.c "$CA/aes.c" "$CA/sha256.c" -lm
# 2. forge.c driver linked against 100% unmodified SBCvole.c
"$CC" $CFLAGS -Dmain=sbcvole_artifact_main_unused -c SBCvole.c -o nomain.o
"$CC" $CFLAGS -o driver forge.c nomain.o SBCtools.c \
    "$CA/aes.c" "$CA/sha256.c" -lm
echo "built driver and pristine (zero-patch unmodified SBCvole.c) for D=$D $MM in $W/work"
