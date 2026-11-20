#!/usr/bin/env bash
set -euo pipefail

KAT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$KAT_DIR/hawk_shim}"

bash "$KAT_DIR/../fetch_third_party.sh" --ensure

if [ -n "${CC:-}" ]; then
    COMPILER="$CC"
elif command -v cc >/dev/null 2>&1; then
    COMPILER="cc"
elif command -v gcc >/dev/null 2>&1; then
    COMPILER="gcc"
elif command -v clang >/dev/null 2>&1; then
    COMPILER="clang"
else
    echo "ERROR: No C compiler found (set CC or install gcc/clang)." >&2
    exit 1
fi

"$COMPILER" -O2 -I"$KAT_DIR/hawk_ref" -o "$OUT" \
    "$KAT_DIR/hawk_shim.c" \
    "$KAT_DIR/hawk_ref/hawk_sign.c" \
    "$KAT_DIR/hawk_ref/hawk_vrfy.c" \
    "$KAT_DIR/hawk_ref/hawk_kgen.c" \
    "$KAT_DIR/hawk_ref/ng_fxp.c" \
    "$KAT_DIR/hawk_ref/ng_hawk.c" \
    "$KAT_DIR/hawk_ref/ng_mp31.c" \
    "$KAT_DIR/hawk_ref/ng_ntru.c" \
    "$KAT_DIR/hawk_ref/ng_poly.c" \
    "$KAT_DIR/hawk_ref/ng_zint31.c" \
    "$KAT_DIR/hawk_ref/sha3.c" \
    -lm
