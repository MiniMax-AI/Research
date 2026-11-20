#!/bin/sh
# Build a ground-truth instance generator for one KEMSET.
# usage: build_gen.sh <ROMAN>  ->  ciphertext_attacks/parameterized_tsm/gen<ROMAN>
set -eu
k=$1
CXX="${CXX:-g++}"
HERE="$(cd "$(dirname "$0")" && pwd)"
bash "$HERE/../../fetch_third_party.sh" --ensure
cd "$HERE/../../reference/native"
SRC="Reference_Implementation/KEMSET-$k"
EXTRA_INC="${NTL_PREFIX:+-I$NTL_PREFIX/include}"
EXTRA_LIB="${NTL_PREFIX:+-L$NTL_PREFIX/lib -Wl,-rpath,$NTL_PREFIX/lib}"
if [ -z "$EXTRA_INC" ] && [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/include/NTL/ZZ.h" ]; then
  EXTRA_INC="-I$CONDA_PREFIX/include"
  EXTRA_LIB="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
fi
for pfx in /opt/homebrew /usr/local /opt/homebrew/opt/openssl /usr/local/opt/openssl /opt/homebrew/opt/gmp /usr/local/opt/gmp; do
  if [ -d "$pfx/include" ] && [ -d "$pfx/lib" ]; then
    EXTRA_INC="$EXTRA_INC -I$pfx/include"
    EXTRA_LIB="$EXTRA_LIB -L$pfx/lib"
  fi
done
INC="-I. -I$SRC/src -I$SRC/src/ffi -I$SRC/lib/hash -I$SRC/lib/rng $EXTRA_INC"
O="$(mktemp -d)"
trap 'rm -rf "$O"' EXIT
for f in kem parsing decoder; do "$CXX" -O2 -fcommon -c "$SRC/src/$f.cpp" $INC -o "$O/$f.o"; done
for f in ffi_field ffi_elt ffi_vec; do "$CXX" -O2 -fcommon -c "$SRC/src/ffi/$f.cpp" $INC -o "$O/$f.o"; done
for f in hash rng; do "$CXX" -O2 -fcommon -c "$SRC/lib/$f/$f.c" $INC -o "$O/$f.o"; done
NTL_LIB=""
if [ -n "${NTL_PREFIX:-}" ] && [ -f "$NTL_PREFIX/lib/libntl.a" ]; then
  NTL_LIB="$NTL_PREFIX/lib/libntl.a"
fi
"$CXX" -O2 -fcommon -o "$HERE/gen$k" "$HERE/gen_instance.cpp" $INC \
    "$O/kem.o" "$O/parsing.o" "$O/decoder.o" "$O/ffi_field.o" "$O/ffi_elt.o" "$O/ffi_vec.o" "$O/hash.o" "$O/rng.o" \
    $EXTRA_LIB ${NTL_LIB:--lntl} -lgmp -lssl -lcrypto
echo "built $HERE/gen$k"
