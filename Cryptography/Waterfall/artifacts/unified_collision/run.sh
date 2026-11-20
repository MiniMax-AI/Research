#!/bin/bash
set -e
cd "$(dirname "$0")"
MODE="${1:---demo}"
case "$MODE" in
  --demo|--full) ;;
  -h|--help)
    echo "Usage: ./run.sh [--demo | --full]"
    exit 0
    ;;
  *)
    echo "Usage: ./run.sh [--demo | --full]" >&2
    exit 1
    ;;
esac

bash ../fetch_third_party.sh --ensure

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [ -z "${CC:-}" ]; then
  if command -v gcc >/dev/null 2>&1; then
    CC=gcc
  else
    CC=cc
  fi
fi

ARCH_FLAGS=""
if printf 'int main(void){return 0;}\n' | "$CC" -O3 -march=native -x c - -o "$TMP/arch_probe" >/dev/null 2>&1; then
  ARCH_FLAGS="-march=native"
fi

COMPAT_FLAGS=""
if ! printf '#include <stdint.h>\n#include <immintrin.h>\nint main(void){return (int)(_pext_u32(5u,3u)+_pdep_u32(5u,3u));}\n' \
     | "$CC" -O3 $ARCH_FLAGS -x c - -o "$TMP/bmi2_probe" >/dev/null 2>&1; then
  mkdir -p "$TMP/compat"
  cat > "$TMP/compat/immintrin.h" << 'EOF'
#ifndef PORTABLE_IMMINTRIN_COMPAT_H
#define PORTABLE_IMMINTRIN_COMPAT_H
#include <stdint.h>
static inline uint32_t _pext_u32(uint32_t val, uint32_t mask) {
    uint32_t res = 0;
    for (uint32_t bb = 1; mask; bb <<= 1) {
        if (val & (mask & (0u - mask))) res |= bb;
        mask &= mask - 1;
    }
    return res;
}
static inline uint32_t _pdep_u32(uint32_t val, uint32_t mask) {
    uint32_t res = 0;
    for (uint32_t bb = 1; mask; bb <<= 1) {
        if (val & bb) res |= mask & (0u - mask);
        mask &= mask - 1;
    }
    return res;
}
#endif
EOF
  COMPAT_FLAGS="-I$TMP/compat"
fi

echo "=================================================================="
echo "Unified Minimal-Length Full-State Collision Verification (Waterfall)"
echo "Mode: $MODE"
echo "=================================================================="

# 1. Build unmodified official NIST reference implementation & independent verifier
"$CC" -O2 -o "$TMP/nist_ver" reference/nist_driver.c reference/Waterfall_ref.c
"$CC" -O2 -o "$TMP/stream_ver" reference/verify_reference.c

# 2. Verify shipped 448-byte (7-block) and 512-byte (8-block) Waterfall-256 collision witnesses
echo "[1] Verifying pre-computed 448-byte (7-block) & 512-byte (8-block) Waterfall-256 witnesses against NIST C reference:"
for sz in 448 512; do
  m1="witnesses/collide_${sz}B_m1.bin"
  m2="witnesses/collide_${sz}B_m2.bin"
  cmp -s "$m1" "$m2" && { echo "FAIL: messages identical"; exit 1; }
  d1=$("$TMP/nist_ver" "$m1" 256)
  d2=$("$TMP/nist_ver" "$m2" 256)
  d_stream=$("$TMP/stream_ver" "$m1")
  [ "$d1" = "$d2" ] && [ "$d1" = "$d_stream" ] || { echo "FAIL: digest mismatch on ${sz}B Waterfall-256"; exit 1; }
  echo "  PASS  ${sz}-byte ($((sz/64))-block) witness verified by NIST C & stream verifier: ${d1:0:32}..."
done

# 3. Verify the dual-excursion and quiet-counter coset-scan ablation solvers
echo "[2] Compiling and verifying ablation collision solvers (solve5.c & attack7.c):"
"$CC" -O3 $ARCH_FLAGS -o "$TMP/ablation_solve5" ablation_solvers/route_a/solve5.c -lm
"$TMP/ablation_solve5" 1 > "$TMP/outA.txt" 2>/dev/null
python3 -c '
import sys
tmp = sys.argv[1]
ls = open(f"{tmp}/outA.txt").read().strip().splitlines()
b1, b2 = bytes.fromhex(ls[0]), bytes.fromhex(ls[1])
assert len(b1) == 2816 and len(b2) == 2816 and b1 != b2
open(f"{tmp}/mA1.bin", "wb").write(b1)
open(f"{tmp}/mA2.bin", "wb").write(b2)
' "$TMP"
dA1=$("$TMP/nist_ver" "$TMP/mA1.bin" 256)
dA2=$("$TMP/nist_ver" "$TMP/mA2.bin" 256)
[ "$dA1" = "$dA2" ] || { echo "FAIL: dual-excursion solve5.c"; exit 1; }
echo "  PASS  Dual-excursion ablation solver (solve5.c):       44-block (2,816 B) collision verified: ${dA1:0:32}..."

"$CC" -O3 $ARCH_FLAGS $COMPAT_FLAGS -o "$TMP/ablation_attack7" ablation_solvers/route_b/attack7.c
"$TMP/ablation_attack7" > "$TMP/outB.txt" 2>/dev/null
python3 -c '
import sys
tmp = sys.argv[1]
ls = open(f"{tmp}/outB.txt").read().strip().splitlines()
b1, b2 = bytes.fromhex(ls[0]), bytes.fromhex(ls[1])
assert len(b1) == 576 and len(b2) == 576 and b1 != b2
open(f"{tmp}/mB1.bin", "wb").write(b1)
open(f"{tmp}/mB2.bin", "wb").write(b2)
' "$TMP"
dB1=$("$TMP/nist_ver" "$TMP/mB1.bin" 256)
dB2=$("$TMP/nist_ver" "$TMP/mB2.bin" 256)
[ "$dB1" = "$dB2" ] || { echo "FAIL: quiet-counter attack7.c"; exit 1; }
echo "  PASS  Quiet-counter coset-scan solver (attack7.c):      9-block (  576 B) collision verified: ${dB1:0:32}..."

# 4. Build and run unified live collision generator for 9, 8, and 7 blocks
echo "[3] Compiling and running unified minimal-length collision solvers (9, 8, and 7 blocks):"
"$CC" -O3 $ARCH_FLAGS $COMPAT_FLAGS -DPATTERN_FILE="\"pattern_8blk.h\"" -o "$TMP/ucol8" unified_collide.c
"$CC" -O3 $ARCH_FLAGS $COMPAT_FLAGS -DPATTERN_FILE="\"pattern_7blk.h\"" -o "$TMP/ucol7" unified_collide.c
"$CC" -O3 $ARCH_FLAGS $COMPAT_FLAGS -DPATTERN_FILE="\"pattern_9blk.h\"" -o "$TMP/ucol9" unified_collide.c

"$TMP/ucol9" 0 256 > "$TMP/out9.txt" 2>/dev/null
python3 -c '
import sys
tmp = sys.argv[1]
ls = open(f"{tmp}/out9.txt").read().strip().splitlines()
b1, b2 = bytes.fromhex(ls[0]), bytes.fromhex(ls[1])
assert len(b1) == 576 and len(b2) == 576 and b1 != b2
open(f"{tmp}/m1_9.bin", "wb").write(b1)
open(f"{tmp}/m2_9.bin", "wb").write(b2)
' "$TMP"
d9_1=$("$TMP/nist_ver" "$TMP/m1_9.bin" 256)
d9_2=$("$TMP/nist_ver" "$TMP/m2_9.bin" 256)
[ "$d9_1" = "$d9_2" ] || { echo "FAIL: live 576B Waterfall-256"; exit 1; }
echo "  PASS  Live 9-block (576 B) solver verified for Waterfall-256: ${d9_1:0:32}..."

for hbits in 224 256 384 512; do
  "$TMP/ucol8" 0 "$hbits" > "$TMP/out8_${hbits}.txt" 2>/dev/null
  python3 -c '
import sys
tmp, hbits = sys.argv[1], sys.argv[2]
ls = open(f"{tmp}/out8_{hbits}.txt").read().strip().splitlines()
b1, b2 = bytes.fromhex(ls[0]), bytes.fromhex(ls[1])
assert len(b1) == 512 and len(b2) == 512 and b1 != b2
open(f"{tmp}/m1.bin", "wb").write(b1)
open(f"{tmp}/m2.bin", "wb").write(b2)
' "$TMP" "$hbits"
  d1=$("$TMP/nist_ver" "$TMP/m1.bin" "$hbits")
  d2=$("$TMP/nist_ver" "$TMP/m2.bin" "$hbits")
  [ "$d1" = "$d2" ] || { echo "FAIL: live 512B Waterfall-${hbits}"; exit 1; }
  echo "  PASS  Live 8-block (512 B) solver verified for Waterfall-${hbits}: ${d1:0:32}..."
done

"$TMP/ucol7" 0 256 > "$TMP/out7.txt" 2>/dev/null
python3 -c '
import sys
tmp = sys.argv[1]
ls = open(f"{tmp}/out7.txt").read().strip().splitlines()
b1, b2 = bytes.fromhex(ls[0]), bytes.fromhex(ls[1])
assert len(b1) == 448 and len(b2) == 448 and b1 != b2
open(f"{tmp}/m1_7.bin", "wb").write(b1)
open(f"{tmp}/m2_7.bin", "wb").write(b2)
' "$TMP"
d7_1=$("$TMP/nist_ver" "$TMP/m1_7.bin" 256)
d7_2=$("$TMP/nist_ver" "$TMP/m2_7.bin" 256)
[ "$d7_1" = "$d7_2" ] || { echo "FAIL: live 448B Waterfall-256"; exit 1; }
echo "  PASS  Live 7-block (448 B) minimal-length solver verified for Waterfall-256: ${d7_1:0:32}..."

if [ "$MODE" = "--full" ]; then
  echo "[4] Running Theorem 2 SAT lower-bound check and 20-seed multi-collision sweep:"
  python3 gen_patterns.py
  for s in $(seq 0 19); do
    "$TMP/ucol8" $((s * 10)) 256 > "$TMP/sweep.txt" 2>/dev/null
    python3 -c '
import sys
tmp = sys.argv[1]
ls = open(f"{tmp}/sweep.txt").read().strip().splitlines()
open(f"{tmp}/s1.bin", "wb").write(bytes.fromhex(ls[0]))
open(f"{tmp}/s2.bin", "wb").write(bytes.fromhex(ls[1]))
assert ls[0] != ls[1]
' "$TMP"
    ds1=$("$TMP/nist_ver" "$TMP/s1.bin" 256)
    ds2=$("$TMP/nist_ver" "$TMP/s2.bin" 256)
    [ "$ds1" = "$ds2" ] || { echo "FAIL on seed $s"; exit 1; }
  done
  echo "  PASS  20/20 independent seeds produced verified 512-byte Waterfall-256 collisions."
fi

echo "RESULT: SUCCESS"
