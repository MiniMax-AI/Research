#!/bin/bash
# verify_s16.sh -- full reproduction of the sparse Block-Wiedemann / Lanczos engine and cost model
# Runtime: ~6 min on 8 cores (excluding the optional LOCKER-I solve, ~3 min).
set -e
case "$(uname -m)" in
  x86_64|amd64) ;;
  *)
    echo "ERROR: verify_s16.sh requires an x86_64 CPU with PCLMULQDQ/SSE2 support (detected: $(uname -m))." >&2
    exit 1
    ;;
esac
HERE="$(cd "$(dirname "$0")" && pwd)"
CLEAN_B=0
if [ -z "${BUILD:-}" ]; then
  B="$(mktemp -d)"
  CLEAN_B=1
else
  B="$BUILD"; mkdir -p "$B"
fi
trap '[ "$CLEAN_B" = 1 ] && rm -rf "$B"' EXIT
export B
CXX=${CXX:-g++}
PY=${PY:-python3}
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
cd "$HERE"

echo "### 0. build"
"$CXX" -O3 -march=native -mpclmul -std=c++17 -funroll-loops          -o "$B/test_gf2n" test_gf2n.cpp
"$CXX" -O3 -march=native -mpclmul -std=c++17 -fopenmp -funroll-loops -o "$B/tsm_bw"    tsm_bw.cpp

echo "### 1. GF(2^n) PCLMULQDQ kernel: exactness vs. bit-level reference + inverse + timing"
"$B/test_gf2n"

echo "### 2. POSITIVE CONTROL: genuine ciphertext -> exact kernel vector"
for sd in 7 11 23; do
  $PY gen_instance.py --m 7 --n 11 --r 3 --d 3 --P 5 --gmod 7,1,0 --Pmod 11,2,0 --out "$B/i.bin" --seed $sd >/dev/null
  "$B/tsm_bw" build "$B/i.bin" "$B/m.mat" 7,1,0 11,2,0 2>/dev/null
  "$B/tsm_bw" solve "$B/m.mat" 11,2,0 11 2>/dev/null
done

echo "### 3. POSITIVE CONTROL, SQUARE SYSTEM (m=20,n=23,r=3,P=5): 80x80, genuine u"
$PY gen_instance.py --m 20 --n 23 --r 3 --d 4 --P 5 --gmod 20,3,0 --Pmod 23,5,0 --out "$B/isq.bin" --seed 501 >/dev/null
"$B/tsm_bw" build "$B/isq.bin" "$B/msq.mat" 20,3,0 23,5,0 2>/dev/null
"$B/tsm_bw" solve "$B/msq.mat" 23,5,0 23 2>/dev/null

echo "### 4. NEGATIVE CONTROL: uniform random u -> 0 solutions"
for sd in 601 602 603; do
  $PY gen_instance.py --m 20 --n 23 --r 3 --d 4 --P 5 --gmod 20,3,0 --Pmod 23,5,0 --out "$B/ineg.bin" --seed $sd --random-u >/dev/null
  "$B/tsm_bw" build "$B/ineg.bin" "$B/mneg.mat" 20,3,0 23,5,0 >/dev/null 2>&1
  "$B/tsm_bw" solve "$B/mneg.mat" 23,5,0 23 2>/dev/null || true   # exit 2 is the EXPECTED rejection
done

echo "### 5. scaled TSM-d1 (r=6) and TSM-d2 (m=15,r=4) calibration instances"
$PY gen_instance.py --m 15 --n 17 --r 4 --d 4 --P 7 --gmod 15,1,0 --Pmod 17,3,0 --out "$B/i15.bin" --seed 202 >/dev/null
"$B/tsm_bw" build "$B/i15.bin" "$B/m15.mat" 15,1,0 17,3,0 2>/dev/null
"$B/tsm_bw" solve "$B/m15.mat" 17,3,0 17 2>/dev/null
$PY gen_instance.py --m 67 --n 83 --r 6 --d 8 --P 9 --gmod 67,5,2,1,0 --Pmod 83,7,4,2,0 --out "$B/i67.bin" --seed 303 >/dev/null
"$B/tsm_bw" build "$B/i67.bin" "$B/m67.mat" 67,5,2,1,0 83,7,4,2,0 2>/dev/null
"$B/tsm_bw" solve "$B/m67.mat" 83,7,4,2,0 83 2>/dev/null

echo "### 6. REAL LOCKER-I (m=71,n=83,r=5,d=7,P=9) end-to-end  [~3 min]"
$PY gen_instance.py --m 71 --n 83 --r 5 --d 7 --P 9 --gmod 71,6,0 --Pmod 83,7,4,2,0 --out "$B/iL1.bin" --seed 101 >/dev/null
"$B/tsm_bw" build "$B/iL1.bin" "$B/mL1.mat" 71,6,0 83,7,4,2,0 2>/dev/null
"$B/tsm_bw" solve "$B/mL1.mat" 83,7,4,2,0 83 2>/dev/null

echo "### 7. REAL LOCKER-III-SCALE MATRIX (m=107,n=89,r=6,P=14): 183612 x 192192, nnz 82258176"
$PY gen_instance.py --m 107 --n 89 --r 6 --d 8 --P 14 --gmod 107,9,7,4,0 --Pmod 89,6,5,3,0 --out "$B/iL3.bin" --seed 404 >/dev/null
"$B/tsm_bw" build "$B/iL3.bin" "$B/mL3.mat" 107,9,7,4,0 89,6,5,3,0
for beta in 1 2 4; do "$B/tsm_bw" bench "$B/mL3.mat" 89,6,5,3,0 89 $beta 3; done
rm -f "$B/mL3.mat"

echo "### 8. exact non-asymptotic cost tables"
$PY cost_model.py > "$B/cost.json"
$PY - <<'PYEOF'
import json,os
rows=json.load(open(os.environ['B']+'/cost.json'))
for r in rows:
    if r['P_min'] is None: continue
    print('%-22s 2^%-4d claimed -> 2^%.1f bit-ops (Lanczos)  margin %3.0f bits  sparse %.2f GB  dense %.0f GB'
          % (r['name'],r['claimed_log2'],r['log2_bitops_block_lanczos'],
             r['claimed_log2']-r['log2_bitops_block_lanczos'],
             r['sparse_bytes']/1e9, r['dense_bytes']/1e9))
PYEOF
echo "### ALL SPARSE KRYLOV CHECKS COMPLETE"
