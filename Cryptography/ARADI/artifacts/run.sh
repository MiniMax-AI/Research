#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# =====================================================================
# Single-Core Verification & Experimental Suite for NSA ARADI (5R–12R)
#
#   ./run.sh --demo         fast end-to-end verification (~30 s, 1 core; default)
#   ./run.sh --higher-only  fast higher-round (6R dist + 9R–12R) check (~6 s, 1 core)
#   ./run.sh --full         full statistical suite across 5R–9R + 10R–12R audit
# =====================================================================

export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

MODE="${1:---demo}"
SEED="${SEED:-20260409}"
PYTHON="${PYTHON:-python3}"
if [ -z "${CC:-}" ]; then
    if command -v gcc >/dev/null 2>&1; then
        CC=gcc
    else
        CC=cc
    fi
fi
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/aradi.XXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT

CFLAGS_OPT=(-std=c99 -O3 -Wall -Wextra)
if printf 'int main(void){return 0;}\n' | "$CC" -x c -std=c99 -O3 -march=native - -o "$SCRATCH/probe" >/dev/null 2>&1; then
    CFLAGS_OPT+=(-march=native)
elif printf 'int main(void){return 0;}\n' | "$CC" -x c -std=c99 -O3 -mcpu=native - -o "$SCRATCH/probe" >/dev/null 2>&1; then
    CFLAGS_OPT+=(-mcpu=native)
fi
rm -f "$SCRATCH/probe"

if [ "$MODE" = "--full" ]; then
    TRIALS_6CP=200
    TRIALS_4CP=10
    TRIALS_6R=10
    TRIALS_7R=10
    TRIALS_8R=10
    TRIALS_9R=20
elif [ "$MODE" = "--demo" ]; then
    TRIALS_6CP=25
    TRIALS_4CP=2
    TRIALS_6R=2
    TRIALS_7R=2
    TRIALS_8R=2
    TRIALS_9R=2
elif [ "$MODE" = "--higher-only" ]; then
    TRIALS_9R=2
else
    echo "Usage: $0 [--demo | --higher-only | --full]" >&2
    exit 2
fi

FAIL=0

if [ "$MODE" != "--higher-only" ]; then
    echo "=== [1] Official NSA Test Vectors, Orthogonal Island Proof & LATINCRYPT 2025 Audit ==="
    "$PYTHON" verify_algebra_and_audit.py || FAIL=1
    echo

    echo "=== [2] 5-Round ARADI Ultra-Low-Data Key Recovery (${TRIALS_6CP}x 6-CP MITM & ${TRIALS_4CP}x 4-CP 256-Bit MK) ==="
    "$PYTHON" attack_5r_lowdata.py --trials-6cp "$TRIALS_6CP" --trials-4cp "$TRIALS_4CP" --seed "$SEED" || FAIL=1
    echo

    echo "=== [3] 6-Round ARADI 48-CP Full 256-Bit Master-Key Recovery (${TRIALS_6R} trials) ==="
    "$PYTHON" attack_6r_48cp.py --trials "$TRIALS_6R" --seed "$SEED" || FAIL=1
    echo

    echo "=== [4] 7-Round ARADI 3,072-CP Full 256-Bit Master-Key Recovery (${TRIALS_7R} trials) ==="
    "$PYTHON" attack_7r_fht.py --trials "$TRIALS_7R" --seed "$SEED" || FAIL=1
    echo

    echo "=== [5] 8-Round ARADI 90,112-CP Full 256-Bit Master-Key Recovery (${TRIALS_8R} trials, 1-core C) ==="
    "$CC" "${CFLAGS_OPT[@]}" attack_8r_mitm.c -o "$SCRATCH/attack_8r_mitm" -lm
    "$SCRATCH/attack_8r_mitm" --trials "$TRIALS_8R" --seed "$SEED" || FAIL=1
    echo
fi

echo "=== [6] Higher-Round Algebraic Theorems, 6R d=16 Distinguisher & 9R–12R Progressive List-Join Audit ==="
"$PYTHON" verify_higher_rounds.py || FAIL=1
echo

echo "=== [7] 9-Round ARADI Nested Joint-WHT Peel Verification (${TRIALS_9R} trials, 1-core C) ==="
"$CC" "${CFLAGS_OPT[@]}" attack_9r_nested_fht.c -o "$SCRATCH/attack_9r_nested_fht" -lm
"$SCRATCH/attack_9r_nested_fht" --trials "$TRIALS_9R" --kbits 0 --cubes 5 --shifts 16 || FAIL=1
echo

if [ "$FAIL" -eq 0 ]; then
    echo "RESULT: SUCCESS"
    exit 0
else
    echo "RESULT: FAIL"
    exit 1
fi
