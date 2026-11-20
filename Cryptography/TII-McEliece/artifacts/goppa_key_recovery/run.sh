#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export TMPDIR="$TMP"
export PYTHONDONTWRITEBYTECODE=1

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

echo "=============================================================================="
echo "TII McEliece (Binary Goppa Code) Key-Recovery Artifact Suite ($MODE)"
echo "=============================================================================="

# Compile native C/C++ helpers into scratch directory $TMP (no writes inside repo).
# Prefer OpenMP-capable GNU gcc/g++ (including Homebrew gcc-*/g++-* on macOS when
# CC/CXX are unset), and fall back cleanly to single-threaded C11/C++17 compilation
# (e.g., Apple Clang without libomp).
compile_c_helper() {
    local src="$1" out="$2"
    local cands=()
    if [ -n "${CC:-}" ]; then
        cands=("$CC")
    else
        cands=("gcc")
        for d in "" "/opt/homebrew/bin/" "/usr/local/bin/"; do
            for v in 16 15 14 13 12 11; do
                if command -v "${d}gcc-${v}" >/dev/null 2>&1; then
                    cands+=("${d}gcc-${v}")
                fi
            done
        done
    fi
    for cc in "${cands[@]}"; do
        if "$cc" -O3 -Wall -march=native -mpopcnt -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cc" -O3 -Wall -march=native -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cc" -O3 -Wall -mcpu=native -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cc" -O2 -Wall -fopenmp "$src" -o "$out" 2>/dev/null; then
            return 0
        fi
    done
    local fb="${CC:-gcc}"
    "$fb" -O3 -Wall -march=native "$src" -o "$out" 2>/dev/null ||
    "$fb" -O3 -Wall -mcpu=native "$src" -o "$out" 2>/dev/null ||
    "$fb" -O2 -Wall "$src" -o "$out"
}

compile_cxx_helper() {
    local src="$1" out="$2"
    local cands=()
    if [ -n "${CXX:-}" ]; then
        cands=("$CXX")
    else
        cands=("g++")
        for d in "" "/opt/homebrew/bin/" "/usr/local/bin/"; do
            for v in 16 15 14 13 12 11; do
                if command -v "${d}g++-${v}" >/dev/null 2>&1; then
                    cands+=("${d}g++-${v}")
                fi
            done
        done
    fi
    for cxx in "${cands[@]}"; do
        if "$cxx" -O3 -std=c++17 -Wall -march=native -mpopcnt -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cxx" -O3 -std=c++17 -Wall -march=native -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cxx" -O3 -std=c++17 -Wall -mcpu=native -fopenmp "$src" -o "$out" 2>/dev/null ||
           "$cxx" -O2 -std=c++17 -Wall -fopenmp "$src" -o "$out" 2>/dev/null; then
            return 0
        fi
    done
    local fb="${CXX:-g++}"
    "$fb" -O3 -std=c++17 -Wall -march=native "$src" -o "$out" 2>/dev/null ||
    "$fb" -O3 -std=c++17 -Wall -mcpu=native "$src" -o "$out" 2>/dev/null ||
    "$fb" -O2 -std=c++17 -Wall "$src" -o "$out"
}

compile_c_helper route_a/minwords_mt.c "$TMP/minwords_mt"
compile_cxx_helper route_a/rankmask.cpp "$TMP/rankmask"

export MINWORDS_MT_BIN="$TMP/minwords_mt"
export RANKMASK_BIN="$TMP/rankmask"
export TII_OUTPUT_DIR="$TMP"
THREADS="${TII_THREADS:-64}"

if command -v sage >/dev/null 2>&1; then
    SAGE_PY="$(sage -c 'import sys; print(sys.executable)')"
    if [ "$MODE" = "--demo" ]; then
        echo "[0] Running live SageMath Route A key recovery on TII-133, 83, 91, 253..."
        "$SAGE_PY" tii_runner.py A 133,83,91,253 "$THREADS"
    else
        echo "[0] Running live SageMath Route A (133,83,91,125,253,254) and Route B (172)..."
        "$SAGE_PY" tii_runner.py A 133,83,91,125,253,254 "$THREADS"
        "$SAGE_PY" tii_runner.py B 172
    fi
else
    echo "sage unavailable; running pure-Python GF(2)/GF(2^m) witness & algebraic verifier on shipped witnesses."
fi

python3 verify_witnesses_pure.py "$MODE" "$TMP"
echo "RESULT: SUCCESS"
