#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

export PYTHONDONTWRITEBYTECODE=1

MODE="${1:---demo}"
SUITE="${2:-}"

case "$MODE" in
    --demo|--full) ;;
    --suite)
        case "$SUITE" in
            ciphertext|sparse|pubkey|window|falsifier|ledger) ;;
            *)
                echo "Unknown suite '$SUITE'. Available: ciphertext | sparse | pubkey | window | falsifier | ledger" >&2
                exit 1
                ;;
        esac
        ;;
    -h|--help)
        echo "Usage: ./run.sh [--demo | --full | --suite <ciphertext|sparse|pubkey|window|falsifier|ledger>]"
        exit 0
        ;;
    *)
        echo "Usage: ./run.sh [--demo | --full | --suite <ciphertext|sparse|pubkey|window|falsifier|ledger>]" >&2
        exit 1
        ;;
esac

./fetch_third_party.sh --ensure

HOST_ARCH="$(uname -m)"
CXX="${CXX:-g++}"
CC="${CC:-gcc}"
export CXX CC

cxx_accepts() {
    printf 'int main() { return 0; }\n' | "$CXX" -x c++ - -o /dev/null "$@" >/dev/null 2>&1
}

OMP_FLAGS=()
if cxx_accepts -std=c++17 -fopenmp; then
    OMP_FLAGS=(-fopenmp)
else
    for omp_pfx in /opt/homebrew/opt/libomp /usr/local/opt/libomp; do
        if [ -d "$omp_pfx" ] && cxx_accepts -std=c++17 -Xpreprocessor -fopenmp "-I$omp_pfx/include" "-L$omp_pfx/lib" -lomp; then
            OMP_FLAGS=(-Xpreprocessor -fopenmp "-I$omp_pfx/include" "-L$omp_pfx/lib" -lomp)
            break
        fi
    done
fi

ARCH_FLAGS=()
if cxx_accepts -std=c++17 -march=native; then
    ARCH_FLAGS+=(-march=native)
fi
if cxx_accepts -std=c++17 -mpclmul; then
    ARCH_FLAGS+=(-mpclmul)
fi

run_ciphertext_schur() {
    echo ">>> [Ciphertext / Schur LOCKER-I] verify_5trials_standalone.py"
    case "$HOST_ARCH" in
        x86_64|amd64) ;;
        *) echo "    (Note: $HOST_ARCH lacks x86_64 PCLMULQDQ; kemsolve.cpp uses scalar carry-less fallback, ~60-130 s/trial)" ;;
    esac
    python3 ciphertext_attacks/schur_locker1/verify_5trials_standalone.py
}

run_cyclotomic_check() {
    echo ">>> [Ciphertext / Cyclotomic Z^n-1] make_phi_variant.py & kemsolve_phi83.cpp build"
    local tmp_dir
    tmp_dir="$(mktemp -d)"
    python3 ciphertext_attacks/cyclotomic_ring/make_phi_variant.py 83 "$tmp_dir/gen_phi83.cpp"
    "$CXX" -O3 -std=c++17 ${ARCH_FLAGS[@]+"${ARCH_FLAGS[@]}"} ${OMP_FLAGS[@]+"${OMP_FLAGS[@]}"} -o "$tmp_dir/kemsolve_phi83" ciphertext_attacks/cyclotomic_ring/kemsolve_phi83.cpp
    rm -rf "$tmp_dir"
    echo "[PASS] Cyclotomic Phi_83 solver derived and compiled cleanly."
}

run_sparse_krylov() {
    local mode="${1:-quick}"
    echo ">>> [Ciphertext / Sparse Krylov] PCLMULQDQ GF(2^n) & CSR Wiedemann engine ($mode)"
    case "$HOST_ARCH" in
        x86_64|amd64)
            if [ "$mode" = "full" ]; then
                bash ciphertext_attacks/sparse_krylov/verify_s16.sh
            else
                local tmp_dir
                tmp_dir="$(mktemp -d)"
                "$CXX" -O3 ${ARCH_FLAGS[@]+"${ARCH_FLAGS[@]}"} -std=c++17 -funroll-loops -o "$tmp_dir/test_gf2n" ciphertext_attacks/sparse_krylov/test_gf2n.cpp
                "$CXX" -O3 ${ARCH_FLAGS[@]+"${ARCH_FLAGS[@]}"} -std=c++17 ${OMP_FLAGS[@]+"${OMP_FLAGS[@]}"} -funroll-loops -Wno-unused-result -o "$tmp_dir/tsm_bw" ciphertext_attacks/sparse_krylov/tsm_bw.cpp
                "$tmp_dir/test_gf2n"
                python3 ciphertext_attacks/sparse_krylov/gen_instance.py --m 20 --n 23 --r 3 --d 4 --P 5 --gmod 20,3,0 --Pmod 23,5,0 --out "$tmp_dir/isq.bin" --seed 501 >/dev/null
                "$tmp_dir/tsm_bw" build "$tmp_dir/isq.bin" "$tmp_dir/msq.mat" 20,3,0 23,5,0 2>/dev/null
                "$tmp_dir/tsm_bw" solve "$tmp_dir/msq.mat" 23,5,0 23 2>/dev/null
                python3 ciphertext_attacks/sparse_krylov/cost_model.py >/dev/null
                rm -rf "$tmp_dir"
                echo "[PASS] Sparse PCLMULQDQ kernel, CSR Wiedemann solver, and cost model verified."
            fi
            ;;
        *)
            echo "[SKIP] Native sparse C++ kernels (gf2n_pclmul.hpp, test_gf2n.cpp, tsm_bw.cpp) require x86_64 PCLMULQDQ/SSE2 (detected: $HOST_ARCH)."
            local tmp_dir
            tmp_dir="$(mktemp -d)"
            python3 ciphertext_attacks/sparse_krylov/gen_instance.py --m 20 --n 23 --r 3 --d 4 --P 5 --gmod 20,3,0 --Pmod 23,5,0 --out "$tmp_dir/isq.bin" --seed 501 >/dev/null
            python3 ciphertext_attacks/sparse_krylov/cost_model.py > "$tmp_dir/cost.json"
            if [ "$mode" = "full" ]; then
                COST_JSON="$tmp_dir/cost.json" python3 - <<'PYEOF'
import json, os
for r in json.load(open(os.environ["COST_JSON"])):
    if r["P_min"] is None:
        continue
    print("%-22s 2^%-4d claimed -> 2^%.1f bit-ops (Lanczos)  margin %3.0f bits  sparse %.2f GB  dense %.0f GB"
          % (r["name"], r["claimed_log2"], r["log2_bitops_block_lanczos"],
             r["claimed_log2"] - r["log2_bitops_block_lanczos"],
             r["sparse_bytes"] / 1e9, r["dense_bytes"] / 1e9))
PYEOF
            fi
            rm -rf "$tmp_dir"
            echo "[PASS] Sparse Krylov Python instance generator and exact gate cost model verified."
            ;;
    esac
}

run_pubkey_direct() {
    echo ">>> [Public Key / Direct Closed-Form] closed_form_rollo_pk_break.py"
    python3 pubkey_attacks/direct_closed_form/closed_form_rollo_pk_break.py
}

run_pubkey_window() {
    local flag="${1---quick}"
    echo ">>> [Public Key / Subspace Window] verify_window_attack.py ${flag}"
    if [ -n "$flag" ]; then
        python3 pubkey_attacks/subspace_window/verify_window_attack.py "$flag"
    else
        python3 pubkey_attacks/subspace_window/verify_window_attack.py
    fi
}

run_falsifier() {
    local mode="${1:-quick}"
    echo ">>> [Barriers / Falsifier] Partial-Frobenius syzygy & Plucker deficit ($mode)"
    if [ "$mode" = "full" ]; then
        local tmp_dir=""
        if pkg-config --exists m4ri 2>/dev/null; then
            tmp_dir="$(mktemp -d)"
            "$CC" -O3 -o "$tmp_dir/gf2rank" benchmarks_and_barriers/falsifier/gf2rank.c $(pkg-config --cflags --libs m4ri)
            GF2RANK="$tmp_dir/gf2rank" python3 benchmarks_and_barriers/falsifier/s19_falsifier.py
            rm -rf "$tmp_dir"
        elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/include/m4ri/m4ri.h" ]; then
            tmp_dir="$(mktemp -d)"
            "$CC" -O3 -o "$tmp_dir/gf2rank" benchmarks_and_barriers/falsifier/gf2rank.c "-I$CONDA_PREFIX/include" "-L$CONDA_PREFIX/lib" "-Wl,-rpath,$CONDA_PREFIX/lib" -lm4ri -lm
            GF2RANK="$tmp_dir/gf2rank" python3 benchmarks_and_barriers/falsifier/s19_falsifier.py
            rm -rf "$tmp_dir"
        else
            python3 benchmarks_and_barriers/falsifier/s19_falsifier.py
        fi
    else
        python3 benchmarks_and_barriers/falsifier/exp5_sigma_equivariance.py
        python3 benchmarks_and_barriers/falsifier/exp10_degree2_obstruction.py
    fi
}

run_ledger() {
    local mode="${1:-quick}"
    echo ">>> [Benchmarks / 27-Set Complexity Table] ($mode)"
    if [ "$mode" = "full" ]; then
        python3 benchmarks_and_barriers/estimators/frobenius_applicability.py
    fi
    python3 benchmarks_and_barriers/estimators/rollo_structural_attack.py
    python3 benchmarks_and_barriers/estimators/compare_all27_master.py
}

echo "=============================================================================="
echo "LOCKER / ROLLO Ideal-LRPC Cryptanalysis Artifact Suite ($MODE ${SUITE})"
echo "=============================================================================="

case "$MODE" in
    --demo)
        run_ciphertext_schur
        echo ""
        run_pubkey_direct
        echo ""
        python3 benchmarks_and_barriers/falsifier/exp5_sigma_equivariance.py
        echo ""
        python3 benchmarks_and_barriers/estimators/compare_all27_master.py
        ;;
    --full)
        run_ciphertext_schur
        echo ""
        run_cyclotomic_check
        echo ""
        run_sparse_krylov quick
        echo ""
        run_pubkey_direct
        echo ""
        run_pubkey_window --quick
        echo ""
        run_falsifier quick
        echo ""
        run_ledger quick
        ;;
    --suite)
        case "$SUITE" in
            ciphertext) run_ciphertext_schur && run_cyclotomic_check ;;
            sparse)     run_sparse_krylov full ;;
            pubkey)     run_pubkey_direct ;;
            window)     run_pubkey_window "" ;;
            falsifier)  run_falsifier full ;;
            ledger)     run_ledger full ;;
            *)
                echo "Unknown suite '$SUITE'. Available: ciphertext | sparse | pubkey | window | falsifier | ledger"
                exit 1
                ;;
        esac
        ;;
    *)
        echo "Usage: ./run.sh [--demo | --full | --suite <ciphertext|sparse|pubkey|window|falsifier|ledger>]"
        exit 1
        ;;
esac

echo "=============================================================================="
echo "RESULT: SUCCESS"
echo "=============================================================================="
