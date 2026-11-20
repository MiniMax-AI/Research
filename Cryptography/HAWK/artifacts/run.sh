#!/bin/bash
set -e
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1

MODE="${1:---demo}"
case "$MODE" in
-h|--help)
    echo "usage: $0 [--demo|--kat|--full]"
    echo "  --demo, --kat  Fast verification of algebraic certificates, security estimates, + 592B secret key (~4 s, default)"
    echo "  --full         Rebuild Lambda_tau & Lambda' from pk.bin, run g6k d=80 tour (if installed),"
    echo "                 verify certificates & security estimates, run tau-tower descent, and verify secret key (~31 min)"
    exit 0
    ;;
--demo|--kat|--full)
    ;;
*)
    echo "usage: $0 [--demo|--kat|--full]" >&2
    exit 2
    ;;
esac

./fetch_third_party.sh --ensure

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export HAWK_WORK="$TMP"

if [ -n "${PY:-}" ]; then
    :
elif command -v sage >/dev/null 2>&1; then
    PY="$(sage -c 'import sys; print(sys.executable)')"
else
    PY="python3"
fi

if ! "$PY" -c "import sage.all" >/dev/null 2>&1; then
    echo "ERROR: SageMath (sage.all) not found via '$PY'." >&2
    echo "Install SageMath 10.x and ensure 'sage' is on PATH (or set PY=/path/to/sage-python)." >&2
    exit 1
fi

SHIM_BIN="$TMP/hawk_shim"
verify_with_shim() {
    local pk_file="$1"
    local sk_file="$2"
    local label="$3"
    local rc sig_out sig_hex

    if [ ! -x "$SHIM_BIN" ]; then
        bash ./kat/build_shim.sh "$SHIM_BIN"
    fi

    rc=0
    sig_out=$("$SHIM_BIN" signverify "$pk_file" "$sk_file" 2>&1) || rc=$?

    if [ "$rc" -eq 0 ]; then
        echo "$sig_out" | sed 's/^/  /'
        echo "$sig_out" | grep -Fq "SIGN(recovered_sk)/VERIFY(original_pk) = PASS" || {
            echo "FAIL: reference signverify did not report PASS" >&2
            exit 1
        }
        sig_hex=$(echo "$sig_out" | awk '/^sig\[:32\]:/ {print $2}')
        echo "  PASS  $label: sig[:32]=$sig_hex"
        return 0
    fi

    echo "$sig_out" | sed 's/^/  /' >&2
    echo "FAIL: reference signverify (exit code $rc)" >&2
    exit 1
}

echo "======================================================================"
echo "HAWK-256 Single-Core Key Recovery Verification Suite"
echo "Mode: $MODE"
echo "======================================================================"

case "$MODE" in
--demo|--kat)
    echo "[1/3] Verifying sublattice certificates, Cholesky pullback, and V_tau..."
    "$PY" -u chain/certificates.py

    echo ""
    echo "[2/3] Verifying GSA block sizes, AGPS20 gate complexities, and scaling bounds..."
    "$PY" -u chain/estimate_security.py

    echo ""
    echo "[3/3] Verifying recovered 592-byte secret key on C reference implementation..."
    verify_with_shim data/pk.bin data/recovered_skd.bin "Reference C signer/verifier acceptance (pk=450B, skd=592B)"
    ;;

--full)
    if ! command -v flatter >/dev/null 2>&1; then
        echo "ERROR: 'flatter' not found on PATH (required by --full for Stages 2 & 4)." >&2
        echo "Install flatter (https://github.com/keeganryan/flatter) or run './run.sh --demo'." >&2
        exit 1
    fi

    echo "[1/5] Rebuilding saturated integer kernel of Lambda_tau and Cholesky-preconditioned Lambda'..."
    "$PY" -u chain/build_lambda.py --full

    echo ""
    echo "[2/5] Running single-thread progressive g6k sieve tour (d = 40 -> 80)..."
    if "$PY" -c "import g6k" 2>/dev/null; then
        "$PY" -u chain/g6k_tour.py "$TMP/Bs.txt" 21 512 1 40,50,60,70,80
    else
        echo "  g6k not installed; using pre-computed reduced basis data/Bred.sobj."
    fi

    echo ""
    echo "[3/5] Verifying algebraic certificates and coordinate pullback v = u * U..."
    "$PY" -u chain/certificates.py

    echo ""
    echo "[4/5] Verifying GSA block sizes, AGPS20 gate complexities, and scaling bounds..."
    "$PY" -u chain/estimate_security.py --simulator

    echo ""
    echo "[5/5] Running tau-tower descent (ell = 9 -> 3), encoding key, and verifying on C reference..."
    "$PY" -u chain/recover.py
    verify_with_shim data/pk.bin "$TMP/recovered_skd.bin" "Full-chain recovered key accepted by C reference"
    ;;
esac

echo "RESULT: SUCCESS"
