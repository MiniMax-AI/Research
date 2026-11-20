#!/bin/bash
# =============================================================================
# run.sh -- Top-level verification runner for all SBC_VOLE / SBC_VOLE^M
#           cryptanalysis and protocol-repair artifacts.
#
#   ./run.sh --demo     Fast end-to-end verification across all 3 suites (~40 s)
#   ./run.sh --full     Extended verification on headline (D, tau) = (13, 10)
#
# First run automatically fetches third-party reference archives if absent;
# subsequent runs execute completely offline using only a C compiler and Python 3.8+.
# =============================================================================
set -e
set -u
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"

MODE="${1:---demo}"
case "$MODE" in
  --demo|--full) ;;
  -h|--help)
    echo "Usage: $0 [--demo | --full]"
    exit 0
    ;;
  *)
    echo "Usage: $0 [--demo | --full]" >&2
    exit 2
    ;;
esac

bash ./fetch_third_party.sh --ensure

rc=0
run_suite () {
  local dir="$1"
  local title="$2"
  echo "======================================================================="
  echo ">>> [$dir] $title ($MODE)"
  echo "======================================================================="
  if bash "./$dir/run.sh" "$MODE"; then
    echo ">>> [$dir] PASS"
  else
    echo ">>> [$dir] FAIL"
    rc=1
  fi
  echo
}

run_suite "key_recovery"       "Single-Signature Secret-Key Recovery & EUF-CMA Forgery"
run_suite "hvzk_distinguisher" "HVZK Zero-Conditional-Entropy & Game 2 Distinguisher"
run_suite "protocol_repair"    "Correlated VOLE-Masking Protocol Repair"

echo "======================================================================="
if [ "$rc" -eq 0 ]; then
  echo "ALL 3 ARTIFACT SUITES PASSED ($MODE)"
  echo "RESULT: SUCCESS"
  exit 0
else
  echo "ONE OR MORE ARTIFACT SUITES FAILED ($MODE)"
  echo "RESULT: FAIL"
  exit 1
fi
