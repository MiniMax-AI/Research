#!/bin/bash
# =====================================================================
# COLM Verification Artifact Suite (Top-Level Runner)
#
# Executes:
#   1. algebra_and_proof_audit/run.sh  -- Pure Python 3 GF(2^128) derivations
#      and security-proof bound audit
#   2. splice_forgery/run.sh           -- End-to-end C splice forgery tests
#      against CAESAR SUPERCOP and J. Cryptol. 2024 implementations
#
# Usage:
#   ./run.sh [--demo|--full]
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:---demo}"
case "$MODE" in
    --demo|--full) ;;
    -h|--help)
        echo "Usage: $0 [--demo|--full]"
        echo "  --demo   Fast verification (~10-15 s, default)"
        echo "  --full   Full verification (~2 min, 1,000 trials/variant + 2,000 controls)"
        exit 0
        ;;
    *)
        echo "Usage: $0 [--demo|--full]" >&2
        exit 1
        ;;
esac

./fetch_third_party.sh --ensure

echo "====================================================================="
echo "  Part I: GF(2^128) Algebraic Derivations & Proof Audit"
echo "====================================================================="
./algebra_and_proof_audit/run.sh
echo

echo "====================================================================="
echo "  Part II: End-to-End Two-Query Splice Forgeries (${MODE})"
echo "====================================================================="
./splice_forgery/run.sh "$MODE"
echo

echo "====================================================================="
echo "  ALL COLM ARTIFACT VERIFICATIONS PASSED (${MODE})"
echo "====================================================================="
