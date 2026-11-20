#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "=== [1] GF(2^128) Splice Index-Set Derivation & Length-Floor Proofs (algebra.py) ==="
python3 -B algebra.py | tee "$TMP/algebra.log"
echo "  PASS  algebra.py: GF(2^128) fade, transplant, and COLM_127 obstruction checks"
echo

echo "=== [2] Case A / B / C Algebraic & Exhaustive Subset Verification (derivation.py) ==="
python3 -B derivation.py | tee "$TMP/derivation.log"
echo "  PASS  derivation.py: Case A, Case B, and Case C (m=4..12) verification"
echo

echo "=== [3] Quantitative Evaluation of Theorem 2 & Lemma 3 (thm2.py) ==="
python3 -B thm2.py | tee "$TMP/thm2.log"
echo "  PASS  thm2.py: J. Cryptol. 2024 (2^-105.92) and 2016 Addendum (2^-107.56) bounds"
echo
echo "RESULT: SUCCESS"
