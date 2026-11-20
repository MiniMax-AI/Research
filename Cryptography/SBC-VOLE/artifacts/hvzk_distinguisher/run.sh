#!/bin/bash
set -e
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
MODE=demo
if [ "${1:-}" = "--full" ]; then MODE=full; fi
out="$(python3 hvzk_sim.py --mode "$MODE" --out "$TMP/hvzk_sim_results.json")"
echo "$out" | grep -v '^RESULT:'
echo "$out" | grep -q '^RESULT: SUCCESS'
echo "RESULT: SUCCESS"
