#!/bin/bash
# =====================================================================
# Two-Query Encryption-Only Splice Forgery on COLM_0
#
# Evaluates all four splice variants (fade, transplant, misuse, prefix)
# against:
#   - 2016 CAESAR SUPERCOP implementations: colm0v1/{ref,aesni4,aesni8}
#     and colm127v1/ref (control)
#   - 2024 Journal of Cryptology specification (Fig. 2 & Fig. 3):
#     AES-COLM_0-128 (128-bit key) and AES-COLM_0-256 (256-bit key)
# together with 16 length-floor and associated-data boundary sweeps and
# one-query baseline checks.
#
#   ./run.sh --demo   fast verification (25 trials per variant)
#   ./run.sh --full   full verification (1000 trials per variant)
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:---demo}"
case "$MODE" in
    --full)
        T=1000
        TS=40
        export NRAND=2000
        ;;
    --demo)
        T=25
        TS=3
        export NRAND=200
        ;;
    -h|--help)
        echo "Usage: $0 [--demo|--full]"
        echo "  --demo   Fast verification (25 trials per variant, default)"
        echo "  --full   Full verification (1,000 trials per variant + 2,000 random controls)"
        exit 0
        ;;
    *)
        echo "Usage: $0 [--demo|--full]" >&2
        exit 1
        ;;
esac

../fetch_third_party.sh --ensure

if [ -z "${CC:-}" ]; then
    if command -v gcc >/dev/null 2>&1; then
        CC=gcc
    else
        CC=cc
    fi
fi
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT
FAIL=0

echo "=== 1. Building CAESAR SUPERCOP and 2024 J. Cryptol. COLM binaries ==="
"$CC" -O2 -w -I. -Icolm0v1_ref -o "$SCRATCH/splice0" splice.c colm0v1_ref/encrypt.c -lm
"$CC" -O2 -w -I. -Icolm127v1_ref -DTAG_STRIDE=127 -DEXPECT_REJECT -o "$SCRATCH/splice127" splice.c colm127v1_ref/encrypt.c -lm
"$CC" -O2 -w -I. -DCOLM_KEY_BITS=128 -o "$SCRATCH/joc128" splice.c joc2024_colm0.c -lm
"$CC" -O2 -w -I. -DCOLM_KEY_BITS=256 -o "$SCRATCH/joc256" splice.c joc2024_colm0.c -lm
"$CC" -O2 -w -I. -o "$SCRATCH/aes_kat" aes_kat.c
"$CC" -O2 -w -I. -Icolm0v1_ref -o "$SCRATCH/onequery" onequery.c colm0v1_ref/encrypt.c -lm
echo "  built: splice0 (colm0v1/ref), splice127 (colm127v1/ref), joc128 (AES-COLM_0-128), joc256 (AES-COLM_0-256), aes_kat, onequery"
echo

echo "=== 2. FIPS-197 Known-Answer Test for AES-128 ==="
"$SCRATCH/aes_kat"
echo

echo "=== 3. Four Two-Query Encryption-Only Splice Variants on COLM_0 (colm0v1/ref) ==="
echo "-------------------------------------------------------------------"
"$SCRATCH/splice0" fade       20261103 "$T" 130 || FAIL=1
echo "-------------------------------------------------------------------"
"$SCRATCH/splice0" transplant 20261103 "$T" 129 || FAIL=1
echo "-------------------------------------------------------------------"
"$SCRATCH/splice0" misuse     20261103 "$T" 131 1 || FAIL=1
echo "-------------------------------------------------------------------"
"$SCRATCH/splice0" prefix     20261103 "$T" 300 100 16 80 129 || FAIL=1
echo

echo "=== 3b. 2024 Journal of Cryptology Specification: AES-COLM_0-128 & AES-COLM_0-256 ==="
echo "  -- AES-COLM_0-128 (J. Cryptol. 2024, Fig. 2 & Fig. 3, 128-bit key) --"
"$SCRATCH/joc128" fade       20261104 "$T" 130 || FAIL=1
"$SCRATCH/joc128" transplant 20261104 "$T" 129 || FAIL=1
echo "  -- AES-COLM_0-256 (J. Cryptol. 2024, Fig. 2 & Fig. 3, 256-bit key) --"
"$SCRATCH/joc256" fade       20261104 "$T" 130 || FAIL=1
"$SCRATCH/joc256" transplant 20261104 "$T" 129 || FAIL=1
"$SCRATCH/joc256" misuse     20261104 "$T" 131 1 || FAIL=1
"$SCRATCH/joc256" prefix     20261104 "$T" 300 100 16 80 129 || FAIL=1
echo

echo "=== 4. Length-Floor and Associated-Data Boundary Sweeps ==="
check_sweep() {
    local label="$1" expect="$2"
    shift 2
    local out rc=0
    out=$(NRAND=10 "$@" 2>&1) || rc=$?
    local got="REJECTED"
    if echo "$out" | grep -q "forgeries ACCEPTED = [1-9]"; then
        got="ACCEPTED"
    fi
    printf "  %-58s %-10s (expected %s)\n" "$label" "$got" "$expect"
    if [ "$got" != "$expect" ]; then
        FAIL=1
    fi
}

check_sweep "fade at l = 130 (exact minimum)"              ACCEPTED "$SCRATCH/splice0" fade 11 "$TS" 130
check_sweep "fade at l = 129 (one block below floor)"      REJECTED "$SCRATCH/splice0" fade 11 "$TS" 129
check_sweep "fade at l = 128"                              REJECTED "$SCRATCH/splice0" fade 11 "$TS" 128
check_sweep "fade at l = 131"                              ACCEPTED "$SCRATCH/splice0" fade 11 "$TS" 131
check_sweep "fade at l = 258"                              ACCEPTED "$SCRATCH/splice0" fade 11 "$TS" 258
check_sweep "transplant at l = 129 (exact minimum)"        ACCEPTED "$SCRATCH/splice0" transplant 11 "$TS" 129
check_sweep "transplant at l = 128 (Condition C fails)"    REJECTED "$SCRATCH/splice0" transplant 11 "$TS" 128
check_sweep "transplant at l = 130"                        ACCEPTED "$SCRATCH/splice0" transplant 11 "$TS" 130
check_sweep "transplant at l = 200"                        ACCEPTED "$SCRATCH/splice0" transplant 11 "$TS" 200
check_sweep "misuse, s = 1, l = 131 (exact minimum)"       ACCEPTED "$SCRATCH/splice0" misuse 11 "$TS" 131 1
check_sweep "misuse, s = 1, l = 130 (one block below)"     REJECTED "$SCRATCH/splice0" misuse 11 "$TS" 130 1
check_sweep "prefix, shared prefix 1..128 blocks"          REJECTED "$SCRATCH/splice0" prefix 11 "$TS" 300 100 16 80 128
check_sweep "prefix, shared prefix 1..129 blocks"          ACCEPTED "$SCRATCH/splice0" prefix 11 "$TS" 300 100 16 80 129
check_sweep "AD lengths |A1|=0, |A2|=0 (empty AD)"         ACCEPTED env AD1=0 AD2=0 "$SCRATCH/splice0" fade 11 "$TS" 130
check_sweep "AD lengths |A1|=1, |A2|=47 (unaligned AD)"    ACCEPTED env AD1=1 AD2=47 "$SCRATCH/splice0" fade 11 "$TS" 130
check_sweep "query 1 complete, query 2 15-byte last block" ACCEPTED env PARTIAL=16 PARTIAL2=15 "$SCRATCH/splice0" fade 11 "$TS" 130
echo

echo "=== 5. Control: COLM_127 (tau = 127 Intermediate Tags) Rejects All Splices ==="
"$SCRATCH/splice127" fade       20261103 "$T" 130 || FAIL=1
"$SCRATCH/splice127" transplant 20261103 "$T" 129 || FAIL=1
"$SCRATCH/splice127" misuse     20261103 "$T" 131 1 || FAIL=1
echo

echo "=== 6. Cross-Implementation Verification (SUPERCOP ref/aesni4/aesni8 & J. Cryptol. 2024) ==="
"$CC" -O2 -w -I. -Icolm0v1_ref \
    -Dcrypto_aead_encrypt=ref_encrypt -Dcrypto_aead_decrypt=ref_decrypt \
    -c colm0v1_ref/encrypt.c -o "$SCRATCH/ref.o"
"$CC" -O2 -w -I. -DCOLM_KEY_BITS=128 \
    -Dcrypto_aead_encrypt=ni_encrypt -Dcrypto_aead_decrypt=ni_decrypt \
    -c joc2024_colm0.c -o "$SCRATCH/joc128.o"
"$CC" -O2 -w -o "$SCRATCH/xcheck_joc128" xcheck.c "$SCRATCH/ref.o" "$SCRATCH/joc128.o"
echo "  -- cross-checking colm0v1/ref against J. Cryptol. 2024 Fig. 2/3 (AES-COLM_0-128) --"
"$SCRATCH/xcheck_joc128" "$T" 20261104 || FAIL=1

ARCH="$(uname -m)"
HAS_AESNI_AVX2=0
if [ "$ARCH" = "x86_64" ] || [ "$ARCH" = "amd64" ]; then
    if grep -qE '\baes\b' /proc/cpuinfo 2>/dev/null && grep -qE '\bavx2\b' /proc/cpuinfo 2>/dev/null; then
        HAS_AESNI_AVX2=1
    elif [ "$(uname -s)" = "Darwin" ] && \
         [ "$(sysctl -n hw.optional.avx2_0 2>/dev/null || echo 0)" = "1" ] && \
         sysctl -n machdep.cpu.features 2>/dev/null | grep -qiw AES; then
        HAS_AESNI_AVX2=1
    fi
fi

if [ "$HAS_AESNI_AVX2" -eq 1 ]; then
    for ni in aesni4 aesni8; do
        src_c=$(ls colm0v1_${ni}/encrypt*.c)
        "$CC" -O2 -w -maes -msse4.1 -mavx2 -I. -Icolm0v1_${ni} \
            -Dcrypto_aead_encrypt=ni_encrypt -Dcrypto_aead_decrypt=ni_decrypt \
            -c "$src_c" -o "$SCRATCH/${ni}.o"
        "$CC" -O2 -w -o "$SCRATCH/xcheck_${ni}" xcheck.c "$SCRATCH/ref.o" "$SCRATCH/${ni}.o"
        echo "  -- cross-checking colm0v1/ref against colm0v1/${ni} --"
        "$SCRATCH/xcheck_${ni}" "$T" 20261103 || FAIL=1
    done
else
    echo "  [skip] colm0v1/aesni4 and colm0v1/aesni8 require x86_64 with AES-NI + AVX2 (detected host: ${ARCH})"
fi
echo

echo "=== 7. One-Query Baseline Strategies (All 8 Strategies Rejected) ==="
"$SCRATCH/onequery" "$T" 130 || FAIL=1
echo

if [ "$FAIL" -eq 0 ]; then
    echo "RESULT: SUCCESS"
    exit 0
else
    echo "RESULT: FAIL"
    exit 1
fi
