#!/bin/bash
# run.sh -- the masked SBC_VOLE repair, implemented on the official ASIACRYPT
# 2025 artifact and costed with the artifact's own serialiser.
#
#   ./run.sh --demo    D=9,  tau=15   (fast; ~15 s)
#   ./run.sh --full    D=13, tau=10   (the paper's headline parameter set)
#
# Automatically fetches third-party reference files on first run if absent;
# needs only a C compiler (gcc or clang) + python3 (stdlib only).
# Prints exactly one RESULT line.
set -e
set -u
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"

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

CC="${CC:-$(command -v gcc >/dev/null 2>&1 && echo gcc || echo cc)}"
if [ "$MODE" = "--full" ]; then D=13; TAU=10; PAPER_BASIC=2978; PAPER_MM=2658
else                          D=9;  TAU=15; PAPER_BASIC=3458; PAPER_MM=2978; fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
SRC="$TMP/build"
mkdir -p "$SRC"
cp -r crypto-algorithms "$TMP/"
python3 make_masked.py "$SRC" || { echo "RESULT: FAIL"; exit 1; }
cp forge_noh0.inc recover_key.py sbcvole.py multi_attack.py pooled_leak.py "$SRC"/

CFL="-O2 -DL$D -DSILENT -DNBREPEATS=1 -DGENERICAES -DGENERICSHA"
AES="../crypto-algorithms/aes.c ../crypto-algorithms/sha256.c"

build() { # $1 = output, rest = extra defines
  local out="$1"; shift
  ( cd "$SRC" && "$CC" $CFL "$@" -I. -o "$out" SBCvole.c SBCtools.c $AES -lm ) 2>&1 \
    | grep -iE "^[^ ]*error" | head -3 || true
  [ -x "$SRC/$out" ] || { echo "BUILD FAILED: $out"; return 1; }
}
run()  { ( cd "$SRC" && ./"$1" 2>&1 ); }
szof() { run "$1" | grep "Compressed size" | sed 's/.*://'; }
corof(){ run "$1" | grep "Return cor"     | sed 's/.*: //'; }

FAIL=0
note() { printf '  PASS  %-52s %s\n' "$1" "$2"; }

echo "== parameter set (D,tau) = ($D,$TAU), artifact asiacrypt-2025-a10"

# ---------------------------------------------------------------- 1. sizes
for MM in "" "-DMINMAXCOMMIT"; do
  TAGv=$([ -n "$MM" ] && echo mm || echo basic)
  build "p_$TAGv" $MM                              || FAIL=1
  build "m_$TAGv" $MM -DMASKB                      || FAIL=1
  P=$(szof "p_$TAGv"); M=$(szof "m_$TAGv"); PC=$(corof "p_$TAGv"); MC=$(corof "m_$TAGv")
  PAPER=$([ -n "$MM" ] && echo "$PAPER_MM" || echo "$PAPER_BASIC")
  note "unmodified artifact, $TAGv: size"  "$P B (paper Table: $PAPER B) verify=$PC"
  note "masked artifact,      $TAGv: size" "$M B   delta = +$((M-P)) B   verify=$MC"
  [ "$P" = "$PAPER" ] || { echo "  !! unmodified size does not match the paper"; FAIL=1; }
  [ "$PC" = "1" ] && [ "$MC" = "1" ] || { echo "  !! a build failed to verify"; FAIL=1; }
  DELTA=$((M-P)); EXPECT=$(( (TAU-1)*32 + ( (2*TAU+7)/8 - (TAU+1+7)/8 ) ))
  note "  predicted delta from the artifact's own bit packing" "+$EXPECT B"
  [ "$DELTA" = "$EXPECT" ] || { echo "  !! delta does not match the packing formula"; FAIL=1; }
done

# ------------------------------------------- 2. Section 4 key-recovery attack on both builds
# 16 independent signatures per build: the artifact's driver reseeds itself
# between iterations (SBCvole.c:1478, randomness[0]++), so each record is a
# distinct instance.
DUMPN=16
for cfg in "dump_plain:" "dump_masked:-DMASKB"; do
  NAME=${cfg%%:*}; DEFS=${cfg#*:}
  ( cd "$SRC" && "$CC" $CFL -DNBREPEATS=$DUMPN -DDUMPATTACK $DEFS -I. -o "$NAME" \
      SBCvole.c SBCtools.c $AES -lm ) 2>&1 | grep -iE "^[^ ]*error" | head -3 || true
  [ -x "$SRC/$NAME" ] || { echo "BUILD FAILED: $NAME"; FAIL=1; }
  rm -f "$SRC/$NAME.txt"
  ( cd "$SRC" && DUMPATTACK="$PWD/$NAME.txt" ./"$NAME" >/dev/null 2>&1 )
done
for NAME in dump_plain dump_masked; do
  echo "-- Section 4's recover_key.py, unmodified, against $NAME:"
  ( cd "$SRC" && python3 multi_attack.py "$NAME.txt" ) | tail -4 | sed 's/^/   /'
done
( cd "$SRC" && python3 multi_attack.py dump_plain.txt )  | grep -q "recovers the key on every record" \
  || { echo "  !! Section 4 key-recovery attack did not reproduce on the unmodified artifact"; FAIL=1; }
( cd "$SRC" && python3 multi_attack.py dump_masked.txt ) | grep -q "fails on every record" \
  || { echo "  !! Section 4 key-recovery attack still produced candidates against the masked build"; FAIL=1; }

# ------------------------- 2b. pooled leakage: 0 bits of key at every m --------
# Deficiency 256 for every m, and 200/200 random keys admitted by the masked
# model against 0/200 by the unmasked one -- on the authors' own data.
cp pooled_leak.py "$SRC"/ 2>/dev/null
echo "-- pooled key leakage from the masked build:"
LRC=0
( cd "$SRC" && python3 pooled_leak.py dump_masked.txt 1 2 4 16 ) > "$SRC/leak.out" 2>&1 || LRC=$?
grep -v '^RESULT:' "$SRC/leak.out" | sed 's/^/   /'
[ "$LRC" = "0" ] || { echo "  !! pooled-leak measurement did not meet its prediction"; FAIL=1; }

# ------------------------------------------------------- 3. pre-challenge h_0 binding check
# dw hashed into h_0 (correct) vs not (mis-bound), with and without the
# cheating prover of forge_noh0.inc, all with a FALSE witness.
build c_ok_forge   -DMASKB -DFORGE_NOH0                        || FAIL=1
build c_bad_nof    -DMASKB -DMASKB_NOH0 -DFALSEWITNESS         || FAIL=1
build c_bad_forge  -DMASKB -DMASKB_NOH0 -DFORGE_NOH0           || FAIL=1
build c_false_only -DFALSEWITNESS                              || FAIL=1
for spec in "c_ok_forge:0" "c_bad_nof:0" "c_bad_forge:1" "c_false_only:0"; do
  b="${spec%%:*}"; want="${spec#*:}"
  C=$(FORGE_NOH0=1 bash -c "cd '$SRC' && ./$b 2>/dev/null" | grep "Return cor" | sed 's/.*: //')
  note "false witness, $b -> verifier says" "$C (want $want)"
  [ "$C" = "$want" ] || { echo "  !! control $b gave $C, wanted $want"; FAIL=1; }
done

echo
if [ "$FAIL" = "0" ]; then
  echo "RESULT: SUCCESS"
  exit 0
else
  echo "RESULT: FAIL"
  exit 1
fi
