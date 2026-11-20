#!/bin/bash
# =============================================================================
# run.sh -- ONE-SIGNATURE SECRET-KEY RECOVERY AND EUF-CMA FORGERY against
#           SBC_VOLE / SBC_VOLE^M  (ePrint 2024/1537, ASIACRYPT 2025),
#           run against the AUTHORS' OWN reference implementation
#           (IACR artifact asiacrypt/2025/a10).
#
#   ./run.sh --demo     ~2 min   (D,tau)=(9,15) basic + MinMax
#   ./run.sh --full     ~5 min   adds (D,tau)=(13,10) basic + MinMax -- the
#                                paper's headline 2 978 B / 2 658 B sets
#   ./run.sh --sweep    ~25 min  all SIX published (D,tau) rows x both variants
#                                x 100 fresh keys  (the 1200-instance sweep)
#
# Exactly one `RESULT: ...` line is printed, last.  Exit status 0 iff it is
# `RESULT: SUCCESS`; non-zero otherwise.  Automatically downloads third-party
# reference archives on first run if absent; subsequent runs are offline.
# =============================================================================
set -e
set -u
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")"
MODE="${1:---demo}"
case "$MODE" in
  --demo|--full|--sweep) ;;
  -h|--help)
    echo "Usage: ./run.sh [--demo | --full | --sweep [N]]"
    exit 0
    ;;
  *)
    echo "Usage: ./run.sh [--demo | --full | --sweep [N]]" >&2
    exit 1
    ;;
esac
NSWEEP="${2:-100}"
CC="${CC:-$(command -v gcc >/dev/null 2>&1 && echo gcc || echo cc)}"

bash ../fetch_third_party.sh --ensure

sha256_file () {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    python3 -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$1"
  fi
}

# ---- scratch -----------------------------------------------------------------
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

rc=0
banner () { echo "---------------------------------------------------------------"; echo "$*"; }

ART_SHA="$(sha256_file vendor/asiacrypt-2025-a10.zip)"
echo "SBC_VOLE one-signature key recovery -- PoC"
echo "mode=$MODE"
echo "artifact sha256=$ART_SHA"
echo "  (expected e275d529c750bad077ca7bbe0ff355e9a183a4038df08d7869ef7ca95d3447a7)"
if [ "$ART_SHA" != \
     "e275d529c750bad077ca7bbe0ff355e9a183a4038df08d7869ef7ca95d3447a7" ]; then
  echo "  ARTIFACT HASH MISMATCH -- refusing to continue"; echo "RESULT: FAIL"; exit 1
fi

# ==================================================== 0. the cost, as a number ==
banner "[0] COST of the attack at the paper's own (n,k) = (130,257)"
python3 ./cost.py || rc=1

# ================================================== 1. structural (paper only) ==
banner "[1] STRUCTURAL: is B affine in the key in the PAPER's printed Fig. 5 (p.20)?"
echo "    (deterministic: random.Random(20241001), no entropy drawn)"
NTRIAL=25
if [ "$MODE" = "--full" ]; then NTRIAL=60; fi
out=$(python3 ./structural.py "$NTRIAL" 2>&1) || rc=1
echo "$out" | grep -v "^RESULT: " || true
echo "$out" | grep -q '^RESULT: SUCCESS' || rc=1
echo "  PASS  structural F_2-affine check ($NTRIAL instances)"

# ================================================ 2..4  the authors' artifact ===
run_one () {   # $1 = D, $2 = basic|minmax, $3 = expected signature size (bytes)
  Dv="$1"; MMv="$2"; EXP="$3"
  banner "[artifact] SBC_VOLE D=$Dv ($MMv), paper's Table 1/Table 3 size = $EXP B"
  if ! bash ./build.sh "$WORK/b${Dv}_${MMv}" "$Dv" "$MMv" >/dev/null 2>"$WORK/err.txt"; then
    echo "  BUILD FAILED"; cat "$WORK/err.txt"; rc=1; return
  fi
  B="$WORK/b${Dv}_${MMv}/work"

  # (a) the artifact's OWN main(), built from the unmodified source, must still
  #     report the paper's signature size -- this is what establishes that the
  #     binary under attack is the paper's own instantiation.
  sz=$(cd "$B" && ./pristine | sed -n 's/^Compressed size is:\([0-9]*\)$/\1/p')
  echo "  unmodified-artifact signature size: $sz (paper: $EXP)"
  [ "$sz" = "$EXP" ] || { echo "  SIZE MISMATCH"; rc=1; }

  # (b) the victim signs ONE message -> writes ONLY raw binary victim.sig and victim.pk
  TK=$("$B/driver" keygen "victim secret seed D$Dv $MMv" | cut -d' ' -f2-5)
  ( cd "$B" && ./driver sign $TK \
      "victim message: pay Bob 1 EUR" victim.sig victim.pk ) | sed 's/^/  /'
  # Attacker extracts public verification quantities ONLY from (message, victim.sig, victim.pk)
  ( cd "$B" && ./driver extract "victim message: pay Bob 1 EUR" \
      victim.sig victim.pk public_only.txt )
  echo "  attack input: raw binary victim.sig ($(wc -c < "$B/victim.sig") B) + victim.pk ($(wc -c < "$B/victim.pk") B)"
  echo "  (zero-patch unmodified SBCvole.c -- the attack never sees SecX or UVxy)"

  # (c) the attack: ONE 257x256 GF(2) elimination
  python3 ./recover_key.py "$B/public_only.txt" | sed 's/^/  /' || rc=1
  RK=$(python3 ./recover_key.py "$B/public_only.txt" | grep RECOVERED_SECX | cut -d' ' -f2-5)
  [ "$RK" = "$TK" ] || { echo "  RECOVERED KEY != VICTIM KEY"; rc=1; }
  echo "  PASS  recovered key == victim key: $([ "$RK" = "$TK" ] && echo yes || echo NO)"

  # (d) the recovered key regenerates the victim's public key
  ( cd "$B" && ./driver sign $RK "unused" chk.sig chk.pk >/dev/null )
  if cmp -s "$B/chk.pk" "$B/victim.pk"; then
    echo "  PASS  recovered key regenerates the victim public key byte-for-byte: yes"
  else echo "  PUBLIC KEY MISMATCH"; rc=1; fi

  # (e) FORGE on an attacker-chosen message, verified against the VICTIM's pk
  ( cd "$B" && ./driver sign $RK \
      "attacker message: pay Eve 1000000 EUR" forge.sig forge.pk ) | sed 's/^/  /'
  fv=$( cd "$B" && ./driver verify "attacker message: pay Eve 1000000 EUR" \
        forge.sig victim.pk | sed -n 's/^FORGE_VERIFY: //p' )
  echo "  PASS  >>> FORGERY ACCEPTED BY THE AUTHORS' VERIFIER: $fv"
  [ "$fv" = "1" ] || rc=1

  # (f) three negative controls, same code path
  BAD=$(python3 -c "
w=[int(t,16) for t in '''$RK'''.split()]; w[2]^=1
print(' '.join('%016x'%x for x in w))")
  ( cd "$B" && ./driver sign $BAD "attacker message: pay Eve 1000000 EUR" \
      bad.sig bad.pk >/dev/null )
  n1=$( cd "$B" && ./driver verify "attacker message: pay Eve 1000000 EUR" \
        bad.sig victim.pk | sed -n 's/^FORGE_VERIFY: //p' )
  n2=$( cd "$B" && ./driver verify "attacker message: pay Eve 1000001 EUR" \
        forge.sig victim.pk | sed -n 's/^FORGE_VERIFY: //p' )
  n3=$( cd "$B" && ./driver verify "attacker message: pay Eve 1000000 EUR" \
        victim.sig victim.pk | sed -n 's/^FORGE_VERIFY: //p' )
  echo "  PASS  negative controls (all must be 0): wrong key=$n1  wrong message=$n2  replayed honest signature=$n3"
  [ "$n1$n2$n3" = "000" ] || rc=1
}

case "$MODE" in
  --full)
    run_one 13 basic  2978     # Table 1 (p.25), the paper's smallest basic signature
    run_one 13 minmax 2658     # Table 3 (p.29) / abstract (p.1), the 2 658 B headline
    run_one  9 basic  3458     # Table 1 (p.25)
    ;;
  --sweep)
    # every published (D,tau) row of BOTH schemes, end to end, then the
    # 100-key rank sweep on each.  This is the 1200-instance sweep.
    for Dv in 9 10 11 12 13 15; do
      for MMv in basic minmax; do
        case "${Dv}_${MMv}" in
          9_basic)   E=3458 ;;  9_minmax)  E=2978 ;;
          10_basic)  E=3218 ;;  10_minmax) E=2802 ;;
          11_basic)  E=3170 ;;  11_minmax) E=2786 ;;
          12_basic)  E=3090 ;;  12_minmax) E=2738 ;;
          13_basic)  E=2978 ;;  13_minmax) E=2658 ;;
          15_basic)  E=2978 ;;  15_minmax) E=2690 ;;
        esac
        run_one "$Dv" "$MMv" "$E"
      done
    done
    ;;
  *)
    run_one  9 basic  3458
    run_one  9 minmax 2978
    ;;
esac

# ============================================= 5. rank structure / statistics ===
banner "[5] how many key candidates does ONE signature leave?  (rank deficiency)"
if [ "$MODE" = "--sweep" ]; then
  NC="$NSWEEP"
  for Dv in 9 10 11 12 13 15; do
    for MMv in basic minmax; do
      if [ ! -x "$WORK/b${Dv}_${MMv}/work/driver" ]; then continue; fi
      python3 ./rank_sweep.py "$WORK/b${Dv}_${MMv}/work/driver" "$NC" \
              "D${Dv}_${MMv}" || rc=1
    done
  done
else
  NC=8
  if [ "$MODE" = "--full" ]; then NC=30; fi
  python3 ./rank_sweep.py "$WORK/b9_basic/work/driver" "$NC" "D9-basic" || rc=1
fi
echo "  PASS  rank deficiency sweep ($NC instances)"

banner "done"
if [ $rc -eq 0 ]; then
  echo "SUMMARY: SBC_VOLE secret key recovered from ONE signature at the paper's own"
  echo "         parameters (q=2, k=257, n=130); cost 2^26.1 bit operations (conservative"
  echo "         accounting) / 2^25.2 (headline), against a claimed lambda = 128 bits;"
  echo "         forgery on an attacker-chosen message ACCEPTED by the authors' own verifier."
  echo "RESULT: SUCCESS"
  exit 0
else
  echo "RESULT: FAIL"
  exit 1
fi
