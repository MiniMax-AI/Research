#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify the official NIST PQC
#                         Round 1 LOCKER C/C++ reference implementation
#                         (KEMSET-I..IX) and CLI bridge.
#
# Usage:
#   ./fetch_third_party.sh            Verify or download third-party files
#   ./fetch_third_party.sh --ensure   Silent fast-path if already present;
#                                     informs user and downloads if missing
#   ./fetch_third_party.sh --force    Re-download all third-party files
#   ./fetch_third_party.sh --help     Show usage help and exit
# =============================================================================
set -euo pipefail

ART_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-}"

usage() {
    cat <<'EOF'
Usage: ./fetch_third_party.sh [--ensure | --force | --help]

Options:
  (no args)   Verify third-party reference files, or download if missing
  --ensure    Silent fast-path (exit 0) if already verified; inform user and download if missing
  --force     Re-download and re-verify all third-party reference files from upstream
  -h, --help  Show this help message and exit
EOF
}

if [ "$#" -gt 1 ]; then
    echo "[LOCKER-ROLLO] ERROR: unexpected extra arguments: $*" >&2
    usage >&2
    exit 2
fi

case "$MODE" in
    ""|--ensure|--force) ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        echo "[LOCKER-ROLLO] ERROR: unknown option: $MODE" >&2
        usage >&2
        exit 2
        ;;
esac

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT INT TERM

python3 - "$ART_DIR" "$MODE" "$TMP_DIR" <<'PY'
from __future__ import annotations

import hashlib
import http.client
import io
import json
import os
import re
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zipfile

art_dir = sys.argv[1]
mode = sys.argv[2]
tmp = sys.argv[3]
native_dir = os.path.join(art_dir, "reference", "native")
manifest_file = os.path.join(native_dir, ".sha256_manifest.json")

BENCH_COMMIT = "d243a0a03d7aa367472411114f140d0aace583cd"
BENCH_BASE = f"https://raw.githubusercontent.com/ethz-spylab/cryptanalysis-benchmark/{BENCH_COMMIT}/algorithms/T2_locker/environment/algorithm_src"

ZIP_URLS = [
    "https://csrc.nist.gov/CSRC/media/Projects/Post-Quantum-Cryptography/documents/round-1/submissions/LOCKER.zip",
    f"{BENCH_BASE}/LOCKER.zip",
]
ZIP_SHA256 = "633418593a7b7df9a42deae38421a5f7f82e8b958eb2e853940a46cb1a58c689"
MANIFEST_SHA256 = "a0a917bc0d6196d3fcf128cc8def4f336a5320d7a14cbe11ce534224f006f3f7"

AUX_FILES = {
    "cli.c": "73a2aa8be003cfbf74d5a73d0e1d747e82d8c1e57b188947c670c2437de50afa",
    "cli_stdin.h": "a9720312d4a24219521c852254584b9dca77ab0a9f1d963e056cd2cbdc6e2abc",
    "kem_c_api.cpp": "3cc245e956565c618215499d29027851e28d39c806b2d4ef1d67dc1884249e1a",
}

MAKEFILE_CONTENT = b"""all: cli

cli:
\tsh ./build_cli.sh
\ttest -x ./cli

clean:
\trm -f cli *.o
"""

BUILD_CLI_CONTENT = b"""#!/bin/sh
set -eu
cd "$(dirname "$0")"
CC="${CC:-gcc}"
CXX="${CXX:-g++}"
if [ -f /usr/include/NTL/ZZ.h ]; then
  NTL_CFLAGS=""
  NTL_LDFLAGS=""
elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/include/NTL/ZZ.h" ]; then
  NTL_CFLAGS="-I$CONDA_PREFIX/include"
  NTL_LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
elif [ -f /opt/homebrew/include/NTL/ZZ.h ]; then
  NTL_CFLAGS="-I/opt/homebrew/include"
  NTL_LDFLAGS="-L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib"
elif [ -f /usr/local/include/NTL/ZZ.h ]; then
  NTL_CFLAGS="-I/usr/local/include"
  NTL_LDFLAGS="-L/usr/local/lib -Wl,-rpath,/usr/local/lib"
elif [ -n "${NTL_PREFIX:-}" ] && [ -f "$NTL_PREFIX/include/NTL/ZZ.h" ]; then
  NTL_CFLAGS="-I$NTL_PREFIX/include"
  NTL_LDFLAGS="-L$NTL_PREFIX/lib -Wl,-rpath,$NTL_PREFIX/lib"
else
  echo "NTL development files were not found" >&2
  exit 1
fi
for pfx in /opt/homebrew/opt/openssl /usr/local/opt/openssl /opt/homebrew/opt/gmp /usr/local/opt/gmp; do
  if [ -d "$pfx/include" ] && [ -d "$pfx/lib" ]; then
    NTL_CFLAGS="$NTL_CFLAGS -I$pfx/include"
    NTL_LDFLAGS="$NTL_LDFLAGS -L$pfx/lib"
  fi
done
SRC="Reference_Implementation/KEMSET-I"
INCLUDES="-I. -I$SRC/src -I$SRC/src/ffi -I$SRC/lib/hash -I$SRC/lib/rng $NTL_CFLAGS"
"$CC" -O2 -fcommon -c cli.c $INCLUDES -o cli.o
"$CXX" -O2 -fcommon -c kem_c_api.cpp $INCLUDES -o kem_c_api.o
"$CXX" -O2 -fcommon -c "$SRC"/src/kem.cpp $INCLUDES -o kem.o
"$CXX" -O2 -fcommon -c "$SRC"/src/parsing.cpp $INCLUDES -o parsing.o
"$CXX" -O2 -fcommon -c "$SRC"/src/decoder.cpp $INCLUDES -o decoder.o
"$CXX" -O2 -fcommon -c "$SRC"/src/ffi/ffi_field.cpp $INCLUDES -o ffi_field.o
"$CXX" -O2 -fcommon -c "$SRC"/src/ffi/ffi_elt.cpp $INCLUDES -o ffi_elt.o
"$CXX" -O2 -fcommon -c "$SRC"/src/ffi/ffi_vec.cpp $INCLUDES -o ffi_vec.o
"$CXX" -O2 -fcommon -c "$SRC"/lib/hash/hash.c $INCLUDES -o hash.o
"$CXX" -O2 -fcommon -c "$SRC"/lib/rng/rng.c $INCLUDES -o rng.o
"$CXX" -O2 -o cli \\
    cli.o kem_c_api.o kem.o parsing.o decoder.o \\
    ffi_field.o ffi_elt.o ffi_vec.o hash.o rng.o \\
    $NTL_LDFLAGS -lntl -lgmp -lssl -lcrypto
chmod 0755 cli
"""


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def all_present() -> bool:
    if not os.path.isfile(manifest_file) or sha256_file(manifest_file) != MANIFEST_SHA256:
        return False
    try:
        with open(manifest_file, "r", encoding="utf-8") as f:
            expected = json.load(f)
    except Exception:
        return False
    if len(expected) < 225:
        return False
    for rel, exp in expected.items():
        p = os.path.join(native_dir, rel)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    return True


def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    curl_bin = shutil.which("curl")
    last_err: Exception | None = None
    if curl_bin:
        dl_path = os.path.join(tmp, f"dl_{os.getpid()}_{abs(hash(url)) & 0xFFFFFFFF:08x}.bin")
        cmd = [
            curl_bin, "-fsSL",
            "--retry", "3",
            "--retry-delay", "2",
            "--connect-timeout", "20",
            "--max-time", str(timeout),
            "-A", "Mozilla/5.0 (compatible; artifact-fetcher/1.0)",
            "-o", dl_path,
            url,
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0 and os.path.isfile(dl_path):
            with open(dl_path, "rb") as f:
                data = f.read()
            os.remove(dl_path)
            return data
        last_err = RuntimeError(f"curl exited {r.returncode}: {r.stderr.strip()}")
    for attempt in range(4):
        if attempt > 0:
            time.sleep(min(2 ** attempt, 8))
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (compatible; artifact-fetcher/1.0)"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, ssl.SSLError, TimeoutError, ConnectionError, http.client.HTTPException, OSError) as exc:
            last_err = exc
    raise RuntimeError(f"failed to fetch {url}: {last_err}")


if mode != "--force" and all_present():
    if mode != "--ensure":
        print("[LOCKER-ROLLO] Official NIST Round 1 LOCKER reference files (reference/native/) are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[LOCKER-ROLLO] Third-party NIST LOCKER reference files (reference/native/) not found locally (or --force requested); downloading from upstream...", file=sys.stderr)

zip_bytes = None
last_err = None
for url in ZIP_URLS:
    print(f"  -> Fetching {url}", file=sys.stderr)
    try:
        data = fetch_bytes(url)
        got = hashlib.sha256(data).hexdigest()
        if got != ZIP_SHA256:
            raise RuntimeError(f"SHA-256 mismatch: expected {ZIP_SHA256}, got {got}")
        zip_bytes = data
        break
    except Exception as exc:
        last_err = exc

if zip_bytes is None:
    raise SystemExit(f"[LOCKER-ROLLO] ERROR: failed to download LOCKER.zip ({last_err})")

aux_data: dict[str, bytes] = {}
print(f"  -> Fetching LOCKER CLI bridge files @ {BENCH_COMMIT[:8]} (cli.c, cli_stdin.h, kem_c_api.cpp)", file=sys.stderr)
for fn, exp in AUX_FILES.items():
    try:
        raw = fetch_bytes(f"{BENCH_BASE}/{fn}")
    except Exception as exc:
        raise SystemExit(f"[LOCKER-ROLLO] ERROR: failed to download {fn} ({exc})")
    got = hashlib.sha256(raw).hexdigest()
    if got != exp:
        raise SystemExit(f"[LOCKER-ROLLO] ERROR: SHA-256 mismatch for {fn}: expected {exp}, got {got}")
    aux_data[fn] = raw


def transform_source(rel_path: str, raw: bytes) -> bytes:
    txt = raw.replace(b"\r\n", b"\n").decode("latin1")
    txt = txt.replace("LOCKER", "KEMSET").replace("locker", "kemset")
    if rel_path.endswith("src/ffi/ffi_vec.cpp"):
        # Fix upstream in-place aliasing bug in ffi_vec_add(se, se, e1, ...):
        # zeroing `o` before `add(o, v1, v2)` clobbers `v1` when `&o == &v1`.
        txt, n = re.subn(
            r"(void\s+ffi_vec_add\s*\([^)]*\)\s*\{\s*)ffi_vec_set_zero\s*\(\s*o\s*,\s*size\s*\)\s*;\s*",
            r"\1",
            txt,
            count=1,
        )
        if n != 1:
            raise SystemExit(f"[LOCKER-ROLLO] ERROR: failed to apply aliasing fix in {rel_path}")
    return txt.encode("latin1")


stage = os.path.join(tmp, "native")
os.makedirs(stage, exist_ok=True)
recorded_hashes: dict[str, str] = {}

with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
    for info in zf.infolist():
        if info.is_dir() or not info.filename.startswith("Reference_Implementation/LOCKER-"):
            continue
        rel = info.filename.replace("LOCKER", "KEMSET").replace("locker", "kemset")
        content = transform_source(rel, zf.read(info.filename))
        dst = os.path.join(stage, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "wb") as f:
            f.write(content)
        recorded_hashes[rel] = hashlib.sha256(content).hexdigest()

for fn, raw in aux_data.items():
    content = transform_source(fn, raw) if fn != "cli_stdin.h" else raw.replace(b"\r\n", b"\n")
    dst = os.path.join(stage, fn)
    with open(dst, "wb") as f:
        f.write(content)
    recorded_hashes[fn] = hashlib.sha256(content).hexdigest()

for fn, content, mode_bits in (
    ("Makefile", MAKEFILE_CONTENT, 0o644),
    ("build_cli.sh", BUILD_CLI_CONTENT, 0o755),
):
    dst = os.path.join(stage, fn)
    with open(dst, "wb") as f:
        f.write(content)
    os.chmod(dst, mode_bits)
    recorded_hashes[fn] = hashlib.sha256(content).hexdigest()

stage_manifest = os.path.join(stage, ".sha256_manifest.json")
with open(stage_manifest, "w", encoding="utf-8") as f:
    json.dump(recorded_hashes, f, indent=2, sort_keys=True)

got_manifest_sha = sha256_file(stage_manifest)
if got_manifest_sha != MANIFEST_SHA256:
    raise SystemExit(
        f"[LOCKER-ROLLO] ERROR: SHA-256 mismatch for generated .sha256_manifest.json: expected {MANIFEST_SHA256}, got {got_manifest_sha}"
    )

os.makedirs(os.path.dirname(native_dir), exist_ok=True)
if os.path.exists(native_dir):
    shutil.rmtree(native_dir)
shutil.move(stage, native_dir)

print("[LOCKER-ROLLO] Official NIST LOCKER reference files (KEMSET-I..IX) downloaded and verified successfully.", file=sys.stderr)
PY
