#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify the official NIST SHA-3
#                         Round 1 Waterfall C reference implementation.
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
    echo "[Waterfall] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[Waterfall] ERROR: unknown option: $MODE" >&2
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
import os
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
ref_dir = os.path.join(art_dir, "unified_collision", "reference")

ZIP_URLS = [
    "https://csrc.nist.rip/groups/ST/hash/sha-3/Round1/documents/Waterfall.zip",
    "https://raw.githubusercontent.com/ethz-spylab/cryptanalysis-benchmark/d243a0a03d7aa367472411114f140d0aace583cd/algorithms/T2_waterfall/environment/algorithm_src/Waterfall.zip",
]
ZIP_SHA256 = "eb6e7e805f6bb737a687992e5f57607cf61955cfd7b198ae3cf8514a8a9a0bf7"

EXPECTED_OUT = {
    "hash.h": "d8c0cc8b79cff24baffb27001ae555f01a3d614bca0656ab2895c082c747a632",
    "Waterfall_ref.c": "a2d6c36d6de05af71551f569f08d2fce9419849bde00b25fede1f22d6153501c",
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def all_present() -> bool:
    for fn, exp in EXPECTED_OUT.items():
        p = os.path.join(ref_dir, fn)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    return True


def fetch_bytes(url: str, timeout: int = 90) -> bytes:
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
        print("[Waterfall] Official NIST SHA-3 Round 1 Waterfall reference files are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[Waterfall] Third-party NIST Waterfall reference files (unified_collision/reference/{Waterfall_ref.c,hash.h}) not found locally (or --force requested); downloading from upstream...", file=sys.stderr)

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
    raise SystemExit(f"[Waterfall] ERROR: failed to download Waterfall.zip ({last_err})")

with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
    wf_h = zf.read("Waterfall/Reference_Implementation/Waterfall.h").replace(b"\r\n", b"\n")
    wf_c = (
        zf.read("Waterfall/Reference_Implementation/Waterfall.c")
        .replace(b"\r\n", b"\n")
        .replace(b'#include "Waterfall.h"', b'#include "hash.h"')
    )

outputs = {"hash.h": wf_h, "Waterfall_ref.c": wf_c}
for fn, exp in EXPECTED_OUT.items():
    got = hashlib.sha256(outputs[fn]).hexdigest()
    if got != exp:
        raise SystemExit(f"[Waterfall] ERROR: SHA-256 mismatch for extracted {fn}: expected {exp}, got {got}")

os.makedirs(ref_dir, exist_ok=True)
for fn, content in outputs.items():
    stage = os.path.join(tmp, fn)
    with open(stage, "wb") as f:
        f.write(content)
    shutil.copy2(stage, os.path.join(ref_dir, fn))

print("[Waterfall] Official NIST Waterfall reference files downloaded and SHA-256 verified successfully.", file=sys.stderr)
PY
