#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify external reference code
#                         for SBC-VOLE cryptanalysis and protocol-repair suites.
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
    echo "[SBC-VOLE] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[SBC-VOLE] ERROR: unknown option: $MODE" >&2
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
import os
import shutil
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request

art_dir = sys.argv[1]
mode = sys.argv[2]
tmp = sys.argv[3]

IACR_URL = "https://artifacts.iacr.org/files/asiacrypt/2025/asiacrypt-2025-a10.zip"
IACR_SHA256 = "e275d529c750bad077ca7bbe0ff355e9a183a4038df08d7869ef7ca95d3447a7"

BCON_COMMIT = "cfbde48414baacf51fc7c74f275190881f037d32"
BCON_BASE = f"https://raw.githubusercontent.com/B-Con/crypto-algorithms/{BCON_COMMIT}"
# Expected SHA-256 after LF normalization
BCON_FILES = {
    "aes.c": "f59be01f535478a16264a82f304802a627362ca27328c35efd9af816ea0e0d45",
    "aes.h": "e0c9a81b9951b558356de49e863341be38680c6cdc702306dd8ce9c0a9e8c032",
    "sha256.c": "4cbc93d30c28cee2a056d6f52cbce6e66ff9450f1d712ee4f59761b9b4486fad",
    "sha256.h": "a946e621775fd8251624e77d19f25407c4b6074ab780027e5797b77c7acf79b5",
}

ZIP_TARGETS = [
    os.path.join(art_dir, "key_recovery", "vendor", "asiacrypt-2025-a10.zip"),
    os.path.join(art_dir, "protocol_repair", "vendor", "asiacrypt-2025-a10.zip"),
]
BCON_DIRS = [
    os.path.join(art_dir, "key_recovery", "vendor", "crypto-algorithms"),
    os.path.join(art_dir, "protocol_repair", "crypto-algorithms"),
]


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def all_present() -> bool:
    for p in ZIP_TARGETS:
        if not os.path.isfile(p) or sha256_file(p) != IACR_SHA256:
            return False
    for d in BCON_DIRS:
        for fn, exp in BCON_FILES.items():
            p = os.path.join(d, fn)
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
    raise SystemExit(f"[SBC-VOLE] ERROR: failed to download {url} ({last_err})")


if mode != "--force" and all_present():
    if mode != "--ensure":
        print("[SBC-VOLE] All third-party reference files are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[SBC-VOLE] Third-party reference files not found locally (or --force requested); downloading from upstream...", file=sys.stderr)

print(f"  -> Fetching {IACR_URL}", file=sys.stderr)
zip_data = fetch_bytes(IACR_URL)
got_zip = hashlib.sha256(zip_data).hexdigest()
if got_zip != IACR_SHA256:
    raise SystemExit(f"[SBC-VOLE] ERROR: SHA-256 mismatch for asiacrypt-2025-a10.zip: expected {IACR_SHA256}, got {got_zip}")
stage_zip = os.path.join(tmp, "asiacrypt-2025-a10.zip")
with open(stage_zip, "wb") as f:
    f.write(zip_data)

stage_bcon = os.path.join(tmp, "crypto-algorithms")
os.makedirs(stage_bcon, exist_ok=True)
print(f"  -> Fetching B-Con/crypto-algorithms @ {BCON_COMMIT[:8]} (aes.c, aes.h, sha256.c, sha256.h)", file=sys.stderr)
for fn, exp in BCON_FILES.items():
    raw = fetch_bytes(f"{BCON_BASE}/{fn}").replace(b"\r\n", b"\n")
    got = hashlib.sha256(raw).hexdigest()
    if got != exp:
        raise SystemExit(f"[SBC-VOLE] ERROR: SHA-256 mismatch for {fn}: expected {exp}, got {got}")
    with open(os.path.join(stage_bcon, fn), "wb") as f:
        f.write(raw)

for p in ZIP_TARGETS:
    os.makedirs(os.path.dirname(p), exist_ok=True)
    shutil.copy2(stage_zip, p)

for d in BCON_DIRS:
    os.makedirs(d, exist_ok=True)
    for fn in BCON_FILES:
        shutil.copy2(os.path.join(stage_bcon, fn), os.path.join(d, fn))

print("[SBC-VOLE] Third-party reference files downloaded and verified successfully.", file=sys.stderr)
PY
