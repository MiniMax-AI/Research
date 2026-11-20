#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify external TII McEliece
#                         Challenge scripts and Track-2 public keys.
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
    echo "[TII-McEliece] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[TII-McEliece] ERROR: unknown option: $MODE" >&2
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
import shutil
import ssl
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request

art_dir = sys.argv[1]
mode = sys.argv[2]
tmp = sys.argv[3]

UPSTREAM_COMMIT = "e502cd59f8b1b3f05bf251a4f66d2b17f5602a63"
UPSTREAM_URL = f"https://codeload.github.com/ElenaKirshanova/tii_decoding_challenge/tar.gz/{UPSTREAM_COMMIT}"
TARBALL_SHA256 = "b1f1bc120c4a052385b5b64be623fe1adef0886b6bd92d8371ead3c8af1c9131"

goppa_dir = os.path.join(art_dir, "goppa_key_recovery")
tii_repo = os.path.join(goppa_dir, "tii_repo")
pk_dir = os.path.join(tii_repo, "public_keyRec")
witness_dir = os.path.join(goppa_dir, "witnesses")

CORE_FILES = {
    "McElieceKeyGen.sage": "cd4bbe934f67783e22ceafce4ec05f7fab6b54966ac055024c35ca8a3ec95673",
    "check_solution_track2.sage": "f6e646beec13f7c032a3b7264f8a3854f05a056bda0f6b06a4ac538c8221758e",
    "key_rec_instances": "b6f0382011f26122caf06037c847587cc6dbcf491a4fdaa6188614c97e284bad",
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


results_json = os.path.join(witness_dir, "results.json")
with open(results_json, "r", encoding="utf-8") as f:
    expected_pks = {
        f"pk_McEliece_{r['sec']}.txt": r["pk_sha256"]
        for r in json.load(f)["results"]
    }


def all_present() -> bool:
    for rel, exp in CORE_FILES.items():
        p = os.path.join(tii_repo, rel)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    for fn, exp in expected_pks.items():
        p = os.path.join(pk_dir, fn)
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
    raise SystemExit(f"[TII-McEliece] ERROR: failed to download {url} ({last_err})")


if mode != "--force" and all_present():
    if mode != "--ensure":
        print("[TII-McEliece] All third-party TII challenge files (94 public keys + scripts) are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[TII-McEliece] Third-party reference files (goppa_key_recovery/tii_repo/) not found locally (or --force requested); downloading from upstream...", file=sys.stderr)
print(f"  -> Fetching ElenaKirshanova/tii_decoding_challenge @ {UPSTREAM_COMMIT[:8]}", file=sys.stderr)

tar_bytes = fetch_bytes(UPSTREAM_URL, timeout=120)

got_tar_sha = hashlib.sha256(tar_bytes).hexdigest()
if got_tar_sha != TARBALL_SHA256:
    raise SystemExit(
        f"[TII-McEliece] ERROR: SHA-256 mismatch for upstream tarball: expected {TARBALL_SHA256}, got {got_tar_sha}"
    )

prefix = f"tii_decoding_challenge-{UPSTREAM_COMMIT}/"
wanted_members = {prefix + k: k for k in CORE_FILES}
for fn in expected_pks:
    wanted_members[f"{prefix}public_keyRec/{fn}"] = f"public_keyRec/{fn}"

stage_repo = os.path.join(tmp, "tii_repo")
os.makedirs(os.path.join(stage_repo, "public_keyRec"), exist_ok=True)

with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tf:
    for member in tf.getmembers():
        rel = wanted_members.get(member.name)
        if rel is None or not member.isfile():
            continue
        fobj = tf.extractfile(member)
        if fobj is None:
            continue
        data = fobj.read()
        out_path = os.path.join(stage_repo, rel)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "wb") as out_f:
            out_f.write(data)

# Verify extracted core files and all 94 Track-2 public keys before installing
for rel, exp in CORE_FILES.items():
    p = os.path.join(stage_repo, rel)
    if not os.path.isfile(p) or sha256_file(p) != exp:
        raise SystemExit(f"[TII-McEliece] ERROR: verification failed for extracted {rel}")
for fn, exp in expected_pks.items():
    p = os.path.join(stage_repo, "public_keyRec", fn)
    if not os.path.isfile(p) or sha256_file(p) != exp:
        raise SystemExit(f"[TII-McEliece] ERROR: verification failed for extracted public_keyRec/{fn}")

if os.path.exists(tii_repo):
    shutil.rmtree(tii_repo)
shutil.move(stage_repo, tii_repo)

print("[TII-McEliece] Extracted and verified 94 official Track-2 public keys + challenge scripts into goppa_key_recovery/tii_repo/.", file=sys.stderr)
PY
