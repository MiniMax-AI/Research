#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Stream and SHA-256-verify CAESAR SUPERCOP reference
#                         implementations for COLM_0 and COLM_127.
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
    echo "[COLM] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[COLM] ERROR: unknown option: $MODE" >&2
        usage >&2
        exit 2
        ;;
esac

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT INT TERM

python3 - "$ART_DIR" "$MODE" "$TMP_DIR" <<'PY'
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request

art_dir = sys.argv[1]
mode = sys.argv[2]
tmp = sys.argv[3]
sf_dir = os.path.join(art_dir, "splice_forgery")

PRIMARY_URL = "https://bench.cr.yp.to/supercop/supercop-20261018.tar.xz"
INDEX_URL = "https://bench.cr.yp.to/supercop.html"

EXPECTED_SHA256 = {
    "colm0v1_ref/crypto_aead.h": "6792d19251f5233abffcd15e5a4ac2be4d83a93901773864024d749bb4c70f2d",
    "colm0v1_ref/encrypt.c": "01d6829978ecd7ec249fb30bbfe996d4e6766bbd1bf6788b0e8b7d0d6ea6cce7",
    "colm0v1_ref/module.h": "c98da79bfde46e1f246f9d7daf507ed75f264625652ec76155e88a522561dc39",
    "colm127v1_ref/crypto_aead.h": "1401715267ac159ed45959fb36ca5ebb721426244107dd532f8c30f86c2e7229",
    "colm127v1_ref/encrypt.c": "52a72e0f4f26a850a5ae26e3eb01a328c3110a5df0e0cb03156ed727dcf131a7",
    "colm127v1_ref/module.h": "1862ad21f40a003cabfcaee8047286143d4f81c26988949f3f32ffbef48c32eb",
    "colm0v1_aesni4/crypto_aead.h": "7191d55c0ad5647d8a38eb7286effb3eee82a5b5cf58ef78e072e5257107d916",
    "colm0v1_aesni4/encrypt-optpar4.c": "8bee4b80e435e00abecd3da5da9f834622076dcc0b1f58dec2ae624b881ad3c3",
    "colm0v1_aesni8/crypto_aead.h": "7191d55c0ad5647d8a38eb7286effb3eee82a5b5cf58ef78e072e5257107d916",
    "colm0v1_aesni8/encrypt-optpar8.c": "1740830bb576ea393551a0e72c0de4bfe3a1cc7acd59554de6ebfe2c58e9cc52",
}


def make_supercop_header(abytes: int, with_protos: bool) -> bytes:
    hdr = (
        f"#define CRYPTO_KEYBYTES 16\n"
        f"#define CRYPTO_NSECBYTES 0\n"
        f"#define CRYPTO_NPUBBYTES 8\n"
        f"#define CRYPTO_ABYTES {abytes}\n"
        f"#define CRYPTO_NOOVERLAP 1\n"
    ).encode("ascii")
    if with_protos:
        hdr += (
            b"int crypto_aead_encrypt(unsigned char *c,unsigned long long *clen,"
            b"const unsigned char *m,unsigned long long mlen,"
            b"const unsigned char *ad,unsigned long long adlen,"
            b"const unsigned char *nsec,const unsigned char *npub,const unsigned char *k);\n"
            b"int crypto_aead_decrypt(unsigned char *m,unsigned long long *mlen,"
            b"unsigned char *nsec,const unsigned char *c,unsigned long long clen,"
            b"const unsigned char *ad,unsigned long long adlen,"
            b"const unsigned char *npub,const unsigned char *k);\n"
        )
    return hdr


UPSTREAM_MAP = {
    "crypto_aead/colm0v1/ref/encrypt.c": "colm0v1_ref/encrypt.c",
    "crypto_aead/colm0v1/ref/module.h": "colm0v1_ref/module.h",
    "crypto_aead/colm127v1/ref/encrypt.c": "colm127v1_ref/encrypt.c",
    "crypto_aead/colm127v1/ref/module.h": "colm127v1_ref/module.h",
    "crypto_aead/colm0v1/aesni4/encrypt-optpar4.c": "colm0v1_aesni4/encrypt-optpar4.c",
    "crypto_aead/colm0v1/aesni8/encrypt-optpar8.c": "colm0v1_aesni8/encrypt-optpar8.c",
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def all_present() -> bool:
    for rel, exp in EXPECTED_SHA256.items():
        p = os.path.join(sf_dir, rel)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    return True


if mode != "--force" and all_present():
    if mode != "--ensure":
        print("[COLM] All third-party CAESAR SUPERCOP COLM files are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[COLM] Third-party CAESAR SUPERCOP files (splice_forgery/colm*) not found locally (or --force requested); streaming from upstream...", file=sys.stderr)


def extract_colm_from_stream(stream) -> dict[str, bytes]:
    out: dict[str, bytes] = {}
    with tarfile.open(fileobj=stream, mode="r|xz") as tf:
        for member in tf:
            parts = member.name.split("/", 1)
            if len(parts) < 2:
                continue
            sub = parts[1]
            if sub in UPSTREAM_MAP and member.isfile():
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                out[UPSTREAM_MAP[sub]] = fobj.read()
                if len(out) == len(UPSTREAM_MAP):
                    break
    return out


def resolve_urls() -> list[str]:
    urls = [PRIMARY_URL]
    curl_bin = shutil.which("curl")
    html = ""
    if curl_bin:
        r = subprocess.run(
            [curl_bin, "-fsSL", "--retry", "2", "--connect-timeout", "15", "--max-time", "30", INDEX_URL],
            capture_output=True,
            text=True,
        )
        if r.returncode == 0:
            html = r.stdout
    if not html:
        try:
            req = urllib.request.Request(INDEX_URL, headers={"User-Agent": "Mozilla/5.0 (compatible; artifact-fetcher/1.0)"})
            with urllib.request.urlopen(req, timeout=20) as r:
                html = r.read().decode("utf-8", "replace")
        except Exception:
            pass
    for m in re.findall(r"(?:https://bench\.cr\.yp\.to/)?supercop/supercop-\d{8}\.tar\.xz", html):
        u = m if m.startswith("http") else f"https://bench.cr.yp.to/{m.lstrip('/')}"
        if u not in urls:
            urls.append(u)
    return urls


collected: dict[str, bytes] = {}
last_err: Exception | None = None
curl_bin = shutil.which("curl")

for url in resolve_urls():
    print(f"  -> Streaming COLM prefix from {url}", file=sys.stderr)
    if curl_bin:
        for attempt in range(3):
            if attempt > 0:
                time.sleep(2 ** attempt)
            proc = subprocess.Popen(
                [
                    curl_bin, "-fsSL",
                    "--retry", "2",
                    "--connect-timeout", "20",
                    "-A", "Mozilla/5.0 (compatible; artifact-fetcher/1.0)",
                    url,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            try:
                assert proc.stdout is not None
                collected = extract_colm_from_stream(proc.stdout)
                if len(collected) == len(UPSTREAM_MAP):
                    break
            except Exception as exc:
                last_err = exc
                collected.clear()
            finally:
                if proc.poll() is None:
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
        if len(collected) == len(UPSTREAM_MAP):
            break

    for attempt in range(3):
        if attempt > 0:
            time.sleep(2 ** attempt)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (compatible; artifact-fetcher/1.0)"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                collected = extract_colm_from_stream(resp)
            if len(collected) == len(UPSTREAM_MAP):
                break
        except Exception as exc:
            last_err = exc
            collected.clear()
    if len(collected) == len(UPSTREAM_MAP):
        break

if len(collected) != len(UPSTREAM_MAP):
    raise SystemExit(f"[COLM] ERROR: failed to stream COLM files from SUPERCOP ({last_err})")

# Synthesize SUPERCOP-generated crypto_aead.h headers
collected["colm0v1_ref/crypto_aead.h"] = make_supercop_header(16, with_protos=True)
collected["colm127v1_ref/crypto_aead.h"] = make_supercop_header(528, with_protos=True)
collected["colm0v1_aesni4/crypto_aead.h"] = make_supercop_header(16, with_protos=False)
collected["colm0v1_aesni8/crypto_aead.h"] = make_supercop_header(16, with_protos=False)

for rel, exp in EXPECTED_SHA256.items():
    got = hashlib.sha256(collected[rel]).hexdigest()
    if got != exp:
        raise SystemExit(f"[COLM] ERROR: SHA-256 mismatch for {rel}: expected {exp}, got {got}")

for rel, content in collected.items():
    p = os.path.join(tmp, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "wb") as f:
        f.write(content)
for subdir in ("colm0v1_ref", "colm127v1_ref", "colm0v1_aesni4", "colm0v1_aesni8"):
    dst_dir = os.path.join(sf_dir, subdir)
    if os.path.exists(dst_dir):
        shutil.rmtree(dst_dir)
    shutil.move(os.path.join(tmp, subdir), dst_dir)

print("[COLM] Third-party CAESAR SUPERCOP COLM files streamed and SHA-256 verified successfully.", file=sys.stderr)
PY
