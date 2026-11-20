#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify the official NIST SHA-3
#                         Round 1 Blender C reference (with 64-bit LP64
#                         <stdint.h> width bindings) and EACirc C++ reference.
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
    echo "[Blender] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[Blender] ERROR: unknown option: $MODE" >&2
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
ref_dir = os.path.join(art_dir, "reference")

ZIP_URLS = [
    "https://web.archive.org/web/20141112163008id_/http://csrc.nist.gov/groups/ST/hash/sha-3/Round1/documents/Blender.zip",
    "https://raw.githubusercontent.com/ethz-spylab/cryptanalysis-benchmark/d243a0a03d7aa367472411114f140d0aace583cd/algorithms/T2_blender/environment/algorithm_src/Blender.zip",
]
ZIP_SHA256 = "ea77dd90acb65a10a76e6746f02dec4986a6f7a32cbdad820f48e7cd1fd2e8e4"

EACIRC_COMMIT = "b92d96ad16679c24a9402b0d33378d1f318db23a"
EACIRC_BASE = f"https://raw.githubusercontent.com/crocs-muni/eacirc-streams/{EACIRC_COMMIT}/streams/hash/sha3/hash_functions/Blender"
EACIRC_FILES = {
    "Blender_sha3.cpp": "60decaaae6e7bc8d26023a83e2a6b2daaaae9b51eda2c0d26390455a3ffc44d2",
    "Blender_sha3.h": "333c211ca35c568c4f17c43ee1fb8942788b85cda22e6c94f7a5eb0c4f39048b",
}

EXPECTED_OUT = {
    "api.h": "b5327cc524ef353e30555843b2c96e7a704d43563d8957db4d4452d1a36b5053",
    "core.h": "e5d6325d394a26076647366cf14bdabe652d733f0c1fa0152d20f6143013cd51",
    "core.c": "3ea0e2f2d3c5c407a7fe9aa8f8fab4c7edd6ded047b0cbba9f39b9e62885578e",
    **EACIRC_FILES,
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
        print("[Blender] Official NIST SHA-3 Round 1 Blender & EACirc reference files are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

print("[Blender] Third-party NIST Blender reference files (reference/{core.c,core.h,api.h,Blender_sha3.*}) not found locally (or --force requested); downloading from upstream...", file=sys.stderr)

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
    raise SystemExit(f"[Blender] ERROR: failed to download Blender.zip ({last_err})")

with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
    api_h = zf.read("Blender/Reference Implementation/SHA3api_ref.h").replace(b"\r\n", b"\n").decode("latin1")
    core_h = zf.read("Blender/Reference Implementation/Blender.h").replace(b"\r\n", b"\n").decode("latin1")
    core_c = zf.read("Blender/Reference Implementation/Blender.c").replace(b"\r\n", b"\n").decode("latin1")

# Apply deterministic 64-bit LP64 <stdint.h> fixed-width bindings & C11 MSVC-intrinsics cleanup
api_h_bytes = (
    api_h.replace("#define SHA3API_H\n", "#define SHA3API_H\n#include <stdint.h>\n", 1)
    .replace("typedef unsigned long long DataLength;", "typedef uint64_t DataLength;")
    .replace("typedef unsigned long DataLength;", "typedef uint32_t DataLength;")
    .replace('#include "Blender.h"', '#include "core.h"')
).encode("latin1")

core_h_bytes = (
    core_h.replace("#define BLENDER_H\n", "#define BLENDER_H\n#include <stdint.h>\n", 1)
    .replace("#define Uint64 unsigned long long", "#define Uint64 uint64_t")
    .replace("unsigned long llw;", "uint32_t llw;")
    .replace("unsigned long hlw;", "uint32_t hlw;")
    .replace("#define Uint32 unsigned long", "#define Uint32 uint32_t")
    .replace("#define Uint16 unsigned short", "#define Uint16 uint16_t")
    .replace("#define Uint8 unsigned char", "#define Uint8 uint8_t")
).encode("latin1")

idx = core_c.index("BitSequence bigBuffer[1000];")
core_c = core_c[:idx] + "HashReturn testhash (void) { return SUCCESS; }\nHashReturn timehash (void) { return SUCCESS; }\n"
core_c_bytes = (
    core_c.replace("#include <stdio.h>\n", "#include <stdio.h>\n#include <stdint.h>\n", 1)
    .replace('#include "SHA3api_ref.h"', '#include "api.h"')
    .replace(
        'char *Author = "Colin Bradbury";\nchar *Algorithm = "Blender";',
        '#ifdef VERBOSE\nchar *Author = "Colin Bradbury";\nchar *Algorithm = "Blender";\n#endif',
    )
    .replace("unsigned long long", "uint64_t")
).encode("latin1")

outputs = {
    "api.h": api_h_bytes,
    "core.h": core_h_bytes,
    "core.c": core_c_bytes,
}

print(f"  -> Fetching crocs-muni/eacirc-streams @ {EACIRC_COMMIT[:8]} (Blender_sha3.cpp, Blender_sha3.h)", file=sys.stderr)
for fn in EACIRC_FILES:
    outputs[fn] = fetch_bytes(f"{EACIRC_BASE}/{fn}").replace(b"\r\n", b"\n")

for fn, exp in EXPECTED_OUT.items():
    got = hashlib.sha256(outputs[fn]).hexdigest()
    if got != exp:
        raise SystemExit(f"[Blender] ERROR: SHA-256 mismatch for {fn}: expected {exp}, got {got}")

os.makedirs(ref_dir, exist_ok=True)
for fn, content in outputs.items():
    stage = os.path.join(tmp, fn)
    with open(stage, "wb") as f:
        f.write(content)
    shutil.copy2(stage, os.path.join(ref_dir, fn))

print("[Blender] Official NIST Blender & EACirc reference files downloaded and SHA-256 verified successfully.", file=sys.stderr)
PY
