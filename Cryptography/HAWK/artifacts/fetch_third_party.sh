#!/usr/bin/env bash
# =============================================================================
# fetch_third_party.sh -- Download and SHA-256-verify external HAWK reference C
#                         implementation (hawk-sign/dev), hawk_recover package
#                         (anthropics/cryptography-research-demo), and AGPS20
#                         classical sieve cost table (jschanck/eprint-2019-1161).
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
    echo "[HAWK] ERROR: unexpected extra arguments: $*" >&2
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
        echo "[HAWK] ERROR: unknown option: $MODE" >&2
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
import tarfile
import time
import urllib.error
import urllib.request

art_dir = sys.argv[1]
mode = sys.argv[2]
tmp = sys.argv[3]

HAWK_REF_COMMIT = "1b9fef52559273fe7b40fe3e22968eaedd3a4c2a"
HAWK_REF_URL = f"https://codeload.github.com/hawk-sign/dev/tar.gz/{HAWK_REF_COMMIT}"
HAWK_REF_TAR_SHA256 = "9ea8a39e57c44a5d3999721d26c2932ffb8ae9871ba02dcb60fb4a7a10c06254"

ANTH_COMMIT = "fa01c398b42bb7d94eadff75a42dfb45da484457"
ANTH_URL = f"https://codeload.github.com/anthropics/cryptography-research-demo/tar.gz/{ANTH_COMMIT}"
ANTH_TAR_SHA256 = "3f318c7124e1fe93e6450e77d3702c3c4555b3cd24470dbc93006ed6b7364fed"

AGPS20_COMMIT = "1ec49e5b5f4e0063e6c470a6c75bf425483987ca"
AGPS20_CSV_URL = f"https://raw.githubusercontent.com/jschanck/eprint-2019-1161/{AGPS20_COMMIT}/data/cost-estimate-list_decoding-classical.csv"
AGPS20_CSV_SHA256 = "5fa1bdee7d904d00cd743a1afb2b74cb32529f734736c444fd3a408552db4cba"

HAWK_REF_EXPECTED = {
    "hawk.h": "84bc9a6c59a7f35a346a866304c65549910c9e246ab08d64fc1357835f00f5b1",
    "hawk_config.h": "bc2db43211d9732f3178d80cbfb687f4cae25319427d3699c5a174c31080076f",
    "hawk_inner.h": "219cb6b38d447c68476edaea4e7b1060464defc86b9a88f6d818324425c2d1af",
    "hawk_kgen.c": "81cf4f3db9e2446d9d4572f2c5df43759d8f39fd55ee38160be0cead4d585eab",
    "hawk_sign.c": "fb60aa13426f3b9fbc0d4c689ddee50a3e3f349cbd215d7085609b4616118d55",
    "hawk_vrfy.c": "b70b86a53c10d7fa2ab75022aab44d3c7452c4143eaa86cd6ae98cc06dac77ae",
    "modq.h": "ea99bbca5a7e2cd7bfe7bf8bcf354c765b8f93d9259e92b4a457a9a534d93a1e",
    "ng_config.h": "15dac9fb98cf951ecccfd2c8fbec9617922a1456d63d305c571af68b1ba6f30c",
    "ng_fxp.c": "0c95f39c01b0505a56e5f8589cd6f8b3d826916d035f757a2bd84d1eec77392c",
    "ng_hawk.c": "f7663bccf91c0f7b2ec143672e4869589b977c2f0651c7b8525403f6cd4888df",
    "ng_inner.h": "9aaa6bee50d15aecdb48e554f7d1dc8b0c20b8fc42615db39b710565d75b21ad",
    "ng_mp31.c": "32e70e2e8d501a5e909bd678018d3c0cb27a6503cdb211aa720e154a50d537be",
    "ng_ntru.c": "e24833b64cbf0b796c1e13031bad64a30b85c3b4c3e74964f47b05e6e1d38b06",
    "ng_poly.c": "7bbede59fbd863235a17f4f3ff7c7fb03caede28d0b843422c0d9138af4ab3d9",
    "ng_zint31.c": "c39d9486d7359714c0c22dff2324cc11ea8bc016d4dae6e62a787e6baf565f66",
    "sha3.c": "f715556701d4ccb3f18661e92596b3ee8e745428a396a5fcc1f4ea342e066223",
    "sha3.h": "d7570e988e67dc02e018429c1669173aaf8d916354c3d2152e5f049278ba77aa",
}

HAWK_SHIM_C_SHA256 = "b71ffe36d711c37a0d803735782f4f912ca73f726420fb0be6b15d44aa5eae45"
SHIM_PY_PATCHED_SHA256 = "6e7c6bd516ca1700118e13a1c0159f7d785c562ee0ade6a824b46ecefbd4ede0"

HAWK_RECOVER_UNMODIFIED = {
    "__init__.py": "99a6d2c7cc078ce555dc7d9c48718c39d18ae23603f0cdf534a3f130ee4eb4ac",
    "cocycle.py": "f2ce99ee09aae97d05d451d63facb3bc13f05fa6f2df981e33550826b7b70b96",
    "decode.py": "5fc3df80ec02d080e33f46c9e49313627fcd0cad8f1650f29e83c59ebcc313bc",
    "encode.py": "d6cdeaeed5511ec5cf428a9939f5493ad9aa1e36de1e531ef20b1a4ad3458003",
    "infra/__init__.py": "0c22c87cfb3783d42bac2f793fa3c3de2872e5edb3111641153295999ea24f83",
    "infra/checkpoint.py": "547f6aa7ac9910126e58398379bd9d5142409ebfc4a14b974efeed650505dd19",
    "infra/sieve_backends/__init__.py": "2e09cb0a12d02679494e7f029f12947b57b373fdf887917a8080d2e2bd10aac5",
    "infra/sieve_backends/bgj.py": "e40276d6acbaaf3afcdb3e8092722788487ba349bf2f0774a855108e4eca341f",
    "intertwiner.py": "e3aa7ffb4c664710e7123825ea45aafbb44ad06f36dfec08e1400f6d42a0006e",
    "params.py": "4ed35a4d008d4a5c913ef6d89650b3c23067db4095adea355928b2582b10eda8",
    "ring.py": "97e9657f638e1b9654f1e583f57f79c966521a3b1c17cb1b888c36b16a5e78a9",
    "svp.py": "fb8fb61542c98b7ac90beb5991c3b6e609a7361ea7cd3985eccdd3107799e333",
    "tower.py": "63a90989305c0ba8477b66f15bf04c72baa3be1bd6b507505a3a7ffce5dec5d7",
    "types.py": "be3c2ed447d90af42613e8027791fba40d087acae3ead7dbb208e3c980329f59",
}

kat_dir = os.path.join(art_dir, "kat")
hawk_ref_dir = os.path.join(kat_dir, "hawk_ref")
hawk_shim_c = os.path.join(kat_dir, "hawk_shim.c")
hawk_recover_dir = os.path.join(art_dir, "hawk_recover")
agps20_dir = os.path.join(art_dir, "agps20")
agps20_csv = os.path.join(agps20_dir, "cost-estimate-list_decoding-classical.csv")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def patch_shim_py(raw: bytes) -> bytes:
    txt = raw.replace(b"\r\n", b"\n").decode("utf-8")
    txt = txt.replace(
        "import hashlib\n",
        "import atexit\nimport hashlib\nimport os\nimport shutil\n",
        1,
    )
    ensure_block = '''_KAT_DIR = Path(__file__).resolve().parents[2] / "kat"
_CACHED_SHIM: Path | None = None


def _ensure_shim() -> Path:
    global _CACHED_SHIM
    prebuilt = _KAT_DIR / "hawk_shim"
    if prebuilt.is_file() and os.access(prebuilt, os.X_OK):
        return prebuilt
    if _CACHED_SHIM is not None and _CACHED_SHIM.is_file():
        return _CACHED_SHIM

    build_script = _KAT_DIR / "build_shim.sh"
    require(build_script.is_file(), f"build_shim.sh missing at {build_script}")
    tmp_dir = Path(tempfile.mkdtemp(prefix="hawk_shim_"))
    atexit.register(shutil.rmtree, tmp_dir, True)
    out = tmp_dir / "hawk_shim"
    r = subprocess.run(
        ["bash", str(build_script), str(out)],
        capture_output=True,
        text=True,
    )
    require(
        r.returncode == 0 and out.is_file(),
        f"failed to compile hawk_shim (rc={r.returncode}): {(r.stdout + r.stderr).strip()}",
    )
    _CACHED_SHIM = out
    return out'''
    if '_SHIM = Path(__file__).parent.parent / "ref_impl" / "hawk_shim"' not in txt:
        raise SystemExit("[HAWK] ERROR: unexpected upstream hawk_recover/infra/shim.py structure")
    txt = txt.replace('_SHIM = Path(__file__).parent.parent / "ref_impl" / "hawk_shim"', ensure_block, 1)
    txt = txt.replace('require(_SHIM.exists(), f"hawk_shim missing at {_SHIM}")', "shim_bin = _ensure_shim()")
    txt = txt.replace("str(_SHIM)", "str(shim_bin)")
    return txt.encode("utf-8")


def ref_present() -> bool:
    for fn, exp in HAWK_REF_EXPECTED.items():
        p = os.path.join(hawk_ref_dir, fn)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    return True


def anth_present() -> bool:
    if not os.path.isfile(hawk_shim_c) or sha256_file(hawk_shim_c) != HAWK_SHIM_C_SHA256:
        return False
    for rel, exp in HAWK_RECOVER_UNMODIFIED.items():
        p = os.path.join(hawk_recover_dir, rel)
        if not os.path.isfile(p) or sha256_file(p) != exp:
            return False
    shim_py = os.path.join(hawk_recover_dir, "infra", "shim.py")
    if not os.path.isfile(shim_py) or sha256_file(shim_py) != SHIM_PY_PATCHED_SHA256:
        return False
    return True


def agps20_present() -> bool:
    return os.path.isfile(agps20_csv) and sha256_file(agps20_csv) == AGPS20_CSV_SHA256


def all_present() -> bool:
    return ref_present() and anth_present() and agps20_present()


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
    raise SystemExit(f"[HAWK] ERROR: failed to download {url} ({last_err})")


if mode != "--force" and all_present():
    if mode != "--ensure":
        print("[HAWK] All third-party HAWK C reference, hawk_recover, and AGPS20 files are already present and SHA-256 verified.", file=sys.stderr)
    sys.exit(0)

need_ref = mode == "--force" or not ref_present()
need_anth = mode == "--force" or not anth_present()
need_agps = mode == "--force" or not agps20_present()

print("[HAWK] Missing third-party reference files (or --force requested); downloading from upstream...", file=sys.stderr)

if need_ref:
    print(f"  -> Fetching hawk-sign/dev @ {HAWK_REF_COMMIT[:8]}", file=sys.stderr)
    ref_tar_bytes = fetch_bytes(HAWK_REF_URL)
    got_ref_sha = hashlib.sha256(ref_tar_bytes).hexdigest()
    if got_ref_sha != HAWK_REF_TAR_SHA256:
        raise SystemExit(f"[HAWK] ERROR: SHA-256 mismatch for hawk-sign/dev tarball: expected {HAWK_REF_TAR_SHA256}, got {got_ref_sha}")

    stage_ref = os.path.join(tmp, "hawk_ref")
    os.makedirs(stage_ref, exist_ok=True)
    ref_prefix = f"dev-{HAWK_REF_COMMIT}/src/"
    with tarfile.open(fileobj=io.BytesIO(ref_tar_bytes), mode="r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile() or not m.name.startswith(ref_prefix):
                continue
            fn = m.name[len(ref_prefix):]
            if fn in HAWK_REF_EXPECTED:
                data = tf.extractfile(m).read()
                got = hashlib.sha256(data).hexdigest()
                if got != HAWK_REF_EXPECTED[fn]:
                    raise SystemExit(f"[HAWK] ERROR: SHA-256 mismatch for {fn}")
                with open(os.path.join(stage_ref, fn), "wb") as f:
                    f.write(data)

    for fn in HAWK_REF_EXPECTED:
        if not os.path.isfile(os.path.join(stage_ref, fn)):
            raise SystemExit(f"[HAWK] ERROR: missing {fn} in hawk-sign/dev tarball")

    os.makedirs(kat_dir, exist_ok=True)
    if os.path.exists(hawk_ref_dir):
        shutil.rmtree(hawk_ref_dir)
    shutil.move(stage_ref, hawk_ref_dir)

if need_anth:
    print(f"  -> Fetching anthropics/cryptography-research-demo @ {ANTH_COMMIT[:8]}", file=sys.stderr)
    anth_tar_bytes = fetch_bytes(ANTH_URL)
    got_anth_sha = hashlib.sha256(anth_tar_bytes).hexdigest()
    if got_anth_sha != ANTH_TAR_SHA256:
        raise SystemExit(f"[HAWK] ERROR: SHA-256 mismatch for cryptography-research-demo tarball: expected {ANTH_TAR_SHA256}, got {got_anth_sha}")

    stage_rec = os.path.join(tmp, "hawk_recover")
    stage_shim_c = os.path.join(tmp, "hawk_shim.c")
    os.makedirs(stage_rec, exist_ok=True)
    anth_prefix = f"cryptography-research-demo-{ANTH_COMMIT}/HAWK/hawk_recover/"
    with tarfile.open(fileobj=io.BytesIO(anth_tar_bytes), mode="r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile() or not m.name.startswith(anth_prefix):
                continue
            rel = m.name[len(anth_prefix):]
            if rel == "ref_impl/hawk_shim.c":
                data = tf.extractfile(m).read()
                got = hashlib.sha256(data).hexdigest()
                if got != HAWK_SHIM_C_SHA256:
                    raise SystemExit("[HAWK] ERROR: SHA-256 mismatch for hawk_shim.c")
                with open(stage_shim_c, "wb") as f:
                    f.write(data)
            elif rel == "infra/shim.py":
                patched = patch_shim_py(tf.extractfile(m).read())
                got_shim_py = hashlib.sha256(patched).hexdigest()
                if got_shim_py != SHIM_PY_PATCHED_SHA256:
                    raise SystemExit(
                        f"[HAWK] ERROR: SHA-256 mismatch for patched infra/shim.py: expected {SHIM_PY_PATCHED_SHA256}, got {got_shim_py}"
                    )
                dst = os.path.join(stage_rec, rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                with open(dst, "wb") as f:
                    f.write(patched)
            elif rel in HAWK_RECOVER_UNMODIFIED:
                data = tf.extractfile(m).read()
                got = hashlib.sha256(data).hexdigest()
                if got != HAWK_RECOVER_UNMODIFIED[rel]:
                    raise SystemExit(f"[HAWK] ERROR: SHA-256 mismatch for hawk_recover/{rel}")
                dst = os.path.join(stage_rec, rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                with open(dst, "wb") as f:
                    f.write(data)

    if not os.path.isfile(stage_shim_c):
        raise SystemExit("[HAWK] ERROR: missing hawk_shim.c in cryptography-research-demo tarball")
    for rel in list(HAWK_RECOVER_UNMODIFIED) + ["infra/shim.py"]:
        if not os.path.isfile(os.path.join(stage_rec, rel)):
            raise SystemExit(f"[HAWK] ERROR: missing hawk_recover/{rel}")

    os.makedirs(kat_dir, exist_ok=True)
    shutil.copy2(stage_shim_c, hawk_shim_c)
    if os.path.exists(hawk_recover_dir):
        shutil.rmtree(hawk_recover_dir)
    shutil.move(stage_rec, hawk_recover_dir)

if need_agps:
    print(f"  -> Fetching jschanck/eprint-2019-1161 @ {AGPS20_COMMIT[:7]} (classical sieve cost table)", file=sys.stderr)
    agps_bytes = fetch_bytes(AGPS20_CSV_URL)
    got_agps_sha = hashlib.sha256(agps_bytes).hexdigest()
    if got_agps_sha != AGPS20_CSV_SHA256:
        raise SystemExit(f"[HAWK] ERROR: SHA-256 mismatch for AGPS20 CSV: expected {AGPS20_CSV_SHA256}, got {got_agps_sha}")
    os.makedirs(agps20_dir, exist_ok=True)
    with open(agps20_csv, "wb") as f:
        f.write(agps_bytes)

print("[HAWK] Third-party HAWK C reference, hawk_recover, and AGPS20 files downloaded and SHA-256 verified successfully.", file=sys.stderr)
PY
