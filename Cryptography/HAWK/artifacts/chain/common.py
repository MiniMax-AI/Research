"""Shared helpers for the HAWK-256 end-to-end key-recovery chain.

Every path in this artefact is relative to the `artifacts/` directory; missing
third-party reference files (`hawk_recover/`, `kat/hawk_ref/`) are fetched and
SHA-256-verified automatically on first run via `fetch_third_party.sh`.
"""
from __future__ import annotations

import gzip
import os
import re
import resource
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))          # artifacts/chain
ROOT = os.path.dirname(HERE)                                # artifacts
DATA = os.path.join(ROOT, "data")
WORK = os.environ.get("HAWK_WORK", os.path.join(ROOT, ".work"))

subprocess.run(["bash", os.path.join(ROOT, "fetch_third_party.sh"), "--ensure"], check=True)

T0 = time.time()
_R0 = resource.getrusage(resource.RUSAGE_SELF)


def cpu() -> float:
    """CPU seconds of this process AND of every subprocess it spawned (flatter)."""
    a = resource.getrusage(resource.RUSAGE_SELF)
    b = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ((a.ru_utime - _R0.ru_utime) + (a.ru_stime - _R0.ru_stime)
            + b.ru_utime + b.ru_stime)


def log(msg: str) -> None:
    print(f"[{time.time()-T0:8.1f}s wall, {cpu():8.1f}s cpu] {msg}", flush=True)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def resolve(name: str) -> str:
    """Find `name` in data/ (preferring a .gz), else in the work directory."""
    for cand in (os.path.join(DATA, name + ".gz"), os.path.join(DATA, name),
                 os.path.join(WORK, name), name):
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(name)


def read_mat(path: str):
    """Read a Sage-style `[[ ... ] ... ]` integer matrix, transparently gunzipping."""
    from sage.all import ZZ, matrix
    txt = _open(path).read()
    rows = re.findall(r"\[([^\[\]]*)\]", txt)
    return matrix(ZZ, [[int(x) for x in r.split()] for r in rows])


def write_mat(path: str, B) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        f.write("[")
        for i in range(B.nrows()):
            f.write("[" + " ".join(str(int(B[i, j])) for j in range(B.ncols())) + "]\n")
        f.write("]\n")


def ensure_work() -> str:
    os.makedirs(WORK, exist_ok=True)
    return WORK
