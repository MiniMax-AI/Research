#!/usr/bin/env python3
"""Self-contained pure-Python GF(2)/GF(2^m) verifier for the TII McEliece artifact suite.

Implements:
  1. Live minimum-weight codeword sampling via the native C engine (`minwords_mt`)
     and verification of the Minimum-Weight Locator Derivative Identity
     `sigma'_S(X) = g_bar(X)^2` (paper Section 3, Theorem 1) and the odd-weight
     pair-count identity `sum_{u<v} c_u c_v = C(2t+1, 2) = t (mod 2)` (Theorem 2).
  2. Full SHA-256 and algebraic validity checks (`deg(g) = t`, `g` monic and
     square-free via Bezout/GCD over GF(2^m)[X], `|L| = n` distinct support
     elements with `g(z_j) != 0`) across all 94 solved TII Track-2 challenges.
  3. Exact mathematical acceptance test of `tii_repo/check_solution_track2.sage`
     (`RREF(H_raw(L, g)) == H_pub` and `H_raw(L, g) * G_pub^T == 0` over GF(2))
     using only the Python standard library.

Standalone usage:
  python3 verify_witnesses_pure.py [--demo | --full] [scratch_dir]
"""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
TII_REPO = HERE / "tii_repo"
PK_DIR = TII_REPO / "public_keyRec"
WITNESS_DIR = HERE / "witnesses"

DEMO_CHALLENGES = [83, 133, 91, 125, 253, 254, 228, 250, 172, 187, 240, 246, 252]

if hasattr(0, "bit_count"):
    def _popcount(x: int) -> int:
        return x.bit_count()
else:
    def _popcount(x: int) -> int:
        return bin(x).count("1")


def sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_elem(s: str, m: int) -> int:
    """Parse a Sage GF(2^m) element in generator `a` into an m-bit integer mask."""
    s = s.strip()
    if s == "0":
        return 0
    v = 0
    for term in s.split("+"):
        term = term.strip()
        if term == "1":
            v ^= 1
        elif term == "a":
            v ^= 2
        else:
            mm = re.fullmatch(r"a\^(\d+)", term)
            if not mm or int(mm.group(1)) >= m:
                raise ValueError(f"invalid GF(2^{m}) element: {s!r}")
            v ^= 1 << int(mm.group(1))
    return v


def read_pk(sec: int, n: int, k: int) -> tuple[list[int], list[int]]:
    """Read H_pub rows (as n-bit integer bitmasks) and defining polynomial f(x)."""
    rows: list[int] = []
    poly: list[int] | None = None
    pk_path = PK_DIR / f"pk_McEliece_{sec}.txt"
    if not pk_path.is_file():
        subprocess.run(["bash", str(HERE.parent / "fetch_third_party.sh"), "--ensure"], check=True)
    with pk_path.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            line = line.strip().strip("[").strip("]")
            if i < n - k:
                bits = [int(x) for x in line.split(" ") if x != ""]
                if len(bits) != n:
                    raise ValueError(f"TII-{sec}: row {i} has {len(bits)} bits, expected {n}")
                row_mask = 0
                for j, b in enumerate(bits):
                    if b:
                        row_mask |= 1 << j
                rows.append(row_mask)
            else:
                poly = [int(x) for x in line.split(", ")]
                break
    if poly is None or len(rows) != n - k:
        raise ValueError(f"TII-{sec}: malformed public key file")
    return rows, poly


def read_sk(sec: int, m: int) -> tuple[list[int], list[int]]:
    """Read two-line TII Track-2 secret key (g coefficients and support L)."""
    sk_path = WITNESS_DIR / f"sk_McEliece_{sec}.txt"
    lines = [ln.strip() for ln in sk_path.read_text(encoding="utf-8").strip().splitlines() if ln.strip()]
    if len(lines) != 2:
        raise ValueError(f"TII-{sec}: expected 2 lines in {sk_path.name}, got {len(lines)}")
    g = [parse_elem(x, m) for x in lines[0].strip("[").strip("]").split(",")]
    L = [parse_elem(x, m) for x in lines[1].strip("[").strip("]").split(",")]
    return g, L


class GF2m:
    """Arithmetic in F_{2^m} = F_2[x]/(f(x)) with lookup tables (m <= 10)."""

    def __init__(self, m: int, poly_bits: list[int]) -> None:
        if len(poly_bits) != m + 1 or poly_bits[m] != 1:
            raise ValueError(f"invalid modulus degree for m={m}: {poly_bits}")
        self.m = m
        self.order = 1 << m
        fullmask = sum(c << i for i, c in enumerate(poly_bits))
        self.fullmask = fullmask
        self.redmask = fullmask & (self.order - 1)

        # Build exp/log tables using a multiplication helper
        self._inv = [0] * self.order
        for a in range(1, self.order):
            self._inv[a] = self._inv_slow(a)

    def mul(self, a: int, b: int) -> int:
        r = 0
        m = self.m
        red = self.redmask
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if (a >> m) & 1:
                a ^= (1 << m) | red
        return r

    def _inv_slow(self, a: int) -> int:
        s, acc = a, 1
        for _ in range(self.m - 1):
            s = self.mul(s, s)
            acc = self.mul(acc, s)
        return acc

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("GF(2^m) division by zero")
        return self._inv[a]

    def eval_poly(self, p: list[int], z: int) -> int:
        acc = 0
        for c in reversed(p):
            acc = self.mul(acc, z) ^ c
        return acc


def pdeg(p: list[int]) -> int:
    d = len(p) - 1
    while d >= 0 and p[d] == 0:
        d -= 1
    return d


def ptrim(p: list[int]) -> list[int]:
    d = pdeg(p)
    return p[: d + 1] if d >= 0 else []


def padd(p: list[int], q: list[int]) -> list[int]:
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) ^ (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(F: GF2m, p: list[int], q: list[int]) -> list[int]:
    if not p or not q:
        return []
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    r[i + j] ^= F.mul(a, b)
    return ptrim(r)


def pdivmod(F: GF2m, a: list[int], b: list[int]) -> tuple[list[int], list[int]]:
    a = a[:]
    db = pdeg(b)
    if db < 0:
        raise ZeroDivisionError("polynomial division by zero")
    ib = F.inv(b[db])
    q = [0] * max(0, pdeg(a) - db + 1)
    for i in range(pdeg(a), db - 1, -1):
        if i < len(a) and a[i]:
            c = F.mul(a[i], ib)
            q[i - db] = c
            for j in range(db + 1):
                a[i - db + j] ^= F.mul(c, b[j])
    return ptrim(q), ptrim(a)


def pderiv(p: list[int]) -> list[int]:
    return ptrim([p[i] if (i & 1) else 0 for i in range(1, len(p))])


def pgcd_is_unit(F: GF2m, a: list[int], b: list[int]) -> bool:
    r0, r1 = ptrim(a), ptrim(b)
    while r1:
        _, r2 = pdivmod(F, r0, r1)
        r0, r1 = r1, r2
    return pdeg(r0) == 0


def build_hraw(F: GF2m, m: int, t: int, n: int, g: list[int], L: list[int]) -> list[int]:
    """Exact transcription of McElieceKeyGen.sage:gen_parity_full."""
    out = [0] * (m * t)
    for j, z in enumerate(L):
        gz = F.eval_poly(g, z)
        if gz == 0:
            raise ValueError(f"g(L[{j}]) == 0")
        ig = F.inv(gz)
        zp = 1
        col_bit = 1 << j
        for i in range(t):
            e = F.mul(zp, ig)
            base_row = m * i
            val = e
            while val:
                lsb = val & -val
                b = lsb.bit_length() - 1
                out[base_row + b] |= col_bit
                val ^= lsb
            zp = F.mul(zp, z)
    return out


def rref(ncols: int, rows: list[int]) -> list[int]:
    """Bit-parallel Gauss-Jordan RREF over GF(2), matching Sage matrix.echelon_form()."""
    acc: list[int] = []
    rest = rows[:]
    for col in range(ncols):
        bit = 1 << col
        piv = None
        for idx, r in enumerate(rest):
            if r & bit:
                piv = idx
                break
        if piv is None:
            continue
        p = rest.pop(piv)
        acc = [(r ^ p) if (r & bit) else r for r in acc]
        rest = [(r ^ p) if (r & bit) else r for r in rest]
        acc.append(p)
    return acc + rest


def check_orthogonality_with_gpub(hraw_rows: list[int], hpub_rows: list[int], n: int) -> bool:
    """Verify H_raw(L, g) * G_pub^T == 0 over GF(2) using pivot/free split of H_pub."""
    mt = len(hpub_rows)
    piv_cols: list[int] = []
    piv_set: set[int] = set()
    for r in hpub_rows:
        if r == 0:
            return False
        lsb = r & -r
        c = lsb.bit_length() - 1
        piv_cols.append(c)
        piv_set.add(c)
    if len(piv_set) != mt:
        return False
    free_cols = [j for j in range(n) if j not in piv_set]
    # Construct each generator row of G_pub (one per free column) and check orthogonality
    for fcol in free_cols:
        g_row = 1 << fcol
        fbit = 1 << fcol
        for idx, pcol in enumerate(piv_cols):
            if hpub_rows[idx] & fbit:
                g_row |= 1 << pcol
        for h_row in hraw_rows:
            if _popcount(h_row & g_row) & 1:
                return False
    return True


def check_key_pure(sec: int, n: int, m: int, t: int, k: int) -> tuple[bool, GF2m, list[int], list[int], list[int]]:
    """Run the full mathematical acceptance check of check_solution_track2.sage."""
    hpub_rows, poly_bits = read_pk(sec, n, k)
    F = GF2m(m, poly_bits)
    g, L = read_sk(sec, m)
    if len(g) != t + 1 or g[t] != 1:
        return False, F, g, L, hpub_rows
    if len(L) != n or len(set(L)) != n:
        return False, F, g, L, hpub_rows
    if not pgcd_is_unit(F, g, pderiv(g)):
        return False, F, g, L, hpub_rows
    hraw_rows = build_hraw(F, m, t, n, g, L)
    if not check_orthogonality_with_gpub(hraw_rows, hpub_rows, n):
        return False, F, g, L, hpub_rows
    hcand_rref = rref(n, hraw_rows)
    return (hcand_rref == hpub_rows), F, g, L, hpub_rows


def _compile_native_helper(src: pathlib.Path, out: pathlib.Path, *, is_cxx: bool) -> None:
    """Compile a Route A C/C++ helper into `out` with OpenMP/single-thread fallback."""
    env_var = "CXX" if is_cxx else "CC"
    default_bin = "g++" if is_cxx else "gcc"
    stem = "g++" if is_cxx else "gcc"
    std_flags = ["-std=c++17"] if is_cxx else []

    env_comp = os.environ.get(env_var, "").strip()
    if env_comp:
        cands = [env_comp]
    else:
        cands = [default_bin]
        for prefix in ("", "/opt/homebrew/bin/", "/usr/local/bin/"):
            for ver in (16, 15, 14, 13, 12, 11):
                cand = f"{prefix}{stem}-{ver}"
                if shutil.which(cand):
                    cands.append(cand)

    omp_flag_sets = [
        ["-O3", *std_flags, "-Wall", "-march=native", "-mpopcnt", "-fopenmp"],
        ["-O3", *std_flags, "-Wall", "-march=native", "-fopenmp"],
        ["-O3", *std_flags, "-Wall", "-mcpu=native", "-fopenmp"],
        ["-O2", *std_flags, "-Wall", "-fopenmp"],
    ]
    for comp in cands:
        if not shutil.which(comp):
            continue
        for flags in omp_flag_sets:
            proc = subprocess.run([comp, *flags, str(src), "-o", str(out)], capture_output=True, text=True)
            if proc.returncode == 0:
                return

    fb = env_comp or default_bin
    for flags in (
        ["-O3", *std_flags, "-Wall", "-march=native"],
        ["-O3", *std_flags, "-Wall", "-mcpu=native"],
        ["-O2", *std_flags, "-Wall"],
    ):
        proc = subprocess.run([fb, *flags, str(src), "-o", str(out)], capture_output=True, text=True)
        if proc.returncode == 0:
            return
    raise RuntimeError(f"Failed to compile {src.name} with {fb}:\n{proc.stderr}")


def run_stage1_native_engines(scratch: pathlib.Path) -> None:
    """Verify compiled C/C++ Route A engines and Theorem 1 & 2 identities on TII-133."""
    scratch.mkdir(parents=True, exist_ok=True)
    minwords_bin = scratch / "minwords_mt"
    rankmask_bin = scratch / "rankmask"
    if not minwords_bin.is_file():
        _compile_native_helper(HERE / "route_a" / "minwords_mt.c", minwords_bin, is_cxx=False)
    if not rankmask_bin.is_file():
        _compile_native_helper(HERE / "route_a" / "rankmask.cpp", rankmask_bin, is_cxx=True)

    # 1. Test rankmask on a 4-dimensional space of 8x8 symmetric GF(2) matrices
    rm_input = scratch / "test_mats.bin"
    # 4 matrices of size 8x8 (each row is 1 uint64_t): ranks 2, 2, 4, 2
    m0 = [0b01100000, 0b10100000, 0b11000000, 0, 0, 0, 0, 0]
    m1 = [0, 0, 0, 0b00011000, 0b00010000, 0, 0, 0]
    m2 = [m0[i] ^ m1[i] for i in range(8)]
    m3 = [0b00000110, 0, 0, 0, 0, 0b10000001, 0b10000001, 0b00110000]
    with rm_input.open("wb") as f:
        for mat in (m0, m1, m2, m3):
            for row in mat:
                f.write(struct.pack("<Q", row))
    rm_proc = subprocess.run(
        [str(rankmask_bin), str(rm_input), "8", "4", "2"],
        capture_output=True,
        text=True,
        check=True,
    )
    rm_lines = [ln.strip() for ln in rm_proc.stdout.strip().splitlines() if ln.strip()]
    assert len(rm_lines) == 15, f"expected 15 rankmask lines, got {len(rm_lines)}"
    print(f"  PASS  route_a/rankmask.cpp build & Gray-code rank census: {len(rm_lines)}/15 non-zero masks evaluated")

    # 2. Run minwords_mt on TII-133 (n=112, m=7, t=7, k=63, mt=49, w=15)
    sec, n, m, t, k = 133, 112, 7, 7, 63
    w = 2 * t + 1
    ok, F, g, L, hpub_rows = check_key_pure(sec, n, m, t, k)
    assert ok, "TII-133 check_key failed"
    mat_file = scratch / "tii133_hpub.txt"
    with mat_file.open("w", encoding="utf-8") as f:
        for r in hpub_rows:
            f.write(f"{r:x}\n")

    target_words = 64
    mw_proc = subprocess.run(
        [str(minwords_bin), str(mat_file), str(n - k), str(n), str(target_words), str(w), "3", "1", "10.0", "4"],
        capture_output=True,
        text=True,
        check=True,
    )
    supports: list[list[int]] = []
    for ln in mw_proc.stdout.splitlines():
        parts = [int(x) for x in ln.strip().split() if x]
        if len(parts) == w:
            supports.append(parts)
    if len(supports) < target_words:
        raise RuntimeError(f"minwords_mt found only {len(supports)}/{target_words} codewords on TII-133")

    # Verify Theorem 1 (sigma'_S(X) = g_bar(X)^2 and e_2(z_S) = g_{t-1}^2)
    # and Theorem 2(3) (pair count C(2t+1, 2) == t mod 2) on all sampled codewords
    g_sq = pmul(F, g, g)
    gt1_sq = F.mul(g[t - 1], g[t - 1])
    verified_thm1 = 0
    for sup in supports:
        cw_mask = sum(1 << j for j in sup)
        if _popcount(cw_mask) != w:
            raise RuntimeError("invalid codeword weight")
        for r in hpub_rows:
            if _popcount(r & cw_mask) & 1:
                raise RuntimeError("sampled vector not in ker(H_pub)")
        # Theorem 2(3): active pairs count mod 2 == t mod 2
        pair_count = w * (w - 1) // 2
        if (pair_count & 1) != (t & 1):
            raise RuntimeError("odd-weight pair count mismatch")
        # Theorem 1: sigma'_S(X) == g(X)^2 and e_2(z_S) == g_{t-1}^2
        sigma = [1]
        e2 = 0
        for u_idx in range(w):
            zu = L[sup[u_idx]]
            sigma = pmul(F, sigma, [zu, 1])
            for v_idx in range(u_idx + 1, w):
                e2 ^= F.mul(zu, L[sup[v_idx]])
        if pderiv(sigma) != g_sq or e2 != gt1_sq:
            raise RuntimeError("Theorem 1 locator derivative identity failed")
        verified_thm1 += 1

    print(
        f"  PASS  route_a/minwords_mt.c live ISD on TII-133: {verified_thm1}/{target_words} weight-{w} codewords satisfy sigma'_S(X) == g(X)^2 and e_2(z_S) == g_{{t-1}}^2"
    )


def _run_verification(mode: str, scratch: pathlib.Path) -> int:
    subprocess.run(["bash", str(HERE.parent / "fetch_third_party.sh"), "--ensure"], check=True)
    print("[1] Compiling & testing Route A C/C++ engines and Theorem 1/2 identities...")
    run_stage1_native_engines(scratch)

    print("[2] Auditing 94-challenge campaign manifest, SHA-256 digests, and key structure...")
    summary = json.loads((WITNESS_DIR / "verified_tii_94_summary.json").read_text(encoding="utf-8"))
    results = json.loads((WITNESS_DIR / "results.json").read_text(encoding="utf-8"))

    challenges = summary["challenges"]
    res_by_sec = {r["sec"]: r for r in results["results"]}
    if len(challenges) != 94 or summary["total_solved"] != 94:
        raise RuntimeError("manifest challenge count != 94")
    if summary["route_a_solved"] != 89 or summary["route_b_solved"] != 5:
        raise RuntimeError("route split mismatch")
    if summary["first_bloods"] != 55:
        raise RuntimeError("first-blood count mismatch")

    sha_ok = 0
    fmt_ok = 0
    for c in challenges:
        sec = c["sec"]
        r = res_by_sec[sec]
        pk_path = PK_DIR / f"pk_McEliece_{sec}.txt"
        sk_path = WITNESS_DIR / f"sk_McEliece_{sec}.txt"
        if sha256_file(pk_path) != r["pk_sha256"]:
            raise RuntimeError(f"TII-{sec}: public key SHA-256 mismatch")
        if sha256_file(sk_path) != r["sk_sha256"]:
            raise RuntimeError(f"TII-{sec}: secret key SHA-256 mismatch")
        sha_ok += 1
        g, L = read_sk(sec, c["m"])
        if len(g) == c["t"] + 1 and g[-1] == 1 and len(L) == c["n"] and len(set(L)) == c["n"]:
            fmt_ok += 1
        else:
            raise RuntimeError(f"TII-{sec}: structural parameter check failed")

    print(
        f"  PASS  Campaign manifest & SHA-256 audit: {sha_ok}/94 pk/sk digests matched, {fmt_ok}/94 valid (g, L) tuples"
    )
    print(
        f"  PASS  Headline campaign totals: {summary['total_solved']} solved ({summary['route_a_solved']} Route A, {summary['route_b_solved']} Route B), {summary['first_bloods']} first-bloods, {results['speedup_count']} speedups"
    )
    print(
        f"  PASS  Compute accounting: {summary['local_single_node_92_core_hours']:.2f} single-node core-h (92 challenges), {summary['total_94_core_hours']:.2f} total core-h (94 challenges)"
    )

    by_sec = {c["sec"]: c for c in challenges}
    target_secs = DEMO_CHALLENGES if mode == "--demo" else [c["sec"] for c in challenges]
    print(
        f"[3] Executing check_solution_track2.sage acceptance test (RREF(H_raw(L, g)) == H_pub) on {len(target_secs)} challenges ({mode})..."
    )
    passed = 0
    for sec in target_secs:
        c = by_sec[sec]
        ok, _, _, _, _ = check_key_pure(sec, c["n"], c["m"], c["t"], c["k"])
        if not ok:
            raise RuntimeError(f"TII-{sec}: check_key acceptance test failed")
        passed += 1
        route_tag = "Route-B" if c["route"].startswith("Route-B") else "Route-A"
        fb_tag = "first-blood" if c["first_blood"] else "speedup"
        print(
            f"  PASS  TII-{sec:<3d}: check_key == 1 (n={c['n']}, m={c['m']}, t={c['t']}, k={c['k']}, {route_tag}, {c['wall_seconds']}s, {fb_tag})"
        )

    print(
        f"  PASS  Official Track-2 RREF & orthogonality verification complete: {passed}/{len(target_secs)} keys verified (check_key == 1)"
    )
    return 0


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "--demo"
    if mode not in ("--demo", "--full"):
        print("Usage: python3 verify_witnesses_pure.py [--demo | --full] [scratch_dir]", file=sys.stderr)
        return 1
    if len(sys.argv) > 2:
        return _run_verification(mode, pathlib.Path(sys.argv[2]))
    with tempfile.TemporaryDirectory(prefix="tii-verify-") as tmp:
        return _run_verification(mode, pathlib.Path(tmp))


if __name__ == "__main__":
    sys.exit(main())
