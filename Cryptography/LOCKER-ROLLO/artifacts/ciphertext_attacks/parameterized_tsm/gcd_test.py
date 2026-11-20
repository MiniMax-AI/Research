#!/usr/bin/env python3
"""End-to-end TSM check on synthetic instances, sweeping gcd(m,n).

For each (m, n, r) it generates a synthetic ideal-LRPC instance, runs the
attack kernel (./tsm) plus the reference recovery (tsm_recover.py)
and verifies that the *recovered shared secret* equals the one the instance was
built with, and that the recovered support equals the true error support E.

The point of the sweep is the role of gcd(m,n).  Tests whether
gcd(m,n) > 1 (where R = GF(2^m)[X]/(f) is no longer a field, and
Q = GF(2^n)[T]/(g) splits) might block the attack.  Here that is tested directly.

usage: gcd_test.py [--pairs-auto] m:n:r[:P] ...
"""
from __future__ import annotations
import json, os, subprocess, sys, tempfile, time
from math import comb, gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def suggested_pairs(m, r):
    for p in range(r + 1, 40):
        if comb(p, r) * 2 ** r - comb(p - 1, r) * m == comb(p - 1, r):
            return p
    for p in range(r + 1, 40):
        nullity = comb(p, r) * 2 ** r - comb(p - 1, r) * m
        if nullity > 0 and comb(p - 1, r) * 2 ** r >= 2 * nullity:
            return p
    for p in range(r + 1, 40):
        if comb(p - 1, r) * m >= comb(p, r) * 2 ** r - 1:
            return p
    raise RuntimeError("no feasible P")


def one(m, n, r, p, seed, threads=2):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inst = tmp / "i.json"
        pr = subprocess.run([sys.executable, str(HERE / "synth.sage"), "--m", str(m),
                             "--n", str(n), "--r", str(r), "--seed", str(seed),
                             "--out", str(inst), "--quiet"], capture_output=True, text=True)
        if pr.returncode != 0:
            return {"m": m, "n": n, "r": r, "P": p, "gcd": gcd(m, n), "ok": False,
                    "stage": "synth", "err": pr.stderr.strip()[-200:]}
        d = json.loads(inst.read_text())
        out = tmp / "p.bin"
        t0 = time.monotonic()
        pr = subprocess.run([str(HERE / "tsm"), "--m", str(m), "--n", str(n), "--r", str(r),
                             "--pairs", str(p), "--ring", ",".join(map(str, d["ring"])),
                             "--field", ",".join(map(str, d["field"])),
                             "--pk", d["public_key"], "--ct", d["encapsulation"],
                             "--out", str(out), "--stats"], capture_output=True, text=True,
                            env={**os.environ, "OMP_NUM_THREADS": str(threads)})
        lin = round(time.monotonic() - t0, 3)
        stats = {}
        for line in pr.stderr.splitlines():
            if line.startswith("{"):
                stats = json.loads(line)
        row = {"m": m, "n": n, "r": r, "P": p, "gcd": int(gcd(m, n)), "seed": seed,
               "lin_rc": pr.returncode, "lin_seconds": lin, "stats": stats}
        if pr.returncode != 0:
            row.update({"ok": False, "stage": "tsm",
                        "err": (pr.stderr.strip().splitlines() or [""])[-1]})
            return row
        rr = subprocess.run([sys.executable, str(HERE / "tsm_recover.py"), str(inst), str(out),
                             "--m", str(m), "--n", str(n), "--r", str(r), "--pairs", str(p),
                             "--ring", ",".join(map(str, d["ring"])),
                             "--field", ",".join(map(str, d["field"])), "--no-verify"],
                            capture_output=True, text=True)
        if rr.returncode != 0:
            row.update({"ok": False, "stage": "recover",
                        "err": (rr.stderr.strip().splitlines() or [""])[-1]})
            return row
        res = json.loads([l for l in rr.stdout.splitlines() if l.strip()][-1])
        got = sorted(int(x, 16) for x in res["support"])
        want = sorted(d["support_echelon"])
        row.update({"ok": res["real_shared_secret"] == d["shared_secret"],
                    "digest_verified": res["verified_digest"],
                    "support_exact": got == want,
                    "digest_matches_published": res.get("matches_published")})
        return row


def main():
    threads = int(os.environ.get("OMP_NUM_THREADS", "2"))
    specs = sys.argv[1:]
    if not specs:
        specs = ["15:11:3", "15:9:3", "15:10:3", "15:12:3", "15:14:3",
                 "21:11:3", "21:9:3", "35:11:4", "35:22:4"]
    rows = []
    for spec in specs:
        parts = [int(x) for x in spec.split(":")]
        m, n, r = parts[0], parts[1], parts[2]
        p = parts[3] if len(parts) > 3 else suggested_pairs(m, r)
        for seed in (1, 2):
            row = one(m, n, r, p, seed, threads)
            rows.append(row)
            print(json.dumps(row, separators=(",", ":")), flush=True)
    json.dump(rows, open(HERE / "gcd_test.json", "w"), indent=1)
    print("SUMMARY ok=%d/%d" % (sum(1 for x in rows if x.get("ok")), len(rows)))


if __name__ == "__main__":
    main()
