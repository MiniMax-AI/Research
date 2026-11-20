#!/usr/bin/env python3
"""Design-point sweep for the TSM linearisation.

For a grid of (m, n, r, P) it builds a synthetic ideal-LRPC instance
(synth.sage), runs the runtime-parameterised attack kernel
(./tsm) and records the *measured* linear-algebra invariants:

    nrows  = C(P-1,r) * m          equations over GF(2^n)      (predicted)
    ncols  = C(P,r)   * 2^r        unknowns  over GF(2^n)      (by construction)
    rank   = nrows                  (full row rank: the prediction under test)
    nullity= ncols - rank           (= C(P,r)2^r - C(P-1,r)m when rank = nrows)
    exit 0 => the Frobenius pinning stage isolated exactly one Plucker line

It also records the analytic crossover
    P* = m*r/(m - 2^r)   (for m > 2^r)
beyond which the linearisation is *overdetermined* (nullity < 0 for every larger
P), which is what limits the method as r grows.

usage: design_sweep.py [--quick] [--out design_sweep.json]
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys, tempfile, time
from math import comb, gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def analytic(m, r, p):
    ncols = comb(p, r) * 2 ** r
    nrows = comb(p - 1, r) * m
    nrel = comb(p - 1, r) * 2 ** r
    return nrows, ncols, ncols - nrows, nrel


def crossover(m, r):
    if m <= 2 ** r:
        return None
    return m * r / (m - 2 ** r)


def run_case(m, n, r, p, seed, threads, timeout=3600):
    with tempfile.TemporaryDirectory() as tmp:
        inst = Path(tmp) / "i.json"
        pr = subprocess.run([sys.executable, str(HERE / "synth.sage"), "--m", str(m),
                             "--n", str(n), "--r", str(r), "--seed", str(seed),
                             "--out", str(inst), "--quiet"],
                            capture_output=True, text=True, timeout=timeout)
        if pr.returncode != 0:
            return {"ok": False, "stage": "synth", "err": pr.stderr.strip()[-300:]}
        d = json.loads(inst.read_text())
        out = Path(tmp) / "p.bin"
        cmd = [str(HERE / "tsm"), "--m", str(m), "--n", str(n), "--r", str(r),
               "--pairs", str(p), "--ring", ",".join(map(str, d["ring"])),
               "--field", ",".join(map(str, d["field"])),
               "--pk", d["public_key"], "--ct", d["encapsulation"],
               "--out", str(out), "--stats"]
        t0 = time.monotonic()
        try:
            pr = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                                env={**os.environ, "OMP_NUM_THREADS": str(threads)})
        except subprocess.TimeoutExpired:
            return {"ok": False, "stage": "timeout", "seconds": timeout}
        wall = round(time.monotonic() - t0, 2)
        stats, note = {}, ""
        for line in pr.stderr.splitlines():
            if line.startswith("{"):
                stats = json.loads(line)
            elif "kernel" in line or "nullity" in line:
                note = line
        row = {"m": m, "n": n, "r": r, "P": p, "gcd": int(gcd(m, n)), "seed": seed,
               "rc": pr.returncode, "seconds": wall, "stats": stats}
        pnrows, pncols, pnull, nrel = analytic(m, r, p)
        row["pred_nrows"], row["pred_ncols"] = pnrows, pncols
        row["pred_nullity"], row["nrel"] = pnull, nrel
        row["pinned"] = (pr.returncode == 0)
        if stats:
            row["nullity_matches"] = (stats.get("nullity") == pnull)
            row["rank_matches"] = (stats.get("rank") == pnrows and stats.get("nrows") == pnrows
                                   and stats.get("ncols") == pncols)
        else:
            row["note"] = (note or pr.stderr.strip().splitlines()[-1:] or [""])[0] \
                if isinstance(note or "", str) else ""
        return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", default=str(HERE / "design_sweep.json"))
    args = ap.parse_args()

    cases = []
    # (1) gcd(m,n) = 1 control vs gcd(m,n) > 1, at a clean design point r=3, P=6, m=15
    cases += [(15, 11, 3, 6, 1), (15, 9, 3, 6, 1), (15, 10, 3, 6, 1), (15, 12, 3, 6, 1)]
    if not args.quick:
        # (2) P-sweep at LOCKER-I (KEMSET-I: m=71, n=83, r=5)
        for p in (7, 8, 9, 10):
            cases.append((71, 83, 5, p, 101 + p))
        # (3) a second r=5 set with a different (m,n) at its own design point
        cases += [(101, 83, 5, 7, 7), (79, 101, 5, 8, 8)]
        # (4) r=6 boundary: is the linearisation overdetermined past P*?
        cases += [(41, 11, 6, 8, 3), (41, 11, 6, 9, 3)]

    rows = []
    for (m, n, r, p, seed) in cases:
        row = run_case(m, n, r, p, seed, args.threads)
        row["crossover_P"] = crossover(m, r)
        rows.append(row)
        print(json.dumps(row, separators=(",", ":")), flush=True)
    json.dump(rows, open(args.out, "w"), indent=1)
    print("SUMMARY rows=%d pinned=%d" % (len(rows), sum(1 for x in rows if x.get("pinned"))))


if __name__ == "__main__":
    main()
