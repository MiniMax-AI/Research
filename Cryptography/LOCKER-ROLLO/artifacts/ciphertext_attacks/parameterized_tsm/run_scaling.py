#!/usr/bin/env python3
"""Scaling harness for the Transposed Support-Minor (TSM) attack.

For each KEMSET parameter set it generates fresh instances with the unmodified
reference implementation (./gen<SET>), runs the attack on the *public*
data only (pk, encapsulation), and checks that the recovered shared secret
equals the one produced by the reference encapsulation.  The recovery is also
self-validated against d = SHA512(m) inside the ciphertext.

usage: run_scaling.py [--sets I,II,IV,VII,VIII] [--trials 1] [--pairs P]
"""
from __future__ import annotations
import argparse, json, os, re, statistics, subprocess, sys, tempfile, time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# KEMSET -> (m, n, d, r, security, ring modulus exponents, GF(2^m) modulus exponents)
PARAMS = {
    "I":    (71, 83, 7, 5, 128, [83, 7, 4, 2, 0],        [71, 6, 0]),
    "II":   (101, 83, 7, 5, 192, [83, 7, 4, 2, 0],       [101, 7, 6, 1, 0]),
    "III":  (107, 89, 8, 6, 256, [89, 38, 0],            [107, 9, 7, 4, 0]),
    "IV":   (79, 101, 7, 5, 128, [101, 7, 6, 1, 0],      [79, 9, 0]),
    "V":    (97, 103, 8, 6, 192, [103, 9, 0],            [97, 6, 0]),
    "VI":   (107, 103, 8, 6, 256, [103, 9, 0],           [107, 9, 7, 4, 0]),
    "VII":  (83, 149, 8, 5, 128, [149, 10, 9, 7, 0],     [83, 7, 4, 2, 0]),
    "VIII": (101, 149, 8, 5, 192, [149, 10, 9, 7, 0],    [101, 7, 6, 1, 0]),
    "IX":   (109, 157, 9, 6, 256, [157, 6, 5, 2, 0],     [109, 5, 4, 2, 0]),
}


def suggested_pairs(m: int, r: int) -> int:
    """Choose the number P of Frobenius shifts.

    Unknowns  : C(P,r) * 2^r over GF(2^n)   (maximal minors of the r x 2P matrix)
    Equations : C(P-1,r) * m  over GF(2^n)  (6x6 minors of [X|U], one per (r+1)-subset
                                             containing pair 0, expanded over GF(2^m))
    nullity   = C(P,r)2^r - C(P-1,r) m      (when the system has full row rank)

    Preferred ("clean") design point: nullity == C(P-1,r), i.e. the GF(2^n)-kernel is
    exactly the space of Frobenius orbits, so the Frobenius relations are consistent and
    isolate a single line.  That is  P*2^r/(P-r) == m+1.
    Otherwise take the smallest P whose Frobenius relations outnumber the kernel
    coordinates by at least a factor 2 (verified empirically for all r=5 sets).
    """
    from math import comb
    for p in range(r + 1, 40):
        if comb(p, r) * 2 ** r - comb(p - 1, r) * m == comb(p - 1, r):
            return p
    for p in range(r + 1, 40):
        nullity = comb(p, r) * 2 ** r - comb(p - 1, r) * m
        nrel = comb(p - 1, r) * 2 ** r
        if nullity > 0 and nrel >= 2 * nullity:
            return p
    raise RuntimeError("no feasible P")


def generate(setname: str):
    exe = HERE / f"gen{setname}"
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            d[k.strip()] = v.strip().removeprefix("0x")
    return d


def attack(setname: str, pk: str, ct: str, pairs: int, scratch: Path, timeout: int):
    m, n, _d, r, _sec, ring, field = PARAMS[setname]
    plk = scratch / "p.bin"
    cmd = [str(HERE / "tsm"), "--m", str(m), "--n", str(n), "--r", str(r),
           "--pairs", str(pairs), "--ring", ",".join(map(str, ring)),
           "--field", ",".join(map(str, field)), "--pk", pk, "--ct", ct,
           "--out", str(plk), "--stats"]
    t0 = time.monotonic()
    pr = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "12")})
    lin_time = time.monotonic() - t0
    stats = {}
    for line in pr.stderr.splitlines():
        if line.startswith("{"):
            stats = json.loads(line)
    if pr.returncode != 0:
        return {"ok": False, "stage": "tsm", "rc": pr.returncode,
                "err": pr.stderr.strip().splitlines()[-1] if pr.stderr.strip() else "",
                "lin_seconds": round(lin_time, 2), "stats": stats}
    t1 = time.monotonic()
    inst = scratch / "inst.json"
    inst.write_text(json.dumps({"public_key": pk, "encapsulation": ct}))
    rr = subprocess.run([sys.executable, str(HERE / "tsm_recover.py"), str(inst), str(plk),
                         "--m", str(m), "--n", str(n), "--r", str(r), "--pairs", str(pairs),
                         "--ring", ",".join(map(str, ring)), "--field", ",".join(map(str, field)),
                         "--no-verify"], capture_output=True, text=True, timeout=timeout)
    rec_time = time.monotonic() - t1
    if rr.returncode != 0:
        return {"ok": False, "stage": "recover", "err": rr.stderr.strip().splitlines()[-1],
                "lin_seconds": round(lin_time, 2), "rec_seconds": round(rec_time, 2),
                "stats": stats}
    res = json.loads([l for l in rr.stdout.splitlines() if l.strip()][-1])
    res.update({"ok": True, "lin_seconds": round(lin_time, 2), "rec_seconds": round(rec_time, 2),
                "stats": stats})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="I,II,IV,VII,VIII")
    ap.add_argument("--trials", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=0, help="override P")
    ap.add_argument("--timeout", type=int, default=3600)
    args = ap.parse_args()
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        for s in args.sets.split(","):
            m, n, d, r, sec, ring, field = PARAMS[s]
            pairs = args.pairs or suggested_pairs(m, r)
            for t in range(args.trials):
                gt = generate(s)
                t0 = time.monotonic()
                try:
                    res = attack(s, gt["pk"], gt["ct"], pairs, scratch, args.timeout)
                except subprocess.TimeoutExpired:
                    res = {"ok": False, "stage": "timeout"}
                wall = time.monotonic() - t0
                good = res.get("ok") and res.get("real_shared_secret") == gt["ss"] \
                    and res.get("verified_digest")
                row = {"set": s, "m": m, "n": n, "d": d, "r": r, "claimed_security": sec,
                       "pairs": pairs, "trial": t, "success": bool(good),
                       "wall_seconds": round(wall, 2),
                       "lin_seconds": res.get("lin_seconds"), "rec_seconds": res.get("rec_seconds"),
                       "stats": res.get("stats", {}),
                       "e2_nonzero": bool(int(gt["e2"], 16) != 0) if "e2" in gt else None,
                       "note": res.get("err", "")}
                rows.append(row)
                print(json.dumps(row, separators=(",", ":")), flush=True)
    ok = sum(1 for r in rows if r["success"])
    print(f"SUMMARY {ok}/{len(rows)} instances broken")
    json.dump(rows, open(HERE / "scaling_results.json", "w"), indent=1)


if __name__ == "__main__":
    main()
