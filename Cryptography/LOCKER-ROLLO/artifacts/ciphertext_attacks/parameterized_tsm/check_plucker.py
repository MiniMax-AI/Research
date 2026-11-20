#!/usr/bin/env python3
"""Check that ./tsm recovers the *true* Pluecker vector.

The Pluecker coordinates of the TSM linearisation are
    p_{(Q,s)} = det [ x_i^{2^{q_j}} ]_{i,j},   x = a if s_j = 0 else b,
i.e. determinants of r x r matrices with entries in GF(2^n) -- they never leave
GF(2^n), so they can be computed from the ground-truth secret (a, b) alone, with
no reference to Q = GF(2^n)[T]/(g).  That makes this test valid even when
gcd(m,n) > 1 and Q is not a field.

It therefore isolates the *algebraic core* of the attack (linearisation +
Frobenius pinning) from the recovery stage, and answers whether gcd(m,n) = 1 is
needed for the core.

usage: check_plucker.py m n r P [seed]
"""
from __future__ import annotations
import itertools, json, os, struct, subprocess, sys, tempfile
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main():
    m, n, r, P = (int(x) for x in sys.argv[1:5])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    from sage.all import GF, PolynomialRing, matrix

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inst = tmp / "i.json"
        pr = subprocess.run([sys.executable, str(HERE / "synth.sage"), "--m", str(m),
                             "--n", str(n), "--r", str(r), "--seed", str(seed),
                             "--out", str(inst), "--quiet"], capture_output=True, text=True)
        assert pr.returncode == 0, pr.stderr[-400:]
        d = json.loads(inst.read_text())
        out = tmp / "p.bin"
        pr = subprocess.run([str(HERE / "tsm"), "--m", str(m), "--n", str(n), "--r", str(r),
                             "--pairs", str(P), "--ring", ",".join(map(str, d["ring"])),
                             "--field", ",".join(map(str, d["field"])),
                             "--pk", d["public_key"], "--ct", d["encapsulation"],
                             "--out", str(out), "--stats"], capture_output=True, text=True,
                            env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "2")})
        stats = {}
        for line in pr.stderr.splitlines():
            if line.startswith("{"):
                stats = json.loads(line)
        if pr.returncode != 0:
            print(json.dumps({"m": m, "n": n, "r": r, "P": P, "gcd": int(gcd(m, n)),
                              "tsm_rc": pr.returncode, "stats": stats,
                              "plucker_matches_truth": None}))
            return

        F2 = PolynomialRing(GF(2), "z")
        zz = F2.gen()
        L = GF(2 ** n, "a", modulus=sum(zz ** int(e) for e in d["ring"]))
        raw = out.read_bytes()
        ncols, rr, pp, nwords = struct.unpack("<IIII", raw[:16])
        body = raw[16:]
        got = []
        assert len(body) == ncols * 32, (len(body), ncols)   # Fe is always 4 words wide
        for i in range(ncols):
            w = struct.unpack_from("<4Q", body, i * 32)[:nwords]
            v = 0
            for k, word in enumerate(w):
                v |= word << (64 * k)
            got.append(L.from_integer(v))

        gen = L.gen()
        def as_field(row):          # binary polynomial -> element of GF(2^n)
            v = L(0)
            for t, bit in enumerate(row):
                if bit:
                    v += gen ** t
            return v
        a = [as_field(row) for row in d["truth"]["a"]]
        b = [as_field(row) for row in d["truth"]["b"]]
        combos = list(itertools.combinations(range(P), r))
        truth = []
        for Q in combos:
            for mask in range(1 << r):
                Mx = matrix(L, r, r,
                            lambda i, j: (b[i] if (mask >> j) & 1 else a[i]) ** (2 ** Q[j]))
                truth.append(Mx.det())
        assert len(truth) == ncols, (len(truth), ncols)

        scale = None
        mism = 0
        nonzero = 0
        for t, gv in zip(truth, got):
            if t:
                nonzero += 1
                s = gv / t
                if scale is None:
                    scale = s
                elif s != scale:
                    mism += 1
        print(json.dumps({"m": m, "n": n, "r": r, "P": P, "gcd": int(gcd(m, n)), "seed": seed,
                          "tsm_rc": 0, "stats": stats, "coords": ncols,
                          "nonzero_true_minors": nonzero,
                          "plucker_matches_truth": mism == 0,
                          "mismatches": mism}, separators=(",", ":")))


if __name__ == "__main__":
    main()
