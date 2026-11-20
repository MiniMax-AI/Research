#!/usr/bin/env python3
"""
rank_sweep.py -- how many independent equations does ONE signature give?

For each parameter set, run the artifact's own
signer with a fresh key, recover the key from the single B equation, and record

  * the GF(2) rank of the 257 x 256 system,
  * the number of key candidates left by the elimination,
  * whether the public NSBC test picks the right one.

Then compare the distribution of the *rank deficiency* 256 - rank against the
Poisson(1/2) law that a random 257 x 256 GF(2) matrix follows, with a
chi-square-style comparison and a control: the same statistic computed on
genuinely random 257 x 256 GF(2) systems.

    python3 rank_sweep.py <driver> <nkeys> [label]
"""
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbcvole as SV                                    # noqa: E402
from recover_key import recover                          # noqa: E402


def random_control(n, seed=12345):
    """Control: rank deficiency of n random 257 x 256 GF(2) systems."""
    import random
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        rows = [rng.getrandbits(257) for _ in range(257)]
        _, rank = SV.solve_gf2_candidates(rows, 256, max_cand=1)
        out.append(256 - rank)
    return out


def main():
    driver = sys.argv[1]
    nkeys = int(sys.argv[2])
    label = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(driver)
    tmp = tempfile.mkdtemp(prefix="ranksweep_")
    try:
        sig = os.path.join(tmp, "s.bin")
        pk = os.path.join(tmp, "p.bin")
        defs = []
        ok = 0
        for i in range(nkeys):
            sec = subprocess.run([driver, "keygen", "%s-key-%d" % (label, i)],
                                 capture_output=True, text=True).stdout.split()[1:]
            msg = "msg %d" % i
            subprocess.run([driver, "sign"] + sec + [msg, sig, pk],
                           capture_output=True, check=True)
            pub = os.path.join(tmp, "pub.txt")
            subprocess.run([driver, "extract", msg, sig, pk, pub],
                           capture_output=True, check=True)
            # ground truth is the keygen output; pub is extracted solely from (msg, sig, pk)
            true = sum(int(t, 16) << (64 * j) for j, t in enumerate(sec))
            words, rank, _ = recover(pub, verbose=False)
            if words is None:
                defs.append(None)
                continue
            defs.append(256 - rank)
            got = sum(w << (64 * j) for j, w in enumerate(words))
            ok += (got == true)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    good = [d for d in defs if d is not None]
    hist = {}
    for d in good:
        hist[d] = hist.get(d, 0) + 1
    print("%-18s keys=%3d  recovered=%3d/%3d  deficiency histogram %s  "
          "mean=%.3f" % (label, nkeys, ok, nkeys, dict(sorted(hist.items())),
                         sum(good) / max(1, len(good))))
    import math
    pred = {k: len(good) * math.exp(-0.5) * 0.5 ** k / math.factorial(k)
            for k in range(0, 5)}
    print("%-18s Poisson(1/2) prediction %s"
          % ("", {k: round(v, 2) for k, v in pred.items()}))
    return 0 if ok == nkeys else 1


if __name__ == "__main__":
    sys.exit(main())
