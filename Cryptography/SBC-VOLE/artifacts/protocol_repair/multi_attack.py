#!/usr/bin/env python3
"""
multi_attack.py -- run the Section 4 single-signature key-recovery attack
against every record in a DUMPATTACK dump, rather than just the first.

    python3 multi_attack.py <dumpfile>

The artifact's benchmark driver appends one record per signing iteration (and
reseeds itself between iterations via `randomness[0]++`, SBCvole.c:1478), so a
run with -DNBREPEATS=N yields N independent signatures.  This script splits them
and calls recover_key.recover() -- unmodified -- on each.

Prints one line per record and a tally.  Exit status 0 iff the outcome is
uniform: every record recovered, or none recovered.  A mixed tally is the
interesting failure and is reported as such.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recover_key as RK                                       # noqa: E402

KEYS = ["SECX", "UVXY", "DELTA", "B", "FX", "FY", "UV"]


def split(path):
    """The artifact dumps SECX/UVXY from its key-generation loop and
    DELTA/B/FX/FY/UV from its verification loop, so the two groups are separate.
    Split the verifier group on DELTA (one per signing iteration) and attach the
    key-generation lines to each record for the correctness check only -- the
    attack itself never reads them."""
    keylines, recs, cur = [], [], None
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] in ("SECX", "UVXY"):
            keylines.append(line)
        elif p[0] == "DELTA":
            if cur is not None:
                recs.append(cur)
            cur = {}
            cur["DELTA"] = [line]
        elif cur is not None:
            cur.setdefault(p[0], []).append(line)
    if cur is not None:
        recs.append(cur)
    for r in recs:
        r["SECX"] = keylines[:2]
    return recs


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    recs = split(sys.argv[1])
    if not recs:
        print("no records in dump")
        return 2
    # A record must carry exactly one of each public field; bail loudly if not.
    bad = [i for i, r in enumerate(recs)
           if any(len(r.get(k, [])) != 1 for k in ("DELTA", "B", "FX", "FY", "UV"))]
    if bad:
        print("ragged records:", bad[:5])
        return 2

    recovered = 0
    matched = 0
    ranks = []
    for i, r in enumerate(recs):
        fd, tmp = tempfile.mkstemp(suffix=".txt")
        with os.fdopen(fd, "w") as fh:
            for k in KEYS:
                for ln in r.get(k, []):
                    fh.write(ln)
        words, rank, true = RK.recover(tmp, verbose=False)
        os.unlink(tmp)
        ranks.append(rank)
        if words is None:
            print("  record %2d: rank %3d -> NO CANDIDATE (system inconsistent)" % (i, rank))
        else:
            recovered += 1
            got = sum(w << (64 * j) for j, w in enumerate(words))
            ok = (true is not None and got == true)
            matched += ok
            print("  record %2d: rank %3d -> KEY RECOVERED, matches the artifact's SecX: %s"
                  % (i, rank, ok))

    n = len(recs)
    print("%d record(s): recovered %d, of which correct %d; rank range %d..%d"
          % (n, recovered, matched, min(ranks), max(ranks)))
    if recovered == 0:
        print("VERDICT: the attack fails on every record")
        return 0
    if matched == n:
        print("VERDICT: the attack recovers the key on every record")
        return 0
    print("VERDICT: MIXED -- this is a failure of the uniformity claim")
    return 1


if __name__ == "__main__":
    sys.exit(main())
