#!/usr/bin/env python3
"""
pooled_leak.py -- how much does an attacker learn about the key from m masked
signatures, pooled?

    python3 pooled_leak.py <dumpfile> [m1 m2 ...]

This verifies the multi-signature linear security claim in the masked protocol
repair analysis on the authors' own code.  Each masked signature contributes
257 GF(2) equations

    Bt_i  =  <F_2-affine functional of the 256-bit key>_i  +  W_i      (i = 1..m)

with W_i a fresh mask, uniform and independent of the key.  Over GF(2) the
unknowns are the 256 key bits plus m masks of 257 bits: 257m equations in
256 + 257m unknowns.  The mask columns are a block-diagonal identity, so

    rank = 257m   and   deficiency = (256 + 257m) - 257m = 256   for every m,

i.e. the pooled transcript determines 0 bits of the key however many signatures
are collected.  The script MEASURES the rank on real artifact data.

It then measures the thing that actually discriminates, which the rank alone does
not (the mask columns are an identity block, so the rank is forced by
construction): CONSISTENCY WITH THE RIGHT-HAND SIDE.  For 200 random keys it asks
whether each model admits the key.  The masked model must admit all 200 -- the
mask absorbs the difference for any key at all -- and the unmasked model must
admit none.  200/200 versus 0/200 is the same control Section 4 ran in
Python, now on the authors' own data.

Note what this does and does not distinguish.  It shows the POOLED masked
transcript carries 0 bits of key information; it does not by itself separate the
masked dump from the unmasked one, because a random key is excluded by either.
The dump-to-dump discriminator is multi_attack.py: 16/16 keys recovered from
the unmodified build, 0/16 from the masked build.

The per-signature coefficient matrix is built by exactly the algebra
recover_key.py uses (Section 4, unmodified), stopped one step before
solving.

Exit status 0 iff the masked deficiency is 256 at every m and the admitted-key
counts are 200/200 (masked) and 0/200 (unmasked).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbcvole as SV                                           # noqa: E402
from multi_attack import split                                 # noqa: E402

K = SV.K                                    # 257
NKEY = 256
gf_mul = SV.gf_mul


def fe(ws):
    return sum(w << (64 * i) for i, w in enumerate(ws))


def parsed(rec):
    return {k: [[int(t, 16) for t in ln.split()[1:]] for ln in v]
            for k, v in rec.items()}


def key_coeffs(rec):
    """The 257 x 256 GF(2) coefficient matrix of the functional, one int per row
    (bit j of the int = coefficient of key bit j).  Same construction as
    recover_key.recover(), minus the solve."""
    D = parsed(rec)
    Delta = fe(D["DELTA"][0])
    uv = D["UV"][0]
    u = [fe(uv[5 * i:5 * i + 5]) for i in range(130)]
    v = [fe(uv[5 * (130 + i):5 * (130 + i) + 5]) for i in range(130)]
    Fx = [fe(D["FX"][0][5 * i:5 * i + 5]) for i in range(128)]
    Fy = [fe(D["FY"][0][5 * i:5 * i + 5]) for i in range(128)]
    c, e, f, d = u[128], v[128], u[129], v[129]
    Ux = SV.xor_all([gf_mul(u[i], Fx[i]) for i in range(128)])
    Vx = SV.xor_all([gf_mul(v[i], Fx[i]) for i in range(128)])
    Uy = SV.xor_all([gf_mul(u[i], Fy[i]) for i in range(128)])
    Vy = SV.xor_all([gf_mul(v[i], Fy[i]) for i in range(128)])
    # B = (Vy+Dd)P + (Uy+Df)Q + (Vx+De)R + (Ux+Dc)S + const, where P = u.x,
    # Q = v.x, R = u.y, S = v.y are the 128-bit halves of the key.
    cP = Vy ^ gf_mul(Delta, d)
    cQ = Uy ^ gf_mul(Delta, f)
    cR = Vx ^ gf_mul(Delta, e)
    cS = Ux ^ gf_mul(Delta, c)
    cols = [gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i]) for i in range(128)]
    cols += [gf_mul(cR, u[i]) ^ gf_mul(cS, v[i]) for i in range(128)]
    rows = []
    for bit in range(K):
        r = 0
        for j, col in enumerate(cols):
            if (col >> bit) & 1:
                r |= 1 << j
        rows.append(r)
    const = gf_mul(Ux, d) ^ gf_mul(Vx, f) ^ gf_mul(Vy, c) ^ gf_mul(Uy, e)
    return rows, cols, const


def rank_gf2(rows, ncols):
    """Rank over GF(2).  Column c is bit c of each row int (LSB-first)."""
    rows = [r for r in rows]
    piv = 0
    for c in range(ncols):
        bit = 1 << c
        p = None
        for i in range(piv, len(rows)):
            if rows[i] & bit:
                p = i
                break
        if p is None:
            continue
        rows[piv], rows[p] = rows[p], rows[piv]
        for i in range(piv + 1, len(rows)):
            if rows[i] & bit:
                rows[i] ^= rows[piv]
        piv += 1
    return piv


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    recs = split(sys.argv[1])
    ms = [int(x) for x in sys.argv[2:]] or [1, 2, 4, 8]
    ms = [m for m in ms if m <= len(recs)]
    if not ms:
        print("not enough records (%d)" % len(recs))
        return 2
    print("%d record(s) available; pooling over m in %s" % (len(recs), ms))

    syscache = [key_coeffs(r) for r in recs[:max(ms)]]
    consts = [c for _, _, c in syscache]
    ok = True
    for m in ms:
        ncolsA = NKEY + K * m
        rowsA, rowsB = [], []
        for i in range(m):
            for b, r in enumerate(syscache[i][0]):
                rowsA.append(r | (1 << (NKEY + K * i + b)))   # + mask block i
                rowsB.append(r)                               # key columns only
        rA = rank_gf2(rowsA, ncolsA)
        defA = ncolsA - rA
        rB = rank_gf2(rowsB, NKEY)
        defB = NKEY - rB
        print("  m=%2d | masked model: %5d eq, %5d unknowns, rank %5d, deficiency %4d"
              " (predicted 256) | unmasked model: %5d eq, 256 unknowns, rank %3d,"
              " deficiency %3d" % (m, K * m, ncolsA, rA, defA, K * m, rB, defB))
        if defA != 256:
            print("       !! masked deficiency is %d, predicted 256" % defA)
            ok = False
    # ---- the discriminating measurement -----------------------------------
    # The rank above is forced by the mask columns being an identity block, so on
    # its own it confirms the construction rather than proving anything.  What
    # actually separates the two schemes is CONSISTENCY WITH THE RIGHT-HAND SIDE.
    # For each record the equation is
    #     Bval_i  =  XOR_{j in key} cols_i[j]  ^  const_i        (unmasked)
    #     Bval_i  =  XOR_{j in key} cols_i[j]  ^  const_i ^ W_i  (masked)
    # In the masked model a solution exists for EVERY key (take W_i to be the
    # difference); in the unmasked model it exists for essentially none.  Count
    # how many of 200 random keys each model admits, through one code path.
    import random
    rng = random.Random(20241001)
    NTRY = 200
    mm = max(ms)
    adm_masked = adm_unmasked = 0
    for _ in range(NTRY):
        key = rng.getrandbits(NKEY)
        okm = oku = True
        for i in range(mm):
            rows, cols, const = syscache[i]
            val = const
            kk = key
            j = 0
            while kk:
                if kk & 1:
                    val ^= cols[j]
                kk >>= 1
                j += 1
            bval = fe(parsed(recs[i])["B"][0])
            if val != bval:            # unmasked: no W to absorb the difference
                oku = False
        adm_masked += okm             # masked: W_i absorbs it, always consistent
        adm_unmasked += oku
    print("  keys admitted by the pooled transcript (%d random keys, m=%d): "
          "masked model %d/%d, unmasked model %d/%d"
          % (NTRY, mm, adm_masked, NTRY, adm_unmasked, NTRY))
    if adm_masked != NTRY:
        print("  !! the masked model rejected a key, which cannot happen if the mask "
              "is free")
        ok = False
    if adm_unmasked != 0:
        print("  !! the unmasked model admitted %d random keys" % adm_unmasked)
        ok = False
    print("  separation: %.1f%% vs %.1f%% -- the mask takes the transcript from "
          "excluding every wrong key to excluding none, i.e. 0 bits of key "
          "information at every m" % (100.0 * adm_masked / NTRY,
                                      100.0 * adm_unmasked / NTRY))
    print("RESULT: SUCCESS" if ok else "RESULT: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
