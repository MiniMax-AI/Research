#!/usr/bin/env python3
"""
recover_key.py -- secret-key recovery from ONE signature of the authors' own
reference implementation (IACR artifact asiacrypt/2025/a10, SBCvole.c).

    python3 recover_key.py <dumpfile>

The input transcript is extracted either from a raw binary signature and public
key via the zero-patch C driver (`forge.c` `extract`) or from a `-DDUMPATTACK`
build of the reference artifact.  Every field the attack consumes -- DELTA, B,
FX, FY, UV -- is public: it is what the verifier CheckVOLESignature computes
(SBCvole.c:713-793) or deserialises from the signature (SBCvole.c:1155-1156).
Optional SECX / UVXY lines (when present in `-DDUMPATTACK` dumps) are used only
to check the answer; the attack itself never reads them.

Prints, on success:
    RECOVERED_SECX <w0> <w1> <w2> <w3>      (hex, the artifact's SecX[4])
    RECOVERED_RANK <rank>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbcvole as SV                                    # noqa: E402

K = SV.K
gf_mul = SV.gf_mul


def load(path):
    d = {}
    for line in open(path):
        p = line.split()
        if p:
            d.setdefault(p[0], []).append([int(t, 16) for t in p[1:]])
    return d


def fe(ws):
    return sum(w << (64 * i) for i, w in enumerate(ws))


def dot(vec, bits):
    s = 0
    for i, b in enumerate(bits):
        if b:
            s ^= vec[i]
    return s


def recover(path, verbose=True, max_cand=4096):
    D = load(path)
    Delta = fe(D["DELTA"][0])
    Bval = fe(D["B"][0])
    uv = D["UV"][0]
    # NOTE: nothing below reads D["SECX"] or D["UVXY"].  The attack consumes only
    # DELTA, B, FX, FY, UV -- all of which CheckVOLESignature computes or reads
    # out of the signature (`forge.c` `extract` writes only these five public
    # lines, with no secret lines present).
    u = [fe(uv[5 * i:5 * i + 5]) for i in range(130)]
    v = [fe(uv[5 * (130 + i):5 * (130 + i) + 5]) for i in range(130)]
    Fx = [fe(D["FX"][0][5 * i:5 * i + 5]) for i in range(128)]
    Fy = [fe(D["FY"][0][5 * i:5 * i + 5]) for i in range(128)]

    c, e, f, d = u[128], v[128], u[129], v[129]
    Ux = SV.xor_all([gf_mul(u[i], Fx[i]) for i in range(128)])
    Vx = SV.xor_all([gf_mul(v[i], Fx[i]) for i in range(128)])
    Uy = SV.xor_all([gf_mul(u[i], Fy[i]) for i in range(128)])
    Vy = SV.xor_all([gf_mul(v[i], Fy[i]) for i in range(128)])

    # B = (Vy+D d) P + (Uy+D f) Q + (Vx+D e) R + (Ux+D c) S + (c Vy + d Ux + e Uy + f Vx)
    cP = Vy ^ gf_mul(Delta, d)
    cQ = Uy ^ gf_mul(Delta, f)
    cR = Vx ^ gf_mul(Delta, e)
    cS = Ux ^ gf_mul(Delta, c)
    rhs = Bval ^ gf_mul(Ux, d) ^ gf_mul(Vx, f) ^ gf_mul(Vy, c) ^ gf_mul(Uy, e)

    cols = [gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i]) for i in range(128)]
    cols += [gf_mul(cR, u[i]) ^ gf_mul(cS, v[i]) for i in range(128)]
    ncol = len(cols)
    rows = []
    for bit in range(K):
        r = 0
        for j, col in enumerate(cols):
            if (col >> bit) & 1:
                r |= 1 << j
        rows.append(r | (((rhs >> bit) & 1) << ncol))

    sols, rank = SV.solve_gf2_candidates(rows, ncol, max_cand=max_cand)
    if verbose:
        print("system: %d GF(2) equations x %d unknowns ; rank = %d ; %d candidate(s)"
              % (K, ncol, rank, len(sols)))
    key = None
    for sol in sols:
        xb = [(sol >> i) & 1 for i in range(128)] + [1, 0]
        yb = [(sol >> (128 + i)) & 1 for i in range(128)] + [0, 1]
        # public test: does this candidate solve NSBC for the PUBLIC key (u,v)?
        if gf_mul(dot(u, xb), dot(v, yb)) == gf_mul(dot(u, yb), dot(v, xb)):
            key = sol
            break
    if key is None:
        return None, rank, None
    words = [(key >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(4)]
    true = None
    if "SECX" in D:
        s = D["SECX"][0]
        xt = [(s[0] >> i) & 1 for i in range(64)] + [(s[1] >> i) & 1 for i in range(64)]
        yt = [(s[2] >> i) & 1 for i in range(64)] + [(s[3] >> i) & 1 for i in range(64)]
        true = sum(b << i for i, b in enumerate(xt + yt))
    return words, rank, true


def main():
    if len(sys.argv) < 2:
        print("usage: recover_key.py <dumpfile>")
        return 2
    words, rank, true = recover(sys.argv[1])
    if words is None:
        print("RECOVERY: FAIL (no candidate passed the public NSBC test)")
        return 1
    print("RECOVERED_RANK %d" % rank)
    print("RECOVERED_SECX %016x %016x %016x %016x" % tuple(words))
    if true is not None:
        got = sum(w << (64 * i) for i, w in enumerate(words))
        print("MATCHES_ARTIFACT_SECRET_KEY %s" % (got == true))
        if got != true:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
