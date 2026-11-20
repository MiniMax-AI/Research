#!/usr/bin/env python3
"""
structural.py -- is the B leak a property of the PAPER, or only of the authors' C?

This script implements Figure 5 (p.20) and Figure 10 (p.24) of ePrint 2024/1537
exactly as printed -- no reference to the artifact -- and establishes three
things.

(1) AFFINENESS.  With the paper's own printed padding

        Fig. 5 step 4:  B = (v.(A_Fy,0,1))(u.(x,1,0)) + (u.(A_Fx,1,0))(v.(y,0,1))
                         - (u.(A_Fy,0,1))(v.(x,1,0)) - (v.(A_Fx,1,0))(u.(y,0,1))

    B is an F_2-AFFINE function of the 256 secret bits of (x,y).  This is tested
    directly, not argued: for random keys z1,z2 the second difference

        B(z1) + B(z2) + B(z1+z2) + B(0)

    must vanish identically if B is affine.  (In F_{2^k} subtraction is XOR, so
    the two minus signs above are pluses; the four quadratic-in-(x,y) terms are
    Delta*S*P, Delta*P*S, Delta*R*Q and Delta*Q*R, and they cancel pairwise.)

(2) CLOSED FORM.  The affine functional is

        B = (Vy* + (1+D) v_129) P + (Uy* + (1+D) u_129) Q
          + (Vx* + (1+D) v_128) R + (Ux* + (1+D) u_128) S
          + (u_128 Vy* + v_129 Ux* + v_128 Uy* + u_129 Vx*)

    with P = u.x, Q = v.x, R = u.y, S = v.y taken over the first n-2 = 128
    coordinates, Ux* = u.F_x(D*), Vx* = v.F_x(D*), Uy* = u.F_y(D*),
    Vy* = v.F_y(D*) -- every one of them computed by the verifier itself
    (Fig. 10 steps 11-13).  Tested against the brute-force evaluation of Fig. 5.

(3) THE ATTACK.  One 257x256 GF(2) elimination on that single equation recovers
    (x,y).  Run for the paper's n = 130, k = 257.

Run:  python3 structural.py [ntrials]
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbcvole as SV                                    # noqa: E402

K = SV.K
N = SV.N_DEFAULT
gf_mul = SV.gf_mul
xor_all = SV.xor_all


def dot(vec, bits):
    """F_2 dot product: `bits` holds 0/1 coefficients."""
    s = 0
    for i, b in enumerate(bits):
        if b:
            s ^= vec[i]
    return s


def fdot(vec, elts):
    """F_{2^k} dot product: `elts` holds field elements (Fig. 5 steps 1-2)."""
    s = 0
    for i, a in enumerate(elts):
        if a:
            s ^= gf_mul(vec[i], a)
    return s


# ---------------------------------------------------------------- Fig. 5, printed
def fig5_B(u, v, Ax, Ay, x, y):
    """Fig. 5 step 4, from Huth and Joux (ePrint 2024/1537, Fig. 5, step 4).

    Ax, Ay are the n-2 = 128-vectors A_Fx, A_Fy of steps 1-2.
    x, y are the n-2 = 128-bit secret vectors; the paper pads them (x,1,0),(y,0,1)
    and pads the A-vectors (A_Fx,1,0),(A_Fy,0,1).
    """
    xa = x + [1, 0]          # (x,1,0)
    ya = y + [0, 1]          # (y,0,1)
    Axa = Ax + [1, 0]        # (A_Fx,1,0)
    Aya = Ay + [0, 1]        # (A_Fy,0,1)
    t1 = gf_mul(fdot(v, Aya), dot(u, xa))
    t2 = gf_mul(fdot(u, Axa), dot(v, ya))
    t3 = gf_mul(fdot(u, Aya), dot(v, xa))
    t4 = gf_mul(fdot(v, Axa), dot(u, ya))
    return t1 ^ t2 ^ t3 ^ t4          # char 2: the two '-' are '+'


def fig10_Phi(u, v, Fx, Fy):
    """Fig. 10 step 13, from Huth and Joux (ePrint 2024/1537, Fig. 10, step 13)."""
    return (gf_mul(fdot(u, Fx + [1, 0]), fdot(v, Fy + [0, 1]))
            ^ gf_mul(fdot(u, Fy + [0, 1]), fdot(v, Fx + [1, 0])))


def closed_form_B(u, v, Fx, Fy, Delta, x, y):
    """The affine functional claimed in the docstring."""
    Ux = xor_all([gf_mul(u[i], Fx[i]) for i in range(N - 2)])
    Vx = xor_all([gf_mul(v[i], Fx[i]) for i in range(N - 2)])
    Uy = xor_all([gf_mul(u[i], Fy[i]) for i in range(N - 2)])
    Vy = xor_all([gf_mul(v[i], Fy[i]) for i in range(N - 2)])
    c, e, f, d = u[128], v[128], u[129], v[129]
    one = 1                                   # the field element 1 of F_{2^k}
    Dp = Delta ^ one
    P = dot(u[:128], x)
    Q = dot(v[:128], x)
    R = dot(u[:128], y)
    S = dot(v[:128], y)
    return (gf_mul(Vy ^ gf_mul(Dp, d), P)
            ^ gf_mul(Uy ^ gf_mul(Dp, f), Q)
            ^ gf_mul(Vx ^ gf_mul(Dp, e), R)
            ^ gf_mul(Ux ^ gf_mul(Dp, c), S)
            ^ gf_mul(c, Vy) ^ gf_mul(d, Ux) ^ gf_mul(e, Uy) ^ gf_mul(f, Vx))


def instance(rng):
    u = [rng.getrandbits(K) for _ in range(N)]
    v = [rng.getrandbits(K) for _ in range(N)]
    Fx = [rng.getrandbits(K) for _ in range(N - 2)]
    Fy = [rng.getrandbits(K) for _ in range(N - 2)]
    Delta = rng.getrandbits(K)
    return u, v, Fx, Fy, Delta


def main():
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = random.Random(20241001)
    print("paper's own parameters: n = %d, k = %d, field F_{2^%d}" % (N, K, K))

    # ---- (1) affineness of the printed Fig. 5 step 4 -------------------------
    bad = 0
    for _ in range(ntr):
        u, v, Fx, Fy, Delta = instance(rng)
        # A_Fx = F_x(Delta*) + Delta*.x  (the VOLE relation of Sec 3.2 / Fig. 10
        # step 11: F_x(Delta*) = sum_m alpha_m x_m^{[[b_m]]} = A_Fx + Delta*.x)
        z1x = [rng.getrandbits(1) for _ in range(N - 2)]
        z1y = [rng.getrandbits(1) for _ in range(N - 2)]
        z2x = [rng.getrandbits(1) for _ in range(N - 2)]
        z2y = [rng.getrandbits(1) for _ in range(N - 2)]
        z0x = [0] * (N - 2)
        z0y = [0] * (N - 2)

        def Bof(zx, zy):
            Ax = [Fx[i] ^ gf_mul(Delta, zx[i]) if zx[i] else Fx[i]
                  for i in range(N - 2)]
            Ay = [Fy[i] ^ gf_mul(Delta, zy[i]) if zy[i] else Fy[i]
                  for i in range(N - 2)]
            return fig5_B(u, v, Ax, Ay, zx, zy)

        s3x = [a ^ b for a, b in zip(z1x, z2x)]
        s3y = [a ^ b for a, b in zip(z1y, z2y)]
        if Bof(z1x, z1y) ^ Bof(z2x, z2y) ^ Bof(s3x, s3y) ^ Bof(z0x, z0y) != 0:
            bad += 1
    print("[1] second difference B(z1)+B(z2)+B(z1+z2)+B(0) over %d random "
          "instances of the PRINTED Fig. 5 step 4: %d non-zero "
          "-> B is F_2-affine in (x,y): %s" % (ntr, bad, bad == 0))

    # ---- (2) the closed form reproduces the printed figure -------------------
    mism = 0
    for _ in range(ntr):
        u, v, Fx, Fy, Delta = instance(rng)
        x = [rng.getrandbits(1) for _ in range(N - 2)]
        y = [rng.getrandbits(1) for _ in range(N - 2)]
        Ax = [Fx[i] ^ gf_mul(Delta, x[i]) if x[i] else Fx[i] for i in range(N - 2)]
        Ay = [Fy[i] ^ gf_mul(Delta, y[i]) if y[i] else Fy[i] for i in range(N - 2)]
        if fig5_B(u, v, Ax, Ay, x, y) != closed_form_B(u, v, Fx, Fy, Delta, x, y):
            mism += 1
    print("[2] closed form vs printed Fig. 5 step 4 over %d instances: %d "
          "mismatches -> closed form correct: %s" % (ntr, mism, mism == 0))

    # ---- (3) key recovery from that single equation --------------------------
    ok = 0
    ranks = {}
    for _ in range(20):
        u, v, Fx, Fy, Delta = instance(rng)
        # plant an honest NSBC key so that a public test exists
        x = [rng.getrandbits(1) for _ in range(N - 2)]
        y = [rng.getrandbits(1) for _ in range(N - 2)]
        # force the NSBC relation (u.x')(v.y') = (u.y')(v.x') by solving for v_128
        # exactly as the artifact's CreateKeyFromSecret does (SBCvole.c:96-99)
        Ux = dot(u, x + [1, 0])
        Uy = dot(u, y + [0, 1])
        Vy = dot(v, y + [0, 1])
        Vx_part = dot(v[:128], x)          # sum_{i<128} v_i x_i
        v[128] = gf_mul(gf_mul(Ux, Vy) ^ gf_mul(Uy, Vx_part), SV.gf_inv(Uy))
        Vx = dot(v, x + [1, 0])
        assert gf_mul(Ux, Vy) == gf_mul(Uy, Vx)

        B = closed_form_B(u, v, Fx, Fy, Delta, x, y)
        UxE = xor_all([gf_mul(u[i], Fx[i]) for i in range(N - 2)])
        VxE = xor_all([gf_mul(v[i], Fx[i]) for i in range(N - 2)])
        UyE = xor_all([gf_mul(u[i], Fy[i]) for i in range(N - 2)])
        VyE = xor_all([gf_mul(v[i], Fy[i]) for i in range(N - 2)])
        c, e, f, d = u[128], v[128], u[129], v[129]
        Dp = Delta ^ 1
        cP = VyE ^ gf_mul(Dp, d)
        cQ = UyE ^ gf_mul(Dp, f)
        cR = VxE ^ gf_mul(Dp, e)
        cS = UxE ^ gf_mul(Dp, c)
        rhs = B ^ gf_mul(c, VyE) ^ gf_mul(d, UxE) ^ gf_mul(e, UyE) ^ gf_mul(f, VxE)
        cols = [gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i]) for i in range(N - 2)]
        cols += [gf_mul(cR, u[i]) ^ gf_mul(cS, v[i]) for i in range(N - 2)]
        ncol = len(cols)
        rows = []
        for bit in range(K):
            r = 0
            for j, col in enumerate(cols):
                if (col >> bit) & 1:
                    r |= 1 << j
            rows.append(r | (((rhs >> bit) & 1) << ncol))
        sols, rank = SV.solve_gf2_candidates(rows, ncol, max_cand=256)
        ranks[rank] = ranks.get(rank, 0) + 1
        true = sum(b << i for i, b in enumerate(x + y))
        found = None
        for sol in sols:
            xb = [(sol >> i) & 1 for i in range(128)] + [1, 0]
            yb = [(sol >> (128 + i)) & 1 for i in range(128)] + [0, 1]
            if gf_mul(dot(u, xb), dot(v, yb)) == gf_mul(dot(u, yb), dot(v, xb)):
                found = sol
                break
        ok += (found == true)
    print("[3] key recovery from the PRINTED Fig. 5 step 4, 20 honest keys, "
          "n=130 k=257: %d/20 recovered exactly" % ok)
    print("    rank histogram of the 257x256 system: %s" % dict(sorted(ranks.items())))

    good = (bad == 0) and (mism == 0) and (ok == 20)
    print("RESULT: %s" % ("SUCCESS" if good else "FAIL"))
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
