# -*- coding: utf-8 -*-
"""Verify the closed-form F_2-affine expression for B (Section 4, Lemma 1).

Usage:  python3 verify_lemma1.py

Evaluates Lemma 1 term by term and checks exact agreement with both
structural.py's closed_form_B and Figure 5 Step 4 across 40 random instances
at the specification parameters (n = 130, k = 257).
"""
import os
import random
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sbcvole as SV
from structural import fig5_B, closed_form_B, dot, fdot, N, K
gf_mul = SV.gf_mul
xor_all = SV.xor_all

def lemma1_closed_form(u, v, Fx, Fy, D, x, y):
    """Exactly as printed in Paper Section 4 (Lemma 1)."""
    Ux = xor_all([gf_mul(u[i], Fx[i]) for i in range(N-2)])   # Ux* = u.F_x(D*)
    Vx = xor_all([gf_mul(v[i], Fx[i]) for i in range(N-2)])   # Vx*
    Uy = xor_all([gf_mul(u[i], Fy[i]) for i in range(N-2)])   # Uy*
    Vy = xor_all([gf_mul(v[i], Fy[i]) for i in range(N-2)])   # Vy*
    c = u[N-2]; e = v[N-2]; f = u[N-1]; d = v[N-1]            # c,e,f,g (d := g) in Paper Section 4 (Lemma 1)
    P = dot(u[:N-2], x); Q = dot(v[:N-2], x)
    R = dot(u[:N-2], y); S = dot(v[:N-2], y)
    one = 1
    return (gf_mul(Vy ^ gf_mul(D ^ one, d), P)
          ^ gf_mul(Uy ^ gf_mul(D ^ one, f), Q)
          ^ gf_mul(Vx ^ gf_mul(D ^ one, e), R)
          ^ gf_mul(Ux ^ gf_mul(D ^ one, c), S)
          ^ gf_mul(c, Vy) ^ gf_mul(d, Ux) ^ gf_mul(e, Uy) ^ gf_mul(f, Vx))

rng = random.Random(777)
bad = 0
for t in range(40):
    u  = [rng.getrandbits(K) for _ in range(N)]
    v  = [rng.getrandbits(K) for _ in range(N)]
    Fx = [rng.getrandbits(K) for _ in range(N-2)]
    Fy = [rng.getrandbits(K) for _ in range(N-2)]
    D  = rng.getrandbits(K)
    x  = [rng.getrandbits(1) for _ in range(N-2)]
    y  = [rng.getrandbits(1) for _ in range(N-2)]
    # printed Fig. 5, with the coordinate-wise VOLE relation A_i = F_i + D*x_i
    Ax = [Fx[i] ^ gf_mul(D, x[i]) if x[i] else Fx[i] for i in range(N-2)]
    Ay = [Fy[i] ^ gf_mul(D, y[i]) if y[i] else Fy[i] for i in range(N-2)]
    printed        = fig5_B(u, v, Ax, Ay, x, y)
    lemma1_val     = lemma1_closed_form(u, v, Fx, Fy, D, x, y)
    structural_val = closed_form_B(u, v, Fx, Fy, D, x, y)
    if not (printed == lemma1_val == structural_val):
        bad += 1
        print("MISMATCH at trial", t)
print("trials=40  mismatches=%d  -> Paper Section 4 (Lemma 1) closed form is %s"
      % (bad, "CORRECT" if bad == 0 else "WRONG"))
sys.exit(1 if bad else 0)
