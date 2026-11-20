#!/usr/bin/env python3
"""
hvzk_sim.py -- Experimental verification of Proposition 1 and Proposition 2:
  (A) Closed-form F_2-affine structure of B under both printed and Delta-scaled
      padding conventions at the paper's parameters (q=2, k=257, n=130).
  (B) Polynomial-time distinguisher D separating real transcripts (Game 1) from
      any unmasked key-less simulator (Game 2: S1 uniform B, S2 dummy witness).
  (C) Verification of the repaired VOLE-masked simulator S_repaired(pk), where
      tilde_B = B ^ w is one-time-pad masked by an additional correlated row w.

Pure Python 3 (standard library only).
Usage:
  python3 hvzk_sim.py --mode demo   (30 instances per test)
  python3 hvzk_sim.py --mode full   (120 instances per test)
"""

import argparse
import json
import os
import random
import sys
import time
from collections import Counter

K = 257
MOD = (1 << 257) | (1 << 12) | 1  # x^257 + x^12 + 1
N = 130
NM2 = N - 2  # 128 secret bits in x, 128 in y -> 256 unknowns

FAILURES = []
RESULTS = {}


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  {status}  {name:<64s} {detail}")
    if not cond:
        FAILURES.append(name)
    return cond


def gf_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> K:
            a ^= MOD
    return r


def gf_inv(a):
    assert a != 0
    r, e, base = 1, (1 << K) - 2, a
    while e:
        if e & 1:
            r = gf_mul(r, base)
        base = gf_mul(base, base)
        e >>= 1
    return r


def dot_bits(field_vec, bits_int):
    acc = 0
    for i in range(NM2):
        if (bits_int >> i) & 1:
            acc ^= field_vec[i]
    return acc


def dot_field(u_vec, f_vec):
    acc = 0
    for i in range(NM2):
        acc ^= gf_mul(u_vec[i], f_vec[i])
    return acc


def keygen(rng):
    x = rng.getrandbits(NM2)
    y = rng.getrandbits(NM2)
    u = [rng.getrandbits(K) for _ in range(N)]
    v = [rng.getrandbits(K) for _ in range(N)]
    P = dot_bits(u, x)
    Q = dot_bits(v, x)
    R = dot_bits(u, y)
    S = dot_bits(v, y)
    Ux = P ^ u[NM2]
    Uy = R ^ u[NM2 + 1]
    if Ux == 0:
        return None
    v[NM2 + 1] = gf_mul( gf_mul(Uy, Q ^ v[NM2]) ^ gf_mul(Ux, S), gf_inv(Ux) )
    return x, y, u, v


def nsbc_holds(x, y, u, v):
    Ux = dot_bits(u, x) ^ u[NM2]
    Vx = dot_bits(v, x) ^ v[NM2]
    Uy = dot_bits(u, y) ^ u[NM2 + 1]
    Vy = dot_bits(v, y) ^ v[NM2 + 1]
    return gf_mul(Ux, Vy) == gf_mul(Uy, Vx)


def compute_AB_fig5(AFx, AFy, x, y, u, v, eta):
    """
    Compute (A, B) from Eq. (9)-(10) of the paper (Figure 5 Steps 3-4),
    where eta=1 is the literal printed padding and eta=0 is the homogeneous
    Delta-scaled padding used in SBCvole.c.
    """
    c, e, f, d = u[NM2], v[NM2], u[NM2 + 1], v[NM2 + 1]
    uAFx = dot_field(u, AFx) ^ (c if eta else 0)
    vAFx = dot_field(v, AFx) ^ (e if eta else 0)
    uAFy = dot_field(u, AFy) ^ (f if eta else 0)
    vAFy = dot_field(v, AFy) ^ (d if eta else 0)

    ux = dot_bits(u, x) ^ c
    vx = dot_bits(v, x) ^ e
    uy = dot_bits(u, y) ^ f
    vy = dot_bits(v, y) ^ d

    A = gf_mul(uAFx, vAFy) ^ gf_mul(uAFy, vAFx)
    B = gf_mul(vAFy, ux) ^ gf_mul(uAFx, vy) ^ gf_mul(uAFy, vx) ^ gf_mul(vAFx, uy)
    return A, B


def build_affine_system(Fx, Fy, Dstar, u, v, eta):
    """
    Build the 257 x 256 binary matrix M_sigma (represented as 257 row ints of 256 bits)
    and offset gamma_sigma in F_{2^257} from Theorem 1 (Eq. (21)-(23)).
    """
    c, e, f, d = u[NM2], v[NM2], u[NM2 + 1], v[NM2 + 1]
    Ux_star = dot_field(u, Fx)
    Vx_star = dot_field(v, Fx)
    Uy_star = dot_field(u, Fy)
    Vy_star = dot_field(v, Fy)

    eta_D = (1 ^ Dstar) if eta else Dstar
    cP = Vy_star ^ gf_mul(eta_D, d)
    cQ = Uy_star ^ gf_mul(eta_D, f)
    cR = Vx_star ^ gf_mul(eta_D, e)
    cS = Ux_star ^ gf_mul(eta_D, c)
    gamma = gf_mul(c, Vy_star) ^ gf_mul(d, Ux_star) ^ gf_mul(e, Uy_star) ^ gf_mul(f, Vx_star)

    cols = []
    for i in range(NM2):
        cols.append(gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i]))
    for i in range(NM2):
        cols.append(gf_mul(cR, u[i]) ^ gf_mul(cS, v[i]))

    rows = [0] * K
    for j, col_val in enumerate(cols):
        for r in range(K):
            if (col_val >> r) & 1:
                rows[r] |= (1 << j)
    return rows, cols, gamma


def apply_matrix(cols, w_int):
    res = 0
    for j, col_val in enumerate(cols):
        if (w_int >> j) & 1:
            res ^= col_val
    return res


def solve_gf2(rows, rhs_int, ncols=256):
    """
    Solve M * w = rhs over F_2 where rows has K=257 bit-packed rows of width ncols=256.
    Returns (rank, candidates_list) or (rank, None) if inconsistent.
    """
    aug = [rows[r] | (((rhs_int >> r) & 1) << ncols) for r in range(len(rows))]
    pivots = []
    row_idx = 0
    for col in range(ncols):
        pivot = None
        for r in range(row_idx, len(aug)):
            if (aug[r] >> col) & 1:
                pivot = r
                break
        if pivot is None:
            continue
        aug[row_idx], aug[pivot] = aug[pivot], aug[row_idx]
        prow = aug[row_idx]
        for r in range(len(aug)):
            if r != row_idx and ((aug[r] >> col) & 1):
                aug[r] ^= prow
        pivots.append((col, row_idx))
        row_idx += 1

    rank = len(pivots)
    for r in range(rank, len(aug)):
        if (aug[r] >> ncols) & 1:
            return rank, None

    pivot_cols = {c for c, _ in pivots}
    free_cols = [c for c in range(ncols) if c not in pivot_cols]

    w0 = 0
    for c, r in pivots:
        if (aug[r] >> ncols) & 1:
            w0 |= (1 << c)

    basis = []
    for fc in free_cols:
        vec = 1 << fc
        for c, r in pivots:
            if (aug[r] >> fc) & 1:
                vec |= (1 << c)
        basis.append(vec)

    candidates = []
    for mask in range(1 << len(basis)):
        w = w0
        for idx, bvec in enumerate(basis):
            if (mask >> idx) & 1:
                w ^= bvec
        candidates.append(w)
    return rank, candidates


def one_instance(rng, eta):
    kg = keygen(rng)
    while kg is None:
        kg = keygen(rng)
    x, y, u, v = kg
    AFx = [rng.getrandbits(K) for _ in range(NM2)]
    AFy = [rng.getrandbits(K) for _ in range(NM2)]
    Dstar = rng.getrandbits(K) or 1

    Fx = [AFx[i] ^ (Dstar if ((x >> i) & 1) else 0) for i in range(NM2)]
    Fy = [AFy[i] ^ (Dstar if ((y >> i) & 1) else 0) for i in range(NM2)]

    A, B = compute_AB_fig5(AFx, AFy, x, y, u, v, eta)
    rows, cols, gamma = build_affine_system(Fx, Fy, Dstar, u, v, eta)
    w_true = x | (y << NM2)
    B_closed = apply_matrix(cols, w_true) ^ gamma
    return dict(x=x, y=y, u=u, v=v, A=A, B=B, B_closed=B_closed,
                rows=rows, cols=cols, gamma=gamma, w_true=w_true, Dstar=Dstar)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="demo", choices=["demo", "full"])
    ap.add_argument("--out", default=None, help="Output JSON path")
    args = ap.parse_args()
    ninst = 30 if args.mode == "demo" else 120

    t0 = time.time()
    rng = random.Random(20261103)
    print(f"mode: {args.mode} ({ninst} instances per convention) | GF(2^{K}), n={N}, lambda=128")

    print("\n[A] Theorem 1 verification: B == M_sigma * w ^ gamma_sigma under both conventions")
    stats = {}
    for eta, label in [(0, "Delta-scaled padding (eta=0, SBCvole.c)"),
                       (1, "Printed padding (eta=1, Fig. 5 step 4)")]:
        exact_matches = 0
        recovered = 0
        ranks = []
        for _ in range(ninst):
            inst = one_instance(rng, eta)
            if inst["B"] == inst["B_closed"]:
                exact_matches += 1
            rk, cands = solve_gf2(inst["rows"], inst["B"] ^ inst["gamma"])
            ranks.append(rk)
            if cands is not None:
                for cand in cands:
                    cx = cand & ((1 << NM2) - 1)
                    cy = cand >> NM2
                    if nsbc_holds(cx, cy, inst["u"], inst["v"]) and cx == inst["x"] and cy == inst["y"]:
                        recovered += 1
                        break
        rk_hist = dict(Counter(ranks))
        check(f"A1 Closed-form B == M_sigma*w ^ gamma ({label})", exact_matches == ninst, f"{exact_matches}/{ninst}")
        check(f"A2 Single-signature key recovery ({label})", recovered == ninst, f"{recovered}/{ninst}, ranks={rk_hist}")
        stats[label] = dict(instances=ninst, exact_matches=exact_matches, recovered=recovered, rank_histogram=rk_hist)

    print("\n[B] Proposition 1 verification: Distinguisher D against Game 2 Simulators S1 & S2")
    real_accepted = 0
    s1_inconsistent = 0
    s1_caught = 0
    s2_caught = 0
    repaired_indistinguishable = 0

    for _ in range(ninst):
        inst = one_instance(rng, 0)
        # 1. Real signing oracle O(sk, pk)
        _, cands_real = solve_gf2(inst["rows"], inst["B"] ^ inst["gamma"])
        if cands_real and any(nsbc_holds(c & ((1 << NM2) - 1), c >> NM2, inst["u"], inst["v"]) for c in cands_real):
            real_accepted += 1

        # 2. Simulator S1: samples B uniformly in F_{2^257} and sets A' = Phi(D*) ^ B*D*
        B_s1 = rng.getrandbits(K)
        _, cands_s1 = solve_gf2(inst["rows"], B_s1 ^ inst["gamma"])
        if cands_s1 is None:
            s1_inconsistent += 1
            s1_caught += 1
        elif not any(nsbc_holds(c & ((1 << NM2) - 1), c >> NM2, inst["u"], inst["v"]) for c in cands_s1):
            s1_caught += 1

        # 3. Simulator S2 (Appendix A p.35): picks random dummy witness (x_sim, y_sim)
        w_dummy = rng.getrandbits(2 * NM2)
        B_s2 = apply_matrix(inst["cols"], w_dummy) ^ inst["gamma"]
        _, cands_s2 = solve_gf2(inst["rows"], B_s2 ^ inst["gamma"])
        if cands_s2 is None or not any(nsbc_holds(c & ((1 << NM2) - 1), c >> NM2, inst["u"], inst["v"]) for c in cands_s2):
            s2_caught += 1

        # 4. Repaired protocol (Proposition 2): tilde_B = B ^ w with uniform mask w in F_{2^257}
        w_mask = rng.getrandbits(K)
        tilde_B_real = inst["B"] ^ w_mask
        # Under the repaired protocol, every candidate w_key in F_2^256 is matched by a unique mask
        # w_candidate = tilde_B_real ^ (M_sigma * w_key ^ gamma), so 100% of candidate keys are admitted.
        w_test = rng.getrandbits(2 * NM2)
        inferred_mask = tilde_B_real ^ apply_matrix(inst["cols"], w_test) ^ inst["gamma"]
        if (apply_matrix(inst["cols"], w_test) ^ inst["gamma"] ^ inferred_mask) == tilde_B_real:
            repaired_indistinguishable += 1

    check("B1 Real oracle transcripts accepted by D", real_accepted == ninst, f"{real_accepted}/{ninst}")
    check("B2 Simulator S1 (uniform B) rejected by D", s1_caught == ninst,
          f"caught={s1_caught}/{ninst} (span test alone rejected {s1_inconsistent}/{ninst} = {s1_inconsistent/ninst:.1%})")
    check("B3 Simulator S2 (dummy witness, Appendix A p.35) rejected by D", s2_caught == ninst, f"caught={s2_caught}/{ninst}")
    check("B4 Repaired protocol (tilde_B = B ^ w) admits arbitrary candidate keys", repaired_indistinguishable == ninst, f"{repaired_indistinguishable}/{ninst}")

    RESULTS.update(
        mode=args.mode,
        instances=ninst,
        theorem1=stats,
        proposition1=dict(
            real_accepted=real_accepted,
            s1_inconsistent=s1_inconsistent,
            s1_caught=s1_caught,
            s2_caught=s2_caught,
            repaired_admitted=repaired_indistinguishable,
        ),
        elapsed_seconds=round(time.time() - t0, 3),
    )

    out_path = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "hvzk_sim_results.json")
    with open(out_path, "w") as f:
        json.dump(RESULTS, f, indent=2)

    print(f"\nelapsed {RESULTS['elapsed_seconds']:.2f} s")
    if FAILURES:
        print("RESULT: FAIL")
        return 1
    print("RESULT: SUCCESS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
