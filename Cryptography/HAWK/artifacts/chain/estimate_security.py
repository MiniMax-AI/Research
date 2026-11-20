#!/usr/bin/env python3
"""Concrete security and AGPS20 gate-complexity estimator for HAWK-256/512/1024.

Reproduces and verifies every theoretical and heuristic security estimate in the
paper (Sections 3.5, 4.1, 6.2, and 6.3 / Table 4):

  1. Lattice invariants & structural Hermite limit (Proposition 2, Theorem 3,
     Proposition 4): exact log2 det Gram(Lambda_tau), log2 det Gram(Lambda'),
     Hermite ratios, and required GSA root Hermite factor delta_0,req = 0.996624 < 1.
  2. Expected BKZ block sizes beta(Lambda_tau) = (92, 205, 432) and
     beta(Lambda') = (58, 156, 351) from the GSA unique-SVP projection equation
     (plus invariance of beta_eff under strict ceiling inequalities, and optional
     `--simulator` verification of the DDGR20 / HAWK BKZ simulator block sizes).
  3. Effective sieve dimensions beta_eff = 8 * floor((beta - d4f(beta)) / 8)
     following the HAWK specification (§5.3) and Straznickas-Weis (Appendix C),
     and classical list-decoding gate counts log2(n) + C_gate(beta_eff) in the
     AGPS20 model (`jschanck/eprint-2019-1161` @ `1ec49e5`), including direct
     spherical-wedge circuit evaluation at d = 40 (self-checked against d = 64)
     alongside two-point boundary affine continuation from d = 64, 72.
  4. Provable (n/2 + 1)-dimensional SVP per-call and total gate bounds, and
     sieve database / runtime scaling from d = 118 down to d = 80.

Usage:
    python3 chain/estimate_security.py               # Fast verification (< 0.5 s)
    python3 chain/estimate_security.py --simulator   # Also run DDGR20 BKZ simulator
"""
from __future__ import annotations

import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from common import ROOT  # auto-runs fetch_third_party.sh --ensure

AGPS20_CSV = os.path.join(ROOT, "agps20", "cost-estimate-list_decoding-classical.csv")

PARAM_SETS = (
    ("HAWK-256", 256, 64),
    ("HAWK-512", 512, 128),
    ("HAWK-1024", 1024, 256),
)


def load_agps20_table(path: str = AGPS20_CSV) -> dict[int, float]:
    """Load tabulated classical BDGL16 list-decoding log2(gates) from AGPS20."""
    table: dict[int, float] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            table[int(row["d"])] = float(row["log_cost"])
    assert 64 in table and 872 in table, "Incomplete AGPS20 classical cost table"
    return table


def agps20_list_decoding_exact(d: int, max_n: int = 8192) -> tuple[int, float]:
    """Evaluate the exact AGPS20 classical BDGL16 list-decoding gate count at sieve
    dimension `d` from first principles (Albrecht-Gheorghiu-Postlethwaite-Schanck,
    ASIACRYPT 2020; `jschanck/eprint-2019-1161` `cost.py:list_decoding`).

    Used to evaluate d = 40 (below the tabulated CSV range d >= 64) without
    ad-hoc extrapolation, and self-checked against the CSV at d = 64.
    """
    from scipy.integrate import quad
    from scipy.optimize import fminbound
    from scipy.special import betainc, gammaln

    word_size = 32
    k_div_n = 1.0 / 3.0
    sphere_ratio = math.exp(
        -0.5 * math.log(math.pi) + gammaln(0.5 * d) - gammaln(0.5 * (d - 1))
    )

    def wedge_scaled(theta1: float, x: float) -> float:
        """W(d, theta1, theta1, x) / sin(theta1)^(d - 2)."""
        if x <= 0.0 or x >= 2.0 * theta1:
            return 0.0
        sin_t1 = math.sin(theta1)
        tan_half_x = math.tan(0.5 * x)

        def integrand(u: float) -> float:
            ratio_t = tan_half_x / math.tan(u)
            if ratio_t >= 1.0:
                return 0.0
            return ((math.sin(u) / sin_t1) ** (d - 2)) * betainc(
                0.5 * (d - 2), 0.5, 1.0 - ratio_t * ratio_t
            )

        return quad(integrand, 0.5 * x, theta1, epsabs=0.0, epsrel=1e-9)[0] * sphere_ratio

    def eval_for_filter(n_code: int, theta1: float) -> tuple[float, bool]:
        k = int(k_div_n * n_code)
        upper = min(math.pi, 2.0 * theta1)
        x_peak = math.acos(math.cos(theta1) ** 2)
        norm_peak = max(
            wedge_scaled(theta1, x_peak) * (math.sin(x_peak) ** (d - 2)), 1e-300
        )

        def wa_norm(x: float) -> float:
            if x <= 0.0 or x >= upper:
                return 0.0
            return wedge_scaled(theta1, x) * ((math.sin(x) ** (d - 2)) / norm_peak)

        lo = max(0.0, x_peak - 8.0 / math.sqrt(d))
        hi = min(upper, x_peak + 8.0 / math.sqrt(d))
        den = quad(wa_norm, lo, hi, points=[x_peak], epsabs=0.0, epsrel=1e-9)[0]
        pf_val = (
            quad(
                lambda x: betainc(n_code - k, k + 1, 1.0 - x / math.pi) * wa_norm(x),
                lo,
                hi,
                points=[x_peak],
                epsabs=0.0,
                epsrel=1e-9,
            )[0]
            / den
        )

        hi_ngr = min(math.pi / 3.0, upper)
        x_ngr = min(math.pi / 3.0, x_peak)
        w_ngr = wedge_scaled(theta1, x_ngr)

        def wa_ngr(x: float) -> float:
            if x <= 0.0 or x >= hi_ngr:
                return 0.0
            return (wedge_scaled(theta1, x) / w_ngr) * (
                (math.sin(x) / math.sin(x_ngr)) ** (d - 2)
            )

        lo_ngr = max(0.0, x_ngr - 8.0 / math.sqrt(d))
        ngr_den = quad(wa_ngr, lo_ngr, hi_ngr, epsabs=0.0, epsrel=1e-9)[0]
        ngr_num = quad(
            lambda x: betainc(n_code - k, k + 1, 1.0 - x / math.pi) * wa_ngr(x),
            lo_ngr,
            hi_ngr,
            epsabs=0.0,
            epsrel=1e-9,
        )[0]
        eta = 1.0 - ngr_num / ngr_den

        log2_c_pi3 = math.log2(0.5 * betainc(0.5 * (d - 1), 0.5, 0.75))
        log2_n_db = 1.0 - math.log2(1.0 - eta) - log2_c_pi3
        log2_w0 = math.log2(wedge_scaled(theta1, math.pi / 3.0)) + (d - 2) * math.log2(
            math.sin(theta1)
        )
        log2_filters = -log2_w0

        m = math.log2(d)
        ip_cost = (d / m) * (word_size ** 2)
        log2_z = log2_filters / m
        z_sub = 2.0 ** log2_z
        cost_tree_iter = 5.0 * word_size + log2_z
        cost_compare_swap = word_size + log2_z
        preprocess_cost = m * z_sub * ip_cost + m * z_sub * log2_z * cost_compare_swap

        log2_c_t1 = math.log2(
            0.5 * betainc(0.5 * (d - 1), 0.5, math.sin(theta1) ** 2)
        )
        insert_cost = preprocess_cost + (2.0 ** (log2_filters + log2_c_t1)) * cost_tree_iter
        query_cost = insert_cost
        log2_bucket_size = log2_filters + log2_n_db + 2.0 * log2_c_t1

        ell = math.ceil(math.log2(n_code))
        look_gates = n_code + (n_code - ell - 1) * 5 + ell
        log2_total = log2_n_db + math.log2(
            insert_cost + query_cost + look_gates * (2.0 ** log2_bucket_size)
        )
        popcounts_dominate = (1.0 / pf_val) > ((word_size ** 2) * d / float(n_code))
        return log2_total, popcounts_dominate

    n_pow2 = 1
    while n_pow2 < d:
        n_pow2 *= 2
    while True:
        n_code = n_pow2 - 1
        theta_opt = float(
            fminbound(
                lambda t: eval_for_filter(n_code, float(t))[0],
                math.pi / 6.0 + 1e-2,
                math.pi / 2.0 - 1e-2,
                xtol=1e-6,
            )
        )
        log_cost, dom = eval_for_filter(n_code, theta_opt)
        if dom or n_pow2 >= max_n:
            return n_code, log_cost
        n_pow2 *= 2


def rhf_cn11(beta: float) -> float:
    """Chen-Nguyen root Hermite factor delta_beta (Paper Section 2.4)."""
    return (((math.pi * beta) ** (1.0 / beta)) * beta / (2.0 * math.pi * math.e)) ** (
        1.0 / (2.0 * (beta - 1.0))
    )


def solve_gsa_beta(n: int, log_hermite_ratio: float, exp_offset: int = 0) -> float:
    """Solve sqrt(beta / n) * (lambda_1 / vol^(1/n)) = delta_beta^(2*beta - n - exp_offset)."""
    def f(b: float) -> float:
        return (
            0.5 * math.log(b / n)
            + log_hermite_ratio
            - (2.0 * b - n - exp_offset) * math.log(rhf_cn11(b))
        )

    lo, hi = 36.01, float(n)
    flo = f(lo)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) * flo <= 0.0:
            hi = mid
        else:
            lo, flo = mid, f(mid)
    return 0.5 * (lo + hi)


def d4f_hawk(beta: int) -> int:
    """Dimensions-for-free d4f(beta) rounded as in [HAWK25, §5.3] and [SW26, App. C]."""
    return round(beta * math.log(4.0 / 3.0) / math.log(beta / (2.0 * math.pi * math.e)))


def beta_eff_hawk(beta: int) -> tuple[int, int]:
    """Return (d4f, beta_eff) with beta_eff lowered to a multiple of 8 ([HAWK25, §5.3])."""
    d4f = d4f_hawk(beta)
    beff = ((beta - d4f) // 8) * 8
    return d4f, beff


def provable_call_upper_bound(n: int) -> float:
    """Upper bound on exact-SVP(n/2+1) oracle calls in Theorem 6.1 ([SW26, §6.1, p.12])."""
    p = n // 2 + 1
    ln_max_b_s2 = 0.5 * n * math.log(2.0)
    it_s2 = (p ** 2) * math.log(p) + p * math.log(2.0 * (p - 1) * ln_max_b_s2) + 3.0 * p
    d_odd = n + 1
    kp = n // 2
    ln_max_b_s3 = 0.5 * (n - 1) * math.log(2.0)
    it_s3 = (
        d_odd * (1.0 + (kp + 1) * math.log(kp + 1))
        + (kp + 1) * math.log(2.0 * kp * ln_max_b_s3)
        + 3.0 * (kp + 1)
    )
    return (2.0 * it_s2 + n) + 2.0 * (2.0 * it_s3 + n)


def run_ddgr20_bkz_simulator() -> dict[int, tuple[float, float]]:
    """Run the Chen-Nguyen / DDGR20 BKZ simulator (as in `hawk-sign/aux` `BKZ_simulator.sage`)
    on both Lambda_tau (n/2 + 1 minimal pairs) and Lambda' (unique minimal pair).
    """
    from scipy.special import gammainc

    rk = [
        0.789527997160000, 0.780003183804613, 0.750872218594458, 0.706520454592593,
        0.696345241018901, 0.660533841808400, 0.626274718790505, 0.581480717333169,
        0.553171463433503, 0.520811087419712, 0.487994338534253, 0.459541470573431,
        0.414638319529319, 0.392811729940846, 0.339090376264829, 0.306561491936042,
        0.276041187709516, 0.236698863270441, 0.196186341673080, 0.161214212092249,
        0.110895134828114, 0.0678261623920553, 0.0272807162335610, -0.0234609979600137,
        -0.0320527224746912, -0.0940331032784437, -0.129109087817554, -0.176965384290173,
        -0.209405754915959, -0.265867993276493, -0.299031324494802, -0.349338597048432,
        -0.380428160303508, -0.427399405474537, -0.474944677694975, -0.530140672818150,
        -0.561625221138784, -0.612008793872032, -0.669011014635905, -0.713766731570930,
        -0.754041787011810, -0.808609696192079, -0.859933249032210, -0.884479963601658,
        -0.886666930030433,
    ]
    sim_c = [0.0] + [rk[-i] - sum(rk[-i:]) / i for i in range(1, 46)]
    for d in range(46, 1050):
        gh2 = ((math.pi * d) ** (1.0 / d)) * d / (2.0 * math.pi * math.e)
        sim_c.append((0.5 * math.log(gh2) - math.log(0.5) / d) / math.log(2.0))

    def predict_beta(d: int, logvol: float, num_targets: int) -> float:
        delta2 = 1.021900
        l = [
            ((d - 1 - 2 * i) * math.log(delta2) + logvol / d) / math.log(2.0)
            for i in range(d)
        ]
        rem_pr = 1.0
        avg_beta = 0.0
        for beta in range(2, d):
            for k in range(d - 1):
                blk = min(beta, d - k)
                f_idx = k + blk
                lma = sum(l[k:f_idx]) / blk + sim_c[blk]
                if lma < l[k]:
                    diff = l[k] - lma
                    l[k] -= diff
                    step = diff / (f_idx - k - 1)
                    for a in range(k + 1, f_idx):
                        l[a] += step
            proba = float(gammainc(0.5 * beta, 0.5 * (2.0 ** (2.0 * l[d - beta]))))
            for j in range(2, int(d / beta + 1)):
                idx = d - j * (beta - 1) - 1
                xt = 2.0 ** (2.0 * l[idx])
                proba = min(proba, float(gammainc(0.5 * (d - idx), 0.5 * xt)))
            for _ in range(num_targets):
                avg_beta += beta * rem_beta_step(rem_pr, proba)
                rem_pr *= 1.0 - proba
            if rem_pr < 0.001:
                avg_beta += beta * rem_pr
                break
        return avg_beta

    def rem_beta_step(rem: float, pr: float) -> float:
        return rem * pr

    out: dict[int, tuple[float, float]] = {}
    for _name, n, _claimed in PARAM_SETS:
        lv_lt = 0.5 * n * math.log(n) + 0.25 * (n - 2) * math.log(2.0)
        lv_lp = 0.5 * n * math.log(n) + 0.25 * (5 * n - 6) * math.log(2.0)
        out[n] = (predict_beta(n, lv_lt, n // 2 + 1), predict_beta(n, lv_lp, 1))
    return out


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    run_sim = "--simulator" in argv

    agps_table = load_agps20_table()

    # Self-check exact AGPS20 spherical-wedge evaluator against tabulated CSV at d = 64,
    # then evaluate d = 40 directly.
    n64, c64_exact = agps20_list_decoding_exact(64)
    assert n64 == 2047 and abs(c64_exact - agps_table[64]) < 1e-4, (
        f"AGPS20 d=64 self-check mismatch: got {c64_exact}, expected {agps_table[64]}"
    )
    n40, c40_exact = agps20_list_decoding_exact(40)
    slope_64_72 = (agps_table[72] - agps_table[64]) / 8.0
    c40_affine = agps_table[64] + slope_64_72 * (40 - 64)

    # 1. Structural Hermite limit & determinant certificates at n = 256 (Sections 3.4, 4.1)
    n256 = 256
    logdet_lt_256 = n256 * math.log2(2 * n256) + (n256 // 2 - 1)
    logdet_lp_256 = logdet_lt_256 + 2 * (n256 - 1)
    hermite_lt_256 = 2.0 ** (-(n256 - 2) / (4.0 * n256))
    hermite_lp_256 = 2.0 ** (-(5 * n256 - 6) / (4.0 * n256))
    rhf_req_256 = hermite_lp_256 ** (1.0 / (n256 - 1))
    assert logdet_lt_256 == 2431.0 and logdet_lp_256 == 2941.0
    assert abs(hermite_lp_256 - 0.4222) < 1e-4
    assert abs(rhf_req_256 - 0.996624) < 1e-6 and rhf_req_256 < 1.0

    print("  [1] Lattice Determinant & Structural Hermite Limit Checks (n = 256):")
    print(f"      log2 det Gram(Lambda_tau) = {logdet_lt_256:.2f}  |  log2 det Gram(Lambda') = {logdet_lp_256:.2f}")
    print(
        f"      Hermite ratio: Lambda_tau = {hermite_lt_256:.4f} (-> 2^(-1/4) = {2**-0.25:.4f}), "
        f"Lambda' = {hermite_lp_256:.4f} (-> 2^(-5/4) = {2**-1.25:.4f})"
    )
    print(f"      Required GSA root Hermite factor at n = 256: delta_0,req = {rhf_req_256:.6f} < 1")

    # 2. GSA block sizes, effective sieve dimensions, and AGPS20 classical gate counts (Table 4)
    print("\n  [2] Heuristic GSA Block Sizes & AGPS20 Gate Complexity (Paper Section 6.3, Table 4):")
    print(
        "      Variant     n    Claim  beta(Lt)->beff  beta(Lp)->beff  Gates(Lt)  Gates(Lp)   DeltaSec"
    )

    expected_paper = {
        256: (92, 72, 58, 40, 50.83, 41.56, -9.28),
        512: (205, 176, 156, 136, 80.79, 69.64, -11.15),
        1024: (432, 392, 351, 312, 146.48, 122.73, -23.74),
    }

    for name, n, claimed in PARAM_SETS:
        log_h_lt = -(n - 2) / (4.0 * n) * math.log(2.0)
        log_h_lp = -(5 * n - 6) / (4.0 * n) * math.log(2.0)

        beta_lt_cont = solve_gsa_beta(n, log_h_lt, exp_offset=0)
        beta_lp_cont = solve_gsa_beta(n, log_h_lp, exp_offset=0)
        beta_lt = round(beta_lt_cont)
        beta_lp = round(beta_lp_cont)

        # Also check ceiling under 2*beta - n - 1: beta_eff must remain invariant
        _, beff_lt_ceil = beta_eff_hawk(math.ceil(solve_gsa_beta(n, log_h_lt, exp_offset=1)))
        _, beff_lp_ceil = beta_eff_hawk(math.ceil(solve_gsa_beta(n, log_h_lp, exp_offset=1)))

        _, beff_lt = beta_eff_hawk(beta_lt)
        _, beff_lp = beta_eff_hawk(beta_lp)
        assert beff_lt == beff_lt_ceil and beff_lp == beff_lp_ceil

        cgate_lt = agps_table[beff_lt]
        cgate_lp = agps_table[beff_lp] if beff_lp in agps_table else c40_exact

        gates_lt = math.log2(n) + cgate_lt
        gates_lp = math.log2(n) + cgate_lp
        delta_sec = gates_lp - gates_lt

        exp_row = expected_paper[n]
        assert (beta_lt, beff_lt, beta_lp, beff_lp) == exp_row[:4], (
            f"Blocksize mismatch at n={n}: got {(beta_lt, beff_lt, beta_lp, beff_lp)}"
        )
        assert abs(round(gates_lt, 2) - exp_row[4]) < 1e-9
        assert abs(round(gates_lp, 2) - exp_row[5]) < 1e-9
        assert abs(round(delta_sec, 2) - exp_row[6]) < 1e-9

        print(
            f"      {name:9s} {n:4d}  {claimed:3d}b    "
            f"{beta_lt:3d} -> {beff_lt:3d}      {beta_lp:3d} -> {beff_lp:3d}     "
            f"2^{gates_lt:6.2f}   2^{gates_lp:6.2f}   {delta_sec:+6.2f} bits"
        )

    print(
        f"      [d=40 cross-checks] Exact AGPS20 (n_code={n40}): 2^{8 + c40_exact:.2f} "
        f"| d=64..72 affine continuation: 2^{8 + c40_affine:.2f} "
        f"| d=64 clamped: 2^{8 + agps_table[64]:.2f}"
    )

    # 3. Provable (n/2 + 1)-dimensional SVP bounds & d = 118 -> 80 sieve scaling
    prov_summary = []
    for name, n, _ in PARAM_SETS:
        p = n // 2 + 1
        _, beff_prov = beta_eff_hawk(p)
        c_call = agps_table[beff_prov]
        c_tot = c_call + math.log2(provable_call_upper_bound(n))
        prov_summary.append(f"{name}: 2^{round(c_call)}/call (2^{math.ceil(c_tot)} total)")

    db_ratio = 2.0 ** (0.2075 * (118 - 80))
    time_ratio = 2.0 ** (0.2345 * (118 - 80))
    assert round(db_ratio) == 236 and round(time_ratio) == 481

    print("\n  [3] Provable SVP(n/2+1) Bounds & Sieve Scaling (Sections 6.2-6.3):")
    print("      " + "  |  ".join(prov_summary))
    print(f"      Sieve drop d = 118 -> 80: database ratio = {db_ratio:.1f}x, predicted sieve work = {time_ratio:.1f}x")

    if run_sim:
        sim = run_ddgr20_bkz_simulator()
        print("\n  [4] DDGR20 / HAWK BKZ Simulator Block Sizes:")
        for name, n, _ in PARAM_SETS:
            b_lt_s, b_lp_s = sim[n]
            print(f"      {name:9s}: beta_sim(Lambda_tau) = {b_lt_s:.2f} -> beta_sim(Lambda') = {b_lp_s:.2f}")
        assert abs(sim[256][0] - 82.58) < 0.02 and abs(sim[256][1] - 60.19) < 0.02
        assert abs(sim[512][0] - 197.68) < 0.02 and abs(sim[512][1] - 158.18) < 0.02
        assert abs(sim[1024][0] - 429.16) < 0.02 and abs(sim[1024][1] - 358.29) < 0.02

    print("  PASS  All GSA block sizes, AGPS20 gate complexities, and scaling bounds verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
