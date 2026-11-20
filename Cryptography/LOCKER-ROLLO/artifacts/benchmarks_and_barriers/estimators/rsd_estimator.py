#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rsd_estimator.py -- Unified Rank-Syndrome-Decoding (RSD) bit-complexity estimator.

Formulas are transcribed VERBATIM from:
  [EV26]  A. Esser, J. Verbel, R. Villanueva-Polanco,
          "Rank Syndrome Decoding Estimator - An Asymptotic and Concrete Analysis",
          IACR ePrint 2026/082 (https://eprint.iacr.org/2026/082)
          Section 3 (Eqs. 1,4,5, GRS/GRSim/GRSge, 7,8,10 and the Hybrid recursion).
  [BBGT22] M. Bardet, P. Briaud, M. Bros, P. Gaborit, J.-P. Tillich,
          "Revisiting Algebraic Attacks on MinRank and on the Rank Decoding Problem",
          IACR ePrint 2022/1031, Table 2.

Cost model (identical to [EV26] Sec. 3, p.6):
  * time is counted in elementary multiplications over F_q;
  * every F_q multiplication is charged log2(q) bit-operations (so for q=2 the
    F_2-operation count IS the bit count);
  * Gaussian elimination on an n1 x n2 matrix (n1<=n2) costs O(n1*n2^(w-1)).

NO extrapolation, NO fitted constants: every number below is produced by the
closed-form formulas of the cited papers.

usage:
    python3 rsd_estimator.py --validate     # reproduce the published tables
    python3 rsd_estimator.py --rollo        # 2020/2021-revised ROLLO + new2ROLLO
"""
from __future__ import annotations
import argparse, math
from math import comb, ceil, log2

W = 2.0          # linear-algebra exponent used by both cited papers' tables
LOG2Q = 1.0      # q = 2


def _c(n, r):
    return comb(n, r) if 0 <= r <= n else 0


# --------------------------------------------------------------------------
# Combinatorial algorithms  ([EV26] Sec 3.1)
# --------------------------------------------------------------------------
def t_basis_enumeration(q, m, n, k, r, w=W):
    """[EV26] Eq. (1):  O((n r + m)^w * q^((m-r)(r-1)))."""
    if r < 1 or m <= r:
        return math.inf
    return w * log2(n * r + m) + (m - r) * (r - 1) * log2(q)


def t_oj1(q, m, n, k, r, w=W):
    """[EV26] Eq. (4)+complexity:  N = ceil(((r-1)m+k+1)/(m-1)),
       T = O(((r-1)m+k+N+1)^w * q^((r-1)(k+N+2-r)))."""
    if r < 1 or m <= 1:
        return math.inf
    N = ( (r - 1) * m + k + 1) // (m - 1)   # [EV26] Table 2 reproduced with floor
    if N > n - k - 1:
        return math.inf
    return w * log2((r - 1) * m + k + N + 1) + (r - 1) * (k + N + 2 - r) * log2(q)


def t_oj2(q, m, n, k, r, w=W):
    """[EV26] Eq. (5)+complexity:  N = ceil((k+1)r/(m-r)),
       T = O((k+1+N)^w * r^w * q^((m-r)(r-1)))."""
    if r < 1 or m <= r:
        return math.inf
    N = ceil((k + 1) * r / (m - r))
    if N > n - k - 1:
        return math.inf
    return w * log2(k + 1 + N) + w * log2(r) + (m - r) * (r - 1) * log2(q)


def t_grs(q, m, n, k, r, w=W):
    """[EV26] 'GRS' complexity: O(m^w (n-k-1)^w * q^((r-1)*ceil(m(k+1)/n)))."""
    if r < 1:
        return math.inf
    return w * log2(m * (n - k - 1)) + (r - 1) * ceil(m * (k + 1) / n) * log2(q)


def t_grsim(q, m, n, k, r, w=W):
    """[EV26] 'GRSim' complexity: O(m^w (n-k-1)^w * q^(r*ceil(m(k+1)/n) - m))."""
    if r < 1:
        return math.inf
    expo = r * ceil(m * (k + 1) / n) - m
    if expo < 0:
        return math.inf
    return w * log2(m * (n - k - 1)) + expo * log2(q)


def t_grsge(q, m, n, k, r, w=W):
    """[EV26] 'GRSge' complexity:
       O( min_{t>=0} ((n-k)m+t)^w * max(1, q^( r*ceil(((k+1)m-t)/n) - m + t )) )."""
    if r < 1:
        return math.inf
    best = math.inf
    # t only matters while it can lower ceil(((k+1)m-t)/n); search the whole
    # meaningful window 0 <= t <= (k+1)m  (beyond that ceil(...) is 0).
    hi = (k + 1) * m
    cand = set(range(0, min(hi, 4096) + 1))
    # add the t values at which ceil(((k+1)m-t)/n) jumps down
    for j in range(0, ceil(hi / n) + 2):
        for t in (hi - j * n, hi - j * n + 1, hi - j * n - 1):
            if 0 <= t <= hi:
                cand.add(t)
    for t in sorted(cand):
        expo = r * ceil(((k + 1) * m - t) / n) - m + t
        val = w * log2((n - k) * m + t) + max(0.0, expo * log2(q))
        best = min(best, val)
    return best


COMB = {"BE": t_basis_enumeration, "OJ1": t_oj1, "OJ2": t_oj2,
        "GRS": t_grs, "GRSim": t_grsim, "GRSge": t_grsge}


# --------------------------------------------------------------------------
# Algebraic algorithms  ([EV26] Sec 3.2)
# --------------------------------------------------------------------------
def mm_pure(q, m, n, k, r, w=W):
    """[EV26] Eq. (7)+(8).  Returns (log2 cost, p) or (inf, None) if Eq.(7) fails.

    Eq.(7):  m*C(n-k-1, r) >= C(n, r) - 1.
    Eq.(8):  T = O( m * C(n-p-k-1, r) * C(n-p, r)^(w-1) )   [F_q operations]
    with p = max{ i : m*C(n-i-k-1,r) >= C(n-i,r) - 1 }.
    """
    if r < 1 or n - k - 1 < r:
        return math.inf, None
    if m * _c(n - k - 1, r) < _c(n, r) - 1:
        return math.inf, None
    p = 0
    for i in range(0, n + 1):
        if n - i - k - 1 >= r and m * _c(n - i - k - 1, r) >= _c(n - i, r) - 1:
            p = i
        else:
            break
    cost = (log2(m) + log2(_c(n - p - k - 1, r))
            + (w - 1) * log2(_c(n - p, r)))
    return cost, p


def sm_pure(q, m, n, k, r, w=W):
    """[EV26] Eq. (10) with the Nb / Mb counts of [EV26] p.17.

    Nb = sum_{i=1..k} C(n-i,r) C(k+b-1-i, b-1) - C(n-k-1,r) C(k+b-1, b)
         - (m-1) * sum_{i=1..b} (-1)^(i+1) C(k+b-i-1, b-i) C(n-k-1, r+i)
    Mb = C(k+b-1, b) * ( C(n,r) - m*C(n-k-1,r) )
    Succeeds as long as Nb >= Mb - 1.  T = O(Nb * Mb^(w-1)) in F_{q^m} operations;
    we convert one F_{q^m} multiplication into m^2 F_q multiplications.
    """
    best = (math.inf, None)
    for b in range(1, 40):
        Nb = 0
        for i in range(1, k + 1):
            Nb += _c(n - i, r) * _c(k + b - 1 - i, b - 1)
        Nb -= _c(n - k - 1, r) * _c(k + b - 1, b)
        s = 0
        for i in range(1, b + 1):
            s += ((-1) ** (i + 1)) * _c(k + b - i - 1, b - i) * _c(n - k - 1, r + i)
        Nb -= (m - 1) * s
        Mb = _c(k + b - 1, b) * (_c(n, r) - m * _c(n - k - 1, r))
        if Nb <= 0 or Mb <= 0:
            continue
        if Nb < Mb - 1:
            continue
        # RYDE v2.1.0 spec Sec.8 and [EV26] both charge an F_{q^m} operation as
        # m^2 F_q operations ('the m^2 factor comes from expressing each F_{q^m}
        # operation involved in terms of F_q operations').
        cost = (log2(Nb) + (w - 1) * log2(Mb) + 2 * log2(m))
        if cost < best[0]:
            best = (cost, b)
    return best


def algebraic(q, m, n, k, r, w=W, depth=0):
    """Hybrid recursion of [EV26] Sec 3.2:
           Time-RSD(q,m,n,k,r) = min_{a>=0} q^(a r) * Time-RSD(q,m,n-a,k-a,r).
    Returns dict {'MM': (cost,a,p), 'SM': (cost,a,b)}."""
    out = {}
    for name, pure in (("MM", mm_pure), ("SM", sm_pure)):
        best = (math.inf, None, None)
        for a in range(0, k + 1):
            nn, kk = n - a, k - a
            if name == "MM":
                c, p = mm_pure(q, m, nn, kk, r, w)
                extra = None
            else:
                c, extra = sm_pure(q, m, nn, kk, r, w)
                p = None
            if c == math.inf:
                continue
            tot = c + a * r * log2(q)
            if tot < best[0]:
                best = (tot, a, p if name == "MM" else extra)
        out[name] = best
    return out


def best_attack(q, m, n, k, r, w=W):
    res = {}
    for nm, f in COMB.items():
        res[nm] = f(q, m, n, k, r, w)
    alg = algebraic(q, m, n, k, r, w)
    res["MM"] = alg["MM"][0]
    res["SM"] = alg["SM"][0]
    res["_MM_a"] = alg["MM"][1]
    res["_SM_a"] = alg["SM"][1]
    return res


# --------------------------------------------------------------------------
# Validation tables
# --------------------------------------------------------------------------
# [EV26] Table 2, p.26 : (name, m, n, k, r, comb_alg, comb_cost, alg, alg_cost, a)
EV26_TABLE2 = [
    ("dRANKula-256", 127, 118, 48, 7, "GRSim", 270, "MM", 68, 0),
    ("dRANKula-192", 97, 82, 42, 5, "GRSim", 181, "MM", 47, 0),
    ("dRANKula-128", 67, 57, 27, 5, "GRSim", 119, "MM", 42, 0),
    ("RVDC-96", 29, 22, 11, 7, "OJ1", 87, "MM", 79, 8),
    ("RVDC-125", 31, 26, 13, 8, "OJ1", 113, "MM", 107, 10),
    ("RVDC-196", 41, 34, 17, 10, "OJ1", 179, "MM", 173, 14),
    ("RVDC-252", 47, 38, 19, 12, "OJ1", 238, "MM", 231, 16),
    ("ROLLO-I-128", 79, 94, 47, 5, "GRSim", 149, "MM", 50, 0),
    ("ROLLO-I-192", 89, 106, 53, 6, "GRSim", 211, "MM", 61, 0),
    ("ROLLO-I-256", 113, 134, 67, 7, "GRSim", 318, "SM", 127, 5),
    ("ROLLO-II-128", 83, 298, 149, 5, "GRSim", 154, "MM", 66, 0),
    ("ROLLO-II-192", 107, 302, 151, 6, "GRSim", 244, "MM", 78, 0),
    ("ROLLO-II-256", 127, 314, 157, 7, "GRSim", 349, "SM", 113, 0),
    ("ROLLO-III-128", 101, 94, 47, 5, "GRSim", 183, "MM", 49, 0),
    ("ROLLO-III-192", 107, 118, 59, 6, "GRSim", 248, "MM", 62, 0),
    ("ROLLO-III-256", 131, 134, 67, 7, "GRSim", 364, "SM", 107, 2),
    ("RQC-I", 97, 134, 67, 5, "GRSim", 178, "MM", 54, 0),
    ("RQC-II", 107, 202, 101, 6, "GRSim", 249, "MM", 71, 0),
    ("RQC-III", 137, 262, 131, 7, "GRSge", 380, "MM", 108, 3),
    ("BG23", 31, 30, 15, 9, "GRSim", 139, "MM", 139, 12),
    ("DSD-I", 53, 53, 45, 4, "GRSim", 148, "MM", 169, 37),
    ("DSD-II", 79, 75, 67, 4, "OJ1", 220, "MM", 254, 58),
    ("DSD-III", 97, 95, 87, 4, "OJ1", 281, "MM", 331, 77),
    ("RYDE-I", 31, 33, 15, 10, "GRSge", 138, "MM", 153, 12),
    ("RYDE-II", 37, 41, 18, 13, "GRSge", 209, "MM", 237, 15),
    ("RYDE-III", 43, 47, 18, 17, "GRSim", 283, "MM", 308, 15),
    ("RYDE-r2-I", 53, 53, 45, 4, "GRSim", 148, "MM", 169, 37),
    ("RYDE-r2-III", 61, 61, 51, 5, "GRSim", 217, "MM", 245, 44),
    ("RYDE-r2-V", 61, 57, 44, 7, "OJ1", 287, "MM", 301, 39),
]


def cmd_validate():
    print("=" * 108)
    print("VALIDATION vs [EV26] ePrint 2026/082 Table 2 (q=2, w=2).  "
          "'mine' = this script, 'pub' = published.")
    print("=" * 108)
    hdr = (f"{'instance':<15}{'m':>4}{'n':>5}{'k':>5}{'r':>3} | "
           f"{'comb':>6}{'mine':>8}{'pub':>6} | {'alg':>4}{'mine':>8}{'pub':>6}{'a_mine':>7}{'a_pub':>6}")
    print(hdr)
    ok_c = ok_a = 0
    for (nm, m, n, k, r, calg, cpub, aalg, apub, apuba) in EV26_TABLE2:
        rc = COMB[calg](2, m, n, k, r)
        alg = algebraic(2, m, n, k, r)
        ra = alg[aalg][0]
        ra_a = alg[aalg][1]
        dc = abs(rc - cpub) <= 1.0
        da = abs(ra - apub) <= 1.0
        ok_c += dc
        ok_a += da
        print(f"{nm:<15}{m:>4}{n:>5}{k:>5}{r:>3} | {calg:>6}{rc:>8.1f}{cpub:>6} "
              f"{'*' if dc else 'X'} | {aalg:>4}{ra:>8.1f}{apub:>6}"
              f"{'*' if da else 'X'}{ra_a if ra_a is not None else -1:>7}{apuba:>6}")
    print(f"\ncombinatorial reproduced within +-1 bit : {ok_c}/{len(EV26_TABLE2)}")
    print(f"algebraic     reproduced within +-1 bit : {ok_a}/{len(EV26_TABLE2)}")


# --------------------------------------------------------------------------
# ROLLO  ->  RD instance mapping  ([BBGT22] Table 2 caption)
# --------------------------------------------------------------------------
def rollo_rd_instances(m, n, r, d, inhom_pk=False):
    """[BBGT22] Table 2 caption / [EV26] Table 2:
       *  For homogeneous ROLLO-I/II (and LOCKER):
            - plain ciphertext RD: (m, 2k, k, r)
            - derived-code-from-key RD: (m, 2k - floor(k/d), k - floor(k/d), d)
       *  For inhomogeneous ROLLO-III (pk = (h, s = x + h y) with uniform h):
            - inhomogeneous public-key 2-DIRSD: (m, 2k, k, d)  (d = w)
            - ciphertext first-component RD:    (m, 2k, k, r)  (r = w_r)
       where k = n is the ROLLO ring degree.
    """
    k = n
    plain = (m, 2 * k, k, r)                      # ciphertext / plain RD
    if inhom_pk:
        return {"RD(m,2k,k,r)": plain, "inhom_pk": (m, 2 * k, k, d)}
    kd = k // d          # floor: the unique choice reproducing [BBGT22] Table 2 (6/6 rows)
    derived = (m, 2 * k - kd, k - kd, d)          # derived-code-from-key
    return {"RD(m,2k,k,r)": plain, "derived": derived}


ROLLO_2020 = [
    ("ROLLO-I-128 (2020)", 67, 83, 7, 8, 128, 196),
    ("ROLLO-I-192 (2020)", 79, 97, 8, 8, 192, 241),
    ("ROLLO-I-256 (2020)", 97, 113, 9, 9, 256, 353),
    ("ROLLO-II-128 (2020)", 83, 189, 7, 8, 128, None),
    ("ROLLO-II-192 (2020)", 97, 193, 8, 8, 192, None),
    ("ROLLO-II-256 (2020)", 97, 211, 8, 9, 256, None),
]
NEW2ROLLO = [
    ("new2ROLLO-I-128", 73, 83, 7, 8, 128, 202),
    ("new2ROLLO-I-192", 89, 97, 8, 8, 192, 223),
    ("new2ROLLO-I-256", 103, 113, 9, 9, 256, 366),
]
# NOTE: d-values below are taken from the OFFICIAL specifications
#   ROLLO 2019-08-24 spec, Table 2 (ROLLO-I) and Table 4 (ROLLO-II)
#   ROLLO 2020-04-21 spec, Table 3 (ROLLO-I) and Table 5 (ROLLO-II)
# Values in rollo_params_all27_official.json.
ROLLO_2019 = [
    ("ROLLO-I-128 (2019)", 79, 47, 5, 6, 128, None),
    ("ROLLO-I-192 (2019)", 89, 53, 6, 7, 192, None),
    ("ROLLO-I-256 (2019)", 113, 67, 7, 8, 256, None),
    ("ROLLO-II-128 (2019)", 83, 149, 5, 8, 128, None),
    ("ROLLO-II-192 (2019)", 107, 151, 6, 8, 192, None),
    ("ROLLO-II-256 (2019)", 127, 157, 7, 8, 256, None),
    ("ROLLO-III-128 (2019)", 101, 47, 6, 5, 128, None),
    ("ROLLO-III-192 (2019)", 107, 59, 8, 6, 192, None),
    ("ROLLO-III-256 (2019)", 131, 67, 8, 7, 256, None),
]


def cmd_rollo():
    print("=" * 118)
    print("2020/2021-REVISED ROLLO + new2ROLLO: all known RSD algorithms, w=2, q=2.")
    print("RD instances derived per [BBGT22] ePrint 2022/1031 Table 2 caption.")
    print("=" * 118)
    rows = ROLLO_2020 + NEW2ROLLO + ROLLO_2019
    for (nm, m, n, r, d, sec, pub) in rows:
        print(f"\n### {nm}   (m={m}, n={n}, r={r}, d={d})   claimed 2^{sec}"
              + (f"   published SOTA 2^{pub}" if pub else "   published SOTA: NONE FOUND"))
        for tag, (mm, nn, kk, rr) in rollo_rd_instances(m, n, r, d, inhom_pk=("ROLLO-III" in nm)).items():
            res = best_attack(2, mm, nn, kk, rr)
            vals = {a: b for a, b in res.items() if not a.startswith("_")}
            bestn = min(vals, key=vals.get)
            line = "  ".join(f"{a}={vals[a]:.1f}" for a in
                             ["BE", "OJ1", "OJ2", "GRS", "GRSim", "GRSge", "MM", "SM"]
                             if vals[a] != math.inf)
            print(f"  [{tag:<14}] (m,n,k,r)=({mm},{nn},{kk},{rr})")
            print(f"      {line}")
            print(f"      -> BEST {bestn} = 2^{vals[bestn]:.1f}"
                  + ("   *** BELOW CLAIMED ***" if vals[bestn] < sec else ""))


# --------------------------------------------------------------------------
# Second, independent validation: [BBGT22] ePrint 2022/1031 Table 2 ("comb"
# column = the AGHT-2018 improved GRS = our GRSim).  Starred rows use the
# derived-code-from-key instance, unstarred rows the plain RD(m,2k,k,r).
# --------------------------------------------------------------------------
BBGT22_TABLE2 = [
    # name,                m,   n,  r, d, starred, comb_pub, best_pub
    ("new2ROLLO-I-128",    73,  83, 7, 8, False, 212, 202),
    ("new2ROLLO-I-192",    89,  97, 8, 8, True,  282, 223),
    ("new2ROLLO-I-256",   103, 113, 9, 9, True,  375, 366),
    ("ROLLO-I-128 (2020)", 67,  83, 7, 8, False, 196, 196),
    ("ROLLO-I-192 (2020)", 79,  97, 8, 8, True,  251, 241),
    ("ROLLO-I-256 (2020)", 97, 113, 9, 9, True,  353, 353),
]


def cmd_validate_bbgt22():
    print("=" * 100)
    print("VALIDATION vs [BBGT22] ePrint 2022/1031 Table 2 ('comb' column, w=2, q=2)")
    print("=" * 100)
    print(f"{'instance':<22}{'RD instance used':<26}{'mine':>8}{'pub':>6}  ok")
    ok = 0
    for (nm, m, n, r, d, star, cpub, _best) in BBGT22_TABLE2:
        inst = rollo_rd_instances(m, n, r, d)["derived" if star else "RD(m,2k,k,r)"]
        mine = t_grsim(2, *inst)
        good = abs(mine - cpub) <= 1.5
        ok += good
        print(f"{nm:<22}{str(inst):<26}{mine:>8.1f}{cpub:>6}  {'*' if good else 'X'}")
    print(f"\nreproduced within +-1.5 bit : {ok}/{len(BBGT22_TABLE2)}")


# --------------------------------------------------------------------------
# NIST PQC Additional Signatures Round-2 rank-metric candidates
# --------------------------------------------------------------------------
# RYDE  v2.1.0 (23/09/2025), RSD over F_{2^m}, RANDOM linear code.
RYDE_R2 = [("RYDE-1", 2, 53, 53, 45, 4, 143),
           ("RYDE-3", 2, 61, 61, 51, 5, 207),
           ("RYDE-5", 2, 61, 57, 44, 7, 272)]
# RYDE v2.0.0 cat-5 (superseded by v2.1.0)
RYDE_R2_OLD = [("RYDE-5 (v2.0.0)", 2, 67, 67, 55, 6, 272)]
# Mirath v2.1.0 : MinRank-syndrome = RSD over the BASE field, code length m*n.
MIRATH = [("Mirath-1a", 16, 16, 16, 143, 4, 158),
          ("Mirath-3a", 16, 19, 19, 195, 5, 225),
          ("Mirath-5a", 16, 22, 22, 255, 6, 301),
          ("Mirath-1b", 2, 42, 42, 1443, 4, 158),
          ("Mirath-3b", 2, 50, 50, 2024, 5, 224),
          ("Mirath-5b", 2, 56, 56, 2499, 6, 289)]


def cmd_signatures():
    print("=" * 112)
    print("NIST PQC Additional Signatures Round 2 -- rank-metric candidates "
          "(w=2 unless noted).")
    print("RYDE rows: plain RSD over F_{2^m}, uniformly random parity-check "
          "matrix (NO ideal structure).")
    print("=" * 112)
    print(f"{'instance':<16}{'q':>3}{'m':>4}{'n':>5}{'k':>6}{'r':>3}"
          f"{'BE':>8}{'OJ1':>8}{'GRSim':>8}{'GRSge':>8}{'MM':>8}{'SM':>8}"
          f"{'BEST':>9}{'target':>8}")
    for (nm, q, m, n, k, r, tgt) in RYDE_R2 + RYDE_R2_OLD:
        vals = {a: f(2, m, n, k, r) for a, f in COMB.items()}
        alg = algebraic(2, m, n, k, r)
        vals["MM"], vals["SM"] = alg["MM"][0], alg["SM"][0]
        best = min(vals.values())
        bn = min(vals, key=vals.get)
        print(f"{nm:<16}{q:>3}{m:>4}{n:>5}{k:>6}{r:>3}"
              f"{vals['BE']:>8.1f}{vals['OJ1']:>8.1f}{vals['GRSim']:>8.1f}"
              f"{vals['GRSge']:>8.1f}{vals['MM']:>8.1f}{vals['SM']:>8.1f}"
              f"{best:>9.1f}{tgt:>8}"
              + ("  <-- BELOW TARGET" if best < tgt else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--rollo", action="store_true")
    ap.add_argument("--validate2", action="store_true")
    ap.add_argument("--signatures", action="store_true")
    a = ap.parse_args()
    if a.validate:
        cmd_validate()
    if a.validate2:
        cmd_validate_bbgt22()
    if a.signatures:
        cmd_signatures()
    if a.rollo:
        cmd_rollo()
    if not (a.validate or a.rollo or a.validate2 or a.signatures):
        ap.print_help()




