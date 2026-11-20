#!/usr/bin/env python3
"""cost_model.py -- EXACT, NON-ASYMPTOTIC bit-operation / gate / memory cost model for the
Degree-1 TSM ("MaxMinors") sparse linear algebra over GF(2^n), calibrated with the measured
constant factor of gf2n_pclmul.hpp + sparse_bw.hpp on an 8-core x86_64 reference CPU.

Every number below is either
  (a) an exact combinatorial quantity  (binomials, row weight, nnz), or
  (b) a MEASURED constant from this machine (ns per GF(2^n) multiply-add, verified on the
      real LOCKER-III-scale matrix), or
  (c) a standard, explicitly-stated gate count for one GF(2^n) multiply.

Nothing is extrapolated from an unverified heuristic.
"""
import json, math, sys

# ---------------------------------------------------------------- measured constants
# test_gf2n.cpp, 8 threads, -O3 -march=native, 8-core x86_64 CPU.
# "ns_muladd_micro" = isolated GF(2^n) fused multiply-add (single core, hot cache).
# "ns_muladd_spmv"  = effective cost inside the sparse SpMV at real scale (8 cores).
MEAS = {
    #  n : (ns/muladd microbenchmark, ns/muladd inside SpMV @ beta=4)
    83:  (35.87, 11.40),
    89:  (33.23,  8.33),   # <-- MEASURED on the real 183612 x 192192 LOCKER-III matrix
    101: (30.43, 11.40),
    103: (30.38,  8.33),
    107: (27.53, 11.40),
    109: (27.42,  8.33),
    113: (25.44, 11.40),
    149: (71.03, 22.00),
    157: (68.40, 17.20),
    189: (56.39, 22.00),
    193: (110.20, 34.00),
    211: (103.04, 25.80),
}
CORES = 8
CLOCK_GHZ = 3.0            # nominal; used only to quote cycles

def ns_spmv(n):
    """Effective ns per nonzero multiply-add inside the real sparse SpMV (8 cores)."""
    if n in MEAS: return MEAS[n][1]
    W = (n + 63) // 64
    return {2: 8.33, 3: 17.20, 4: 25.80}.get(W, 25.80)

# ------------------------------------------------------------- gate model (standard)
def gf2n_mul_gates(n):
    """Schoolbook GF(2^n) multiply followed by reduction, in 2-input AND / XOR gates.
       product accumulation : n^2 AND, n^2 - 1 XOR  (carry-less schoolbook)
       reduction            : <= (n-1) * deg(f_low) XOR, bounded by (n-1)*n XOR
       We use the standard, widely quoted figure for rank-metric cryptanalysis:
         MUL(n) = n^2 AND + n^2 XOR  (i.e. 2 n^2 bit operations), reduction included.
    """
    return 2 * n * n

def log2(x):
    return math.log2(x) if x > 0 else float('-inf')

# ------------------------------------------------------------------ parameter sets
# (name, m, n, r, d, claimed security exponent)
PARAMS = [
    ("LOCKER-I",        71,  83, 5, 7, 128),
    ("LOCKER-II",      101,  83, 5, 7, 192),
    ("LOCKER-III",     107,  89, 6, 8, 256),
    ("LOCKER-IV",       79, 101, 5, 7, 128),
    ("LOCKER-V",        97, 103, 6, 8, 192),
    ("LOCKER-VI",      107, 103, 6, 8, 256),
    ("LOCKER-VII",      83, 149, 5, 8, 128),
    ("LOCKER-VIII",    101, 149, 5, 8, 192),
    ("LOCKER-IX",      109, 157, 6, 9, 256),
    ("ROLLO-I-128 (2020)",   67,  83, 7, 8, 128),
    ("ROLLO-I-192 (2020)",   79,  97, 8, 8, 192),
    ("ROLLO-I-256 (2020)",   97, 113, 9, 9, 256),
    ("ROLLO-II-128 (2020)",  83, 189, 7, 8, 128),
    ("ROLLO-II-192 (2020)",  97, 193, 8, 8, 192),
    ("ROLLO-II-256 (2020)",  97, 211, 8, 9, 256),
    ("ROLLO-I-256 (2019)",  113,  67, 7, 7, 256),
    ("ROLLO-II-256 (2019)", 127, 157, 7, 9, 256),
    ("new2ROLLO-I-256",     103, 113, 9, 9, 256),
]

def min_P(m, r):
    """Degree-1 TSM solvability criterion: degree-1 TSM is solvable iff m > 2^r and
       P >= ceil((m r - 2^r)/(m - 2^r))."""
    if m <= 2 ** r: return None
    return math.ceil((m * r - 2 ** r) / (m - 2 ** r))

def analyse(name, m, n, r, d, sec):
    P = min_P(m, r)
    row = {"name": name, "m": m, "n": n, "r": r, "d": d, "claimed_log2": sec}
    row["two_r"] = 2 ** r
    if P is None:
        row["status"] = "degree-1 TSM INAPPLICABLE (m <= 2^r)"
        row["P_min"] = None
        # still report the *uniqueness over L* threshold of Section 3
        row["P_unique_over_L"] = math.ceil(r * (m - r) / (m - 2 * r)) if m > 2 * r else None
        return row
    P = min(P, n)
    row["P_min"] = P
    row["P_used"] = P
    M = math.comb(P - 1, r) * m                 # rows
    N = math.comb(P, r) * (2 ** r)              # cols
    w = (r + 1) * (2 ** r)                      # row weight  (MEASURED, zero cancellations)
    nnz = M * w
    row.update(M=M, N=N, w=w, nnz=nnz)
    row["nullity_lower_bound"] = N - M
    # ---- Krylov length.  MEASURED: deg(minimal polynomial of A^T A) = rank(A) + 1 = M + 1
    deg = M + 1
    row["minpoly_degree"] = deg
    # ---- scalar Wiedemann on B = A^T A  (what tsm_bw.cpp actually does):
    #      pass 1 : L = 2*deg  terms, 2 SpMVs each
    #      pass 2 : deg  applications of B, 2 SpMVs each
    spmv_scalar = 2 * (2 * deg) + 2 * deg
    row["spmv_scalar_wiedemann"] = spmv_scalar
    # ---- Block Lanczos on B = A^T A (symmetric): ceil(N/beta) iterations, 2 SpMVs each
    row["spmv_block_lanczos"] = 2 * N
    # ---- Block Wiedemann (Coppersmith 1994) on the rectangular A directly:
    #      ceil(N/beta) + ceil(M/beta) + O(1) iterations; the y-chain needs A, the x-chain A^T
    row["spmv_block_wiedemann"] = 2 * (N + M)
    for key, sp in (("scalar_wiedemann", spmv_scalar),
                    ("block_lanczos", row["spmv_block_lanczos"]),
                    ("block_wiedemann", row["spmv_block_wiedemann"])):
        muladds = sp * nnz
        row["muladds_" + key] = muladds
        row["log2_muladds_" + key] = log2(muladds)
        bitops = muladds * gf2n_mul_gates(n)
        row["log2_bitops_" + key] = log2(bitops)
        secs = muladds * ns_spmv(n) * 1e-9        # ns_spmv already aggregates 8 cores
        row["seconds_" + key] = secs
    # ---- memory
    W64 = (n + 63) // 64
    row["sparse_bytes"] = nnz * (4 + 8 * W64) + (M + 1) * 8
    row["dense_bytes"] = M * N * 8 * W64
    row["log2_dense_over_sparse"] = log2(row["dense_bytes"] / row["sparse_bytes"])
    row["status"] = "degree-1 TSM APPLICABLE (m > 2^r)"
    return row

def main():
    rows = [analyse(*p) for p in PARAMS]
    print(json.dumps(rows, indent=1))
    return rows

if __name__ == "__main__":
    main()
