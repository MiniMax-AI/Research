#!/usr/bin/env python3
"""
exp2_rank_truth_table.py -- EXACT GF(2) rank truth table for the TSM support-minor
(Plucker) linearisation, TRUE Frobenius  H_q = rho^q(H_0),  U_q = rho^q(U_0)
versus the RANDOM model used by the degree-2 bi-monomial model.

Linearisation (identical to parameterized_tsm/tsm.cpp):
  columns : (Q, mask),  Q in C([P],r),  mask in {0,1}^r        -> C(P,r)*2^r unknowns over L
  rows    : for every (r+1)-subset S' of [P] containing 0, the Q-determinant
              det[ C_{q_1} ; ... ; C_{q_r} ; U ] = 0,  C_{q,k} = w_{2q,k} + H_q w_{2q+1,k}
            expanded along the last row and then over the 2^r masks
              => |C([P] containing 0, r+1)| * m  rows over L.
  Frobenius relations (valid ONLY in the TRUE system):
              p_{(Q+1,mask)} = sigma( p_{(Q,mask)} ),   sigma = squaring in L = GF(2^n)
            (GF(2)-linear, NOT L-linear).

Everything is flattened to an exact GF(2) matrix and ranked with M4RI.
"""
import sys, os, json, itertools, random, struct, subprocess, argparse, time, tempfile, shutil, atexit
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rollo_gf2_engine import RingQ
from exp1_census import make_instance, canon

_RANKBIN_CACHE = None
_RANKBIN_TRIED = False

def _get_rankbin():
    global _RANKBIN_CACHE, _RANKBIN_TRIED
    env_bin = os.environ.get("GF2RANK")
    if env_bin and os.path.isfile(env_bin) and os.access(env_bin, os.X_OK):
        return env_bin
    if _RANKBIN_TRIED:
        return _RANKBIN_CACHE
    _RANKBIN_TRIED = True
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "gf2rank.c")
    if not os.path.isfile(src):
        return None
    cc = os.environ.get("CC", "gcc")
    cflags_libs = None
    try:
        pc = subprocess.run(["pkg-config", "--cflags", "--libs", "m4ri"],
                            capture_output=True, text=True)
        if pc.returncode == 0:
            cflags_libs = pc.stdout.split()
    except OSError:
        pass
    if cflags_libs is None:
        conda = os.environ.get("CONDA_PREFIX", "")
        for pfx in ([conda] if conda else []) + ["/opt/homebrew", "/usr/local"]:
            if os.path.isfile(os.path.join(pfx, "include", "m4ri", "m4ri.h")):
                cflags_libs = [f"-I{pfx}/include", f"-L{pfx}/lib", f"-Wl,-rpath,{pfx}/lib", "-lm4ri", "-lm"]
                break
    if cflags_libs is None:
        return None
    tmpdir = tempfile.mkdtemp(prefix="gf2rank_")
    atexit.register(shutil.rmtree, tmpdir, True)
    out_bin = os.path.join(tmpdir, "gf2rank")
    res = subprocess.run([cc, "-O3", "-o", out_bin, src, *cflags_libs],
                         capture_output=True, text=True)
    if res.returncode == 0 and os.path.isfile(out_bin):
        _RANKBIN_CACHE = out_bin
    return _RANKBIN_CACHE

# ---------------------------------------------------------------- helpers
def mul_tables(R):
    """POW[a][j] = a * X^j  in L = GF(2^n); BLK[a][i] = sum_j bit_i(a X^j) << j."""
    L, n = R.L, R.n
    P = L.mod
    xpow = [0] * n
    x = 1
    for j in range(n):
        xpow[j] = x
        x = L.mul(x, 2) if False else (x << 1) ^ (P if (x >> (n - 1)) & 1 else 0)
        x &= (1 << n) - 1
    size = 1 << n
    POW = np.zeros((size, n), dtype=np.uint32)
    for a in range(size):
        for j in range(n):
            POW[a, j] = L.mul(a, xpow[j]) if a else 0
    BLK = np.zeros((size, n), dtype=np.uint32)
    for i in range(n):
        BLK[:, i] = np.bitwise_or.reduce(((POW >> i) & 1) << np.arange(n, dtype=np.uint32), axis=1)
    return POW, BLK

def sigma_blk(R):
    """SIGBLK[i] = sum_j bit_i(X^(2j) mod P) << j   (matrix of z -> z^2 on L)."""
    L, n = R.L, R.n
    P = L.mod
    out = []
    vals = []
    x = 1
    for j in range(n):
        vals.append(x)
        x = (x << 1) ^ (P if (x >> (n - 1)) & 1 else 0)
        x &= (1 << n) - 1
    sq = []
    for j in range(n):
        y = vals[j]
        # X^(2j) mod P
        r = 1
        for _ in range(2 * j):
            r = (r << 1) ^ (P if (r >> (n - 1)) & 1 else 0)
            r &= (1 << n) - 1
        sq.append(r)
    for i in range(n):
        v = 0
        for j in range(n):
            if (sq[j] >> i) & 1: v |= 1 << j
        out.append(v)
    return out

def build_laplace_L(R, Hs, Us, P, r):
    """Return (fives, sixes, A) with A a numpy uint32 array (nrows_L, ncols) of L-elements."""
    m = R.m
    fives = list(itertools.combinations(range(P), r))
    fidx = {t: i for i, t in enumerate(fives)}
    sixes = [s for s in itertools.combinations(range(P), r + 1) if s[0] == 0]
    ncols = len(fives) * (1 << r)
    A = np.zeros((len(sixes) * m, ncols), dtype=np.uint32)
    for eq, sel in enumerate(sixes):
        for omit in range(r + 1):
            pairs = tuple(x for i, x in enumerate(sel) if i != omit)
            block = fidx[pairs]
            for sm in range(1 << r):
                coef = Us[sel[omit]]
                for k in range(r):
                    if (sm >> k) & 1: coef = R.qmul(coef, Hs[pairs[k]])
                col = block * (1 << r) + sm
                for d in range(m):
                    A[eq * m + d, col] = coef[d]
    return fives, sixes, A

def pack_bits(b):
    """(nrows, nbits) uint8 array of 0/1 -> list of python ints, logical bit k = integer bit k."""
    packed = np.packbits(b, axis=1, bitorder="little")
    return [int.from_bytes(row.tobytes(), "little") for row in packed]

def flatten(A, BLK, n):
    """A (nrows_L, ncols) of L-elements -> list of GF(2) row ints of width ncols*n.
       L-entry a becomes the n x n GF(2) block M_a with (M_a)_{ij} = bit_i(a*X^j)."""
    nrows, ncols = A.shape
    out = []
    ar = np.arange(n, dtype=np.uint32)
    for i in range(n):
        blk = BLK[A, i]                      # (nrows, ncols) uint32, each < 2^n
        b = ((blk[:, :, None] >> ar) & 1).astype(np.uint8).reshape(nrows, ncols * n)
        out.extend(pack_bits(b))
    return out

def gf2_rank_py(rows, ncols):
    """pure-python GF(2) rank (bitset rows); used when the M4RI helper is unavailable."""
    rows = list(rows); rk = 0
    for c in range(ncols):
        p = -1
        for i in range(rk, len(rows)):
            if (rows[i] >> c) & 1: p = i; break
        if p < 0: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        pr = rows[rk]
        for i in range(rk + 1, len(rows)):
            if (rows[i] >> c) & 1: rows[i] ^= pr
        rk += 1
        if rk == len(rows): break
    return rk

def gf2_rank(rows, ncols, tag=""):
    rankbin = _get_rankbin()
    if not rankbin:
        return gf2_rank_py(rows, ncols)
    wpr = (ncols + 63) // 64
    with tempfile.NamedTemporaryFile(prefix=f"gf2mat_{tag}_", suffix=".bin", delete=False) as f:
        path = f.name
        f.write(struct.pack("<QQ", len(rows), ncols))
        for r in rows:
            f.write(r.to_bytes(wpr * 8, "little"))
    try:
        out = subprocess.check_output([rankbin, path]).decode().split()
        return int(out[0])
    finally:
        if os.path.exists(path):
            os.unlink(path)

def pack_rows(list_of_ints, ncols):
    return list_of_ints

# ---------------------------------------------------------------- true Plucker vector
def true_plucker(R, alphas, betas, P, r):
    """p[(Q,mask)] = det of the r x r matrix with rows w_{2q+bit}."""
    L, m = R.L, R.m
    ws = []
    a = list(alphas); b = list(betas)
    for q in range(P):
        ws.append((2 * q, a[:]))
        ws.append((2 * q + 1, b[:]))
        a = [L.fmul(x, x) for x in a]
        b = [L.fmul(x, x) for x in b]
    W = {i: v for i, v in ws}
    fives = list(itertools.combinations(range(P), r))
    p = {}
    for Q in fives:
        for sm in range(1 << r):
            idx = [2 * Q[k] + ((sm >> k) & 1) for k in range(r)]
            M = [W[i] for i in idx]
            p[(Q, sm)] = det_L(M, L)
    return p

def det_L(M, L):
    M = [row[:] for row in M]
    n = len(M)
    d = 1
    for c in range(n):
        piv = -1
        for i in range(c, n):
            if M[i][c]: piv = i; break
        if piv < 0: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]
        d = L.fmul(d, M[c][c])
        iv = L.finv(M[c][c])
        for i in range(c + 1, n):
            if M[i][c]:
                f = L.fmul(M[i][c], iv)
                M[i] = [x ^ L.fmul(f, y) for x, y in zip(M[i], M[c])]
    return d

# ---------------------------------------------------------------- main
def run(m, n, r, d, P, seed, verbose=True):
    R, inst = make_instance(m, n, r, d, seed)
    L = R.L
    H0, U0 = R.Phi(inst["h"]), R.Phi(inst["u"])
    # true alpha,beta
    E = inst["E"]
    # recover alpha_k, beta_k  (e_i = sum_k b_k a_k, a_k in L)
    def recover(e):
        cols = []
        for j in range(n):
            v = e[j]; rs = E[:]; tag = [1 << k for k in range(r)]; rk = 0
            for c in range(m):
                p = -1
                for i in range(rk, r):
                    if (rs[i] >> c) & 1: p = i; break
                if p < 0: continue
                rs[rk], rs[p] = rs[p], rs[rk]; tag[rk], tag[p] = tag[p], tag[rk]
                for i in range(r):
                    if i != rk and ((rs[i] >> c) & 1): rs[i] ^= rs[rk]; tag[i] ^= tag[rk]
                rk += 1
            t, s = v, 0
            for i in range(rk):
                c = (rs[i] & -rs[i]).bit_length() - 1
                if (t >> c) & 1: t ^= rs[i]; s ^= tag[i]
            assert t == 0
            cols.append(s)
        al = [0] * r
        for j in range(n):
            for k in range(r):
                if (cols[j] >> k) & 1: al[k] ^= (1 << j)
        return al
    alphas, betas = recover(inst["e1"]), recover(inst["e2"])
    Hs, Us = [H0], [U0]
    for q in range(1, P):
        Hs.append(R.rho(Hs[-1])); Us.append(R.rho(Us[-1]))
    # sanity: U_q = sum_k b~_k (rho^q alpha_k + H_q rho^q beta_k)
    for q in [0, min(1, P - 1), P - 1]:
        acc = [0] * m
        aq = list(alphas); bq = list(betas)
        for _ in range(q):
            aq = [L.fmul(x, x) for x in aq]; bq = [L.fmul(x, x) for x in bq]
        for k in range(r):
            t1 = R.K_to_Q(E[k])
            acc = [x ^ y for x, y in zip(acc, R.qmul(t1, [aq[k]] + [0] * (m - 1)))]
            acc = [x ^ y for x, y in zip(acc, R.qmul(t1, R.qmul(Hs[q], [bq[k]] + [0] * (m - 1))))]
        assert acc == Us[q], f"TSM identity failed at q={q}"

    POW, BLK = mul_tables(R)
    SIG = sigma_blk(R)
    res = {}

    # ---------- TRUE model
    fives, sixes, A = build_laplace_L(R, Hs, Us, P, r)
    ncolsL = A.shape[1]
    ptrue = true_plucker(R, alphas, betas, P, r)
    pvec = np.zeros(ncolsL, dtype=np.uint32)
    for ci, Q in enumerate(fives):
        for sm in range(1 << r):
            pvec[ci * (1 << r) + sm] = ptrue[(Q, sm)]
    # POSITIVE CONTROL: the true Plucker vector must satisfy A p = 0
    bad = 0
    nzc = np.nonzero(pvec)[0]
    for i in range(A.shape[0]):
        s = 0
        row = A[i]
        for c in nzc:
            if row[c]: s ^= L.mul(int(row[c]), int(pvec[c]))
        if s: bad += 1
    res["true_laplace_residual_rows"] = int(bad)
    # Frobenius relations on the true p
    fidx = {t: i for i, t in enumerate(fives)}
    frob_bad = 0
    nrel = 0
    for Q in itertools.combinations(range(P - 1), r):
        Q1 = tuple(x + 1 for x in Q)
        if Q1 not in fidx: continue
        for sm in range(1 << r):
            nrel += 1
            if ptrue[(Q1, sm)] != L.fmul(ptrue[(Q, sm)], ptrue[(Q, sm)]): frob_bad += 1
    res["true_frobenius_violations"] = int(frob_bad)
    res["n_frobenius_relations"] = int(nrel)

    rows_L = flatten(A, BLK, n)
    ncols_gf2 = ncolsL * n
    rk = gf2_rank(rows_L, ncols_gf2, "true")
    res["TRUE_laplace_only"] = dict(rows=len(rows_L), cols=ncols_gf2, rank=rk,
                                    nullity=ncols_gf2 - rk, nullity_in_L_units=(ncols_gf2 - rk) / n)
    # + Frobenius relations
    frows = []
    for Q in itertools.combinations(range(P - 1), r):
        Q1 = tuple(x + 1 for x in Q)
        if Q1 not in fidx: continue
        c0 = fidx[Q] * (1 << r); c1 = fidx[Q1] * (1 << r)
        for sm in range(1 << r):
            for i in range(n):
                v = (1 << i) << ((c1 + sm) * n)
                v ^= SIG[i] << ((c0 + sm) * n)
                frows.append(v)
    rk2 = gf2_rank(rows_L + frows, ncols_gf2, "truefrob")
    res["TRUE_laplace_plus_frobenius"] = dict(rows=len(rows_L) + len(frows), cols=ncols_gf2, rank=rk2,
                                              nullity=ncols_gf2 - rk2, nullity_in_L_units=(ncols_gf2 - rk2) / n)

    # ---------- RANDOM (probe) model
    rng = random.Random(seed ^ 0x5eed)
    Hr, Ur = [], []
    for q in range(P):
        Hq = [rng.randrange(1 << n) for _ in range(m)]
        aq = [rng.randrange(1 << n) for _ in range(r)]
        bq = [rng.randrange(1 << n) for _ in range(r)]
        Uq = [0] * m
        for k in range(r):
            t1 = R.K_to_Q(E[k])
            Uq = [x ^ y for x, y in zip(Uq, R.qmul(t1, [aq[k]] + [0] * (m - 1)))]
            Uq = [x ^ y for x, y in zip(Uq, R.qmul(t1, R.qmul(Hq, [bq[k]] + [0] * (m - 1))))]
        Hr.append(Hq); Ur.append(Uq)
    _, _, Ar = build_laplace_L(R, Hr, Ur, P, r)
    rows_R = flatten(Ar, BLK, n)
    rkr = gf2_rank(rows_R, ncols_gf2, "rand")
    res["RANDOM_laplace_only"] = dict(rows=len(rows_R), cols=ncols_gf2, rank=rkr,
                                      nullity=ncols_gf2 - rkr, nullity_in_L_units=(ncols_gf2 - rkr) / n)
    rkr2 = gf2_rank(rows_R + frows, ncols_gf2, "randfrob")
    res["RANDOM_laplace_plus_frobenius"] = dict(rows=len(rows_R) + len(frows), cols=ncols_gf2, rank=rkr2,
                                                nullity=ncols_gf2 - rkr2, nullity_in_L_units=(ncols_gf2 - rkr2) / n)
    res["P"] = P; res["m"] = m; res["n"] = n; res["r"] = r; res["seed"] = seed
    res["ncols_L"] = int(ncolsL); res["nrows_L"] = int(A.shape[0])
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--n", type=int, default=13)
    ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--P", type=str, default="3,4,5,6,7,8")
    ap.add_argument("--seed", type=int, default=777)
    ap.add_argument("--json-out", type=str, default="")
    args = ap.parse_args()
    allres = []
    print(f"=== EXACT GF(2) RANK TRUTH TABLE  m={args.m} n={args.n} r={args.r} (2r={2*args.r}, 2^r={2**args.r}) ===")
    import math
    print(f"{'P':>3} {'colsL':>7} {'rowsL':>7} | {'TRUElap':>8} {'TRUEfrob':>9} {'PREDICT':>8} | {'RANDlap':>8} {'RANDfrob':>9} | ctrl")
    for P in [int(x) for x in args.P.split(",")]:
        t0 = time.time()
        res = run(args.m, args.n, args.r, args.d, P, args.seed)
        allres.append(res)
        r_, m_, n_, P_ = args.r, args.m, args.n, P
        pred = n_ * ((1 << r_) * math.comb(P_ - 1, r_ - 1) - m_ * math.comb(P_ - 2, r_ - 1))
        res["PREDICTED_true_nullity"] = pred
        obs = res['TRUE_laplace_plus_frobenius']['nullity']
        flag = "OK" if obs == max(pred, 0) or (pred < 0 and obs == 0) else ("MISMATCH" if pred >= 0 else "?")
        print(f"{P:>3} {res['ncols_L']:>7} {res['nrows_L']:>7} | "
              f"{res['TRUE_laplace_only']['nullity']:>8} {obs:>9} {pred:>8} | "
              f"{res['RANDOM_laplace_only']['nullity']:>8} {res['RANDOM_laplace_plus_frobenius']['nullity']:>9} | "
              f"resid={res['true_laplace_residual_rows']} frobviol={res['true_frobenius_violations']} {flag} [{time.time()-t0:.1f}s]")
        if args.json_out:
            json.dump(allres, open(args.json_out, "w"), indent=1)

if __name__ == "__main__":
    main()
