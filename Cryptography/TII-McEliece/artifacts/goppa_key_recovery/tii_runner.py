#!/usr/bin/env python3
"""
Unified TII Official Track-2 Key-Recovery Runner & Strict Verifier
Supports:
  - Route A: Shortened Primal Code + Multi-threaded Bit-Parallel ISD (LB3 / Stern4 / Stern6)
             + Public Weight-Form Quotient + Gram Dual Linear Completion
  - Route B: Shortened Dual Code Higher-Order Vanishing (HOV p=3) + HOVER Frobenius Single-Direction
             + PGL2-Free Flag Multiplier Completion
  - Strict Official TII Verification:
      1) Distinct support elements L in GF(2^m) (defined by the official polynomial f(x) in pk_McEliece_<sec>.txt)
      2) Square-free degree-t Goppa polynomial g(x) (and monic irreducible check) with g(L_i) != 0
      3) Reconstructed parity-check matrix H_cand = gen_parity_full(L, g, echelonForm=True) is
         100% BIT-FOR-BIT IDENTICAL to the official TII challenge matrix H_pub!
      4) Also executes TII's official `sage check_solution_track2.sage <sec> <sk_file>` to confirm return code 1.
"""
import atexit
import json
import math
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from math import comb

sys.dont_write_bytecode = True

from sage.all import GF, Matrix, PolynomialRing, matrix, vector

_CAND_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else ""
if not os.path.isdir(os.path.join(_CAND_DIR, "route_a")) and sys.argv:
    _CAND_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))
if not os.path.isdir(os.path.join(_CAND_DIR, "route_a")):
    _CAND_DIR = os.getcwd()
BASE_DIR = _CAND_DIR
TII_REPO = f"{BASE_DIR}/tii_repo"
WITNESS_DIR = os.environ.get("TII_OUTPUT_DIR", f"{BASE_DIR}/witnesses")


def _ensure_third_party():
    if not os.path.isfile(os.path.join(TII_REPO, "key_rec_instances")):
        fetcher = os.path.join(BASE_DIR, "..", "fetch_third_party.sh")
        subprocess.run(["bash", fetcher, "--ensure"], check=True)


def load_tii_instance_meta(sec):
    _ensure_third_party()
    toys = [
        {"n": 32, "k": 22, "t": 2, "m": 5, "bitcomplexity": 22.0},
        {"n": 28, "k": 18, "t": 2, "m": 5, "bitcomplexity": 39.0},
    ]
    for inst in toys:
        if int(round(inst["bitcomplexity"])) == sec:
            return inst
    with open(f"{TII_REPO}/key_rec_instances") as f:
        for line in f:
            if line.strip():
                inst = json.loads(line)
                if int(round(inst["bitcomplexity"])) == sec:
                    return inst
    raise ValueError(f"Unknown TII challenge security level: {sec}")


def load_tii_public_key(sec):
    meta = load_tii_instance_meta(sec)
    n, k, t, m = meta["n"], meta["k"], meta["t"], meta["m"]
    pk_path = f"{TII_REPO}/public_keyRec/pk_McEliece_{sec}.txt"
    rows = []
    flist = None
    with open(pk_path, "r") as f:
        for i, line in enumerate(f):
            if i < n - k:
                line = line.strip("[").strip("]\n")
                rows.append([int(el) for el in line.split(" ") if el != ""])
            else:
                line = line.strip("[").strip("]\n")
                flist = [int(el) for el in line.split(", ")]
                break
    H_pub = matrix(GF(2), n - k, n, rows)
    GFx = PolynomialRing(GF(2), "x")
    modulus_poly = GFx(flist)
    ff = GF(2**m, name="a", modulus=modulus_poly)
    return meta, H_pub, flist, ff


def tii_gen_parity_full(L, g, m, n, ff, echelonForm=True):
    """Exact replica of McElieceKeyGen.sage:gen_parity_full."""
    t = g.degree()
    V = matrix.vandermonde(L).transpose()
    V1 = V[0:t]
    Gdiag = matrix.diagonal([ff.one() / g(L[i]) for i in range(n)])
    Hbar = V1 * Gdiag
    _, _, to_V = ff.vector_space(GF(2), map=True)
    H = matrix(GF(2), m * t, n)
    for i in range(t):
        for j in range(n):
            tmp_vec = to_V(Hbar[i, j])
            for b in range(m):
                H[m * i + b, j] = tmp_vec[b]
    if echelonForm:
        H = H.echelon_form()
    return H


def save_and_verify_tii_solution(sec, meta, H_pub, flist, ff, L_elems, g_poly, elapsed, route_name, extra_info=None):
    n, k, t, m = meta["n"], meta["k"], meta["t"], meta["m"]
    # Normalize g_poly to monic
    if g_poly.degree() != t:
        raise RuntimeError(f"Invalid degree of g: {g_poly.degree()} != {t}")
    g_poly = g_poly / g_poly.leading_coefficient()

    # Check distinctness and non-vanishing
    if len(set(L_elems)) != n:
        raise RuntimeError(f"Support elements not distinct: {len(set(L_elems))} != {n}")
    if g_poly.gcd(g_poly.derivative()) != 1:
        raise RuntimeError("Goppa polynomial is not square-free")
    if any(g_poly(x) == 0 for x in L_elems):
        raise RuntimeError("Goppa polynomial vanishes on support")

    # Reconstruct RREF parity-check matrix and check bit-for-bit equality with H_pub
    H_cand = tii_gen_parity_full(L_elems, g_poly, m, n, ff, echelonForm=True)
    bit_equal = (H_cand == H_pub)
    if not bit_equal:
        raise RuntimeError("H_cand != H_pub (RREF bit-for-bit check failed!)")

    # Write TII official Track-2 witness file:
    # Line 0: [g_0, g_1, ..., g_t] in variable 'a'
    # Line 1: [L_0, L_1, ..., L_{n-1}] in variable 'a'
    os.makedirs(WITNESS_DIR, exist_ok=True)
    sk_txt_path = f"{WITNESS_DIR}/sk_McEliece_{sec}.txt"
    g_list = [g_poly[i] for i in range(t + 1)]
    with open(sk_txt_path, "w") as f:
        f.write(str(g_list) + "\n")
        f.write(str(list(L_elems)) + "\n")

    # Run official TII check_solution_track2.sage:check_key
    sage_bin = os.environ.get("SAGE_BIN") or shutil.which("sage") or "sage"
    chk = subprocess.run(
        [sage_bin, "-c", f'load("check_solution_track2.sage"); print(check_key({sec}, "{sk_txt_path}"))'],
        cwd=TII_REPO,
        capture_output=True,
        text=True,
    )
    official_out = chk.stdout.strip().splitlines()[-1].strip() if chk.stdout.strip() else ""
    if official_out != "1":
        raise RuntimeError(f"Official check_solution_track2.sage failed: stdout={official_out!r}, stderr={chk.stderr!r}")

    rss_divisor = (1024.0 * 1024.0) if sys.platform == "darwin" else 1024.0
    rss_mb = max(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
    ) / rss_divisor

    phys_cores = 1 if route_name.startswith("Route-B") else 64
    core_sec = round(float(elapsed) * phys_cores, 2)
    core_hr = round(core_sec / 3600.0, 4)

    record = {
        "sec": sec,
        "n": n,
        "m": m,
        "t": t,
        "k": k,
        "mt": m * t,
        "bitcomplexity": meta["bitcomplexity"],
        "route": route_name,
        "physical_cores": phys_cores,
        "wall_seconds": round(elapsed, 3),
        "core_seconds": core_sec,
        "core_hours": core_hr,
        "peak_rss_mb": round(rss_mb, 1),
        "rref_bit_identical": bool(bit_equal),
        "official_checker_output": int(official_out),
        "g_is_irreducible": bool(g_poly.is_irreducible()),
        "witness_file": f"witnesses/sk_McEliece_{sec}.txt",
        "g_coeffs_int": [int(g_poly[i]._integer_representation()) for i in range(t + 1)],
        "support_first10_int": [int(L_elems[i]._integer_representation()) for i in range(min(10, n))],
        "extra": extra_info or {},
    }
    json_path = f"{WITNESS_DIR}/verified_tii_{sec}.json"
    with open(json_path, "w") as f:
        json.dump(record, f, indent=2)
    print(f"[SUCCESS] TII-{sec} ({n},{m},{t}) SOLVED & VERIFIED in {elapsed:.2f}s ({core_hr:.4f} core-h on {phys_cores}c) | RSS={rss_mb:.1f} MB | RREF_EQ=True | TII_CHECK=1", flush=True)
    return record

# ===========================================================================
# Route A: Shortened Primal Code + Multi-Threaded Bit-Parallel ISD
# ===========================================================================
sys.path.insert(0, f"{BASE_DIR}/route_a")
import goppa_attack

MINWORDS_MT_BIN = os.environ.get("MINWORDS_MT_BIN", f"{BASE_DIR}/route_a/minwords_mt")
RANKMASK_BIN = os.environ.get("RANKMASK_BIN", f"{BASE_DIR}/route_a/rankmask")
_SCRATCH_BIN_DIR = None


def _compile_helper_into_scratch(src, out, is_cxx=False):
    env_var = "CXX" if is_cxx else "CC"
    default_bin = "g++" if is_cxx else "gcc"
    stem = "g++" if is_cxx else "gcc"
    std_flags = ["-std=c++17"] if is_cxx else []
    env_comp = os.environ.get(env_var, "").strip()
    if env_comp:
        cands = [env_comp]
    else:
        cands = [default_bin]
        for prefix in ("", "/opt/homebrew/bin/", "/usr/local/bin/"):
            for ver in (16, 15, 14, 13, 12, 11):
                cand = f"{prefix}{stem}-{ver}"
                if shutil.which(cand):
                    cands.append(cand)
    omp_flag_sets = [
        ["-O3", *std_flags, "-Wall", "-march=native", "-mpopcnt", "-fopenmp"],
        ["-O3", *std_flags, "-Wall", "-march=native", "-fopenmp"],
        ["-O3", *std_flags, "-Wall", "-mcpu=native", "-fopenmp"],
        ["-O2", *std_flags, "-Wall", "-fopenmp"],
    ]
    for comp in cands:
        if not shutil.which(comp):
            continue
        for flags in omp_flag_sets:
            if subprocess.run([comp, *flags, src, "-o", out], capture_output=True).returncode == 0:
                return
    fb = env_comp or default_bin
    for flags in (
        ["-O3", *std_flags, "-Wall", "-march=native"],
        ["-O3", *std_flags, "-Wall", "-mcpu=native"],
        ["-O2", *std_flags, "-Wall"],
    ):
        proc = subprocess.run([fb, *flags, src, "-o", out], capture_output=True, text=True)
        if proc.returncode == 0:
            return
    raise RuntimeError(f"Failed to compile {src} with {fb}: {proc.stderr}")


def _ensure_route_a_binaries():
    global MINWORDS_MT_BIN, RANKMASK_BIN, _SCRATCH_BIN_DIR
    if os.path.isfile(MINWORDS_MT_BIN) and os.path.isfile(RANKMASK_BIN):
        return
    if _SCRATCH_BIN_DIR is None:
        _SCRATCH_BIN_DIR = tempfile.mkdtemp(prefix="tii-bins-")
        atexit.register(shutil.rmtree, _SCRATCH_BIN_DIR, ignore_errors=True)
    if not os.path.isfile(MINWORDS_MT_BIN):
        mw_out = os.path.join(_SCRATCH_BIN_DIR, "minwords_mt")
        _compile_helper_into_scratch(f"{BASE_DIR}/route_a/minwords_mt.c", mw_out, is_cxx=False)
        MINWORDS_MT_BIN = mw_out
    if not os.path.isfile(RANKMASK_BIN):
        rm_out = os.path.join(_SCRATCH_BIN_DIR, "rankmask")
        _compile_helper_into_scratch(f"{BASE_DIR}/route_a/rankmask.cpp", rm_out, is_cxx=True)
        RANKMASK_BIN = rm_out


def sample_words_mt(HA, target, weight, mode=0, seed=1, max_seconds=600.0, threads=64):
    """Invoke the OpenMP lock-free bit-parallel ISD sampler."""
    import tempfile
    rows, cols = HA.nrows(), HA.ncols()
    fd, path = tempfile.mkstemp(prefix="tii-ha-", suffix=".txt")
    os.close(fd)
    try:
        goppa_attack._write_parity(HA, path)
        cmd = [
            MINWORDS_MT_BIN,
            path,
            str(rows),
            str(cols),
            str(target),
            str(weight),
            str(mode),
            str(seed),
            repr(float(max_seconds)),
            str(threads),
        ]
        done = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        if os.path.exists(path):
            os.unlink(path)
    words = []
    seen = set()
    for line in done.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != weight:
            continue
        sup = tuple(int(x) for x in parts)
        if sup not in seen:
            seen.add(sup)
            words.append(sup)
    return words


def optimal_shortening_sizes(n, m, t):
    """Compute optimal shortening window ladder [n'_1, n'_2, ..., n] for (n, m, t)."""
    mt = m * t
    w = 2 * t + 1
    q = (2 * m + 1) if t == 5 else (m + 1)
    cands = []
    margin_ladder = (2.1, 2.6, 3.2) if mt >= 100 else (1.42, 1.65, 2.0, 2.6)
    for margin in margin_ladder:
        for np_val in range(mt + 2 * m + 6, n + 1):
            d = np_val - mt - 1
            F_d = 1 + d + d * (d - 1) // 2
            w_need = max(1, F_d - q)
            supply = comb(np_val, w) / (2**mt)
            if supply >= margin * w_need and np_val >= 2 * mt + 2 * m - 4:
                if np_val not in cands:
                    cands.append(np_val)
                break
    if n not in cands:
        cands.append(n)
    return cands


def solve_route_a(sec, threads=64, max_isd_seconds=600.0, mode=0):
    _ensure_route_a_binaries()
    t0 = time.time()
    meta, H_pub, flist, ff = load_tii_public_key(sec)
    n, k, t, m = meta["n"], meta["k"], meta["t"], meta["m"]
    mt = m * t
    assert H_pub.nrows() == mt and H_pub.ncols() == n

    sizes = optimal_shortening_sizes(n, m, t)
    print(f"[Route A] Starting TII-{sec} (n={n}, m={m}, t={t}, k={k}, mt={mt}) | sizes={sizes} | threads={threads}", flush=True)

    import random
    rng = random.Random(20261103 ^ sec)

    for nprime in sizes:
        margins = (280, 600, 1200) if nprime < n else (16, 64, 280, 600)
        for margin in margins:
            for attempt_idx in range(2 if nprime < n else 1):
                if nprime < n:
                    # Deterministic full-rank shortening selection
                    A = None
                    for _ in range(100):
                        cand_A = sorted(rng.sample(range(n), nprime))
                        if H_pub.matrix_from_columns(cand_A).rank() == mt:
                            A = cand_A
                            break
                    if A is None:
                        continue
                else:
                    A = list(range(n))

                HA = H_pub.matrix_from_columns(A)
                if HA.rank() != mt:
                    continue
                kp = nprime - mt
                aff = kp - 1
                feats = 1 + aff + aff * (aff - 1) // 2
                q_est = (2 * m + 1) if t == 5 else (m + 1)
                need = feats - q_est + margin
                t_isd0 = time.time()
                words = sample_words_mt(
                    HA,
                    need,
                    weight=2 * t + 1,
                    mode=mode,
                    seed=rng.randrange(1, 1 << 40),
                    max_seconds=max_isd_seconds,
                    threads=threads,
                )
                # If within reach (>= 25% collected), union additional runs until saturation!
                if feats * 25 // 100 <= len(words) < feats - q_est + 12:
                    word_set = set(words)
                    for extra_run in range(10):
                        if len(word_set) >= feats - q_est + 16:
                            break
                        more = sample_words_mt(
                            HA,
                            need,
                            weight=2 * t + 1,
                            mode=mode,
                            seed=rng.randrange(1, 1 << 40),
                            max_seconds=max_isd_seconds,
                            threads=threads,
                        )
                        word_set.update(more)
                        print(f"  [n'={nprime}] union run {extra_run+1}: total distinct words = {len(word_set)}/{feats - q_est}", flush=True)
                    words = list(word_set)
                t_isd = time.time() - t_isd0
                print(f"  [n'={nprime}] sampled {len(words)}/{need} weight-{2*t+1} words in {t_isd:.2f}s", flush=True)
                min_words_needed = (feats - q_est - 22) if nprime == n else (feats - q_est - 2)
                if len(words) < min_words_needed:
                    break  # nprime population too small, move to larger nprime

                try:
                    GA_ = HA.right_kernel().basis_matrix()
                    E, relations, aff_dim, pair_index, _ = goppa_attack.affine_polars(GA_, words)
                    q_dim = relations.nrows()
                    print(f"  [n'={nprime}] relation quotient dim = {q_dim} (expected {q_est})", flush=True)
                    if not (m <= q_dim <= 80):
                        continue

                    G0 = E.transpose() * GA_
                    mats = goppa_attack.polar_forms(relations, aff_dim, pair_index)
                    t_form = goppa_attack.chart_weight_form(G0)
                    mask = goppa_attack.express_in_basis(mats, t_form)
                    divided = None
                    if mask:
                        drop = (mask & -mask).bit_length() - 1
                        divided = mats[:drop] + mats[drop + 1:]

                    components = None
                    if divided is not None and len(divided) == m:
                        components = [divided]
                        print(f"  [n'={nprime}] quotient mod T_pub has exact dim {m} -> direct lift!", flush=True)
                    elif divided is not None:
                        try:
                            components, nlow, cut = goppa_attack.low_rank_components(
                                divided, aff_dim, nprime=nprime, m=m, finder=RANKMASK_BIN,
                                restricted=(len(divided) > 16), target=m, log=lambda s: print("   ", s, flush=True)
                            )
                        except Exception:
                            components = None
                    if components is None:
                        if len(mats) == m + 1:
                            # Even t or T_pub inside mats: low_rank_components on mats (size m+1) finds the unique high-rank direction!
                            components, nlow, cut = goppa_attack.low_rank_components(
                                mats, aff_dim, nprime=nprime, m=m, finder=RANKMASK_BIN,
                                restricted=False, target=m + 1, log=lambda s: print("   ", s, flush=True)
                            )
                        else:
                            components, nlow, cut = goppa_attack.low_rank_components(
                                mats, aff_dim, nprime=nprime, m=m, finder=RANKMASK_BIN,
                                restricted=(len(mats) > 17), target=m + 1, log=lambda s: print("   ", s, flush=True)
                            )

                    pivots = G0.pivots()
                    U = G0[:, pivots].inverse()
                    Gs = U * G0
                    pivot_set = set(pivots)
                    nonpivots = [j for j in range(nprime) if j not in pivot_set]

                    for which, raw in enumerate(components):
                        Bs = [U * B * U.transpose() for B in raw]
                        try:
                            W, blocks, I = goppa_attack.lift_all(Bs, Gs, nonpivots, pivots, m=m, log=lambda s: print("   ", s, flush=True))
                        except Exception as exc:
                            print(f"   lift_all failed: {exc!r}", flush=True)
                            continue

                        projections = range(1 << m) if len(Bs) == m + 1 else range(1)
                        for proj in projections:
                            if len(Bs) == m + 1:
                                Wp = matrix(GF(2), [W.row(i) + (((proj >> i) & 1) * W.row(m)) for i in range(m)])
                                bp = [[[int(blocks[tp][i][j]) ^ ((proj >> tp & 1) * int(blocks[m][i][j]))
                                        for j in range(m)] for i in range(m)] for tp in range(m)]
                            else:
                                Wp, bp = W, blocks[:m]
                            try:
                                labels_A = goppa_attack.make_field(Wp, bp, I, ff, m=m)
                            except Exception:
                                continue
                            if len(set(labels_A)) != nprime:
                                continue
                            poly = goppa_attack.goppa_from_word(ff, labels_A, words[0], t=t)
                            if poly.degree() != t or poly.gcd(poly.derivative()) != 1:
                                continue
                            if any(poly(x) == 0 for x in labels_A):
                                continue
                            if goppa_attack.parity_matrix(ff, labels_A, poly, m, t).row_space() != HA.row_space():
                                continue
                            labels_full = goppa_attack.extend_support(H_pub, ff, A, labels_A, poly, m=m, t=t)
                            elapsed = time.time() - t0
                            return save_and_verify_tii_solution(
                                sec, meta, H_pub, flist, ff, labels_full, poly, elapsed,
                                route_name=f"Route-A (n'={nprime}, ISD={t_isd:.2f}s, threads={threads})",
                                extra_info={"nprime": nprime, "words_sampled": len(words), "isd_seconds": round(t_isd, 3), "quotient_dim": q_dim},
                            )
                except Exception as exc:
                    print(f"  [n'={nprime}] attempt failed: {exc!r}", flush=True)
                    continue
    raise RuntimeError(f"Route A failed to solve TII-{sec}")


# ===========================================================================
# Route B: Shortened Dual Code HOV (p=3) + HOVER Frobenius Single-Direction
# ===========================================================================
def solve_route_b(sec):
    sys.path.insert(0, f"{BASE_DIR}/route_b")
    from hovlib import hov_params
    from hover_direction import directions_from_dual
    from hover_fibers import extract_fibers
    from hover_recover import conic_params, first_factors, multiplier_square, goppa_from_square, key_from_lh, complete_shortened
    import random as rnd

    t0 = time.time()
    meta, H_pub, flist, ff = load_tii_public_key(sec)
    n, k, t, m = meta["n"], meta["k"], meta["t"], meta["m"]
    e, sstar, kappa, a, b = hov_params(m, t, 3)
    print(f"[Route B HOVER] Starting TII-{sec} (n={n}, m={m}, t={t}, k={k}) | p=3, s*={sstar}, kappa={kappa}", flush=True)

    rng = rnd.Random(20261103 ^ sec)
    for attempt in range(4):
        short_cols = list(range(sstar)) if attempt == 0 else sorted(rng.sample(range(n), sstar))
        st = directions_from_dual(H_pub, m, t, n, ff, short_cols=short_cols, rng=rng, verbose=True)
        if st is None:
            continue
        cols = st["cats"][0].column_space().basis_matrix() * st["M"].change_ring(ff)
        FF, fib = extract_fibers(cols, ff, verbose=True)
        if FF is None:
            continue
        lams = first_factors(FF, fib, ff)
        params = conic_params(lams, ff, rng=rng, verbose=True)
        if params is None:
            continue
        ms = None
        for retry in range(4):
            ms = multiplier_square(H_pub, params, st["keep"], ff, t, rng=rnd.Random(97531 + 7919 * retry), verbose=True)
            if ms is not None:
                break
        if ms is None:
            continue
        S, beta, w2 = ms
        fit = goppa_from_square(beta, w2, ff, t, verbose=True)
        if fit is None:
            continue
        l, h = fit
        cands = key_from_lh(l, h, params, ff, t, verbose=True)
        for (ghat, alpha_keep, b_val, d_val) in cands:
            alpha_full = complete_shortened(H_pub, ghat, alpha_keep, st["keep"], ff, verbose=True)
            if alpha_full is None:
                continue
            try:
                elapsed = time.time() - t0
                return save_and_verify_tii_solution(
                    sec, meta, H_pub, flist, ff, alpha_full, ghat, elapsed,
                    route_name=f"Route-B (HOVER p=3, s*={sstar}, kappa={kappa})",
                    extra_info={"sstar": sstar, "kappa": kappa},
                )
            except Exception as exc:
                print(f"  Candidate verification failed: {exc!r}", flush=True)
                continue
    raise RuntimeError(f"Route B failed to solve TII-{sec}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: $(sage -c 'import sys; print(sys.executable)') tii_runner.py <route: A|B> <sec1,sec2,...> [threads] [mode] [max_isd_sec]")
        sys.exit(1)
    route = sys.argv[1].upper()
    secs = [int(x) for x in sys.argv[2].split(",")]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    mode = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    max_isd_sec = float(sys.argv[5]) if len(sys.argv) > 5 else 600.0

    for sec in secs:
        if route == "A":
            solve_route_a(sec, threads=threads, max_isd_seconds=max_isd_sec, mode=mode)
        elif route == "B":
            solve_route_b(sec)
        else:
            raise ValueError(f"Unknown route {route}")
