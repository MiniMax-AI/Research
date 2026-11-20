#!/usr/bin/env python3
"""Zero-dependency end-to-end verifier for LOCKER-I (m=71, n=83, d=7, r=5).
Compiles kemsolve.cpp with g++ only (no NTL/Sage/M4RI needed), generates 5 fresh
spec-authentic LOCKER-I keypairs & ciphertexts from /dev/urandom, and verifies 100% recovery.
"""
import hashlib, os, secrets, subprocess, tempfile, time
from pathlib import Path

M, N, D, R = 71, 83, 7, 5
MOD_T = (1 << 71) | (1 << 6) | 1

def gf2m_mul(a: int, b: int) -> int:
    res = 0
    while b:
        lsb = b & -b
        res ^= a << (lsb.bit_length() - 1)
        b ^= lsb
    while res >= (1 << 71):
        hi = res >> 71
        res = (res & ((1 << 71) - 1)) ^ hi ^ (hi << 6)
    return res

def gf2m_inv(a: int) -> int:
    u, v, g1, g2 = a, MOD_T, 1, 0
    while u != 1:
        j = u.bit_length() - v.bit_length()
        if j < 0:
            u, v, g1, g2, j = v, u, g2, g1, -j
        u ^= v << j
        g1 ^= g2 << j
    while g1 >= (1 << 71):
        hi = g1 >> 71
        g1 = (g1 & ((1 << 71) - 1)) ^ hi ^ (hi << 6)
    return g1

def poly_mul_modP(p1: list[int], p2: list[int]) -> list[int]:
    tmp = [0] * (2 * N - 1)
    for i, a in enumerate(p1):
        if a:
            for j, b in enumerate(p2):
                if b: tmp[i + j] ^= gf2m_mul(a, b)
    for deg in range(2 * N - 2, N - 1, -1):
        c = tmp[deg]
        if c:
            sh = deg - N
            tmp[sh] ^= c; tmp[sh + 2] ^= c; tmp[sh + 4] ^= c; tmp[sh + 7] ^= c
    return tmp[:N]

def poly_inv_modP(a: list[int]):
    def trim(p):
        while len(p) > 1 and p[-1] == 0: p.pop()
        return p
    def p_divmod(num, den):
        num, den = trim(num[:]), trim(den[:])
        if len(num) < len(den) or den == [0]: return [0], num
        q = [0] * (len(num) - len(den) + 1)
        inv_lc = gf2m_inv(den[-1])
        while len(num) >= len(den) and num != [0]:
            sh = len(num) - len(den)
            sc = gf2m_mul(num[-1], inv_lc)
            q[sh] ^= sc
            for i, c in enumerate(den):
                if c: num[i + sh] ^= gf2m_mul(sc, c)
            num = trim(num)
        return trim(q), num
    def p_mul(p1, p2):
        out = [0] * (len(p1) + len(p2) - 1)
        for i, x in enumerate(p1):
            if x:
                for j, y in enumerate(p2):
                    if y: out[i + j] ^= gf2m_mul(x, y)
        return trim(out)
    def p_add(p1, p2):
        return trim([(p1[i] if i < len(p1) else 0) ^ (p2[i] if i < len(p2) else 0) for i in range(max(len(p1), len(p2)))])

    old_r, r = [1, 0, 1, 0, 1, 0, 0, 1] + [0] * (N - 8) + [1], trim(a[:])
    old_t, t = [0], [1]
    while r != [0]:
        q, rem = p_divmod(old_r, r)
        old_r, r = r, rem
        old_t, t = t, p_add(old_t, p_mul(q, t))
    if len(old_r) != 1 or old_r[0] == 0: return None
    inv_c = gf2m_inv(old_r[0])
    return [gf2m_mul(c, inv_c) for c in old_t] + [0] * (N - len(old_t))

def rref_gf2(vecs: list[int]) -> list[int]:
    r, rk = vecs[:], 0
    for col in range(M):
        piv = next((i for i in range(rk, len(r)) if (r[i] >> col) & 1), -1)
        if piv < 0: continue
        r[rk], r[piv] = r[piv], r[rk]
        for i in range(len(r)):
            if i != rk and ((r[i] >> col) & 1): r[i] ^= r[rk]
        rk += 1
    return r[:rk]

def sample_subspace(dim: int) -> list[int]:
    while True:
        rb = rref_gf2([int.from_bytes(secrets.token_bytes(9), "little") & ((1 << M) - 1) for _ in range(dim)])
        if len(rb) == dim: return rb

def sample_vec(basis: list[int]) -> list[int]:
    dim = len(basis)
    while True:
        v = [0] * N
        for j in range(N):
            msk = secrets.randbelow(1 << dim)
            for k in range(dim):
                if (msk >> k) & 1: v[j] ^= basis[k]
        if len(rref_gf2(v)) == dim: return v

def encode_poly(p: list[int]) -> bytes:
    return b"".join(x.to_bytes(9, "little") for x in p)

def build_kemsolve(cpp: Path, bin_p: Path) -> None:
    cxx = os.environ.get("CXX", "g++")
    probe_src = "int main() { return 0; }\n"

    def accepts(flags: list[str]) -> bool:
        res = subprocess.run(
            [cxx, "-x", "c++", "-", "-o", os.devnull, *flags],
            input=probe_src, capture_output=True, text=True,
        )
        return res.returncode == 0

    omp_flags: list[str] = []
    if accepts(["-std=c++17", "-fopenmp"]):
        omp_flags = ["-fopenmp"]
    else:
        for pfx in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
            cand = ["-Xpreprocessor", "-fopenmp", f"-I{pfx}/include", f"-L{pfx}/lib", "-lomp"]
            if Path(pfx).is_dir() and accepts(["-std=c++17", *cand]):
                omp_flags = cand
                break

    arch_flags: list[str] = []
    if accepts(["-std=c++17", "-march=native"]):
        arch_flags.append("-march=native")
    if accepts(["-std=c++17", "-mpclmul"]):
        arch_flags.append("-mpclmul")

    subprocess.run(
        [cxx, "-O3", "-std=c++17", *arch_flags, *omp_flags, "-o", str(bin_p), str(cpp)],
        check=True,
    )


def main():
    here = Path(__file__).resolve().parent
    cpp = here / "kemsolve.cpp"
    n_trials = int(os.environ.get("TRIALS", "5"))
    with tempfile.TemporaryDirectory() as td:
        bin_p = Path(td) / "kemsolve"
        build_kemsolve(cpp, bin_p)
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "12"), BLK="16", MAXREL="256")
        for trial in range(n_trials):
            while True:
                F = sample_subspace(D)
                x = sample_vec(F)
                xi = poly_inv_modP(x)
                if xi is not None:
                    h = poly_mul_modP(xi, sample_vec(F))
                    break
            E = sample_subspace(R)
            e1, e2 = sample_vec(E), sample_vec(E)
            u = [a ^ b for a, b in zip(e1, poly_mul_modP(e2, h))]
            assert len(rref_gf2(u)) == M
            m_seed = secrets.token_bytes(64)
            u_b = encode_poly(u)
            v_b = bytes(a ^ b for a, b in zip(m_seed, hashlib.sha512(encode_poly(E)).digest()))
            d_b = hashlib.sha512(m_seed).digest()
            ss_true = hashlib.sha512(m_seed + u_b + v_b).hexdigest()
            t0 = time.monotonic()
            proc = subprocess.run([str(bin_p), encode_poly(h).hex(), (u_b + v_b + d_b).hex()], env=env, capture_output=True, text=True, check=True)
            dt = time.monotonic() - t0
            ss_rec = proc.stdout.strip().splitlines()[-1].strip()
            assert ss_rec == ss_true
            print(f"Trial {trial}: rank_GF2(u)={M}, kemsolve={dt:.4f}s, ss_match=True (ss={ss_rec[:24]}...)")
        print(f"ALL {n_trials} FRESH RANDOM LOCKER-I TRIALS VERIFIED 100%!")

if __name__ == "__main__":
    main()
