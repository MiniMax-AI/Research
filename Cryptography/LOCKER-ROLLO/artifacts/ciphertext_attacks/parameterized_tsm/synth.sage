#!/usr/bin/env sage
"""Synthetic ideal-LRPC instance generator.

Produces exactly the same wire encoding as the reference implementation for
*arbitrary* (m, n, r), so that `./tsm` can be exercised on parameter
regimes the shipped KEMSET table does not contain (in particular gcd(m,n) > 1,
where R = GF(2^m)[X]/(f) is no longer a field).

Model (identical to the KEM):
    R  = GF(2^m)[X]/(f),  deg f = n,
    h  = pk                      (a public ring element)
    e1 = sum_k E_k a_k,  e2 = sum_k E_k b_k,  E = <E_1..E_r>_GF(2) subset GF(2^m),
    a_k, b_k in GF(2)[X]/(f)     (binary polynomials),
    u  = e1 + e2 h               (the encapsulation syndrome)

usage: synth.sage --m 71 --n 83 --r 5 [--seed 1] [--out inst.json]
"""
import argparse, hashlib, json, sys
from sage.all import GF, PolynomialRing, ZZ, gcd, matrix, set_random_seed


def poly_exponents(p):
    return sorted([int(e) for e in p.exponents()], reverse=True)


def build(m, n, seed):
    set_random_seed(seed)
    F2 = GF(2)
    P = PolynomialRing(F2, "x")
    x = P.gen()
    g = P.irreducible_element(m)          # GF(2^m) modulus (T-direction)
    f = P.irreducible_element(n)          # ring modulus   (X-direction)
    L = GF(2 ** n, "a", modulus=f)        # GF(2^n) = GF(2)[X]/(f)
    K = GF(2 ** m, "t", modulus=g)        # GF(2^m) = GF(2)[T]/(g)
    return F2, P, x, g, f, L, K


def encode(vec, m, n):
    """Reference ffi_vec encoding: n coefficients, ceil(m/8) bytes each, LE."""
    width = (m + 7) // 8
    out = bytearray()
    for c in vec:
        v = int(c.to_integer()) if hasattr(c, "to_integer") else int(c)
        out += int(v).to_bytes(width, "little")
    assert len(out) == n * width
    return bytes(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    m, n, r = args.m, args.n, args.r
    F2, P, x, g, f, L, K = build(m, n, args.seed)

    # R = GF(2^m)[X]/(f) -- a field iff gcd(m,n) = 1
    R = PolynomialRing(K, "X").quotient(f, "Xb")

    def to_R(coeffs):                     # coeffs: list of n elements of K
        return R(coeffs)

    def random_R():
        return R([K.random_element() for _ in range(n)])

    def coeffs_of(elem):
        c = list(elem.lift().list())
        return c + [K(0)] * (n - len(c))

    # error support E (r-dimensional over GF(2) inside GF(2^m))
    while True:
        E = [K.random_element() for _ in range(r)]
        Mw = matrix(GF(2), r, m, lambda i, j: (E[i].to_integer() >> j) & 1)
        if Mw.rank() == r:
            break

    # binary coefficient polynomials a_k, b_k in GF(2)[X]/(f)
    def rand_bin():
        return [F2.random_element() for _ in range(n)]

    a = [rand_bin() for _ in range(r)]
    b = [rand_bin() for _ in range(r)]
    A = [to_R([K(0)] * n) for _ in range(r)]
    B = [to_R([K(0)] * n) for _ in range(r)]
    for k in range(r):
        A[k] = sum((E[k] * a[k][j]) * R.gen() ** j for j in range(n) if a[k][j])
        B[k] = sum((E[k] * b[k][j]) * R.gen() ** j for j in range(n) if b[k][j])
    e1 = sum(A); e2 = sum(B)
    h = random_R()
    u = e1 + e2 * h

    # canonical (RREF) basis of E, exactly as the KEM hashes it
    rows = [int(v.to_integer()) for v in E]
    rank = 0
    for column in range(m):
        piv = next((i for i in range(rank, len(rows)) if (rows[i] >> column) & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(rank + 1, len(rows)):
            if (rows[i] >> column) & 1:
                rows[i] ^= rows[rank]
        rank += 1
        if rank == len(rows):
            break
    rows = rows[:rank]
    for i in range(1, rank):
        pv = (rows[i] & -rows[i]).bit_length() - 1
        for j in range(i):
            if (rows[j] >> pv) & 1:
                rows[j] ^= rows[i]
    width = (m + 7) // 8
    mask = hashlib.sha512(b"".join(v.to_bytes(width, "little") for v in rows)).digest()
    msg = hashlib.sha512(bytes([0x5a]) * 40).digest()[:64]
    v_tail = bytes(a ^ b for a, b in zip(msg, mask))
    d_tail = hashlib.sha512(msg).digest()

    pk = encode(coeffs_of(h), m, n)
    ct = encode(coeffs_of(u), m, n) + v_tail + d_tail
    res = {
        "m": m, "n": n, "r": r, "gcd": int(gcd(m, n)),
        "ring": poly_exponents(f), "field": poly_exponents(g),
        "public_key": pk.hex(), "encapsulation": ct.hex(),
        "shared_secret": hashlib.sha512(msg + ct[:n * width] + v_tail).hexdigest(),
        "support_echelon": [int(v) for v in rows],
        "e2_is_zero": bool(e2 == 0),
        "truth": {
            "E": [int(v.to_integer()) for v in E],
            "a": [[int(v) for v in ak] for ak in a],
            "b": [[int(v) for v in bk] for bk in b],
        },
    }
    text = json.dumps(res, separators=(",", ":"))
    if args.out:
        open(args.out, "w").write(text)
    if not args.quiet:
        print(text)


if __name__ == "__main__":
    main()
