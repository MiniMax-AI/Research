#!/usr/bin/env python3
"""Parameterised recovery stage of the Transposed Support-Minor attack.

Input : instance JSON ({"public_key":hex,"encapsulation":hex[, "shared_secret":hex]}),
        Pluecker vector produced by ./tsm.
Output: JSON {"real_shared_secret":hex, "support":[...], "verified_digest":bool}

The recovered support is validated against d = SHA512(m) taken from the
ciphertext, so a successful run is a proof (up to 2^-512) that the rank-r
error support of the encapsulation was recovered.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, struct, sys

from sage.all import GF, Matrix, PolynomialRing, vector


def build_fields(n, ring_low, m, field_low):
    F2 = PolynomialRing(GF(2), "x0")
    x0 = F2.gen()
    ring_mod = x0 ** n
    for e in ring_low:
        ring_mod += x0 ** e
    L = GF(2 ** n, "x", modulus=ring_mod)
    PR = PolynomialRing(L, "T")
    T = PR.gen()
    field_mod = T ** m
    for e in field_low:
        field_mod += T ** e
    Q = PR.quotient(field_mod, "z")
    return L, Q


def parse_coeffs(value, n, mbytes):
    raw = bytes.fromhex(value.removeprefix("0x"))
    return [int.from_bytes(raw[mbytes * i:mbytes * i + mbytes], "little") for i in range(n)]


def transpose_to_q(values, L, Q, m):
    return Q([L.from_integer(sum(((v >> bit) & 1) << j for j, v in enumerate(values)))
              for bit in range(m)])


def rho(value, Q):
    return Q([c * c for c in value.lift().list()])


def transformed(value, L, Q, m, n, mbytes, pairs):
    cols = [transpose_to_q(parse_coeffs(value, n, mbytes), L, Q, m)]
    for _ in range(1, pairs):
        cols.append(rho(cols[-1], Q))
    return cols


def coordinate_keys(pairs, r):
    keys = []
    for combo in itertools.combinations(range(pairs), r):
        for mask in range(1 << r):
            keys.append(tuple((p, (mask >> k) & 1) for k, p in enumerate(combo)))
    return keys


def key_for(selection):
    return tuple(sorted(selection.items()))


def read_pluecker(path, L, nwords):
    raw = open(path, "rb").read()
    ncols, r, pairs, nw = struct.unpack_from("<IIII", raw, 0)
    body = raw[16:]
    if nw != nwords:
        raise ValueError(f"word-size mismatch {nw} != {nwords}")
    if len(body) != ncols * 4 * 8:
        raise ValueError("bad Pluecker length")
    vals = []
    for i in range(ncols):
        words = struct.unpack_from("<4Q", body, i * 32)
        vals.append(L.from_integer(sum(words[k] << (64 * k) for k in range(nwords))))
    return vals, r, pairs


def reconstruct_generator(values, L, r, pairs, verify=True):
    keys = coordinate_keys(pairs, r)
    if len(keys) != len(values):
        raise ValueError("Pluecker length does not match (pairs, r)")
    lookup = {key: values[i] for i, key in enumerate(keys)}
    chart = {p: 1 for p in range(r)}
    chart_minor = lookup[key_for(chart)]
    if not chart_minor:
        raise ValueError("the fixed B_0..B_{r-1} chart is singular")

    G = Matrix(L, r, 2 * pairs)
    for i in range(r):
        G[i, 2 * i + 1] = 1
    for pair in range(r, pairs):
        for side in (0, 1):
            for row in range(r):
                sel = dict(chart); del sel[row]; sel[pair] = side
                G[row, 2 * pair + side] = lookup[key_for(sel)] / chart_minor
    for i in range(r):
        sel = dict(chart); sel[i] = 0
        G[i, 2 * i] = lookup[key_for(sel)] / chart_minor
    for i in range(r):
        cands = [(p, s) for p in range(r, pairs) for s in (0, 1) if G[i, 2 * p + s]]
        if not cands:
            raise ValueError("no nonsingular alternate chart")
        pair, side = cands[0]
        v = G.column(2 * pair + side)
        for row in range(r):
            if row == i:
                continue
            sel = {j: 1 for j in range(r) if j not in (i, row)}
            sel[i] = 0; sel[pair] = side
            det = lookup[key_for(sel)] / chart_minor
            G[row, 2 * i] = (det + G[i, 2 * i] * v[row]) / v[i]
    if verify:
        for key, claimed in lookup.items():
            cols = [2 * p + s for p, s in key]
            if G.matrix_from_columns(cols).det() * chart_minor != claimed:
                raise ValueError("reconstructed generator fails a Pluecker minor")
    return G


def solve_error_transform(G, H, U, Q, r, pairs):
    L = Q.base_ring()
    z = Q.gen()
    m = Q.modulus().degree()

    def q_to_vec(elem):
        c = list(elem.lift().list())
        return c + [L.zero()] * (m - len(c))

    eqs = Matrix(Q, pairs, r,
                 [Q(G[k, 2 * q]) + H[q] * Q(G[k, 2 * q + 1])
                  for q in range(pairs) for k in range(r)])
    A_L = Matrix(L, pairs * m, r * m)
    for q in range(pairs):
        for k in range(r):
            cur_pow = eqs[q, k]
            for j in range(m):
                col_vec = q_to_vec(cur_pow)
                for i in range(m):
                    A_L[q * m + i, k * m + j] = col_vec[i]
                if j + 1 < m:
                    cur_pow = cur_pow * z
    u_L = vector(L, [c for q in range(pairs) for c in q_to_vec(U[q])])
    try:
        c_L = A_L.solve_right(u_L)
    except ValueError as exc:
        raise ValueError("could not solve the transformed error equations") from exc
    coeffs = vector(Q, [Q([c_L[k * m + j] for j in range(m)]) for k in range(r)])
    if eqs * coeffs != vector(Q, U):
        raise ValueError("could not solve the transformed error equations")
    B = sum((coeffs[k] * Q(G[k, 1]) for k in range(r)), Q.zero())
    cur = B
    for q in range(pairs):
        expect = sum((coeffs[k] * Q(G[k, 2 * q + 1]) for k in range(r)), Q.zero())
        if cur != expect:
            raise ValueError("recovered B columns fail the Frobenius check")
        cur = rho(cur, Q)
    return B


def canonical_binary_basis(elements, m):
    rows = elements[:]
    rank = 0
    for col in range(m):
        piv = next((i for i in range(rank, len(rows)) if (rows[i] >> col) & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(rank + 1, len(rows)):
            if (rows[i] >> col) & 1:
                rows[i] ^= rows[rank]
        rank += 1
        if rank == len(rows):
            break
    rows = rows[:rank]
    for i in range(1, rank):
        piv = (rows[i] & -rows[i]).bit_length() - 1
        for j in range(i):
            if (rows[j] >> piv) & 1:
                rows[j] ^= rows[i]
    return rows


def support_from_transpose(B, L, m, n, r):
    coeffs = list(B.lift().list())
    coeffs += [L.zero()] * (m - len(coeffs))
    ints = [int(c.to_integer()) for c in coeffs]
    error = [sum(((ints[bit] >> j) & 1) << bit for bit in range(m)) for j in range(n)]
    supp = canonical_binary_basis(error, m)
    if len(supp) != r:
        raise ValueError(f"recovered error has rank {len(supp)}, expected {r}")
    return supp


def recover_secret(ct_bytes, supp, n, mbytes):
    u_end = n * mbytes
    u, v, d = ct_bytes[:u_end], ct_bytes[u_end:u_end + 64], ct_bytes[u_end + 64:u_end + 128]
    encoded = b"".join(x.to_bytes(mbytes, "little") for x in supp)
    mask = hashlib.sha512(encoded).digest()
    msg = bytes(a ^ b for a, b in zip(v, mask))
    verified = hashlib.sha512(msg).digest() == d
    return hashlib.sha512(msg + u + v).hexdigest(), verified


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("instance")
    ap.add_argument("pluecker")
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--pairs", type=int, required=True)
    ap.add_argument("--ring", required=True, help="n,e1,e2,... exponents of the ring modulus")
    ap.add_argument("--field", required=True, help="m,e1,e2,... exponents of the GF(2^m) modulus")
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args()

    ring = [int(x) for x in args.ring.split(",")]
    field = [int(x) for x in args.field.split(",")]
    n, m, r, pairs = args.n, args.m, args.r, args.pairs
    mbytes = (m + 7) // 8
    nwords = (n + 63) // 64
    L, Q = build_fields(n, ring[1:], m, field[1:])
    inst = json.load(open(args.instance))
    inst = inst.get("result", inst)
    values, r_file, pairs_file = read_pluecker(args.pluecker, L, nwords)
    if (r_file, pairs_file) != (r, pairs):
        raise ValueError("Pluecker header disagrees with the parameters")
    G = reconstruct_generator(values, L, r, pairs, verify=not args.no_verify)
    H = transformed(inst["public_key"], L, Q, m, n, mbytes, pairs)
    U = transformed(inst["encapsulation"], L, Q, m, n, mbytes, pairs)
    B = solve_error_transform(G, H, U, Q, r, pairs)
    supp = support_from_transpose(B, L, m, n, r)
    ss, verified = recover_secret(bytes.fromhex(inst["encapsulation"].removeprefix("0x")),
                                  supp, n, mbytes)
    out = {"real_shared_secret": ss, "verified_digest": bool(verified),
           "support": [hex(x) for x in supp]}
    shown = inst.get("shared_secret")
    if shown:
        out["matches_published"] = ss == shown.removeprefix("0x").lower()
    print(json.dumps(out, separators=(",", ":")))


if __name__ == "__main__":
    main()
