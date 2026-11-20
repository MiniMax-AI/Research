#!/usr/bin/env python3
"""Recover the real KEM shared secret from a solved Pluecker vector."""

import argparse
import hashlib
import itertools
import json
import struct

from sage.all import GF, Matrix, PolynomialRing, vector


M = 71
N = 83
R = 5
PAIRS = 9


def fields():
    binary = PolynomialRing(GF(2), "x0")
    x0 = binary.gen()
    L = GF(2**N, "x", modulus=x0**N + x0**7 + x0**4 + x0**2 + 1)
    polynomials = PolynomialRing(L, "T")
    T = polynomials.gen()
    Q = polynomials.quotient(T**M + T**6 + 1, "z")
    return L, Q


def parse_coefficients(value):
    raw = bytes.fromhex(value.removeprefix("0x"))
    return [int.from_bytes(raw[9*i:9*i+9], "little") for i in range(N)]


def transpose_to_q(values, L, Q):
    return Q([L.from_integer(sum(((value >> bit) & 1) << j
                                 for j, value in enumerate(values)))
              for bit in range(M)])


def rho(value, Q):
    return Q([coefficient * coefficient for coefficient in value.lift().list()])


def transformed_columns(value, L, Q):
    columns = [transpose_to_q(parse_coefficients(value), L, Q)]
    for _ in range(1, PAIRS):
        columns.append(rho(columns[-1], Q))
    return columns


def coordinate_keys():
    keys = []
    for pairs in itertools.combinations(range(PAIRS), R):
        for mask in range(1 << R):
            keys.append(tuple((pair, (mask >> k) & 1)
                              for k, pair in enumerate(pairs)))
    return keys


def key_for(selection):
    return tuple(sorted(selection.items()))


def read_pluecker(path, L):
    raw = open(path, "rb").read()
    expected = len(coordinate_keys()) * 16
    if len(raw) != expected:
        raise ValueError(f"bad Pluecker vector length {len(raw)} (expected {expected})")
    return [L.from_integer(lo | (hi << 64))
            for lo, hi in struct.iter_unpack("<QQ", raw)]


def reconstruct_generator(values, L):
    keys = coordinate_keys()
    lookup = {key: values[index] for index, key in enumerate(keys)}
    chart = {pair: 1 for pair in range(R)}  # B_0,...,B_4
    chart_minor = lookup[key_for(chart)]
    if not chart_minor:
        raise ValueError("the fixed B_0,...,B_4 chart is singular")

    # Columns are A_0,B_0,A_1,B_1,... .  Normalize the B chart to I_5.
    G = Matrix(L, R, 2 * PAIRS)
    for i in range(R):
        G[i, 2*i + 1] = 1

    # Every column belonging to an external pair is obtained by replacing
    # one chart column in a maximal minor.
    for pair in range(R, PAIRS):
        for side in (0, 1):
            for row in range(R):
                selection = dict(chart)
                del selection[row]
                selection[pair] = side
                G[row, 2*pair + side] = lookup[key_for(selection)] / chart_minor

    # The diagonal entry of A_i is another one-column replacement.
    for i in range(R):
        selection = dict(chart)
        selection[i] = 0
        G[i, 2*i] = lookup[key_for(selection)] / chart_minor

    # For an off-diagonal A_i[k], combine A_i with a known external column v.
    # The corresponding normalized minor is A_i[i] v[k] + A_i[k] v[i].
    for i in range(R):
        candidates = [(pair, side) for pair in range(R, PAIRS)
                      for side in (0, 1) if G[i, 2*pair + side]]
        if not candidates:
            raise ValueError("could not find a nonsingular alternate chart")
        pair, side = candidates[0]
        v = G.column(2*pair + side)
        for row in range(R):
            if row == i:
                continue
            selection = {j: 1 for j in range(R) if j not in (i, row)}
            selection[i] = 0
            selection[pair] = side
            determinant = lookup[key_for(selection)] / chart_minor
            G[row, 2*i] = (determinant + G[i, 2*i] * v[row]) / v[i]

    # This catches chart/orientation mistakes before any hash-derived answer.
    for key, claimed in lookup.items():
        columns = [2*pair + side for pair, side in key]
        actual = G.matrix_from_columns(columns).det() * chart_minor
        if actual != claimed:
            raise ValueError("reconstructed generator fails a Pluecker minor")
    return G


def solve_error_transform(instance, G, L, Q):
    H = transformed_columns(instance["public_key"], L, Q)
    U = transformed_columns(instance["encapsulation"], L, Q)
    equations = Matrix(Q, PAIRS, R,
                       [Q(G[k, 2*q]) + H[q] * Q(G[k, 2*q + 1])
                        for q in range(PAIRS) for k in range(R)])
    coefficients = None
    for rows in itertools.combinations(range(PAIRS), R):
        square_system = equations.matrix_from_rows(rows)
        if square_system.det():
            coefficients = square_system.solve_right(vector(Q, [U[q] for q in rows]))
            break
    if coefficients is None or equations * coefficients != vector(Q, U):
        raise ValueError("could not solve the transformed error equations")
    B = sum((coefficients[k] * Q(G[k, 1]) for k in range(R)), Q.zero())
    # In this chart B_0 is enough; its Frobenius orbit must reproduce all B_q.
    current = B
    for q in range(PAIRS):
        expected = sum((coefficients[k] * Q(G[k, 2*q + 1]) for k in range(R)),
                       Q.zero())
        if current != expected:
            raise ValueError("recovered B columns fail the Frobenius check")
        current = rho(current, Q)
    return B


def canonical_binary_basis(elements):
    rows = elements[:]
    rank = 0
    for column in range(M):
        pivot = next((i for i in range(rank, len(rows))
                      if (rows[i] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank + 1, len(rows)):
            if (rows[i] >> column) & 1:
                rows[i] ^= rows[rank]
        rank += 1
        if rank == len(rows):
            break
    rows = rows[:rank]
    for i in range(1, rank):
        pivot = (rows[i] & -rows[i]).bit_length() - 1
        for j in range(i):
            if (rows[j] >> pivot) & 1:
                rows[j] ^= rows[i]
    return rows


def support_from_transpose(B, L):
    coefficients = list(B.lift().list())
    coefficients += [L.zero()] * (M - len(coefficients))
    integers = [int(value.to_integer()) for value in coefficients]
    error = [sum(((integers[bit] >> j) & 1) << bit for bit in range(M))
             for j in range(N)]
    support = canonical_binary_basis(error)
    if len(support) != R:
        raise ValueError(f"recovered error has rank {len(support)}, expected {R}")
    return support


def recover_secret(instance, support):
    ciphertext = bytes.fromhex(instance["encapsulation"].removeprefix("0x"))
    u_end = N * 9
    if len(ciphertext) != u_end + 128:
        raise ValueError("unexpected ciphertext length")
    u = ciphertext[:u_end]
    v = ciphertext[u_end:u_end + 64]
    digest = ciphertext[u_end + 64:]
    encoded_support = b"".join(value.to_bytes(9, "little") for value in support)
    mask = hashlib.sha512(encoded_support).digest()
    message = bytes(a ^ b for a, b in zip(v, mask))
    if hashlib.sha512(message).digest() != digest:
        raise ValueError("recovered support does not satisfy the ciphertext digest")
    return hashlib.sha512(message + u + v).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("instance")
    parser.add_argument("pluecker")
    parser.add_argument("--show-support", action="store_true")
    args = parser.parse_args()
    instance = json.load(open(args.instance))
    instance = instance.get("result", instance)
    L, Q = fields()
    values = read_pluecker(args.pluecker, L)
    generator = reconstruct_generator(values, L)
    transformed_error = solve_error_transform(instance, generator, L, Q)
    support = support_from_transpose(transformed_error, L)
    result = {"real_shared_secret": recover_secret(instance, support)}
    if args.show_support:
        result["support"] = [hex(value) for value in support]
    print(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
