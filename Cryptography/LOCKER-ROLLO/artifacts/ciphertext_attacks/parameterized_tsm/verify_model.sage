#!/usr/bin/env python3
"""Verify the TSM model against ground truth (independent re-derivation).

Ground truth comes from ./gen<SET>, which prints the internal values of
the *unmodified* reference implementation: E, e1, e2, pk, ct, ss, m.

Checks performed
  1. ideal-LRPC relation            u = e1 + e2*h   inside R = GF(2^m)[X]/(f)
  2. transposition                  Phi : R -> Q = GF(2^n)[T]/(g)  is a ring
                                    homomorphism on these elements
  3. support decomposition          e_i = sum_k E_k c_k,  c_k in GF(2)^n
  4. the P Frobenius equations      U_q = sum_k E~_k (A_{q,k} + H_q B_{q,k})
  5. the TRUE Pluecker vector of G, its chart minor, and (optionally) equality
     up to scale with the vector produced by ./tsm

usage: verify_model.sage GT.txt --m .. --n .. --r .. --pairs .. --ring .. --field ..
                          [--pluecker FILE]
"""
from __future__ import annotations
import argparse, os, struct, sys
from sage.all import GF, Matrix, PolynomialRing, vector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tsm_recover import coordinate_keys, key_for


def load_gt(path):
    d = {}
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            d[k.strip()] = v.strip().removeprefix("0x")
    return d


def coeffs(hx, count, mbytes):
    raw = bytes.fromhex(hx)
    return [int.from_bytes(raw[mbytes * i:mbytes * i + mbytes], "little") for i in range(count)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gt")
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--pairs", type=int, required=True)
    ap.add_argument("--ring", required=True)
    ap.add_argument("--field", required=True)
    ap.add_argument("--pluecker")
    args = ap.parse_args()
    ring = [int(x) for x in args.ring.split(",")]
    field = [int(x) for x in args.field.split(",")]
    m, n, r, P = args.m, args.n, args.r, args.pairs
    mbytes = (m + 7) // 8
    gt = load_gt(args.gt)

    F2 = PolynomialRing(GF(2), "z")
    z = F2.gen()
    fmod2 = z ** n
    for e in ring[1:]:
        fmod2 += z ** e
    gmod2 = z ** m
    for e in field[1:]:
        gmod2 += z ** e
    if not (fmod2.is_irreducible() and gmod2.is_irreducible()):
        print("!! a modulus is reducible: the transposition is not a field embedding")
    L = GF(2 ** n, "a", modulus=fmod2)        # GF(2^n): the "horizontal" field
    Fm = GF(2 ** m, "b", modulus=gmod2)       # GF(2^m): the "vertical" field

    RX = PolynomialRing(Fm, "X"); X = RX.gen()
    R = RX.quotient(RX(fmod2), "Xb")          # R = GF(2^m)[X]/(f)
    PR = PolynomialRing(L, "T")
    Q = PR.quotient(PR(gmod2), "q")           # Q = GF(2^n)[T]/(g)

    E = [Fm.from_integer(v) for v in coeffs(gt["E_raw"], r, mbytes)]
    e1 = R([Fm.from_integer(v) for v in coeffs(gt["e1"], n, mbytes)])
    e2 = R([Fm.from_integer(v) for v in coeffs(gt["e2"], n, mbytes)])
    h = R([Fm.from_integer(v) for v in coeffs(gt["pk"], n, mbytes)])
    u = R([Fm.from_integer(v) for v in coeffs(gt["ct"][:2 * n * mbytes], n, mbytes)])
    print("1. u == e1 + e2*h              :", e1 + e2 * h == u)
    print("   e2 != 0, h != 0             :", e2 != 0, h != 0)
    print("   dim_{GF(2)} supp(E)         :",
          Matrix(GF(2), [[(int(e.to_integer()) >> b) & 1 for b in range(m)] for e in E]).rank())

    # ---- transposition Phi -------------------------------------------------
    def phi(elt):
        cols = [int(c.to_integer()) for c in elt.lift().list()]
        cols += [0] * (n - len(cols))
        return Q([L.from_integer(sum(((v >> bit) & 1) << j for j, v in enumerate(cols)))
                  for bit in range(m)])
    print("2. Phi multiplicative          :", phi(e1 * h) == phi(e1) * phi(h),
          "| additive:", phi(e1 + e2) == phi(e1) + phi(e2))
    print("   Phi(u) == Phi(e1)+Phi(e2)Phi(h):", phi(u) == phi(e1) + phi(e2) * phi(h))

    # ---- support decomposition --------------------------------------------
    Emat = Matrix(GF(2), r, m, [[(int(e.to_integer()) >> b) & 1 for b in range(m)] for e in E])
    def decompose(elt):
        cols = elt.lift().list()
        cols += [Fm.zero()] * (n - len(cols))
        acc = [0] * r
        for j, c in enumerate(cols):
            v = vector(GF(2), [(int(c.to_integer()) >> b) & 1 for b in range(m)])
            if v.is_zero():
                continue
            sol = Emat.solve_left(v)
            for k in range(r):
                if sol[k]:
                    acc[k] |= 1 << j
        return [L.from_integer(x) for x in acc]
    a = decompose(e1)
    b = decompose(e2)
    def rebuild(cv):
        out = R.zero()
        for k in range(r):
            cs = [Fm.zero()] * n
            for j in range(n):
                if (int(cv[k].to_integer()) >> j) & 1:
                    cs[j] = E[k]
            out += R(cs)
        return out
    print("3. e1 = sum E_k a_k            :", rebuild(a) == e1,
          "| e2 = sum E_k b_k:", rebuild(b) == e2)
    moore = Matrix(L, r, r, [[b[k] ** (2 ** q) for q in range(r)] for k in range(r)]).det()
    print("   Moore det of (b_1..b_r)     :", "NONZERO" if moore else "ZERO")

    # ---- Frobenius equations ----------------------------------------------
    Et = [Q([L.one() if (int(e.to_integer()) >> bit) & 1 else L.zero() for bit in range(m)])
          for e in E]
    print("   Phi(E_k) == E~_k            :", all(phi(R([E[k]])) == Et[k] for k in range(r)))
    Hq, Uq = phi(h), phi(u)
    ok = True
    for q in range(P):
        lhs = sum((Et[k] * (Q(a[k] ** (2 ** q)) + Hq * Q(b[k] ** (2 ** q))) for k in range(r)), Q.zero())
        ok = ok and (lhs == Uq)
        Hq = Q([c * c for c in Hq.lift().list()])
        Uq = Q([c * c for c in Uq.lift().list()])
    print(f"4. all {P} Frobenius equations hold:", ok)

    # ---- true Pluecker vector ---------------------------------------------
    G = Matrix(L, r, 2 * P)
    for k in range(r):
        ca, cb = a[k], b[k]
        for q in range(P):
            G[k, 2 * q] = ca
            G[k, 2 * q + 1] = cb
            ca, cb = ca * ca, cb * cb
    keys = coordinate_keys(P, r)
    pl = {key_for(dict(key)): G.matrix_from_columns([2 * p + s for p, s in key]).det() for key in keys}
    chart = key_for({p: 1 for p in range(r)})
    nz = sum(1 for v in pl.values() if v)
    print(f"5. chart (B_0..B_{r-1}) minor    :", "NONZERO" if pl[chart] else "ZERO",
          f"| nonzero maximal minors: {nz}/{len(pl)}")

    if args.pluecker:
        raw = open(args.pluecker, "rb").read()
        ncols, rr, pp, nw = struct.unpack_from("<IIII", raw, 0)
        body = raw[16:]
        got = [L.from_integer(sum(struct.unpack_from("<4Q", body, i * 32)[k] << (64 * k)
                                 for k in range(nw))) for i in range(ncols)]
        scale, bad = None, 0
        for key, val in zip(keys, got):
            t = pl[key_for(dict(key))]
            if t and val:
                if scale is None:
                    scale = val / t
                elif val != scale * t:
                    bad += 1
            elif bool(t) != bool(val):
                bad += 1
        print("   tsm Pluecker == true (up to scale):", bad == 0, f"(mismatches={bad})")


if __name__ == "__main__":
    main()
