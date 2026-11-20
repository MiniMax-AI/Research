#!/usr/bin/env python3
"""Analytic design-point theory for the TSM linearisation.

The TSM system over GF(2^n) has

    ncols   = C(P,r) * 2^r                unknowns  (maximal minors of the r x 2P matrix G)
    nrows   = C(P-1,r) * m                equations (one Q-equation per (r+1)-subset
                                                     containing pair 0, expanded over GF(2^n))
    nullity = C(P-1,r) * ( 2^r * P/(P-r) - m )        when the system has full row rank
    nrel    = C(P-1,r) * 2^r              Frobenius relations p_{Q+1} = sigma(p_Q)

Consequences, all of which are checked numerically by design_sweep.py:

  (D1) nullity > 0  <=>  2^r * P/(P-r) > m.  For m > 2^r this bounds P:
       P < P* = m*r/(m - 2^r).  For m <= 2^r every P >= r+1 has a kernel.
  (D2) the "clean" design point nullity == C(P-1,r) (kernel exactly the space of
       Frobenius orbits) is  2^r * P/(P-r) = m+1, i.e.  P = r(m+1)/(m+1-2^r).
       It exists only if m+1 > 2^r, and only when that quotient is an integer.
  (D3) the Frobenius pinning stage is a GF(2) linear system with
           rows = nrel * n,  cols = nullity * n,
       so its *memory* is  nrel * nullity * n^2 bits.  This, not the GF(2^n)
       elimination, is what limits r.
"""
from math import comb

SETS = {  # name: (m, n, d, r, claimed security)
    "I": (71, 83, 7, 5, 128), "II": (101, 83, 7, 5, 192), "III": (107, 89, 8, 6, 256),
    "IV": (79, 101, 7, 5, 128), "V": (97, 103, 8, 6, 192), "VI": (107, 103, 8, 6, 256),
    "VII": (83, 149, 8, 5, 128), "VIII": (101, 149, 8, 5, 192), "IX": (109, 157, 9, 6, 256),
}


def feasible(m, r):
    """All P with nullity > 0, plus the clean point if it exists."""
    out = []
    for p in range(r + 1, 60):
        ncols = comb(p, r) * 2 ** r
        nrows = comb(p - 1, r) * m
        if ncols - nrows > 0:
            out.append(p)
        else:
            break
    clean = None
    if m + 1 > 2 ** r:
        num, den = r * (m + 1), (m + 1) - 2 ** r
        if num % den == 0 and num // den > r:
            clean = num // den
    return out, clean


def best(m, n, r):
    """Largest feasible P (most Frobenius relations per kernel bit)."""
    feas, clean = feasible(m, r)
    if not feas:
        return None
    p = clean if clean else max(feas)
    ncols = comb(p, r) * 2 ** r
    nrows = comb(p - 1, r) * m
    nrel = comb(p - 1, r) * 2 ** r
    nullity = ncols - nrows
    return {"P": p, "clean": clean is not None, "ncols": ncols, "nrows": nrows,
            "nullity": nullity, "nrel": nrel,
            "elim_field_mults": nrows * ncols * n,
            "pin_rows": nrel * n, "pin_cols": nullity * n,
            "pin_bytes": nrel * nullity * n * n // 8}


def min_memory(m, n, r, ratio=2.0):
    """Smallest-memory feasible design point.

    The Frobenius pinning matrix has nrel*n rows and nullity*n columns, so its size
    is nrel*nullity*n^2 bits.  Taking the *largest* feasible P maximises the
    redundancy nrel/nullity but blows the memory up quadratically; the useful
    choice is the smallest P whose redundancy nrel/nullity is still >= `ratio`.
    """
    bestrow = None
    for p in range(r + 1, 60):
        ncols = comb(p, r) * 2 ** r
        nrows = comb(p - 1, r) * m
        nullity = ncols - nrows
        if nullity <= 0:
            break
        nrel = comb(p - 1, r) * 2 ** r
        if nrel < ratio * nullity:
            continue
        gb = nrel * nullity * n * n / 8 / 2 ** 30
        if bestrow is None or gb < bestrow["pin_gib"]:
            bestrow = {"P": p, "nullity": nullity, "nrel": nrel,
                       "redundancy": nrel / nullity, "pin_gib": gb,
                       "elim_field_mults": nrows * ncols * n}
    return bestrow


if __name__ == "__main__":
    hdr = ("set", "m", "n", "r", "sec", "P", "clean", "ncols", "nrows", "nullity",
           "nrel", "elim(1e9 mult)", "pin rows x cols", "pin GiB")
    print("%-5s %4s %4s %2s %4s %3s %5s %9s %9s %8s %8s %14s %22s %9s" % hdr)
    for name, (m, n, d, r, sec) in SETS.items():
        b = best(m, n, r)
        if b is None:
            print("%-5s %4d %4d %2d %4d  -- no feasible P (nullity < 0 for every P > r)" %
                  (name, m, n, r, sec))
            continue
        print("%-5s %4d %4d %2d %4d %3d %5s %9d %9d %8d %8d %14.2f %22s %9.1f" % (
            name, m, n, r, sec, b["P"], "yes" if b["clean"] else "no",
            b["ncols"], b["nrows"], b["nullity"], b["nrel"],
            b["elim_field_mults"] / 1e9,
            "%d x %d" % (b["pin_rows"], b["pin_cols"]),
            b["pin_bytes"] / 2 ** 30))
    print()
    print("minimum-memory design point (smallest P with nrel >= 2*nullity):")
    print("  %-5s %3s %10s %10s %10s %14s %10s" %
          ("set", "P", "nullity", "nrel", "redundancy", "elim(1e9 mult)", "pin GiB"))
    for name, (m, n, d, r, sec) in SETS.items():
        b = min_memory(m, n, r)
        if b is None:
            print("  %-5s  -- none" % name)
            continue
        print("  %-5s %3d %10d %10d %10.2f %14.2f %10.1f" %
              (name, b["P"], b["nullity"], b["nrel"], b["redundancy"],
               b["elim_field_mults"] / 1e9, b["pin_gib"]))
    print()
    print("crossover P* = m*r/(m-2^r)  (largest P with a nontrivial GF(2^n) kernel):")
    for name, (m, n, d, r, sec) in SETS.items():
        if m > 2 ** r:
            print("  %-5s m=%3d r=%d  P* = %.2f" % (name, m, r, m * r / (m - 2 ** r)))
        else:
            print("  %-5s m=%3d r=%d  m <= 2^r: every P >= r+1 has a kernel" % (name, m, r))
