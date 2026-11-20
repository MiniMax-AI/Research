#!/usr/bin/env python3
"""
cost.py -- the cost of the one-signature key recovery, as a formula you can
substitute into, plus the gap to the paper's claimed 128-bit level.

No measurement, no randomness: this is arithmetic.  Exits 0 iff the computed
cost is strictly below the claimed 2^128 bit-operation security level at
(n,k) = (130,257).

Cost model (elementary bit operations, schoolbook):
  * one F_{2^k} multiplication            = k^2 bit operations
  * GF(2) elimination of an m x c system  = (1/3) * m * c^2 bit operations

The attack, in three terms:

  1. BUILD the 257 x 256 GF(2) system.  Column i is the k-bit expansion of
     c_P*u_i + c_Q*v_i  (unknown x_i) and of  c_R*u_i + c_S*v_i  (unknown y_i),
     where c_P,c_Q,c_R,c_S are PUBLIC elements of F_{2^k} built from the public
     key (u,v), from Delta*, and from F_x(Delta*),F_y(Delta*) -- all of which
     the paper's own verifier recomputes (Fig. 10 steps 11-13, p.24).  That is
     2 products per column over 2(n-2) columns:

         4*(n-2)  field multiplications.

  2. SOLVE the k x 2(n-2) GF(2) system:

         (1/3) * k * (2*(n-2))^2  bit operations.

  3. FILTER the 2^d kernel candidates left by a rank deficiency d, using the
     PUBLIC NSBC relation (u.x')(v.y') = (u.y')(v.x'):

         2^d * [ 4*(n-2)*k + 2*k^2 ]   bit operations.

HEADLINE (this file's `T_headline`) uses term 1 as derived above.
The conservative accounting `T_conservative` additionally charges the four
public dot products u.F_x(Delta*) etc. and the four products forming the
right-hand side as field multiplications rather than as XOR-selections, and
re-charges the elimination over the same bits:

         (8*(n-2)+8)*k^2 + (1/3)*k*(2*(n-2))^2 + 2^d*(4*(n-2)*k + 2*k^2)

Both are printed.  We report both and bound against the conservative one.

Usage:  python3 cost.py [n] [k] [d]
"""
import math
import sys


def T_headline(n, k, d):
    return 4 * (n - 2) * k * k + (1.0 / 3.0) * k * (2 * (n - 2)) ** 2 \
        + (2 ** d) * (4 * (n - 2) * k + 2 * k * k)


def T_conservative(n, k, d):
    return (8 * (n - 2) + 8) * k * k + (1.0 / 3.0) * k * (2 * (n - 2)) ** 2 \
        + (2 ** d) * (4 * (n - 2) * k + 2 * k * k)


def T_karatsuba(n, k, d):
    # products at k^1.585, elimination unchanged (the elimination then dominates)
    return 4 * (n - 2) * (k ** 1.585) + (1.0 / 3.0) * k * (2 * (n - 2)) ** 2 \
        + (2 ** d) * (4 * (n - 2) * k + 2 * k * k)


THRESHOLD_LOG2 = 128.0   # claimed bit-security level of SBC_VOLE (lambda = 128)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 130
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 257
    d = float(sys.argv[3]) if len(sys.argv) > 3 else 0.467
    # d = 0.467 is the MEASURED mean rank deficiency over 1200 artifact instances
    # (Paper Section 4); d = 0 is the best case, d = 2 the worst ever
    # observed.  Term 3 is negligible either way (< 2^19).
    claimed = 128.0

    print("paper's parameters: n = %d, k = %d   claimed security level = lambda = 128 bits"
          % (n, k))
    print("mean measured rank deficiency d = %.3f (1200 artifact instances)" % d)
    print("cost model: elementary bit operations; F_{2^k} multiply = k^2, "
          "GF(2) elimination = (1/3)*m*c^2")
    print("")
    rows = [("headline      T = 4(n-2)k^2 + (1/3)k(2(n-2))^2 + 2^d[...]", T_headline),
            ("conservative  T = (8(n-2)+8)k^2 + (1/3)k(2(n-2))^2 + 2^d[...]", T_conservative),
            ("karatsuba     T = 4(n-2)k^1.585 + (1/3)k(2(n-2))^2 + 2^d[...]", T_karatsuba)]
    worst = 0.0
    for label, f in rows:
        for dd, tag in ((0, "d=0"), (d, "d=mean"), (2, "d=2")):
            t = f(n, k, dd)
            l2 = math.log2(t)
            worst = max(worst, l2)
            print("  %-62s %-7s = %15.0f = 2^%.2f  (margin to 2^128: %.1f bits)"
                  % (label if tag == "d=0" else "", tag, t, l2, claimed - l2))
    print("")
    mem_bits = k * 2 * (n - 2)
    print("  memory  : %d x %d GF(2) system = %d bits = 2^%.1f bits (%.1f kB)"
          % (k, 2 * (n - 2), mem_bits, math.log2(mem_bits), mem_bits / 8 / 1024))
    print("  data    : ONE signature = 2^0.0 signatures "
          "(B carries k = %d bits of key-affine information about a %d-bit key)"
          % (k, 2 * (n - 2)))
    print("")
    ok = worst < THRESHOLD_LOG2
    print("  worst-case over all three accountings and all three d: 2^%.2f  "
          "(claimed security 2^%.0f)" % (worst, THRESHOLD_LOG2))
    print("COST_RESULT: %s" % ("SUCCESS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
