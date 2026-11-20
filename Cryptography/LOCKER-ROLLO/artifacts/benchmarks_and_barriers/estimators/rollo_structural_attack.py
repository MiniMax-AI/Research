#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rollo_structural_attack.py
==========================
Reproduces the *structural attack on ideal LRPC codes* of the ROLLO 2020-04-21
specification, Section 6.2, page 28 (verbatim from the PDF):

    "the dual of an ideal LRPC code contains much more codeword of weight d than
     a random code. ... For all 1 <= i <= n, Supp(h_i) = F  =>  C^perp contains
     q^n codewords of the same support.  Thus, we have considered an attack in

         O( (n m)^w  *  q^( d*ceil(m/2) - m - n ) )

     for the choice of the parameters of ROLLO-I and ROLLO-II."

The `- n` term is the q^n gain from the ideal-LRPC dual structure; it is invisible
to any random-code estimator (e.g. cryptographic_estimators.RankSDEstimator).

Cross-validated against ePrint 2023/1387 (Song, Zhang, Huang, Wu,
"Blockwise Rank Decoding Problem and LRPC Codes"), which states verbatim:
  * "the original parameters only achieve the 180-bit security"           (Lake  = ROLLO-I-192)
  * "the original parameters only achieve about 103-, 142-, 173-bit security" (Locker = ROLLO-II)

Run:  python3 rollo_structural_attack.py
"""
from math import ceil, log2

# (name, m, n, r, d, claimed_security_bits, published_by_2023_1387_or_None)
SETS = [
    ("ROLLO-I-128  (2020)",  67,  83, 7, 8, 128, None),
    ("ROLLO-I-192  (2020)",  79,  97, 8, 8, 192, 180),
    ("ROLLO-I-256  (2020)",  97, 113, 9, 9, 256, None),
    ("ROLLO-II-128 (2020)",  83, 189, 7, 8, 128, 103),
    ("ROLLO-II-192 (2020)",  97, 193, 8, 8, 192, 142),
    ("ROLLO-II-256 (2020)",  97, 211, 8, 9, 256, 173),
    ("new2ROLLO-I-128",      73,  83, 7, 8, 128, None),
    ("new2ROLLO-I-192",      89,  97, 8, 8, 192, None),
    ("new2ROLLO-I-256",     103, 113, 9, 9, 256, None),
    ("ROLLO-I-256  (2019)", 113,  67, 7, 8, 256, None),
    ("ROLLO-II-256 (2019)", 127, 157, 7, 8, 256, None),
]


def structural_log2(m, n, d, w):
    """log2 of  O((n m)^w * 2^(d*ceil(m/2) - m - n))  for q = 2."""
    return d * ceil(m / 2) - m - n + w * log2(n * m)


if __name__ == "__main__":
    print("=" * 100)
    print("ROLLO 2020-04-21 spec Sec. 6.2 structural attack on ideal LRPC codes, q = 2")
    print("=" * 100)
    print(f"{'set':<21}{'m':>4}{'n':>5}{'d':>3}{'claim':>6}{'exp':>6}"
          f"{'w=2':>8}{'w=2.81':>8}{'2023/1387':>11}  verdict(w=2.81)")
    ok = tot = 0
    for nm, m, n, r, d, sec, pub in SETS:
        e = d * ceil(m / 2) - m - n
        c2, c281 = structural_log2(m, n, d, 2.0), structural_log2(m, n, d, 2.81)
        v = "BELOW CLAIMED" if c281 < sec else "above"
        mark = ""
        if pub is not None:
            tot += 1
            good = abs(c281 - pub) <= 0.5
            ok += good
            mark = f"{pub}{'*' if good else 'X'}"
        print(f"{nm:<21}{m:>4}{n:>5}{d:>3}{sec:>6}{e:>6}{c2:>8.1f}{c281:>8.1f}"
              f"{mark:>11}  {v} ({c281 - sec:+.1f})")
    print(f"\nreproduced ePrint 2023/1387 published values within 0.5 bit : {ok}/{tot}")
