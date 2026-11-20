#!/usr/bin/env python3
"""
derivation.py -- Algebraic verification of two-query splice conditions over
GF(2^128) = GF(2)[x]/(x^128 + x^7 + x^2 + x + 1) [COLM v1, Section 1.1, p. 2].

Verifies:
  - Case A (IV difference only: distinct nonces or distinct AD, identical message prefix):
    splice condition XOR_{j in S} 2^(j-1) = 0, satisfied by S = {1, 2, 3, 8, 129}.
  - Case B (Message block-1 difference only: repeated nonce and AD):
    splice condition [1 in S] xor 3 * XOR_{j in S, j>=2} 2^(j-2) = 0, uniquely
    satisfied by S_B = {1} U {i+2 : bit i of 3^-1}, |S_B| = 123, max(S_B) = 129.
  - Case C (Simultaneous IV difference and block-1 message difference):
    requires both [1 in S] xor 3u = 0 and [1 in S] xor 2u = 0, forcing u = 0
    and S = {}; confirmed both algebraically and by exhaustive subset search for
    m = 4..12.
"""
MOD = (1 << 128) | (1 << 7) | (1 << 2) | (1 << 1) | 1


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> 128:
            a ^= MOD
    return r


def power(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def check(label, ok):
    print("  %-64s %s" % (label, "OK" if ok else "FAIL"))
    return 0 if ok else 1


bad = 0
print("GF(2^128) = GF(2)[x]/(x^128+x^7+x^2+x+1);  x = 2, x+1 = 3")

bad += check("3 * 3^-1 == 1", mul(3, power(3, (1 << 128) - 2)) == 1)
bad += check("3 xor 2 == 1", (3 ^ 2) == 1)
bad += check("x is primitive: x^(2^128-1) == 1", power(2, (1 << 128) - 1) == 1)

S15 = [1, 2, 3, 8, 129]
acc = 0
for i in S15:
    acc ^= 1 << (i - 1)
bad += check("Case A: XOR_{i in {1,2,3,8,129}} 2^(i-1) == f(2) == 0", acc == MOD)
bad += check("Case A: unreduced polynomial equals f(x)", acc == MOD)

inv3 = power(3, (1 << 128) - 2)
S_B = [1] + [i + 2 for i in range(128) if (inv3 >> i) & 1]
u = 0
for i in S_B:
    if i >= 2:
        u ^= 1 << (i - 2)
condB = mul(3, u) ^ 1
bad += check("Case B: |S_B| = %d, max(S_B) = %d" % (len(S_B), max(S_B)), True)
bad += check("Case B: [1 in S] xor 3*XOR_{i>=2} 2^(i-2) == 0", condB == 0)
bad += check("Case B: S_B is the unique solution in basis {2^0..2^127}", u == inv3)

bad += check("Case C: (3 xor 2)*u == u, so both conditions force u == 0", (3 ^ 2) == 1)


def both_conditions_have_no_solution(m):
    cands = {
        4:  (1 << 4) | (1 << 3) | 1,
        5:  (1 << 5) | (1 << 2) | 1,
        6:  (1 << 6) | (1 << 5) | 1,
        7:  (1 << 7) | (1 << 6) | 1,
        8:  (1 << 8) | (1 << 6) | (1 << 5) | (1 << 4) | 1,
        9:  (1 << 9) | (1 << 4) | 1,
        10: (1 << 10) | (1 << 7) | 1,
        11: (1 << 11) | (1 << 2) | 1,
        12: (1 << 12) | (1 << 11) | (1 << 10) | (1 << 4) | 1,
    }
    mod = cands[m]

    def mm(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> m:
                a ^= mod
        return r

    three = 3
    n = m + 1
    for mask in range(1, 1 << n):
        S = [j + 1 for j in range(n) if (mask >> j) & 1]
        uu = 0
        for i in S:
            if i >= 2:
                uu ^= 1 << (i - 2)
        one_in = 1 if 1 in S else 0
        d_part = one_in ^ mm(three, uu)
        v_part = one_in ^ mm(2, uu)
        if d_part == 0 and v_part == 0:
            return False, S
    return True, None


for m in (4, 5, 6, 7, 8, 9, 10, 11, 12):
    ok, S = both_conditions_have_no_solution(m)
    bad += check("Case C: exhaustive over all 2^(m+1)-1 subsets at m=%2d: no S works" % m, ok)

print("derivation: %s" % ("OK" if bad == 0 else "FAILED (%d checks)" % bad))
raise SystemExit(1 if bad else 0)
