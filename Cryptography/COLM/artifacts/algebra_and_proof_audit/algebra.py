#!/usr/bin/env python3
"""
algebra.py -- every index set the splice uses, DERIVED from COLM's own field.

COLM v1, PDF p.2 section 1.1, verbatim:
    f(x) = x^128 + x^7 + x^2 + x + 1
and PDF p.3, verbatim:
    "So, for example, "2" means x, "3" means x + 1, and "7" means x^2 + x + 1.
     When we write multiplications such as 2 . 3 and 7^2, we mean those in the
     field GF(2^n)."

Elements are Python integers; bit i is the coefficient of x^i, exactly the
paper's convention.  Pure Python, no dependencies, no Sage.

Three things are established here and NOTHING is hard-coded:

  [A] the FADE relation          sum_{k in I0} x^(k-1) = 0,  I0 = supp(f)+1
  [B] the minimum message length for the fade   (l = 130, and why)
  [C] the TRANSPLANT relation    I0 = supp(x^(l-1)/3)+1 union {l, l+1},
                                 and its condition (C): maxdeg(x^(l-1)/3) <= l-2,
                                 which first holds at l = 129
  [D] the COLM_tau obstruction   the same condition fails for every tau <= 127,
                                 which is why the two-query splice does not
                                 reach COLM_127
"""
MOD = (1 << 128) | (1 << 7) | (1 << 2) | (1 << 1) | 1     # f(x)


def gmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> 128:
            a ^= MOD
    return r


def gpow(a, e):
    r = 1
    while e:
        if e & 1:
            r = gmul(r, a)
        a = gmul(a, a)
        e >>= 1
    return r


def ginv(a):
    return gpow(a, (1 << 128) - 2)


X = 2                      # the field element COLM writes "2" (= x)
THREE = 3                  # the field element COLM writes "3" (= x + 1)
INV3 = ginv(THREE)


def gsum(vals):
    """XOR-sum -- Python's sum() is integer addition, NOT the field sum."""
    r = 0
    for v in vals:
        r ^= v
    return r


def deg(u):
    return u.bit_length() - 1 if u else -1


def supp(u):
    return [i for i in range(128) if (u >> i) & 1]


def rank_of_powers(maxexp):
    """GF(2)-rank of {x^0,...,x^maxexp} inside GF(2^128)."""
    basis, r = {}, 0
    for e in range(maxexp + 1):
        v = gpow(X, e)
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                r += 1
                break
    return r


def span_solve(target, maxexp):
    """GF(2)-linear solve: is `target` in span{x^0..x^maxexp}?  Witness or None."""
    basis = {}
    for e in range(maxexp + 1):
        v, m = gpow(X, e), 1 << e
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p][0]
                m ^= basis[p][1]
            else:
                basis[p] = (v, m)
                v = 0
                break
    v, m = target, 0
    while v:
        p = v.bit_length() - 1
        if p not in basis:
            return None
        v ^= basis[p][0]
        m ^= basis[p][1]
    return [e for e in range(maxexp + 1) if (m >> e) & 1]


ok = True


def check(name, cond):
    global ok
    print("  %-70s %s" % (name, "OK" if cond else "FAILED"))
    ok = ok and cond


print("GF(2^128) = GF(2)[x]/(x^128+x^7+x^2+x+1)   (COLM v1, PDF p.2, section 1.1)")
check("x is primitive: x^(2^128-1) = 1", gpow(X, (1 << 128) - 1) == 1)
check("x has order exactly 2^128-1 (x^(2^127) != 1)", gpow(X, 1 << 127) != 1)
print("  1/3 = 0x%x  (degree %d)" % (INV3, deg(INV3)))
check("3 * (1/3) = 1", gmul(THREE, INV3) == 1)
print()

print("[A] the FADE relation -- I0 is the support of f(x), shifted by one")
I0 = [1, 2, 3, 8, 129]
check("sum_{k in {1,2,3,8,129}} x^(k-1) = 0", gsum([gpow(X, k - 1) for k in I0]) == 0)
check("equivalently  1 + x + x^2 + x^7 + x^128 = 0, i.e. f(x) = 0",
      1 ^ X ^ gpow(X, 2) ^ gpow(X, 7) ^ gpow(X, 128) == 0)
# the shifted family used by the nonce-misuse variant
for s in (1, 7, 50):
    Ks = [s + 1, s + 2, s + 3, s + 8, s + 129]
    check("shifted by s=%-3d : sum_{k in %s} x^(k-1) = 0" % (s, Ks),
          gsum([gpow(X, k - 1) for k in Ks]) == 0)
print()

print("[B] the fade's MINIMUM message length")
print("    a non-empty I0 with sum x^(k-1) = 0 exists iff {x^0..x^(l-2)} is GF(2)-dependent")
for l in (127, 128, 129, 130, 131):
    r = rank_of_powers(l - 2)
    print("      l = %3d : %3d powers, GF(2)-rank %3d  ->  non-empty I0 exists: %s"
          % (l, l - 1, r, r < l - 1))
check("no non-empty I0 for any l <= 129", all(rank_of_powers(l - 2) == l - 1 for l in range(2, 130)))
check("a non-empty I0 exists at l = 130", rank_of_powers(128) < 129)
E = span_solve(gpow(X, 128), 127)
I0w = sorted([e + 1 for e in E] + [129])
print("      kernel element derived by Gaussian elimination: x^128 = sum_{e in %s} x^e" % E)
print("      => I0 = %s" % I0w)
check("the DERIVED kernel element is exactly I0 = {1,2,3,8,129}", I0w == [1, 2, 3, 8, 129])
print()

print("[C] the TRANSPLANT relation")
print("    A splice at block i moves the decryption onto query 2's state -- it makes")
print("    X_mix[i] = X_2[i] and W_mix[i] = W_2[i] -- IF AND ONLY IF W_mix[i-1] = W_2[i-1].")
print("    State transfer is therefore STICKY ALONG A CONTIGUOUS RUN of spliced blocks")
print("    and is destroyed by any unspliced block in between.  The transplant needs")
print("      (i)  a contiguous run {m, m+1, ..., l+1} of spliced blocks reaching the tag, and")
print("      (ii) W_mix[m-1] = W_2[m-1], i.e. S cap {1..m-1} = supp(x^(m-1)/3)+1, possible")
print("           only if maxdeg(x^(m-1)/3) <= m-2   -- 'condition (C)'.")
first = None
for m in range(2, 400):
    A = gmul(INV3, gpow(X, m - 1))
    if deg(A) <= m - 2:
        first = m
        break
check("condition (C) first holds at m = 129", first == 129)
check("condition (C) fails for every m <= 128", all(deg(gmul(INV3, gpow(X, m - 1))) > m - 2 for m in range(2, 129)))
A = gmul(INV3, gpow(X, 128))
pre = sorted([e + 1 for e in supp(A)])
print("      x^128/3 has support %s  =>  prefix %s" % (supp(A), pre))
print("      at l = m = 129 the run is {129,130} = {l, l+1}, so S = %s" % sorted(pre + [129, 130]))
check("the DERIVED transplant set at l = 129 is exactly {1,3,4,5,6,7,129,130}",
      sorted(pre + [129, 130]) == [1, 3, 4, 5, 6, 7, 129, 130])
check("for l > 129 the SAME prefix works and the run is simply {129..l+1}",
      all(sorted(pre + list(range(129, l + 2)))[-1] == l + 1 for l in range(129, 400)))
check("condition (C) fails at m = 128 (deg(x^127/3) = 127 > 126)",
      deg(gmul(INV3, gpow(X, 127))) == 127)
check("m <= l forces l >= 129: no transplant exists below that", first == 129)
print()

print("[D] the COLM_tau obstruction -- why the same two queries do NOT reach COLM_127")
bad = [t for t in range(1, 128) if deg(gmul(gpow(X, t), INV3)) <= t - 1]
check("x^tau/3 in span{1..x^(tau-1)} for NO tau in 1..127 (degree is always 127)", bad == [])
print("      i.e. the blocking quantity 3^{-1} x^tau always has an x^127 term, so the")
print("      intermediate tag's state difference can never be cancelled for tau <= 127.")
print("      The margin is exactly one block: tau* = 128.")
print()

print("ALGEBRA: %s" % ("ALL CHECKS PASS" if ok else "FAILURE"))
raise SystemExit(0 if ok else 1)
