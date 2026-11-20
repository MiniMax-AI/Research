#!/usr/bin/env python3
"""gen_instance.py -- generate a spec-authentic LOCKER/ROLLO Ideal-LRPC instance and
emit the transposed Frobenius data (H_q, U_q)_{q<P} that defines the Degree-1 TSM
("MaxMinors") sparse linear system over L = GF(2^n).

Ideal-LRPC algebraic model:
    K = GF(2^m) = GF(2)[T]/(g),  R = K[X]/(P),  L = GF(2^n) = GF(2)[X]/(P)
    public key  h = x^{-1} y,  supp(x), supp(y) subset F, dim_{F2} F = d
    ciphertext  u = e1 + e2 h, supp(e1) = supp(e2) = E, dim_{F2} E = r
Transposition Phi : R -> Q = L[T]/(g) transposes the n x m bit matrix:
    Phi(v)[t] = sum_{i<n} bit_t(v_i) * X^i  in L.
Partial Frobenius rho squares each L-coefficient; H_q = rho^q(H_0), U_q = rho^q(U_0).

Output: magic 'TSMI', int32 (m,n,r,P,W), then 2P blocks of m*W little-endian uint64.
"""
import argparse, random, struct, sys

def gf_mul(a, b, m, mod):
    r = 0
    while b:
        lsb = b & -b
        r ^= a << (lsb.bit_length() - 1)
        b ^= lsb
    while r.bit_length() > m:
        p = r.bit_length() - 1
        r ^= 1 << p
        r ^= mod << (p - m)
    return r

def gf_inv(a, m, mod):
    # a^(2^m - 2)
    r = 1; base = a; e = (1 << m) - 2
    while e:
        if e & 1: r = gf_mul(r, base, m, mod)
        base = gf_mul(base, base, m, mod)
        e >>= 1
    return r

def gf_sqr(a, m, mod): return gf_mul(a, a, m, mod)

def poly_mul_R(p1, p2, n, Pmod, m, gmod):
    tmp = [0] * (2 * n)
    for i, a in enumerate(p1):
        if a:
            for j, b in enumerate(p2):
                if b: tmp[i + j] ^= gf_mul(a, b, m, gmod)
    for i in range(2 * n - 1, n - 1, -1):
        c = tmp[i]
        if c:
            tmp[i] = 0
            for j in range(n):
                if (Pmod >> j) & 1: tmp[i - n + j] ^= c
    return tmp[:n]

def poly_inv_R(p, n, Pmod, m, gmod):
    def deg(a):
        for i in range(len(a) - 1, -1, -1):
            if a[i]: return i
        return -1
    def trim(a):
        a = a[:]
        while len(a) > 1 and a[-1] == 0: a.pop()
        return a
    def padd(a, b):
        L = max(len(a), len(b)); r = [0] * L
        for i in range(len(a)): r[i] ^= a[i]
        for i in range(len(b)): r[i] ^= b[i]
        return trim(r)
    def pmul(a, b):
        if deg(a) < 0 or deg(b) < 0: return [0]
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y: r[i + j] ^= gf_mul(x, y, m, gmod)
        return trim(r)
    def pdivmod(a, b):
        a = trim(a[:]); q = [0] * max(1, len(a))
        binv = gf_inv(b[-1], m, gmod)
        while deg(a) >= deg(b) and deg(a) >= 0:
            s = deg(a) - deg(b)
            c = gf_mul(a[deg(a)], binv, m, gmod)
            q[s] = c
            for i in range(len(b)):
                if b[i]: a[s + i] ^= gf_mul(c, b[i], m, gmod)
            a = trim(a)
        return trim(q), a
    modpoly = [0] * (n + 1); modpoly[n] = 1
    for j in range(n):
        if (Pmod >> j) & 1: modpoly[j] = 1
    r0, r1 = trim(modpoly), trim(p[:])
    s0, s1 = [0], [1]
    while deg(r1) > 0:
        q, rr = pdivmod(r0, r1)
        r0, r1 = r1, rr
        s0, s1 = s1, padd(s0, pmul(q, s1))
    il = gf_inv(r1[0], m, gmod)
    s1 = [gf_mul(c, il, m, gmod) for c in s1]
    out = [0] * n
    for i in range(min(n, len(s1))): out[i] = s1[i]
    return out

def gf2n_sqr(a, n, Pmod):
    r = 0; i = 0; b = a
    while b:
        if b & 1: r ^= 1 << (2 * i)
        b >>= 1; i += 1
    for q in range(2 * n - 1, n - 1, -1):
        if (r >> q) & 1:
            r ^= 1 << q
            r ^= Pmod << (q - n)
    return r

def main():
    ap = argparse.ArgumentParser()
    for k in ('m', 'n', 'r', 'd', 'P'):
        ap.add_argument('--' + k, type=int, required=True)
    ap.add_argument('--gmod', type=str, required=True)
    ap.add_argument('--Pmod', type=str, required=True)
    ap.add_argument('--out', type=str, required=True)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--random-u', action='store_true',
                    help='NEGATIVE CONTROL: replace the ciphertext by a uniform random u in R')
    a = ap.parse_args()
    m, n, r, d, P = a.m, a.n, a.r, a.d, a.P
    gmod = 0
    for e in a.gmod.split(','):
        if int(e) < m: gmod |= 1 << int(e)
    Pmod = 0
    for e in a.Pmod.split(','):
        if int(e) < n: Pmod |= 1 << int(e)
    W = (n + 63) // 64
    rng = random.Random(a.seed)
    rnd = lambda k: rng.getrandbits(k)

    def span(basis):
        v = 0
        for b in basis:
            if rng.getrandbits(1): v ^= b
        return v

    Fb = [rnd(m) for _ in range(d)]
    while True:
        x = [span(Fb) for _ in range(n)]
        if any(x): break
    y = [span(Fb) for _ in range(n)]
    xinv = poly_inv_R(x, n, Pmod, m, gmod)
    h = poly_mul_R(xinv, y, n, Pmod, m, gmod)
    assert poly_mul_R(x, h, n, Pmod, m, gmod) == y, "keygen inverse check failed"

    Eb = [rnd(m) for _ in range(r)]
    e1 = [span(Eb) for _ in range(n)]
    e2 = [span(Eb) for _ in range(n)]
    if a.random_u:
        u = [rnd(m) for _ in range(n)]          # uniform syndrome: NOT of the form e1 + e2 h
    else:
        u = poly_mul_R(e2, h, n, Pmod, m, gmod)
        u = [u[i] ^ e1[i] for i in range(n)]

    def transpose(v):
        out = [0] * m
        for t in range(m):
            s = 0
            for i in range(n):
                if (v[i] >> t) & 1: s ^= 1 << i
            out[t] = s
        return out
    H0, U0 = transpose(h), transpose(u)
    H = [H0]; U = [U0]
    for q in range(1, P):
        H.append([gf2n_sqr(c, n, Pmod) for c in H[-1]])
        U.append([gf2n_sqr(c, n, Pmod) for c in U[-1]])

    with open(a.out, 'wb') as f:
        f.write(b'TSMI')
        f.write(struct.pack('<5i', m, n, r, P, W))
        for arr in (H, U):
            for el in arr:
                for c in el:
                    f.write(struct.pack('<%dQ' % W, *[(c >> (64 * k)) & ((1 << 64) - 1) for k in range(W)]))
    # also dump the secret support E for end-to-end verification
    with open(a.out + '.secret', 'w') as f:
        f.write('%d %d %d %d\n' % (m, n, r, d))
        for b in Eb: f.write('%x\n' % b)
    print(f"OK m={m} n={n} r={r} d={d} P={P} W={W} -> {a.out}")

if __name__ == '__main__':
    main()
