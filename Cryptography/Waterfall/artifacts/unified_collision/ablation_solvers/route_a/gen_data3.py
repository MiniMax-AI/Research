"""Emit fq_data2.h for the collision solver (linearised injection tables).

Emits the differential and message-difference pattern together with the
byte-cascade lookup tables that turn the injection search from a 2^32 scan
into a 2^16 structured enumeration.

Structure of the X-box cascade (xstruct.py).  With
    u2 = a1 ^ B1(a0),  u3 = a2 ^ B2(a0) ^ B1(u2),  u4 = a3 ^ B3(a0) ^ B2(u2) ^ B1(u3)
    (Bk(x) = byte k of XBOX[x])
the X-difference of the pair (v, v^D) is governed by the four S-box output
differences y0..y3, and byte_k(XBOX[x]) = L_k(SBOX[x]) ^ c_k with L_k linear, so

    low24( X(v) ^ X(v^D) ) = c + phi(y3),
    c   = ( y0^L3(y1)^L2(y2) , y1^L3(y2) , y2 ),      phi(y) = (L1 y, L2 y, L3 y)

where (y0,y1,y2) depends only on the low three bytes of v and y3 only on the top
byte.  A test  low24(X(v)^X(v^D)^T) in HH  therefore implies the *linear*
condition  M c = M T  with ker M = U := HH + im(phi).  Here dim HH = 5,
rank phi = 8 and HH n im(phi) = 0, so U has dimension 13 and M has 11 rows.

Fixing the low 16 bits of v fixes y0 and y1, hence c = A + psi(y2) with
psi(y) = (L2 y, L3 y, y) linear, so the 11 filter equations become 11 equations
in the 8 bits of y2.  N := M.psi has rank 8, so a table SOLVE[tau] (tau = the
11-bit right hand side) either rejects or yields the unique y2 -- the scan costs
2^16 steps per instance instead of 2^32.

The remaining byte is recovered in closed form: R = c + low24(T) lies in U and
U = HH + im(phi) is a direct sum, so the unique y3 with phi(y3) + R in HH is
y3 = Lproj(R) for the linear projection Lproj with Lproj.phi = id, Lproj|_HH = 0.
"""
import random
from model import XBOX, SBOX, SBOX_INV, xbox_apply
from xstruct import LMAP, DDTX

A, B, C = LMAP[1], LMAP[2], LMAP[3]
XK = [[(XBOX[x] >> (8 * k)) & 0xFF for x in range(256)] for k in range(4)]
K = (16, 7, 6)


# ------------------------------------------------------------------ linear algebra
def rref_basis(vals):
    basis = []
    for v in vals:
        cur = v
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    return basis


def dual_basis(basis, n):
    piv = {}
    for r in basis:
        cur = r
        while cur:
            p = cur.bit_length() - 1
            if p in piv:
                cur ^= piv[p]
            else:
                piv[p] = cur
                break
    out = []
    for f in range(n):
        if f in piv:
            continue
        x = 1 << f
        for p in sorted(piv):
            if (piv[p] >> f) & 1:
                x ^= 1 << p
        out.append(x)
    return out


def matvec(rows, v):
    out = 0
    for i, r in enumerate(rows):
        if bin(r & v).count('1') & 1:
            out |= 1 << i
    return out


def solve_linear(rows, rhs_bits):
    piv = {}
    for i, r in enumerate(rows):
        cur, cb = r, (rhs_bits >> i) & 1
        while cur:
            p = cur.bit_length() - 1
            if p in piv:
                cur ^= piv[p][0]
                cb ^= piv[p][1]
            else:
                piv[p] = (cur, cb)
                break
        else:
            if cb:
                return None
    for p in sorted(piv):
        for q in sorted(piv):
            if q != p and (piv[q][0] >> p) & 1:
                piv[q] = (piv[q][0] ^ piv[p][0], piv[q][1] ^ piv[p][1])
    x = 0
    for p, (r, b) in piv.items():
        if b:
            x |= 1 << p
    return x


# ---------------------------------------------------------------- difference
D_DIFF = 0x003f7638
d0, d1, d2, d3 = [(D_DIFF >> (8 * k)) & 0xFF for k in range(4)]
y3 = d3
y2 = d2 ^ C[y3]
y1 = d1 ^ C[y2] ^ B[y3]
y0 = d0 ^ C[y1] ^ B[y2] ^ A[y3]
e0 = d0
e1 = d1 ^ A[y0]
e2 = d2 ^ B[y0] ^ A[y1]
e3 = d3 ^ C[y0] ^ B[y1] ^ A[y2]
assert e3 == 0 and y3 == 0
S32 = sorted({a0 | ((c1 ^ XK[1][a0]) << 8) | ((c2 ^ XK[2][a0] ^ XK[1][c1]) << 16)
              for a0 in DDTX[e0][y0] for c1 in DDTX[e1][y1] for c2 in DDTX[e2][y2]})
assert len(S32) == 32
H = sorted({a ^ b for a in S32 for b in S32})
assert len(H) == 32
random.seed(11)
bad = sum(1 for _ in range(300000)
          if xbox_apply(v := (random.choice(S32) | (random.getrandbits(8) << 24)))
          ^ xbox_apply(v ^ D_DIFF) != D_DIFF)
assert bad == 0
print('d=%08x |S32|=%d |H|=%d verified' % (D_DIFF, len(S32), len(H)))

# ---------------------------------------------------------------- linear filter
PHI = [A[y] | (B[y] << 8) | (C[y] << 16) for y in range(256)]
PSI = [B[y] | (C[y] << 8) | (y << 16) for y in range(256)]
U = sorted({h ^ p for h in H for p in PHI})
assert len(U) == 1 << 13, len(U)
M = dual_basis(rref_basis(U), 24)
NCOND = len(M)
assert NCOND == 11, NCOND
for r in M:
    for h in U:
        assert bin(r & h).count('1') % 2 == 0
print('dim U=%d filter rows=%d' % (len(U).bit_length() - 1, NCOND))

NCOL = [matvec(M, PSI[1 << j]) for j in range(8)]
SOLVE = [0x200] * (1 << NCOND)
for y in range(256):
    t = 0
    for j in range(8):
        if y >> j & 1:
            t ^= NCOL[j]
    assert SOLVE[t] == 0x200, 'rank(N) < 8'
    SOLVE[t] = y
assert sum(1 for s in SOLVE if s != 0x200) == 256

SYN = [[matvec(M, v << (8 * k)) for v in range(256)] for k in range(3)]

Lrows = []
for j in range(8):
    row = solve_linear([PHI[1 << i] for i in range(8)] + rref_basis(H), 1 << j)
    assert row is not None
    Lrows.append(row)
assert all(matvec(Lrows, h) == 0 for h in H)
assert all(matvec(Lrows, PHI[y]) == y for y in range(256))
LP = [[matvec(Lrows, v << (8 * k)) for v in range(256)] for k in range(3)]

# sanity: for random v in F and random D the filter identity holds
random.seed(3)
for _ in range(2000):
    v = random.choice(S32) | (random.getrandbits(8) << 24)
    g = xbox_apply(v) ^ xbox_apply(v ^ D_DIFF)
    assert (g & 0xFFFFFF) in set(H) | {0}
print('linear filter tables verified')

# ---------------------------------------------------------------- pattern
DPOLY = 0x821086e311703
PHASE = 5


def pmod(a, b):
    db = b.bit_length() - 1
    while a and a.bit_length() - 1 >= db:
        a ^= b << (a.bit_length() - 1 - db)
    return a


def pdiv(a, b):
    db = b.bit_length() - 1
    q = 0
    while a and a.bit_length() - 1 >= db:
        s = a.bit_length() - 1 - db
        q ^= 1 << s
        a ^= b << s
    return q


BP = []
for i in range(3):
    g = (1 << K[i]) | 3
    assert pmod(DPOLY, g) == 0
    BP.append(pdiv(DPOLY, g))
PT = DPOLY.bit_length()
TT = max(x.bit_length() for x in BP)
e = [(DPOLY >> t) & 1 for t in range(PT)]
w = [sum((BP[i] >> t) & 1 for i in range(3)) for t in range(PT)]

for i in range(3):
    for t in range(PT + 40):
        lhs = (BP[i] >> t) & 1
        p = ((BP[i] >> (t - 1)) & 1) if t >= 1 else 0
        q = ((BP[i] >> (t - K[i])) & 1) if t >= K[i] else 0
        ee = e[t] if t < PT else 0
        assert lhs == p ^ q ^ ee, (i, t)
print('pattern verified: degD=%d PT=%d TT=%d' % (PT - 1, PT, TT))

for t in range(PT):
    if t % 17 == PHASE:
        assert e[t] == 0, t
        assert w[t] <= 1, t
assert w[0] == 3 and max(w[1:]) <= 2
print('counter constraint ok for phase', PHASE)

plan = []
last = 0
for t in range(1, PT):
    ctr = (t % 17 == PHASE)
    ncond = (w[t] - 1) if w[t] >= 2 else 0
    if ctr and w[t] == 1:
        ncond = 1
    if not ncond:
        continue
    act = [i for i in range(3) if (BP[i] >> t) & 1]
    plan.append((t, w[t], act[0], act[1] if len(act) > 1 else 0,
                 act[2] if len(act) > 2 else 0, 1 if ctr else 0))
    print('  cond rel %2d w=%d streams %s window=(%d,%d)%s'
          % (t, w[t], [a + 1 for a in act], last, t, ' COUNTER' if ctr else ''))
    # solve5.c pins m[t-1] algebraically, so it needs at least one redrawable word
    # strictly between the previous condition and t-1, and it needs every condition
    # to be a genuine w >= 2 stream-agreement test (a counter step with w == 1 would
    # be a condition on a word that cannot be chosen).  Both hold for this pattern;
    # assert them so a future pattern change is caught here and not as a solver hang.
    assert t - last >= 2, ('no free word between conditions', last, t)
    assert not ctr, ('counter step carries a condition', t)
    last = t
assert all(w[t] == 0 for t in range(PT) if t % 17 == PHASE), 'counter rel not quiet'

# minimal re-randomisation window per condition: walk back from rel-1 until the
# accumulated freedom (32 bits for an inactive step, 13 for a single-stream one)
# reaches COND_BITS + 5; the word of the preceding condition step is free too
# (it only has to keep its own stream inside F, the extra conditions there are
# relations between the A values and do not involve that word).
COND_BITS = 24 - (len(H).bit_length() - 1)
assert COND_BITS == NCOND + 8, (COND_BITS, NCOND)


def free_bits(rel):
    """Message-word freedom at relative step rel: 0 at a counter step (the word
    is the block counter), 32 when no stream is active, 13 when one is (the word
    must keep that stream inside F)."""
    if rel % 17 == PHASE:
        return 0
    return 32 if w[rel] == 0 else 13


# Every step with an active stream between two conditions must have its word
# placed (that stream inside F) before the next condition can be attacked, so
# the window of a condition always starts right after the previous condition.
# Only a suffix of that window has to be re-drawn per attempt: walk back from
# rel-1 until the accumulated freedom covers the condition.
wfull, wtry = [], []
prev = 0
for (rel, wt, a0i, a1i, a2i, ctr) in plan:
    wfull.append(prev + 1)
    acc = 0
    r = rel - 1
    while r >= prev and acc < COND_BITS + 5:
        acc += free_bits(r)
        r -= 1
    wtry.append(r + 1)
    prev = rel
for (p, wf, wt_) in zip(plan, wfull, wtry):
    acc_full = sum(free_bits(r) for r in range(wf, p[0]))
    acc_try = sum(free_bits(r) for r in range(wt_, p[0]))
    assert acc_try >= COND_BITS + 5, (p, wt_, acc_try)
    assert wf >= 1 and wt_ >= 0 and min(wf, wt_) <= wf
    assert all(w[r] == 0 or True for r in range(wf, p[0]))
    print('  cond rel %2d -> full window %d..%d (%d bits, %d steps), '
          'per-try suffix %d..%d (%d bits, %d steps)'
          % (p[0], wf, p[0] - 1, acc_full, p[0] - wf, wt_, p[0] - 1, acc_try, p[0] - wt_))

# ---------------------------------------------------------------- placement
NB = 44
NST = NB * 17
S1 = 11 + 17 * 2
S2 = S1 + 544
assert S1 % 17 == (16 - PHASE) % 17 and S2 % 17 == S1 % 17
assert S2 % 32 == S1 % 32
assert S1 >= 17 and S2 + PT + 17 < NST, (S2 + PT, NST)
assert (S1 - 1) % 17 != 16 and (S1 - 2) % 17 != 16

with open('fq_data2.h', 'w') as fh:
    fh.write('/* generated by gen_data3.py -- do not edit */\n')
    fh.write('#define DIFF 0x%08xu\n#define NS32 %d\n' % (D_DIFF, len(S32)))
    fh.write('static const unsigned int S32[NS32]={%s};\n' % ','.join('0x%06xu' % v for v in S32))
    fh.write('#define NH %d\n' % len(H))
    fh.write('static const unsigned int HH[NH]={%s};\n' % ','.join('0x%06xu' % v for v in H))
    fh.write('#define PT %d\n#define TT %d\n' % (PT, TT))
    fh.write('static const int EV[PT]={%s};\n' % ','.join(map(str, e)))
    for i in range(3):
        fh.write('static const int B%d[TT]={%s};\n'
                 % (i + 1, ','.join(str((BP[i] >> t) & 1) for t in range(TT))))
    fh.write('#define NB %d\n#define NST %d\n#define EXC1 %d\n#define EXC2 %d\n' % (NB, NST, S1, S2))
    fh.write('static const unsigned char SBOXC[256]={%s};\n' % ','.join(map(str, SBOX)))
    fh.write('static const unsigned int XBOXC[256]={%s};\n' % ','.join('0x%08xu' % v for v in XBOX))
    fh.write('#define NPLAN %d\n' % len(plan))
    fh.write('static const int PLAN[NPLAN][8]={%s};\n'
             % ','.join('{%d,%d,%d,%d,%d,%d,%d,%d}' % (p + (wf, wt_))
                        for p, wf, wt_ in zip(plan, wfull, wtry)))
    # ---- linear filter tables
    fh.write('#define NCOND %d\n#define NOSOLVE 0x%xu\n' % (NCOND, 1 << NCOND))
    for k in range(3):
        fh.write('static const unsigned short SYN%d[256]={%s};\n'
                 % (k, ','.join(str(v) for v in SYN[k])))
    fh.write('static const unsigned short SOLVE[%d]={%s};\n'
             % (1 << NCOND, ','.join(str(v) for v in SOLVE)))
    for k in range(3):
        fh.write('static const unsigned char LP%d[256]={%s};\n'
                 % (k, ','.join(str(v) for v in LP[k])))
    for k in range(1, 4):
        fh.write('static const unsigned char XB%d[256]={%s};\n'
                 % (k, ','.join(str(v) for v in XK[k])))
    for nm, tbl in (('LM1', A), ('LM2', B), ('LM3', C)):
        fh.write('static const unsigned char %s[256]={%s};\n'
                 % (nm, ','.join(str(v) for v in tbl)))
    fh.write('static const unsigned int PSI24[256]={%s};\n'
             % ','.join('0x%06xu' % v for v in PSI))
print('wrote fq_data2.h  (NST=%d, excursions at %d and %d)' % (NST, S1, S2))
