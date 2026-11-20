#!/usr/bin/env python3
"""
sbcvole.py -- reference implementation of SBC_VOLE (Huth & Joux, ePrint 2024/1537,
Asiacrypt 2025), Figures 3/5/9 (sign) and Figure 10 (verify), together with the
one-signature linear secret-key-recovery attack of Section 4 of the paper.

Field F_{2^k}, k = 257, modulus x^257 + x^12 + 1 (irreducibility checked here).
Secret key: x, y in F_2^{n-2} with n = 130, i.e. 256 bits.

usage:  ./sbcvole.py --demo    (D=5,  tau=3,  ~2 s)
        ./sbcvole.py --full    (D=13, tau=10, the paper's headline parameter set)
        ./sbcvole.py --D 9 --tau 15
"""
import hashlib, os, sys, time, random, argparse, copy

K = 257
MOD = (1 << 257) | (1 << 12) | 1          # x^257 + x^12 + 1
N_DEFAULT = 130


# ---------------------------------------------------------------- F_{2^257} ---
def gf_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> K:
            a ^= MOD
    return r


def gf_inv(a):
    assert a
    r, e, base = 1, (1 << K) - 2, a
    while e:
        if e & 1:
            r = gf_mul(r, base)
        base = gf_mul(base, base)
        e >>= 1
    return r


def modulus_is_irreducible():
    """x^(2^257) == x (mod MOD); 257 is prime so this is iff irreducible."""
    r = 2
    for _ in range(K):
        r = gf_mul(r, r)
    return r == 2


# ------------------------------------------------------------ hashes / PRGs ---
def H2(salt, j, level, idx, node):
    return hashlib.sha256(b"H2|" + salt + b"|%d|%d|%d|" % (j, level, idx) + node).digest()[:16]


def expand_leaf(salt, j, idx, seed):
    o = hashlib.shake_256(b"leaf|" + salt + b"|%d|%d|" % (j, idx) + seed).digest(81)
    return (int.from_bytes(o[0:16], "big"),                       # y-share, 128 b
            int.from_bytes(o[16:49], "big") & ((1 << 256) - 1),   # z-share, in F_{2^257}
            int.from_bytes(o[49:81], "big"))                       # com, 256 b


def H(tag, *parts):
    return hashlib.sha256(tag + b"|" + b"|".join(parts)).digest()


def derive_mu_alpha(h0, tau, D, n):
    need = ((n - 2) + tau * D) * 33
    buf = hashlib.shake_256(b"PRG0|" + h0).digest(need)
    p, out = 0, []
    for _ in range((n - 2) + tau * D):
        out.append(int.from_bytes(buf[p:p + 33], "big") & ((1 << K) - 1))
        p += 33
    return out[:n - 2], out[n - 2:]


def derive_ichall(h1, tau, D):
    buf = hashlib.shake_256(b"PRG1|" + h1).digest(4 * tau)
    return [int.from_bytes(buf[4 * t:4 * t + 4], "big") % (1 << D) for t in range(tau)]


def bs(*vals):
    return b"|".join(b"%d" % v for v in vals)


# ------------------------------------------------------------------- cGGM ----
def _children(salt, j, level, idx, T):
    L = H2(salt, j, level, idx, T)
    return L, bytes(a ^ b for a, b in zip(T, L))


def cggm_leaves(salt, j, s0, s1, depth):
    """Correlated GGM, 2^depth leaves; XOR of all leaves == s0 ^ s1.
    level l holds 2^l nodes (l = 1..depth); seeds sit at level 1."""
    nodes, level = [s0, s1], 1
    while level < depth:
        nxt = []
        for idx, T in enumerate(nodes):
            L, R = _children(salt, j, level, idx, T)
            nxt.append(L); nxt.append(R)
        nodes, level = nxt, level + 1
    return nodes


def cggm_subtree(salt, j, T, level, base, depth):
    """Leaves of the subtree rooted at node T sitting at (level, base)."""
    nodes = [T]
    while level < depth:
        nxt = []
        for off, t in enumerate(nodes):
            L, R = _children(salt, j, level, base + off, t)
            nxt.append(L); nxt.append(R)
        nodes, level, base = nxt, level + 1, base * 2
    return nodes


def cggm_puncture(salt, j, s0, s1, depth, ileaf):
    """Punctured key: the `depth` siblings of the path nodes, level 1..depth."""
    key, level = [], 1
    nodes = [s0, s1]
    path = ileaf >> (depth - 1)          # index of the path node at level 1
    key.append((1, path ^ 1, nodes[path ^ 1]))
    while level < depth:
        nxt = []
        for idx, T in enumerate(nodes):
            L, R = _children(salt, j, level, idx, T)
            nxt.append(L); nxt.append(R)
        nodes, level = nxt, level + 1
        path = ileaf >> (depth - level)
        if level < depth:
            key.append((level, path ^ 1, nodes[path ^ 1]))
        else:
            key.append((level, path ^ 1, nodes[path ^ 1]))
    return key


def cggm_recover(salt, j, key, depth, ileaf):
    """All 2^depth leaves except ileaf; that entry is None."""
    out = [None] * (1 << depth)
    for (level, idx, T) in key:
        base = idx << (depth - level)
        sub = cggm_subtree(salt, j, T, level, idx, depth)
        out[base:base + len(sub)] = sub
    out[ileaf] = None
    return out


def xor_all(vals):
    r = 0
    for v in vals:
        r ^= v
    return r


def folds(vals, depth):
    """[(s0_i, s1_i) for i in 0..depth-1]; s0_i ^ s1_i == XOR of all leaves.
    vals[k] = 0 stands for an unknown leaf (it contributes nothing)."""
    res = [None] * depth
    cur = list(vals)
    for bit in range(depth - 1, -1, -1):
        h = len(cur) // 2
        lo, hi = cur[:h], cur[h:]
        res[bit] = (xor_all(lo), xor_all(hi))
        cur = [a ^ b for a, b in zip(lo, hi)]
    return res


# ------------------------------------------------------------------ scheme ---
class Params:
    def __init__(self, n=N_DEFAULT, D=5, tau=3):
        self.n, self.k, self.D, self.tau = n, K, D, tau


def vec_dot(vec, bits):
    """sum_i vec_i * bits_i over F_{2^k}, bits_i in {0,1}."""
    s = 0
    for i, b in enumerate(bits):
        if b:
            s ^= vec[i]
    return s


def keygen(P, rnd):
    n = P.n
    while True:
        x = [rnd.getrandbits(1) for _ in range(n - 2)]
        y = [rnd.getrandbits(1) for _ in range(n - 2)]
        u = [rnd.getrandbits(K) for _ in range(n)]
        v = [rnd.getrandbits(K) for _ in range(n - 1)]
        xp, yp = x + [1, 0], y + [0, 1]
        uxp, uyp = vec_dot(u, xp), vec_dot(u, yp)
        vxp = vec_dot(v + [0], xp)
        if uxp == 0:
            continue
        v_last = gf_mul(gf_mul(uyp, vxp), gf_inv(uxp)) ^ vec_dot(v + [0], y)
        v.append(v_last)
        if vec_dot(u, xp) == 0 or vec_dot(v, xp) == 0:
            continue
        if gf_mul(vec_dot(u, xp), vec_dot(v, yp)) != gf_mul(vec_dot(u, yp), vec_dot(v, xp)):
            raise AssertionError("planting failed")
        return dict(x=x, y=y, u=u, v=v)


def sign(P, key, msg, salt=None, minmax=False, rnd=None):
    n, D, tau = P.n, P.D, P.tau
    x, y, u, v = key['x'], key['y'], key['u'], key['v']
    xv = sum(x[i] << i for i in range(n - 2))
    yv = sum(y[i] << i for i in range(n - 2))
    salt = salt or os.urandom(32)
    Rx = (rnd or random).randbytes(16) if rnd else os.urandom(16)
    Rxx = bytes(a ^ b for a, b in zip(Rx, xv.to_bytes(16, "big")))

    pre = [int.from_bytes(t, "big") for t in cggm_leaves(salt, 1000, Rx, Rxx, tau)]
    assert xor_all(pre) == xv
    pre_x0 = [0] * tau
    pre_x1 = [0] * tau
    for idx, val in enumerate(pre):
        for j in range(tau):
            if (idx >> j) & 1:
                pre_x1[j] ^= val
            else:
                pre_x0[j] ^= val

    trees = []
    for j in range(tau):
        lv = cggm_leaves(salt, j, pre_x0[j].to_bytes(16, "big"),
                         pre_x1[j].to_bytes(16, "big"), D)
        li = [int.from_bytes(t, "big") for t in lv]
        ys, zs, cs = [], [], []
        for idx, sd in enumerate(lv):
            a, b, c = expand_leaf(salt, j, idx, sd)
            ys.append(a); zs.append(b); cs.append(c)
        trees.append(dict(root=(pre_x0[j], pre_x1[j]), seeds=lv, li=li,
                          ys=ys, zs=zs, cs=cs, Xj=xor_all(cs)))
    for j in range(tau):
        assert xor_all(trees[j]['li']) == xv

    fx = [folds(trees[j]['li'], D) for j in range(tau)]
    fy = [folds(trees[j]['ys'], D) for j in range(tau)]
    fz = [folds(trees[j]['zs'], D) for j in range(tau)]
    fc = [folds(trees[j]['cs'], D) for j in range(tau)]

    Yj = [xor_all(trees[j]['ys']) for j in range(tau)]
    Zj = [xor_all(trees[j]['zs']) for j in range(tau)]
    z = Zj[0]
    dy = [yv ^ Yj[j] for j in range(tau)]
    dz = [z ^ Zj[j] for j in range(tau)]

    if minmax:
        # Sec 6.1: implicit commitment = (global min, global max) of the 2*lambda-bit
        # leaf commitments; abort if a challenge index hits either extremum.
        flat = [(trees[j]['cs'][i], j, i) for j in range(tau) for i in range(1 << D)]
        mn = min(flat); mx = max(flat)
        hcom = H(b"H3", salt, b"%d" % mn[0], b"%d" % mx[0], bs(*dy), bs(*dz[1:]))
        h0 = H(b"H0", msg, hcom, salt, bs(*dy), bs(*dz[1:]))
    else:
        hcom = H(b"H3", salt,
                 bs(*[fc[j][i][b] for j in range(tau) for i in range(D) for b in (0, 1)]),
                 bs(*[trees[j]['Xj'] for j in range(tau)]))
        h0 = H(b"H0", msg, hcom, salt, bs(*dy), bs(*dz[1:]),
               bs(*[trees[j]['Xj'] for j in range(tau)]))
    mu, al = derive_mu_alpha(h0, tau, D, n)

    # A_Fx , A_Fy  (Fig. 5 steps 1-2)
    AFx = [0] * (n - 2)
    AFy = [0] * (n - 2)
    for j in range(tau):
        for i in range(D):
            a = al[j * D + i]
            s0x = fx[j][i][0]
            s0y = fy[j][i][0]
            for t in range(n - 2):
                if (s0x >> t) & 1:
                    AFx[t] ^= a
                if (s0y >> t) & 1:
                    AFy[t] ^= a

    # Phi(Delta) = (u.X(D))(v.Y(D)) + (u.Y(D))(v.X(D)),
    #   X(D) = (A_Fx + D*x, D, 0),  Y(D) = (A_Fy + D*y, 0, D)
    xp, yp = x + [1, 0], y + [0, 1]
    uA0 = sum((gf_mul(u[t], AFx[t]) for t in range(n - 2)), 0) ^ 0
    uA0 = xor_all([gf_mul(u[t], AFx[t]) for t in range(n - 2)])
    vA0 = xor_all([gf_mul(v[t], AFx[t]) for t in range(n - 2)])
    uB0 = xor_all([gf_mul(u[t], AFy[t]) for t in range(n - 2)])
    vB0 = xor_all([gf_mul(v[t], AFy[t]) for t in range(n - 2)])
    uxp, vxp = vec_dot(u, xp), vec_dot(v, xp)
    uyp, vyp = vec_dot(u, yp), vec_dot(v, yp)
    Aval = gf_mul(uA0, vB0) ^ gf_mul(uB0, vA0)
    Bval = (gf_mul(uA0, vyp) ^ gf_mul(uxp, vB0)
            ^ gf_mul(uB0, vxp) ^ gf_mul(uyp, vA0))
    Cval = gf_mul(uxp, vyp) ^ gf_mul(uyp, vxp)          # == 0 by NSBC
    assert Cval == 0, "NSBC not satisfied"

    bval = z ^ vec_dot(mu, y)
    aj = []
    for j in range(tau):
        AFyj = [0] * (n - 2)
        AFzj = 0
        for i in range(D):
            a = al[j * D + i]
            s0y = fy[j][i][0]
            for t in range(n - 2):
                if (s0y >> t) & 1:
                    AFyj[t] ^= a
            AFzj ^= gf_mul(a, fz[j][i][0])
        aj.append(AFzj ^ xor_all([gf_mul(AFyj[t], mu[t]) for t in range(n - 2)]))

    h1 = H(b"H1", msg, hcom, salt, bs(*dy), bs(*dz[1:]), h0,
           b"%d" % Aval, b"%d" % Bval, bs(*aj), b"%d" % bval)
    ich = derive_ichall(h1, tau, D)
    if minmax:
        hit = (mn[1], ich[mn[1]]) == (mn[1], mn[2]) or (mx[1], ich[mx[1]]) == (mx[1], mx[2])
        if hit:
            return None                      # signer aborts and restarts (Sec 6.1)
    keys = [cggm_puncture(salt, j, trees[j]['root'][0].to_bytes(16, "big"),
                          trees[j]['root'][1].to_bytes(16, "big"), D, ich[j])
            for j in range(tau)]
    ip = sum((ich[j] >> (D - 1)) << j for j in range(tau))
    prekey = cggm_puncture(salt, 1000, Rx, Rxx, tau, ip)

    sig = dict(h1=h1, salt=salt, dy=dy, dz=dz[1:], B=Bval, b=bval, keys=keys,
               prekey=prekey, Xj=(None if minmax else [trees[j]['Xj'] for j in range(tau)]),
               minmax=minmax, msg=msg)
    return sig


def verify(P, pk, sig):
    """Figure 10.  Returns (accepted, public view)."""
    n, D, tau = P.n, P.D, P.tau
    u, v = pk['u'], pk['v']
    ich = derive_ichall(sig['h1'], tau, D)

    KL, KY, KZ, KC = [], [], [], []
    for j in range(tau):
        lv = cggm_recover(sig['salt'], j, sig['keys'][j], D, ich[j])
        KL.append([0 if t is None else int.from_bytes(t, "big") for t in lv])
        yy, zz, cc = [], [], []
        for idx, sd in enumerate(lv):
            if sd is None:
                yy.append(0); zz.append(0); cc.append(0)
            else:
                a, b, c = expand_leaf(sig['salt'], j, idx, sd)
                yy.append(a); zz.append(b); cc.append(c)
        KY.append(yy); KZ.append(zz); KC.append(cc)

    # steps 3-5: the com shares, using X_j
    if sig.get('minmax'):
        # verifier recomputes M'min / M'max over the forest with leaf i*_j deleted
        rem = [KC[j][i] for j in range(tau) for i in range(1 << D) if i != ich[j]]
        hcom = H(b"H3", sig['salt'], b"%d" % min(rem), b"%d" % max(rem),
                 bs(*sig['dy']), bs(*sig['dz']))
        dz = [0] + list(sig['dz'])
        h0 = H(b"H0", sig['msg'], hcom, sig['salt'], bs(*sig['dy']), bs(*dz[1:]))
    else:
        fc = []
        for j in range(tau):
            pr = folds(KC[j], D)
            row = []
            for i in range(D):
                b = ((ich[j] >> i) & 1) ^ 1
                known, other = pr[i][b], sig['Xj'][j] ^ pr[i][b]
                row.append((known, other) if b == 0 else (other, known))
            fc.append(row)
        hcom = H(b"H3", sig['salt'],
                 bs(*[fc[j][i][b] for j in range(tau) for i in range(D) for b in (0, 1)]),
                 bs(*sig['Xj']))
        dz = [0] + list(sig['dz'])
        h0 = H(b"H0", sig['msg'], hcom, sig['salt'], bs(*sig['dy']), bs(*dz[1:]), bs(*sig['Xj']))
    mu, al = derive_mu_alpha(h0, tau, D, n)

    # steps 9-12
    Fx = [0] * (n - 2)
    Fy = [0] * (n - 2)
    Delta = 0
    ajp = []
    for j in range(tau):
        Dj = 0
        for i in range(D):
            b = ((ich[j] >> i) & 1) ^ 1
            Dj ^= gf_mul(al[j * D + i], b)
        Delta ^= Dj
        fxj = folds(KL[j], D)
        fyj = folds(KY[j], D)
        fzj = folds(KZ[j], D)
        Fyj = [0] * (n - 2)
        Zj = 0
        for i in range(D):
            b = ((ich[j] >> i) & 1) ^ 1
            a = al[j * D + i]
            sx, sy, sz = fxj[i][b], fyj[i][b], fzj[i][b]
            if b:
                sy ^= sig['dy'][j]
                sz ^= dz[j]
            for t in range(n - 2):
                if (sx >> t) & 1:
                    Fx[t] ^= a
                if (sy >> t) & 1:
                    Fyj[t] ^= a
            Zj ^= gf_mul(a, sz)
        for t in range(n - 2):
            Fy[t] ^= Fyj[t]
        ajp.append(Zj ^ xor_all([gf_mul(Fyj[t], mu[t]) for t in range(n - 2)])
                   ^ gf_mul(sig['b'], Dj))

    # step 13-14
    c, e, f, d = u[n - 2], v[n - 2], u[n - 1], v[n - 1]
    Ux = xor_all([gf_mul(u[t], Fx[t]) for t in range(n - 2)])
    Vx = xor_all([gf_mul(v[t], Fx[t]) for t in range(n - 2)])
    Uy = xor_all([gf_mul(u[t], Fy[t]) for t in range(n - 2)])
    Vy = xor_all([gf_mul(v[t], Fy[t]) for t in range(n - 2)])
    uX = Ux ^ gf_mul(Delta, c)
    vX = Vx ^ gf_mul(Delta, e)
    uY = Uy ^ gf_mul(Delta, f)
    vY = Vy ^ gf_mul(Delta, d)
    Phi = gf_mul(uX, vY) ^ gf_mul(uY, vX)
    Aprime = Phi ^ gf_mul(sig['B'], Delta)
    h1p = H(b"H1", sig['msg'], hcom, sig['salt'], bs(*sig['dy']), bs(*dz[1:]), h0,
            b"%d" % Aprime, b"%d" % sig['B'], bs(*ajp), b"%d" % sig['b'])
    view = dict(Ux=Ux, Vx=Vx, Uy=Uy, Vy=Vy, Delta=Delta, Fx=Fx, Fy=Fy, ich=ich)
    return h1p == sig['h1'], view


# ------------------------------------------------------------------ attack ---
def build_system(P, pk, sig, view, conv=0):
    """
    B = (Vy + (D+1)*d)*P + (Uy + (D+1)*f)*Q + (Vx + (D+1)*e)*R + (Ux + (D+1)*c)*S
        + (c*Vy + d*Ux + e*Uy + f*Vx)
    with Ux = u.(F_x(D*),0,0), Vx = v.(F_x(D*),0,0),
         Uy = u.(F_y(D*),0,0), Vy = v.(F_y(D*),0,0)   -- all verifier-computable.
      P = sum_{i<n-2} u_i x_i,  Q = sum_{i<n-2} v_i x_i
      R = sum_{i<n-2} u_i y_i,  S = sum_{i<n-2} v_i y_i
    => one F_{2^257}-affine equation in the 256 F_2-unknowns (x, y),
       i.e. 257 linear equations over F_2 in 256 unknowns.
    """
    n = P.n
    u, v = pk['u'], pk['v']
    Ux, Vx, Uy, Vy, Dl = view['Ux'], view['Vx'], view['Uy'], view['Vy'], view['Delta']
    c, e, f, d = u[n - 2], v[n - 2], u[n - 1], v[n - 1]
    if conv == 0:        # X(D)=(A_Fx+Dx, D, 0), Y(D)=(A_Fy+Dy, 0, D)  -- the only
        t = Dl           #   convention under which an honest signature verifies
    elif conv == 1:      # X(D)=(A_Fx+Dx, 1, 0), Y(D)=(A_Fy+Dy, 0, 1), true D-coeff
        t = 0
    else:                # Figure 5 exactly as printed (A_F padded with (1,0)/(0,1))
        t = Dl ^ 1
    cP = Vy ^ gf_mul(t, d)
    cQ = Uy ^ gf_mul(t, f)
    cR = Vx ^ gf_mul(t, e)
    cS = Ux ^ gf_mul(t, c)
    rhs = sig['B'] ^ gf_mul(Ux, d) ^ gf_mul(c, Vy) ^ gf_mul(Uy, e) ^ gf_mul(f, Vx)
    if conv == 1:
        rhs = sig['B']
    cols = [gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i]) for i in range(n - 2)]
    cols += [gf_mul(cR, u[i]) ^ gf_mul(cS, v[i]) for i in range(n - 2)]
    ncol = len(cols)
    rows = []
    for bit in range(K):
        r = 0
        for j, col in enumerate(cols):
            if (col >> bit) & 1:
                r |= 1 << j
        rows.append(r | (((rhs >> bit) & 1) << ncol))
    return rows, ncol


def solve_gf2(rows, ncol):
    piv = {}
    for r in rows:
        cur = r
        while True:
            low = cur & ((1 << ncol) - 1)
            if low == 0:
                break
            p = low.bit_length() - 1
            if p in piv:
                cur ^= piv[p]
            else:
                piv[p] = cur
                break
    rank = len(piv)
    if rank < ncol:
        return None, rank
    sol = 0
    for p in sorted(piv):
        r = piv[p]
        val = (r >> ncol) & 1
        val ^= bin((r & ((1 << ncol) - 1) & ~(1 << p)) & sol).count("1") & 1
        if val:
            sol |= 1 << p
    return sol, rank


def solve_gf2_candidates(rows, ncol, max_cand=4):
    """All solutions of the affine system when the kernel has dimension <= log2(max_cand).
    Returns (list_of_solutions, rank)."""
    piv = {}
    for r in rows:
        cur = r
        while True:
            low = cur & ((1 << ncol) - 1)
            if low == 0:
                if (cur >> ncol) & 1:          # 0 = 1, inconsistent
                    return [], len(piv)
                break
            p = low.bit_length() - 1
            if p in piv:
                cur ^= piv[p]
            else:
                piv[p] = cur
                break
    rank = len(piv)
    nfree = ncol - rank
    if (1 << nfree) > max_cand:
        return [], rank
    free = [j for j in range(ncol) if j not in piv]
    # back-substitute for each assignment of the free variables
    sols = []
    for mask in range(1 << nfree):
        sol = 0
        for t, j in enumerate(free):
            if (mask >> t) & 1:
                sol |= 1 << j
        for p in sorted(piv):
            r = piv[p]
            val = (r >> ncol) & 1
            val ^= bin((r & ((1 << ncol) - 1) & ~(1 << p)) & sol).count("1") & 1
            if val:
                sol |= 1 << p
        sols.append(sol)
    return sols, rank


def attack_one(P, pk, sig, view, convs=(0, 1, 2)):
    """ONE signature.  The 257x256 GF(2) system is over-determined by one bit, so
    its kernel has dimension 0 (probability ~0.58) or 1; in the latter case there are
    exactly 2 candidate keys and the NSBC relation for the public key decides between
    them at zero cost.  Hence a single signature always suffices."""
    for conv in convs:
        rows, ncol = build_system(P, pk, sig, view, conv)
        sols, rank = solve_gf2_candidates(rows, ncol)
        for sol in sols:
            if nsbc_ok(P, pk, sol):
                return sol, rank, conv, len(sols)
    return None, -1, -1, 0


def attack(P, pk, sigs, views, convs=(0, 1, 2)):
    """Try every padding convention; accept a solution only if it solves NSBC
    for the public key (a check that uses no secret)."""
    best = None
    for conv in convs:
        rows, ncol = None, None
        for sg, vw in zip(sigs, views):
            r, nc = build_system(P, pk, sg, vw, conv)
            rows = r if rows is None else rows + r
            ncol = nc
        sol, rank = solve_gf2(rows, ncol)
        if sol is not None and nsbc_ok(P, pk, sol):
            return sol, rank, conv
        best = best or (sol, rank, conv)
    return best


# -------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--D", type=int)
    ap.add_argument("--tau", type=int)
    ap.add_argument("--minmax", action="store_true")
    a = ap.parse_args()
    D, tau = (13, 10) if a.full else ((a.D, a.tau) if a.D else (5, 3))

    print("F_{2^257} modulus x^257+x^12+1 irreducible:", modulus_is_irreducible())
    P = Params(D=D, tau=tau)
    mm = a.minmax
    print("parameters n=%d k=%d D=%d tau=%d  variant=%s  (paper: n=130, k=257)"
          % (P.n, P.k, D, tau, "SBC_VOLE^M (MinMax)" if mm else "SBC_VOLE (basic)"))
    rnd = random.Random(int(os.environ.get("SEED", "20241001")))

    t = time.time(); key = keygen(P, rnd); print("keygen: NSBC solution planted+checked [%.2fs]" % (time.time() - t))
    pk = dict(u=key['u'], v=key['v'])
    xv = sum(key['x'][i] << i for i in range(P.n - 2))
    yv = sum(key['y'][i] << i for i in range(P.n - 2))

    tries = 0
    t = time.time()
    while True:
        tries += 1
        sig = sign(P, key, b"hello world", minmax=mm)
        if sig is not None:
            break
        if tries > 500:
            print("RESULT: FAIL"); return
    ts = time.time() - t
    print("sign [%.2fs]%s" % (ts, ("  (MinMax aborts before success: %d)" % (tries - 1)) if mm else ""))
    t = time.time(); ok, view = verify(P, pk, sig); tv = time.time() - t
    print("verify: %s  [%.2fs]" % (ok, tv))
    if not ok:
        print("reference implementation does not verify -- aborting")
        print("RESULT: FAIL"); return

    # ---- the attack -------------------------------------------------------
    def hit(s):
        return (s is not None and (s & ((1 << (P.n - 2)) - 1)) == xv
                and (s >> (P.n - 2)) == yv)

    t = time.time()
    sol1, rank1, conv1, ncand = attack_one(P, pk, sig, view)
    t1 = time.time() - t
    rows, ncol = build_system(P, pk, sig, view, 0)
    print("attack/1 signature: ONE F_{2^257} equation in the 4 field elements "
          "P=u.x, Q=v.x, R=u.y, S=v.y")
    print("attack/1 signature: = %d GF(2) equations (its 257 bit-slices) in the %d "
          "unknown bits of (x,y)" % (P.k, ncol))
    print("attack/1 signature: rank = %d, %d candidate key(s) after kernel enumeration; "
          "NSBC[u,v] picks one -> key recovered: %s   [%.3fs]"
          % (rank1, ncand, hit(sol1), t1))

    sig2 = sign(P, key, b"second message", minmax=mm, rnd=rnd)
    ok2, view2 = verify(P, pk, sig2)
    print("second, independent signature verifies:", ok2)

    rows2, _ = build_system(P, pk, sig2, view2, 0)
    t = time.time(); sol2, rank2 = solve_gf2(rows + rows2, ncol); t2 = time.time() - t
    print("attack/2 signatures (no NSBC test needed): rank = %d (need %d) -> key recovered: %s"
          "   [%.3fs]" % (rank2, ncol, hit(sol2), t2))
    print("attack: recovered key solves NSBC[u,v] for the PUBLIC key: %s"
          % (sol1 is not None and nsbc_ok(P, pk, sol1)))

    # ---- controls ---------------------------------------------------------
    bad = copy.deepcopy(sig); bad['B'] ^= (1 << 77)
    sb, rkb, _, _ = attack_one(P, pk, bad, view)
    ctrl1 = hit(sb)
    print("control 1 (one bit of B flipped in signature 1) still recovers the key: %s" % ctrl1)

    sig3 = sign(P, key, b"third message", minmax=mm)
    ok3, view3 = verify(P, pk, sig3)
    ctrl2 = (predict_B(P, pk, view3, sol2) == sig3['B'])
    print("control 2: third signature verifies: %s ; its B is predicted by the recovered key: %s"
          % (ok3, ctrl2))

    wrong = sol2 ^ (1 << 3)
    ctrl3 = (predict_B(P, pk, view3, wrong) == sig3['B'])
    print("control 3 (recovered key with one bit flipped) reproduces B:", ctrl3)

    import math
    ops = ncol ** 3 // 3
    build = 2 * ncol * (P.k ** 2)
    print("cost: system build 2*%d*%d^2 = %d bit-ops ~ 2^%.1f ; elimination (1/3)*%d^3 = %d "
          "~ 2^%.1f" % (ncol, P.k, build, math.log2(build), ncol, ops, math.log2(ops)))
    print("cost: TOTAL ~ 2^%.1f bit-ops ; DATA = 1 signature" % math.log2(build + ops))
    print("measured wall clock for the whole attack: %.3f s" % (t1 + t2))

    if hit(sol2) and not ctrl1 and ctrl2 and not ctrl3:
        print("RESULT: SUCCESS")
    else:
        print("RESULT: FAIL")


def nsbc_ok(P, pk, sol):
    """Public check: does the candidate (x,y) solve NSBC[u,v]?"""
    n = P.n
    xb = [(sol >> i) & 1 for i in range(n - 2)] + [1, 0]
    yb = [(sol >> (n - 2 + i)) & 1 for i in range(n - 2)] + [0, 1]
    return (gf_mul(vec_dot(pk['u'], xb), vec_dot(pk['v'], yb))
            == gf_mul(vec_dot(pk['u'], yb), vec_dot(pk['v'], xb)))


def predict_B(P, pk, view, sol):
    n = P.n
    u, v = pk['u'], pk['v']
    xb = [(sol >> i) & 1 for i in range(n - 2)]
    yb = [(sol >> (n - 2 + i)) & 1 for i in range(n - 2)]
    c, e, f, d = u[n - 2], v[n - 2], u[n - 1], v[n - 1]
    Dl = view['Delta']
    cP = view['Vy'] ^ gf_mul(Dl, d)
    cQ = view['Uy'] ^ gf_mul(Dl, f)
    cR = view['Vx'] ^ gf_mul(Dl, e)
    cS = view['Ux'] ^ gf_mul(Dl, c)
    r = (gf_mul(view['Ux'], d) ^ gf_mul(c, view['Vy'])
         ^ gf_mul(view['Uy'], e) ^ gf_mul(f, view['Vx']))
    for i in range(n - 2):
        if xb[i]:
            r ^= gf_mul(cP, u[i]) ^ gf_mul(cQ, v[i])
        if yb[i]:
            r ^= gf_mul(cR, u[i]) ^ gf_mul(cS, v[i])
    return r


if __name__ == "__main__":
    main()
