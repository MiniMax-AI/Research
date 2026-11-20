"""Label extraction from one HOV direction (HOVER Prop. 2), fibre method.

  V      = 9-dim order-2 derivative space of the direction h  (column space of Cat(h))
  c_j    = N_h * (column j of RREF(H_shortened))  in V   -- rank-one tensors
  c_j    = nu(alpha_j) (x) w_j    (after a row-coordinate change)

The 3-dimensional linear fibres  F_j = nu(alpha_j)(x)W  are recovered from the tangent
spaces of the 15-dimensional quadratic ideal.  With V = F_1 + F_2 + F_3 and the
identifications F_i -> F_1 read off a fourth fibre, each rank-one tensor yields its
first-factor coordinates (lambda_1:lambda_2:lambda_3) up to a diagonal rescaling, i.e.
the label alpha_j up to a Moebius transformation.  The residual Moebius is fixed by the
Frobenius action alpha -> alpha^2, which is computable because c_j^[2] is again a
rank-one tensor of the same V.
"""
from __future__ import annotations
from sage.all import Matrix, vector, PolynomialRing


def ideal_quadrics(cols, Fq, dim=9):
    idx = [(i, j) for i in range(dim) for j in range(i, dim)]
    Ev = Matrix(Fq, len(cols), len(idx))
    for r, c in enumerate(cols):
        for k, (i, j) in enumerate(idx):
            Ev[r, k] = c[i] * c[j]
    return Ev.right_kernel().basis_matrix(), idx


def _polar(qv, idx, dim, Fq):
    Bm = Matrix(Fq, dim, dim)
    for k, (i, j) in enumerate(idx):
        if qv[k] and i != j:
            Bm[i, j] = qv[k]
            Bm[j, i] = qv[k]
    return Bm


def tangent(Ib, idx, c, dim, Fq):
    rows = [_polar(qv, idx, dim, Fq) * c for qv in Ib.rows()]
    return Matrix(Fq, rows).right_kernel()


def fibre_of(Tb, Ib, idx, dim, Fq):
    """3-dim fibre inside the 4-dim tangent space Tb.

    Every ideal quadric restricted to Tb is a product L * (linear form), where
    ker L is the fibre.  The polar form of q = L*l has image span{L,l}, so intersecting
    the polar images of two restricted quadrics yields L by pure linear algebra.
    """
    k = Tb.nrows()
    if k != 4:
        return None
    rows = Tb.rows()
    def restrict(qv):
        """polar matrix of q restricted to Tb (k x k)."""
        Bm = _polar(qv, idx, dim, Fq)
        return Matrix(Fq, k, k, lambda a, b: rows[a] * Bm * rows[b])
    imgs = []
    for qv in Ib.rows():
        Pm = restrict(qv)
        if Pm.rank() == 2:
            imgs.append(Pm.column_space())
            if len(imgs) == 2:
                break
    if len(imgs) < 2:
        return None
    inter = imgs[0].intersection(imgs[1])
    if inter.dimension() != 1:
        for qv in Ib.rows():
            Pm = restrict(qv)
            if Pm.rank() == 2:
                it = imgs[0].intersection(Pm.column_space())
                if it.dimension() == 1:
                    inter = it
                    break
    if inter.dimension() != 1:
        return None
    L = inter.basis()[0]                      # linear form on Tb (kills the fibre)
    K = Matrix(Fq, 1, k, list(L)).right_kernel().basis_matrix()
    return K * Tb


class FibreFrame:
    """V = F1 + F2 + F3 with identifications F_i -> F_1 taken from a 4th fibre."""

    def __init__(self, cols, Fq, dim=9, seeds=(0, 1, 2, 3)):
        self.Fq, self.dim = Fq, dim
        self.Ib, self.idx = ideal_quadrics(cols, Fq, dim)
        Fs = []
        for j in seeds:
            T = tangent(self.Ib, self.idx, cols[j], dim, Fq)
            if T.dimension() != 4:
                raise ValueError("tangent space dimension %d" % T.dimension())
            F = fibre_of(T.basis_matrix(), self.Ib, self.idx, dim, Fq)
            if F is None or F.nrows() != 3:
                raise ValueError("fibre extraction failed")
            Fs.append(F)
        self.F1, self.F2, self.F3, self.F4 = Fs
        base = self.F1.stack(self.F2).stack(self.F3)
        if base.rank() != 9:
            raise ValueError("fibres not transverse")
        self.base = base                      # 9 x 9, rows = basis of V
        self.BinvT = base.transpose().inverse()   # x -> coordinates w.r.t. those rows
        # projections:  x -> component in F_i  (as a vector of the 9-dim space)
        self.proj = []
        for i in range(3):
            sel = Matrix(Fq, 3, 9)
            for r in range(3):
                sel[r, 3 * i + r] = Fq.one()
            self.proj.append(sel * self.BinvT)     # 3 x 9 -> coordinates in F_i block
        # identifications psi_{i->1}: F_1 -> F_i  from the graph of F4 in F1 + Fi
        self.psi = [None]
        for i in (1, 2):
            p1 = self.proj[0] * self.F4.transpose()      # 3 x 3 : F4 -> F1 coords
            pi = self.proj[i] * self.F4.transpose()      # 3 x 3 : F4 -> Fi coords
            # for x in F1 (coords u) the element of F4 mapping to it is p1^{-1}u,
            # and its Fi-component is pi * p1^{-1} u
            self.psi.append(pi * p1.inverse())          # 3x3, F1-coords -> Fi-coords

    def lambdas(self, x):
        """First-factor coordinates (1 : r2 : r3) of a rank-one tensor x in V."""
        Fq = self.Fq
        u = [self.proj[i] * x for i in range(3)]         # 3-vectors in block coords
        # identify u_i back to F1 coordinates:  psi[i] maps F1->Fi, so invert
        v = [u[0]]
        for i in (1, 2):
            v.append(self.psi[i].solve_right(u[i]))
        # all v_i are scalar multiples of v[0]; take ratios
        base = None
        for c in v[0]:
            if c:
                base = c
                break
        if base is None:
            return None
        r = [Fq.one()]
        for i in (1, 2):
            rr = None
            for a, b in zip(v[0], v[i]):
                if a:
                    rr = b / a
                    break
            if rr is None:
                return None
            r.append(rr)
        return vector(Fq, r)


# ---------------------------------------------------------------- conic handling
class ConicParam:
    """Standard-form parametrisation of a nondegenerate conic in P^2 over F_q.

    Builds a basis (e0,e1,e2) in which the conic is  x1^2 + x0*x2 = 0,  so that a point
    (x0,x1,x2) on it corresponds to the parameter t = x1/x0 of the Veronese (1,t,t^2).
    """

    def __init__(self, cq, p0, p2, Fq):
        c00, c01, c02, c11, c12, c22 = list(cq)
        self.Fq = Fq
        self.Bm = Matrix(Fq, 3, 3, [[0, c01, c02], [c01, 0, c12], [c02, c12, 0]])
        self.cq = (c00, c01, c02, c11, c12, c22)
        q = self.value
        e0 = vector(Fq, p0)
        e2 = vector(Fq, p2)
        if q(e0) != 0 or q(e2) != 0:
            raise ValueError("base points not on the conic")
        s = e0 * self.Bm * e2
        if s == 0:
            raise ValueError("degenerate choice of conic base points")
        e2 = e2 / s
        # e1: B(e0,.)=0 and B(e2,.)=0, scaled to q(e1)=1
        A = Matrix(Fq, 2, 3, [list(e0 * self.Bm), list(e2 * self.Bm)])
        K = A.right_kernel()
        if K.dimension() != 1:
            raise ValueError("conic degenerate")
        e1 = K.basis()[0]
        v = q(e1)
        if v == 0:
            raise ValueError("conic degenerate (q(e1)=0)")
        e1 = e1 / v.sqrt()
        self.basis = Matrix(Fq, 3, 3, [list(e0), list(e1), list(e2)])
        self.binv = self.basis.transpose().inverse()   # point -> (x0,x1,x2)

    def value(self, x):
        c00, c01, c02, c11, c12, c22 = self.cq
        return (c00 * x[0] ** 2 + c01 * x[0] * x[1] + c02 * x[0] * x[2]
                + c11 * x[1] ** 2 + c12 * x[1] * x[2] + c22 * x[2] ** 2)

    def t(self, P):
        x = self.binv * vector(self.Fq, P)
        if x[0] == 0:
            return None                      # t = infinity
        return x[1] / x[0]


def fit_frobenius(ts, tps, Fq):
    """Fit  f(t) = (A t^2 + B)/(C t^2 + D)  from pairs (t, f(t)); returns (A,B,C,D)."""
    rows = []
    for t, tp in zip(ts, tps):
        if t is None or tp is None:
            continue
        rows.append([t ** 2, Fq.one(), tp * t ** 2, tp])
    if len(rows) < 5:
        return None
    K = Matrix(Fq, rows).right_kernel()
    if K.dimension() < 1:
        return None
    v = K.basis()[0]
    return tuple(v)


def moebius_from_fixed(roots, perm):
    """M with M(0)=r0, M(1)=r1, M(inf)=r2 given as a 2x2 matrix [[a,b],[c,d]]."""
    Fq = roots[0].parent()
    r0, r1, r2 = [roots[i] for i in perm]
    # M(t) = (a t + b)/(c t + d):  M(r0)=0 -> a r0 + b = 0 ; M(r1)=1 ; M(r2)=inf -> c r2+d=0
    a = Fq.one()
    b = r0                                   # char 2
    c = Fq.one()
    d = r2
    # enforce M(r1) = 1:  (a r1 + b) = (c r1 + d)  -> scale (a,b) by (c r1 + d)/(a r1 + b)
    num = c * r1 + d
    den = a * r1 + b
    if den == 0:
        return None
    s = num / den
    return Matrix(Fq, 2, 2, [[a * s, b * s], [c, d]])


def moebius_apply(M, t):
    Fq = M.base_ring()
    if t is None:
        return M[0, 0] / M[1, 0] if M[1, 0] != 0 else None
    den = M[1, 0] * t + M[1, 1]
    if den == 0:
        return None
    return (M[0, 0] * t + M[0, 1]) / den


def moebius_inv_apply(M, t):
    return moebius_apply(M.inverse(), t)


def conic_point_outside(cq, avoid, Fq, dim=3, tries=200):
    """A point of the conic q(x)=0 (coefficients cq over the pairs i<=j) that is not
    projectively equal to any point in `avoid`."""
    c00, c01, c02, c11, c12, c22 = list(cq)
    def qval(x):
        return (c00 * x[0] ** 2 + c01 * x[0] * x[1] + c02 * x[0] * x[2]
                + c11 * x[1] ** 2 + c12 * x[1] * x[2] + c22 * x[2] ** 2)
    Bm = Matrix(Fq, 3, 3, [[0, c01, c02], [c01, 0, c12], [c02, c12, 0]])
    av = [vector(Fq, a) for a in avoid if a is not None]
    import random as _r
    rng = _r.Random(12345)
    for _ in range(tries):
        p0 = av[rng.randrange(len(av))]
        v = vector(Fq, [Fq.random_element() for _ in range(3)])
        qv = qval(v)
        if qv == 0:
            continue
        lam = (p0 * Bm * v) / qv
        p = p0 + lam * v
        if qval(p) != 0:
            continue
        ok = True
        for a in av:
            # proportional?
            sc = None
            for i in range(3):
                if p[i] != 0:
                    sc = a[i] / p[i]
                    break
            if sc is None:
                continue
            if a == sc * p:
                ok = False
                break
        if ok:
            return p
    return None
