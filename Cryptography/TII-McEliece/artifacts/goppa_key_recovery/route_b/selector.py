"""HOVER selector (ePrint 2026/1778, Alg. 1): from a binary HOV kernel to the
Frobenius-ordered hidden directions.

  Cat(f)[r,S] = coefficient of x^S in df/dx_r,  |S| = p-1  (eq. 6)
  beta(z,lambda) = (sum_i lambda_i M_i)^T z   (eq. 7)
  T            = row-compressed linearisation of beta
  E_beta       = {B in End_F2(L) : (I_kappa ⊗ B)(ker T) ⊆ ker T}   (eq. 8/11)
"""
from __future__ import annotations
from itertools import combinations
from sage.all import GF, Matrix, vector, PolynomialRing, identity_matrix


def cat_matrix(coeffs, kappa, F):
    """coeffs: dict {(i,j,l): c}.  Returns kappa x C(kappa,2) matrix Cat(f)."""
    pairs = list(combinations(range(kappa), 2))
    pidx = {p: k for k, p in enumerate(pairs)}
    M = Matrix(F, kappa, len(pairs))
    for (i, j, l), c in coeffs.items():
        if c == 0:
            continue
        M[i, pidx[(j, l)]] += c
        M[j, pidx[(i, l)]] += c
        M[l, pidx[(i, j)]] += c
    return M


def rows_to_forms(K, triples, F):
    """Binary kernel basis rows -> list of coefficient dicts over F."""
    out = []
    for row in K:
        d = {}
        for c, tr in enumerate(triples):
            if row[c]:
                d[tr] = F.one()
        out.append(d)
    return out


class Selector:
    def __init__(self, Kbasis, triples, kappa, F2=GF(2)):
        """Kbasis: rows (binary) of the HOV kernel; triples: column monomial list."""
        import numpy as np
        self.F2 = F2
        self.kappa = kappa
        self.triples = triples
        self.d = len(Kbasis)
        self.forms = rows_to_forms(Kbasis, triples, F2)
        self.tindex = {tr: k for k, tr in enumerate(triples)}
        self.N = kappa * (kappa - 1) // 2
        pairs = list(combinations(range(kappa), 2))
        self.pidx = {p: k for k, p in enumerate(pairs)}
        # catalecticants of the kernel basis, built with numpy scatter (eq. 6)
        self.M = []
        for f in self.forms:
            Marr = np.zeros((kappa, self.N), dtype=np.uint8)
            rows, cols = [], []
            for (i, j, l) in f:
                rows += (i, j, l)
                cols += (self.pidx[(j, l)], self.pidx[(i, l)], self.pidx[(i, j)])
            if rows:
                Marr[np.array(rows), np.array(cols)] ^= 1
            self.M.append(Marr)
        self.trip = np.array(triples, dtype=np.int64)
        # trip_pairs[c, k] = index of the pair {triple c} \ {triple c[k]}
        self.trip_pairs = np.array(
            [[self.pidx[(tr[1], tr[2])], self.pidx[(tr[0], tr[2])],
              self.pidx[(tr[0], tr[1])]] for tr in triples], dtype=np.int64)

    def build_T(self):
        """T_e in F2^{N x kappa*d}, row-compressed to full row rank."""
        import numpy as np
        from bitsolve import row_basis_indices
        F2 = self.F2
        kap, d, N = self.kappa, self.d, self.N
        Te = np.zeros((N, kap * d), dtype=np.uint8)
        for i, f in enumerate(self.forms):
            if not f:
                continue
            cols = np.fromiter((self.tindex[tr] for tr in f), dtype=np.int64,
                               count=len(f))
            rs = self.trip[cols]                       # (len(f), 3) vertices
            prs = self.trip_pairs[cols]                # (len(f), 3) complement pairs
            Te[prs.ravel(), rs.ravel() * d + i] = 1
        piv = row_basis_indices(Te)                     # independent ROWS of Te
        self.rho = len(piv)
        self.T = Matrix(F2, self.rho, kap * d, Te[np.array(piv, dtype=np.int64)]
                        .flatten().tolist())
        return self.T

    def algebra(self):
        """E_beta as a list of basis matrices (d x d) over F2, via eq. (11).

        Equation (s,z):  sum_{r,a,b} T[s,(r,a)] B[a,b] z_r[b] = 0, i.e. the coefficient
        of the unknown B[a,b] is  sum_r T[s,(r,a)] z[r*d+b] = (U Z_z)[(s,a), b].
        """
        import numpy as np
        from bitsolve import row_basis_indices
        F2 = self.F2
        d, kap = self.d, self.kappa
        T = self.T
        rho = T.nrows()
        Tnp = np.array(T.rows(), dtype=np.uint8).reshape(rho, kap, d)
        U = np.ascontiguousarray(Tnp.transpose(0, 2, 1)).reshape(rho * d, kap)
        Z = T.right_kernel().basis_matrix()          # nz x (kap*d)
        nz = Z.nrows()
        Zb = np.array(Z.rows(), dtype=np.uint8).reshape(nz, kap, d)
        # coef[z][(s,a)][b]
        coef = (np.einsum('xr,zrb->zxb', U, Zb, dtype=np.int64) & 1).astype(np.uint8)
        Earr = np.ascontiguousarray(coef.reshape(nz * rho, d * d))
        # only the row space of Earr matters for its right kernel, so compress first:
        # turning all 5.8e6 entries into a Sage matrix costs ~1.8 s, this costs ~0.25 s.
        idx = row_basis_indices(Earr)
        E = Matrix(F2, len(idx), d * d,
                   Earr[np.array(idx, dtype=np.int64)].flatten().tolist())
        K = E.right_kernel().basis_matrix()
        self.Edim = K.nrows()
        self.Ebasis = [Matrix(F2, d, d, list(row)) for row in K.rows()]
        return self.Ebasis

    def field_generator(self, rng):
        """Random element of E_beta with irreducible degree-d char. polynomial."""
        F2 = self.F2
        d = self.d
        Rx = PolynomialRing(F2, 'x')
        tries = 0
        while True:
            tries += 1
            co = [rng(0) for _ in range(self.Edim)]
            if not any(co):
                continue
            B = sum((c * M for c, M in zip(co, self.Ebasis) if c), Matrix(F2, d, d))
            cp = B.charpoly(Rx.gen())
            if cp.degree() == d and cp.is_irreducible():
                return B, cp, tries


def eigenlines(B, cp, Fq):
    """Split B over Fq: returns list of (eigenvalue, eigenvector) with v^[2] next."""
    d = B.nrows()
    Rx = PolynomialRing(Fq, 'x')
    roots = cp.change_ring(Fq).roots()
    out = []
    for lam, mult in roots:
        K = (B.change_ring(Fq) - lam * identity_matrix(Fq, d)).right_kernel()
        assert K.dimension() == 1, "eigenspace not 1-dimensional"
        out.append((lam, K.basis()[0]))
    # order by Frobenius: squaring an eigenvector squares its eigenvalue
    byval = {lam: v for lam, v in out}
    lam0 = sorted(byval)[0]
    order = []
    lam = lam0
    for _ in range(d):
        order.append((lam, byval[lam]))
        lam = lam ** 2
    assert len(set(l for l, _ in order)) == d, "eigenvalues do not form one Frobenius cycle"
    return order
