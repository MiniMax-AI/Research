"""Stages 1-3 of Route B (Frobenius-reduced HOVER): shortened-dual HOV kernel,
incremental E_beta stabilization, and single-direction Frobenius rank selection.

Applies two exact reductions from the Frobenius structure of the extension field:

  (F1)  The m = 10 Frobenius-ordered directions satisfy  h_{r+1} = mu_r * h_r^{[2]}
        (entrywise), hence  supp(h_{r+1}) = supp(h_r)  and
        Cat(h_{r+1}) = mu_r * Cat(h_r)^{[2]} ,  so  rank Cat(h_{r+1}) = rank Cat(h_r).
        The selector's rank condition  sum_r rank Cat(h_r) == rho  therefore collapses to
        the single-rank test  m * rank Cat(h_0) == rho , and the m-1 = 9 conjugate
        catalecticants -- none of which is used downstream (hover_fibers consumes
        cats[0] only) -- need not be built.  Cost of the stage drops by a factor m.

  (F2)  E_beta has dimension m, so the stacked system of eq. (11) has rank exactly
        d^2 - m = rho; each kernel vector z already contributes rho rows, so the stack
        can be grown one z at a time and stopped as soon as its kernel has dimension m
        (which then *equals* E_beta).  This removes the factor nz = C(kappa,2) - rho
        from the dominant step of the algebra computation.
"""
import sys, time
import random as rnd
from sage.all import *
from hovlib import hov_params, shorten, hov_kernel
from selector import Selector, rows_to_forms, cat_matrix, eigenlines
from selincr import algebra_incremental


def directions_from_dual(H, m, t, n, Fq, short_cols=None, p=3, rng=None, verbose=True):
    """Returns dict with the HOV kernel, the single Frobenius direction h_0 used
    downstream, M = RREF(Hs) and the kept columns."""
    e, sstar, kappa, a, b = hov_params(m, t, p)
    if short_cols is None:
        short_cols = list(range(sstar))
    keep = [c for c in range(n) if c not in set(short_cols)]
    Hs = shorten(H, short_cols)
    t0 = time.time()
    Kb, tri, P = hov_kernel(Hs, p)
    t_k = time.time() - t0
    if verbose:
        print(f"  HOV kernel dim {Kb.nrows()} ({t_k:.1f}s)  [m={m}]", flush=True)
    if Kb.nrows() != m:
        return None
    t1 = time.time()
    sel = Selector(Kb.rows(), tri, kappa)
    sel.build_T()
    t_T = time.time() - t1
    if verbose:
        print(f"  rho={sel.rho} (expect {m*p*p})  [build_T {t_T:.1f}s]", flush=True)
    t1 = time.time()
    Eb, used = algebra_incremental(sel, m)
    t_E = time.time() - t1
    if Eb is None:
        print(f"  E_beta failed after {used} z-blocks", flush=True)
        return None
    if verbose:
        print(f"  dim E_beta={sel.Edim} (expect {m}) via {used} z-block(s), "
              f"{sel.Echeck} extra z re-checked  [{t_E:.1f}s]", flush=True)
    B, cp, tries = sel.field_generator(lambda i: (rng or rnd).randrange(2))
    order = eigenlines(B, cp, Fq)
    # (F1): only the first Frobenius direction is materialised
    t1 = time.time()
    binforms = rows_to_forms(Kb.rows(), tri, Fq)

    def combine(vec):
        d = {}
        for c, f in zip(vec, binforms):
            if c:
                for k, val in f.items():
                    d[k] = d.get(k, Fq.zero()) + c * val
        return {k: v for k, v in d.items() if v}

    h0 = combine(order[0][1])
    cat0 = cat_matrix(h0, kappa, Fq)
    r0 = cat0.rank()
    t_h = time.time() - t1
    if verbose:
        print(f"  rank Cat(h_0)={r0}, m*rank={m*r0} (rho={sel.rho})  "
              f"[direction {t_h:.1f}s]", flush=True)
    if m * r0 != sel.rho:          # == sum_r rank Cat(h_r) by (F1)
        return None
    M = Hs.echelon_form()
    return dict(sel=sel, hs=[h0], cats=[cat0], ranks=[r0], kappa=kappa, e=e,
                sstar=sstar, a=a, b=b, M=M, keep=keep, short_cols=list(short_cols),
                Kb=Kb, tri=tri, p=p)
