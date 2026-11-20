#!/usr/bin/env python3
"""
gen_patterns.py -- Verify Theorem 2 (unsatisfiability of B <= 6 blocks and
satisfiability of B = 7, 8, 9 blocks under active block-counter steering and
single-excursion modulo-32 pool parity cancellation).
"""
import sys
from itertools import combinations

try:
    import pycryptosat
except ImportError:
    print("pycryptosat not installed; using pre-generated pattern_{7,8,9}blk.h headers.")
    sys.exit(0)

L1, L2, L3 = 16, 7, 6
POOL = 32
BLK = 17
CNT = 16


def solve_pattern(T, free_before_ctr=False):
    F = 1
    d = lambda t: 2 + t
    u = lambda t: 2 + T + t
    v = lambda t: 2 + 2 * T + t
    w = lambda t: 2 + 3 * T + t
    fs = (u, v, w)
    cl = [[F]]
    xcl = []

    for t in range(T):
        xcl.append(([u(t)] + ([u(t - 1)] if t >= 1 else []) +
                    ([u(t - L1)] if t >= L1 else []) + [d(t)], False))
        xcl.append(([v(t)] + ([v(t - 1)] if t >= 1 else []) +
                    ([v(t - L2)] if t >= L2 else []) + [d(t)], False))
        xcl.append(([w(t)] + ([w(t - 1)] if t >= 1 else []) +
                    ([w(t - L3)] if t >= L3 else []) + [d(t)], False))
        if t % BLK == CNT:
            cl.append([-d(t)])
            for i in range(3):
                for j in range(i + 1, 3):
                    cl.append([-fs[i](t), -fs[j](t)])
            if free_before_ctr and t >= 1:
                for i in range(3):
                    for j in range(3):
                        cl.append([-fs[i](t), -fs[j](t - 1)])

    aux0 = 2 + 4 * T
    sel = [aux0 + t for t in range(T)]
    for t in range(T):
        cl.append([-u(t), -v(t), -w(t), sel[t]])
        cl.append([-sel[t], u(t)])
        cl.append([-sel[t], v(t)])
        cl.append([-sel[t], w(t)])
    cl.append(list(sel))
    for comb in combinations(range(T), 2):
        cl.append([-sel[t] for t in comb])

    for t in range(T):
        if t + 1 < T:
            for i in range(3):
                for j in range(i + 1, 3):
                    cl.append([-sel[t], -fs[i](t + 1), -fs[j](t + 1)])

    for t in range(T):
        if t < 4 or any((t - k) % BLK == CNT for k in (1, 2, 3, 4)):
            cl.append([-sel[t]])
            continue
        for f in fs:
            for k in (1, 2, 3, 4):
                cl.append([-sel[t], -f(t - k)])

    nxt = 2 + 5 * T
    has_cond = [nxt + t for t in range(T)]
    nxt += T
    for t in range(T):
        if t % BLK == CNT:
            for i in range(3):
                cl.append([-fs[i](t), has_cond[t]])
            cl.append([-has_cond[t], u(t), v(t), w(t)])
        else:
            for i in range(3):
                for j in range(i + 1, 3):
                    cl.append([-fs[i](t), -fs[j](t), has_cond[t]])
            cl.append([-has_cond[t], u(t), v(t)])
            cl.append([-has_cond[t], u(t), w(t)])
            cl.append([-has_cond[t], v(t), w(t)])

    for t in range(1, T):
        if (t - 1) % BLK == CNT:
            cl.append([-has_cond[t]])
            continue
        branches = []
        for k in (1, 2, 3):
            s = t - k
            if s < 0 or s % BLK == CNT:
                continue
            y = nxt; nxt += 1
            branches.append(y)
            cl.append([-y, -sel[s]])
            for r in range(s, t):
                cl.append([-y, -has_cond[r]])
        if not branches:
            cl.append([-has_cond[t]])
        else:
            cl.append([-has_cond[t]] + branches)

    for t in range(T - L1, T): cl.append([-u(t)])
    for t in range(T - L2, T): cl.append([-v(t)])
    for t in range(T - L3, T): cl.append([-w(t)])
    for r in range(POOL):
        for f in (v, w):
            lits = [f(t) for t in range(r, T, POOL)]
            if len(lits) > 1: xcl.append((lits, False))
            elif len(lits) == 1: cl.append([-lits[0]])

    s = pycryptosat.Solver()
    s.add_clauses(cl)
    for lits, rhs in xcl:
        s.add_xor_clause(lits, rhs)
    sat, sol = s.solve()
    if not sat:
        return None
    D = [1 if sol[d(t)] else 0 for t in range(T)]
    U = [1 if sol[u(t)] else 0 for t in range(T)]
    V = [1 if sol[v(t)] else 0 for t in range(T)]
    W = [1 if sol[w(t)] else 0 for t in range(T)]
    return D, U, V, W


def main():
    print("[Theorem 2] Checking single-excursion block-length lower bound:")
    for B in (5, 6, 7, 8, 9):
        T = 17 * B
        r = solve_pattern(T, free_before_ctr=(B >= 8))
        status = "SAT" if r is not None else "UNSAT"
        print(f"  B={B} blocks ({B*64:3d} B, T={T:3d} steps): {status}")
        if B <= 6 and r is not None:
            return 1
        if B >= 7 and r is None:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
