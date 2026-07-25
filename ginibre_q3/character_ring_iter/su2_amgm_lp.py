#!/usr/bin/env python3
"""Small dense two-phase simplex, and the AM-GM allocation LP built on it.

scipy is not available in this environment, so this is a self-contained
solver.  It runs in floating point and is used only to *propose* an
allocation.  Proposals are rounded to a fixed rational denominator and then
re-checked in interval arithmetic by su2_amgm_certify; no result from this
file is accepted as proof evidence.

Bland's rule is used throughout, so the method terminates rather than cycling.
Speed is not the concern here: the largest instance has fifteen terms.
"""

from __future__ import annotations

import numpy as np

TOL = 1e-9


def _simplex(c, A, b, n_real):
    """maximize c.x subject to A x <= b, x >= 0, via two-phase simplex.

    b may contain negative entries.  Returns (status, x).
    """
    m, n = A.shape
    A = A.copy().astype(float)
    b = b.copy().astype(float)

    # Introduce slacks FIRST, turning A x <= b into the equality
    # A x + s = b with s >= 0.  Only then may a row be negated, because
    # negating an equality is sound whereas negating an inequality reverses
    # its direction and silently discards the constraint.
    T = np.zeros((m + 1, n + m + m + 1))
    T[:m, :n] = A
    T[:m, n:n + m] = np.eye(m)
    T[:m, -1] = b
    for i in range(m):
        if T[i, -1] < 0:
            T[i, :n + m] *= -1.0
            T[i, -1] *= -1.0
    # Artificials give the starting basis.
    T[:m, n + m:n + 2 * m] = np.eye(m)
    # objective: minimise sum of artificials  ->  maximise -sum
    T[m, n + m:n + 2 * m] = 1.0
    basis = list(range(n + m, n + 2 * m))
    for i in range(m):
        T[m, :] -= T[i, :]

    def pivot(T, basis, allowed):
        # Dantzig's rule (steepest reduced cost) is far faster in practice but
        # can cycle on degenerate vertices, which these transportation-shaped
        # problems produce in quantity.  Fall back to Bland's rule, which is
        # slow but provably terminating, once an instance looks degenerate.
        it, bland_after = 0, 4 * (T.shape[1] + T.shape[0])
        cap = 200 * (T.shape[1] + T.shape[0])
        while True:
            it += 1
            if it > cap:
                return False           # abandon a pathological instance
            col = -1
            if it < bland_after:
                cand = [j for j in allowed if T[-1, j] < -TOL]
                if cand:
                    col = min(cand, key=lambda j: T[-1, j])
            else:
                for j in allowed:
                    if T[-1, j] < -TOL:
                        col = j
                        break
            if col < 0:
                return True
            row, best = -1, None
            for i in range(len(basis)):
                if T[i, col] > TOL:
                    r = T[i, -1] / T[i, col]
                    if best is None or r < best - TOL or (
                        abs(r - best) <= TOL and basis[i] < basis[row]
                    ):
                        best, row = r, i
            if row < 0:
                return False           # unbounded
            T[row, :] /= T[row, col]
            for i in range(T.shape[0]):
                if i != row and abs(T[i, col]) > 0:
                    T[i, :] -= T[i, col] * T[row, :]
            basis[row] = col

    allowed_p1 = list(range(n + 2 * m))
    if not pivot(T, basis, allowed_p1):
        return "unbounded", None
    if T[-1, -1] > TOL:
        return "infeasible", None

    # Drive any artificial still basic out of the basis.  Leaving one in place
    # lets Phase II pivot around a row it cannot represent, which silently
    # returns a point violating that row's constraint.
    drop = []
    for i in range(len(basis)):
        if basis[i] >= n + m:
            piv = -1
            for j in range(n + m):
                if abs(T[i, j]) > TOL:
                    piv = j
                    break
            if piv < 0:
                drop.append(i)          # redundant row
                continue
            T[i, :] /= T[i, piv]
            for r in range(T.shape[0]):
                if r != i and abs(T[r, piv]) > 0:
                    T[r, :] -= T[r, piv] * T[i, :]
            basis[i] = piv
    if drop:
        keep = [i for i in range(len(basis)) if i not in drop]
        T = np.vstack([T[keep, :], T[-1:, :]])
        basis = [basis[i] for i in keep]

    # Phase II: drop artificials, restore the real objective.
    T2 = np.hstack([T[:, :n + m], T[:, -1:]])
    T2[-1, :] = 0.0
    for j in range(n):
        T2[-1, j] = -c[j]
    for i, bi in enumerate(basis):
        if bi < n + m and abs(T2[-1, bi]) > 0:
            T2[-1, :] -= T2[-1, bi] * T2[i, :]
    if not pivot(T2, basis, list(range(n + m))):
        return "unbounded", None

    x = np.zeros(n + m)
    for i, bi in enumerate(basis):
        if bi < n + m:
            x[bi] = T2[i, -1]
    return "optimal", x[:n_real]


def propose_allocation(log_lam, log_c, signs, free_count):
    """Propose an AM-GM allocation maximising the worst-case slack.

    log_lam[t][l], log_c[t] describe the terms; signs[t] is +1 or -1.

    Variables are alpha[p][x] for positive p and negative x, plus a scalar
    margin delta.  Returns (delta, alpha) or (None, None) when infeasible.
    """
    pos = [t for t in range(len(signs)) if signs[t] > 0]
    neg = [t for t in range(len(signs)) if signs[t] < 0]
    if not neg:
        return None, None
    if not pos:
        return None, None
    P, X, L = len(pos), len(neg), free_count
    nv = P * X + 1                      # alpha entries then delta
    dj = P * X

    def a_idx(pi, xi):
        return pi * X + xi

    rows, rhs = [], []

    # (N) sum_p alpha[p][x] = 1, as two inequalities.
    for xi in range(X):
        r = np.zeros(nv)
        for pi in range(P):
            r[a_idx(pi, xi)] = 1.0
        rows.append(r.copy()); rhs.append(1.0)
        rows.append(-r); rhs.append(-1.0)

    # (C) capacity: sum_x alpha[p][x] <= 1.
    for pi in range(P):
        r = np.zeros(nv)
        for xi in range(X):
            r[a_idx(pi, xi)] = 1.0
        rows.append(r); rhs.append(1.0)

    # (D) domination per free coordinate, and (K) the constant condition.
    #     -sum_p alpha[p][x] g[p] + delta <= -h[x]
    for xi, x in enumerate(neg):
        for l in range(L):
            r = np.zeros(nv)
            for pi, p in enumerate(pos):
                r[a_idx(pi, xi)] = -log_lam[p][l]
            r[dj] = 1.0
            rows.append(r); rhs.append(-log_lam[x][l])
        r = np.zeros(nv)
        for pi, p in enumerate(pos):
            r[a_idx(pi, xi)] = -log_c[p]
        r[dj] = 1.0
        rows.append(r); rhs.append(-log_c[x])

    c = np.zeros(nv)
    c[dj] = 1.0                          # maximise the margin
    Am, bm = np.array(rows), np.array(rhs)
    status, sol = _simplex(c, Am, bm, nv)
    if status != "optimal" or sol is None:
        return None, None

    # Never trust the solver.  A proposal that violates its own constraints is
    # discarded here rather than being rounded and handed to the verifier,
    # where it would look like a certificate failure instead of a solver bug.
    if np.any(sol < -1e-7) or np.max(Am.dot(sol) - bm) > 1e-6:
        return None, None

    delta = sol[dj]
    alpha = [[float(sol[a_idx(pi, xi)]) for xi in range(X)] for pi in range(P)]
    return delta, (pos, neg, alpha)
