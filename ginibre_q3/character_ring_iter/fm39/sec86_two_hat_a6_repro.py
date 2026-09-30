"""FM-SEC86 (luna_max_neptune): level-1 two-hat words phi_1(h_u h_v hat S_P hat S_Q h_1^a)
for all labels at suffix a = 3..6, by positive edge-orbit certificates of the
label generating function N_a / D_4.  Run from this directory."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product

import numpy as np
import sympy as sp
from scipy.optimize import linprog
from n0check import h, shat, mul, phi

t1, t2, u, v = sp.symbols("t1 t2 u v")
X = (t1, t2, u, v)
EDGES = list(combinations(range(4), 2))
PERMS = [(0, 1, 2, 3), (1, 0, 2, 3),
         (0, 1, 3, 2), (1, 0, 3, 2)]

def catmom(m):
    if m % 2:
        return sp.Integer(0)
    n = m // 2
    return sp.binomial(2*n, n) // (n + 1)

def edge_den(ds):
    return sp.prod(1 - ds[i]*ds[j]
                   for i in range(len(ds))
                   for j in range(i + 1, len(ds)))

# J_m(S), for every subset of the four grading variables.
SUBSETS = []
for k in range(5):
    SUBSETS.extend(combinations(X, k))

J = {(m, ()): catmom(m) for m in range(9)}
for S in SUBSETS[1:]:
    n = len(S)
    if n == 1:
        J[(0, S)] = sp.Integer(1)
    elif n == 4:
        J[(0, S)] = sp.cancel((1 - sp.prod(S)) / edge_den(S))
    else:
        J[(0, S)] = sp.cancel(1 / edge_den(S))

for m in range(8):
    for S in SUBSETS[1:]:
        d = S[0]
        J[(m + 1, S)] = sp.cancel(
            ((1 + d**2)*J[(m, S)] - J[(m, S[1:])]) / d
        )

def j(m, ds):
    key = tuple(sorted(ds, key=X.index))
    return J[(m, key)]

D4 = sp.prod(1 - X[i]*X[j]
             for i in range(4)
             for j in range(i + 1, 4))

def exact_numerator(a):
    # Expand (x-y)^2 (x+y)^a and integrate each x/y monomial.
    integral_sum = 0
    for sides in product((0, 1), repeat=2):
        dx = (u, v) + tuple(t for t, side in zip((t1, t2), sides)
                            if side == 0)
        dy = (u, v) + tuple(t for t, side in zip((t1, t2), sides)
                            if side == 1)
        for k in range(a + 1):
            b = sp.binomial(a, k)
            integral_sum += b * (
                j(k + 2, dx)*j(a - k, dy)
                + j(k, dx)*j(a - k + 2, dy)
                - 2*j(k + 1, dx)*j(a - k + 1, dy)
            )
    rational = sp.cancel(integral_sum * D4 / 2)
    num, den = sp.fraction(rational)
    assert den == 1
    poly = sp.Poly(sp.expand(num), *X)
    out = {}
    for mon, coeff in poly.terms():
        assert coeff.q == 1
        out[mon] = int(coeff)
    return out

def compositions(total, n=4, prefix=()):
    if n == 1:
        yield prefix + (total,)
        return
    for k in range(total + 1):
        yield from compositions(total - k, n - 1, prefix + (k,))

def denominator_series(degree):
    series = {(0, 0, 0, 0): 1}
    for i, j0 in EDGES:
        nxt = defaultdict(int)
        for mon, c in series.items():
            for k in range((degree - sum(mon)) // 2 + 1):
                q = list(mon)
                q[i] += k
                q[j0] += k
                nxt[tuple(q)] += c
        series = dict(nxt)
    return series

def check_direct(a, N, degree=10):
    den = denominator_series(degree)
    quotient_series = defaultdict(int)
    for m, c in N.items():
        for d, b in den.items():
            q = tuple(x0 + y0 for x0, y0 in zip(m, d))
            if sum(q) <= degree:
                quotient_series[q] += c*b

    checks = 0
    for total in range(degree + 1):
        for P, Q, uu, vv in compositions(total):
            word = mul(mul(h(uu), h(vv)), mul(shat(P), shat(Q)))
            direct = phi(1, word, a)
            assert direct == quotient_series.get((P, Q, uu, vv), 0)
            checks += 1
    print("a =", a, "direct Catalan checks =", checks)

def transform(term, perm):
    m, F = term
    mt = tuple(m[perm.index(i)] for i in range(4))
    edge_id = {e: i for i, e in enumerate(EDGES)}
    Ft = tuple(sorted(
        edge_id[tuple(sorted((perm[i], perm[j])))]
        for k in F
        for i, j in [EDGES[k]]
    ))
    return mt, Ft

def expand_term(m, F):
    out = defaultdict(int)
    for mask in range(1 << len(F)):
        q = list(m)
        sign = 1
        for bit, edge in enumerate(F):
            if mask >> bit & 1:
                i, j = EDGES[edge]
                q[i] += 1
                q[j] += 1
                sign = -sign
        out[tuple(q)] += sign
    return {k: z for k, z in out.items() if z}

def orbit_certificate(N):
    degree = max(map(sum, N))
    rows = list(m for d in range(degree + 1)
                for m in compositions(d))
    row_id = {m: i for i, m in enumerate(rows)}
    target = np.array([N.get(m, 0) for m in rows], dtype=float)

    terms, columns, seen = [], [], set()
    for fmask in range(1 << len(EDGES)):
        F = tuple(i for i in range(len(EDGES)) if fmask >> i & 1)
        for d in range(degree - 2*len(F) + 1):
            for m in compositions(d):
                orbit = sorted({transform((m, F), p) for p in PERMS})
                if orbit[0] in seen:
                    continue
                seen.add(orbit[0])

                col = defaultdict(Fraction)
                for mm, FF in orbit:
                    for mon, coeff in expand_term(mm, FF).items():
                        col[mon] += Fraction(coeff, len(orbit))
                terms.append((m, F))
                columns.append(col)

    A = np.zeros((len(rows), len(columns)))
    for j, col in enumerate(columns):
        for mon, coeff in col.items():
            A[row_id[mon], j] = float(coeff)

    result = linprog(np.ones(len(columns)), A_eq=A, b_eq=target,
                     bounds=(0, None), method="highs")
    assert result.success, result.message

    chosen = [(i, x) for i, x in enumerate(result.x) if x > 1e-8]
    rational = [(i, Fraction(float(x)).limit_denominator(100000))
                for i, x in chosen]

    # Exact re-expansion; this is the certificate check.
    check = defaultdict(Fraction)
    for i, coeff in rational:
        for mon, value in columns[i].items():
            check[mon] += coeff*value
    assert all(check.get(m, Fraction(0)) == Fraction(N.get(m, 0))
               for m in set(check) | set(N))

    print("degree", degree, "numerator monomials", len(N),
          "orbit terms", len(rational), "exact re-expansion: OK")
    print("certificate rows:")
    print([(terms[i][0], terms[i][1], str(c)) for i, c in rational])

for a in (3, 4, 5, 6):
    N = exact_numerator(a)
    print("a =", a, "degree =", max(map(sum, N)),
          "numerator monomials =", len(N))
    orbit_certificate(N)
    if a in (3, 4):
        check_direct(a, N)
