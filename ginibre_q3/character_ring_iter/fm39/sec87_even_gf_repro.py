"""FM-SEC87 (luna_max_neptune): EVEN-form label generating functions at fixed length; numerator formula; L<=5 edge-cone certificates."""
from itertools import combinations, product
from collections import defaultdict

def comps(total, m, prefix=()):
    if m == 1:
        yield prefix + (total,)
        return
    for k in range(total + 1):
        yield from comps(total - k, m - 1, prefix + (k,))

def fusion_zero(labels):
    # Multiplicity of the trivial SU(2) representation.
    state = {0: 1}
    for n in labels:
        nxt = defaultdict(int)
        for j, mult in state.items():
            for k in range(abs(j - n), j + n + 1, 2):
                nxt[k] += mult
        state = nxt
    return state.get(0, 0)

def poly_add(P, Q, scale=1):
    out = defaultdict(int, P)
    for a, c in Q.items():
        out[a] += scale * c
    return {a: c for a, c in out.items() if c}

def poly_mul(P, Q, cap=None):
    out = defaultdict(int)
    for a, ca in P.items():
        for b, cb in Q.items():
            z = tuple(x + y for x, y in zip(a, b))
            if cap is None or sum(z) <= cap:
                out[z] += ca * cb
    return {a: c for a, c in out.items() if c}

def edge_factor(m, i, j):
    e = [0] * m
    e[i] = e[j] = 1
    return {(0,) * m: 1, tuple(e): -1}

def K_from_fusion(m):
    if m == 0:
        return {(): 1}
    if m <= 3:
        return {(0,) * m: 1}

    # The multigraded Pluecker numerator has total degree at most m(m-3).
    bound = m * (m - 3)
    H = {}
    for degree in range(bound + 1):
        for a in comps(degree, m):
            value = fusion_zero(a)
            if value:
                H[a] = value

    P = H
    for i, j in combinations(range(m), 2):
        P = poly_mul(P, edge_factor(m, i, j), cap=bound)
    return P

def embed(P, positions, L):
    out = {}
    for a, c in P.items():
        b = [0] * L
        for k, i in enumerate(positions):
            b[i] = a[k]
        out[tuple(b)] = c
    return out

def numerator(L, t):
    K = {m: K_from_fusion(m) for m in range(L + 1)}
    T = set(range(t))
    out = {}
    for mask in range(1 << L):
        A = tuple(i for i in range(L) if (mask >> i) & 1)
        B = tuple(i for i in range(L) if not ((mask >> i) & 1))
        term = poly_mul(embed(K[len(A)], A, L), embed(K[len(B)], B, L))
        for i in A:
            for j in B:
                term = poly_mul(term, edge_factor(L, i, j))
        sign = (-1) ** len(T.intersection(B))
        out = poly_add(out, term, sign)
    return out

def inverse_denominator_box(L, bound):
    # Coefficients of prod_(i<j) (1-x_i*x_j)^(-1), truncated to the box.
    state = {(0,) * L: 1}
    for i, j in combinations(range(L), 2):
        nxt = defaultdict(int)
        for a, c in state.items():
            k = 0
            while a[i] + k <= bound and a[j] + k <= bound:
                b = list(a)
                b[i] += k
                b[j] += k
                nxt[tuple(b)] += c
                k += 1
        state = {a: c for a, c in nxt.items() if c}
    return state

def direct_value(labels, T, cache):
    total = 0
    L = len(labels)
    for mask in range(1 << L):
        A = tuple(sorted(labels[i] for i in range(L) if (mask >> i) & 1))
        B = tuple(sorted(labels[i] for i in range(L) if not ((mask >> i) & 1)))
        if A not in cache:
            cache[A] = fusion_zero(A)
        if B not in cache:
            cache[B] = fusion_zero(B)
        sign = (-1) ** sum(i in T for i in range(L) if not ((mask >> i) & 1))
        total += sign * cache[A] * cache[B]
    return total

def from_numerator(labels, N, denominator):
    total = 0
    for a, c in N.items():
        b = tuple(n - k for n, k in zip(labels, a))
        if min(b) >= 0:
            total += c * denominator.get(b, 0)
    return total

for L, box in ((4, 3), (5, 3), (6, 2)):
    denominator = inverse_denominator_box(L, box)
    for t in range(0, L + 1, 2):
        T = set(range(t))
        N = numerator(L, t)
        cache = {}
        mismatches = negatives = extra_zeros = 0
        least = None
        for labels in product(range(box + 1), repeat=L):
            direct = direct_value(labels, T, cache)
            rational = from_numerator(labels, N, denominator)
            mismatches += (direct != rational)
            negatives += (direct < 0)
            least = direct if least is None else min(least, direct)
            total = sum(labels)
            largest = max(labels)
            expected_zero = (
                total % 2 == 1
                or largest > total - largest
                or any(labels[i] == 0 for i in T)
            )
            extra_zeros += (not expected_zero and direct == 0)
        print(L, t, box, mismatches, negatives, extra_zeros, least)

for T in ({0, 1}, {0, 1, 2, 3}):
    cache = {}
    for labels in (
        (0, 0, 0, 0), (1, 1, 1, 1), (2, 2, 2, 2),
        (2, 2, 2, 3), (5, 1, 1, 1)
    ):
        print(tuple(sorted(T)), labels, direct_value(labels, T, cache))


# ---- block ----
from itertools import combinations, permutations
from collections import defaultdict
from fractions import Fraction
from math import lcm
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

def stabilizer(L, t):
    out = []
    for p in permutations(range(t)):
        for q in permutations(range(t, L)):
            g = list(range(L))
            for i, j in zip(range(t), p):
                g[i] = j
            for i, j in zip(range(t, L), q):
                g[i] = j
            out.append(tuple(g))
    return out

def transform_m(m, g):
    out = [0] * len(m)
    for i in range(len(m)):
        out[g[i]] = m[i]
    return tuple(out)

def transform_F(F, g):
    return tuple(sorted(tuple(sorted((g[i], g[j]))) for i, j in F))

def transform_pair(m, F, g):
    return transform_m(m, g), transform_F(F, g)

def canon_m(m, group):
    return min(transform_m(m, g) for g in group)

def all_leq(M, L):
    for degree in range(M + 1):
        yield from comps(degree, L)

def certify(L, t, M):
    N = numerator(L, t)
    G = stabilizer(L, t)
    edges = list(combinations(range(L), 2))

    # Canonical edge-set representatives under the stabilizer.
    F_reps = set()
    for mask in range(1 << len(edges)):
        F = tuple(edges[k] for k in range(len(edges)) if (mask >> k) & 1)
        F_reps.add(min(transform_F(F, g) for g in G))

    # Canonical pairs (m,F); Aut(F) acts on the monomial shift.
    reps = set()
    for F in F_reps:
        if 2 * len(F) > M:
            continue
        aut = [g for g in G if transform_F(F, g) == F]
        for m in all_leq(M - 2 * len(F), L):
            reps.add((canon_m(m, aut), F))
    reps = sorted(reps)

    rows = sorted({canon_m(a, G) for a in all_leq(M, L)})
    row_index = {a: i for i, a in enumerate(rows)}

    rr, cc, vv = [], [], []
    for j, (m, F) in enumerate(reps):
        orbit = {transform_pair(m, F, g) for g in G}
        P = defaultdict(int)
        for mm, FF in orbit:
            term = {mm: 1}
            for i, k in FF:
                term = poly_mul(term, edge_factor(L, i, k), cap=M)
            for a, c in term.items():
                P[a] += c

        # Orbit-invariant columns: use the coefficient at each canonical row.
        for a, c in P.items():
            if c and a == canon_m(a, G):
                rr.append(row_index[a])
                cc.append(j)
                vv.append(float(c))

    A = coo_matrix((vv, (rr, cc)), shape=(len(rows), len(reps))).tocsr()
    b = np.array([N.get(a, 0) for a in rows], dtype=float)
    result = linprog(
        np.ones(len(reps)), A_eq=A, b_eq=b,
        bounds=(0, None), method="highs"
    )
    assert result.success, result.message

    cert = [
        (Fraction(float(result.x[j])).limit_denominator(1_000_000), reps[j])
        for j in range(len(reps)) if result.x[j] > 1e-8
    ]

    # Exact expansion of the orbit certificate.
    total = {}
    for coefficient, (m, F) in cert:
        orbit = {transform_pair(m, F, g) for g in G}
        for mm, FF in orbit:
            term = {mm: Fraction(1)}
            for i, k in FF:
                term = poly_mul(term, edge_factor(L, i, k))
            total = poly_add(total, term, coefficient)

    assert total == N
    denominator = lcm(*(c.denominator for c, _ in cert))
    print("certificate", L, t, M, len(cert), denominator, "exact")
    for coefficient, (m, F) in cert:
        print(coefficient, m, F)

for args in ((4, 0, 10), (4, 2, 10), (4, 4, 10),
             (5, 0, 14), (5, 2, 16), (5, 4, 16)):
    certify(*args)


# ---- block ----
from math import comb
L, M = 6, 22
print(sum(comb(L * (L - 1) // 2, f) * comb(M - 2*f + L, L)
          for f in range(M // 2 + 1)))
