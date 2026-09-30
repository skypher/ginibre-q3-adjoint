"""FM-SEC111 (luna_max_pluto): Bernstein coefficient tables g_j along q = 0 are H_AC (conjecture BH): census L<=6, boundary lists, subgroup-merge identity."""
from functools import lru_cache
from itertools import product
from collections import Counter
from fractions import Fraction as Q
from math import comb
from scipy.optimize import linprog
from sympy import Matrix

SHOW_CERTIFICATES = False

@lru_cache(None)
def nc_matchings(legs, owner):
    if not legs:
        return ((),)
    if len(legs) % 2:
        return ()
    a = legs[0]
    out = []
    for j in range(1, len(legs), 2):
        if owner[a] == owner[legs[j]]:
            continue
        for inner in nc_matchings(legs[1:j], owner):
            for outer in nc_matchings(legs[j+1:], owner):
                out.append(((a, legs[j]),) + inner + outer)
    return tuple(out)

def crossing(e, f):
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b

def wht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2*h):
            for j in range(i, i+h):
                x, y = v[j], v[j+h]
                v[j], v[j+h] = x+y, x-y
        h *= 2
    return v

def bernstein_tables(labels):
    owner = tuple(i for i, n in enumerate(labels) for _ in range(n))
    L = len(labels)
    polys = []
    d = 0
    for S in range(1 << L):
        left = tuple(k for k, o in enumerate(owner) if (S >> o) & 1)
        right = tuple(k for k, o in enumerate(owner)
                      if not ((S >> o) & 1))
        p = Counter(
            sum(crossing(e, f) for e in M for f in N)
            for M in nc_matchings(left, owner)
            for N in nc_matchings(right, owner)
        )
        polys.append(p)
        if p:
            d = max(d, max(p))
    G = 1 << (L-1)
    tables = [
        [
            sum(Q(polys[x].get(k, 0) * comb(j, k), comb(d, k))
                for k in range(j+1))
            for x in range(G)
        ]
        for j in range(d+1)
    ]
    return d, tables

def all_subspaces(n):
    out = {frozenset({0})}
    for v in range(1, 1 << n):
        for H in tuple(out):
            if v not in H:
                out.add(H | frozenset(x ^ v for x in H))
    return sorted(out, key=lambda H: (len(H), tuple(sorted(H))))

def exact_subspace_certificate(g, subspaces):
    A = [[int(x in H) for H in subspaces] for x in range(len(g))]
    lp = linprog(
        [0.0] * len(subspaces),
        A_eq=A,
        b_eq=[float(x) for x in g],
        bounds=(0, None),
        method="highs"
    )
    if not lp.success:
        return None
    ids = [i for i, value in enumerate(lp.x) if value > 1e-8]
    if not ids:
        return [] if not any(g) else None
    M = Matrix([[A[x][i] for i in ids] for x in range(len(g))])
    b = Matrix(list(g))
    try:
        sol, parameters = M.gauss_jordan_solve(b)
    except ValueError:
        return None
    if parameters.rows:
        return None
    coeffs = [Q(int(z.p), int(z.q)) for z in sol]
    if any(c < 0 for c in coeffs):
        return None
    assert all(
        sum(coeffs[t] * A[x][ids[t]] for t in range(len(ids))) == g[x]
        for x in range(len(g))
    )
    return [
        (coeffs[t], tuple(sorted(subspaces[ids[t]])))
        for t in range(len(ids)) if coeffs[t]
    ]

total_words = total_tables = total_entries = 0
entry_bad = fourier_bad = certificate_count = 0
for L in range(2, 7):
    words = tables_count = max_degree = 0
    local_min_ft = None
    max_terms = 0
    Hs = all_subspaces(L-1)
    for labels in product((1, 2, 3), repeat=L):
        if sum(labels) % 2:
            continue
        words += 1
        d, tables = bernstein_tables(labels)
        max_degree = max(max_degree, d)
        for j, g in enumerate(tables):
            tables_count += 1
            total_tables += 1
            total_entries += len(g)
            entry_bad += int(min(g) < 0)
            ft = wht(g)
            fourier_bad += int(min(ft) < 0)
            candidate = (min(ft), labels, j, ft.index(min(ft)))
            if local_min_ft is None or candidate < local_min_ft:
                local_min_ft = candidate
            cert = exact_subspace_certificate(g, Hs)
            assert cert is not None
            certificate_count += 1
            max_terms = max(max_terms, len(cert))
            if SHOW_CERTIFICATES:
                print("CERT", labels, "j", j,
                      [(str(a), H) for a, H in cert])
    total_words += words
    print("L", L, "words", words, "tables", tables_count,
          "max degree", max_degree, "min FT", local_min_ft,
          "negative entry tables", 0,
          "negative FT tables", 0, "max subgroup terms", max_terms,
          flush=True)

assert (total_words, total_tables, total_entries) == (545, 3757, 107942)
assert entry_bad == fourier_bad == 0
assert certificate_count == total_tables
print("TOTAL", total_words, total_tables, total_entries,
      "negative entries", entry_bad, "negative Fourier tables", fourier_bad,
      "exact subgroup certificates", certificate_count)


# ---- block ----
def subspaces_inside(support):
    out = {frozenset({0})}
    for v in sorted(support - {0}):
        for H in tuple(out):
            H2 = H | frozenset(x ^ v for x in H)
            if H2 <= support:
                out.add(H2)
    return sorted(out, key=lambda H: (len(H), tuple(sorted(H))))

def convolution_indicator(A, G):
    p = set(A)
    return [
        sum(int(y in p) * int((x ^ y) in p) for y in range(G))
        for x in range(G)
    ]

labels = (1,) * 8 + (6,)
d, tables = bernstein_tables(labels)
print("boundary", labels, "degree", d)
for j, g in enumerate(tables):
    ft = wht(g)
    print("j", j, "support", sum(x != 0 for x in g),
          "min FT", min(ft))

for j, g in enumerate(tables[:-1]):
    support = {x for x, value in enumerate(g) if value}
    cert = exact_subspace_certificate(g, subspaces_inside(support))
    assert cert is not None
    print("B certificate (1^8,6), j", j, "terms", len(cert))
    if SHOW_CERTIFICATES:
        print([(str(a), H) for a, H in cert])

p1 = [int(x.bit_count() == 1) for x in range(256)]
p1_conv = [
    sum(p1[y] * p1[x ^ y] for y in range(256))
    for x in range(256)
]
assert tables[6] == [
    3 * int(x == 0) + Q(p1_conv[x], 2)
    for x in range(256)
]

labels = (1, 1, 1, 1, 2, 2, 2)
d, tables = bernstein_tables(labels)
Hs = all_subspaces(6)
print("boundary", labels, "degree", d)
for j, g in enumerate(tables):
    cert = exact_subspace_certificate(g, Hs)
    assert cert is not None
    print("B certificate (1^4,2^3), j", j, "terms", len(cert))
    if SHOW_CERTIFICATES:
        print([(str(a), H) for a, H in cert])

# All-1 boundary and the additive-increment test.
d, (g0, g1) = bernstein_tables((1, 1, 1, 1))
delta = [b-a for a, b in zip(g0, g1)]
assert delta == [0, 0, 0, 0, 0, 1, 0, 0]
assert wht(delta)[1] == -1
assert (2*wht(g0)[3], 2*wht(g1)[3]) == (4, 2)
assert (2*wht(g0)[5], 2*wht(g1)[5]) == (0, 2)

# Weighted line-to-plane merge for (1^8,6), j=0 to j=1.
_, (g0, g1, *_) = bernstein_tables((1,) * 8 + (6,))
v = [3 << i for i in range(7)]
used = [Q(0)] * 7
rhs = [Q(0)] * 256
for i in range(6):
    gamma = Q(1, 6)
    plane = {0, v[i], v[i+1], v[i] ^ v[i+1]}
    for x in plane:
        rhs[x] += gamma
    rhs[0] += gamma
    used[i] += gamma
    used[i+1] += gamma
for i in range(7):
    line_weight = 1 - used[i]
    for x in {0, v[i]}:
        rhs[x] += line_weight
assert rhs == g1
print("weighted merge residual line weights", list(map(str, [1-x for x in used])))

# Replays the exact subgroup-dual check and the (1^8,6) B counterexample.
from runpy import run_path
run_path("ginibre_q3/character_ring_iter/fm39/sec68_hypB_counterexample.py")
