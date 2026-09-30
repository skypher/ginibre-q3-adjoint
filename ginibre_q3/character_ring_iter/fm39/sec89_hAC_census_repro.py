"""FM-SEC89 (luna_max_mars): H_AC dictionary census; dictionary failures at (1,1,1,1,2,2) and (2^6,4);
exact H_AC certificate at (2^6,4) with subspace-orbit indicators; gap-two insertion law; path-factor tables."""
from fractions import Fraction as Q
from itertools import permutations
from math import comb
import random

def add_row(row, n):
    out = {}
    for j, v in row.items():
        for k in range(abs(j - n), j + n + 1, 2):
            out[k] = out.get(k, 0) + v
    return out

def subset_rows(labels):
    rows = [{0: 1}]
    for n in labels:
        old = rows[:]
        for row in old:
            rows.append(add_row(row, n))
    return rows

def target(labels):
    d = len(labels) - 1
    rows = subset_rows(labels[:-1])
    full = (1 << d) - 1
    return [
        rows[x].get(0, 0) * rows[full ^ x].get(labels[-1], 0)
        for x in range(1 << d)
    ]

def square(p):
    out = [Q(0)] * len(p)
    for x, a in enumerate(p):
        if a:
            for y, b in enumerate(p):
                if b:
                    out[x ^ y] += a * b
    return out

def path_factor(labels, j, k):
    rows = subset_rows(labels)
    d = len(labels) - 1
    full = (1 << len(labels)) - 1
    p = []
    for x in range(1 << d):
        a = rows[x].get(j, 0) * rows[full ^ x].get(k, 0)
        b = rows[x].get(k, 0) * rows[full ^ x].get(j, 0)
        p.append(a if j == k else a + b)
    return p

# Exact path-factor re-expansions, including the first profile-sphere failure.
boundary = [
    ((1,) * 8 + (6,), {(0, 14): Q(3), (7, 7): Q(1, 2)}),
    ((1,) * 10 + (2,), {
        (0, 12): Q(1567, 92), (1, 11): Q(157, 92),
        (2, 10): Q(33, 92), (3, 9): Q(59, 460),
        (4, 8): Q(41, 460),
    }),
    ((1,) * 6 + (2,), {
        (0, 8): Q(2), (1, 7): Q(1, 2), (2, 6): Q(1, 4),
    }),
    ((3, 3, 3, 3, 2), {(0, 14): Q(7), (3, 11): Q(1, 2)}),
    ((1, 1, 1, 1, 2, 2), {
        (0, 8): Q(2), (1, 7): Q(1, 2), (2, 6): Q(1, 4),
    }),
]
for labels, terms in boundary:
    rhs = [Q(0)] * (1 << (len(labels) - 1))
    for (j, k), a in terms.items():
        rhs = [
            v + a * w for v, w in zip(rhs, square(path_factor(labels, j, k)))
        ]
    assert rhs == list(map(Q, target(labels)))

# Profile-sphere obstruction at (1^4, 2^2).
labels = (1, 1, 1, 1, 2, 2)
totals = tuple(labels.count(v) for v in sorted(set(labels)))

def profile(x):
    return tuple(
        sum(bool(x >> i & 1) for i, v in enumerate(labels[:-1]) if v == label)
        for label in sorted(set(labels))
    )

def canonical_profile(w):
    complement = tuple(totals[i] - w[i] for i in range(len(totals)))
    return min(w, complement)

x = (1 << 0) | (1 << 1) | (1 << 4)
f = target(labels)
profiles = sorted({
    canonical_profile(profile(y))
    for y in range(1 << (len(labels) - 1))
})
assert profile(x) == (2, 1) and f[x] == 1
for w in profiles:
    p = [
        int(canonical_profile(profile(y)) == w)
        for y in range(1 << (len(labels) - 1))
    ]
    assert square(p)[x] == 0

# Farkas witness for spheres + all endpoint paths + averaged point-pair factors
# at (2^6, 4). Rows are the per-mask values for weights 0 through 6.
labels = (2, 2, 2, 2, 2, 2, 4)
d = 6
T = sum(labels)
ff = target(labels)
reps = [(1 << k) - 1 for k in range(d + 1)]
b = [Q(ff[x]) for x in reps]
columns = []

for w in range(4):
    p = [int(x.bit_count() in {w, d - w}) for x in range(1 << d)]
    c = square(p)
    columns.append([c[x] for x in reps])

for j in range(T + 1):
    for k in range(j, T + 1):
        p = path_factor(labels, j, k)
        if any(p):
            c = square(p)
            col = [c[x] for x in reps]
            if any(col):
                columns.append(col)

for w in (2, 3, 4):
    columns.append([
        Q(2) if k == 0 else Q(2, comb(d, w)) if k == w else Q(0)
        for k in range(d + 1)
    ])

Y = (2270, 61188, -6810, -45400, 3405, 0, 1135)

def dot(u, v):
    return sum(Q(a) * Q(b) for a, b in zip(u, v))

assert dot(Y, b) == -30645
assert all(dot(Y, col) >= 0 for col in columns)

# H_AC certificate for (2^6, 4) using two subspace orbits.
def span(gens):
    H = {0}
    for g in gens:
        H |= {x ^ g for x in tuple(H)}
    return H

def perm_mask(x, pi):
    return sum(1 << pi[i] for i in range(6) if x >> i & 1)

def subspace_orbit(gens):
    H = span(gens)
    return {
        tuple(sorted(perm_mask(x, pi) for x in H))
        for pi in permutations(range(6))
    }

O1 = subspace_orbit((3, 5, 9))
O2 = subspace_orbit((7, 25, 42))
assert (len(O1), len(O2)) == (15, 30)

rhs = [Q(15 if x == 0 else 0) for x in range(1 << d)]
for orbit, a in ((O1, Q(1, 8)), (O2, Q(1, 24))):
    for H in orbit:
        p = [int(x in H) for x in range(1 << d)]
        rhs = [u + a * v for u, v in zip(rhs, square(p))]
assert rhs == list(map(Q, ff))

# Exact insertion-rule check on 200 random lists, with labels <= 6.
rng = random.Random(8901)
seen = 0
max_inserted = 0
while seen < 200:
    m = rng.randint(2, 11)
    mu = tuple(rng.randint(1, 6) for _ in range(m))
    n = sum(mu) - 2
    if n < 1:
        continue

    p = [Q(0)] * (1 << m)
    for i, z in enumerate(mu):
        if z == 1:
            p[1 << i] = Q(1)

    t = mu.count(1)
    residual = Q(m - 1) - Q(t, 2)
    rhs = [v / 2 for v in square(p)]
    rhs[0] += residual
    assert rhs == list(map(Q, target(mu + (n,))))

    seen += 1
    max_inserted = max(max_inserted, n)

assert seen == 200 and max_inserted == 47
