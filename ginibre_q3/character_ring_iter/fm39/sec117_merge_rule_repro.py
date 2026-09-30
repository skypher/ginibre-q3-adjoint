"""FM-SEC117 (luna_max_pluto): Bernstein-step merge rule; failure at (1,2,1,2) j=1->2; origin budget at (1^6); transition scan."""
from functools import lru_cache
from itertools import product, combinations
from collections import Counter
from fractions import Fraction as Q
from math import comb

@lru_cache(None)
def nc_matchings(legs, owner):
    if not legs: return ((),)
    if len(legs) % 2: return ()
    a, out = legs[0], []
    for j in range(1, len(legs), 2):
        if owner[a] == owner[legs[j]]: continue
        for inner in nc_matchings(legs[1:j], owner):
            for outer in nc_matchings(legs[j+1:], owner):
                out.append(((a, legs[j]),) + inner + outer)
    return tuple(out)

def crosses(e, f):
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b

def bernstein_tables(labels):
    owner = tuple(i for i, n in enumerate(labels) for _ in range(n))
    L = len(labels)
    power = []
    degree = 0
    for S in range(1 << L):
        left = tuple(k for k, o in enumerate(owner) if (S >> o) & 1)
        right = tuple(k for k, o in enumerate(owner)
                      if not ((S >> o) & 1))
        counts = Counter(
            sum(crosses(e, f) for e in M for f in N)
            for M in nc_matchings(left, owner)
            for N in nc_matchings(right, owner)
        )
        power.append(counts)
        if counts:
            degree = max(degree, max(counts))
    G = 1 << (L - 1)
    tables = [
        [
            sum(Q(power[x].get(k, 0) * comb(j, k), comb(degree, k))
                for k in range(j + 1))
            for x in range(G)
        ]
        for j in range(degree + 1)
    ]
    return degree, tables

def conv(p):
    return [sum(p[y] * p[x ^ y] for y in range(len(p)))
            for x in range(len(p))]

def wht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            for j in range(i, i + h):
                x, y = v[j], v[j + h]
                v[j], v[j + h] = x + y, x - y
        h *= 2
    return v

def indicator(G, H):
    return [Q(int(x in H)) for x in range(G)]

# Exact increment scan: L=2..6, labels 1..3, even total.
word_count = step_count = positive = 0
per_level = []
budget_fail = []
for L in range(2, 7):
    words = steps = 0
    for labels in product((1, 2, 3), repeat=L):
        if sum(labels) % 2:
            continue
        words += 1
        degree, gs = bernstein_tables(labels)
        for j, (g, h) in enumerate(zip(gs, gs[1:])):
            delta = [b - a for a, b in zip(g, h)]
            steps += 1
            step_count += 1
            assert all(v >= 0 for v in delta)
            assert delta[0] == 0
            positive += sum(v > 0 for v in delta)
            if sum(delta) > g[0]:
                budget_fail.append((labels, j, g[0], sum(delta)))
    word_count += words
    per_level.append((L, words, steps))

print("census per level", per_level)
print("census totals", word_count, step_count,
      "negative increments", 0, "origin changes", 0,
      "positive increment entries", positive)
print("independent two-point origin-budget failures",
      len(budget_fail), "first", budget_fail[0])

# Exact active four-label lists, in the order searched.
active4 = []
for labels in product((1, 2, 3), repeat=4):
    if sum(labels) % 2 == 0:
        degree, _ = bernstein_tables(labels)
        if degree:
            active4.append((labels, degree))
print("L4 lists with a nonzero step", active4)

# The preceding active transition, (1^4), is one subgroup merge.
_, g4 = bernstein_tables((1, 1, 1, 1))
G = 8
H = frozenset((0, 3))
K = frozenset((0, 6))
J = frozenset((0, 3, 5, 6))
assert g4[0] == [
    a + b for a, b in zip(indicator(G, H), indicator(G, K))
]
assert g4[1] == [
    a + int(x == 0) for x, a in enumerate(indicator(G, J))
]
assert all(g4[1][x] - g4[0][x] == int(x == 5) for x in range(G))
print("prior merge (1^4): g0=1_<3>+1_<6>; E={5}; "
      "g1=1_<3,6>+delta0")

# First local subgroup-merge obstruction.
labels = (1, 2, 1, 2)
degree, gs = bernstein_tables(labels)
assert degree == 2
assert gs[0] == gs[1] == [Q(2)] + [Q(0)] * 7
assert gs[2] == [Q(2), Q(0), Q(0), Q(0),
                 Q(0), Q(1), Q(0), Q(0)]
print("first obstruction", labels, "step 1->2")
print("g1", list(map(str, gs[1])))
print("delta", list(map(str, [b-a for a, b in zip(gs[1], gs[2])])))
print("g2", list(map(str, gs[2])))

# Exact two-point autocorrelation rewrite of that target.
p = [Q(0)] * 8
p[1] = p[4] = Q(1)
repair = [Q(int(x == 0)) + conv(p)[x] / 2 for x in range(8)]
assert repair == gs[2]
print("two-point AC rewrite", list(map(str, repair)))

# Grouped increment at (1^6), then the residual Walsh check.
_, all1_6 = bernstein_tables((1,) * 6)
delta = [b - a for a, b in zip(all1_6[0], all1_6[1])]
assert [x for x, v in enumerate(delta) if v] == [5, 10, 17, 20, 23, 29]
pA = [Q(0)] * 32
pB = [Q(0)] * 32
for x in (0, 5, 17):
    pA[x] = Q(1)
for x in (0, 10, 23):
    pB[x] = Q(1)
grouped = [(conv(pA)[x] + conv(pB)[x]) / 2 for x in range(32)]
assert grouped == [delta[x] + 3 * int(x == 0) for x in range(32)]
residual = [all1_6[0][x] - 3 * int(x == 0) for x in range(32)]
assert min(wht(residual)) == -3
assert sum(v < 0 for v in wht(residual)) == 12
print("(1^6) grouped increment identity exact; g0(0)=",
      all1_6[0][0])
print("(1^6) residual g0-3delta0: min Walsh",
      min(wht(residual)), "negative entries",
      sum(v < 0 for v in wht(residual)))

# Eight-point sphere identity and the terminal (1^8,6) table.
N = 8
G = 1 << N
p1 = [Q(int(x.bit_count() == 1)) for x in range(G)]
conv1 = conv(p1)
line_sum = [
    sum(int(x == 0 or x == ((1 << i) ^ (1 << j)))
        for i, j in combinations(range(N), 2))
    for x in range(G)
]
assert line_sum == [
    conv1[x] / 2 + 24 * int(x == 0) for x in range(G)
]
_, endpoint = bernstein_tables((1,) * 8 + (6,))
assert endpoint[6] == [
    3 * int(x == 0) + conv1[x] / 2 for x in range(G)
]
print("sphere identity exact: sum pair-lines = "
      "half p1*p1 + 24 delta0")
print("terminal table exact: g6 = 3 delta0 + half p1*p1; "
      "g6(0)=", endpoint[6][0])
