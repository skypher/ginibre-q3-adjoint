from collections import defaultdict
from itertools import combinations
import random

FAMS = ("distinguished", "smallest", "duplicate", "opposite",
        "largest-smallest", "max-gap")

def canon(B):
    return tuple(sorted(B, key=lambda z: (abs(z), z)))

def entry(C, A, B):
    # Coefficient of U_A(x) U_B(y) in the product for C.
    rem = sum(abs(z) for z in C)
    d = {(0, 0): 1}
    for z in sorted(C, key=lambda q: abs(q), reverse=True):
        n, eps = abs(z), (1 if z > 0 else -1)
        rem -= n
        nd = defaultdict(int)
        for (s, t), v in d.items():
            if abs(t - B) <= rem:
                for x in range(abs(s - n), s + n + 1, 2):
                    if abs(x - A) <= rem:
                        nd[x, t] += v
            if abs(s - A) <= rem:
                for y in range(abs(t - n), t + n + 1, 2):
                    if abs(y - B) <= rem:
                        nd[s, y] += eps * v
        d = nd
    return d.get((A, B), 0)

def class_pairs(L):
    vals = []
    for z in L:
        if z not in vals:
            vals.append(z)
    return [(u, v) for i, u in enumerate(vals) for v in vals[i:]
            if u != v or L.count(u) >= 2]

def analyze(B, signed_p):
    B = canon(B)
    p = abs(signed_p)
    sigma = -1 if sum(z < 0 for z in B) % 2 else 1
    pv = sigma * p
    assert signed_p == pv
    assert len({abs(z) for z in B}) == len(B) or all(
        len({z > 0 for z in B if abs(z) == n}) == 1
        for n in {abs(z) for z in B})
    W = sum(abs(z) for z in B)
    delta = (W - p) // 2
    assert p > max(abs(z) for z in B) and p >= 6
    assert (W - p) % 2 == 0 and delta >= 8
    assert max(abs(z) for z in B) <= delta
    assert sum(abs(z) >= 3 for z in B) >= 2

    L = B + (pv,)
    pairs = class_pairs(L)
    maxgap = max(abs(abs(u) - abs(v)) for u, v in pairs)
    lo, hi = min(map(abs, L)), max(map(abs, L))
    rows = []
    for u, v in pairs:
        C = list(L)
        C.remove(u)
        C.remove(v)
        D = (1 if v > 0 else -1) * entry(C, abs(u), abs(v))
        if D < 0:
            continue
        flags = set()
        if u == pv or v == pv:
            flags.add("distinguished")
        if min(abs(u), abs(v)) == lo:
            flags.add("smallest")
        if u == v:
            flags.add("duplicate")
        if (u > 0) != (v > 0):
            flags.add("opposite")
        if min(abs(u), abs(v)) == lo and max(abs(u), abs(v)) == hi:
            flags.add("largest-smallest")
        if abs(abs(u) - abs(v)) == maxgap:
            flags.add("max-gap")
        rows.append((u, v, D, flags))
    return L, rows

def draw(rng, kind):
    for _ in range(10000):
        if kind == "small":
            ns = [rng.randint(1, 7) for _ in range(rng.randint(8, 15))]
        elif kind == "mixed":
            count = rng.randint(5, 9)
            ns = ([rng.randint(1, 4) for _ in range(count - 2)]
                  + [rng.randint(10, 30), rng.randint(10, 35)])
        elif kind == "medium":
            ns = [rng.randint(5, 35) for _ in range(rng.randint(4, 8))]
        else:
            ns = rng.sample(range(2, 36), rng.randint(5, 8))

        # One sign per absolute label, so B is pair-free.
        signs = {n: (-1 if rng.randrange(2) else 1) for n in set(ns)}
        B = canon([signs[n] * n for n in ns])
        W, mx = sum(ns), max(ns)
        choices = [p for p in range(mx + 1, min(80, W - 16) + 1)
                   if (W - p) % 2 == 0 and (W - p) // 2 >= 8
                   and mx <= (W - p) // 2]
        if 40 <= W <= 200 and choices and sum(n >= 3 for n in ns) >= 2:
            p = rng.choice(choices)
            sigma = -1 if sum(z < 0 for z in B) % 2 else 1
            return B, sigma * p
    raise RuntimeError(kind)

rng = random.Random(168168)
cases = []
for kind in ("small", "mixed", "medium", "distinct"):
    made = []
    while len(made) < 4:
        case = draw(rng, kind)
        if case not in made:
            made.append(case)
    cases.extend((kind, B, p) for B, p in made)

B124 = tuple(map(int, """
25 7 -1 -1 -1 -1 -1 2 -1 -1 -1 -1 -1 -1 -1 2 -1 -1 -1 2
-1 -1 -1 -1 2 -1 -1 -1 2 -1 -1 -1 -1 2 -1 -1 2 2 2 2 -1
2 -1 2 -1 2 -1 2 -1 2 2 -1 -1 2 -1 2 -1 -1 -1 -1 -1 2 2
-1 2 2 2 2 2
""".split()))
B144 = tuple(map(int, """
-14 -3 -3 2 -3 -1 2 2 2 2 2 -3 2 2 -3 -3 2 2 2 2 2 2
-3 -3 2 2 2 -3 2 2 -3 -3 2 2 2 2 -3 -3 -3 -3 2 2 -3
-3 -3 2 -3 2 -3 2 2 2 2 2 2 -3
""".split()))
cases += [
    ("killer124", B124, 54),
    ("killer144", B144, -44),
    ("p-miss1", (-2, 9, 11, 14, -15, -17), -30),
    ("p-miss2", (-1, 2, 6, 8, 9, 12, -14), 16),
    ("opposite-miss", (5, 12, 13, 18, 28, 29, 30, 34), 37),
    ("extreme-miss", (-2, 7, 10, 13, -14, 15, -16, -24), 43),
]

subsets = [s for k in range(1, len(FAMS) + 1)
           for s in combinations(FAMS, k)]
stats, pooled = {}, []
for name, B, p in cases:
    group = name if name in ("small", "mixed", "medium", "distinct") else "fixtures"
    L, rows = analyze(B, p)
    assert rows, (name, "expected a flip descent")
    pooled.append((group, L, p, rows))
    st = stats.setdefault(group, {"n": 0, "miss": {f: 0 for f in FAMS},
                                  "cover": {s: 0 for s in subsets}})
    st["n"] += 1
    for f in FAMS:
        if not any(f in row[3] for row in rows):
            st["miss"][f] += 1
    for s in subsets:
        if any(any(f in row[3] for f in s) for row in rows):
            st["cover"][s] += 1
    if name in ("p-miss1", "p-miss2", "opposite-miss", "extreme-miss"):
        print("WITNESS", name, L, p,
              [(u, v, D) for u, v, D, _ in rows])

def minimum_unions(cover, n):
    valid = [s for s in subsets if cover[s] == n]
    size = min(map(len, valid))
    return [s for s in valid if len(s) == size]

for group, st in stats.items():
    print("GROUP", group, "n", st["n"], "misses", st["miss"],
          "minimum unions", minimum_unions(st["cover"], st["n"]))

pooled_miss = {
    f: sum(not any(f in row[3] for row in rows)
           for _, _, _, rows in pooled)
    for f in FAMS
}
pooled_cover = {
    s: sum(any(any(f in row[3] for f in s) for row in rows)
           for _, _, _, rows in pooled)
    for s in subsets
}
print("POOLED", len(pooled), "misses", pooled_miss,
      "minimum unions", minimum_unions(pooled_cover, len(pooled)))
assert pooled_miss == {
    "distinguished": 2, "smallest": 2, "duplicate": 10,
    "opposite": 1, "largest-smallest": 8, "max-gap": 8
}
