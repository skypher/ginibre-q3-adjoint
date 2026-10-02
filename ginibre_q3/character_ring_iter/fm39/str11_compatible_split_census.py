from itertools import product
from datetime import datetime
from math import prod
import random

def stamp(s):
    print(datetime.now().strftime("%H:%M:%S"), s, flush=True)

def step(T, n, eps):
    out = {}
    for (i, j), v in T.items():
        for ii in range(abs(i-n), i+n+1, 2):
            out[(ii, j)] = out.get((ii, j), 0) + v
        for jj in range(abs(j-n), j+n+1, 2):
            out[(i, jj)] = out.get((i, jj), 0) + eps*v
    return {k: v for k, v in out.items() if v}

def build_tables(counts, types):
    zero = (0,) * len(counts)
    memo = {zero: {(0, 0): 1}}

    def get(v):
        if v not in memo:
            t = next(i for i, x in enumerate(v) if x)
            w = list(v)
            w[t] -= 1
            memo[v] = step(get(tuple(w)), *types[t])
        return memo[v]

    for v in product(*(range(c+1) for c in counts)):
        get(v)
    return memo

def compatible(A, B):
    return all(A[k] * B[k] >= 0 for k in A.keys() & B.keys())

def weight(v, types):
    return sum(c*n for c, (n, e) in zip(v, types))

def factor_count(v):
    return sum(v)

stamp("bounded census begin: labels 1..5, lengths 2..10")
profiles = profiles_with_split = profiles_no_split = 0
candidate_oriented = compatible_ordered = compatible_unordered = 0
max_peel_fail = min_peel_fail = 0
shape = {"peel": 0, "balanced_factors": 0,
         "balanced_weight": 0, "near_half_counts": 0}
per_length = {}

for counts0 in product(range(11), repeat=5):
    L = sum(counts0)
    if not 2 <= L <= 10:
        continue
    support = [i for i, c in enumerate(counts0) if c]
    for signs in product((1, -1), repeat=len(support)):
        minus = sum(counts0[i] for i, e in zip(support, signs) if e == -1)
        if minus % 2:
            continue
        types = tuple((i+1, e) for i, e in zip(support, signs))
        counts = tuple(counts0[i] for i in support)
        profiles += 1
        per_length[L] = per_length.get(L, 0) + 1
        tabs = build_tables(counts, types)
        local = []

        for a in product(*(range(c+1) for c in counts)):
            if not any(a) or a == counts:
                continue
            b = tuple(c-x for c, x in zip(counts, a))
            candidate_oriented += 1
            if compatible(tabs[a], tabs[b]):
                compatible_ordered += 1
                local.append((a, b))
                if a <= b:
                    compatible_unordered += 1
                    na, nb = factor_count(a), factor_count(b)
                    wa, wb = weight(a, types), weight(b, types)
                    if min(na, nb) == 1:
                        shape["peel"] += 1
                    if abs(na-nb) <= 1:
                        shape["balanced_factors"] += 1
                    if abs(wa-wb) <= max(n for n, e in types):
                        shape["balanced_weight"] += 1
                    if all(abs(x-y) <= 1 for x, y in zip(a, b)):
                        shape["near_half_counts"] += 1

        if local:
            profiles_with_split += 1
        else:
            profiles_no_split += 1

        for index, kind in ((len(types)-1, "max"), (0, "min")):
            a = [0] * len(counts)
            a[index] = 1
            a = tuple(a)
            b = tuple(c-x for c, x in zip(counts, a))
            if not compatible(tabs[a], tabs[b]):
                if kind == "max":
                    max_peel_fail += 1
                else:
                    min_peel_fail += 1

        if profiles % 3000 == 0:
            stamp(f"bounded profiles={profiles}; with_split={profiles_with_split}")

stamp(f"bounded done profiles={profiles}, with_split={profiles_with_split}, "
      f"no_split={profiles_no_split}")
print("by_length", per_length)
print("oriented cuts", candidate_oriented,
      "compatible ordered", compatible_ordered,
      "compatible unordered", compatible_unordered)
print("compatible shape counts", shape)
print("max/min peel failures", max_peel_fail, min_peel_fail)

assert profiles == profiles_with_split == 19018
assert profiles_no_split == max_peel_fail == min_peel_fail == 0
assert candidate_oriented == 1524902
assert compatible_ordered == 1143612
assert compatible_unordered == 572647
assert shape == {"peel": 72490, "balanced_factors": 129685,
                 "balanced_weight": 163754, "near_half_counts": 51575}

stamp("random census begin: 300 seeded profiles")
rng = random.Random(20261002)
nsamp = 300
random_no_split = random_max_fail = random_min_fail = 0
weights = []

for z in range(nsamp):
    L = rng.randint(3, 10)
    raw = [rng.randint(1, 18) for _ in range(L)]
    counts_by_label = {}
    for n in raw:
        counts_by_label[n] = counts_by_label.get(n, 0) + 1
    labels = sorted(counts_by_label)
    eps = [rng.choice((1, -1)) for _ in labels]
    if sum(counts_by_label[n] for n, e in zip(labels, eps) if e == -1) % 2:
        idx = next(i for i, n in enumerate(labels) if counts_by_label[n] % 2)
        eps[idx] *= -1

    types = tuple(zip(labels, eps))
    counts = tuple(counts_by_label[n] for n in labels)
    W = sum(n*c for n, c in zip(labels, counts))
    weights.append(W)
    tabs = build_tables(counts, types)
    found = False

    for a in product(*(range(c+1) for c in counts)):
        if not any(a) or a == counts:
            continue
        b = tuple(c-x for c, x in zip(counts, a))
        if compatible(tabs[a], tabs[b]):
            found = True
            break
    if not found:
        random_no_split += 1

    for idx, kind in ((len(types)-1, "max"), (0, "min")):
        a = [0] * len(counts)
        a[idx] = 1
        a = tuple(a)
        b = tuple(c-x for c, x in zip(counts, a))
        if not compatible(tabs[a], tabs[b]):
            if kind == "max":
                random_max_fail += 1
            else:
                random_min_fail += 1

    if z % 50 == 0:
        stamp(f"random {z}/300; no_split={random_no_split}; "
              f"max_peel_fail={random_max_fail}")

stamp(f"random done samples={nsamp}, W_range=[{min(weights)},{max(weights)}], "
      f"meanW={sum(weights)/nsamp:.3f}, no_split={random_no_split}, "
      f"max/min peel failures={random_max_fail}/{random_min_fail}")

assert (min(weights), max(weights), sum(weights),
        random_no_split, random_max_fail, random_min_fail) == (10, 125, 18617, 0, 0, 0)
print("ALL ASSERTIONS PASS")
