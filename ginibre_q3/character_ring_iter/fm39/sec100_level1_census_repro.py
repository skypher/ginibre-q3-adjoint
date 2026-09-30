"""FM-SEC100 (luna_max_uranus): level-1 census P_T vs N_T (length <= 9, labels <= 4);
component weights; local smoothing obstruction at (1,2,1,2); consumer checks."""
from functools import lru_cache
from itertools import combinations_with_replacement, combinations
from collections import Counter, defaultdict
from math import factorial, comb, prod

@lru_cache(None)
def fusion(xs):
    # SU(2) tensor-product multiplicities, indexed by highest weight.
    dp = {0: 1}
    for n in xs:
        nxt = defaultdict(int)
        for j, multiplicity in dp.items():
            for k in range(abs(j - n), j + n + 1, 2):
                nxt[k] += multiplicity
        dp = dict(nxt)
    return dp.get(0, 0)

def split_counts(labels, i=0, j=1):
    positive = negative = 0
    for mask in range(1 << len(labels)):
        left = tuple(labels[k] for k in range(len(labels))
                     if (mask >> k) & 1)
        right = tuple(labels[k] for k in range(len(labels))
                      if not ((mask >> k) & 1))
        count = fusion(left) * fusion(right)
        same_side = bool((mask >> i) & 1) == bool((mask >> j) & 1)
        if same_side:
            positive += count
        else:
            negative += count
    return positive, negative

rows = []
positive_profiles = zero_profiles = negative_profiles = 0
total_signatures = total_ordered_profiles = 0

for length in range(2, 10):
    signatures = {}
    for labels in combinations_with_replacement(range(1, 5), length):
        if sum(labels) % 2:
            continue
        for i, j in combinations(range(length), 2):
            pair = tuple(sorted((labels[i], labels[j])))
            rest = tuple(labels[k] for k in range(length)
                         if k not in (i, j))
            signature = (pair, rest)
            if signature not in signatures:
                signatures[signature] = split_counts(pair + rest)

    profile_count = positive_total = negative_total = 0
    for (pair, rest), (P, N) in signatures.items():
        position_choices = comb(length, 2) * (
            2 if pair[0] != pair[1] else 1
        )
        rest_orders = factorial(len(rest)) // prod(
            factorial(n) for n in Counter(rest).values()
        )
        multiplicity = position_choices * rest_orders

        profile_count += multiplicity
        positive_total += multiplicity * P
        negative_total += multiplicity * N
        if P > N:
            positive_profiles += multiplicity
        elif P == N:
            zero_profiles += multiplicity
        else:
            negative_profiles += multiplicity

    rows.append((length, len(signatures), profile_count,
                 positive_total, negative_total))
    total_signatures += len(signatures)
    total_ordered_profiles += profile_count

for row in rows:
    print(*row)
print("totals", total_signatures, total_ordered_profiles,
      sum(row[3] for row in rows), sum(row[4] for row in rows))
print("profile signs", positive_profiles, zero_profiles, negative_profiles)
for labels, i, j in [
    ((1, 1, 1, 1), 0, 1),
    ((1, 5, 2, 2), 0, 1),
    ((1, 2, 1, 2), 0, 1),
    ((1, 2, 1, 2), 0, 2),
]:
    print("boundary", labels, (i, j), split_counts(labels, i, j))


# ---- block ----
from itertools import combinations, product
from collections import Counter

def all_matchings(owners):
    def rec(left):
        if not left:
            yield ()
            return
        a = left[0]
        for ix in range(1, len(left)):
            b = left[ix]
            if owners[a] == owners[b]:
                continue
            rest = left[1:ix] + left[ix + 1:]
            for tail in rec(rest):
                yield ((a, b),) + tail
    yield from rec(tuple(range(len(owners))))

def crosses(edge1, edge2):
    a, b = sorted(edge1)
    c, d = sorted(edge2)
    return a < c < b < d or c < a < d < b

def component_weight(labels, T, matching):
    # Return None if the matching is outside the component model.
    n = len(labels)
    eps = [-1 if i in T else 1 for i in range(n)]
    owners = tuple(i for i, size in enumerate(labels) for _ in range(size))
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for a, b in matching:
        union(owners[a], owners[b])

    roots = sorted({find(i) for i in range(n)})
    root_id = {root: k for k, root in enumerate(roots)}
    block_component = [root_id[find(i)] for i in range(n)]
    chord_component = [block_component[owners[a]] for a, b in matching]

    # Each block-connected component must be internally noncrossing.
    for x, y in combinations(range(len(matching)), 2):
        if (chord_component[x] == chord_component[y]
                and crosses(matching[x], matching[y])):
            return None

    # Crossing graph on block-connected components.
    adjacency = [set() for _ in roots]
    for x, y in combinations(range(len(matching)), 2):
        if (chord_component[x] != chord_component[y]
                and crosses(matching[x], matching[y])):
            u, v = chord_component[x], chord_component[y]
            adjacency[u].add(v)
            adjacency[v].add(u)

    color = {}
    factors = []
    for start in range(len(adjacency)):
        if start in color:
            continue
        color[start] = 0
        stack = [start]
        vertices = []
        while stack:
            v = stack.pop()
            vertices.append(v)
            for w in adjacency[v]:
                if w not in color:
                    color[w] = 1 - color[v]
                    stack.append(w)
                elif color[w] == color[v]:
                    return None

        vertices = set(vertices)
        side = [1, 1]
        for i in range(n):
            v = block_component[i]
            if v in vertices:
                side[color[v]] *= eps[i]
        factors.append(side[0] + side[1])

    weight = 1
    for factor in factors:
        weight *= factor
    return weight

def matching_summary(labels, T):
    owners = tuple(i for i, size in enumerate(labels) for _ in range(size))
    states = []
    for matching in all_matchings(owners):
        weight = component_weight(labels, T, matching)
        if weight is not None:
            states.append((matching, weight))
    counts = Counter("positive" if w > 0 else
                     "negative" if w < 0 else "zero"
                     for _, w in states)
    return counts, sum(w for _, w in states), states

for labels, T in [
    ((1, 1, 1, 1), {0, 1}),
    ((1, 5, 2, 2), {0, 1}),
    ((1, 2, 1, 2), {0, 1}),
    ((1, 2, 1, 2), {0, 2}),
]:
    counts, total, states = matching_summary(labels, T)
    print("grouped", labels, tuple(sorted(T)), dict(counts), total)
    if labels == (1, 1, 1, 1):
        print("states", states)

def smooth_crossing_pair(matching, x, y):
    a, b = sorted(matching[x])
    c, d = sorted(matching[y])
    if not (a < c < b < d):
        return []
    other = [matching[k] for k in range(len(matching))
             if k not in (x, y)]
    return [
        tuple(sorted(other + [(a, c), (b, d)])),
        tuple(sorted(other + [(a, d), (c, b)])),
    ]

def test_profile(labels, T):
    owners = tuple(i for i, size in enumerate(labels) for _ in range(size))
    negatives = []
    for matching in all_matchings(owners):
        weight = component_weight(labels, T, matching)
        if weight is not None and weight < 0:
            negatives.append((matching, weight))

    targets = {}
    for matching, weight in negatives:
        candidates = []
        for x, y in combinations(range(len(matching)), 2):
            for target in smooth_crossing_pair(matching, x, y):
                target_weight = component_weight(labels, T, target)
                if target_weight is not None and target_weight > 0:
                    candidates.append(target)
        if not candidates:
            return "undefined", matching, weight
        target = min(candidates)
        if target in targets:
            return "collision", targets[target], (matching, weight), target
        targets[target] = (matching, weight)
    return "ok", len(negatives), len(targets)

for endpoint_count in range(4, 7):
    found = False
    for length in range(2, endpoint_count + 1):
        for labels in product(range(1, 5), repeat=length):
            if sum(labels) != endpoint_count:
                continue
            for T_tuple in combinations(range(length), 2):
                result = test_profile(labels, set(T_tuple))
                if result[0] != "ok":
                    print("first issue", endpoint_count, length, labels,
                          T_tuple, result)
                    if result[0] == "undefined":
                        matching = result[1]
                        owners = tuple(i for i, size in enumerate(labels)
                                       for _ in range(size))
                        for x, y in combinations(range(len(matching)), 2):
                            targets = smooth_crossing_pair(matching, x, y)
                            if targets:
                                print("resolutions", matching[x], matching[y],
                                      targets,
                                      [component_weight(labels, set(T_tuple), z)
                                       for z in targets])
                    found = True
                    break
            if found:
                break
        if found:
            break
        print("passed endpoint count", endpoint_count)
    if found:
        break


# ---- block ----
from math import comb
from functools import lru_cache
from sympy import symbols, Poly, expand, Rational

x, y = symbols("x y")

@lru_cache(None)
def semicircle_moment(n):
    if n % 2:
        return Rational(0)
    m = n // 2
    return Rational(comb(2*m, m), m + 1)

def average(poly):
    return sum(
        coefficient * semicircle_moment(i) * semicircle_moment(j)
        for (i, j), coefficient in Poly(expand(poly), x, y).terms()
    )

weyl_denominator = average((x - y)**2)
hat_s2 = x**2 + y**2 - 2

for k in range(1, 7):
    value = average((x - y)**2 * hat_s2**k * (x + y)**2)
    print("hatS2 power", k, value / weyl_denominator)

U = [1, x]
for n in range(2, 5):
    U.append(expand(x * U[-1] - U[-2]))
# main-agent fix: the second factor must be in y (the printed version used x twice and gave 25)
Uy = [expand(u).subs(x, y) if hasattr(expand(u), "subs") else u for u in U]
h4 = expand(sum(U[k] * Uy[4-k] for k in range(5)))
value = average((x-y)**2 * h4 * hat_s2**2)
print("h4 hatS2^2", value / weyl_denominator)
