"""FM-SEC105 (luna_max_venus): groupings G1 (straightening targets), G2 (partitions) fail at (1^4); G3 (NC partitions) factorization and screen."""
from itertools import combinations, product
from collections import defaultdict, Counter
from fractions import Fraction as Q

def matchings(owner):
    def rec(legs):
        if not legs:
            yield ()
            return
        a = legs[0]
        for j in range(1, len(legs)):
            b = legs[j]
            if owner[a] != owner[b]:
                for tail in rec(legs[1:j] + legs[j+1:]):
                    yield tuple(sorted(((a, b),) + tail))
    yield from rec(list(range(len(owner))))

def cross(e, f):
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b

def span(gens):
    W = {0}
    for g in gens:
        W |= {x ^ g for x in list(W)}
    return W

def coset(M, L, owner):
    G = 1 << (L - 1)
    par = list(range(L))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a, b in M:
        par[find(owner[a])] = find(owner[b])

    comps = defaultdict(list)
    for e in M:
        comps[find(owner[e[0]])].append(e)
    nodes = list(comps)

    if any(cross(e, f) for v in nodes
           for e, f in combinations(comps[v], 2)):
        return [0] * G, None

    adj = {v: set() for v in nodes}
    for v, w in combinations(nodes, 2):
        if any(cross(e, f) for e in comps[v] for f in comps[w]):
            adj[v].add(w)
            adj[w].add(v)

    masks = {
        v: sum(1 << i for i in range(L) if find(i) == v)
        for v in nodes
    }
    color, gens, A, part = {}, [], 0, []

    for root in nodes:
        if root in color:
            continue
        color[root] = 0
        todo = [root]
        blocks, g = [], 0
        while todo:
            v = todo.pop()
            blocks += [i for i in range(L) if masks[v] >> i & 1]
            g ^= masks[v]
            for w in adj[v]:
                if w not in color:
                    color[w] = 1 - color[v]
                    todo.append(w)
                elif color[w] == color[v]:
                    return [0] * G, None
            if color[v]:
                A ^= masks[v]
        gens.append(g)
        part.append(tuple(sorted(blocks)))

    W = span(gens)
    return [int((s ^ A) in W) for s in range(G)], tuple(sorted(part))

def children(M, owner):
    crossing = next(((e, f) for e, f in combinations(M, 2)
                     if cross(e, f)), None)
    if crossing is None:
        return []
    (a, b), (c, d) = crossing
    if c < a:
        a, b, c, d = c, d, a, b
    out = []
    for raw in (((a, c), (b, d)), ((a, d), (b, c))):
        child = tuple(sorted(
            tuple(sorted(e))
            for e in tuple(e for e in M if e not in ((a, b), (c, d)))
            + raw
        ))
        if all(owner[i] != owner[j] for i, j in child):
            out.append(child)
    return out

def leaves(M, owner, memo):
    if M not in memo:
        cs = children(M, owner)
        if not cs:
            assert not any(cross(e, f) for e, f in combinations(M, 2))
            memo[M] = Counter({M: 1})
        else:
            out = Counter()
            for C in cs:
                out.update(leaves(C, owner, memo))
            memo[M] = out
    return memo[M]

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

def groups(labels, kind, weight):
    L = len(labels)
    owner = sum(([i] * k for i, k in enumerate(labels)), [])
    G = 1 << (L - 1)
    memo = {}
    out = defaultdict(lambda: [Q(0)] * G)

    for M in matchings(owner):
        v, _ = coset(M, L, owner)
        if not any(v):       # only admissible component-model terms
            continue
        d = leaves(M, owner, memo)
        tau = sum(d.values())
        for N, c in d.items():
            key = N if kind == "G1" else coset(N, L, owner)[1]
            wt = Q(c) if weight == "raw" else Q(c, tau)
            out[key] = [x + wt*y for x, y in zip(out[key], v)]
    return out

for L in (2, 3):
    for labels in product((1, 2, 3), repeat=L):
        if sum(labels) % 2:
            continue
        for kind in ("G1", "G2"):
            for weight in ("raw", "normalized"):
                assert all(min(wht(v)) >= 0
                           for v in groups(labels, kind, weight).values())
print("G1/G2: all lists of lengths 2,3, labels 1..3 pass; four variants")

labels = (1, 1, 1, 1)
for kind in ("G1", "G2"):
    for weight in ("raw", "normalized"):
        out = groups(labels, kind, weight)
        key = ((0, 1), (2, 3))
        z = wht(out[key])
        print(kind, weight, "key", key, "group", out[key],
              "WHT", z, "WHT[1]", z[1])


# ---- block ----
from itertools import combinations_with_replacement

def fusion(labels):
    d = {0: 1}
    for n in labels:
        e = {}
        for x, v in d.items():
            for y in range(abs(x-n), x+n+1, 2):
                e[y] = e.get(y, 0) + v
        d = e
    return d.get(0, 0)

def f_profile(labels):
    L = len(labels)
    return [
        fusion([labels[i] for i in range(L) if S >> i & 1]) *
        fusion([labels[i] for i in range(L) if not (S >> i & 1)])
        for S in range(1 << L)
    ]

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

profiles = entries = 0
minimum = None
where = None
for L in range(1, 8):
    for labels in combinations_with_replacement((1, 2, 3), L):
        z = wht(f_profile(labels))
        profiles += 1
        entries += len(z)
        for chi, value in enumerate(z):
            if minimum is None or value < minimum:
                minimum = value
                where = (labels, chi, value)
            assert value >= 0
print("G3 exhaustive local factors:", profiles, "profiles;",
      entries, "WHT entries; minimum", minimum, "at", where)

for labels in ((1, 1, 1, 1, 2, 2, 2), (1, 2, 1, 2)):
    z = wht(f_profile(labels))
    print("boundary profile", labels, "minimum", min(z),
          "negative entries", sum(x < 0 for x in z))

# Local subprofiles of (1^8,6): (1^k) and (1^k,6), 0 <= k <= 8.
extra = [(1,)*8] + [(1,)*k + (6,) for k in range(9)]
for labels in extra:
    z = wht(f_profile(labels))
    assert min(z) >= 0
    print("special local profile", labels, "min", min(z),
          "max", max(z), "negative entries", sum(x < 0 for x in z))
