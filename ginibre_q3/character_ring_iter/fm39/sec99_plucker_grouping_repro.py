"""FM-SEC99 (luna_max_venus): component-to-coset identity f = sum_M 1_(A_M+W_M); recursive Pluecker family test
(failure at (1^6)); 15-subgroup H_AC certificate at (1^6); (1^8,6) check."""
from itertools import combinations
from collections import defaultdict

def matchings(legs, owner):
    if not legs:
        yield ()
        return
    a = legs[0]
    for j in range(1, len(legs)):
        b = legs[j]
        if owner[a] != owner[b]:
            for tail in matchings(legs[1:j] + legs[j+1:], owner):
                yield tuple(sorted(((a, b),) + tail))

def cross(e, f):
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b

def n_cross(M):
    return sum(cross(e, f) for e, f in combinations(M, 2))

def span(gens):
    W = {0}
    for g in gens:
        W |= {x ^ g for x in list(W)}
    return W

def coset_term(M, L, owner, G):
    par = list(range(L))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a, b in M:
        par[find(owner[a])] = find(owner[b])

    comp = defaultdict(list)
    for e in M:
        comp[find(owner[e[0]])].append(e)
    nodes = list(comp)

    if any(cross(e, f) for v in nodes
           for e, f in combinations(comp[v], 2)):
        return [0] * G

    adj = {v: set() for v in nodes}
    for v, w in combinations(nodes, 2):
        if any(cross(e, f) for e in comp[v] for f in comp[w]):
            adj[v].add(w)
            adj[w].add(v)

    blockmask = {
        v: sum(1 << i for i in range(L) if find(i) == v)
        for v in nodes
    }
    color = {}
    gens = []
    A = 0

    for root in nodes:
        if root in color:
            continue
        color[root] = 0
        todo = [root]
        cluster = []
        while todo:
            v = todo.pop()
            cluster.append(v)
            for w in adj[v]:
                if w not in color:
                    color[w] = 1 - color[v]
                    todo.append(w)
                elif color[w] == color[v]:
                    return [0] * G

        g = 0
        for v in cluster:
            g ^= blockmask[v]
            if color[v]:
                A ^= blockmask[v]
        gens.append(g)

    W = span(gens)
    return [int((s ^ A) in W) for s in range(G)]

def fusion(labels):
    d = {0: 1}
    for n in labels:
        e = {}
        for x, v in d.items():
            for y in range(abs(x-n), x+n+1, 2):
                e[y] = e.get(y, 0) + v
        d = e
    return d.get(0, 0)

def atlas(labels):
    owner = []
    for i, n in enumerate(labels):
        owner += [i] * n
    L = len(labels)
    G = 1 << (L - 1)
    Ms = list(matchings(list(range(len(owner))), owner))

    total = [0] * G
    for M in Ms:
        v = coset_term(M, L, owner, G)
        total = [x+y for x, y in zip(total, v)]

    direct = []
    for S in range(G):
        a = fusion([labels[i] for i in range(L) if (S >> i) & 1])
        b = fusion([labels[i] for i in range(L)
                    if not ((S >> i) & 1)])
        direct.append(a*b)

    assert total == direct
    return owner, Ms, total

def smooth_children(M, owner):
    crossings = sorted((e, f) for e, f in combinations(M, 2)
                       if cross(e, f))
    if not crossings:
        return []
    (a, b), (c, d) = crossings[0]
    if c < a:
        a, b, c, d = c, d, a, b

    out = []
    for raw in [((a, c), (b, d)), ((a, d), (b, c))]:
        raw = tuple(tuple(sorted(e)) for e in raw)
        child = tuple(sorted(
            tuple(e for e in M if e not in ((a, b), (c, d))) + raw
        ))
        if all(owner[x] != owner[y] for x, y in child):
            out.append(child)
    return out

def descendants(M, owner, memo):
    if M not in memo:
        out = [M]
        for child in smooth_children(M, owner):
            assert n_cross(child) < n_cross(M)
            out += descendants(child, owner, memo)
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

boundary = [
    (1, 1, 1, 1),
    (1,)*6,
    (1,)*8 + (6,),
    (1, 1, 1, 1, 2, 2, 2),
    (1, 2, 1, 2)
]
for labels in boundary:
    _, Ms, v = atlas(labels)
    print("atlas PASS", labels, "matchings", len(Ms),
          "G-size", len(v), "f(0)", v[0])

for labels in [(1, 1, 1, 1), (1,)*6, (1, 2, 1, 2),
               (1, 1, 1, 1, 2, 2, 2)]:
    owner, Ms, _ = atlas(labels)
    L = len(labels)
    G = 1 << (L-1)
    memo = {}
    seen = set()
    negative = []
    count = 0

    for M in sorted(Ms):
        if not n_cross(M):
            continue
        family = tuple(descendants(M, owner, memo))
        if family in seen:
            continue
        seen.add(family)

        value = [0] * G
        for D in family:
            term = coset_term(D, L, owner, G)
            value = [x+y for x, y in zip(value, term)]
        if not any(value):
            continue

        count += 1
        transform = wht(value)
        if min(transform) < 0:
            negative.append((M, family, value, transform))

    print("rooted families", labels, "count", count,
          "negative Fourier", len(negative))
    if labels == (1,)*6:
        M, family, value, transform = negative[0]
        assert M == ((0, 2), (1, 4), (3, 5))
        assert transform[2] == -1
        print("witness root", M, "family size", len(family),
              "Fourier[2]", transform[2])
        print("family", sorted(family))
        print("family values", value)


# ---- block ----
from fractions import Fraction as F

terms = [
 (F(9,10),(0,15)),
 (F(3,20),(0,3,9,10)),
 (F(3,10),(0,3,24,27)),
 (F(3,10),(0,3,29,30)),
 (F(3,10),(0,5,18,23)),
 (F(3,20),(0,6,10,12)),
 (F(3,10),(0,6,17,23)),
 (F(3,10),(0,9,20,29)),
 (F(3,10),(0,10,17,27)),
 (F(3,40),(0,3,5,6,17,18,20,23)),
 (F(3,40),(0,3,5,6,24,27,29,30)),
 (F(9,40),(0,5,9,12,17,20,24,29)),
 (F(9,40),(0,5,9,12,18,23,27,30)),
 (F(3,10),(0,6,10,12,18,20,24,30)),
 (F(11,10),(0,3,5,6,9,10,12,15,
             17,18,20,23,24,27,29,30))
]

def is_subspace(W):
    S = set(W)
    return 0 in S and all((x ^ y) in S for x in S for y in S)

cat = [1, 1, 2, 5]
for c, W in terms:
    assert c > 0 and is_subspace(W)

for s in range(32):
    k = s.bit_count()
    target = cat[k//2] * cat[(6-k)//2] if k % 2 == 0 else 0
    got = sum(c for c, W in terms if s in W)
    assert got == target

print("15-subgroup regrouping: exact at all 32 quotient points")


# ---- block ----
def fusion(labels):
    d = {0: 1}
    for n in labels:
        e = {}
        for x, v in d.items():
            for y in range(abs(x-n), x+n+1, 2):
                e[y] = e.get(y, 0) + v
        d = e
    return d.get(0, 0)

labels = [1]*8 + [6]
for S in range(256):
    left = fusion([labels[i] for i in range(8) if (S >> i) & 1])
    right = fusion([labels[i] for i in range(8)
                    if not ((S >> i) & 1)] + [6])
    target = left * right
    expected = (7 if S == 0 else 0) + (1 if S.bit_count() == 2 else 0)
    assert target == expected

print("(1^8,6): f = 7 delta_0 + sum_{|S|=2} delta_S")
print("p=weight-one indicator: p*p(0)=8, p*p(weight 2)=2")
print("therefore f=p*p/2+3 delta_0")
