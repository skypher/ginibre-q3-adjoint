"""FM-SEC94 (luna_max_jupiter): component-model charging; single-smoothing failure at (1,2,1,2); recursive Pluecker repair on boundary cases."""
from functools import lru_cache
from itertools import combinations, product

def cross(e, f):
    a, b = sorted(e)
    c, d = sorted(f)
    return a < c < b < d or c < a < d < b

def smooth(M, e, f):
    a, b, c, d = sorted(tuple(e) + tuple(f))
    for p, q in (((a, b), (c, d)), ((a, d), (b, c))):
        yield tuple(sorted(
            (set(M) - {tuple(sorted(e)), tuple(sorted(f))})
            | {tuple(sorted(p)), tuple(sorted(q))}
        ))

def owners(labels):
    return tuple(i for i, n in enumerate(labels) for _ in range(n))

@lru_cache(None)
def matchings(labels):
    own = owners(labels)

    @lru_cache(None)
    def rec(rem):
        if not rem:
            return ((),)
        a = rem[0]
        out = []
        for j in range(1, len(rem)):
            b = rem[j]
            if own[a] == own[b]:
                continue
            for tail in rec(rem[1:j] + rem[j + 1:]):
                out.append(tuple(sorted(((a, b),) + tail)))
        return tuple(out)

    return rec(tuple(range(len(own))))

def raw_components(labels, M):
    own = owners(labels)
    parent = list(range(len(labels)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[x] = y

    for a, b in M:
        union(own[a], own[b])

    out = {}
    for e in M:
        out.setdefault(find(own[e[0]]), []).append(e)
    return out

def component_data(labels, M):
    own = owners(labels)
    edges = raw_components(labels, M)

    if any(own[a] == own[b] for a, b in M):
        return None
    if any(cross(e, f) for es in edges.values()
           for e, f in combinations(es, 2)):
        return None

    blocks = {
        r: {own[p] for e in es for p in e}
        for r, es in edges.items()
    }
    roots = list(edges)
    adj = {r: set() for r in roots}

    for i, r in enumerate(roots):
        for s in roots[i + 1:]:
            if any(cross(e, f) for e in edges[r] for f in edges[s]):
                adj[r].add(s)
                adj[s].add(r)

    return edges, blocks, adj

def config_weight(labels, T, M):
    data = component_data(labels, M)
    if data is None:
        return None

    edges, blocks, adj = data
    color = {}
    negative_factors = []
    weight = 1

    for r in edges:
        if r in color:
            continue

        color[r] = 0
        todo = [r]
        cc = []

        while todo:
            u = todo.pop()
            cc.append(u)
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    todo.append(v)
                elif color[v] == color[u]:
                    return None

        side = [set(), set()]
        for u in cc:
            side[color[u]] |= blocks[u]

        eps = lambda B: (-1) ** sum((T >> i) & 1 for i in B)
        factor = eps(side[0]) + eps(side[1])
        weight *= factor

        if factor < 0:
            negative_factors.append((
                {u for u in cc if color[u] == 0},
                {u for u in cc if color[u] == 1},
            ))

    return weight, negative_factors, edges

def negative_crossings(data):
    _, negative_factors, edges = data
    out = []
    for side0, side1 in negative_factors:
        for u in side0:
            for v in side1:
                for e in edges[u]:
                    for f in edges[v]:
                        if cross(e, f):
                            out.append((tuple(sorted(e)), tuple(sorted(f))))
    return sorted(set(out))

def one_step_targets(labels, T, M, data, cache):
    out = set()
    for e, f in negative_crossings(data):
        for N in smooth(M, e, f):
            if N not in cache:
                cache[N] = config_weight(labels, T, N)
            z = cache[N]
            if z is not None and z[0] > 0:
                out.add(N)
    return out

@lru_cache(None)
def canonical_descendants(M):
    out = {M}
    pairs = sorted(
        (tuple(sorted(e)), tuple(sorted(f)))
        for e, f in combinations(M, 2) if cross(e, f)
    )
    if pairs:
        for N in smooth(M, *pairs[0]):
            out.update(canonical_descendants(N))
    return frozenset(out)

def seeded_reachable(M, data):
    out = set()
    pairs = negative_crossings(data)
    if pairs:
        e, f = pairs[0]
        for N in smooth(M, e, f):
            out.update(canonical_descendants(N))
    return out

def maxflow(demand, supply, edges):
    n, m = len(demand), len(supply)
    N = n + m + 2
    source, sink = N - 2, N - 1
    cap = [[0] * N for _ in range(N)]
    adj = [[] for _ in range(N)]

    def add(u, v, c):
        if v not in adj[u]:
            adj[u].append(v)
            adj[v].append(u)
        cap[u][v] += c

    for i, x in enumerate(demand):
        add(source, i, x)
    for j, x in enumerate(supply):
        add(n + j, sink, x)
    for i, js in enumerate(edges):
        for j in js:
            add(i, n + j, sum(demand))

    flow = 0
    while True:
        parent = [-1] * N
        parent[source] = source
        queue = [source]

        for u in queue:
            for v in adj[u]:
                if parent[v] < 0 and cap[u][v]:
                    parent[v] = u
                    queue.append(v)

        if parent[sink] < 0:
            return flow

        amount = 10**9
        v = sink
        while v != source:
            amount = min(amount, cap[parent[v]][v])
            v = parent[v]

        v = sink
        while v != source:
            u = parent[v]
            cap[u][v] -= amount
            cap[v][u] += amount
            v = u
        flow += amount

def fusion(labels):
    state = {0: 1}
    for n in labels:
        nxt = {}
        for j, c in state.items():
            for k in range(abs(j - n), j + n + 1, 2):
                nxt[k] = nxt.get(k, 0) + c
        state = nxt
    return state.get(0, 0)

def direct(labels, T):
    L = len(labels)
    full = (1 << L) - 1
    total = 0
    for S in range(1 << L):
        a = fusion(tuple(labels[i] for i in range(L) if (S >> i) & 1))
        b = fusion(tuple(labels[i] for i in range(L)
                         if ((full ^ S) >> i) & 1))
        total += (-1) ** ((S & T).bit_count()) * a * b
    return total

def audit(labels, Tset):
    T = sum(1 << i for i in Tset)
    valid = {}
    Ms = matchings(labels)

    for M in Ms:
        z = config_weight(labels, T, M)
        if z is not None:
            valid[M] = z

    assert sum(z[0] for z in valid.values()) == direct(labels, T)

    neg = [M for M, z in valid.items() if z[0] < 0]
    pos = [M for M, z in valid.items() if z[0] > 0]
    pos_index = {M: j for j, M in enumerate(pos)}
    cache = dict(valid)
    one_edges = []
    recursive_edges = []

    for M in neg:
        one_edges.append({
            pos_index[P]
            for P in one_step_targets(labels, T, M, valid[M], cache)
            if P in pos_index
        })
        recursive_edges.append({
            pos_index[P]
            for P in seeded_reachable(M, valid[M])
            if P in pos_index
        })

    demand = [-valid[M][0] for M in neg]
    supply = [valid[M][0] for M in pos]
    print(
        "CASE", labels, "T", Tset,
        "matchings", len(Ms),
        "valid", len(valid),
        "negative", len(neg),
        "positive", len(pos),
        "demand", sum(demand),
        "capacity", sum(supply),
        "one_step_isolated", sum(not e for e in one_edges),
        "recursive_flow", maxflow(demand, supply, recursive_edges),
        "F", direct(labels, T),
    )

    if labels == (1, 2, 1, 2):
        M = neg[0]
        print("  source", M, "weight", valid[M][0])
        for e, f in negative_crossings(valid[M]):
            for N in smooth(M, e, f):
                groups = raw_components(labels, N)
                internal = [
                    (x, y)
                    for es in groups.values()
                    for x, y in combinations(es, 2)
                    if cross(x, y)
                ]
                print("  smoothing", e, f, "->", N,
                      "internal_crossings", internal,
                      "model_weight", config_weight(labels, T, N))
        print("  positive recursive targets", [
            (P, valid[P][0])
            for P in pos if P in seeded_reachable(M, valid[M])
        ])

for labels, T in [
    ((1, 1, 1, 1), [0, 1]),
    ((1, 1, 2, 2), [0, 2]),
    ((2, 2, 2, 2, 2, 2), [4, 5]),
    ((1, 2, 1, 2), [0, 1]),
]:
    audit(labels, T)

first = None
for labels in product(range(1, 5), repeat=4):
    if sum(labels) % 2:
        continue
    for T in range(1, 16):
        if T.bit_count() % 2:
            continue
        cache = {}
        for M in matchings(labels):
            z = config_weight(labels, T, M)
            if z is not None and z[0] < 0:
                if not one_step_targets(labels, T, M, z, cache):
                    first = (labels, T, M, z[0])
                    break
        if first:
            break
    if first:
        break
print("FIRST_ORDERED_LENGTH4_ONE_STEP_FAILURE", first)

labels = (1, 2, 1, 2)
T = 0b0011
M = ((0, 3), (1, 5), (2, 4))
M1 = ((0, 5), (1, 3), (2, 4))
M2 = ((0, 5), (1, 4), (2, 3))
print("two_step_path",
      config_weight(labels, T, M)[0],
      config_weight(labels, T, M1),
      config_weight(labels, T, M2)[0])
