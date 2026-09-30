"""FM-SEC103 (luna_max_jupiter): terminal straightening charge; Hall/max-flow census; first failure (1,1,1,1,2), T={0,2}."""
from functools import lru_cache
from itertools import combinations, product
from collections import defaultdict, deque

def cross(e, f):
    a, b = sorted(e); c, d = sorted(f)
    return a < c < b < d or c < a < d < b

def smooth(M, e, f):
    a, b, c, d = sorted(tuple(e) + tuple(f))
    base = set(M) - {tuple(sorted(e)), tuple(sorted(f))}
    return (tuple(sorted(base | {(a,b),(c,d)})),
            tuple(sorted(base | {(a,d),(b,c)})))  # main-agent fix: closing parenthesis

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
            for tail in rec(rem[1:j] + rem[j+1:]):
                out.append(tuple(sorted(((a,b),) + tail)))
        return tuple(out)

    return rec(tuple(range(sum(labels))))

@lru_cache(None)
def straighten(M):
    ef = next(((e,f) for e,f in combinations(M,2) if cross(e,f)), None)
    if ef is None:
        return ((M,1),)
    out = defaultdict(int)
    for child in smooth(M, *ef):
        for N, c in straighten(child):
            out[N] += c
    return tuple(sorted((N,c) for N,c in out.items() if c))

@lru_cache(None)
def factors(labels, M):
    own = owners(labels)
    parent = list(range(len(labels)))

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    for a,b in M:
        if own[a] == own[b]:
            return None
        x, y = find(own[a]), find(own[b])
        if x != y:
            parent[x] = y

    groups = defaultdict(list)
    for e in M:
        groups[find(own[e[0]])].append(e)
    E = list(groups.values())

    if any(cross(e,f) for G in E for e,f in combinations(G,2)):
        return None

    masks = []
    for G in E:
        mask = 0
        for e in G:
            for v in e:
                mask |= 1 << own[v]
        masks.append(mask)

    adj = [set() for _ in E]
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            if any(cross(e,f) for e in E[i] for f in E[j]):
                adj[i].add(j)
                adj[j].add(i)

    color, out = {}, []
    for s in range(len(E)):
        if s in color:
            continue
        color[s] = 0
        todo = [s]
        side = [0,0]
        while todo:
            u = todo.pop()
            side[color[u]] |= masks[u]
            for v in adj[u]:
                if v not in color:
                    color[v] = 1-color[u]
                    todo.append(v)
                elif color[v] == color[u]:
                    return None
        out.append(tuple(side))
    return tuple(out)

def weight(F, T):
    if F is None:
        return None
    ans = 1
    for x,y in F:
        ex = -1 if (x&T).bit_count() % 2 else 1
        ey = -1 if (y&T).bit_count() % 2 else 1
        ans *= ex + ey
    return ans

def maxflow(demand, capacity, edges):
    n, m = len(demand), len(capacity)
    s, t = n+m, n+m+1
    g = [[] for _ in range(t+1)]

    def add(u, v, c):
        g[u].append([v,c,len(g[v])])
        g[v].append([u,0,len(g[u])-1])

    inf = sum(demand)
    for i,x in enumerate(demand):
        add(s, i, x)
    for j,x in enumerate(capacity):
        add(n+j, t, x)
    for i, js in enumerate(edges):
        for j in js:
            add(i, n+j, inf)

    ans = 0
    while True:
        lev = [-1] * len(g)
        lev[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v,c,_ in g[u]:
                if c and lev[v] < 0:
                    lev[v] = lev[u] + 1
                    q.append(v)
        if lev[t] < 0:
            return ans

        it = [0] * len(g)
        def dfs(u, f):
            if u == t:
                return f
            while it[u] < len(g[u]):
                e = g[u][it[u]]
                v, c, r = e
                if c and lev[v] == lev[u] + 1:
                    z = dfs(v, min(f,c))
                    if z:
                        e[1] -= z
                        g[v][r][1] += z
                        return z
                it[u] += 1
            return 0

        while True:
            z = dfs(s, inf)
            if not z:
                break
            ans += z

def analyze(labels, T, Ms=None, FF=None):
    Ms = matchings(labels) if Ms is None else Ms
    FF = {M:factors(labels,M) for M in Ms} if FF is None else FF
    neg, pos, dn, cp = [], [], [], []
    allw = {}

    for M in Ms:
        z = weight(FF[M], T)
        allw[M] = z
        if z is not None and z < 0:
            neg.append(M)
            dn.append(-z)
        elif z is not None and z > 0:
            pos.append(M)
            cp.append(z)

    ix = {M:j for j,M in enumerate(pos)}
    edges = []
    for M in neg:
        edges.append({
            ix[N] for N,c in straighten(M)
            if c > 0 and N in ix
        })

    fl = maxflow(dn, cp, edges) if neg else 0
    return {
        'valid': sum(x is not None for x in allw.values()),
        'neg': len(neg), 'pos': len(pos),
        'demand': sum(dn), 'capacity': sum(cp), 'flow': fl,
        'F': sum(x or 0 for x in allw.values()),
        'negM': neg, 'posM': pos, 'edges': edges, 'weights': allw
    }

profiles = []
for L in range(2,5):
    for ns in product(range(1,5), repeat=L):
        S = sum(ns)
        if S % 2 == 0 and max(ns) <= S//2:
            profiles.append(ns)
for ns in product(range(1,5), repeat=5):
    S = sum(ns)
    if S <= 10 and S % 2 == 0 and max(ns) <= S//2:
        profiles.append(ns)
profiles.sort(key=lambda ns:(len(ns),sum(ns),ns))

failures = []
for ns in profiles:
    Ms = matchings(ns)
    FF = {M:factors(ns,M) for M in Ms}
    for T in range(1 << len(ns)):
        if T.bit_count() % 2:
            continue
        z = analyze(ns, T, Ms, FF)
        if z['flow'] < z['demand']:
            failures.append(
                (ns,T,z['flow'],z['demand'],z['neg'],z['pos'],
                 z['capacity']-z['demand'])
            )

print('CENSUS profiles', len(profiles), 'failures')
for row in failures:
    print(' ', row)

boundaries = [
    ((1,1,1,1),3), ((1,1,2,2),5), ((2,2,2,2,2,2),48),
    ((1,2,1,2),3), ((1,1,1,1,2),3),
    ((1,1,1,1,1,1,1,1,6),15)
]
for ns,T in boundaries:
    z = analyze(ns,T)
    print('BOUNDARY',ns,'T',T,'valid/neg/pos',
          z['valid'],z['neg'],z['pos'],
          'demand/capacity/flow/F',
          z['demand'],z['capacity'],z['flow'],z['F'])

ns, T = (1,1,1,1,2), 5
M = ((0,3),(1,5),(2,4))
z = analyze(ns,T)
print('WITNESS',M,'weight',z['weights'][M],
      'terminal expansion',straighten(M))
print('TERMINAL WEIGHTS',
      [(N,c,weight(factors(ns,N),T)) for N,c in straighten(M)])

first = []
for e,f in combinations(M,2):
    if cross(e,f):
        out = defaultdict(int)
        for child in smooth(M,e,f):
            for N,c in straighten(child):
                out[N] += c
        first.append(tuple(sorted((N,c) for N,c in out.items() if c)))
print('FIRST-CROSSING CHOICES',len(first),
      'distinct expansions',len(set(first)))

P = ((0,5),(1,3),(2,4))
print('POSITIVE INTERMEDIATE',P,'weight',weight(factors(ns,P),T))
