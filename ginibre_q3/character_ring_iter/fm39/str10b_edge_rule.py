from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod
import random

def progress(s):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

@lru_cache(None)
def paths(ns):
    if not ns:
        return ((0, ()),)
    return tuple((c, p+(c,)) for a, p in paths(ns[:-1]) for c in cg(a, ns[-1]))

def clean(p):
    return {e: c for e, c in p.items() if c}

def lower(p, ns):
    o = defaultdict(Q)
    for e, c in p.items():
        for i, n in enumerate(ns):
            if e[i] < n:
                f = list(e)
                f[i] += 1
                o[tuple(f)] += c*(n-e[i])
    return clean(o)

def raise_op(p):
    o = defaultdict(Q)
    for e, c in p.items():
        for i, h in enumerate(e):
            if h:
                f = list(e)
                f[i] -= 1
                o[tuple(f)] += c*h
    return clean(o)

def inner(p, q, ns):
    if len(p) > len(q):
        p, q = q, p
    return sum(c*q.get(e, 0)/prod(comb(n, h) for n, h in zip(ns, e))
               for e, c in p.items())

@lru_cache(None)
def cg_basis(ns):
    if not ns:
        return ((0, (), ({(): Q(1)},)),)
    out = []
    n = ns[-1]
    for a, path, old in cg_basis(ns[:-1]):
        for c in cg(a, n):
            j = (a+n-c)//2
            top = defaultdict(Q)
            for h in range(j+1):
                for e, z in old[h].items():
                    top[e+(j-h,)] += ((-1)**h)*comb(j, h)*z
            top = clean(top)
            assert not raise_op(top)
            states = [top]
            for _ in range(c):
                nxt = lower(states[-1], ns)
                divisor = c-len(states)+1
                states.append({e: z/Q(divisor) for e, z in nxt.items()})
            assert not lower(states[-1], ns)
            out.append((c, path+(c,), tuple(states)))
    assert sum(c+1 for c, _, _ in out) == prod(n+1 for n in ns)
    return tuple(out)

@lru_cache(None)
def half_copies(word):
    ns = tuple(map(abs, word))
    out = defaultdict(lambda: [[], []])
    for mask in range(1 << len(word)):
        ix = [i for i in range(len(word)) if mask >> i & 1]
        iy = [i for i in range(len(word)) if not (mask >> i & 1)]
        parity = sum(word[i] < 0 for i in iy) & 1
        for a, px, _ in cg_basis(tuple(ns[i] for i in ix)):
            for b, py, _ in cg_basis(tuple(ns[i] for i in iy)):
                out[(a, b)][parity].append((mask, px, py))
    for e, o in out.values():
        e.sort()
        o.sort()
    return {k: (tuple(e), tuple(o)) for k, (e, o) in out.items()}

def unmatched(copies):
    out = {}
    for ch, (even, odd) in copies.items():
        m = min(len(even), len(odd))
        kept = {}
        if len(even) > m:
            kept[0] = tuple((i, even[i]) for i in range(m, len(even)))
        if len(odd) > m:
            kept[1] = tuple((i, odd[i]) for i in range(m, len(odd)))
        if kept:
            out[ch] = kept
    return out

def tensor(p, q):
    o = defaultdict(Q)
    for x, a in p.items():
        for y, b in q.items():
            o[tuple(i+j for i, j in zip(x, y))] += a*b
    return clean(o)

def embed(p, ix, total):
    out = {}
    for e, c in p.items():
        v = [0]*total
        for i, h in zip(ix, e):
            v[i] = h
        out[tuple(v)] = c
    return out

def states_for(ns, path):
    return next(s for _, p, s in cg_basis(ns) if p == path)

def pair_invariant(ns, left, right, lpath, rpath, spin):
    out = defaultdict(Q)
    L = states_for(tuple(ns[i] for i in left), lpath)
    R = states_for(tuple(ns[i] for i in right), rpath)
    total = len(ns)
    for h in range(spin+1):
        p = tensor(embed(L[h], left, total), embed(R[spin-h], right, total))
        factor = (-1)**h*comb(spin, h)
        for e, c in p.items():
            out[e] += factor*c
    return clean(out)

def harmonic_vector(word, A, B, ch, ca, cb):
    ns = tuple(map(abs, word))
    maska, pax, pay = ca
    maskb, pbx, pby = cb
    ax = tuple(A[j] for j in range(len(A)) if maska >> j & 1)
    ay = tuple(A[j] for j in range(len(A)) if not (maska >> j & 1))
    bx = tuple(B[j] for j in range(len(B)) if maskb >> j & 1)
    by = tuple(B[j] for j in range(len(B)) if not (maskb >> j & 1))
    xp = pair_invariant(ns, ax, bx, pax, pbx, ch[0])
    yp = pair_invariant(ns, ay, by, pay, pby, ch[1])
    v = tensor(xp, yp)
    assert v and not raise_op(v) and not lower(v, ns)
    ym = sum(1 << i for i in ay+by)
    return ym, v

def harmonic_basis(word, parity):
    A = tuple(range(0, len(word), 2))
    B = tuple(range(1, len(word), 2))
    wa = tuple(word[i] for i in A)
    wb = tuple(word[i] for i in B)
    raw_a = half_copies(wa)
    raw_b = half_copies(wb)
    ca = unmatched(raw_a)
    cb = unmatched(raw_b)
    out = []
    for ch in sorted(ca.keys() & cb.keys()):
        for pa, ra in ca[ch].items():
            for pb, rb in cb[ch].items():
                if pa ^ pb != parity:
                    continue
                for ia, xa in ra:
                    for ib, xb in rb:
                        mask, v = harmonic_vector(word, A, B, ch, xa, xb)
                        out.append((ch, (ia, ib), mask, v))
    return out

def cut_weight(S, T):
    common = S & T
    val = 1
    i = 0
    while common:
        if common & 1:
            val *= i+2
        common >>= 1
        i += 1
    return val

@lru_cache(None)
def gate(word, S, T):
    n = len(word)
    C = tuple(i for i in range(n) if S >> i & 1 and T >> i & 1)
    D = tuple(i for i in range(n) if not ((S | T) >> i & 1))
    R = tuple(i for i in range(n) if S >> i & 1 and not (T >> i & 1))
    Qs = tuple(i for i in range(n) if T >> i & 1 and not (S >> i & 1))
    Fs = []
    for I in (C, D, R, Qs):
        Fs.append({c for c, _ in paths(tuple(abs(word[i]) for i in I))})
    return bool(set.intersection(*Fs))

def rational_rank(M):
    a = [[Q(x) for x in row] for row in M]
    if not a:
        return 0
    r = 0
    for j in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][j]
        a[r] = [v/z for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                z = a[i][j]
                a[i] = [v-z*w for v, w in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

class Dinic:
    def __init__(self, n):
        self.g = [[] for _ in range(n)]

    def add(self, u, v, c):
        a = [v, c, len(self.g[v])]
        b = [u, 0, len(self.g[u])]
        self.g[u].append(a)
        self.g[v].append(b)

    def flow(self, s, t):
        total = 0
        n = len(self.g)
        while True:
            level = [-1]*n
            level[s] = 0
            q = deque([s])
            while q:
                u = q.popleft()
                for v, c, rev in self.g[u]:
                    if c and level[v] < 0:
                        level[v] = level[u]+1
                        q.append(v)
            if level[t] < 0:
                return total
            it = [0]*n

            def dfs(u, f):
                if u == t:
                    return f
                while it[u] < len(self.g[u]):
                    e = self.g[u][it[u]]
                    v, c, rev = e
                    if c and level[v] == level[u]+1:
                        z = dfs(v, min(f, c))
                        if z:
                            e[1] -= z
                            self.g[v][rev][1] += z
                            return z
                    it[u] += 1
                return 0

            while True:
                z = dfs(s, 10**100)
                if not z:
                    break
                total += z

def capflow(sc, tc, G):
    L = sorted(sc)
    R = sorted(tc)
    s = 0
    lo = 1
    ro = lo+len(L)
    t = ro+len(R)
    d = Dinic(t+1)
    inf = sum(sc.values())
    for i, a in enumerate(L):
        d.add(s, lo+i, sc[a])
        for b in G.get(a, ()):
            if b in tc:
                d.add(lo+i, ro+R.index(b), inf)
    for j, b in enumerate(R):
        d.add(ro+j, t, tc[b])
    return d.flow(s, t)

def graph_data(word):
    E = harmonic_basis(word, 0)
    O = harmonic_basis(word, 1)
    ns = tuple(map(abs, word))
    M = [[cut_weight(S, T)*inner(vs, vt, ns)
          for _, _, T, vt in O] for _, _, S, vs in E]
    sc = dict(sorted(Counter(ch for ch, _, _, _ in O).items()))
    tc = dict(sorted(Counter(ch for ch, _, _, _ in E).items()))
    G = defaultdict(set)
    compatible = zeros = 0
    for i, (beta, _, T, vt) in enumerate(E):
        for j, (alpha, _, S, vs) in enumerate(O):
            ok = gate(word, S, T)
            if ok:
                compatible += 1
            if M[i][j]:
                assert ok
                G[alpha].add(beta)
            elif ok:
                zeros += 1
    G = {a: tuple(sorted(G[a])) for a in sorted(sc)}
    rk = rational_rank(M)
    fl = capflow(sc, tc, G)
    one = all(sc[a] <= sum(tc.get(b, 0) for b in G[a]) for a in sc)
    rad = max((max(abs(a[0]-b[0]), abs(a[1]-b[1]))
               for a in sc for b in G[a]), default=0)
    return E, O, M, sc, tc, G, rk, fl, one, rad, compatible, zeros

def dec(d):
    return {(int(a[0]), int(a[1])):
            {(int(b[0]), int(b[1])) for b in s.split()}
            for a, s in d.items()}

reps = [
    (-1,-1,-1,-2,-2,-2,3),
    (-1,-2,-2,-2,3,3,3),
    (1,1,1,1,1,-2,-2,3),
    (1,1,1,-2,-2,3,3,3),
    (1,-2,-2,3,3,3,3,3)
]
partners = [
    (1,1,1,-2,-2,-2,-3),
    (1,-2,-2,-2,-3,-3,-3),
    (-1,-1,-1,-1,-1,-2,-2,-3),
    (-1,-1,-1,-2,-2,-3,-3,-3),
    (-1,-2,-2,-3,-3,-3,-3,-3)
]
G5 = [
    dec({'12':'01 03 10 23 30 32 41',
         '21':'01 03 10 14 23 30 32'}),
    dec({'23':'01 03 05 10 12 21 25 30 34 43 50 52',
         '32':'01 03 05 10 12 21 25 30 34 43 50 52'}),
    dec({'23':'03 05 10 12 14 21 30 41 50',
         '32':'01 03 05 10 12 14 21 30 41 50'}),
    dec({'34':'01 03 05 07 10 12 16 21 30 32 50 61 70',
         '43':'01 03 05 07 10 12 16 21 30 32 50 61 70'}),
    dec({'12':'01 03 10 23 30 32 34 43 45 54',
         '21':'01 03 10 23 30 32 34 43 45 54',
         '16':'05 07 23 27 32 34 36 43 45 50 54 63 70 72',
         '61':'05 07 23 27 32 34 36 43 45 50 54 63 70 72'})
]
extras = [
    (1,-2,-2,3,-4,-4),
    (1,2,2,3,-4,-4),
    (-1,2,2,-4,-4,-5),
    (1,-2,-4,-4,-4,5),
    (2,-3,-3,4,5,5),
    (-2,-3,4,5,5,5)
]
Gx = [
    dec({'12':'01 03 10 30','21':'01 03 10 30'}),
    dec({'12':'01 03 10 14 30 34 41 43','21':'01 03 10 14 30 34 41 43'}),
    dec({'12':'01 03 10 30 34 43','21':'01 03 10 30 34 43'}),
    dec({'12':'01 03 10 30 34 43','21':'01 03 10 30 34 43'}),
    dec({'24':'02 04 06 15 20 35 37 40 51 53 60 73',
         '42':'02 04 06 15 20 35 37 40 51 53 60 73'}),
    dec({'23':'01 03 05 10 25 30 50 52 56 65',
         '32':'01 03 05 10 25 30 50 52 56 65'})
]
expected = {w:g for pair,g in zip(zip(reps, partners), G5) for w in pair}
expected.update(dict(zip(extras, Gx)))

more = [
    (-1,2,3,3,4,-5),(-1,-2,-3,-4,-5,-5),(-1,3,3,4,4,-5),
    (1,1,1,-2,4,-5),(-1,-2,3,3,3,4),(1,-2,3,3,-4,5),
    (2,-3,-3,-3,-4,5),(-1,-3,-3,-4,-4,-5),(1,2,2,-4,-4,5),
    (1,-2,-2,3,3,5),(1,-2,3,3,4,-5),(1,-4,-4,5,5,5),
    (-1,2,3,-4,-5,-5),(-1,3,3,-4,-4,-5),(1,1,2,-3,4,-5),
    (1,1,-2,3,4,-5),(-1,-2,-3,-3,-4,-5),(-2,-3,-4,-5,-5,-5),
    (2,3,3,4,-5,-5),(1,3,3,-4,-4,5),(1,1,1,-2,-3,4),
    (-2,-3,-3,-3,4,5),(-1,-3,4,4,5,5),(1,3,3,3,-4,-4),
    (2,3,-4,-5,-5,-5),(2,3,3,3,-4,-5),(-1,-1,-1,-2,3,4),
    (-1,-1,2,3,4,5),(1,-2,3,-4,5,5),(1,1,1,-2,3,-4)
]
assert len(expected) == 16 and len(more) == 30

for idx, w in enumerate(list(reps)+list(partners)+extras+more, 1):
    progress(f'graph {idx}/46 start {w}')
    E,O,M,sc,tc,G,rk,fl,one,rad,ng,nz = graph_data(w)
    assert rk == fl == len(O) and one
    if w in expected:
        assert {a:set(v) for a,v in G.items()} == expected[w]
    if w == reps[0]:
        assert (ng, nz) == (262, 110)
    progress(f'graph {idx}/46 H=({len(E)},{len(O)}) rank=flow={fl} '
             f'one_vertex={one} edges={sum(map(len,G.values()))} '
             f'radius={rad} gate={ng} compatible_zero={nz}')
progress('PASS 46 exact support graphs; first 16 matched recorded adjacency')

w = reps[0]
E,O,M,sc,tc,G,rk,fl,one,rad,ng,nz = graph_data(w)
sid = next(i for i,x in enumerate(O) if x[0] == (1,2) and i == 1)
tid = next(i for i,x in enumerate(E) if x[0] == (0,1) and i == 0)
assert O[sid][2] == 13 and E[tid][2] == 127
assert M[tid][sid] == 0 and cut_weight(13,127) == 40 and gate(w,13,127)
progress('gate-zero witness: source (1,2), copy 1, mask 13; target (0,1), '
         'copy 0, mask 127; weight 40, overlap 0, common-spin gate true')

def channel_table(word):
    d = {(0,0,0):1}
    for signed in word:
        n = abs(signed)
        q = defaultdict(int)
        for (a,b,p), v in d.items():
            for c in cg(a,n):
                q[(c,b,p)] += v
            for c in cg(b,n):
                q[(a,c,p ^ (signed < 0))] += v
        d = q
    even = defaultdict(int)
    odd = defaultdict(int)
    for (a,b,p), v in d.items():
        (odd if p else even)[(a,b)] += v
    return {k:even[k]-odd[k] for k in even.keys() | odd.keys()}

def dimensions(word):
    A = channel_table(word[::2])
    B = channel_table(word[1::2])
    vals = [A.get(k,0)*B.get(k,0) for k in A.keys() | B.keys()]
    return (sum(max(v,0) for v in vals),
            sum(max(-v,0) for v in vals),
            sum(v > 0 for v in vals),
            sum(v < 0 for v in vals))

for k in range(1,21):
    sg = 1 if k % 2 == 0 else -1
    w = tuple([-i for i in range(1,k+1)] + [sg*(k+1)])
    hp, hm, np, nm = dimensions(w)
    assert hp >= hm
    progress(f'run k={k} plus={hp} minus={hm} gap={hp-hm} '
             f'channel_counts=({np},{nm})')

special = [
    ('sign-min',(-1,2,3,4,-5,6,7,8)),
    ('F1',(1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8))
]
for name, w in special:
    hp, hm, np, nm = dimensions(w)
    assert hp >= hm
    progress(f'{name} W={sum(map(abs,w))} plus={hp} minus={hm} '
             f'gap={hp-hm} channel_counts=({np},{nm})')

for t in range(1,13):
    w = (1,)*(2*t)+(-2,4,-6)
    hp, hm, np, nm = dimensions(w)
    assert hp >= hm
    progress(f'suffix t={t} factors={len(w)} plus={hp} minus={hm} gap={hp-hm}')

rng = random.Random(20261002)
samples = []
seen = set()
draws = 0
while len(samples) < 100:
    draws += 1
    L = rng.randint(4,12)
    nums = [rng.randint(1,14) for _ in range(L)]
    ct = Counter(nums)
    sg = {n:(-1 if rng.randrange(2) else 1) for n in ct}
    if sum(ct[n] for n in ct if sg[n] < 0) % 2:
        continue
    w = tuple(sg[n]*n for n in sorted(ct) for _ in range(ct[n]))
    if sum(map(abs,w)) > 80 or w in seen:
        continue
    seen.add(w)
    hp, hm, np, nm = dimensions(w)
    if hm:
        samples.append((hp-hm,w,hp,hm))
assert len(samples) == 100
assert min(samples) == (40,(1,-2,-3,4,10,12),42,2)
progress(f'conditioned random pair-free W<=80 count=100 draws={draws} '
         f'minimum={min(samples)}')

for label, w in [
    ('run-k6',tuple([-i for i in range(1,7)]+[7])),
    ('suffix-t3',(1,)*6+(-2,4,-6))
]:
    progress(f'{label} full graph start')
    E,O,M,sc,tc,G,rk,fl,one,rad,ng,nz = graph_data(w)
    assert rk == fl == len(O) and one
    progress(f'{label} H=({len(E)},{len(O)}) rank=flow={fl} '
             f'source_classes={len(sc)} edges={sum(map(len,G.values()))} radius={rad}')
progress('PASS')
