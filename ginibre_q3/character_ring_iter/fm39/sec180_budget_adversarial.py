from collections import defaultdict, Counter
from fractions import Fraction
from functools import lru_cache
from datetime import datetime, timezone
import random

def stamp(s):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def mult(F, z, cap=None):
    n = abs(z)
    eps = 1 if z > 0 else -1
    d = len(F)-1
    D = d+n if cap is None else min(d+n, cap)
    G = [[0]*(D+1) for _ in range(D+1)]
    for ax in (0, 1):
        for v in range(min(d, D)+1):
            last = d-v
            P = [0]*(last+3)
            for u in range(last+1):
                P[u+2] = P[u] + (F[u][v] if ax == 0 else F[v][u])
            for t in range(D-v+1):
                lo = abs(t-n)
                hi = min(last, t+n)
                if (hi-lo) % 2:
                    hi -= 1
                if lo <= hi:
                    q = P[hi+2]-P[lo]
                    if ax == 0:
                        G[t][v] += q
                    else:
                        G[v][t] += eps*q
    return G

@lru_cache(maxsize=4096)
def row(w, cap):
    w = tuple(sorted(w, key=lambda z: (-abs(z), z)))
    rem = sum(map(abs, w))
    F = [[1]]
    for z in w:
        rem -= abs(z)
        F = mult(F, z, None if cap < 0 else rem+cap)
    return tuple(tuple(x) for x in F)

def cut(C):
    A = []
    B = []
    wa = wb = 0
    for z in sorted(C, key=lambda z: -abs(z)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    return tuple(A), tuple(B)

def ratio_min(inc, deninc):
    n = d = 0
    best = None
    for t, (x, y) in enumerate(zip(inc, deninc)):
        n += x
        d += y
        if d:
            q = Fraction(n, d)
            if best is None or (q, t) < (best[0], best[1]):
                best = (q, t, n, d)
    return best

@lru_cache(maxsize=20000)
def profile(w, i, j):
    u, v = w[i], w[j]
    C = tuple(z for k, z in enumerate(w) if k not in (i, j))
    A, B = cut(C)
    cap = sum(map(abs, A))
    X = row(A, -1)
    Y = row(B+(u, v), cap)
    Z = row(B+(-u, -v), cap)
    P = []
    M = []
    L = []
    for t in range(cap+1):
        p = m = ell = 0
        for a in range(t+1):
            if a < len(X) and t-a < len(X[a]):
                x = X[a][t-a]
                y = Y[a][t-a]
                z = Z[a][t-a]
                assert (y+z) % 2 == 0 and (y-z) % 2 == 0
                p += x*((y+z)//2)
                m += x*((y-z)//2)
                ell += x*y
        P.append(p)
        M.append(m)
        L.append(ell)
    Q = [p+min(m, 0) for p, m in zip(P, M)]
    plain = ratio_min(L, [max(x, 0) for x in L])
    b_pure = ratio_min(Q, [max(x, 0) for x in P])
    b_net = ratio_min(Q, [max(x, 0) for x in Q])
    neg = sum(max(-x, 0) for x in M)
    pos = sum(max(x, 0) for x in P)
    minPi = min([0]+list(__import__('itertools').accumulate(L)))
    minB = min([0]+list(__import__('itertools').accumulate(Q)))
    return (A, B, cap, tuple(P), tuple(M), tuple(L), tuple(Q),
            plain, b_pure, b_net,
            Fraction(neg, pos) if pos else None, neg, pos, minPi, minB)

def top_pair(w):
    B = w[:-1]
    pairs = [(i, j) for i in range(len(B)) for j in range(i+1, len(B))
             if (abs(B[i])+abs(B[j])) % 2 == 0]
    return max(pairs, key=lambda ij:
               (abs(B[ij[0]])+abs(B[ij[1]]),
                max(abs(B[ij[0]]), abs(B[ij[1]])))) if pairs else None

def canonical(w):
    return tuple(sorted(w[:-1], key=lambda z: (abs(z), z)))+(w[-1],)

def valid(w):
    if len(w) < 4 or len(w) > 28 or sum(map(abs, w)) > 120 or any(z == 0 for z in w):
        return False
    signs = {}
    for z in w:
        n = abs(z)
        e = 1 if z > 0 else -1
        if n in signs and signs[n] != e:
            return False
        signs[n] = e
    return sum(z < 0 for z in w) % 2 == 0 and top_pair(w) is not None

def mutate(w, rng):
    classes = sorted(set(map(abs, w)))
    for _ in range(20):
        B = list(w[:-1])
        p = w[-1]
        v = list(w)
        op = rng.randrange(6)
        if op == 0 and classes:
            n = rng.choice(classes)
            nn = n+rng.choice((-2, -1, 1, 2))
            if nn < 1:
                continue
            v = [(nn if z > 0 else -nn) if abs(z) == n else z for z in v]
        elif op == 1 and len(classes) >= 2:
            a, b = rng.sample(classes, 2)
            ma = sum(abs(z) == a for z in v)
            mb = sum(abs(z) == b for z in v)
            if (ma+mb) % 2:
                continue
            v = [-z if abs(z) in (a, b) else z for z in v]
        elif op == 2:
            plus = [i for i, z in enumerate(B) if z > 0]
            if plus and rng.random() < 0.55:
                B.pop(rng.choice(plus))
            else:
                n = rng.randint(1, 30)
                if any(abs(z) == n and z < 0 for z in v):
                    continue
                B.append(n)
            v = B+[p]
        elif op == 3:
            opts = [n for n in classes
                    if any(abs(z) == n and z < 0 for z in v)
                    and sum(z == -n for z in B) >= 2]
            if opts and rng.random() < 0.5:
                n = rng.choice(opts)
                B.remove(-n)
                B.remove(-n)
            else:
                n = rng.randint(1, 30)
                if any(abs(z) == n and z > 0 for z in v):
                    continue
                B.extend([-n, -n])
            v = B+[p]
        elif op == 4:
            if len(B) < 2:
                continue
            ii, jj = rng.sample(range(len(B)), 2)
            x, y = B[ii], B[jj]
            if (x > 0) != (y > 0):
                continue
            znew = abs(x)+abs(y)
            B = [z for k, z in enumerate(B) if k not in (ii, jj)]+[znew]
            v = B+[p]
        else:
            opts = [i for i, z in enumerate(B) if z > 0 and abs(z) >= 2]
            if not opts:
                continue
            i = rng.choice(opts)
            n = B.pop(i)
            a = rng.randint(1, n-1)
            b = n-a
            B.extend([a, b] if rng.random() < 0.5 else [-a, -b])
            v = B+[p]
        v = canonical(v)
        if valid(v):
            return v
    return w

def direct(w):
    d = {(0, 0): 1}
    for z in w:
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (a, b), c in d.items():
            for x in cg(a, n):
                q[(x, b)] += c
            for y in cg(b, n):
                q[(a, y)] += eps*c
        d = dict(q)
    return d

@lru_cache(maxsize=10000)
def phi(w):
    return direct(w).get((0, 0), 0)

def is_noflip(w):
    p = phi(w)
    for i in range(len(w)):
        for j in range(i+1, len(w)):
            z = list(w)
            z[i] = -z[i]
            z[j] = -z[j]
            if p-phi(tuple(z)) >= 0:
                return False
    return True

def verify(w, i, j, pr):
    A, B, cap, P, M, L, Q, *_ = pr
    X = direct(A)
    Y = direct(B+(w[i], w[j]))
    Z = direct(B+(-w[i], -w[j]))
    FX = row(A, -1)
    FY = row(B+(w[i], w[j]), cap)
    FZ = row(B+(-w[i], -w[j]), cap)
    for U, D, c in ((FX, X, -1), (FY, Y, cap), (FZ, Z, cap)):
        for a in range(len(U)):
            for b in range(len(U[a])):
                if c < 0 or a+b <= c:
                    assert U[a][b] == D.get((a, b), 0)
    p2 = []
    m2 = []
    l2 = []
    for t in range(cap+1):
        p = m = ell = 0
        for a in range(t+1):
            x = X.get((a, t-a), 0)
            y = Y.get((a, t-a), 0)
            z = Z.get((a, t-a), 0)
            p += x*((y+z)//2)
            m += x*((y-z)//2)
            ell += x*y
        p2.append(p)
        m2.append(m)
        l2.append(ell)
    assert tuple(p2) == P and tuple(m2) == M and tuple(l2) == L
    assert sum(L) == phi(w)

def score(w, obj):
    ij = top_pair(w)
    if ij is None:
        return None
    q = profile(w, *ij)
    m = {'plain': q[7], 'budget': q[8],
         'budget_net': q[9], 'mixed': q[10]}[obj]
    if m is None or (isinstance(m, tuple) and not m): return None
    return m[0] if isinstance(m, tuple) else m

def better(a, b, obj):
    return a > b if obj == 'mixed' else a < b

def pair_best(w):
    out = None
    for i in range(len(w)):
        for j in range(i+1, len(w)):
            q = profile(w, i, j)
            m = q[8]
            if m is None:
                continue
            if out is None or m[0] > out[0]:
                out = (m[0], i, j, q)
    return out

# Census records sampled with seed 180 and embedded for a self-contained run.
census = [
    (-1,4,5,6,6,7,8,-13),(-1,-3,-4,-7,10,12,13),
    (-1,-2,-3,-3,-4,-5,-5,-5,-6,-14),(1,2,-3,4,4,-5,7,10,14),
    (-2,3,-4,5,-9,10,-11),(-1,-5,10,11,12,13),
    (1,2,-3,6,-7,9,-11,-15),(1,-3,4,-5,7,-9,10,-11),
    (1,-2,3,7,-8,9,-10,-10),(-2,-3,-5,-8,9,12,15),
    (1,2,-3,4,6,11,-13,14),(-1,-2,-3,-4,-6,-7,-8,-9,22),
    (-1,4,5,8,-9,-12,-13),(-1,-4,6,7,9,10,11),
    (-1,-2,-3,-3,-4,-5,-5,8,9,-22),(1,2,4,-5,8,-9,10,13)
]
tight = (-1,-3,-3,-4,-5,-5,-6,-7)
F1 = (1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)
rng = random.Random(180)
seeds = [canonical(tight), canonical(F1)]+[canonical(w) for w in census]
for k in range(5, 21):
    s = -1 if k % 2 else 1
    for p in range(k+1, k+11):
        w = tuple(-n for n in range(1, k+1))+(s*p,)
        if valid(w):
            seeds.append(canonical(w))

def randword(r):
    for _ in range(100):
        q = r.randint(4, 9)
        labs = r.sample(range(1, 36), q)
        cnt = {n:r.randint(1, 3) for n in labs}
        sg = {n:(-1 if r.randrange(2) else 1) for n in labs}
        if sum(cnt[n] for n in labs if sg[n] < 0) % 2:
            sg[labs[0]] *= -1
        z = [sg[n]*n for n in labs for _ in range(cnt[n])]
        if sum(map(abs, z)) > 120:
            continue
        pi = r.randrange(len(z))
        p = z.pop(pi)
        w = canonical(tuple(z)+(p,))
        if valid(w):
            return w
    return canonical(tight)

seeds += [randword(rng) for _ in range(40)]
seeds = list(dict.fromkeys(seeds))
rng = random.Random(180180)
pool = {}

def record(w):
    w = canonical(w)
    if not valid(w):
        return None
    if w not in pool:
        ij = top_pair(w)
        if ij is None:
            return None
        pool[w] = (ij, profile(w, *ij))
    return pool[w]

for w in seeds:
    record(w)

objs = ('plain', 'budget', 'budget_net', 'mixed')
best = {}
bins = {}
for obj in objs:
    eligible = [w for w in pool if score(w, obj) is not None]
    bw = (max(eligible, key=lambda w:score(w, obj)) if obj == 'mixed'
          else min(eligible, key=lambda w:score(w, obj)))
    best[obj] = (score(bw, obj), bw)
    for restart in range(8):
        cur = seeds[(restart*7+len(obj)) % len(seeds)]
        cv = score(cur, obj)
        for step in range(12):
            cand = []
            for _ in range(2):
                nw = mutate(cur, rng)
                record(nw)
                nv = score(nw, obj)
                if nv is not None:
                    cand.append((nv, nw))
            if not cand:
                continue
            nv, nw = max(cand) if obj == 'mixed' else min(cand)
            if cv is None or better(nv, cv, obj) or rng.randrange(7) == 0:
                cur, cv = nw, nv
            if better(nv, best[obj][0], obj):
                best[obj] = (nv, nw)
            W = sum(map(abs, nw))
            b = min(5, (W-1)//20)
            old = bins.get((obj, b))
            if old is None or better(nv, old[0], obj):
                bins[(obj, b)] = (nv, nw)
        stamp(f'top hill {obj} restart={restart+1}/8 pool={len(pool)} best={best[obj]}')

for _ in range(80):
    record(randword(rng))
for w in list(pool):
    for obj in objs:
        nv = score(w, obj)
        if nv is None:
            continue
        W = sum(map(abs, w))
        b = min(5, (W-1)//20)
        old = bins.get((obj, b))
        if old is None or better(nv, old[0], obj):
            bins[(obj, b)] = (nv, w)
stamp(f'top search pool={len(pool)}')

dpool = list(dict.fromkeys(
    seeds[:24]+[w for _, w in best.values()]+[w for _, w in bins.values()]
))
dpool += [randword(rng) for _ in range(28)]
drows = []
for k, w in enumerate(dict.fromkeys(dpool)):
    q = pair_best(w)
    if q is not None:
        drows.append((q[0], w, q[1], q[2], q[3]))
    if (k+1) % 10 == 0:
        stamp(f'all-pair checks {k+1}/{len(dpool)}')
dbest = min(drows, key=lambda x:x[0])
stamp(f'all-pair candidates={len(drows)} best={dbest[0]} '
      f'word={dbest[1]} pair={(dbest[1][dbest[2]], dbest[1][dbest[3]])}')

shown = {}
for obj, (v, w) in best.items():
    ij, p = pool[w]
    verify(w, *ij, p)
    shown[obj] = (v, w, ij, p, is_noflip(w))
for w in [x[1] for x in bins.values()]+[dbest[1]]:
    if w not in pool:
        record(w)
    ij, p = pool[w]
    verify(w, *ij, p)
verify(dbest[1], dbest[2], dbest[3], dbest[4])

for w, (ij, p) in pool.items():
    if p[13] < 0 or p[14] < 0:
        verify(w, *ij, p)
        print('DIRECTLY_CONFIRMED_NEGATIVE', w, p[13], p[14])
    assert p[13] >= 0 and p[14] >= 0

for obj, (v, w, ij, p, nf) in shown.items():
    print('EXTREME', obj, 'score', v, 'word', w, 'W', sum(map(abs, w)),
          'TopPair', (w[ij[0]], w[ij[1]]), 'noflip', nf,
          'plain', p[7], 'budget_pure', p[8], 'budget_net', p[9],
          'mixed', p[10], 'negative_mixed', p[11],
          'positive_pure', p[12], 'minPi', p[13], 'minB', p[14])
for key, (v, w) in sorted(bins.items()):
    print('BIN', key, 'score', v, 'word', w, 'W', sum(map(abs, w)))
for w in dict.fromkeys(dpool):
    q = pair_best(w)
    if q is not None and q[1] == dbest[2] and w == dbest[1]:
        assert q[3][14] >= 0
print('ALL_PAIR_EXTREME', dbest[0], dbest[1], 'W', sum(map(abs, dbest[1])),
      'pair', (dbest[1][dbest[2]], dbest[1][dbest[3]]),
      'profile', dbest[4][8], 'raw_min_B', dbest[4][14],
      'noflip', is_noflip(dbest[1]))

ij = top_pair(tight)
q = profile(tight, *ij)
verify(tight, *ij, q)
assert q[10] == Fraction(131, 406)
print('TIGHT', tight, 'TopPair', (tight[ij[0]], tight[ij[1]]),
      'mix', q[10], q[11], q[12], 'plain', q[7],
      'budget', q[8], 'noflip', is_noflip(tight))
print('CENSUS_SEED_NOFLIP', sum(is_noflip(w) for w in census), '/', len(census))

zero = (-1,-2,-3,-4,-5,-6)
zi = top_pair(zero)
zp = profile(zero, *zi)
assert all(x == 0 for seq in zp[3:7] for x in seq)
stamp(f'zero-mass ratio excluded: {zero}, TopPair={(zero[zi[0]],zero[zi[1]])}')

print('PASS')
