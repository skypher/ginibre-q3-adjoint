from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
import random, sys

if any(a in ('-h', '--help') for a in sys.argv[1:]):
    print('Usage: python3 -u - [no args]; exact witness/F1 and 30 seeded pair-free random split checks.')
    raise SystemExit(0)

def stamp(s):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def table(w):
    d = {(0, 0): 1}
    for z in sorted(w, key=lambda q: (-abs(q), q)):
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), v in d.items():
            for c in cg(r, n):
                q[c, s] += v
            for c in cg(s, n):
                q[r, c] += eps*v
        d = {k: v for k, v in q.items() if v}
    return d

def laurent_irreps(w):
    d = {(0, 0): 1}
    for z in w:
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), v in d.items():
            for a in range(-n, n+1, 2):
                q[r+a, s] += v
            for b in range(-n, n+1, 2):
                q[r, s+b] += eps*v
        d = dict(q)
    N = sum(map(abs, w))
    out = {}
    for r in range(N, -1, -1):
        for s in range(N, -1, -1):
            v = (d.get((r, s), 0) - d.get((r+2, s), 0)
                 - d.get((r, s+2), 0) + d.get((r+2, s+2), 0))
            if v:
                out[r, s] = v
    return out

def canon(w):
    return tuple(sorted(w, key=lambda z: (abs(z), z)))

def state(A, B):
    A, B = canon(A), canon(B)
    return (A, B) if A <= B else (B, A)

def all_splits(w):
    c = Counter(w)
    typ = tuple(sorted(c, key=lambda z: (abs(z), z)))
    seen = set()
    for alloc in product(*(range(c[z]+1) for z in typ)):
        A = tuple(z for z, k in zip(typ, alloc) for _ in range(k))
        B = tuple(z for z, k in zip(typ, alloc) for _ in range(c[z]-k))
        if A and B:
            seen.add(state(A, B))
    return sorted(seen)

def profile(st):
    A, B = st
    X, Y = table(A), table(B)
    layer = defaultdict(int)
    for (r, s), v in X.items():
        layer[r+s] += v*Y.get((r, s), 0)
    run = mass = 0
    neg = []
    best = None
    for t, val in sorted(layer.items()):
        run += val
        mass += max(val, 0)
        if run < 0:
            neg.append((t, run))
        if mass:
            q = (Fraction(run, mass), t, run, mass)
            if best is None or q[:2] < best[:2]:
                best = q
    return tuple(sorted(layer.items())), neg, best

def independent(st):
    X, Y = laurent_irreps(st[0]), laurent_irreps(st[1])
    assert X == table(st[0]) and Y == table(st[1])
    layer = defaultdict(int)
    for key, v in X.items():
        layer[sum(key)] += v*Y.get(key, 0)
    return tuple(sorted(layer.items()))

def balance(st):
    A, B = st
    ca, cb = sum(z < 0 for z in A), sum(z < 0 for z in B)
    wa, wb = sum(-z for z in A if z < 0), sum(-z for z in B if z < 0)
    ta, tb = sum(map(abs, A)), sum(map(abs, B))
    return (abs(ca-cb), abs(wa-wb), abs(ta-tb),
            tuple(sorted((ca, cb))), tuple(sorted((wa, wb))))

witnesses = {
    'W60a': (1,-2,-4,-4,-4,1,3,3,3,3,3,3,5,5,5,5,6),
    'W62': (1,-2,-4,-4,-4,1,3,3,3,3,3,5,5,5,5,5,6),
    'W60b': (1,-2,-4,-4,-4,1,1,2,3,3,3,3,3,4,4,5,5,8)
}
F1 = (1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)
summaries, all_fails = {}, []

for name, w in [*witnesses.items(), ('F1', F1)]:
    stamp('start '+name)
    splits = all_splits(w)
    fails, best, phi = [], None, None
    for st in splits:
        lay, neg, b = profile(st)
        phi = sum(v for _, v in lay)
        if neg:
            assert lay == independent(st), ('Laurent cross-check failed', name, st)
            fails.append((st, neg[0], balance(st)))
        if b is not None and (best is None or b[:2] < best[0][:2]):
            best = (b, st)
    summaries[name] = (len(splits), fails, best, phi)
    print('SUMMARY', name, 'factors=', len(w), 'weight=', sum(map(abs, w)),
          'splits=', len(splits), 'failures=', len(fails), 'Phi=', phi,
          'minimum=', best, flush=True)
    for st, (T, Pi), bal in fails:
        print('FAIL', name, 'A=', st[0], 'B=', st[1], 'T=', T,
              'Pi=', Pi, 'balance=', bal, flush=True)
        all_fails.append((name, st, bal))
    if name == 'F1':
        assert len(splits) == 599 and not fails and phi == 8150742
        assert best[1] == state((-2,-4,-4,-4),
                               tuple(z for z in F1 if z not in (-2,-4,-4,-4)))
        assert best[0] == (Fraction(2756159,19046176),12,5512318,38092352)
        f1tight = best[1]
        print('F1_TIGHT_BALANCE', balance(f1tight),
              'layers=', profile(f1tight)[0], flush=True)

assert [len(summaries[k][1]) for k in ('W60a','W62','W60b','F1')] == [2,2,3,0]
assert all_fails[6][2] == balance(f1tight) == (4,14,28,(0,4),(0,14))

rng = random.Random(182182)
random_splits = random_fails = random_empty = 0
random_best = None
for ix in range(30):
    for tries in range(1000):
        k = rng.randint(5, 9)
        labs = sorted(rng.sample(range(1, 13), k))
        sig = {n: (-1 if rng.randrange(2) else 1) for n in labs}
        mult = {n: rng.randint(1, 2) for n in labs}
        if sum(mult[n] for n in labs if sig[n] < 0) % 2:
            sig[labs[0]] *= -1
        w = tuple(sig[n]*n for n in labs for _ in range(mult[n]))
        if len(w) <= 12 and sum(map(abs, w)) <= 100:
            break
    stamp(f'random list {ix+1}/30 factors={len(w)} weight={sum(map(abs,w))}')
    for st in all_splits(w):
        random_splits += 1
        lay, neg, b = profile(st)
        if neg:
            assert lay == independent(st)
            random_fails += 1
            print('RANDOM_FAIL', ix+1, 'word=', w, 'split=', st,
                  'first=', neg[0], flush=True)
        if b is None:
            random_empty += 1
        elif random_best is None or b[:2] < random_best[0][:2]:
            random_best = (b, w, st)

assert (random_splits, random_fails, random_empty) == (6308,0,4506)
assert random_best == (
    (Fraction(1775,1877),19,3550,3754),
    (2,2,4,6,6,-7,-8,11),
    ((-8,11),(2,2,4,6,6,-7))
)
print('RANDOM_SUMMARY seed=182182 lists=30 distinct_splits=6308 failures=0',
      'empty_mass_splits=4506 minimum=', random_best, flush=True)
