from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from datetime import datetime, timezone
import random, sys

if any(x in ('-h', '--help') for x in sys.argv[1:]):
    print('Exact height-prefix search for signed SU(2) character splits; run with python3 -u - [no arguments].')
    raise SystemExit(0)

def stamp(msg):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), msg, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

@lru_cache(maxsize=50000)
def table(word):
    d = {(0, 0): 1}
    for z in sorted(word, key=lambda q: (-abs(q), q)):
        n = abs(z); eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), v in d.items():
            for c in cg(r, n): q[c, s] += v
            for c in cg(s, n): q[r, c] += eps*v
        d = {k: v for k, v in q.items() if v}
    return d

def cword(w):
    return tuple(sorted(w, key=lambda z: (abs(z), z)))

def cstate(A, B):
    a, b = cword(A), cword(B)
    return (a, b) if a <= b else (b, a)

def whole(st):
    return st[0] + st[1]

def pairfree(w):
    signs = {}
    for z in w:
        n = abs(z); e = 1 if z > 0 else -1
        if n in signs and signs[n] != e: return False
        signs[n] = e
    return True

def valid(st, need_pairfree):
    A, B = st; w = A+B
    return (bool(A) and bool(B) and 2 <= len(w) <= 18 and
            sum(map(abs, w)) <= 120 and sum(z < 0 for z in w) % 2 == 0 and
            (not need_pairfree or pairfree(w)))

def height_profile(st):
    A, B = st
    X, Y = table(A), table(B)
    layer = defaultdict(int)
    for (r, s), x in X.items():
        v = x*Y.get((r, s), 0)
        if v: layer[r+s] += v
    run = mass = 0; rawmin = 0; rawT = None; best = None
    for t in sorted(layer):
        run += layer[t]
        mass += max(layer[t], 0)
        if run < rawmin: rawmin, rawT = run, t
        if mass:
            q = (Fraction(run, mass), t, run, mass)
            if best is None or q[:2] < best[:2]: best = q
    return best, rawmin, rawT, tuple(sorted(layer.items()))

def all_splits(w):
    C = Counter(w); types = tuple(sorted(C, key=lambda z: (abs(z), z)))
    seen = set()
    for alloc in product(*(range(C[z]+1) for z in types)):
        A = tuple(z for z, k in zip(types, alloc) for _ in range(k))
        B = tuple(z for z, k in zip(types, alloc) for _ in range(C[z]-k))
        if not A or not B: continue
        seen.add(cstate(A, B))
    return sorted(seen)

def random_split(w, rng):
    for _ in range(20):
        A = []; B = []
        for z in w:
            (A if rng.randrange(2) else B).append(z)
        if A and B: return cstate(A, B)
    return cstate((w[0],), w[1:])

def random_word(rng, need_pairfree):
    for _ in range(300):
        k = rng.randint(4, 10)
        labels = rng.sample(range(1, 25), k)
        if need_pairfree:
            mult = {n: rng.randint(1, 3) for n in labels}
            sign = {n: (-1 if rng.randrange(2) else 1) for n in labels}
            if sum(mult[n] for n in labels if sign[n] < 0) % 2:
                odd = next(n for n in labels if mult[n] % 2)
                sign[odd] *= -1
            w = tuple(sign[n]*n for n in labels for _ in range(mult[n]))
        else:
            w = tuple(((-1 if rng.randrange(2) else 1)*n)
                      for n in labels for _ in range(rng.randint(1, 2)))
            if sum(z < 0 for z in w) % 2: w = (-w[0],)+w[1:]
        if len(w) <= 18 and sum(map(abs, w)) <= 120 and sum(z < 0 for z in w) % 2 == 0:
            if not need_pairfree or pairfree(w): return cword(w)
    return (-1, -3, -5, -7)

def randparts(n, k, rng):
    cuts = sorted(rng.sample(range(1, n), k-1))
    return tuple(cuts[i] - (cuts[i-1] if i else 0) for i in range(k-1)) + (n-cuts[-1],)

def mutate(st, rng, need_pairfree):
    entries = [(z, 0) for z in st[0]] + [(z, 1) for z in st[1]]
    out = list(entries); op = rng.randrange(8)
    if op == 0 and len(out) >= 2:
        i = rng.randrange(len(out)); z, side = out[i]
        out[i] = (z, 1-side)
    elif op == 1:
        i = rng.randrange(len(out)); z, side = out[i]
        n = abs(z) + rng.choice((-2, -1, 1, 2))
        if n < 1: return st
        out[i] = ((n if z > 0 else -n), side)
    elif op == 2:
        n = abs(rng.choice(out)[0]); shift = rng.choice((-2, -1, 1, 2))
        if n+shift < 1: return st
        sg = rng.choice((-1, 1))
        out = [((z+shift if z > 0 else z-shift), s)
               if abs(z) == n and (1 if z > 0 else -1) == sg else (z, s)
               for z, s in out]
    elif op == 3 and len(out) >= 2:
        i, j = rng.sample(range(len(out)), 2)
        out[i] = (-out[i][0], out[i][1]); out[j] = (-out[j][0], out[j][1])
    elif op == 4:
        plus = [i for i, (z, _) in enumerate(out) if z > 0]
        if plus and rng.randrange(2): out.pop(rng.choice(plus))
        else: out.append((rng.randint(1, 24), rng.randrange(2)))
    elif op == 5:
        neg = [i for i, (z, _) in enumerate(out) if z < 0]
        if len(neg) >= 2 and rng.randrange(2):
            i, j = rng.sample(neg, 2)
            for k in sorted((i, j), reverse=True): out.pop(k)
        else:
            out.extend(((-rng.randint(1, 18), rng.randrange(2)),
                        (-rng.randint(1, 18), rng.randrange(2))))
    elif op == 6:
        plus = [i for i, (z, _) in enumerate(out) if z > 0]
        neg = [i for i, (z, _) in enumerate(out) if z < 0]
        choices = ([tuple(rng.sample(plus, 2))] if len(plus) >= 2 else [])
        if len(neg) >= 3: choices.append(tuple(rng.sample(neg, 3)))
        if not choices: return st
        inds = choices[rng.randrange(len(choices))]; z0, side = out[inds[0]]
        newz = sum(abs(out[i][0]) for i in inds) * (1 if z0 > 0 else -1)
        for i in sorted(inds, reverse=True): out.pop(i)
        out.append((newz, side))
    else:
        plus = [i for i, (z, _) in enumerate(out) if z >= 2]
        neg = [i for i, (z, _) in enumerate(out) if z <= -3]
        choices = [('p', i) for i in plus] + [('n', i) for i in neg]
        if not choices: return st
        typ, i = choices[rng.randrange(len(choices))]; z, side = out.pop(i); n = abs(z)
        k = 2 if typ == 'p' else 3
        parts = randparts(n, k, rng)
        sg = 1 if typ == 'p' else -1
        out.extend((sg*x, side) for x in parts)
    A = [z for z, s in out if s == 0]; B = [z for z, s in out if s == 1]
    q = cstate(A, B)
    return q if valid(q, need_pairfree) else st

# Independent evaluator: Laurent monomials followed by SU(2) character-basis inversion.
def laurent(word):
    d = {(0, 0): 1}
    for z in word:
        n = abs(z); eps = 1 if z > 0 else -1; q = defaultdict(int)
        for (r, s), v in d.items():
            for a in range(-n, n+1, 2): q[r+a, s] += v
            for b in range(-n, n+1, 2): q[r, s+b] += eps*v
        d = dict(q)
    return d

def laurent_irreps(word):
    L = laurent(word); N = sum(map(abs, word)); out = {}
    for r in range(N, -1, -1):
        for s in range(N, -1, -1):
            v = (L.get((r, s), 0) - L.get((r+2, s), 0)
                 - L.get((r, s+2), 0) + L.get((r+2, s+2), 0))
            if v: out[r, s] = v
    return out

def independent_profile(st):
    A, B = st; X = laurent_irreps(A); Y = laurent_irreps(B)
    assert X == table(A) and Y == table(B)
    layer = defaultdict(int)
    for (r, s), x in X.items():
        v = x*Y.get((r, s), 0)
        if v: layer[r+s] += v
    run = mass = 0; rawmin = 0; rawT = None; best = None
    for t in sorted(layer):
        run += layer[t]; mass += max(layer[t], 0)
        if run < rawmin: rawmin, rawT = run, t
        if mass:
            q = (Fraction(run, mass), t, run, mass)
            if best is None or q[:2] < best[:2]: best = q
    return best, rawmin, rawT, tuple(sorted(layer.items())), X, Y

def verify_independent(st, primary):
    other = independent_profile(st)
    assert primary == other[:4], (st, primary, other[:4])
    return other

F1 = (1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)
tight = (-1,-3,-3,-4,-5,-5,-6,-7)
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
rng = random.Random(181)
seed_states = set()
for w in (F1, tight, *census):
    ww = cword(w)
    if ww == cword(F1) or ww == cword(tight): seed_states.update(all_splits(ww))
    elif len(ww) <= 8: seed_states.update(all_splits(ww))
    else:
        for _ in range(6): seed_states.add(random_split(ww, rng))
for k in range(3, 15):
    tri = k*(k+1)//2
    sigp = 1 if k % 2 == 0 else -1
    for p in range(k+1, min(k+4, 120-tri)+1):
        w = tuple(-i for i in range(1, k+1)) + (sigp*p,)
        if len(w) <= 9: seed_states.update(all_splits(w))
        else:
            for _ in range(5): seed_states.add(random_split(w, rng))
for _ in range(48):
    w = random_word(rng, True)
    for __ in range(4): seed_states.add(random_split(w, rng))
for _ in range(42):
    w = random_word(rng, False)
    for __ in range(4): seed_states.add(random_split(w, rng))
seed_states.update(all_splits((-1,-1,1,1,2)))
for _ in range(30):
    base = random_word(rng, True)
    n, m = rng.sample(range(1, 18), 2)
    w = base + (n, -n, m, -m)
    if len(w) <= 18 and sum(map(abs, w)) <= 120 and sum(z < 0 for z in w) % 2 == 0:
        for __ in range(4): seed_states.add(random_split(w, rng))

stamp(f'seed states={len(seed_states)} F1_splits={len(all_splits(F1))}')

@lru_cache(maxsize=100000)
def eval_state(st):
    return height_profile(st)

failures = set()
def note(st, prof, domain):
    if prof[1] < 0:
        key = (st, prof[1], prof[2])
        if key not in failures:
            alt = verify_independent(st, prof)
            failures.add(key)
            stamp(f'FAIL domain={domain} A={st[0]} B={st[1]} W={sum(map(abs, whole(st)))} T={prof[2]} Pi={prof[1]} ratio={prof[0]} Phi={sum(v for _,v in prof[3])} layers={prof[3]} independent_Pi={alt[1]} independent_layers={alt[3]}')

results = {}
for name, need_pf in (('pair_free', True), ('pairs_allowed', False)):
    states = [st for st in seed_states if valid(st, need_pf)]
    seen = set(); best = None; bins = {}; lows = (0, None)
    for ix, st in enumerate(states):
        p = eval_state(st); note(st, p, name); seen.add(st)
        if p[0] is not None:
            row = (p[0], st, p)
            if best is None or row[0][:2] < best[0][:2]: best = row
            W = sum(map(abs, whole(st))); b = min(5, (W-1)//20)
            if b not in bins or row[0][:2] < bins[b][0][:2]: bins[b] = row
        if p[1] < lows[0]: lows = (p[1], (st, p[2]))
        if (ix+1) % 250 == 0: stamp(f'{name} seed profiles={ix+1}/{len(states)}')
    if not states: raise RuntimeError('no valid seed states')
    elite = sorted((eval_state(st)[0], st) for st in states if eval_state(st)[0] is not None)[:40]
    for restart in range(12):
        cur = elite[restart % len(elite)][1] if restart < len(elite) else elite[rng.randrange(len(elite))][1]
        cv = eval_state(cur)[0]
        for step in range(12):
            proposals = []
            for _ in range(3):
                nw = mutate(cur, rng, need_pf)
                pp = eval_state(nw); note(nw, pp, name); seen.add(nw)
                if pp[0] is not None: proposals.append((pp[0], nw, pp))
                if pp[1] < lows[0]: lows = (pp[1], (nw, pp[2]))
            if not proposals: continue
            nv, nw, np = min(proposals, key=lambda x: (x[0][0], x[0][1]))
            if best is None or nv[:2] < best[0][:2]: best = (nv, nw, np)
            W = sum(map(abs, whole(nw))); b = min(5, (W-1)//20)
            if b not in bins or nv[:2] < bins[b][0][:2]: bins[b] = (nv, nw, np)
            if cv is None or nv[:2] <= cv[:2] or rng.randrange(5) == 0:
                cur, cv = nw, nv
        stamp(f'{name} hill restart={restart+1}/12 visited={len(seen)} best={None if best is None else best[0]}')
    bins = {}
    for st in seen:
        p = eval_state(st)
        if p[0] is None: continue
        W = sum(map(abs, whole(st))); b = min(5, (W-1)//20)
        row = (p[0], st, p)
        if b not in bins or row[0][:2] < bins[b][0][:2]: bins[b] = row
    best = min(((eval_state(st)[0], st, eval_state(st)) for st in seen
                if eval_state(st)[0] is not None), key=lambda x:(x[0][0],x[0][1]))
    results[name] = (best, bins, len(seen), lows)
    stamp(f'{name} visited_total={len(seen)} minimum={best[0]} lowest_raw_prefix={lows[0]}')

f1rows = []
for st in all_splits(F1):
    p = eval_state(st); note(st, p, 'F1_all_splits')
    f1rows.append((p[0], st, p))
f1best = min(f1rows, key=lambda x:(x[0][0],x[0][1]))
f1tight = cstate((-4,-4,-4,-2), tuple(z for z in F1 if z not in (-4,-4,-4,-2)))
f1tp = eval_state(f1tight)
assert len(f1rows) == 599
assert f1best[1] == f1tight
assert f1best[0] == (Fraction(2756159,19046176), 12, 5512318, 38092352)
assert f1tp[1] == 0
verify_independent(f1tight, f1tp)
stamp(f'F1 all_splits={len(f1rows)} minimum={f1best[0]} A={f1best[1][0]} B={f1best[1][1]}')

for name, (best, bins, count, lows) in results.items():
    st = best[1]; p = best[2]; verify_independent(st, p)
    pairfactor = not pairfree(whole(st))
    print('EXTREME', name, 'score', best[0], 'W', sum(map(abs, whole(st))),
          'A', st[0], 'B', st[1], 'pair_bearing', pairfactor,
          'raw_min', p[1], 'raw_T', p[2], 'layers', p[3])
    print('W_BINS', name)
    for b, row in sorted(bins.items()):
        print('BIN', b, 'score', row[0], 'W', sum(map(abs, whole(row[1]))),
              'A', row[1][0], 'B', row[1][1], 'raw_min', row[2][1])
    print('SEARCH_COUNT', name, count, 'LOWEST_RAW_PREFIX', lows)

witnesses = [
    (cstate((1,-2,-4,-4,-4), (1,)+(3,)*6+(5,)*4+(6,)), -24695910,
     ((1,206613818),(3,863409696),(5,452240438),(7,-552023508),
      (9,-663275890),(11,-331660464),(13,70846556),(15,74390904))),
    (cstate((1,-2,-4,-4,-4), (1,)+(3,)*5+(5,)*5+(6,)), -26969676,
     ((1,255355504),(3,1071276562),(5,574413842),(7,-674538760),
      (9,-831357746),(11,-422119078),(13,87231344),(15,95065210))),
    (cstate((1,-2,-4,-4,-4), (1,1,2)+(3,)*5+(4,4)+(5,5)+(8,)), -61656134,
     ((1,488785802),(3,2056983050),(5,1076172248),(7,-1301389998),
      (9,-1568859336),(11,-813347900),(13,112226582),(15,142752892)))
]
for st, expected_pi, expected_layers in witnesses:
    p = eval_state(st); alt = verify_independent(st, p)
    assert p[1] == expected_pi and p[2] == 11 and p[3] == expected_layers
    print('WITNESS', 'A', st[0], 'B', st[1], 'W', sum(map(abs,whole(st))),
          'factors', len(whole(st)), 'pair_free', pairfree(whole(st)),
          'score', p[0], 'Pi11', p[1], 'Phi', sum(v for _,v in p[3]),
          'layers', p[3], 'independent_Pi11', alt[1],
          'independent_layer_match', alt[3] == p[3])
expected_keys = {(st, pi, 11) for st, pi, _ in witnesses}
assert {(st, pi, T) for st, pi, T in failures} == expected_keys
assert results['pair_free'][2] == 4274 and results['pairs_allowed'][2] == 4571
print('PAIR_SCOPE_TEST', eval_state(cstate((-1,-1),(1,1,2))))
print('FAILURE_COUNT', len(failures), 'EXPECTED_FAILURES_REPRODUCED')
print('PASS: exact search rerun and all failure profiles cross-checked')
