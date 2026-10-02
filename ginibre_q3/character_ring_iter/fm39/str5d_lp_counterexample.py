from collections import defaultdict, Counter
from functools import lru_cache
from itertools import product
from datetime import datetime, timezone

def stamp(msg):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC'),
          msg, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def table(word):
    d = {(0, 0): 1}
    for z in sorted(word, key=lambda z: (-abs(z), z)):
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (a, b), v in d.items():
            for c in cg(a, n):
                q[c, b] += v
            for c in cg(b, n):
                q[a, c] += eps*v
        d = {k: v for k, v in q.items() if v}
    return d

def interior_cut(C):
    A, B = [], []
    wa = wb = 0
    for z in sorted(C, key=lambda z: (-abs(z), z)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    key = lambda z: (-abs(z), z)
    return tuple(sorted(A, key=key)), tuple(sorted(B, key=key))

def pair_types(L):
    cnt = Counter(L)
    vals = sorted(cnt, key=lambda z: (abs(z), z))
    out = []
    for i, u in enumerate(vals):
        for v in vals[i:]:
            if u == v and cnt[u] < 2:
                continue
            out.append((u, v))
    return out

def profile(L, u, v):
    C = list(L)
    C.remove(u)
    C.remove(v)
    A, B = interior_cut(C)
    key = lambda z: (-abs(z), z)
    Bp = tuple(sorted(B + (u, v), key=key))
    fa, fb = table(A), table(Bp)
    layers = defaultdict(int)
    for (r, s), x in fa.items():
        y = fb.get((r, s), 0)
        if y:
            layers[r+s] += x*y
    return A, Bp, dict(sorted(layers.items()))

def words(N):
    for m1 in range(N+1):
        for m2 in range(N-m1+1):
            m3 = N-m1-m2
            mult = (m1, m2, m3)
            present = [i+1 for i, m in enumerate(mult) if m]
            for signs in product((1, -1), repeat=len(present)):
                smap = dict(zip(present, signs))
                if sum(m for i, m in enumerate(mult, 1)
                       if m and smap[i] < 0) % 2:
                    continue
                L = tuple(sorted(
                    (smap[i+1]*(i+1)
                     for i, m in enumerate(mult)
                     for _ in range(m)),
                    key=lambda z: (-abs(z), z)))
                yield L

def laurent(word):
    d = {(0, 0): 1}
    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        terms = [((n-2*j, 0), 1) for j in range(n+1)]
        terms += [((0, n-2*j), eps) for j in range(n+1)]
        q = defaultdict(int)
        for (a, b), v in d.items():
            for (c, e), w in terms:
                q[a+c, b+e] += v*w
        d = {k: v for k, v in q.items() if v}
    return d

@lru_cache(None)
def irreps_laurent(word):
    q = laurent(word)
    W = sum(map(abs, word))
    out = {}
    for r in range(W+1):
        for s in range(W+1):
            value = (q.get((r, s), 0) - q.get((r+2, s), 0)
                     - q.get((r, s+2), 0) + q.get((r+2, s+2), 0))
            if value:
                out[r, s] = value
    return out

nwords = ncuts = 0
bad = []
for N in range(2, 21):
    here = cuts = 0
    bad_here = []
    for L in words(N):
        here += 1
        nwords += 1
        records = []
        for u, v in pair_types(L):
            cuts += 1
            ncuts += 1
            records.append((u, v, *profile(L, u, v)))
        if not any(all(x >= 0 for x in rec[4].values())
                   for rec in records):
            bad.append((N, L, records))
            bad_here.append(L)
    stamp(f'N={N} lists={here} pair-class cuts={cuts} '
          f'LP failures={len(bad_here)}')

stamp(f'TOTAL lists={nwords} pair-class cuts={ncuts} failures={len(bad)}')
expected_bad = {
    tuple([-3]*14 + [2]*3 + [1]*2),
    tuple([3]*14 + [2]*3 + [-1]*2),
    tuple([-3]*10 + [2]*7 + [1]*2),
    tuple([3]*10 + [2]*7 + [-1]*2),
}
assert (nwords, ncuts, len(bad)) == (6537, 32277, 4)
assert {L for _, L, _ in bad} == expected_bad
assert all(N == 19 for N, _, _ in bad)

first = tuple([-3]*14 + [2]*3 + [1]*2)
assert bad[0][1] == first
odd = {
    1: 75653592, 3: 229970792, 5: 334590060, 7: 189230932,
    9: 420087042, 11: 178076150, 13: 25842642, 15: 57436610,
    17: 11653672, 19: -661444, 21: 1045396, 23: 196284,
}
even = {
    0: 66733568, 2: 145606664, 4: 133708424, 6: 441442388,
    8: 341014342, 10: 115025058, 12: 199032588, 14: 62534496,
    16: 4928552, 18: 10719802, 20: 2396306, 22: -20460,
}

for N, L, records in bad:
    phis = []
    for u, v, A, Bp, lay in records:
        assert any(x < 0 for x in lay.values())
        assert table(A) == irreps_laurent(A)
        assert table(Bp) == irreps_laurent(Bp)
        phis.append(sum(lay.values()))
    assert len(set(phis)) == 1
    stamp(f'BAD LIST {L}: pair classes={len(records)}, '
          f'Phi={phis[0]}, Laurent cross-check=PASS')
    if L == first:
        assert len(records) == 6
        assert all(lay == (odd if (u, v) in {
            (1, 1), (1, 2), (1, -3), (2, 2)
        } else even) for u, v, A, Bp, lay in records)
        assert phis[0] == 1523121728
        for u, v, A, Bp, lay in records:
            neg = tuple((t, z) for t, z in lay.items() if z < 0)
            stamp(f'  pair=({u},{v}) A={A} Bprime={Bp} '
                  f'negative={neg} profile={tuple(lay.items())}')

stamp('EXACT LP COUNTEREXAMPLE AND FINITE CENSUS VERIFIED')
