import gzip, re, pathlib, datetime
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction

def stamp(msg):
    print(datetime.datetime.now(datetime.timezone.utc).strftime(
        '%Y-%m-%d %H:%M:%S UTC'), msg, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

# Direct evaluator in the U_r(x) U_s(y) character basis.
@lru_cache(None)
def char_table(word):
    d = {(0, 0): 1}
    for z in sorted(word, key=lambda z: (-abs(z), z)):
        n, eps = abs(z), (1 if z > 0 else -1)
        q = defaultdict(int)
        for (r, s), v in d.items():
            for x in cg(r, n):
                q[x, s] += v
            for y in cg(s, n):
                q[r, y] += eps*v
        d = {k: v for k, v in q.items() if v}
    return d

def canonical(word):
    return tuple(sorted(word, key=lambda z: (-abs(z), z)))

def interior_cut(word):
    a, b, wa, wb = [], [], 0, 0
    for z in canonical(word):
        if wa <= wb:
            a.append(z)
            wa += abs(z)
        else:
            b.append(z)
            wb += abs(z)
    if wa > wb:
        a, b = b, a
    return canonical(a), canonical(b)

def layer_from_tables(x, y):
    q = defaultdict(int)
    for (r, s), v in x.items():
        if (r, s) in y:
            q[r+s] += v*y[r, s]
    return tuple(sorted(q.items()))

@lru_cache(None)
def direct_profile(a, b):
    return layer_from_tables(char_table(a), char_table(b))

def pair_sides(word, i, j):
    c = tuple(z for k, z in enumerate(word) if k not in (i, j))
    a, b = interior_cut(c)
    return a, canonical(b + (word[i], word[j]))

# Independent evaluator: Laurent expansion followed by highest-weight
# inversion in both coordinates.
@lru_cache(None)
def laurent_table(word):
    d = {(0, 0): 1}
    for z in word:
        n, eps = abs(z), (1 if z > 0 else -1)
        q = defaultdict(int)
        for (r, s), v in d.items():
            for x in range(-n, n+1, 2):
                q[r+x, s] += v
            for y in range(-n, n+1, 2):
                q[r, s+y] += eps*v
        d = dict(q)
    return d

@lru_cache(None)
def laurent_chars(word):
    d, N, out = laurent_table(word), sum(map(abs, word)), {}
    for r in range(N, -1, -1):
        for s in range(N, -1, -1):
            v = (d.get((r, s), 0) - d.get((r+2, s), 0)
                 - d.get((r, s+2), 0) + d.get((r+2, s+2), 0))
            if v:
                out[r, s] = v
    return out

def independent_profile(a, b):
    return layer_from_tables(laurent_chars(a), laurent_chars(b))

def pair_metrics(word, crosscheck=False):
    word = canonical(word)
    good, zero, bad = [], 0, []
    best = None
    for i in range(len(word)):
        for j in range(i+1, len(word)):
            a, b = pair_sides(word, i, j)
            p = direct_profile(a, b)
            if crosscheck:
                q = independent_profile(a, b)
                assert p == q, (word, i, j, p, q)
            vals = [v for t, v in p]
            mn = min(vals, default=0)
            mass = sum(v for v in vals if v > 0)
            if mn < 0:
                bad.append((i, j, word[i], word[j], p))
            else:
                good.append((i, j, word[i], word[j], mn, mass, p))
                if mass == 0:
                    zero += 1
                ratio = Fraction(mn, mass) if mass else Fraction(0)
                if best is None or ratio > best[0]:
                    best = (ratio, i, j, mn, mass, p)
    return word, good, zero, bad, best

pair_bearing = (17, -14, 14, 14, 13, -11, 11, 11, -10, 9, 7, 4, -3, 2)
pair_free = (21, 17, -16, -16, -13, -13, -7, -7, -6, -6, 4, -3, -3)

for name, word, expected_n, expected_best in (
    ('pair-bearing W=140', pair_bearing, 91,
     Fraction(41830, 6840525217)),
    ('pair-free W=132', pair_free, 78,
     Fraction(86643, 482670154)),
):
    L, good, zero, bad, best = pair_metrics(word, crosscheck=True)
    assert len(good)+len(bad) == expected_n and not zero
    assert len(good) > 0 and best[0] == expected_best
    assert all(v >= 0 for t, v in best[5])
    bad_classes = sorted({tuple(sorted((x[2], x[3]))) for x in bad})
    stamp(f'{name}: W={sum(map(abs,L))} factors={len(L)} '
          f'good={len(good)}/{expected_n} best={best[0]} '
          f'pair={(L[best[1]],L[best[2]])} min={best[3]} mass={best[4]} '
          f'bad-label-classes={bad_classes}')
    print('BEST_LAYERS', name, best[5], flush=True)

# Complete W <= 40 no-flip census, checking every pair and layer.
path = pathlib.Path(
    'ginibre_q3/character_ring_iter/fm39/sec166_census_w40_noflip.log.gz')
rows = []
for line in gzip.decompress(path.read_bytes()).decode().splitlines():
    if line.startswith('NOFLIP'):
        p = int(re.search(r'p=(-?\d+)', line).group(1))
        b = [int(x) for x in re.search(
            r'B=([-\d ]+?)\s+phi', line).group(1).split()]
        rows.append(tuple(b + [p]))
assert len(rows) == 5430
stamp('census loaded: 5430 rows')

no_pair, zero_rows, ranked = [], 0, []
for ix, word in enumerate(rows):
    L = word
    n_good = n_zero = 0
    best = None
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            a, b = pair_sides(L, i, j)
            p = direct_profile(a, b)
            vals = [v for t, v in p]
            mn = min(vals, default=0)
            mass = sum(v for v in vals if v > 0)
            if mn >= 0:
                n_good += 1
                if mass == 0:
                    n_zero += 1
                else:
                    r = Fraction(mn, mass)
                    if best is None or r > best[0]:
                        best = (r, i, j, mn, mass, p)
    if n_good == 0:
        no_pair.append(L)
    if n_zero:
        zero_rows += 1
    ranked.append((best[0] if best else Fraction(0),
                   L, n_good, n_zero, best))
    if (ix+1) % 1000 == 0:
        stamp(f'census progress {ix+1}/5430; '
              f'direct-table cache={char_table.cache_info().currsize}')

assert not no_pair and zero_rows == 0
minimum = min(ranked, key=lambda x: x[0])
expected_word = (-1, -2, -3, -3, -4, -4, -5, -5, -6, -7, 10)
assert minimum[1] == expected_word
assert minimum[0] == Fraction(304, 41763)
assert minimum[2] == 46 and minimum[4][1:5] == (8, 10, 608, 83526)
stamp(f'census complete: no all-pair failures; zero-profile rows={zero_rows}; '
      f'tight row={minimum[1]} best={minimum[0]} pair=(-6,10) '
      f'min={minimum[4][3]} mass={minimum[4][4]}')
print('CENSUS_TIGHT_LAYERS', minimum[4][5], flush=True)
print('PASS: direct census and independent Laurent checks agree.', flush=True)
