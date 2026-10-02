from collections import defaultdict
from datetime import datetime, timezone

def stamp(message):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
          message, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def side(word):
    # Multiplicities split by parity of negative factors sent to y.
    d = [{(0, 0): 1}, {}]
    for z in word:
        n = abs(z)
        nxt = [defaultdict(int), defaultdict(int)]
        for parity, row in enumerate(d):
            for (r, s), value in row.items():
                for rr in cg(r, n):
                    nxt[parity][rr, s] += value
                yp = parity ^ int(z < 0)
                for ss in cg(s, n):
                    nxt[yp][r, ss] += value
        d = [dict(row) for row in nxt]
    even, odd = d
    f = {key: even.get(key, 0) - odd.get(key, 0)
         for key in even.keys() | odd.keys()}
    return {key: value for key, value in f.items() if value}

def signed_table(word):
    d = {(0, 0): 1}
    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        nxt = defaultdict(int)
        for (r, s), value in d.items():
            for rr in cg(r, n):
                nxt[rr, s] += value
            for ss in cg(s, n):
                nxt[r, ss] += eps * value
        d = {key: value for key, value in nxt.items() if value}
    return d

def pair_indices(word, labels):
    indices = []
    for label in labels:
        indices.append(next(i for i, z in enumerate(word)
                            if z == label and i not in indices))
    return tuple(indices)

def interior(word, pair):
    rest = sorted((i for i in range(len(word)) if i not in pair),
                  key=lambda i: (-abs(word[i]), word[i], i))
    A, B = [], []
    wa = wb = 0
    for i in rest:
        if wa <= wb:
            A.append(i)
            wa += abs(word[i])
        else:
            B.append(i)
            wb += abs(word[i])
    if wa > wb:
        A, B = B, A
    return A, B + list(pair)

def test_pair(word, labels):
    pair = pair_indices(word, labels)
    A, B = interior(word, pair)
    fa = side(tuple(word[i] for i in A))
    fb = side(tuple(word[i] for i in B))
    layer = defaultdict(int)
    even_h = defaultdict(int)
    odd_h = defaultdict(int)
    for rs in fa.keys() | fb.keys():
        t = sum(rs)
        product = fa.get(rs, 0) * fb.get(rs, 0)
        layer[t] += product
        if product >= 0:
            even_h[t] += product
        else:
            odd_h[t] -= product
    layer = {t: v for t, v in layer.items() if v}
    even_h = {t: v for t, v in even_h.items() if v}
    odd_h = {t: v for t, v in odd_h.items() if v}
    heights = layer.keys() | even_h.keys() | odd_h.keys()
    assert all(layer.get(t, 0) == even_h.get(t, 0) - odd_h.get(t, 0)
               for t in heights)
    assert all(v >= 0 for v in layer.values())
    assert all(even_h.get(t, 0) >= odd_h.get(t, 0) for t in heights)
    assert sum(layer.values()) == signed_table(word).get((0, 0), 0)
    # Independent generic entries give rank O_t whenever E_t >= O_t.
    assert all(min(even_h.get(t, 0), odd_h.get(t, 0)) == odd_h.get(t, 0)
               for t in heights)
    return (tuple(word[i] for i in pair), sum(even_h.values()),
            sum(odd_h.values()), tuple(sorted(odd_h.items())))

corrected = [
    ((-1,-1,-1,-2,-2,-2,3), (-1,-1)),
    ((1,1,1,-2,-2,-2,-3), (1,1)),
    ((-1,-2,-2,-2,3,3,3), (-1,-2)),
    ((1,-2,-2,-2,-3,-3,-3), (1,-2)),
    ((1,1,1,1,1,-2,-2,3), (1,1)),
    ((-1,-1,-1,-1,-1,-2,-2,-3), (-1,-1)),
    ((1,1,1,-2,-2,3,3,3), (1,1)),
    ((-1,-1,-1,-2,-2,-3,-3,-3), (-1,-1)),
    ((1,-2,-2,3,3,3,3,3), (1,-2)),
    ((-1,-2,-2,-3,-3,-3,-3,-3), (-1,-2)),
]
sec181 = [
    ((1,1,-2)+(3,)*6+(-4,)*3+(5,)*4+(6,), (1,3)),
    ((1,1,-2)+(3,)*4+(-4,)*3+(5,)*5+(6,), (1,1)),
    ((1,1,1,-2,2)+(3,)*5+(-4,)*3+(4,)*2+(5,)*2+(8,), (1,3)),
]
runs = [(tuple(-i for i in range(1,k+1)) +
         ((-1 if k % 2 else 1)*(k+1),), (-1,-2))
        for k in range(1,10)]
f1 = (((1,1,-2)+(3,)*4+(-4,)*3+(5,)*4+(8,), (1,-2)))
checks = [
    ('STR8b-corrected-10', corrected,
     [(80,0),(80,0),(158,0),(158,0),(118,22),(118,22),
      (152,34),(152,34),(408,0),(408,0)]),
    ('SEC181-three-witnesses', sec181,
     [(150956716,30415166),(0,0),(227384652,34061312)]),
    ('all-minus-runs-k1..9', runs,
     [(0,0),(2,0),(4,0),(0,0),(0,0),(186,2),(1108,128),(0,0),(0,0)]),
    ('F1', [f1], [(9460028,1309286)]),
]
for name, rows, expected in checks:
    stamp('START '+name+' cases='+str(len(rows)))
    got = []
    for i, (word, labels) in enumerate(rows):
        result = test_pair(word, labels)
        got.append(result[1:3])
        print(name, i, 'pair=', result[0], 'H_even(d0)=', result[1],
              'H_odd(d0)=', result[2], 'odd_by_height=', result[3],
              flush=True)
    assert got == expected, (name, got, expected)
    stamp('PASS '+name)
stamp('PASS exact layer Euler checks and generic full-column-rank capacity checks')
