from collections import defaultdict
from fractions import Fraction
from datetime import datetime, timezone
from itertools import combinations

def stamp(s):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def table(w):
    d = {(0, 0): 1}
    for z in sorted(w, key=lambda z: (-abs(z), z)):
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), c in d.items():
            for t in cg(r, n):
                q[t, s] += c
            for t in cg(s, n):
                q[r, t] += eps*c
        d = {k: v for k, v in q.items() if v}
    return d

def cut(C, mode):
    A, B = [], []
    wa = wb = 0
    key = (lambda z: abs(z)) if mode == 'incA' else (lambda z: -abs(z))
    for z in sorted(C, key=key):
        if wa > wb or (wa == wb and mode == 'tieB'):
            B.append(z)
            wb += abs(z)
        else:
            A.append(z)
            wa += abs(z)
    if wa > wb:
        A, B = B, A
    return tuple(A), tuple(B)

def layer_profile(w, i, j, mode):
    C = tuple(z for k, z in enumerate(w) if k not in (i, j))
    A, B = cut(C, mode)
    cap = sum(map(abs, A))
    X = table(A)
    Y = table(B+(w[i], w[j]))
    L = [
        sum(X.get((r, h-r), 0)*Y.get((r, h-r), 0)
            for r in range(h+1))
        for h in range(cap+1)
    ]
    run = mass = raw = 0
    rawt = None
    best = None
    for h, x in enumerate(L):
        run += x
        mass += max(x, 0)
        if run < raw:
            raw, rawt = run, h
        if mass:
            q = Fraction(run, mass)
            if best is None or (q, h) < best[:2]:
                best = (q, h, run, mass)
    return A, B, tuple(L), best, raw, rawt

def top_pair(w):
    pairs = [
        (i, j)
        for i in range(len(w)-1)
        for j in range(i+1, len(w)-1)
        if (abs(w[i])+abs(w[j])) % 2 == 0
    ]
    return max(
        pairs,
        key=lambda ij: (
            abs(w[ij[0]])+abs(w[ij[1]]),
            max(abs(w[ij[0]]), abs(w[ij[1]]))
        )
    )

def verify_top(name, w, mode, pair, expected):
    i, j = top_pair(w)
    assert tuple(sorted((w[i], w[j]))) == tuple(sorted(pair))
    A, B, L, best, raw, rawt = layer_profile(w, i, j, mode)
    assert (best[0], best[1], best[2], best[3]) == expected
    assert raw >= 0
    stamp(f'TOP {name} mode={mode} W={sum(map(abs,w))} '
          f'factors={len(w)} pair={w[i]},{w[j]} ratio={best[0]} '
          f'T={best[1]} prefix={best[2]} mass={best[3]} rawmin={raw}')

def verify_route(name, w, mode, pair):
    i, j = next(
        (i, j) for i, j in combinations(range(len(w)), 2)
        if tuple(sorted((w[i], w[j]))) == tuple(sorted(pair))
    )
    _, _, _, best, raw, rawt = layer_profile(w, i, j, mode)
    assert best[0] == 1 and raw >= 0
    stamp(f'BESTPAIR {name} mode={mode} W={sum(map(abs,w))} '
          f'pair={w[i]},{w[j]} ratio={best[0]} T={best[1]} rawmin={raw}')

verify_top(
    'default-pairfree',
    (-1,)*6+(-2,)*3+(-3,)*3+(-6,-7),
    'decA', (-2,-6),
    (Fraction(36598,41141), 13, 73196, 82282)
)
verify_top(
    'default-pairs-allowed',
    (1,1,1,-2,2,3,3,3,3,3,-4,-4,4,4,4,5,5,-8),
    'decA', (5,5),
    (Fraction(172952362,176691733), 25, 345904724, 353383466)
)
verify_top(
    'tieB-pairfree',
    (-1,)*6+(-2,)*3+(-3,)*3+(-4,-7),
    'tieB', (-2,-4),
    (Fraction(2941,4392), 11, 64702, 96624)
)
verify_top(
    'tieB-pairs-allowed',
    (1,1,1,-2,2,3,3,3,3,3,-4,4,4,5,5,8),
    'tieB', (5,5),
    (Fraction(8000385,11348714), 21, 16000770, 22697428)
)
verify_top(
    'increasing-pairfree',
    (-7,-15,-20,-26,62,10),
    'incA', (-26,62),
    (Fraction(1), 8, 6, 6)
)
verify_top(
    'increasing-pairs-allowed',
    (-23,24,-28,28,-35,-22),
    'incA', (-23,-35),
    (Fraction(1), 6, 294, 294)
)

oldpf = (1,1,1,1,2,3,3,3,3,-4,-4,5,5,5,5,6)
oldpairs = (-1,-1,1,1,1,-2,2,-3,3,3,3,3,4,4,4,5,5,-7,-7,8)

for mode, expected in [
    ('decA', (Fraction(20159257,33537232),21,40318514,67074464)),
    ('tieB', (Fraction(20159257,20869940),21,40318514,41739880)),
    ('incA', (Fraction(1),1,116880,116880))
]:
    i, j = top_pair(oldpf)
    q = layer_profile(oldpf, i, j, mode)
    assert q[3][:4] == expected and q[4] >= 0
    stamp(f'RETAINED old pairfree mode={mode} ratio={q[3][0]} '
          f'T={q[3][1]} rawmin={q[4]}')

for mode, expected in [
    ('decA', (Fraction(361867715,369393743),27,723735430,738787486)),
    ('tieB', (Fraction(361867715,361921374),27,723735430,723842748)),
    ('incA', (Fraction(1),2,45544796,45544796))
]:
    i, j = top_pair(oldpairs)
    q = layer_profile(oldpairs, i, j, mode)
    assert q[3][:4] == expected and q[4] >= 0
    stamp(f'RETAINED old pairs-allowed mode={mode} ratio={q[3][0]} '
          f'T={q[3][1]} rawmin={q[4]}')

pfbest = (2,2,8,8,-9,-9,-13,-13,14,-17,21,-26,14)
pairbest = (-7,-7,-21,-23,24,-28,28,-22)
incpairbest = (-23,24,-28,28,-35,-22)

verify_route('pairfree', pfbest, 'decA', (2,2))
verify_route('pairfree', pfbest, 'tieB', (2,2))
verify_route('pairfree', pfbest, 'incA', (-9,21))
verify_route('pairs-allowed', pairbest, 'decA', (-21,-28))
verify_route('pairs-allowed', pairbest, 'tieB', (-21,-28))
verify_route('pairs-allowed', incpairbest, 'incA', (-28,28))

stamp('PASS exact sparse verifier: all listed plain prefixes are nonnegative')
