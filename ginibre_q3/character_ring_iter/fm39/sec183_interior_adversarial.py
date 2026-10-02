from collections import defaultdict
from fractions import Fraction
from datetime import datetime, timezone
from itertools import combinations

def stamp(s):
    print(datetime.now(timezone.utc).isoformat(timespec='seconds'), s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def dense_mul(F, z, cap=None):
    n = abs(z)
    sg = 1 if z > 0 else -1
    d = len(F)-1
    D = d+n if cap is None else min(d+n, cap)
    G = [[0]*(D+1) for _ in range(D+1)]
    for ax in (0, 1):
        for v in range(min(d, D)+1):
            last = d-v
            ps = [0]*(last+3)
            for u in range(last+1):
                ps[u+2] = ps[u] + (F[u][v] if ax == 0 else F[v][u])
            for t in range(D-v+1):
                lo = abs(t-n)
                hi = min(last, t+n)
                if (hi-lo) % 2:
                    hi -= 1
                if lo <= hi:
                    q = ps[hi+2]-ps[lo]
                    if ax == 0:
                        G[t][v] += q
                    else:
                        G[v][t] += sg*q
    return G

def dense_table(w, cap=None):
    w = tuple(sorted(w, key=lambda z: (-abs(z), z)))
    rem = sum(map(abs, w))
    F = [[1]]
    for z in w:
        rem -= abs(z)
        F = dense_mul(F, z, None if cap is None else rem+cap)
    return F

def sparse_table(w):
    d = {(0, 0): 1}
    for z in sorted(w, key=lambda z: (-abs(z), z)):
        n = abs(z)
        sg = 1 if z > 0 else -1
        q = defaultdict(int)
        for (r, s), v in d.items():
            for t in cg(r, n):
                q[t, s] += v
            for t in cg(s, n):
                q[r, t] += sg*v
        d = {k: v for k, v in q.items() if v}
    return d

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

def layers_from_tables(A, B, u, v):
    cap = sum(map(abs, A))
    X = dense_table(A)
    Y = dense_table(B+(u, v), cap)
    Z = dense_table(B+(-u, -v), cap)
    P, M, L, Q = [], [], [], []
    for h in range(cap+1):
        p = m = ell = 0
        for r in range(h+1):
            s = h-r
            x = X[r][s] if r < len(X) and s < len(X[r]) else 0
            y = Y[r][s] if r < len(Y) and s < len(Y[r]) else 0
            z = Z[r][s] if r < len(Z) and s < len(Z[r]) else 0
            assert (y+z) % 2 == 0 and (y-z) % 2 == 0
            p += x*((y+z)//2)
            m += x*((y-z)//2)
            ell += x*y
        P.append(p)
        M.append(m)
        L.append(ell)
        Q.append(p+min(m, 0))
    return P, M, L, Q

def direct_layers(A, B, u, v):
    cap = sum(map(abs, A))
    X = sparse_table(A)
    Y = sparse_table(B+(u, v))
    Z = sparse_table(B+(-u, -v))
    P, M, L, Q = [], [], [], []
    for h in range(cap+1):
        p = m = ell = 0
        for r in range(h+1):
            s = h-r
            x = X.get((r, s), 0)
            y = Y.get((r, s), 0)
            z = Z.get((r, s), 0)
            assert (y+z) % 2 == 0 and (y-z) % 2 == 0
            p += x*((y+z)//2)
            m += x*((y-z)//2)
            ell += x*y
        P.append(p)
        M.append(m)
        L.append(ell)
        Q.append(p+min(m, 0))
    return P, M, L, Q

def prefix(v):
    out = []
    s = 0
    for x in v:
        s += x
        out.append(s)
    return out

def min_ratio(v, den):
    s = d = 0
    best = None
    for t, (x, y) in enumerate(zip(v, den)):
        s += x
        d += y
        if d:
            q = Fraction(s, d)
            if best is None or (q, t) < (best[0], best[1]):
                best = (q, t, s, d)
    return best

def top_pair(w):
    pairs = [(i, j) for i in range(len(w)-1)
             for j in range(i+1, len(w)-1)
             if (abs(w[i])+abs(w[j])) % 2 == 0]
    return max(pairs, key=lambda ij:
               (abs(w[ij[0]])+abs(w[ij[1]]),
                max(abs(w[ij[0]]), abs(w[ij[1]]))))

def verify_top(w, labpair, plain_expected=None, budget_expected=None):
    i, j = top_pair(w)
    assert tuple(sorted((w[i], w[j]))) == tuple(sorted(labpair))
    C = tuple(z for k, z in enumerate(w) if k not in (i, j))
    A, B = cut(C)
    P, M, L, Q = layers_from_tables(A, B, w[i], w[j])
    if plain_expected is not None:
        q = min_ratio(L, [max(x, 0) for x in L])
        assert (q[0], q[1]) == plain_expected, (q, plain_expected)
    if budget_expected is not None:
        q = min_ratio(Q, [max(x, 0) for x in P])
        assert (q[0], q[1]) == budget_expected, (q, budget_expected)
    return A, B, P, M, L, Q

w = (-1,)*6 + (-2,)*3 + (-3,)*3 + (-4, -7)
A, B, P, M, L, Q = verify_top(
    w, (-2, -4),
    (Fraction(32465, 35756), 13),
    (Fraction(-10145, 57268), 13)
)
assert (P, M, L, Q) == direct_layers(A, B, -2, -4)
assert (A, B) == ((-7, -2, -2, -1, -1),
                   (-3, -3, -3, -1, -1, -1, -1))
assert sum(L) == 64930 and sparse_table(w).get((0, 0), 0) == 64930
assert prefix(Q)[13] == -20290
assert sum(max(x, 0) for x in P[:14]) == 114536
assert P == [0,4488,0,50864,0,59184,0,-21852,0,-18374,0,-4304,0,-1248]
assert M == [0,-3200,0,-40100,0,-44420,0,47692,0,37230,0,298,0,-1328]
assert prefix(L) == [0,1288,1288,12052,12052,26816,26816,52656,
                     52656,71512,71512,67506,67506,64930]
assert prefix(Q) == [0,1288,1288,12052,12052,26816,26816,4964,
                     4964,-13410,-13410,-17714,-17714,-20290]
stamp('Budget counterexample verified by dense and sparse exact evaluators')
print('A=', A, 'B=', B)
print('P=', P)
print('M=', M)
print('L=', L)
print('Q=', Q)
print('Pi prefixes=', prefix(L))
print('budget prefixes=', prefix(Q))

W1 = (1,-2,-4,-4,-4,1)+(3,)*6+(5,)*4+(6,)
W2 = (1,-2,-4,-4,-4,1)+(3,)*5+(5,)*5+(6,)
W3 = (1,-2,-4,-4,-4,1,1,2)+(3,)*5+(4,4)+(5,5)+(8,)
for name, ww, expected in [
    ('W1', W1, {(-4,-2)}),
    ('W2', W2, set()),
    ('W3', W3, set())
]:
    plainbad = []
    budbad = set()
    for i, j in combinations(range(len(ww)), 2):
        C = tuple(z for k, z in enumerate(ww) if k not in (i, j))
        aa, bb = cut(C)
        p, m, ell, q = layers_from_tables(aa, bb, ww[i], ww[j])
        if min(prefix(ell), default=0) < 0:
            plainbad.append((ww[i], ww[j]))
        if min(prefix(q), default=0) < 0:
            budbad.add(tuple(sorted((ww[i], ww[j]))))
    assert not plainbad and budbad == expected, (name, plainbad, budbad)
    stamp(f'{name} all-pair audit: factors={len(ww)} '
          f'W={sum(map(abs, ww))} plain_failures=0 '
          f'budget_bad_classes={sorted(budbad)}')

records = [
    ('PF plain overall',
     (1,1,1,1,2,3,3,3,3,-4,-4,5,5,5,5,6), (5,5), 'top',
     (Fraction(20159257,33537232),21), None),
    ('pairs plain overall',
     (-1,-1,1,1,1,-2,2,-3,3,3,3,3,4,4,4,5,5,-7,-7,8), (-7,-7), 'top',
     (Fraction(361867715,369393743),27), None),
    ('pairs budget overall',
     (1,1,1,-2,2,3,3,3,3,3,-4,4,4,5,5,12,8), (-4,12), 'top',
     None, (Fraction(50171483,122034656),24)),
    ('PF best-pair sample',
     (-4,-4,8,-13,14,15,15,-19,-19,28,-1), (-4,8), 'all',
     (Fraction(1),None), None),
    ('pairs best-pair sample',
     (-10,10,11,-23,20), (-10,20), 'all',
     (Fraction(1),None), None),
]
for name, ww, pair, kind, pexp, bexp in records:
    if kind == 'top':
        i, j = top_pair(ww)
    else:
        i, j = next((i, j) for i, j in combinations(range(len(ww)), 2)
                    if tuple(sorted((ww[i], ww[j]))) == tuple(sorted(pair)))
    assert tuple(sorted((ww[i], ww[j]))) == tuple(sorted(pair))
    C = tuple(z for k, z in enumerate(ww) if k not in (i, j))
    aa, bb = cut(C)
    P0, M0, L0, Q0 = layers_from_tables(aa, bb, ww[i], ww[j])
    if pexp is not None:
        q = min_ratio(L0, [max(x, 0) for x in L0])
        assert q[0] == pexp[0] and (pexp[1] is None or q[1] == pexp[1]), (name, q)
    if bexp is not None:
        q = min_ratio(Q0, [max(x, 0) for x in P0])
        assert q[0] == bexp[0] and q[1] == bexp[1], (name, q)
    stamp(f'EXTREME CHECK {name}: W={sum(map(abs, ww))} '
          f'pair={ww[i]},{ww[j]} '
          f'plain={None if pexp is None else min_ratio(L0,[max(x,0) for x in L0])} '
          f'budget={None if bexp is None else min_ratio(Q0,[max(x,0) for x in P0])}')
stamp('PASS: exact recurrence, independent sparse evaluator, witness audits, '
      'and reported global profiles agree')
