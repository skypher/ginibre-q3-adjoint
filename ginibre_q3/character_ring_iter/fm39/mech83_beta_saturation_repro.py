from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
from fractions import Fraction as Q
from random import Random

@lru_cache(None)
def fusion(ns):
    if not ns:
        return {0: 1}
    out = {}
    n = ns[-1]
    for j, v in fusion(ns[:-1]).items():
        for k in range(abs(j-n), j+n+1, 2):
            out[k] = out.get(k, 0) + v
    return out

@lru_cache(None)
def inv(ns):
    if not ns:
        return 1
    if sum(ns) % 2 or 2*max(ns) > sum(ns):
        return 0
    l = len(ns)
    if l == 1:
        return int(ns[0] == 0)
    d = sum(ns)//2
    ans = 0
    for mask in range(1 << l):
        q = d - sum(ns[i]+1 for i in range(l) if mask >> i & 1)
        if q >= 0:
            ans += (-1)**mask.bit_count() * comb(q+l-2, l-2)
    return ans

def cat(j):
    return comb(2*j, j)//(j+1)

def uchar(n):
    return {n-2*j: (-1)**j*comb(n-j, j)
            for j in range(n//2+1)}

def epoly(word, cap=None):
    if cap is None:
        cap = 2*sum(n for n, e in word)
    P = {(0, 0): 1}
    for n, e in word:
        G = {(2*j, 0): 1 for j in range(min(n, cap//2)+1)}
        if n <= cap:
            for j, c in uchar(n).items():
                G[n, j] = G.get((n, j), 0) + e*c
        T = {}
        for (i, j), v in P.items():
            for (k, l), w in G.items():
                if i+k <= cap:
                    T[i+k, j+l] = T.get((i+k, j+l), 0) + v*w
        P = {ij: v for ij, v in T.items() if v}
    out = [0]*(cap//2+1)
    for (i, j), v in P.items():
        if j % 2 == 0:
            assert i % 2 == 0
            out[i//2] += v*cat(j//2)
    return out

def direct(word):
    ns = tuple(n for n, e in word)
    L = len(ns)
    full = (1 << L)-1
    m = [inv(tuple(sorted(ns[i] for i in range(L) if s >> i & 1)))
         for s in range(1 << L)]
    minus = sum((1 << i) for i, (n, e) in enumerate(word) if e < 0)
    return sum((-1)**((s & minus).bit_count())*m[s]*m[full ^ s]
               for s in range(1 << L))

checks = 0
for L in range(2, 8):
    for ns in combinations_with_replacement(range(1, 11), L):
        a = fusion(ns)
        v = a.get(0, 0)
        if v:
            assert all(a.get(2*j, 0) >= v for j in range(ns[-1]+1))
            checks += 1
assert checks == 9640
print("fusion-channel lists:", checks)

for ns, expected in [
    (tuple(range(64, 72)), 998871442),
    (tuple(range(128, 137)), 3116813934958),
    ((64,)*8, 772636800),
]:
    L = len(ns)
    full = (1 << L)-1
    m = [inv(tuple(ns[i] for i in range(L) if s >> i & 1))
         for s in range(1 << L)]
    f = [m[s]*m[full ^ s] for s in range(1 << L)]
    vals = f[:]
    for j in range(L):
        for s in range(1 << L):
            if not s >> j & 1:
                t = s | (1 << j)
                vals[s], vals[t] = vals[s]+vals[t], vals[s]-vals[t]
    alpha = sum((Q(1, ns[s.bit_length()-1]+1)
                 for s in range(1, 1 << (L-1)) if f[s]), Q(0))
    assert alpha <= Q(2**(L-1)-L-1, ns[0]+1)
    top = Q(0)
    for neg, value in enumerate(vals):
        if neg.bit_count() % 2:
            assert value == 0
            continue
        assert abs(value-2*m[full]) <= 2*m[full]*alpha
        beta = Q(0)
        for s in range(1, 1 << (L-1)):
            if (s & neg).bit_count() % 2 and f[s]:
                denominator = ns[s.bit_length()-1]+1
                assert m[full] >= denominator*f[s]
                beta += Q(1, denominator)
        assert beta <= 1
        assert value >= 2*m[full]*(1-beta)
        top = max(top, beta)
    minimum = min(v for s, v in enumerate(vals) if s.bit_count() % 2 == 0)
    assert minimum == expected
    print("all signs:", ns, "minimum", minimum, "max beta", top)

rng = Random(8302)
moments = chambers = 0
for case in range(240):
    bg = [(rng.randrange(1, 7), rng.choice((-1, 1)))
          for _ in range(rng.randrange(0, 8))]
    P = epoly(bg)
    D = len(P)-1
    assert P == P[::-1] and sum(P) > 0
    for j in range(7):
        assert sum(c*Q(2*i-D, 2)**(2*j) for i, c in enumerate(P)) >= 0
        moments += 1
    for h in range(2, 9):
        d = max(D-h+2, (D+1)//2-1, 0)
        value = sum(c*comb(d-i+h-2, h-2)
                    for i, c in enumerate(P) if i <= d)
        lower = Q(sum(P))
        for j in range(1, h-1):
            lower *= Q(2*d-D+2*j, 2*j)
        assert value >= lower >= 0 and value > 0
        chambers += 1
assert moments == chambers == 1680
print("central moments and chamber bounds:", moments, chambers)

examples = [
    ([(1,1)]*2+[(1,-1)]*2+[(3,-1),(4,1)], 2, 11, 204),
    ([(3,1),(3,-1),(4,-1)], 3, 9, 380),
    ([(4,-1)]*3, 6, 8, 6438),
]
for bg, h, d, expected in examples:
    P = epoly(bg)
    D = len(P)-1
    v = sum(c*comb(d-i+h-2, h-2)
            for i, c in enumerate(P) if i <= d)
    rest = bg + [(d+1, 1)]*h
    p = D+h*(d+1)-2*d
    eps = (-1)**sum(e < 0 for n, e in rest)
    G = epoly(rest, 2*d)
    assert G[d]-G[d-1] == v == expected
    assert direct(rest+[(p,eps)]) == 2*v
print("three chamber values:", [2*x[3] for x in examples])

def escape(r):
    n, p = r+1, 6*r+6
    ns = (3,)*(2*r)+(4,)+(n,n,p)
    origin = fusion(ns).get(0, 0)
    bad = good = 0
    for a in range(r+1):
        for b in range(r+1):
            for c in range(2):
                for d in range(3):
                    if a+b+c+d == 0:
                        continue
                    S = (3,)*(a+b)+(4,)*c+(n,)*d
                    T = (3,)*(2*r-a-b)+(4,)*(1-c)+(n,)*(2-d)+(p,)
                    v = (comb(r,a)*comb(r,b)*comb(2,d) *
                         fusion(S).get(0,0)*fusion(T).get(0,0))
                    if b % 2:
                        bad += v
                    else:
                        good += v
    rest = [(3,1)]*r+[(3,-1)]*r+[(4,1)]+[(n,1)]*2
    G = epoly(rest, 2*r)
    assert origin+good-bad == G[r]-G[r-1] > 0
    assert Q(r*r,4) > 1
    S = (p+1)**2+2*(r+2)**2
    W = S+32*r+25
    assert 2*r+3 < 380*(r+1)**2
    assert W < 2**21*(p+1)**2
    assert 2*r+1 < 384*S
    return origin, bad, good, 2*(origin+good-bad)

assert escape(7) == (158984,144011,127806,285558)
assert escape(9) == (6136051,13846446,13039836,10658882)
assert escape(11) == (241442340,1402703027,1364907896,407294418)
assert escape(9)[1] > escape(9)[0]
print("escape and absolute-budget obstruction:", escape(9))

def box(d,g):
    numerator = (-1)**g+sum(2**j*comb(g+j,j) for j in range(2*d+1))
    assert numerator % 2**(2*d+1) == 0
    return numerator//2**(2*d+1)

for d in range(1,6):
    for g in range(15):
        assert box(d,g) == sum(comb(g-2*b+2*d-1,2*d-1)
                              for b in range(g//2+1))
C7 = 3200*8**3+760*8**2
assert C7 == 1687040 and len(str(box(7,C7-1))) == 76
print("finite profile count at distance 7:", box(7,C7-1))
print("PASS")