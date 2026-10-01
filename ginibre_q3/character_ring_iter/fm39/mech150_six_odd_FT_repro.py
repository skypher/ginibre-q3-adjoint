import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
from math import comb
from random import Random

ap = argparse.ArgumentParser(description="FM-MECH150 exact verifier; memory only")
ap.add_argument("--samples", type=int, default=300)
args = ap.parse_args()

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

@lru_cache(maxsize=100000)
def mu(ns, p=0):
    if p < 0 or (sum(ns)-p)%2 or p > sum(ns):
        return 0
    if not ns:
        return int(p == 0)
    if len(ns) == 1:
        return int(ns[0] == p)
    r, s = len(ns), (sum(ns)-p)//2
    return sum((-1)**K.bit_count() *
               C(s-sum(ns[i]+1 for i in range(r) if K>>i&1)+r-2, r-2)
               for K in range(1<<r))

def char_table(word):
    out = {(0, 0): 1}
    for x in word:
        n, eps = abs(x), 1 if x > 0 else -1
        nxt = defaultdict(int)
        for (a, b), v in out.items():
            for j in range(abs(a-n), a+n+1, 2):
                nxt[j, b] += v
            for j in range(abs(b-n), b+n+1, 2):
                nxt[a, j] += eps*v
        out = {ab: v for ab, v in nxt.items() if v}
    return out

def phi(word):
    return char_table(tuple(word)).get((0, 0), 0)

def fwht(v):
    v = list(v)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2*h):
            for j in range(i, i+h):
                a, b = v[j], v[j+h]
                v[j], v[j+h] = a+b, a-b
        h *= 2
    return v

PAIRS = tuple(combinations(range(6), 2))

def data(ns, p):
    M = mu(ns, p)
    A, W = {}, {}
    for i, j in PAIRS:
        rest = tuple(ns[k] for k in range(6) if k not in (i, j))
        A[i, j] = mu(rest, p) if ns[i] == ns[j] else 0
        W[i, j] = mu(rest, 0) if ns[i]+ns[j] >= p else 0
    h = mu(ns[:4], p)
    return M, A, W, h

def graph_value(info, K):
    M, A, W, h = info
    eps = tuple(-1 if K>>i&1 else 1 for i in range(6))
    s = (-1)**K.bit_count()
    return M+sum((A[i,j]+s*W[i,j])*eps[i]*eps[j] for i,j in PAIRS)

def graph_flips(info, K):
    M, A, W, h = info
    x = tuple(-1 if K>>i&1 else 1 for i in range(6))
    s = (-1)**K.bit_count()
    qa = sum(A[i,j]*x[i]*x[j] for i,j in PAIRS)
    qw = sum(W[i,j]*x[i]*x[j] for i,j in PAIRS)
    av = [sum(A[min(i,j),max(i,j)]*x[i]*x[j]
              for j in range(6) if j != i) for i in range(6)]
    wv = [sum(W[min(i,j),max(i,j)]*x[i]*x[j]
              for j in range(6) if j != i) for i in range(6)]
    single = [av[i]+s*(qw-wv[i]) for i in range(6)]
    double = [av[i]+av[j]-2*A[i,j]*x[i]*x[j]
              +s*(wv[i]+wv[j]-2*W[i,j]*x[i]*x[j])
              for i,j in PAIRS]
    return single+double

def check_profile(ns, p, check_signs=False):
    a,b,c,d,e,f = ns
    assert a < b < c < d <= e <= f
    assert all(n%2 for n in ns) and p%2 == 0 and p >= f
    delta = (sum(ns)-p)//2
    assert f <= delta
    r, t = (a+b+c+d-p)//2, (f-e)//2
    assert r >= t >= 0
    info = data(ns, p)
    M, A, W, h = info
    P2, T = sum(A.values()), sum(W.values())
    assert M >= (e+1)*h
    if r <= 3:
        assert h >= (1,3,5,7)[r]
        cap = a+1+(2,7,17,27)[r]
        if r == 0:
            assert e == f and h == 1 and P2 >= h
    else:
        lower = 15 if a >= 5 else 14 if a == 3 else 8 if b == 3 else 9
        assert h >= lower
        cap = 10*(a+1)+4*(b+1)+(c+1)
    assert T <= cap
    assert e*h+P2-cap >= 4
    assert M+P2-T-h >= 4
    if check_signs:
        weights = []
        for S in range(64):
            left = tuple(ns[i] for i in range(6) if S>>i&1)
            right = tuple(ns[i] for i in range(6) if not S>>i&1)
            weights.append(mu(left,0)*mu(right,p))
        values = fwht(weights)
        for K in range(64):
            assert values[K] == graph_value(info,K)
            ds = graph_flips(info,K)
            flips = [(K^(1<<i)) for i in range(6)]
            flips += [K^(1<<i)^(1<<j) for i,j in PAIRS]
            assert all(2*z == values[K]-values[J] for z,J in zip(ds,flips))
            if all(ns[i] != ns[j] or ((K>>i)&1) == ((K>>j)&1)
                   for i,j in PAIRS):
                assert values[K] >= h+4
    return r

counts = [0]*5
profiles = 0
odd = tuple(range(1,18,2))
for bot in combinations(odd,4):
    for e in odd:
        if e < bot[-1]:
            continue
        for f in odd:
            if f < e:
                continue
            ns = bot+(e,f)
            for p in range(f+1, sum(ns)-2*f+1, 2):
                r = check_profile(ns,p,check_signs=True)
                counts[min(r,4)] += 1
                profiles += 1
print("exhaustive unsigned profiles:", profiles, "defect bins:", counts)

rng = Random(150)
accepted = 0
while accepted < args.samples:
    bot = tuple(sorted(rng.sample(range(1,400,2),4)))
    e = rng.randrange(bot[-1],602,2)
    f = rng.randrange(e,802,2)
    ns = bot+(e,f)
    lo, hi = f+1, sum(ns)-2*f
    if lo > hi:
        continue
    p = rng.randrange(lo,hi+1,2)
    check_profile(ns,p)
    accepted += 1
print("additional exact label profiles:", accepted)

examples = [
    ((-1,-3,-5,7,-9,-11),12),
    ((-1,-3,-5,7,-11,-13),14),
    ((-1,3,5,7,9,11),14),
    ((-1,3,5,7,9,9),12),
]
for B,p in examples:
    for reflect in (1,-1):
        B = tuple(reflect*x for x in B)
        ns = tuple(map(abs,B))
        K = sum(1<<i for i,x in enumerate(B) if x < 0)
        info = data(ns,p)
        value = graph_value(info,K)
        ds = graph_flips(info,K)
        assert max(ds) < 0
        sig = (-1)**K.bit_count()
        Lam = B+(sig*p,)
        tab = char_table(B)
        assert value == tab.get((p,0),0)
        assert 2*value == phi(Lam)
        child = char_table(B[:4]).get((p,0),0)
        assert child == info[3] and value-child >= 4
        if reflect == 1:
            prs = [(i,6) for i in range(6)]+list(PAIRS)
            for (i,j), predicted in zip(prs,ds):
                rest = tuple(x for k,x in enumerate(Lam) if k not in (i,j))
                table = char_table(rest)
                actual = (1 if Lam[j]>0 else -1)*table.get(
                    (abs(Lam[i]),abs(Lam[j])),0)
                assert actual == predicted
    print("no-flip example:", ns, p, "g =", value,
          "child =", child, "flip range =", (min(ds),max(ds)))

B,p = (7,9,19,-35,37,-45),52
ns = tuple(map(abs,B))
K = sum(1<<i for i,x in enumerate(B) if x < 0)
info = data(ns,p)
assert graph_value(info,K) == 48987
assert max(graph_flips(info,K)) < 0
assert info[3] == 52

# Normalization of the Walsh pair sum.
Lam = (1,1,1,1)
full = 15
base = phi(Lam)
total = 0
for i,j in combinations(range(4),2):
    flip = tuple(-x if k in (i,j) else x for k,x in enumerate(Lam))
    total += base-phi(flip)
walsh = 0
for S in range(16):
    left = tuple(1 for i in range(4) if S>>i&1)
    right = tuple(1 for i in range(4) if not S>>i&1)
    walsh += S.bit_count()*(4-S.bit_count())*mu(left,0)*mu(right,0)
assert base == 10 and total == 48 and walsh == 24
assert total == 2*walsh

# The selector is not asserted outside the theorem's sector.
B = (-1,)*4+(-2,)+(-3,)*8+(-4,)
g = char_table(B).get((6,0),0)
rest = list(B)
rest.remove(-2)
rest.remove(-4)
child = char_table(tuple(rest)).get((6,0),0)
assert (g, child) == (996550, 1077706)
assert -char_table(B[2:]+(6,)).get((1,1),0) == 240992
print("outside-sector TopPair counterexample:", g, child)
print("PASS")
