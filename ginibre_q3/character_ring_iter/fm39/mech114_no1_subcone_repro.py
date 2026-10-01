import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
from random import Random

ap = argparse.ArgumentParser(
    description="FM-MECH114 exact certificates; no files written.")
ap.parse_args()

def grow(B, p):
    st = {(0, 0): 1}
    remain = sum(abs(z) for z in B)
    for z in B:
        n, eps = abs(z), 1 if z > 0 else -1
        remain -= n
        out = defaultdict(int)
        for (i, j), v in st.items():
            if j <= remain:
                for h in range(abs(i-n), i+n+1, 2):
                    if abs(h-p) <= remain:
                        out[h, j] += v
            if abs(i-p) <= remain:
                for h in range(abs(j-n), j+n+1, 2):
                    if h <= remain:
                        out[i, h] += eps*v
        st = {ij: v for ij, v in out.items() if v}
    return st.get((p, 0), 0)

def even(B):
    B = list(B)
    if sum(z < 0 for z in B) % 2 or sum(map(abs, B)) % 2:
        return 0
    p = max(map(abs, B))
    B.pop(next(i for i, z in enumerate(B) if abs(z) == p))
    return 2*grow(B, p)

def row(N, A):
    c = [Q(1)]
    for j in range(N):
        c.append(
            (-A*c[j]-(N-j+1)*(c[j-1] if j else 0))/(j+1))
    return c

def F(c, n, m, d):
    C = lambda j: c[j] if 0 <= j < len(c) else 0
    D = lambda j: C(j)**2-C(j-1)*C(j+1)
    J = lambda j, h: C(j)*C(j-h)-C(j+1)*C(j-h-1)
    return (
        sum(D(d-j)-D(d-n-1-j) for j in range(m+1))
        - J(d,n)+J(d-m-1,n)-J(d,m)+J(d-n-1,m)
        + C(d-m)*C(d-n)-C(d-m-1)*C(d-n-1)
        - C(d+1)*C(d-n-m-1)+C(d)*C(d-n-m-2)
    )

# Requested formula versus independent SU(2) multiplication.
rng = Random(114)
checks = 0
while checks < 120:
    a, b = rng.randrange(37), rng.randrange(13)
    n, m = rng.randrange(3,15), rng.randrange(3,15)
    N = a+2*b
    lo, hi = max(8,n,m), (N+min(n,m))//2
    if lo > hi:
        continue
    d = rng.randrange(lo,hi+1)
    p = N+n+m-2*d
    if p in (n,m) and (a+b) % 2 == 0:
        continue
    g = grow([-1]*a+[-2]*b+[-n,-m],p)
    assert F(row(N,a),n,m,d) == g >= 0
    checks += 1
assert F(row(66,64),3,3,8) == 2501223044263151142
print("K=2 FORMULA:",checks,"independent comparisons")

# Parity fusion, every compatible core sign, including r=0.
fusion_checks = 0
for core in combinations_with_replacement(range(3,10),3):
    odd = [i for i,n in enumerate(core) if n % 2]
    if len(odd) != 2:
        continue
    i,j = odd
    k = next(h for h in range(3) if h not in odd)
    A,B,C = core[i],core[j],core[k]
    for b in range(6):
        for ep in product((-1,1),repeat=3):
            if (-1)**b*ep[0]*ep[1]*ep[2] != 1:
                continue
            base = [-2]*b
            lhs = even(base+[ep[h]*core[h] for h in range(3)])
            tau = ep[i]*ep[j]
            rhs = 0
            for r in range(abs(A-B),A+B+1,2):
                term = (
                    (1+tau)*even(base+[ep[k]*C]) if r == 0
                    else even(base+[tau*r,ep[k]*C]))
                assert term >= 0
                rhs += term
            assert lhs == rhs >= 0
            fusion_checks += 1
assert fusion_checks == 720
assert even([-2]*6+[-7,-7,8]) == 3390
print("PARITY FUSION:",fusion_checks,"identities")

# Left-boundary continuation: exact negative normalized interval.
c = row(500,Q(997,2))
q = F(c,3,3,250)/c[250]**2
assert -Q(11,1000) < q < -Q(1,100)
for b in (0,1,2):
    assert F(row(500,500-2*b),3,3,250) > 0

# Reciprocal interior continuation.
def centered(N,A,eps):
    h = N//2
    c = [Q(0)]*(N+1)
    if N % 2:
        c[h],c[h+1],j0 = Q(1),Q(eps),h
    elif eps == -1:
        c[h],c[h-1],j0 = Q(0),Q(1),h-1
    else:
        c[h],c[h-1],j0 = Q(1),-A/Q(2*(h+1)),h-1
    for j in range(j0,0,-1):
        c[j-1] = (-A*c[j]-(j+1)*c[j+1])/(N-j+1)
    for j in range(h+1):
        c[N-j] = eps*c[j]
    for j in range(1,N):
        assert (j+1)*c[j+1] == -A*c[j]-(N-j+1)*c[j-1]
    assert all(c[N-j] == eps*c[j] for j in range(N+1))
    return c

c = centered(66,Q(196,3),-1)
q = F(c,4,8,11)/c[11]**2
assert -Q(127,500) < q < -Q(253,1000)
assert c[1] != -Q(196,3)*c[0]
print("CONTINUOUS RELAXATIONS: two exact negative witnesses")

# Every pair-free list with delta=8, k=3, maximum 6 or 7.
counts, minus = Counter(), Counter()
minima = {}
for p in (6,7):
    for low in combinations_with_replacement(range(3,p+1),3):
        core = low+(p,)
        t = 2*p+16-sum(core)
        for b in range(t//2+1):
            a = t-2*b
            ct = Counter([1]*a+[2]*b+list(core))
            types = sorted(ct)
            for ep in product((-1,1),repeat=len(types)):
                if sum(ct[n] for n,e in zip(types,ep) if e < 0) % 2:
                    continue
                W = tuple(
                    e*n for n,e in zip(types,ep)
                    for _ in range(ct[n]))
                assert sum(map(abs,W)) == 2*p+16
                assert sum(abs(z)>=3 for z in W) == 4
                v = even(W)
                assert v > 0,(W,v)
                counts[p] += 1
                minus[p] += all(z < 0 for z in W)
                if p not in minima or v < minima[p][0]:
                    minima[p] = (v,W)

assert dict(counts) == {6:1180,7:2254}
assert dict(minus) == {6:50,7:84}
assert minima[6][0] == 140 and minima[7][0] == 54

examples = [
    ((-1,)*6+(-5,-5,-6,-6),4146),
    ((-1,)*5+(-2,-5,-5,-5,-6),2518),
    ((-1,)*4+(2,-5,-5,-6,-6),1446),
    ((-2,-2,-5,-7,-7,-7),114),
]
for W,v in examples:
    assert even(W) == v

print("FRONTIER:",dict(counts),"ALL MINUS:",dict(minus))
print("MINIMA:",minima)
print("ALL CHECKS PASS")