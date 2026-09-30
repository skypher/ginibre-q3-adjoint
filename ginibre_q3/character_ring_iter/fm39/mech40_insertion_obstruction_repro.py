"""FM-MECH40 (astra_max_ceres): insertion preserving the parent autocorrelation is impossible (Gram obstruction; (2^9)->(2^10) at q^4, (2^15) at q = 0);
channel-reorganized certificates on the census."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb
from pathlib import Path

argparse.ArgumentParser(
    description="FM-MECH40: exact insertion obstruction and finite controls"
).parse_args()

# Load only exact fusion functions and subspace enumeration.
path = Path("ginibre_q3/character_ring_iter/fm39/sec113_hACq_L8_repro.py")
names = {"trim", "padd", "psub", "pmul", "pshift", "qbinom",
         "qfactorial", "linearization", "moment", "subspaces"}
tree = ast.parse(path.read_text())
env = {"lru_cache": lru_cache, "combinations": combinations}
exec(compile(ast.Module(body=[
    node for node in tree.body
    if isinstance(node, ast.FunctionDef) and node.name in names
], type_ignores=[]), str(path), "exec"), env)
add, mul, lin = (env[s] for s in ("padd", "pmul", "linearization"))

def moment(labels):
    return env["moment"](tuple(sorted(labels)))

def coef(p, k):
    return p[k] if k < len(p) else 0

def scale(p, k):
    return tuple(k*x for x in p)

def feasible(columns, target):
    """Exact Phase I; return the number of positive certificate terms."""
    assert min(target) >= 0
    active = [i for i, x in enumerate(target) if x]
    cols = sorted({
        tuple(c[i] for i in active) for c in columns
        if all(not x or target[i] for i, x in enumerate(c))
    })
    rhs = [F(target[i]) for i in active]
    nr, nc = len(rhs), len(cols)
    rows = [[F(cols[j][i]) for j in range(nc)]
            + [F(i == j) for j in range(nr)] for i in range(nr)]
    basis = [nc+i for i in range(nr)]
    red = [sum(row[j] for row in rows) if j < nc else F(0)
           for j in range(nc+nr)]
    objective = -sum(rhs)
    for _ in range(10000):
        en = next((j for j, x in enumerate(red) if x > 0), None)
        if en is None:
            break
        possible = [i for i in range(nr) if rows[i][en] > 0]
        assert possible
        le = min(possible, key=lambda i: (rhs[i]/rows[i][en], basis[i]))
        pivot = rows[le][en]
        rows[le] = [x/pivot for x in rows[le]]
        rhs[le] /= pivot
        for i in range(nr):
            if i == le:
                continue
            v = rows[i][en]
            if v:
                rows[i] = [x-v*y for x, y in zip(rows[i], rows[le])]
                rhs[i] -= v*rhs[le]
        v = red[en]
        objective += v*rhs[le]
        red = [x-v*y for x, y in zip(red, rows[le])]
        basis[le] = en
    else:
        raise AssertionError("Phase I iteration bound")
    assert objective == 0
    weights = [F(0)]*nc
    for i, j in enumerate(basis):
        if j < nc:
            weights[j] = rhs[i]
    assert min(weights) >= 0
    assert all(sum(weights[j]*cols[j][i] for j in range(nc)) == target[k]
               for i, k in enumerate(active))
    return sum(x > 0 for x in weights)

def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def kraw(m, t, s):
    return sum((-1)**j*binom(t, j)*binom(m-t, s-j)
               for j in range(s+1))

R = [moment((2,)*j) for j in range(11)]
first = None
for m in range(3, 9):
    failures = []
    for t in range(m+1):
        P = A = C = (0,)
        for s in range(m+1):
            w = kraw(m, t, s)
            P = add(P, scale(mul(R[s], R[m-s+1]), w))
            A = add(A, scale(mul(R[s], R[m-s+2]), w))
            C = add(C, scale(mul(R[s+1], R[m-s+1]), w))
        for k in range(max(map(len, (P, A, C)))):
            p, a, c = (coef(v, k) for v in (P, A, C))
            gap = 4*p*(a-p)-c*c
            if gap < 0:
                failures.append((k, t, p, a, c, gap))
    if m < 8:
        assert not failures
    else:
        assert len(failures) == 5
        first = min(failures)
assert first == (4, 0, 1020600, 4412408, 3788512, -506106194944)
print("First repeated-2 failure:", (8, first))

# Independent q=0 recurrence for Riordan moments.
state, R0 = {0: 1}, [1]
for _ in range(17):
    nxt = {}
    for h, v in state.items():
        for j in (h+2, h, h-2):
            if j == h and h == 0:
                continue
            if j < 0:
                continue
            nxt[j] = nxt.get(j, 0)+v
    state = nxt
    R0.append(state.get(0, 0))
for m in range(3, 16):
    p = sum(comb(m, s)*R0[s]*R0[m-s+1] for s in range(m+1))
    a = sum(comb(m, s)*R0[s]*R0[m-s+2] for s in range(m+1))
    c = sum(comb(m, s)*R0[s+1]*R0[m-s+1] for s in range(m+1))
    gap = 4*p*(a-p)-c*c
    assert (gap < 0) == (m == 15)
assert (p, a, c, gap) == (
    296028930, 758870302, 740559330, -370376798481060)
print("q=0 failure:", (15, p, a, c, gap))

codes = {
    9: [(80,44), (108,146,200,103,163,66), (160,109,251,180,110,124),
        (96,68,155), (150,), (253,), (137,189,194), (237,), (), (84,11)],
    10: [(224,197,285,461), (322,332,469,148,124,506),
         (14,165,101,195,430,140), (461,337,142,452,455),
         (342,322,242,370), (74,250,224,463,423), (328,229,477,216),
         (27,133,24), (447,368,325,378,467), (416,273),
         (404,238,149,404,123), (73,), (297,33,510), (258,94,122,348),
         (115,), (), (330,243), (447,72), (83,17), (34,), (479,304)]
}

def span(generators):
    H = {0}
    for v in generators:
        H |= {x ^ v for x in H}
    return H

entries = 0
for M, generators in codes.items():
    sizes = [1] + [comb(M,k)//(2 if 2*k == M else 1)
                   for k in range(1, M//2+1)]
    columns = []
    for gs in generators:
        hist = [0]*len(sizes)
        for x in span(gs):
            hist[min(x.bit_count(), M-x.bit_count())] += 1
        columns.append(hist)
    polynomials = [mul(R[k], R[M-k]) for k in range(len(sizes))]
    for degree in range(len(R[M])):
        feasible(columns, [sizes[k]*coef(p,degree)
                           for k,p in enumerate(polynomials)])
    entries += (1 << (M-1))*len(R[M])
    print("H_AC_q control:", M, len(R[M]), "coefficient certificates")
assert entries == 33024

def channel_data(mu, a, n):
    count = Counter(mu)
    types = sorted(count)
    profiles = list(product(*(range(count[t]+1) for t in types)))
    sizes, lefts, rights = [], [], []
    for prof in profiles:
        lefts.append(tuple(t for t,k in zip(types,prof) for _ in range(k)))
        rights.append(tuple(t for t,k in zip(types,prof)
                            for _ in range(count[t]-k)))
        size = 1
        for t,k in zip(types,prof):
            size *= comb(count[t],k)
        sizes.append(size)
    channels = []
    for j in range(n-a,n+a+1,2):
        ell = lin(a,n,(a+n-j)//2)
        channels.append([mul(ell, mul(moment(left),moment(right+(j,))))
                         for left,right in zip(lefts,rights)])
    cross = [mul(moment(left+(a,)),moment(right+(n,)))
             for left,right in zip(lefts,rights)]
    for i,(left,right) in enumerate(zip(lefts,rights)):
        total = (0,)
        for channel in channels:
            total = add(total,channel[i])
        assert total == mul(moment(left),moment(right+(a,n)))
    return profiles, sizes, channels, cross

spaces = {d:list(env["subspaces"](d)) for d in range(1,7)}
instances = certificates = terms = 0
for L in range(2,8):
    for labels in combinations_with_replacement(range(1,4),L):
        n = labels[-1]
        D = (sum(labels)-2*n)//2
        if sum(labels)%2 or not 0 <= D <= 4:
            continue
        for a in sorted(set(labels[:-1])):
            if a >= D:
                continue
            mu = list(labels[:-1])
            mu.remove(a)
            mu = tuple(mu)
            profiles,sizes,channels,cross = channel_data(mu,a,n)
            r, nc, m = len(profiles), len(channels), len(mu)
            index = {p:i for i,p in enumerate(profiles)}
            types = sorted(set(mu))
            classes = [index[tuple(sum(1 for i,t0 in enumerate(mu)
                                      if t0 == t and mask >> i & 1)
                                   for t in types)]
                       for mask in range(1 << m)]
            atoms = set()
            for H in spaces[m+1]:
                hist = [0]*(2*r)
                for x in H:
                    hist[r*(x >> m)+classes[x & ((1 << m)-1)]] += 1
                atoms.add(tuple(hist))
            columns = []
            for h in atoms:
                for j in range(nc):
                    col = [0]*((nc+1)*r)
                    col[j*r:(j+1)*r] = h[:r]
                    col[nc*r:] = h[r:]
                    columns.append(col)
            polys = sum(channels,[]) + cross
            for degree in range(max(map(len,polys))):
                target = [sizes[i%r]*coef(p,degree)
                          for i,p in enumerate(polys)]
                if any(target):
                    terms += feasible(columns,target)
                    certificates += 1
            instances += 1
assert (instances,certificates,terms) == (47,760,4473)
print("Small-label channel certificates:", instances,certificates,terms)

boundaries = [
    ((1,)*3+(2,)*2,1,2), ((1,)*3+(2,),1,2),
    ((1,)*6+(2,),2,4), ((1,)*7,1,6), ((1,)*12,1,3)]
orbits = boundary_entries = 0
for mu,a,n in boundaries:
    prof,sizes,channels,cross = channel_data(mu,a,n)
    orbits += len(prof)
    boundary_entries += sum(sizes)
assert (orbits,boundary_entries) == (55,4400)
print("Boundary channel identities:", orbits,boundary_entries)

# Exact CP allocation at the obstructing q^4 grade, using all channels.
allocation_atoms = [
    (0,2,4,9),(0,4,2,4),(0,8,1,3),(0,12,0,5),
    (1,0,0,1),(1,0,4,0),(1,1,3,2),(1,3,2,0),
    (1,4,0,3),(1,5,0,3),(1,5,1,3),(1,6,0,3),(1,12,0,3),
    (2,2,2,4),(2,5,0,2),(2,5,0,3),(2,6,3,4),(2,7,3,0),
    (2,7,5,2),(2,7,5,3),(2,8,1,6),(2,10,0,4),
    (2,16,1,0),(2,20,5,0)]
columns = []
for channel,ci,a,n in allocation_atoms:
    hist = [0]*18
    for x in span(codes[10][ci]):
        xn = x >> n & 1
        b = (x >> a & 1) ^ xn
        s = sum((x >> k & 1) ^ xn for k in range(10) if k not in (a,n))
        hist[9*b+s] += 1
    col = [0]*36
    col[9*channel:9*(channel+1)] = hist[:9]
    col[27:] = hist[9:]
    columns.append(col)
prof,sizes,channels,cross = channel_data((2,)*8,2,2)
polys = sum(channels,[]) + cross
target = [sizes[i%9]*coef(p,4) for i,p in enumerate(polys)]
assert feasible(columns,target) == 24
print("Witness q^4 reorganized-channel certificate: 24 terms")
print("All assertions passed.")