#!/usr/bin/env python3
"""FM-MECH143: exact averaged fusion, a positive sector, and obstructions.
Run from /home/yang/q3adjoint.  No repository modules or file writes.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import comb
from random import Random

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--samples", type=int, default=80)
args = ap.parse_args()

@lru_cache(None)
def cg(a, b):
    return tuple(range(abs(a-b), a+b+1, 2))

def step(F, z):
    n, eps = abs(z), 1 if z > 0 else -1
    G = defaultdict(int)
    for (i, j), v in F.items():
        for k in cg(i, n):
            G[k, j] += v
        for k in cg(j, n):
            G[i, k] += eps*v
    return {ij: v for ij, v in G.items() if v}

@lru_cache(maxsize=1024)
def row(B):
    if not B:
        return {(0, 0): 1}
    return step(row(B[:-1]), B[-1])

def expectation(B):
    return row(tuple(sorted(B))).get((0, 0), 0)

def fusion(labels):
    F = {0: 1}
    for n in labels:
        G = defaultdict(int)
        for i, v in F.items():
            for j in cg(i, n):
                G[j] += v
        F = dict(G)
    return F

def fused_value(B, mask):
    I = tuple(i for i in range(len(B)) if mask >> i & 1)
    C = tuple(B[i] for i in range(len(B)) if not (mask >> i & 1))
    sigma = (-1)**sum(B[i] < 0 for i in I)
    G = row(tuple(sorted(C)))
    return sum(v*(G.get((j, 0), 0) + sigma*G.get((0, j), 0))
               for j, v in fusion(tuple(abs(B[i]) for i in I)).items())

def residual(B, p):
    assert p >= 6 and p >= max(map(abs, B))
    assert all(-z not in B for z in B)
    W = sum(map(abs, B))
    assert (W-p) % 2 == 0
    delta = (W-p)//2
    assert delta >= 8 and max(map(abs, B)) <= delta
    assert sum(abs(z) >= 3 for z in B) >= 2
    return W, delta

# Independent monomial/Catalan expectation.
@lru_cache(None)
def U(n):
    return tuple((n-2*j, (-1)**j*comb(n-j, j))
                 for j in range(n//2+1))

@lru_cache(None)
def moment(n):
    return 0 if n % 2 else comb(n, n//2)//(n//2+1)

def direct(B):
    F = {(0, 0): 1}
    for z in B:
        eps = 1 if z > 0 else -1
        G = defaultdict(int)
        for (i, j), v in F.items():
            for k, w in U(abs(z)):
                G[i+k, j] += v*w
                G[i, j+k] += eps*v*w
        F = {ij: v for ij, v in G.items() if v}
    return sum(v*moment(i)*moment(j) for (i, j), v in F.items())

def in_positive_antisymmetric_cone(F):
    return all(i != j and F.get((j, i), 0) == -v
               and (i < j or v >= 0) for (i, j), v in F.items())

# Exactly two odd labels, arbitrary even-label background.
rng = Random(143)
for t in range(args.samples):
    C = tuple(rng.choice((-6, -4, -2, 2, 4, 6))
              for _ in range(rng.randrange(6)))
    a, b = rng.choice((1, 3, 5, 7)), rng.choice((1, 3, 5, 7))
    eps = rng.choice((-1, 1))
    eta = eps*(-1)**sum(z < 0 for z in C)
    B = C + (eps*a, eta*b)
    mask = (1 << len(C)) | (1 << (len(C)+1))
    assert expectation(B) == fused_value(B, mask)
    if t % 13 == 0:
        assert direct(B) == expectation(B)
print("two-odd fusion identities:", args.samples)

# Four minus factors: a,b odd; c,d even; d >= a+b+sum(extra).
# Arbitrarily many +2 factors are allowed.
channels = 0
for t in range(args.samples):
    a, b = rng.choice((1, 3, 5, 7)), rng.choice((1, 3, 5, 7))
    c = rng.choice((2, 4, 6, 8, 10))
    extra = tuple(rng.choice((4, 6)) for _ in range(rng.randrange(3)))
    d = a+b+sum(extra)+2*rng.randrange(3)
    power = rng.randrange(5)
    B = (-a, -b, -c, -d) + extra + (2,)*power
    total = 0
    for j in cg(a, b):
        H = row((-d,) + (() if j == 0 else (j,)) + extra + (2,)*power)
        scale = 2 if j == 0 else 1
        assert in_positive_antisymmetric_cone(H)
        assert all((i-jj) % 2 == 0 for i, jj in H)
        total += 2*scale*H.get((c, 0), 0)
        channels += 1
    assert expectation(B) == total >= 0
    if t % 19 == 0:
        assert direct(B) == total
print("four-minus sector:", args.samples, "profiles;", channels, "channels")

example = (-3, -5, -6, -8) + (2,)*5
assert residual((-3, -5, -6)+(2,)*5, 8) == (24, 8)
assert expectation(example) == fused_value(example, 3)
print("residual sector example:", expectation(example))

# All sign values for labels 1,...,8, by a separate subset-fusion/Walsh method.
labels = tuple(range(1, 9))
L, size = len(labels), 1 << len(labels)
mult = [0]*size
mult[0] = 1
for S in range(1, size):
    mult[S] = fusion(tuple(labels[i] for i in range(L) if S >> i & 1)).get(0, 0)
values = [mult[S]*mult[(size-1)^S] for S in range(size)]
h = 1
while h < size:
    for s in range(0, size, 2*h):
        for j in range(s, s+h):
            u, v = values[j], values[j+h]
            values[j], values[j+h] = u+v, u-v
    h *= 2
for K in range(size):
    B = tuple(-n if K >> i & 1 else n for i, n in enumerate(labels))
    assert values[K] == expectation(B)

B = (-1, 2, 3, 4, -5, 6, 7, 8)
K = (1 << 0) | (1 << 4)
assert residual(B[:-1], 8) == (28, 10)
assert values[K] == expectation(B) == direct(B) == 956
even_signs = [J for J in range(size) if J.bit_count() % 2 == 0]
assert min(values[J] for J in even_signs) == 956
minima = [J for J in even_signs if values[J] == 956]
pair_remainders = []
for i, j in combinations(range(L), 2):
    J = K ^ (1 << i) ^ (1 << j)
    flipped = tuple(-n if J >> z & 1 else n for z, n in enumerate(labels))
    assert direct(flipped) == values[J]
    delta = (values[K]-values[J])//2
    assert values[K]-values[J] == 2*delta
    assert delta < 0
    pair_remainders.append(delta)

block_values = []
for I in range(size):
    if I.bit_count() < 2:
        continue
    E = I
    acc = 0
    count = 0
    while True:
        if E.bit_count() % 2 == 0:
            acc += values[K ^ E]
            count += 1
        if E == 0:
            break
        E = (E-1) & I
    assert count == 1 << (I.bit_count()-1)
    Q = fused_value(B, I)
    assert acc == count*Q and Q > 956
    block_values.append(Q)

print("Walsh minima:", [tuple(labels[i] for i in range(L) if J >> i & 1)
                         for J in minima])
print("pair remainders:", len(pair_remainders),
      "range", min(pair_remainders), max(pair_remainders))
print("block fusions:", len(block_values),
      "range", min(block_values), max(block_values))

odd = (0, 2, 4, 6)
Q = [fused_value(B, (1 << odd[0]) | (1 << i)) for i in odd[1:]]
H = fused_value(B, sum(1 << i for i in odd))
assert sum(Q)-2*H == 956
print("four-odd budget:", Q, "single fused value", H)

# (D) is not refuted by the fusion obstruction.
background = B[:-1]
removals = [(i,) for i, z in enumerate(background) if z % 2 == 0]
removals += [(i, j) for i, j in combinations(range(len(background)), 2)
             if (background[i]+background[j]) % 2 == 0]
working = []
for I in removals:
    child = tuple(z for i, z in enumerate(background) if i not in I)
    value = 2*row(tuple(sorted(child))).get((8, 0), 0)
    if value <= 956:
        working.append((value, tuple(background[i] for i in I)))
assert working
print("working removal:", min(working))

# Two controls against stronger reductions.
C = (-3, -5, 6) + (2,)*4
F = C+(1, 7)
assert residual(C+(1,), 7) == (23, 8)
assert expectation(F) == direct(F) == 1468
assert fused_value(F, (1 << len(C)) | (1 << (len(C)+1))) == 2542
print("four-odd pair-fusion control: 1468 < 2542")

old = (-2, -2) + (-3,)*6
new = (-1, -3) + (-3,)*6
assert residual(old, 6) == (22, 8)
assert expectation(old+(6,)) == direct(old+(6,)) == 4052
assert expectation(new+(6,)) == direct(new+(6,)) == 8064
print("minus-spreading control: g_6 = 2026 < 4032")
print("PASS")
