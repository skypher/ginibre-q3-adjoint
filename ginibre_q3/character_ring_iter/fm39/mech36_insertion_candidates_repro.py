"""FM-MECH36 (astra_max_ceres): insertion-stable candidates for H_AC: common-spin Gram, translation lift, scalar allocation (all fail), subset-dependent allocation repair; bounded screen."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement
from math import comb

argparse.ArgumentParser(
    description="FM-MECH36 exact obstruction and repair verifier"
).parse_args()

def fusion_rows(mu):
    rows = [{0: 1}]
    for label in mu:
        for old in rows[:]:
            new = {}
            for j, a in old.items():
                for k in range(abs(j-label), j+label+1, 2):
                    new[k] = new.get(k, 0) + a
            rows.append(new)
    return rows

def wh(a):
    a = list(a)
    step = 1
    while step < len(a):
        for x in range(0, len(a), 2*step):
            for y in range(x, x+step):
                a[y], a[y+step] = a[y]+a[y+step], a[y]-a[y+step]
        step *= 2
    return a

def data(mu, h, n):
    r = fusion_rows(mu)
    last = len(r)-1
    F = {
        j: [r[S].get(0, 0)*r[last ^ S].get(j, 0)
            for S in range(last+1)]
        for j in range(abs(h-n), h+n+1, 2)
    }
    C = [r[S].get(h, 0)*r[last ^ S].get(n, 0)
         for S in range(last+1)]
    return F, C

def scalar_caps(F, C):
    ch = wh(C)
    return {
        j: min(Q(a, abs(c)) for a, c in zip(wh(f), ch) if c)
        for j, f in F.items()
    }

def translation_budget(F, C):
    support = {x for x, a in enumerate(C) if a}
    budget = 0
    for f in F.values():
        supp = {x for x, a in enumerate(f) if a}
        usable = supp and any(
            all(x ^ s in support for x in supp)
            for s in range(len(C))
        )
        if usable:
            budget += sum(f)
    return budget

def square(p, N):
    a = [Q(0)]*N
    for x, b in p.items():
        for y, c in p.items():
            a[x ^ y] += b*c
    return a

def span(basis):
    H = [0]
    for b in basis:
        assert b not in H
        H += [x ^ b for x in H]
    return H

def orbit(x, widths):
    ans = []
    for width in widths:
        ans.append((x & ((1 << width)-1)).bit_count())
        x >>= width
    return tuple(ans)

def code_average(widths, terms):
    N = 1 << sum(widths)
    out = [Q(0)]*N
    for a, basis in terms:
        assert a >= 0
        counts = Counter(orbit(x, widths) for x in span(basis))
        for x in range(N):
            key = orbit(x, widths)
            size = 1
            for w, k in zip(widths, key):
                size *= comb(w, k)
            out[x] += a*Q(counts.get(key, 0), size)
    return out

# Common-spin Gram obstruction: positive odd external spins.
r = fusion_rows((2,))[-1]
spins = (1, 3, 5)
M = [[sum(r.get(j, 0)
          for j in range(abs(h-n), h+n+1, 2))
      for n in spins] for h in spins]
z = (1, -1, 1)
assert M == [[1, 1, 0], [1, 1, 1], [0, 1, 1]]
assert sum(z[i]*M[i][j]*z[j]
           for i in range(3) for j in range(3)) == -1
print("common-spin Gram: quadratic witness -1")

# Translation obstruction and asymmetric channel-1 CP repair.
F, C = data((1, 1, 1), 1, 2)
assert F[1] == [2, 0, 0, 1, 0, 1, 1, 0]
assert F[3] == [1, 0, 0, 0, 0, 0, 0, 0]
assert C == [0, 1, 1, 0, 1, 0, 0, 0]
assert translation_budget(F, C) == 1 and sum(C) == 3
p = {1: Q(1), 2: Q(1), 4: Q(1), 8: Q(1)}
assert [a/2 for a in square(p, 16)] == F[1] + C
print("translation lift: budget 1 < 3; channel-1 CP repair PASS")

# Scalar-allocation obstruction.
F, C = data((1,)*5 + (2,), 1, 2)
caps = scalar_caps(F, C)
assert caps == {1: Q(7, 15), 3: Q(1, 2)}
assert sum(caps.values()) == Q(29, 30)
assert (wh(F[1])[32], wh(C)[32],
        wh(F[3])[0], wh(C)[0]) == (14, 30, 60, 120)
print("scalar allocation: 7/15 + 1/2 = 29/30 < 1")

# Nonuniform CP repair.
# Bits 0..4: old labels 1; bit 5: old label 2;
# bit 6: inserted label 1.
# Each term is its coefficient times the S5-average of 1_span(basis).
cert = {
    1: [
        (Q(13,15),   (24,47)),
        (Q(77,60),   (6,10,18,35)),
        (Q(16,45),   (24,104)),
        (Q(14,15),   (12,18,37,71)),
        (Q(127,90),  (6,10,17,97)),
        (Q(29,40),   (6,10,18,35,65)),
        (Q(1),       (24,72)),
        (Q(89,60),   (12,20,38,68)),
        (Q(113,120), (6,10,17,35,66)),
    ],
    3: [
        (Q(2,3), (95,)),
        (Q(10,3),(80,)),
        (Q(5,6), (12,18,38,65)),
        (Q(5,6), (6,10,18,65)),
        (Q(5,3), (12,18,38,68)),
        (Q(5,3), (6,10,18,66)),
        (Q(1),   (3,5,9,17,65)),
    ]
}
C1 = [(Q(7,9) if S & 32 else Q(1,3))*a
      for S, a in enumerate(C)]
for j, cross in [
    (1, C1),
    (3, [Q(a)-b for a, b in zip(C, C1)])
]:
    assert code_average((5,1,1), cert[j]) == F[j] + cross
print("nonuniform channel CP repair: 16 exact code terms PASS")

# Required boundaries.
boundaries = [
    ((1,)*7, 1, 6),
    ((2,)*3, 2, 2),
    ((1,)*6 + (2,), 2, 4),
]
children = []
for mu, h, n in boundaries:
    F, C = data(mu, h, n)
    A = [sum(f[S] for f in F.values()) for S in range(len(C))]
    assert min(wh(A+C)) >= 0
    children.append(A+C)
    print("boundary", mu, "h,n", (h,n),
          "cross mass", sum(C),
          "translation budget", translation_budget(F,C))

# Explicit H_AC checks at all three boundaries.
p = {1 << i: Q(1) for i in range(8)}
g = [a/2 for a in square(p, 256)]
g[0] += 3
assert g == children[0]

p = {1 << i: Q(1) for i in range(4)}
g = [a/2 for a in square(p, 16)]
for i, j, k in combinations(range(4), 3):
    p = {1 << i: Q(1), (1 << j) ^ (1 << k): Q(1)}
    g = [a+b/2 for a, b in zip(g, square(p, 16))]
assert g == children[1]

known = [
    (Q(45,16), (3,5,9,48)),
    (Q(35,16), (3,5,9,48,81)),
    (Q(105,8), (3,5,9,81)),
    (Q(5,8),   (3,5,24,40)),
    (Q(5),     (3,5,192)),
    (Q(5),     (3,12,53,81)),
    (Q(15,2),  (3,12,69)),
    (Q(15,4),  (3,12,197)),
]
assert code_average((6,2), known) == children[2]

# Exact bounded screen.
checks = scalar_fail = translation_fail = 0
for L in range(2, 8):
    for mu in combinations_with_replacement(range(1,4), L):
        for h, n in [(1,2), (2,2), (2,3)]:
            F, C = data(mu, h, n)
            if not any(C):
                continue
            A = [sum(f[S] for f in F.values())
                 for S in range(len(C))]
            assert min(wh(A+C)) >= 0
            checks += 1
            scalar_fail += sum(scalar_caps(F,C).values()) < 1
            translation_fail += translation_budget(F,C) < sum(C)

assert (checks, scalar_fail, translation_fail) == (167, 10, 112)
print("exact screen:", checks, "instances;",
      scalar_fail, "scalar failures;",
      translation_fail, "translation-budget failures")
print("ALL ASSERTIONS PASS")
