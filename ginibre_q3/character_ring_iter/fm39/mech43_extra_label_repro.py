"""FM-MECH43 (astra_max_ceres): {1,2} plus one label: B certificates (8 tables), IBP elimination (2)-(3), consumer reduction (4)-(5), counterexample to the moment-only route, 637 checks."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb

import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

argparse.ArgumentParser(
    description="FM-MECH43: exact B certificates and one-label-3 reduction"
).parse_args()


@lru_cache(None)
def row(a, b):
    out = {0: 1}
    for n in (1,)*a + (2,)*b:
        nxt = {}
        for j, c in out.items():
            for k in range(abs(j-n), j+n+1, 2):
                nxt[k] = nxt.get(k, 0) + c
        out = nxt
    return out


def spaces(d):
    for k in range(d+1):
        for piv in combinations(range(d), k):
            free = [
                (i, j) for i, p in enumerate(piv)
                for j in range(p+1, d) if j not in piv
            ]
            for mask in range(1 << len(free)):
                basis = [1 << p for p in piv]
                for bit, (i, j) in enumerate(free):
                    if mask >> bit & 1:
                        basis[i] |= 1 << j
                members = [0]
                for v in basis:
                    members += [x ^ v for x in members]
                yield members


all_spaces = list(spaces(7))
assert len(all_spaces) == 29212

cases = [
    (1,7,3), (3,5,3), (5,3,3), (7,1,3),
    (2,6,4), (4,4,4), (6,2,4), (8,0,4),
]
for a, b, n in cases:
    profiles = [
        (i, j) for i in range(0, a+1, 2) for j in range(b+1)
    ]
    ix = {p: k for k, p in enumerate(profiles)}
    target = [
        row(i,j).get(0,0)*row(a-i,b-j).get(n,0)
        for i, j in profiles
    ]
    sizes = [comb(a,i)*comb(b,j) for i,j in profiles]
    rhs = [v*s for v,s in zip(target, sizes)]

    basis = [
        (1 << i) ^ (1 << (a-1)) for i in range(a-1)
    ] + [1 << (a+j) for j in range(b)]

    classes = []
    for x in range(128):
        mask = 0
        for j, v in enumerate(basis):
            if x >> j & 1:
                mask ^= v
        classes.append(ix[
            ((mask & ((1 << a)-1)).bit_count(),
             (mask >> a).bit_count())
        ])

    dictionary = set()
    for members in all_spaces:
        hist = [0]*len(profiles)
        for x in members:
            hist[classes[x]] += 1
        if all(not hist[k] or target[k]
               for k in range(len(profiles))):
            dictionary.add(tuple(hist))

    cols = sorted(dictionary)
    mat = np.array(cols, dtype=float).T
    res = linprog(
        np.ones(len(cols)), A_eq=mat, b_eq=rhs,
        bounds=(0, None), method="highs"
    )
    assert res.success
    supp = [j for j,v in enumerate(res.x) if v > 1e-8]
    exact = Matrix([
        [cols[j][i] for j in supp] for i in range(len(rhs))
    ])
    sol, params = exact.gauss_jordan_solve(Matrix(rhs))
    assert params.rows == 0
    assert all(v >= 0 for v in sol)
    assert exact*sol == Matrix(rhs)
    print("B PASS", (a,b,n), "columns", len(cols), flush=True)


@lru_cache(None)
def catmoment(k):
    return 0 if k % 2 else comb(k, k//2)//(k//2+1)


def mul(p, q):
    out = {}
    for (i,j), v in p.items():
        for (k,l), w in q.items():
            key = (i+k, j+l)
            out[key] = out.get(key, 0) + v*w
    return {key:v for key,v in out.items() if v}


def mom(p, dx=0, dy=0):
    return sum(
        v*catmoment(i+dx)*catmoment(j+dy)
        for (i,j),v in p.items()
    )


Z = {(2,0):1, (0,2):1, (0,0):-2}
identities = inequalities = 0

for A in range(0, 13, 2):
    for E in range(0, 13, 2):
        p = {(0,0):1}
        for _ in range(A):
            p = mul(p, {(1,0):1, (0,1):1})
        for _ in range(E):
            p = mul(p, {(1,0):1, (0,1):-1})
        Ms, Qs = [], []
        for b in range(14):
            Ms.append(mom(p))
            Qs.append(mom(p, 1, 1))
            p = mul(p, Z)

        prefix = Q(0)
        for b in range(13):
            rhs = 4*(A-E)*Ms[b]
            if b:
                rhs += 12*b*Qs[b-1]
            assert (A+E+2*b+6)*Qs[b] == rhs
            prefix = (
                4*(A-E)*Ms[b] + 12*b*prefix
            ) / Q(A+E+2*b+6)
            assert prefix == Qs[b]
            identities += 1
            assert Ms[b+1] >= abs(Qs[b])
            inequalities += 1

print("Pearson and prefix identities:", identities)
print("Actual semicircle inequalities:", inequalities, "failures: 0")

q = Q(0)
for b in range(5):
    q = (8*2**b + 12*b*q)/Q(12+2*b)
assert q == Q(3872,105)
assert Q(32)-q == Q(-512,105)
print("Abstract positive-moment control: margin", Q(32)-q)
print("PASS")
