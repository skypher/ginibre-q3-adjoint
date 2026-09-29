#!/usr/bin/env python3
"""Verifier for Theorem T1-6: T1 (FM3, r = 2, H-only) for every word with at most six factors.

N_6 := (sum_kappa T1(kappa) u^kappa) * prod_{i<j<=6} (1 - u_i u_j).
 (a) N_6 = 5 + p_2 - 2 e_2 + 3 e_4 - e_1 e_5 - 4 e_6 + e_2 e_6 - e_6^2.
     Checked here by computing T1 exactly (C2 Weyl character formula) for every 6-part
     composition of size <= D (default 16; the note records D = 26, which covers the u-degree
     bound 26 coming from the K-polynomial of the 8-vertex Pfaffian ring, a = -2n).
 (b) N_6 = (1/12) [ sum over Hamiltonian 6-cycles C of prod_{e in C} (1 - u_e)
                  + sum over vertices v and 5-cycles C on the other five vertices of
                    u_v^2 prod_{e in C} (1 - u_e) ],
     checked exactly with sympy.  Hence
       12 T1(kappa) = sum_C G_{K6 - C}(kappa) + sum_(v,C) G_{K6 - C}(kappa - 2 e_v) >= 0,
     where G_X(kappa) counts multigraphs with edges in X and degree sequence kappa.
Usage: verify_t1_six_factors.py [D]
"""
import itertools, sys
from collections import defaultdict
import sympy as sp

s = 6
u = sp.symbols('u1:7')
E = [(i, j) for i in range(s) for j in range(i + 1, s)]
e = lambda k: sum(sp.prod([u[t] for t in S]) for S in itertools.combinations(range(s), k))
p2 = sum(x**2 for x in u)
N6 = sp.expand(5 + p2 - 2 * e(2) + 3 * e(4) - e(1) * e(5) - 4 * e(6) + e(2) * e(6) - e(6)**2)

# (b) cycle decomposition
def cycles(vertices):
    vs = list(vertices); first = vs[0]; out = set()
    for perm in itertools.permutations(vs[1:]):
        if perm[0] > perm[-1]:
            continue
        cyc = (first,) + perm
        out.add(frozenset(tuple(sorted((cyc[i], cyc[(i + 1) % len(cyc)]))) for i in range(len(cyc))))
    return out
w = lambda F: sp.prod([1 - u[i] * u[j] for i, j in F])
total = sum(w(C) for C in cycles(range(6)))
for v in range(6):
    total += sum(u[v]**2 * w(C) for C in cycles([x for x in range(6) if x != v]))
assert len(cycles(range(6))) == 60
assert sp.expand(total / 12 - N6) == 0
print("(b) N_6 = (1/12)[Hamiltonian 6-cycles + (vertex, 5-cycle) terms]: exact identity OK")

# (a) N_6 from exact T1 values
D = int(sys.argv[1]) if len(sys.argv) > 1 else 16
WEYL = []
for perm in [(0, 1), (1, 0)]:
    for s0 in (1, -1):
        for s1 in (1, -1):
            WEYL.append((perm, (s0, s1), (1 if perm == (0, 1) else -1) * s0 * s1))
def hk(k):
    P = defaultdict(int)
    for a in range(k + 1):
        for i in range(-a, a + 1, 2):
            for j in range(-(k - a), k - a + 1, 2):
                P[(i, j)] += 1
    return P
def mul(A, B):
    C = defaultdict(int)
    for (i, j), x in A.items():
        for (k, l), y in B.items():
            C[(i + k, j + l)] += x * y
    return C
def mult(P, a, b):
    v = (a + 2, b + 1)
    return sum(det * P.get((sg[0] * v[perm[0]] - 2, sg[1] * v[perm[1]] - 1), 0) for perm, sg, det in WEYL)
cache = {(): {(0, 0): 1}}
def wpoly(lam):
    if lam not in cache:
        cache[lam] = mul(wpoly(lam[:-1]), hk(lam[-1]))
    return cache[lam]
def parts(n, m, l):
    if n == 0:
        yield (); return
    if l == 0:
        return
    for k in range(min(n, m), 0, -1):
        for r in parts(n - k, k, l - 1):
            yield (k,) + r
T1 = {}
for n in range(0, D + 1, 2):
    for lam in parts(n, n, s):
        P = wpoly(lam)
        T1[lam] = 5 * mult(P, 0, 0) - 3 * mult(P, 1, 1) + mult(P, 2, 0)
ser = {}
for comp in itertools.product(range(D + 1), repeat=s):
    if sum(comp) > D or sum(comp) % 2:
        continue
    v = T1.get(tuple(sorted([c for c in comp if c], reverse=True)), 0)
    if v:
        ser[comp] = v
for (i, j) in E:
    new = dict(ser)
    for comp, v in ser.items():
        c2 = list(comp); c2[i] += 1; c2[j] += 1
        if sum(c2) <= D:
            c2 = tuple(c2); new[c2] = new.get(c2, 0) - v
    ser = {k: v for k, v in new.items() if v}
Npoly = sp.Poly(N6, *u).as_dict()
target = {k: int(v) for k, v in Npoly.items() if sum(k) <= D}
assert ser == target, "N_6 mismatch"
print(f"(a) N_6 formula matches exact T1 series through degree {D}")
print("Theorem T1-6 verified (given the degree bound for D < 26).")
