#!/usr/bin/env python3
"""Verifier for Theorem T1-5: FM3 at level r = 2 (T1) for every word with at most
five symmetric-power factors.

Claim.  For kappa in N^5, T1(kappa) = 5 m_(0,0) - 3 m_(1,1) + m_(2,0) >= 0 for
M = Sym^k1 W (x) ... (x) Sym^k5 W, W = C^4, G = Sp(4).

Proof being checked.
 (1) T1 = 8 m0(kappa) + 4 m0(kappa,2) - 3 m0(kappa,1,1) (W (x) W = 1 + U + Sym^2 W).
 (2) First/second fundamental theorems for Sp4: invariants of n vectors are C[omega]/I_6,
     I_6 = ideal of 6 x 6 Pfaffians.  I_6 = 0 for n <= 5; I_6 = (Pf) for n = 6; for n = 7,
     Buchsbaum--Eisenbud: 0 -> R(-(2^7)) -> (+)_l R(-(1^7)-e_l) -> (+)_k R(-(1^7)+e_k) -> R.
     Extracting the coefficients of t^2 (n = 6) and t1 t2 (n = 7) gives
       sum_kappa T1(kappa) u^kappa = N / prod_{i<j<=5} (1 - u_i u_j),
       N = 5 + p_2 - 2 e_2 + 3 e_4 - e_1 e_5.
 (3) N = sum_i c_i u^(m_i) prod_{e in F_i} (1 - u_e), with 56 rational c_i > 0 (CERT5 below),
     so N / prod_e (1 - u_e) = sum_i c_i u^(m_i) / prod_{e not in F_i} (1 - u_e) has
     nonnegative coefficients.
The script checks (3) exactly and checks (2) against a direct Weyl-character computation.
"""
import itertools
from collections import defaultdict
from fractions import Fraction
import sympy as sp

CERT5 = [
    ('2/133', ((0, 1), (0, 2), (1, 3)), (0, 0, 0, 0, 2)),
    ('5/133', ((0, 1), (0, 2), (1, 4)), (0, 0, 0, 2, 0)),
    ('40/133', ((0, 1), (0, 3), (1, 2)), (0, 0, 0, 0, 2)),
    ('6/133', ((0, 1), (0, 3), (1, 4)), (0, 0, 2, 0, 0)),
    ('34/133', ((0, 1), (0, 3), (2, 3)), (0, 0, 0, 0, 2)),
    ('29/133', ((0, 1), (0, 3), (3, 4)), (0, 0, 2, 0, 0)),
    ('1/7', ((0, 1), (0, 4), (1, 2)), (0, 0, 0, 2, 0)),
    ('10/133', ((0, 1), (0, 4), (2, 4)), (0, 0, 0, 2, 0)),
    ('8/133', ((0, 1), (0, 4), (3, 4)), (0, 0, 2, 0, 0)),
    ('18/133', ((0, 1), (1, 3), (2, 3)), (0, 0, 0, 0, 2)),
    ('17/133', ((0, 1), (1, 4), (2, 4)), (0, 0, 0, 2, 0)),
    ('31/133', ((0, 1), (1, 4), (3, 4)), (0, 0, 2, 0, 0)),
    ('13/133', ((0, 2), (0, 3), (1, 2)), (0, 0, 0, 0, 2)),
    ('16/133', ((0, 2), (0, 3), (1, 3)), (0, 0, 0, 0, 2)),
    ('8/133', ((0, 2), (0, 3), (2, 4)), (0, 2, 0, 0, 0)),
    ('5/19', ((0, 2), (0, 4), (3, 4)), (0, 2, 0, 0, 0)),
    ('10/133', ((0, 2), (1, 2), (1, 3)), (0, 0, 0, 0, 2)),
    ('51/133', ((0, 2), (1, 2), (1, 4)), (0, 0, 0, 2, 0)),
    ('32/133', ((0, 3), (0, 4), (1, 4)), (0, 0, 2, 0, 0)),
    ('9/133', ((0, 3), (0, 4), (2, 3)), (0, 2, 0, 0, 0)),
    ('81/133', ((0, 3), (2, 3), (2, 4)), (0, 2, 0, 0, 0)),
    ('31/133', ((0, 4), (1, 2), (2, 4)), (0, 0, 0, 2, 0)),
    ('27/133', ((0, 4), (1, 3), (3, 4)), (0, 0, 2, 0, 0)),
    ('3/19', ((1, 2), (1, 3), (2, 4)), (2, 0, 0, 0, 0)),
    ('36/133', ((1, 2), (1, 3), (3, 4)), (2, 0, 0, 0, 0)),
    ('1/133', ((1, 2), (1, 4), (3, 4)), (2, 0, 0, 0, 0)),
    ('1/133', ((1, 3), (1, 4), (2, 3)), (2, 0, 0, 0, 0)),
    ('8/133', ((1, 3), (1, 4), (2, 4)), (2, 0, 0, 0, 0)),
    ('51/133', ((1, 4), (2, 3), (2, 4)), (2, 0, 0, 0, 0)),
    ('15/133', ((1, 4), (2, 3), (3, 4)), (2, 0, 0, 0, 0)),
    ('22/133', ((0, 1), (0, 2), (1, 3), (3, 4)), (0, 0, 0, 0, 0)),
    ('36/133', ((0, 1), (0, 2), (1, 4), (3, 4)), (0, 0, 0, 0, 0)),
    ('43/133', ((0, 1), (0, 3), (1, 2), (3, 4)), (0, 0, 0, 0, 0)),
    ('6/133', ((0, 1), (0, 3), (2, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('17/133', ((0, 1), (0, 3), (2, 4), (3, 4)), (0, 0, 0, 0, 0)),
    ('18/133', ((0, 1), (0, 4), (1, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('58/133', ((0, 1), (0, 4), (2, 3), (3, 4)), (0, 0, 0, 0, 0)),
    ('58/133', ((0, 1), (1, 3), (2, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('8/133', ((0, 1), (1, 4), (2, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('5/19', ((0, 2), (0, 3), (1, 2), (1, 4)), (0, 0, 0, 0, 0)),
    ('32/133', ((0, 2), (0, 3), (1, 4), (2, 4)), (0, 0, 0, 0, 0)),
    ('29/133', ((0, 2), (0, 4), (1, 2), (3, 4)), (0, 0, 0, 0, 0)),
    ('8/133', ((0, 2), (0, 4), (1, 3), (2, 3)), (0, 0, 0, 0, 0)),
    ('43/133', ((0, 2), (0, 4), (1, 4), (2, 3)), (0, 0, 0, 0, 0)),
    ('10/133', ((0, 2), (1, 2), (1, 4), (3, 4)), (0, 0, 0, 0, 0)),
    ('16/133', ((0, 2), (1, 3), (1, 4), (2, 3)), (0, 0, 0, 0, 0)),
    ('5/19', ((0, 2), (1, 3), (2, 4), (3, 4)), (0, 0, 0, 0, 0)),
    ('5/133', ((0, 3), (0, 4), (1, 2), (1, 3)), (0, 0, 0, 0, 0)),
    ('32/133', ((0, 3), (0, 4), (1, 2), (1, 4)), (0, 0, 0, 0, 0)),
    ('15/133', ((0, 3), (0, 4), (1, 2), (2, 3)), (0, 0, 0, 0, 0)),
    ('46/133', ((0, 3), (1, 2), (1, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('12/133', ((0, 3), (1, 2), (1, 4), (2, 3)), (0, 0, 0, 0, 0)),
    ('23/133', ((0, 3), (1, 4), (2, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('23/133', ((0, 4), (1, 2), (1, 3), (2, 4)), (0, 0, 0, 0, 0)),
    ('16/133', ((0, 4), (1, 2), (1, 3), (3, 4)), (0, 0, 0, 0, 0)),
    ('1/7', ((0, 4), (1, 3), (1, 4), (2, 3)), (0, 0, 0, 0, 0)),
]

nv = 5
u = sp.symbols('u1:6')
E = [(i, j) for i in range(nv) for j in range(i + 1, nv)]
e1 = sum(u); p2 = sum(x**2 for x in u); e2 = sum(u[i] * u[j] for i, j in E)
e4 = sum(sp.prod([u[t] for t in S]) for S in itertools.combinations(range(nv), 4))
e5 = sp.prod(u)
N = sp.expand(5 + p2 - 2 * e2 + 3 * e4 - e1 * e5)
total = 0
for c, F, m in CERT5:
    c = sp.Rational(c)
    assert c > 0
    total += c * sp.prod([u[t]**m[t] for t in range(nv)]) * sp.prod([1 - u[i] * u[j] for i, j in F])
assert sp.expand(total - N) == 0
print(f"(3) exact identity with {len(CERT5)} positive rational terms: OK")

# (2) series coefficients vs direct Sp4 multiplicities
Npoly = sp.Poly(N, *u).as_dict()
def m0_free(k):  # loopless multigraphs on 5 vertices with degree sequence k
    if min(k) < 0 or sum(k) % 2:
        return 0
    deg = list(k)
    def rec(idx):
        if idx == len(E):
            return 1 if not any(deg) else 0
        i, j = E[idx]; tot = 0
        for x in range(min(deg[i], deg[j]) + 1):
            deg[i] -= x; deg[j] -= x
            if all(deg[v] == 0 or any(v in E[t] for t in range(idx + 1, len(E))) for v in (i, j)):
                tot += rec(idx + 1)
            deg[i] += x; deg[j] += x
        return tot
    return rec(0)
def series(k):
    return sum(int(c) * m0_free(tuple(x - y for x, y in zip(k, m))) for m, c in Npoly.items())
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
def T1_direct(k):
    P = defaultdict(int); P[(0, 0)] = 1
    for x in k:
        P = mul(P, hk(x))
    return 5 * mult(P, 0, 0) - 3 * mult(P, 1, 1) + mult(P, 2, 0)
bad = n = 0
for k in itertools.product(range(5), repeat=5):
    if sum(k) % 2 or sum(k) > 12:
        continue
    n += 1
    if series(k) != T1_direct(k):
        bad += 1
        print("MISMATCH", k)
print(f"(2) generating function vs direct Weyl-character T1: {n} compositions, {bad} mismatches")
assert bad == 0
print("Theorem T1-5 verified.")
