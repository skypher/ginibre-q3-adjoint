#!/usr/bin/env python3
"""Verifier for Theorem T1-4: FM3 at level r = 2 (T1) for every word with at most
four symmetric-power factors.

Claim.  For kappa in N^4 put M = Sym^k1 W (x) ... (x) Sym^k4 W, W = C^4, G = Sp(4), and
T1(kappa) = 5 m_(0,0) - 3 m_(1,1) + m_(2,0).  Then T1(kappa) >= 0.

Proof being checked.
 (1) T1 = 8 m0(kappa) + 4 m0(kappa,2) - 3 m0(kappa,1,1), with m0 the Sp4-invariant count,
     because W (x) W = 1 + U + Sym^2 W and Sym^2 W is self-dual.
 (2) By the first and second fundamental theorems for Sp4, the invariants of n vectors are
     C[omega_ij] with no relations for n <= 5 and C[omega_ij]/(Pf) for n = 6.  So
       sum_kappa m0(kappa) u^kappa          = Pi,
       sum_kappa m0(kappa,2) u^kappa        = h_2 Pi,
       sum_kappa m0(kappa,1,1) u^kappa      = (h_1^2 + 1 - e_4) Pi,
     with Pi = prod_{i<j<=4} (1 - u_i u_j)^(-1).  Hence
       sum_kappa T1(kappa) u^kappa = N Pi,   N = 5 + p_2 - 2 e_2 + 3 e_4.
 (3) N = sum of nine terms  u^m prod_{e in F} (1 - u_e)  with coefficient 1, so N Pi is a sum of
     monomials times products of geometric series: every coefficient is >= 0.
This script checks (3) exactly and checks (2) against a direct Weyl-character computation.
"""
import itertools
from collections import defaultdict
import sympy as sp

u = sp.symbols('u1:5')
E = [(i, j) for i in range(4) for j in range(i + 1, 4)]
N = 5 + sum(x**2 for x in u) - 2 * sum(u[i] * u[j] for i, j in E) + 3 * u[0] * u[1] * u[2] * u[3]
w = lambda i, j: 1 - u[i - 1] * u[j - 1]
TERMS = [
    w(1, 2),
    u[0]**2 * w(3, 4),
    u[2]**2 * w(1, 2) * w(2, 4),
    w(1, 3) * w(1, 4),
    u[3]**2 * w(1, 3) * w(2, 3),
    u[1]**2 * w(1, 4) * w(3, 4),
    w(1, 2) * w(2, 4) * w(3, 4),
    w(1, 3) * w(2, 3) * w(2, 4),
    w(1, 4) * w(2, 3) * w(3, 4),
]
assert sp.expand(sum(TERMS) - N) == 0
print("(3) exact identity N = sum of 9 manifestly positive terms: OK")

# (2): series coefficients vs direct Sp4 multiplicities (C2 Weyl character formula)
def m0_free(k):  # loopless multigraphs on 4 vertices with degree sequence k
    a, b, c, d = k
    if min(k) < 0 or (a + b + c + d) % 2:
        return 0
    T = (a + b + c - d) // 2
    n = 0
    for x12 in range(0, T + 1):
        for x13 in range(0, T - x12 + 1):
            x23 = T - x12 - x13
            if a - x12 - x13 >= 0 and b - x12 - x23 >= 0 and c - x13 - x23 >= 0:
                n += 1
    return n

Npoly = sp.Poly(N, *u).as_dict()
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
for k in itertools.product(range(7), repeat=4):
    if sum(k) % 2 or sum(k) > 14:
        continue
    n += 1
    if series(k) != T1_direct(k):
        bad += 1
        print("MISMATCH", k, series(k), T1_direct(k))
print(f"(2) generating function vs direct Weyl-character T1: {n} compositions, {bad} mismatches")
assert bad == 0
print("Theorem T1-4 verified.")
