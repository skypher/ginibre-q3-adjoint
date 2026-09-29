#!/usr/bin/env python3
"""Verifier: FM3 at level r = 1 for every word hat S_p hat S_q h_a h_b (all p, q, a, b >= 0).

G(t1,t2,u1,u2) = sum phi_1(hat S_p hat S_q h_a h_b) t1^p t2^q u1^a u2^b.
 (1) hat S_p = 2 h_p + 2 h_(p-2) - h_(p-1) h_1, so sum_p hat S_p t^p = (2 + 2t^2) H(t) - t H(t) h_1(z) with an extra
     degree-1 vertex z.  With the Sp4 invariant rings on <= 6 vertices (free for n <= 5, C[omega]/(Pf) for n = 6):
       G * prod_{edges of K_4}(1 - x_i x_j) = N,
       N = (2+2t1^2)(2+2t2^2) - t1 (2+2t2^2) h1 - t2 (2+2t1^2) h1 + t1 t2 (h1^2 + 1 - e4),
     h1 = t1 + t2 + u1 + u2, e4 = t1 t2 u1 u2.  Checked against exact phi_1 values (Weyl character formula).
 (2) N = sum of 10 terms x^m prod_{e in F}(1 - x_e) with coefficient 1, so G has nonnegative coefficients.
"""
import itertools
from collections import defaultdict
import sympy as sp
t1, t2, u1, u2 = X = sp.symbols('t1 t2 u1 u2')
h1 = t1 + t2 + u1 + u2; e4 = t1 * t2 * u1 * u2
N = sp.expand((2 + 2*t1**2) * (2 + 2*t2**2) - t1 * (2 + 2*t2**2) * h1 - t2 * (2 + 2*t1**2) * h1 + t1 * t2 * (h1**2 + 1 - e4))
w = lambda i, j: 1 - X[i] * X[j]
CERT = [X[1]**2 * w(0, 1), X[0]**2 * X[1]**2 * w(2, 3), w(0, 1) * w(1, 2), w(0, 1) * w(1, 3),
        X[0]**2 * X[1]**2 * w(0, 2) * w(2, 3), X[0]**2 * w(1, 3) * w(2, 3), w(0, 1) * w(0, 2) * w(0, 3),
        X[0]**2 * w(0, 1) * w(1, 2) * w(2, 3), X[1]**2 * w(0, 2) * w(0, 3) * w(2, 3),
        w(0, 2) * w(0, 3) * w(1, 2) * w(1, 3)]
assert sp.expand(sum(CERT) - N) == 0
print("(2) N = sum of 10 manifestly positive terms: OK")
# (1) against exact phi_1
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
def Sh(m):
    P = defaultdict(int)
    if m == 0:
        P[(0, 0)] = 2; return P
    for i in range(-m, m + 1, 2):
        P[(i, 0)] += 1; P[(0, i)] += 1
    return P
def mul(A, B):
    C = defaultdict(int)
    for (i, j), x in A.items():
        for (k, l), y in B.items():
            C[(i + k, j + l)] += x * y
    return C
def inv(P):
    v = (2, 1)
    return sum(det * P.get((sg[0] * v[perm[0]] - 2, sg[1] * v[perm[1]] - 1), 0) for perm, sg, det in WEYL)
D = 14
ser = {}
for c in itertools.product(range(D + 1), repeat=4):
    if sum(c) > D or sum(c) % 2:
        continue
    v = inv(mul(mul(Sh(c[0]), Sh(c[1])), mul(hk(c[2]), hk(c[3]))))
    if v:
        ser[c] = v
for (i, j) in itertools.combinations(range(4), 2):
    new = dict(ser)
    for c, v in ser.items():
        c2 = list(c); c2[i] += 1; c2[j] += 1
        if sum(c2) <= D:
            c2 = tuple(c2); new[c2] = new.get(c2, 0) - v
    ser = {k: v for k, v in new.items() if v}
target = {k: int(v) for k, v in sp.Poly(N, *X).as_dict().items()}
assert ser == target, "numerator mismatch"
print(f"(1) closed-form numerator matches exact phi_1 series through degree {D}")
print("Level-1 two-label, two-factor sector verified.")
