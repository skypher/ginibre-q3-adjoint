# FM-SEC64 (luna_max_saturn) printed code: strict-OL scans, band-boundary tuples, and the no-consecutive-flat resultant.
from math import comb
from itertools import product

def coeffs(a, e):
    c = [0] * (a + e + 1)
    for u in range(e + 1):
        for v in range(a + 1):
            c[u + v] += (-1)**u * comb(e, u) * comb(a, v)
    return c

def C(c, k):
    return c[k] if 0 <= k < len(c) else 0

def D(c, k):
    return C(c, k)**2 - C(c, k - 1) * C(c, k + 1)

steps, flat, zero_psi = 0, [], []
branch = {"min<=3": [0, 0], "high_p": [0, 0]}

for a, e in product(range(61), repeat=2):
    N, c = a + e, coeffs(a, e)
    for k in range(N // 2 + 1, N + 1):
        steps += 1
        delta = D(c, k) - D(c, k + 1)
        psi = (C(c, k), C(c, k - 1) + C(c, k + 1))
        if psi != (0, 0) and delta == 0:
            flat.append((a, e, k))
        if psi == (0, 0):
            zero_psi.append((a, e, k))
        if min(a, e) <= 3:
            branch["min<=3"][0] += 1
            branch["min<=3"][1] += delta == 0 and psi != (0, 0)
        p = 2 * k - N
        if 3 * p >= N - 2:
            branch["high_p"][0] += 1
            branch["high_p"][1] += delta == 0 and psi != (0, 0)

print("upper-half steps:", steps)
print("nonzero-psi zero deltas:", flat)
print("zero psi:", zero_psi)
print("branch [checked, nonzero-psi zero deltas]:", branch)

# ---- part ----
from math import comb, isqrt
from fractions import Fraction

boundary = []
for e in range(101):
    for N in range(2 * e, 1501):
        d, M = N - 2 * e, N + 3
        A = d * d - M * M - 1
        disc = A * A - 4 * M * M
        if disc < 0:
            continue
        q = isqrt(disc)
        if q * q != disc:
            continue
        for num in set((-A + q, -A - q)):
            if num % 2:
                continue
            y = num // 2
            u = isqrt(y) if y >= 0 else -1
            if u * u != y or u < 2:
                continue
            p = u - 1
            if (N + p) % 2 or 3 * p >= N - 2:
                continue
            k, n = (N + p) // 2, (N - p) // 2
            F = (p + 1)**2 * d*d - 4*p*(p + 2)*(n + 1)*(k + 2)
            assert F == 0
            a = N - e
            ck = sum((-1)**r * comb(e, r) * comb(a, k-r)
                     for r in range(e + 1) if 0 <= k-r <= a)
            ckp1 = sum((-1)**r * comb(e, r) * comb(a, k+1-r)
                       for r in range(e + 1) if 0 <= k+1-r <= a)
            actual = None if ck == 0 else Fraction(ckp1, ck)
            double_root = Fraction((p + 1) * d, 2 * p * (k + 2))
            boundary.append((e, N, p, double_root, actual,
                             actual == double_root))

print("F=0 tuples in box:", len(boundary))
for row in boundary:
    print(row)

# A continuing family of valid p=1 boundary tuples.
d, v = 45, 26
for _ in range(2):
    N = 2 * v - 3
    e = (N - d) // 2
    assert d*d == 3 * (v*v - 1) and (N - d) % 2 == 0
    print((e, N, 1))
    d, v = 7*d + 12*v, 4*d + 7*v

# ---- part ----
import sympy as s

p, n, d, r = s.symbols("p n d r")
k, N = n + p, 2*n + p
u = (d - (k + 1)*r) / (n + 1)
v = (d*r - n) / (k + 2)
w = (d*v - (n - 1)*r) / (k + 3)
b0, b1, b2 = u + r, 1 + v, r + w

P0 = s.factor((n + 1)*(k + 2)*(b1 - r*b0))
P1 = s.factor((k + 2)**2*(k + 3)*(r*b2 - v*b1))
Q = p*(k + 2)*r**2 - (p + 1)*d*r + (p + 2)*(n + 1)
resultant = s.factor(s.resultant(P0, P1, r))
target = ((k + 2)**2 * (p + 2)**2 * (N + 2 - d) *
          (N + 4 - d) * (N + 2 + d) * (N + 4 + d))

print("P0-Q:", s.factor(P0 - Q))
print("resultant-target:", s.factor(resultant - target))
print("resultant:", resultant)
