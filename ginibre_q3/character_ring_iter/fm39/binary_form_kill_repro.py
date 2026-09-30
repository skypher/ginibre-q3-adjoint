
import argparse
from fractions import Fraction as Q
from math import comb
from functools import lru_cache

parser = argparse.ArgumentParser(
    description="FM-MECH27 counterexamples and endpoint identity.")
parser.parse_args()

def coeff(a, e):
    N = a + e
    return [
        sum((-1)**b * comb(e, b) * comb(a, k-b)
            for b in range(max(0, k-a), min(e, k)+1))
        for k in range(N+1)
    ]

def extend(N, d, x, p, q):
    y = {x: Q(p), x+1: Q(q)}
    for k in range(x+1, N):
        y[k+1] = (d*y[k]-(N-k+1)*y[k-1])/(k+1)
    for k in range(x, 0, -1):
        y[k-1] = (d*y[k]-(k+1)*y[k+1])/(N-k+1)
    return [y[k] for k in range(N+1)]

def Phi(y, A, B, C):
    N = len(y)-1
    l = (N+A-B-C)//2
    t = l+B+1
    def z(k):
        return y[k] if 0 <= k <= N else 0
    def D(k):
        return z(k)**2-z(k-1)*z(k+1)
    def P(k):
        return (sum(D(j) for j in range(k, k+C+1))
                -z(k)*z(k+C)+z(k-1)*z(k+C+1))
    def ac(k):
        return z(k)-z(k+C)
    def bc(k):
        return z(k-1)-z(k+C+1)
    return P(l)-P(t)-ac(l)*bc(t)+bc(l)*ac(t)

@lru_cache(None)
def ballot(n, p):
    if n < p or (n-p) % 2:
        return 0
    return (p+1)*comb(n, (n-p)//2)//((n+p)//2+1)

def direct(r, a, u, v, w):
    e = 2*r-3
    N = a+e
    c = coeff(a, e)
    A, B, C = u+1, v+1, w+1
    def cg(p, q):
        return range(abs(p-q), p+q+1, 2)
    def W(p, q):
        return sum(c[k]*ballot(N-k, p)*ballot(k, q)
                   for k in range(N+1))
    return (
        sum(W(l, 0) for d in cg(A, B) for l in cg(d, C))
        -sum(W(d, C) for d in cg(A, B))
        -sum(W(d, B) for d in cg(A, C))
        -sum(W(d, A) for d in cg(B, C))
    )

checks = indefinite = 0
for r in range(2, 10):
    e = 2*r-3
    for a in range(31):
        N, d = a+e, a-e
        if N < 2:
            continue
        A, B, C = N+4, 3, 3
        for p, q in [(1, 0), (0, 1), (1, 1), (d, 4)]:
            y = extend(N, d, N-1, p, q)
            expected = Q(p*p)-Q(d*p*q, 2)+Q((N+2)*q*q, 2)
            assert Phi(y, A, B, C) == expected
        expected = Q(d*d+N+2, 2)
        assert Phi(coeff(a, e), A, B, C) == expected
        assert direct(r, a, N+3, 2, 2) == expected
        checks += 1
        indefinite += d*d > 8*(N+2)
print("edge checks:", checks, "indefinite:", indefinite)

for r, a, u, v, w, x, p, q in [
    (2, 12, 16, 2, 2, 12, 11, 4),
    (7, 0, 14, 2, 2, 10, -11, 4),
    (2, 36, 26, 2, 2, 29, 3, 2),
]:
    e = 2*r-3
    N, d = a+e, a-e
    A, B, C = u+1, v+1, w+1
    y = extend(N, d, x, p, q)
    assert all((k+1)*y[k+1] == d*y[k]-(N-k+1)*y[k-1]
               for k in range(1, N))
    negative = Phi(y, A, B, C)
    actual = direct(r, a, u, v, w)
    assert negative < 0
    assert actual == Phi(coeff(a, e), A, B, C) > 0
    print((r, a, u, v, w), negative, actual)

values = [
    Phi(extend(37, 35, 29, *pq), 27, 3, 3)
    for pq in [(1, 0), (0, 1), (1, 1)]
]
f10, f01, f11 = values
scale = 658003434303600
print("interior numerator:",
      tuple(scale*z for z in (f10, f11-f10-f01, f01)))

y = extend(14, 12, 11, 5, 3)
assert y[11]**2-y[10]*y[12] == 7
assert Phi(y, 14, 3, 3) == -Q(563, 8281)
print("single Newton:", Phi(y, 14, 3, 3), "D_alpha:", 7)

for r, a in [(2, 0), (2, 2), (2, 4),
             (3, 0), (3, 2), (3, 4)]:
    e = 2*r-3
    N = a+e
    x = min(N-1, N//2)
    f10, f01, f11 = [
        Phi(extend(N, a-e, x, *pq), 3, 3, 3)
        for pq in [(1, 0), (0, 1), (1, 1)]
    ]
    assert f10 > 0 and 4*f10*f01 > (f11-f10-f01)**2
    actual = direct(r, a, 2, 2, 2)
    assert actual == Phi(coeff(a, e), 3, 3, 3)
    print("h2^3:", (r, a), "PD", actual)
