"""FM-SEC119 (luna_max_mercury): d = 3 table, append recurrence, radial E/A/H obstruction at (1^7), n = 1."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb

parser = argparse.ArgumentParser(
    description="Exact d=3 all-one radial-factor census."
)
args = parser.parse_args()
ZERO = (F(0),)

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    return trim([
        (a[i] if i < len(a) else F(0))
        + (b[i] if i < len(b) else F(0))
        for i in range(max(len(a), len(b)))
    ])

def scale(p, c):
    return trim([F(c)*x for x in p])

def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)

def qint(n):
    return (F(1),)*n if n > 0 else ZERO

@lru_cache(None)
def qbinom(n, k):
    if k < 0 or k > n:
        return ZERO
    if k == 0 or k == n:
        return (F(1),)
    return add(qbinom(n-1, k), (F(0),)*(n-k)+qbinom(n-1, k-1))

@lru_cache(None)
def qfactorial(n):
    p = (F(1),)
    for j in range(1, n+1):
        p = mul(p, qint(j))
    return p

@lru_cache(None)
def linearization(a, b, k):
    return mul(mul(qbinom(a, k), qbinom(b, k)), qfactorial(k))

@lru_cache(None)
def row(labels):
    current = {0: (F(1),)}
    for b in labels:
        nxt = {}
        for degree, p in current.items():
            for k in range(min(degree, b)+1):
                new_degree = degree+b-2*k
                term = mul(p, linearization(degree, b, k))
                nxt[new_degree] = add(nxt.get(new_degree, ZERO), term)
        current = nxt
    return current

def moment(labels):
    return row(tuple(sorted(labels))).get(0, ZERO)

def radial_profile(n, d=3):
    N = n+2*d
    out = {}
    for r in range(N+1):
        if r % 2:
            out[r] = ZERO
        else:
            out[r] = mul(
                moment((1,)*r),
                row((1,)*(N-r)).get(n, ZERO)
            )
    return out

def sphere_coefficients(profile, N, d=3):
    # For the j-subset sphere p_j:
    # (p_j*p_j)(S)=C(2k,k)C(N-2k,j-k) when |S|=2k.
    b = [ZERO]*(d+1)
    for k in range(d, -1, -1):
        value = scale(profile[2*k], F(1, comb(2*k, k)))
        for j in range(k+1, d+1):
            value = add(value, scale(b[j], -comb(N-2*k, j-k)))
        b[k] = value
    return b

patterns = {
    (1, 1, 1, 1, 1, 1): (5, 6, 3, 1),
    (1, 1, 1, 1, 2): (3, 5, 3, 1),
    (1, 1, 2, 2): (2, 4, 3, 1),
    (2, 2, 2): (1, 3, 3, 1),
    (1, 1, 1, 3): (1, 2, 2, 1),
    (1, 2, 3): (1, 2, 2, 1),
    (3, 3): (1, 2, 2, 1),
}
for labels, expected in patterns.items():
    value = moment(labels)
    assert value == tuple(F(x) for x in expected), (labels, value)
    print("low-shell moment", labels, "=", value, flush=True)

first_failure = None
for N in range(7, 15):
    n = N-6
    profile = radial_profile(n)
    b = sphere_coefficients(profile, N)
    required = F(3*(N-4), 10)*profile[6][0]
    gap = profile[4][0]-required
    print(
        "all-one census", "N=", N, "n=", n,
        "g4(q=0)=", profile[4][0], "g6(q=0)=", profile[6][0],
        "radial lower bound=", required, "gap=", gap,
        "A-sphere coefficient=", b[2],
        flush=True
    )
    if gap < 0 and first_failure is None:
        first_failure = (N, n, profile[4][0], profile[6][0], gap, b[2])

assert first_failure == (
    7, 1, F(4), F(5), F(-1, 2),
    (F(-1, 12), F(-7, 30), F(-17, 60), F(-3, 20))
), first_failure
print("first radial-family failure:", first_failure, flush=True)
