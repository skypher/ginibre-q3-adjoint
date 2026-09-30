"""FM-SEC116 (luna_max_mars): (1^N,n) under q: d = 2 sphere certificates for all n; d = 3 coordinates; N <= 20 census; (8,2,3) two-point repair."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb

parser = argparse.ArgumentParser(description='Exact q-Hermite radial formulas, N<=20 census, and a two-point repair.')
parser.parse_args()
ZERO, ONE = (Q(0),), (Q(1),)

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    return trim([
        (a[i] if i < len(a) else Q(0))
        + (b[i] if i < len(b) else Q(0))
        for i in range(max(len(a), len(b)))
    ])

def sub(a, b):
    return trim([
        (a[i] if i < len(a) else Q(0))
        - (b[i] if i < len(b) else Q(0))
        for i in range(max(len(a), len(b)))
    ])

def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)

def scale(a, c):
    return trim([Q(c) * x for x in a])

def qint(k):
    return (Q(1),) * k if k else ZERO

@lru_cache(None)
def power_expansion(r):
    cur = {0: ONE}
    for _ in range(r):
        nxt = {}
        for k, p in cur.items():
            nxt[k + 1] = add(nxt.get(k + 1, ZERO), p)
            if k:
                nxt[k - 1] = add(nxt.get(k - 1, ZERO), mul(p, qint(k)))
        cur = nxt
    return cur

def target_profile(n, d):
    N = n + 2 * d
    out = []
    for r in range(N + 1):
        if r % 2:
            out.append(ZERO)
        else:
            out.append(mul(
                power_expansion(r).get(0, ZERO),
                power_expansion(N - r).get(n, ZERO)
            ))
    return out

def triangular_coordinates(profile, n, d):
    b = []
    for t in range(d + 1):
        j = d - t
        v = scale(profile[2 * j], Q(1, comb(2 * j, j)))
        for h in range(t):
            v = sub(v, scale(b[h], comb(n + 2 * t, t - h)))
        b.append(v)
    return b

def sphere_coefficients(n, d):
    return triangular_coordinates(target_profile(n, d), n, d)

def qsumints(m):
    out = ZERO
    for k in range(1, m + 1):
        out = add(out, qint(k))
    return out

def M2_formula(m):
    out = ZERO
    for k in range(1, m + 1):
        for x in range(1, k + 2):
            out = add(out, mul(qint(x), qint(k)))
    return out

def M3_formula(m):
    out = ZERO
    for k in range(1, m + 1):
        for h in range(1, k + 2):
            for x in range(1, h + 2):
                out = add(out, mul(mul(qint(x), qint(h)), qint(k)))
    return out

C2 = (Q(2), Q(1))
C3 = (Q(5), Q(6), Q(3), Q(1))

# Check d=2 and d=3 triangular formulas against direct q-Hermite expansion.
for n in range(1, 11):
    m = n + 1
    S = qsumints(m)
    M2 = M2_formula(m)
    b = sphere_coefficients(n, 2)
    expected = [
        scale(C2, Q(1, 6)),
        sub(scale(S, Q(1, 2)), scale(C2, Q(m + 1, 6))),
        add(
            sub(M2, scale(S, Q(m + 3, 2))),
            scale(C2, Q(m * (m + 3), 12))
        )
    ]
    assert b == expected

for n in range(2, 21):
    m = n + 1
    S = qsumints(m)
    M2 = M2_formula(m)
    M3 = M3_formula(m)
    b = sphere_coefficients(n, 3)
    expected = [
        scale(C3, Q(1, 20)),
        sub(
            scale(mul(C2, S), Q(1, 6)),
            scale(C3, Q(m + 1, 20))
        ),
        add(
            sub(
                scale(M2, Q(1, 2)),
                scale(mul(C2, S), Q(m + 3, 6))
            ),
            scale(C3, Q(m * (m + 3), 40))
        ),
        add(
            add(
                add(M3, scale(M2, Q(-(m + 5), 2))),
                scale(mul(C2, S), Q((m + 5) * (m + 2), 12))
            ),
            scale(C3, Q(-(m + 5) * m * (m + 1), 120))
        )
    ]
    assert b == expected
    if n >= 3:
        low = (
            Q(m * (m + 3), 24),
            Q(m * (9 * m + 17), 60),
            Q(39 * m * m + 37 * m - 60, 120),
            Q(63 * m * m - 11 * m - 180, 120)
        )
        assert tuple(b[2][r] for r in range(4)) == low

# Exact sphere-only coefficient census, d<=n+1 and N<=20.
negative = []
cases = 0
for n in range(1, 21):
    for d in range(1, min(n + 1, (20 - n) // 2) + 1):
        cases += 1
        b = sphere_coefficients(n, d)
        bad = [
            (h, k, c)
            for h, p in enumerate(b)
            for k, c in enumerate(p)
            if c < 0
        ]
        if bad:
            negative.append((n + 2 * d, n, d, bad))

assert cases == 69
assert [(N, n, d) for N, n, d, _ in negative] == [
    (8, 2, 3), (11, 3, 4), (14, 4, 5), (17, 5, 6), (20, 6, 7)
]
assert all(all(h == 1 for h, _, _ in bad) for _, _, _, bad in negative)
print(
    'sphere census:', cases,
    'cases; coefficientwise failures:',
    [(N, n, d, bad[0][:2]) for N, n, d, bad in negative]
)

# Exact repair of (N,n,d)=(8,2,3) by all two-point factors with |v|=6.
def pair_orbit_coordinates(N, d):
    n = N - 2 * d
    profile = [ZERO] * (N + 1)
    profile[0] = (Q(2 * comb(N, 2 * d)),)
    profile[2 * d] = (Q(2),)
    return triangular_coordinates(profile, n, d)

N, n, d = 8, 2, 3
b = sphere_coefficients(n, d)
beta = pair_orbit_coordinates(N, d)
alpha = (Q(0), Q(1, 12), Q(0), Q(1, 12))
residual = [
    sub(bh, mul(alpha, (v[0],)))
    for bh, v in zip(b, beta)
]
assert residual == [
    (Q(1, 4), Q(7, 24), Q(3, 20), Q(1, 24)),
    (Q(0), Q(0), Q(1, 15)),
    (Q(3, 4), Q(17, 8), Q(67, 20), Q(23, 8), Q(3, 2), Q(1, 2)),
    (Q(8), Q(26), Q(764, 15), Q(57), Q(56), Q(41),
     Q(24), Q(11), Q(4), Q(1))
]
assert all(x >= 0 for p in residual for x in p)

def xor_square(support, S):
    support_set = set(support)
    return sum(1 for x in support if (x ^ S) in support_set)

masks = range(1 << N)
supports = {
    k: [x for x in masks if x.bit_count() == k]
    for k in range(d + 1)
}
target = target_profile(n, d)
for S in masks:
    value = ZERO
    for h in range(d + 1):
        value = add(
            value,
            scale(residual[h], xor_square(supports[d - h], S))
        )
    pair_sum = sum(
        xor_square([0, v], S)
        for v in masks
        if v.bit_count() == 2 * d
    )
    value = add(value, mul(alpha, (Q(pair_sum),)))
    assert value == target[S.bit_count()]

print(
    'two-point repair: coefficientwise nonnegative; exact re-expansion on',
    1 << N, 'masks PASS'
)