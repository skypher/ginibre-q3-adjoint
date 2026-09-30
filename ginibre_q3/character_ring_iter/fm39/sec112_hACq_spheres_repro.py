"""FM-SEC112 (luna_max_mars): q-deformed sphere certificates for (1^N,n): closed forms, (1^8,2) failure of sphere-only factors and exact repair."""
from functools import lru_cache
from fractions import Fraction as Q
from math import comb

ZERO, ONE = (Q(0),), (Q(1),)

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])

def sub(a, b):
    return trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])

def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)

def scale(a, c):
    return trim([c * x for x in a])

def monomial(k, c):
    return (Q(0),) * k + (Q(c),)

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

def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def divide_by_one_minus_q(p):
    out, total = [], Q(0)
    for x in p:
        total += x
        out.append(total)
    assert total == 0
    return trim(out)

def touchard_riordan(j):
    p = ZERO
    for k in range(j + 1):
        m = binom(2*j, j-k) - binom(2*j, j-k-1)
        p = add(p, monomial(k*(k+1)//2, (-1)**k * m))
    for _ in range(j):
        p = divide_by_one_minus_q(p)
    return p

def target_profile(n, d):
    N = n + 2*d
    out = []
    for s in range(N + 1):
        if s % 2:
            out.append(ZERO)
        else:
            left = power_expansion(s).get(0, ZERO)
            right = power_expansion(N-s).get(n, ZERO)
            out.append(mul(left, right))
    return out

def sphere_coefficients(n, d):
    f = target_profile(n, d)
    b = []
    for t in range(d + 1):
        j = d - t
        v = scale(f[2*j], Q(1, comb(2*j, j)))
        for h in range(t):
            v = sub(v, scale(b[h], comb(n + 2*t, t-h)))
        b.append(v)
    return b

def pair_orbit_coefficients(N, d):
    n = N - 2*d
    f = [ZERO] * (N + 1)
    f[0] = (Q(2 * comb(N, 2*d)),)
    f[2*d] = (Q(2),)
    b = []
    for t in range(d + 1):
        j = d - t
        v = scale(f[2*j], Q(1, comb(2*j, j)))
        for h in range(t):
            v = sub(v, scale(b[h], comb(n + 2*t, t-h)))
        b.append(v)
    return b

for j in range(9):
    assert touchard_riordan(j) == power_expansion(2*j).get(0, ZERO)

first = None
for n in range(1, 9):
    for d in range(n + 2):
        if n + 2*d > 8:
            continue
        for h, p in enumerate(sphere_coefficients(n, d)):
            for k, c in enumerate(p):
                if c < 0:
                    candidate = (n + 2*d, n, d, h, k, c)
                    if first is None or candidate[:3] < first[:3]:
                        first = candidate
assert first == (8, 2, 3, 1, 1, Q(-1, 30))

higher_index_negative = None
for n in range(1, 13):
    for d in range(n + 2):
        for h, p in enumerate(sphere_coefficients(n, d)):
            if h >= 2 and any(c < 0 for c in p):
                higher_index_negative = (n, d, h)
                break
assert higher_index_negative is None

N, n, d = 8, 2, 3
b = sphere_coefficients(n, d)
assert b == [
    (Q(1,4), Q(3,10), Q(3,20), Q(1,20)),
    (Q(0), Q(-1,30), Q(1,15), Q(-1,30)),
    (Q(3,4), Q(11,5), Q(67,20), Q(59,20), Q(3,2), Q(1,2)),
    (Q(8), Q(458,15), Q(764,15), Q(923,15),
     Q(56), Q(41), Q(24), Q(11), Q(4), Q(1)),
]
v = pair_orbit_coefficients(N, d)
assert tuple(p[0] for p in v) == (
    Q(1,10), Q(-2,5), Q(9,10), Q(272,5)
)

alpha = (Q(0), Q(1,12), Q(0), Q(1,12))
residual = [sub(bh, scale(alpha, vh[0])) for bh, vh in zip(b, v)]
assert residual == [
    (Q(1,4), Q(7,24), Q(3,20), Q(1,24)),
    (Q(0), Q(0), Q(1,15)),
    (Q(3,4), Q(17,8), Q(67,20), Q(23,8), Q(3,2), Q(1,2)),
    (Q(8), Q(26), Q(764,15), Q(57), Q(56), Q(41),
     Q(24), Q(11), Q(4), Q(1)),
]
assert all(c >= 0 for p in residual for c in p)

def convolution(support, S):
    support_set = set(support)
    return sum(1 for x in support if (x ^ S) in support_set)

masks = range(1 << N)
sphere_conv = []
for k in range(d + 1):
    support = [x for x in masks if x.bit_count() == k]
    sphere_conv.append([convolution(support, S) for S in masks])

weight6 = [x for x in masks if x.bit_count() == 2*d]
pair_conv = [[convolution([0, v0], S) for S in masks] for v0 in weight6]
f = target_profile(n, d)
for S in masks:
    rhs = ZERO
    for h in range(d + 1):
        rhs = add(rhs, scale(residual[h], sphere_conv[d-h][S]))
    for pc in pair_conv:
        rhs = add(rhs, mul(alpha, (Q(pc[S]),)))
    assert rhs == f[S.bit_count()]

print("Touchard-Riordan / q-Hermite check j=0..8: PASS")
print("first negative sphere coefficient by N<=8:", first)
print("h>=2 coefficients nonnegative in scan n<=12, d<=n+1:",
      higher_index_negative is None)
print("b_h at (N,n,d)=(8,2,3):", b)
print("two-point orbit coordinates:", tuple(p[0] for p in v))
print("residual polynomials:", residual)
print("exact re-expansion on all 256 masks: PASS")
