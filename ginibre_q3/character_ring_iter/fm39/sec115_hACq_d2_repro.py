"""FM-SEC115 (luna_max_mercury): independent proof and census of H_AC_q at distance two."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    return trim([
        (a[i] if i < len(a) else 0)
        + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ])

def scale(p, c):
    return trim([c*x for x in p])

def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)

def sub(a, b):
    return add(a, scale(b, -1))

def qi(n):
    return (F(1),) * n if n > 0 else (F(0),)

@lru_cache(None)
def qbinom(n, k):
    if k < 0 or k > n:
        return (F(0),)
    if k == 0 or k == n:
        return (F(1),)
    return add(
        qbinom(n-1, k),
        (F(0),) * (n-k) + qbinom(n-1, k-1)
    )

@lru_cache(None)
def qfactorial(n):
    p = (F(1),)
    for j in range(1, n+1):
        p = mul(p, qi(j))
    return p

@lru_cache(None)
def linearization(a, b, k):
    return mul(mul(qbinom(a, k), qbinom(b, k)), qfactorial(k))

@lru_cache(None)
def moment_sorted(labels):
    # Coefficient of H_0 in a product of continuous q-Hermite polynomials.
    current = {0: (F(1),)}
    for n in labels:
        nxt = {}
        for d, p in current.items():
            for k in range(min(d, n) + 1):
                degree = d + n - 2*k
                term = mul(p, linearization(d, n, k))
                nxt[degree] = add(nxt.get(degree, (F(0),)), term)
        current = nxt
    return current.get(0, (F(0),))

def moment(labels):
    return moment_sorted(tuple(sorted(labels)))

def one_chord(mu):
    total = 0
    out = (F(0),)
    for a in mu:
        out = add(out, mul(qi(total), qi(a)))
        total += a
    return out

def two_chord(mu):
    prefix = []
    total = 0
    for a in mu:
        prefix.append(total)
        total += a
    out = (F(0),)
    for s, a in zip(prefix, mu):
        term = mul(mul(qbinom(s, 2), qbinom(a, 2)), qfactorial(2))
        out = add(out, term)
    for i, j in combinations(range(len(mu)), 2):
        term = (F(1),)
        for x in (prefix[i], mu[i], prefix[j]-2, mu[j]):
            term = mul(term, qi(x))
        out = add(out, term)
    return out

def map_add(a, b):
    out = dict(a)
    for x, p in b.items():
        out[x] = add(out.get(x, (F(0),)), p)
    return out

def map_scale(a, c):
    return {x: mul(p, c) for x, p in a.items()}

def square(p):
    out = {}
    for x, a in p.items():
        for y, b in p.items():
            out[x ^ y] = add(
                out.get(x ^ y, (F(0),)),
                mul(a, b)
            )
    return out

def point(x, c=F(1)):
    return {x: (c,)}

def delete_two_ones(mu):
    out = []
    removed = 0
    for a in mu:
        if a == 1 and removed < 2:
            removed += 1
        else:
            out.append(a)
    return tuple(out)

def certificate(mu):
    L = len(mu)
    t = mu.count(1)
    v = mu.count(2)
    w = L - t - v

    ones = [1 << i for i, a in enumerate(mu) if a == 1]
    twos = [1 << i for i, a in enumerate(mu) if a == 2]
    E = {x: (F(1),) for x in ones}
    A = {i ^ j: (F(1),) for i, j in combinations(ones, 2)}
    B = {x: (F(1),) for x in twos}
    out = {}

    def take(c, p):
        nonlocal out
        assert all(x >= 0 for x in c)
        out = map_add(out, map_scale(square(p), c))

    if t <= 1:
        take((F(1, 2), F(1, 2)), B)
    elif t <= 3:
        take((F(1, 2), F(1, 2)), map_add(A, B))
        alpha = one_chord(delete_two_ones(mu))
        e = scale(
            sub(alpha, scale((F(1), F(1)), F(t-2))),
            F(1, 2)
        )
        take(e, E)
    else:
        alpha = one_chord(delete_two_ones(mu))
        e = scale(
            sub(alpha, scale((F(2), F(1)), F(t-2, 3))),
            F(1, 2)
        )
        if v == 0:
            take((F(1, 3), F(1, 6)), A)
        elif v <= 3:
            for z in twos:
                P = map_add(A, point(z, F(3*v, 2)))
                Q = map_add(A, point(z, F(3*v)))
                take((F(1, 3*v),), P)
                take((F(0), F(1, 6*v)), Q)
            take((F(1, 2), F(1, 2)), B)
        else:
            for j, k in combinations(twos, 2):
                R = map_add(
                    A,
                    map_add(point(j, F(v-1)), point(k, F(v-1)))
                )
                take((F(v, 4*(v-1)*comb(v, 2)),), R)
            take((F(v-4, 12*(v-1)),), A)
            for z in twos:
                Q = map_add(A, point(z, F(3*v)))
                take((F(0), F(1, 6*v)), Q)
            take((F(0), F(1, 2)), B)
        take(e, E)

    rho = sub(two_chord(mu), out.get(0, (F(0),)))
    assert all(c >= 0 for c in rho), (mu, rho)
    out = map_add(out, {0: rho})
    return out, rho

def ordered_lists(total, prefix=()):
    if total == 0:
        yield prefix
        return
    for a in range(1, min(4, total) + 1):
        yield from ordered_lists(total-a, prefix + (a,))

count = 0
entries = 0
min_residual = None
max_residual_degree = 0

for total in range(5, 13):
    for mu in ordered_lists(total):
        n = total - 4
        g, rho = certificate(mu)
        nfact = qfactorial(n)
        for mask in range(1 << len(mu)):
            left = tuple(
                sorted(mu[i] for i in range(len(mu)) if mask >> i & 1)
            )
            right = (n,) + tuple(
                sorted(mu[i] for i in range(len(mu)) if not (mask >> i & 1))
            )
            expected = mul(moment(left), moment(right))
            actual = mul(g.get(mask, (F(0),)), nfact)
            assert actual == expected, (mu, mask, actual, expected)
            entries += 1
        count += 1
        min_residual = min(min_residual, min(rho)) if min_residual is not None else min(rho)
        max_residual_degree = max(max_residual_degree, len(rho)-1)

samples = (
    (1, 2, 2),
    (1, 1, 1, 1, 1),
    (1, 1, 1, 1, 2),
    (1, 1, 1, 1, 1, 1, 1, 2, 2),
)
for mu in samples:
    alpha = one_chord(delete_two_ones(mu)) if mu.count(1) >= 2 else None
    print("sample", mu, "D=", two_chord(mu), "alpha=", alpha,
          "rho=", certificate(mu)[1])

assert (count, entries, min_residual, max_residual_degree) == (
    3080, 496584, F(0), 17
)
print("PASS:", count, "ordered lists;", entries, "subset entries")
print("minimum residual coefficient:", min_residual)
print("maximum residual degree:", max_residual_degree)