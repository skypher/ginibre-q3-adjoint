"""FM-SEC110 (luna_max_mercury): H_AC_q at distance one (uniform q-certificate), sphere factor failure at (1^8,2), distance-two tests."""
import argparse
from functools import lru_cache
from fractions import Fraction as Q
from itertools import combinations, product, permutations
from math import comb

parser = argparse.ArgumentParser(description="Exact H_AC_q family checks")
parser.parse_args()

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])

def scale(a, c):
    return trim([c*x for x in a])

def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)

def shift(a, k):
    return (0,)*k + a

@lru_cache(None)
def qbinom(n, k):
    if k < 0 or k > n:
        return (0,)
    if k == 0 or k == n:
        return (1,)
    return add(qbinom(n-1, k),
               shift(qbinom(n-1, k-1), n-k))

@lru_cache(None)
def qfactorial(n):
    out = (1,)
    for j in range(1, n+1):
        out = mul(out, (1,)*j)       # [j]_q
    return out

@lru_cache(None)
def linearization(a, b, k):
    return mul(mul(qbinom(a, k), qbinom(b, k)), qfactorial(k))

@lru_cache(None)
def moment(labels):
    # Coefficient of H_0 using exact q-Hermite linearization.
    state = {0: (1,)}
    for n in labels:
        nxt = {}
        for d, p in state.items():
            for k in range(min(d, n)+1):
                e = d+n-2*k
                nxt[e] = add(nxt.get(e, (0,)),
                             mul(p, linearization(d, n, k)))
        state = nxt
    return state.get(0, (0,))

def table(mu, n):
    # Coordinates are mu; n is the final inserted label.
    out = []
    for S in range(1 << len(mu)):
        left = tuple(mu[i] for i in range(len(mu)) if S >> i & 1)
        right = tuple(mu[i] for i in range(len(mu))
                      if not (S >> i & 1)) + (n,)
        out.append(mul(moment(left), moment(right)))
    return out

def xor_square(p):
    out = [Q(0)] * len(p)
    for x, a in enumerate(p):
        for y, b in enumerate(p):
            out[x ^ y] += a*b
    return out

def a_mu(mu):
    M = sum(mu)
    return trim([M-d-sum(max(a-d, 0) for a in mu)
                 for d in range(1, M)])

def endpoint_gap_poly(mu):
    blocks = []
    for i, size in enumerate(mu):
        blocks.extend([i]*size)
    out = [0]*max(1, len(blocks)-1)
    for a in range(len(blocks)):
        for b in range(a+1, len(blocks)):
            if blocks[a] != blocks[b]:
                out[b-a-1] += 1
    return trim(out)

def d1_certificate(mu):
    n = sum(mu)-2
    f = table(mu, n)
    t = mu.count(1)
    p = [Q(0)]*(1 << len(mu))
    for i, label in enumerate(mu):
        if label == 1:
            p[1 << i] = Q(1)
    pp = xor_square(p)
    fac = qfactorial(n)
    delta_coefficient = mul(fac, add(a_mu(mu), (Q(-t, 2),)))
    assert all(c >= 0 for c in fac + delta_coefficient)
    rhs = []
    for x in range(len(f)):
        v = scale(fac, Q(pp[x], 2))
        if x == 0:
            v = add(v, delta_coefficient)
        rhs.append(v)
    assert f == rhs
    assert a_mu(mu) == endpoint_gap_poly(mu)
    return f

# Census driver: all ordered mu with entries <=4 and length <=5.
d1_count = 0
for m in range(1, 6):
    for mu in product(range(1, 5), repeat=m):
        if sum(mu) > 2:
            d1_certificate(mu)
            d1_count += 1
    print("d1 length", m, "completed", flush=True)
print("d1 exact census:", d1_count, "ordered lists PASS")

def sphere_coefficients(N, n):
    k = (N-n)//2
    values = [mul(moment((1,)*(2*j)),
                  moment((1,)*(N-2*j)+(n,))) for j in range(k+1)]
    coeff = [(Q(0),)]*(k+1)
    for j in range(k, -1, -1):
        v = values[j]
        for ell in range(j+1, k+1):
            factor = comb(2*j, j)*comb(N-2*j, ell-j)
            v = add(v, scale(coeff[ell], -factor))
        coeff[j] = scale(v, Q(1, comb(2*j, j)))
    return coeff

def sphere_square(N, ell):
    p = [Q(int(mask.bit_count() == ell))
         for mask in range(1 << N)]
    return xor_square(p)

# Find the first coefficientwise failure of the original sphere factors.
first = None
sphere_cases = 0
for N in range(1, 13):
    for n in range(1, N+1):
        if (N+n) % 2 or 3*n < N-2:
            continue
        sphere_cases += 1
        coeff = sphere_coefficients(N, n)
        if any(z < 0 for p in coeff for z in p):
            first = (N, n, coeff)
            break
    if first:
        break

assert first is not None
N, n, sphere_c = first
assert (N, n) == (8, 2)
assert sphere_c[2] == (Q(0), -Q(1,30), Q(1,30),
                       Q(1,30), -Q(1,30))
assert sum(c * Q(1,2)**i for i, c in enumerate(sphere_c[2])) \
       == -Q(1,160)

full = table((1,)*N, n)
reexp = [(0,)]*(1 << N)
for ell, c in enumerate(sphere_c):
    atom = sphere_square(N, ell)
    reexp = [add(v, scale(c, atom[x])) for x, v in enumerate(reexp)]
assert reexp == full
print("sphere census:", sphere_cases, "cases to first failure")
print("sphere radius-2 coefficient:", sphere_c[2])
print("sphere q=0 coefficients by radius:",
      tuple(p[0] for p in sphere_c))

# Boundary (1^8,6).
mu = (1,)*8
f86 = d1_certificate(mu)
fac6 = qfactorial(6)
A8 = a_mu(mu)
a86 = mul(fac6, add(A8, (Q(-4),)))
b86 = scale(fac6, Q(1,2))
p1 = [Q(int(x.bit_count() == 1)) for x in range(256)]
p7 = [Q(int(x.bit_count() == 7)) for x in range(256)]
ss1, ss7 = xor_square(p1), xor_square(p7)
assert ss1 == ss7
assert f86 == [add(scale(b86, ss1[x]), a86 if x == 0 else (0,))
               for x in range(256)]
print("(1^8,6): A(q) =", A8)
print("(1^8,6): b(q) =", b86)
print("(1^8,6): a(q) =", a86)

# Boundary (1,5,2,2), with 5 distinguished.
f1522 = table((1,2,2), 5)
assert f1522 == [qfactorial(5)] + [(0,)]*7
print("(1,5,2,2): [5]_q! =", qfactorial(5))

# d=2 passing example: mu=(1,2,2), n=1.
mu, n = (1,2,2), 1
f = table(mu, n)
B = [Q(0)]*(1 << len(mu))
for i, label in enumerate(mu):
    if label == 2:
        B[1 << i] = Q(1)
BB = xor_square(B)
c = scale(qfactorial(2), Q(1,2))
rho = add(f[0], scale(c, -BB[0]))
rhs = [add(scale(c, BB[x]), rho if x == 0 else (0,))
       for x in range(len(f))]
assert rhs == f and all(z >= 0 for z in c+rho)
print("d2 mu=(1,2,2): B-square =", c, "delta =", rho)

# d=2 passing example: mu=1^5, n=1.
mu, n = (1,)*5, 1
f = table(mu, n)
A = [Q(0)]*(1 << len(mu))
for i, j in combinations(range(5), 2):
    A[(1 << i) | (1 << j)] = Q(1)
AA = xor_square(A)
c = scale(moment((1,)*4), Q(1,6))
rho = add(f[0], scale(c, -AA[0]))
rhs = [add(scale(c, AA[x]), rho if x == 0 else (0,))
       for x in range(len(f))]
assert rhs == f and all(z >= 0 for z in c+rho)
print("d2 mu=1^5: A-square =", c, "delta =", rho)

# d=2 fixed-factor obstruction: mu=(1^4,2), n=2.
mu, n = (1,)*4 + (2,), 2
f = table(mu, n)
mask1111, mask112 = 15, 19
assert f[mask1111] == (2,3,1) and f[mask112] == (1,2,1)

A = [Q(0)]*32
E = [Q(0)]*32
B = [Q(0)]*32
for i in range(4):
    E[1 << i] = Q(1)
for i, j in combinations(range(4), 2):
    A[(1 << i) | (1 << j)] = Q(1)
B[1 << 4] = Q(1)
P = A[:]
P[1 << 4] = Q(3,2)
D = [Q(int(x == 0)) for x in range(32)]

# q=0 weights from the FM-MECH34 factors.
atoms = [(Q(1,3), P), (Q(1,2), B), (Q(1,3), E),
         (Q(17,12), D)]
base = [(0,)]*32
for weight, p in atoms:
    sq = xor_square(p)
    base = [add(v, scale((z,), weight)) for v, z in zip(base, sq)]
assert [p[0] for p in base] == [p[0] for p in f]
assert f[mask1111] != scale(f[mask112], Q(2))
print("d2 fixed-factor obstruction:",
      "f(1111) =", f[mask1111], "f(112) =", f[mask112])

# Direct inversion-enumerator check for [n]_q!.
def inversion_poly(n):
    out = [0]*(n*(n-1)//2 + 1)
    for pi in permutations(range(n)):
        inv = sum(pi[i] > pi[j]
                  for i in range(n) for j in range(i+1, n))
        out[inv] += 1
    return trim(out)

for n in range(1, 7):
    assert inversion_poly(n) == qfactorial(n)
print("Mahonian check: [n]_q! enumerates inversions for n<=6")
