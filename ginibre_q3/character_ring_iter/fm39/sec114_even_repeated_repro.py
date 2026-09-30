"""FM-SEC114 (luna_max_uranus): repeated even labels (2^M): profiles, rank-one radial obstruction at M=5, B certificate at M=5."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb

# q=2: mu_n = [U_0] U_2^n, computed by exact SU(2) fusion.
def moments(nmax):
    row = {0: 1}
    out = [1]
    for _ in range(nmax):
        new = {}
        for j, multiplicity in row.items():
            for ell in range(abs(j - 2), j + 3, 2):
                new[ell] = new.get(ell, 0) + multiplicity
        row = new
        out.append(row.get(0, 0))
    return out

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def kraw(M, t, k):
    return sum((-1)**j * choose(t, j) * choose(M-t, k-j)
               for j in range(k+1))

def dual_kraw(M, t, k):
    return sum((-1)**j * choose(k, j) * choose(M-k, t-j)
               for j in range(t+1))

def radial_data(M, mu):
    values = [mu[k] * mu[M-k] for k in range(M+1)]
    # Characters on G are the even-weight T; the quotient Fourier sum
    # is half the full-cube sum.
    spectrum = {}
    for t in range(0, M+1, 2):
        z = Q(sum(kraw(M, t, k) * values[k]
                  for k in range(M+1)), 2)
        assert z.denominator == 1
        spectrum[t] = z.numerator
    return values, spectrum

mu = moments(11)
assert mu == [1, 0, 1, 1, 3, 6, 15, 36, 91, 232, 603, 1585]
for M in range(2, 6):
    vals, spectrum = radial_data(M, mu)
    print("M", M, "f-profile", vals[:M//2+1], "Fourier", spectrum)

# Exact principal-root profiles for M=2,3,4.
# Represent A+B*sqrt(2)+C*sqrt(6).
def add3(x, y):
    return tuple(a+b for a, b in zip(x, y))

def scale3(c, x):
    return tuple(c*a for a in x)

def root3(n):
    table = {
        1: (Q(1), Q(0), Q(0)),
        2: (Q(0), Q(1), Q(0)),
        6: (Q(0), Q(0), Q(1)),
    }
    return table[n]

def root_profile(M, spectrum):
    den = 1 << (M-1)
    ans = []
    for k in range(M//2+1):
        z = (Q(0), Q(0), Q(0))
        for t, h in spectrum.items():
            z = add3(z, scale3(dual_kraw(M, t, k), root3(h)))
        ans.append(scale3(Q(1, den), z))
    return ans

assert root_profile(2, radial_data(2, mu)[1]) == [
    (Q(1), Q(0), Q(0)), (Q(0), Q(0), Q(0))]
assert root_profile(3, radial_data(3, mu)[1]) == [
    (Q(1), Q(0), Q(0)), (Q(0), Q(0), Q(0))]
p4 = root_profile(4, radial_data(4, mu)[1])
assert p4 == [
    (Q(0), Q(3, 4), Q(1, 4)),
    (Q(0), Q(0), Q(0)),
    (Q(0), Q(-1, 4), Q(1, 4)),
]
# The last nonzero entry is (sqrt(6)-sqrt(2))/4 > 0.
assert 6 > 2
print("principal radial square root through M=4:", p4)

# Exact sign of a+b*sqrt(2), using rational comparisons.
def sign_a_bsqrt2(a, b):
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    if a > 0 and b < 0:
        return (a*a > 2*b*b) - (a*a < 2*b*b)
    return (2*b*b > a*a) - (2*b*b < a*a)

# M=5: enumerate every radial Fourier square-root sign choice.
failures = []
for e0, e2, e4 in product((-1, 1), repeat=3):
    # Fourier spectrum (16,4,8), hence roots (4,2,2*sqrt(2)).
    ps = [
        (Q(4*e0 + 20*e2, 16), Q(10*e4, 16)),
        (Q(4*e0 + 4*e2, 16), Q(-6*e4, 16)),
        (Q(4*e0 - 4*e2, 16), Q(2*e4, 16)),
    ]
    signs = tuple(sign_a_bsqrt2(a, b) for a, b in ps)
    if min(signs) < 0:
        failures.append(((e0, e2, e4), signs))
assert len(failures) == 8
# For the principal signs, p(1)=1/2-3*sqrt(2)/8<0 iff 16<18.
assert 16 < 18
print("M=5 radial one-square sign patterns:", failures)

# M=5 B certificate. Normalize quotient representatives to last bit 0,
# then average edge-line and triangle-plane subspaces.
def canon(v, M):
    full = (1 << M) - 1
    if (v >> (M-1)) & 1:
        v ^= full
    return v & ((1 << (M-1)) - 1)

def span(gens):
    H = {0}
    for g in gens:
        H |= {x ^ g for x in tuple(H)}
    return H

M = 5
edge_subspaces = [
    span([canon((1 << i) | (1 << j), M)])
    for i, j in combinations(range(M), 2)
]
triangle_subspaces = []
for i, j, k in combinations(range(M), 3):
    triangle_subspaces.append(span([
        canon((1 << i) | (1 << j), M),
        canon((1 << j) | (1 << k), M),
    ]))

for S in range(1 << (M-1)):
    edge_avg = Q(sum(S in H for H in edge_subspaces),
                 len(edge_subspaces))
    tri_avg = Q(sum(S in H for H in triangle_subspaces),
                len(triangle_subspaces))
    target = mu[S.bit_count()] * mu[M-S.bit_count()]
    assert 4*edge_avg + 2*tri_avg == target
print("M=5 B certificate verified on 16 quotient elements")

M = 11
profile = [mu[k] * mu[M-k] for k in range(M//2+1)]
assert profile == [1585, 0, 232, 91, 108, 90]
print("M=11 quotient radial profile:", profile)
