"""FM-SEC90 (luna_max_pluto): braided (q,s) family; boundary witnesses, 41x41 half-plane grid on 544 profiles,
q = 0 Bernstein screen, and the (1^8,2,2) Walsh-Hadamard screen."""
import itertools
from collections import defaultdict
from fractions import Fraction
from math import comb

def matchings(legs, owner):
    if not legs:
        yield ()
        return
    a = legs[0]
    for j in range(1, len(legs)):
        b = legs[j]
        if owner[a] == owner[b]:
            continue
        rest = legs[1:j] + legs[j+1:]
        for m in matchings(rest, owner):
            yield ((a, b),) + m

def crossing(e, f):
    a, b = e
    c, d = f
    return (a < c < b < d) or (c < a < d < b)

def Fqs(labels, T):
    owner = tuple(i for i, n in enumerate(labels) for _ in range(n))
    poly = defaultdict(int)
    for S in range(1 << len(labels)):
        sign = (-1) ** ((S & T).bit_count())
        left = tuple(k for k in range(len(owner))
                     if (S >> owner[k]) & 1)
        right = tuple(k for k in range(len(owner))
                      if not ((S >> owner[k]) & 1))
        if len(left) % 2:
            continue
        for M1 in matchings(left, owner):
            c1 = sum(crossing(e, f)
                     for e, f in itertools.combinations(M1, 2))
            for M2 in matchings(right, owner):
                c2 = sum(crossing(e, f)
                         for e, f in itertools.combinations(M2, 2))
                c12 = sum(crossing(e, f) for e in M1 for f in M2)
                poly[(c1 + c2, c12)] += sign
    return {k: v for k, v in poly.items() if v}

def value(poly, q, s):
    return sum(Fraction(v) * q**i * s**j
               for (i, j), v in poly.items())

def q0_slice(poly):
    out = defaultdict(int)
    for (i, j), v in poly.items():
        if i == 0:
            out[j] += v
    return {j: v for j, v in sorted(out.items()) if v}

def derivative(poly, axis):
    out = defaultdict(int)
    for (i, j), v in poly.items():
        if axis == "q" and i:
            out[(i-1, j)] += i * v
        if axis == "s" and j:
            out[(i, j-1)] += j * v
    return {k: v for k, v in out.items() if v}

def scaled_numerator(poly, qi, si):
    # q=qi/20, s=si/20; denominator 20^D, D=max(i+j).
    D = max((i+j for i, j in poly), default=0)
    numerator = sum(v * qi**i * si**j * 20**(D-i-j)
                    for (i, j), v in poly.items())
    return numerator, D

named = [
    ((1, 1, 1, 1), 0b0011),
    ((1, 1, 1, 1), 0b0101),
    ((1, 5, 2, 2), 0b0110),
]
for labels, T in named:
    print("named", labels, T, Fqs(labels, T))

long_boundary = Fqs((1,)*8 + (6,), 0b001111)
print("q=0 slice for (1^8,6):", q0_slice(long_boundary))

p01 = Fqs((1, 1, 1, 1), 0b0011)
print("derivatives for (1^4), T={0,1}:",
      "dq", derivative(p01, "q"), "ds", derivative(p01, "s"))

profiles = []
for L in range(2, 7):
    for labels in itertools.combinations_with_replacement(range(1, 4), L):
        if sum(labels) % 2 == 0 and sum(labels) <= 12:
            for T in range(1 << L):
                if T.bit_count() % 2 == 0:
                    profiles.append((labels, T, Fqs(labels, T)))

negative_total = 0
negative_region = 0
negative_cone = 0
first_near_boundary = None
min_region = None
min_q0_axis = None
min_boundary = None

for labels, T, poly in profiles:
    for qi in range(-20, 21):
        for si in range(-20, 21):
            num, D = scaled_numerator(poly, qi, si)
            if num < 0:
                negative_total += 1
            if qi == 0 and si == -1 and num < 0 and first_near_boundary is None:
                first_near_boundary = (
                    labels, T, Fraction(num, 20**D))
            if qi >= 0 and si >= -qi:
                if num < 0:
                    negative_region += 1
                v = Fraction(num, 20**D)
                if min_region is None or v < min_region[0]:
                    min_region = (v, labels, T, qi, si)
            if qi >= 0 and abs(si) <= qi and num < 0:
                negative_cone += 1
            if qi == 0 and si >= 0:
                v = Fraction(num, 20**D)
                if min_q0_axis is None or v < min_q0_axis[0]:
                    min_q0_axis = (v, labels, T, si)
            if qi >= 0 and si == -qi:
                v = Fraction(num, 20**D)
                if min_boundary is None or v < min_boundary[0]:
                    min_boundary = (v, labels, T, qi)

print("profiles", len(profiles))
print("full-square negative profile-point pairs", negative_total)
print("negative pairs in q>=0, s>=-q", negative_region)
print("negative pairs in 0<=q, |s|<=q", negative_cone)
print("first at (q,s)=(0,-1/20)", first_near_boundary)
print("minimum in q>=0, s>=-q grid", min_region)
print("minimum on q=0, 0<=s<=1 grid", min_q0_axis)
print("minimum on s=-q, 0<=q<=1 grid", min_boundary)

bernstein_failures = []
bernstein_minimum = None
for labels, T, poly in profiles:
    a = defaultdict(int)
    for (i, j), v in poly.items():
        if i == 0:
            a[j] += v
    d = max(a, default=0)
    if d:
        b = [
            sum(Fraction(a[j] * comb(k, j), comb(d, j))
                for j in range(k+1))
            for k in range(d+1)
        ]
    else:
        b = [Fraction(a[0])]
    if any(x < 0 for x in b):
        bernstein_failures.append((labels, T, b))
    for k, x in enumerate(b):
        if bernstein_minimum is None or x < bernstein_minimum[0]:
            bernstein_minimum = (x, labels, T, d, k)

print("q=0 Bernstein failures", len(bernstein_failures))
print("minimum Bernstein coefficient", bernstein_minimum)
print("Bernstein coefficients for (1^8,6), T={0,1,2,3}:")
a = q0_slice(long_boundary)
d = max(a)
print([
    sum(Fraction(a.get(j, 0) * comb(k, j), comb(d, j))
        for j in range(k+1))
    for k in range(d+1)
])


# ---- block ----
import itertools
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache

def all_Fqs(labels):
    owner = tuple(i for i, n in enumerate(labels) for _ in range(n))
    L = len(labels)

    @lru_cache(None)
    def matchings(legs):
        if not legs:
            return ((),)
        a = legs[0]
        out = []
        for j in range(1, len(legs)):
            b = legs[j]
            if owner[a] != owner[b]:
                out.extend(((a, b),) + m
                           for m in matchings(legs[1:j] + legs[j+1:]))
        return tuple(out)

    def crossing(e, f):
        a, b = e
        c, d = f
        return (a < c < b < d) or (c < a < d < b)

    # A[S] stores unsigned pairing counts by (same-colour, cross-colour)
    # crossing exponents.
    A = [defaultdict(int) for _ in range(1 << L)]
    for S in range(1 << L):
        left = tuple(k for k in range(len(owner))
                     if (S >> owner[k]) & 1)
        right = tuple(k for k in range(len(owner))
                      if not ((S >> owner[k]) & 1))
        if len(left) % 2:
            continue
        for M1 in matchings(left):
            c1 = sum(crossing(e, f)
                     for e, f in itertools.combinations(M1, 2))
            for M2 in matchings(right):
                c2 = sum(crossing(e, f)
                         for e, f in itertools.combinations(M2, 2))
                c12 = sum(crossing(e, f) for e in M1 for f in M2)
                A[S][(c1+c2, c12)] += 1

    exponents = set().union(*(set(d) for d in A))
    F = [defaultdict(int) for _ in range(1 << L)]
    for exponent in exponents:
        vals = [A[S].get(exponent, 0) for S in range(1 << L)]
        h = 1
        while h < (1 << L):
            for i in range(0, 1 << L, 2*h):
                for j in range(i, i+h):
                    x, y = vals[j], vals[j+h]
                    vals[j], vals[j+h] = x+y, x-y
            h *= 2
        for T, v in enumerate(vals):
            if v:
                F[T][exponent] = v
    return [dict(p) for p in F]

def scaled_numerator(poly, qi, si):
    D = max((i+j for i, j in poly), default=0)
    num = sum(v * qi**i * si**j * 20**(D-i-j)
              for (i, j), v in poly.items())
    return num, D

labels = (1,)*8 + (2, 2)
all_polys = all_Fqs(labels)
even_T = [T for T in range(1 << len(labels))
          if T.bit_count() % 2 == 0]
first_negative = None
minimum_region = None
minimum_q0_axis = None
boundary_negative_count = 0

for T in even_T:
    poly = all_polys[T]
    for qi in range(21):
        for si in range(-qi, 21):
            num, D = scaled_numerator(poly, qi, si)
            val = Fraction(num, 20**D)
            if num < 0 and first_negative is None:
                first_negative = (T, qi, si, val)
            if si == -qi and num < 0:
                boundary_negative_count += 1
            if minimum_region is None or val < minimum_region[0]:
                minimum_region = (val, T, qi, si)
    for si in range(21):
        num, D = scaled_numerator(poly, 0, si)
        val = Fraction(num, 20**D)
        if minimum_q0_axis is None or val < minimum_q0_axis[0]:
            minimum_q0_axis = (val, T, si)

region_points = sum(21 + qi for qi in range(21))
print("labels", labels, "even sign masks", len(even_T))
print("grid points per mask", region_points)
print("negative in q>=0, s>=-q", first_negative)
print("negative boundary values", boundary_negative_count)
print("minimum on region grid", minimum_region)
print("minimum on q=0, 0<=s<=1", minimum_q0_axis)
