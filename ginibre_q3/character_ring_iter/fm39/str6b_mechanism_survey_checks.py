#!/usr/bin/env python3
"""Exact SU(2) doubled-character checks for FM-STR6b."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

parser = argparse.ArgumentParser(description="Exact SU(2) FM3 finite-case verifier")
parser.add_argument("--help-only", action="store_true", help="print this help and exit")
args = parser.parse_args()
if args.help_only:
    raise SystemExit(0)

@lru_cache(None)
def cg(a, n):
    return tuple(range(abs(a - n), a + n + 1, 2))

def expansion(signed_labels):
    # Coefficients of prod_i (V_{n_i} x 1 + eps_i 1 x V_{n_i}).
    d = {(0, 0): 1}
    for label in signed_labels:
        n = abs(label)
        eps = -1 if label < 0 else 1
        nd = defaultdict(int)
        for (a, b), coeff in d.items():
            for c in cg(a, n):
                nd[c, b] += coeff
            for c in cg(b, n):
                nd[a, c] += eps * coeff
        d = {key: val for key, val in nd.items() if val}
    return d

def phi(signed_labels):
    return expansion(signed_labels).get((0, 0), 0)

# Each unordered signed multiset appears once: labels 1..4, lengths 1..8.
tested = zero_count = negative_count = 0
minimum = None
minimum_examples = []
minimum_ties = 0
for length in range(1, 9):
    for word_types in combinations_with_replacement(range(8), length):
        minus_count = sum(t < 4 for t in word_types)
        if minus_count & 1:
            continue
        word = [-(t + 1) if t < 4 else (t - 3) for t in word_types]
        value = phi(word)
        tested += 1
        if value < 0:
            negative_count += 1
            print("NEGATIVE", word, value)
        if value == 0:
            zero_count += 1
        if minimum is None or value < minimum:
            minimum = value
            minimum_examples = [word]
            minimum_ties = 1
        elif value == minimum:
            minimum_ties += 1
            if len(minimum_examples) < 8:
                minimum_examples.append(word)

hard_a = [-1, 2, 3, 4, -5, 6, 7, 8]
hard_b = list(range(-1, -18, -1)) + [-17]
assert len(hard_a) == 8 and sum(x < 0 for x in hard_a) == 2
assert len(hard_b) == 18 and sum(x < 0 for x in hard_b) == 18
hard_values = (phi(hard_a), phi(hard_b))

# Exact polynomial integration against independent semicircle laws.
# Even moments are Catalan numbers; odd moments vanish.
def padd(P, Q):
    R = dict(P)
    for ij, c in Q.items():
        R[ij] = R.get(ij, 0) + c
        if not R[ij]:
            del R[ij]
    return R

def pscale(P, c):
    return {ij: c*v for ij, v in P.items() if c*v}

def pmul(P, Q):
    R = defaultdict(int)
    for (i, j), a in P.items():
        for (k, ell), b in Q.items():
            R[i+k, j+ell] += a*b
    return {ij: c for ij, c in R.items() if c}

def ppow(P, n):
    R = {(0, 0): 1}
    for _ in range(n):
        R = pmul(R, P)
    return R

def psub(P, Q):
    return padd(P, pscale(Q, -1))

def catalan_moment(k):
    if k & 1:
        return 0
    j = k // 2
    return comb(2*j, j) // (j+1)

def semicircle_product_integral(P):
    return sum(c*catalan_moment(i)*catalan_moment(j)
               for (i, j), c in P.items())

x = {(1, 0): 1}
y = {(0, 1): 1}
one = {(0, 0): 1}
xmy = psub(x, y)
x2_y2 = padd(pmul(x, x), pmul(y, y))
A = ppow(padd(x, y), 2)                    # a+2b=(x+y)^2
B = ppow(psub(x2_y2, pscale(one, 2)), 2)  # (a-2)^2
S2 = psub(x2_y2, pscale(one, 2))          # a-2

def tilted_expectation(P, r):
    weight = ppow(xmy, 2*r)
    return Fraction(semicircle_product_integral(pmul(weight, P)),
                    semicircle_product_integral(weight))

cov_r1 = (tilted_expectation(pmul(A, B), 1)
          - tilted_expectation(A, 1)*tilted_expectation(B, 1))
cov_r2 = (tilted_expectation(pmul(A, S2), 2)
          - tilted_expectation(A, 2)*tilted_expectation(S2, 2))
assert cov_r1 == -1
assert cov_r2 == Fraction(-7, 25)

# A two-factor signed character product is not positive definite on SU(2)^2.
pd_test = expansion([-1, -1])
assert pd_test.get((1, 1), 0) == -2
assert pd_test.get((0, 0), 0) == 2

assert tested == 6469 and negative_count == 0 and minimum == 0
print("small-screen: labels=1..4, factor-count=1..8, unordered signed multisets")
print("cases=", tested, "zeros=", zero_count, "negative=", negative_count)
print("minimum=", minimum, "minimum_ties=", minimum_ties)
print("minimum_examples=", minimum_examples)
print("hard_a=", hard_a, "Phi=", hard_values[0])
print("hard_b=", hard_b, "Phi=", hard_values[1])
print("tilted-Haar covariance checks:", cov_r1, cov_r2)
print("PD control for (-1,-1): coeff(1,1)=", pd_test[(1, 1)],
      "coeff(0,0)=", pd_test[(0, 0)])
print("PASS: exact integer fusion and rational moment arithmetic")
