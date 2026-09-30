# FM-MECH31 (astra_max_ceres) printed verifier: Hypothesis B for repeated odd labels; obstruction and K_(2,3) repair.
import argparse
from itertools import combinations_with_replacement
from fractions import Fraction as Q
from math import comb

parser = argparse.ArgumentParser(
    description="FM-MECH31 theorem and obstruction verifier")
parser.parse_args()

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def cat(h):
    return comb(2*h, h)//(h+1)

checks = 0
for L in range(2, 7):
    for degrees in combinations_with_replacement(range(1, 5), L):
        D = sum(degrees)
        if D % 2:
            continue
        for alpha, beta in (
                (Q(0), Q(1)), (Q(1, 3), Q(2, 3)),
                (Q(1), Q(2)), (Q(2), Q(2))):
            weights = []
            for A in range(1, 1 << L):
                value = Q(1, 2)
                for i, d in enumerate(degrees):
                    value *= (
                        beta**d-alpha**d if A >> i & 1
                        else alpha**d)
                assert value >= 0
                weights.append(value)

            for S in range(1 << L):
                s = sum(d for i, d in enumerate(degrees)
                        if S >> i & 1)
                left = Q(0) if s % 2 else (
                    alpha**s*beta**(D-s)
                    + alpha**(D-s)*beta**s)/2
                right = Q(0)
                if s % 2 == 0:
                    right = alpha**D
                    for A, w in enumerate(weights, 1):
                        if (S & A) == 0 or (S & A) == A:
                            right += w
                assert left == right
                checks += 1
print("exact heterogeneous monomial identities", checks)

# Column b>0 sums 1_(E intersect C_A) over every |A|=b.
# Column 0 is 1_E.
certificates = {
    2: {0: Q(1)},
    4: {0: Q(1, 2), 2: Q(1, 4)},
    6: {0: Q(25, 14), 4: Q(3, 14)},
    8: {0: Q(287, 74), 4: Q(9, 148), 6: Q(31, 148)},
    10: {0: Q(4956, 503), 6: Q(74, 503), 8: Q(14, 503)},
    12: {
        0: Q(2299, 96), 3: Q(13, 1920),
        5: Q(25, 384), 10: Q(5, 6)
    }
}
for L, coeff in certificates.items():
    assert all(v >= 0 for v in coeff.values())
    for s in range(L+1):
        value = 0 if s % 2 else sum(
            v*(1 if b == 0 else choose(s, b)+choose(L-s, b))
            for b, v in coeff.items())
        target = 0 if s % 2 else cat(s//2)*cat((L-s)//2)
        assert value == target
print("exact all-ones certificates, even lengths 2..12",
      len(certificates))

M = 6
f = [
    0 if S.bit_count() % 2 else
    cat(S.bit_count()//2)
    * (cat((M-S.bit_count())//2+1)
       - cat((M-S.bit_count())//2))
    for S in range(1 << M)
]
y = [
    1 if S.bit_count() == 2 else
    -2 if S.bit_count() == 4 else
    15 if S.bit_count() == 6 else 0
    for S in range(1 << M)
]
for A in range(1 << (M+1)):
    col = [
        int(S.bit_count() % 2 == 0
            and ((S & A) == 0 or (S & A) == A))
        for S in range(1 << M)
    ]
    assert sum(v*w for v, w in zip(y, col)) >= 0
assert sum(v*w for v, w in zip(y, f)) == -15

def span(basis):
    H = {0}
    for v in basis:
        H |= {x ^ v for x in H}
    return H

def enumerator(H):
    return [
        sum(S.bit_count() == j for S in H)
        for j in range(M+1)
    ]

E5 = span([(1 << i)^1 for i in range(1, 5)])
K23 = span([0b001111, 0b110011])
e, k = enumerator(E5), enumerator(K23)
assert e == [1, 0, 10, 0, 5, 0, 0]
assert k == [1, 0, 0, 0, 3, 0, 0]

for j in range(M+1):
    target = 0 if j % 2 else (
        cat(j//2)*(cat((M-j)//2+1)-cat((M-j)//2)))
    value = (
        2*int(j == 0)
        + Q(9, 2)*Q(e[j], comb(M, j))
        + Q(5, 2)*Q(k[j], comb(M, j)))
    assert value == target

assert sum(y[S] for S in K23) == -6
print("restricted-cone separator pairing", -15)
print("full B certificate verified: "
      "2*zero + (9/2)*avg(E5) + (5/2)*avg(K2,3)")
print("separator sum on K2,3", -6)
