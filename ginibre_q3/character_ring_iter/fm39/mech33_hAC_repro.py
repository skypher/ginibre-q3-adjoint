# FM-MECH33 (astra_max_ceres) printed audit: hypothesis H_AC (autocorrelation / complete positivity) certificates and insertion law.
import argparse, random
from math import comb, factorial
from fractions import Fraction as Q

parser = argparse.ArgumentParser()
parser.parse_args()

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def fusion_rows(labels):
    rows = [{0: 1}]
    for n in labels:
        old = rows[:]
        for row in old:
            out = {}
            for j, v in row.items():
                for q in range(abs(j-n), j+n+1, 2):
                    out[q] = out.get(q, 0) + v
            rows.append(out)
    return rows

def table(labels):
    rows = fusion_rows(labels[:-1])
    n, full = labels[-1], len(rows)-1
    return [rows[s].get(0, 0) * rows[full ^ s].get(n, 0)
            for s in range(full+1)]

def coef(n, k, h):
    if h == 0:
        return Q(1, k+1)
    if h == 1:
        return Q(n-k+1, k*(k+1))
    numerator = (n+1)*(n+2*h)*factorial(k-h)
    for i in range(h-2):
        numerator *= n+k+3+i
    return Q(numerator, factorial(k+1))

# Each term is a positive coefficient times the permutation average
# of an even binary code defined by the listed parity-check masks.
codes = {
    10: [
        (Q(1485,32), (1,735)),
        (Q(51,16), (7,280)),
        (Q(525,32), (7,594)),
        (Q(45,4), (7,55,881,618)),
        (Q(51,4), (31,962,743,766,130,710))
    ],
    12: [
        (Q(14157,128), (1,4088)),
        (Q(94545,2368), (7,1997,2557,3050)),
        (Q(126225,1184), (7,1592,3582,60,4029,192)),
        (Q(146421,4736), (31,)),
        (Q(10395,2368), (31,3972,4090)),
        (Q(10791,2368), (31,4075,2063))
    ]
}

radial = {}
for N, terms in codes.items():
    v = [Q(0)]*(N//2+1)
    for weight, checks in terms:
        assert weight > 0
        for x in range(1 << N):
            if (x.bit_count() % 2 == 0 and
                all((x & h).bit_count() % 2 == 0 for h in checks)):
                v[x.bit_count()//2] += weight
    radial[N] = [v[j]/C(N, 2*j) for j in range(N//2+1)]

words = entries = spheres = zero = inner = 0
for N in range(2, 13, 2):
    for n in range(1, N+1):
        f = table((1,)*N + (n,))
        if n % 2:
            assert not any(f)
            zero += 1
        else:
            k = (N-n)//2
            if n >= k-1:
                b = [coef(n, k, h) for h in range(k+1)]
                assert min(b) >= 0
                expected = [
                    sum(b[k-l]*C(2*j,j)*C(N-2*j,l-j)
                        for l in range(k+1))
                    for j in range(N//2+1)
                ]
                spheres += 1
            else:
                assert n == 2 and N in codes
                expected = radial[N]
                inner += 1
            for x, v in enumerate(f):
                assert v == (
                    0 if x.bit_count() % 2
                    else expected[x.bit_count()//2]
                )
        words += 1
        entries += len(f)

rng = random.Random(3303)
insertion = 0
for _ in range(200):
    mu = tuple(rng.randrange(1,8)
               for _ in range(rng.randrange(0,7)))
    h, n = rng.randrange(1,8), rng.randrange(1,8)
    rows = fusion_rows(mu)
    full = len(rows)-1
    A = [
        rows[s].get(0,0) *
        sum(rows[full ^ s].get(j,0)
            for j in range(abs(h-n), h+n+1, 2))
        for s in range(full+1)
    ]
    B = [
        rows[s].get(h,0)*rows[full ^ s].get(n,0)
        for s in range(full+1)
    ]
    assert table(mu+(h,n)) == A+B
    insertion += len(A)+len(B)

N = 8
f = table((1,)*N+(6,))
p = {1 << i for i in range(N)}
for x, v in enumerate(f):
    convolution = sum((y ^ x) in p for y in p)
    assert v == Q(convolution,2) + (Q(N-2,2) if x == 0 else 0)

print("boundary membership:", words, "tables;", entries, "entries")
print("sphere / code / zero:", spheres, inner, zero)
print("insertion law:", insertion, "entries in 200 cases")
print("(1^8,6): PASS")
