import argparse
from fractions import Fraction as Q
from math import comb, factorial
from time import monotonic
import sympy as sp

argparse.ArgumentParser(
    description='Exact FM-SEC126 verifier; screens total degree <= 60.'
).parse_args()

u, v = sp.symbols('u v')
Zs = (u + v)/2 - 2
Ps = (u - v)/4
Qplus = sp.expand(Zs + Ps)
Qminus = sp.expand(Zs - Ps)
assert sp.expand(Qplus - ((3*u + v)/4 - 2)) == 0
assert sp.expand(Qminus - ((u + 3*v)/4 - 2)) == 0


def cat(n):
    return comb(2*n, n)//(n+1)


def mu(m, k):
    return Q(2*factorial(2*m)*factorial(2*m+1)
             *factorial(2*k)*factorial(2*k+1),
             factorial(m)**2*factorial(k)**2
             *factorial(m+k+1)*factorial(m+k+2))


def mu_by_catalan_sum(m, k):
    total = Q(0)
    for i in range(2*m+1):
        for j in range(2*k+1):
            px = 2*m - i + 2*k - j
            py = i + j
            if px % 2 == 0 and py % 2 == 0:
                total += (Q(comb(2*m, i)*comb(2*k, j)*(-1)**j)
                          * cat(px//2)*cat(py//2))
    return total


for m in range(7):
    for k in range(7):
        assert mu(m, k) == mu_by_catalan_sum(m, k)

D = 30  # A+E+2b+2al+2ga <= 60, A=2m, E=2k.
vals = {}
for base_sum in range(D+1):
    for m in range(base_sum+1):
        k = base_sum - m
        vals[(m, k, 0, 0, 0)] = mu(m, k)

started = monotonic()
for t in range(1, D+1):
    made = 0
    for b in range(t+1):
        for al in range(t-b+1):
            ga = t - b - al
            for m in range(D-t+1):
                for k in range(D-t-m+1):
                    if b:
                        value = (vals[(m+1,k,b-1,al,ga)]
                                 + vals[(m,k+1,b-1,al,ga)])/2
                        value -= 2*vals[(m,k,b-1,al,ga)]
                    elif al:
                        value = (3*vals[(m+1,k,b,al-1,ga)]
                                 + vals[(m,k+1,b,al-1,ga)])/4
                        value -= 2*vals[(m,k,b,al-1,ga)]
                    else:
                        value = (vals[(m+1,k,b,al,ga-1)]
                                 + 3*vals[(m,k+1,b,al,ga-1)])/4
                        value -= 2*vals[(m,k,b,al,ga-1)]
                    vals[(m,k,b,al,ga)] = value
                    made += 1
    print(f'layer exponent_total={t} states={made}', flush=True)

negatives = []
zeros = []
inside = outside = negative_outside = 0
min_ratio = None
min_ties = []
smallest = {}
max_degree_states = 0
slice2_count = slice2_zero = slice2_negative = 0
slice2_min = None
slice2_ties = []

for (m,k,b,al,ga), value in vals.items():
    degree = m+k+b+al+ga
    if degree == D:
        max_degree_states += 1
    A, E = 2*m, 2*k
    key = (A,E,b,al,ga)

    if value < 0:
        negatives.append((key, value))
        if not (al <= E and ga <= A):
            negative_outside += 1
    elif value == 0:
        zeros.append(key)
    else:
        ratio = value/mu(m,k)
        if min_ratio is None or ratio < min_ratio:
            min_ratio = ratio
            min_ties = [(key, value)]
        elif ratio == min_ratio:
            min_ties.append((key, value))

        if ratio in smallest:
            smallest[ratio][0] += 1
            if len(smallest[ratio][1]) < 6:
                smallest[ratio][1].append((key, value))
        elif len(smallest) < 8 or ratio < max(smallest):
            smallest[ratio] = [1, [(key, value)]]
            if len(smallest) > 8:
                del smallest[max(smallest)]

    if al <= E and ga <= A:
        inside += 1
    else:
        outside += 1

    if al + ga == 2:
        slice2_count += 1
        if value == 0:
            slice2_zero += 1
        elif value < 0:
            slice2_negative += 1
        else:
            ratio2 = value/mu(m,k)
            if slice2_min is None or ratio2 < slice2_min:
                slice2_min = ratio2
                slice2_ties = [(key, value)]
            elif ratio2 == slice2_min:
                slice2_ties.append((key, value))

assert len(vals) == 324632
assert max_degree_states == 46376
assert not negatives
assert len(zeros) == 8
assert min_ratio == Q(1,5)
assert set(key for key,_ in min_ties) == {
    (4,0,0,0,1), (0,4,0,1,0)
}
assert (slice2_count,slice2_zero,slice2_negative,slice2_min) == (
    13485, 1, 0, Q(4,7)
)
assert set(key for key,_ in slice2_ties) == {
    (6,0,0,0,2), (0,6,0,2,0)
}

q1 = vals[(0,0,0,1,0)]
q2 = vals[(0,0,0,2,0)]
q3 = vals[(0,0,0,3,0)]
r1 = vals[(0,0,0,0,1)]
r2 = vals[(0,0,0,0,2)]
r3 = vals[(0,0,0,0,3)]
assert (q1,q2,q3) == (Q(0),Q(3),Q(8))
assert (r1,r2,r3) == (Q(0),Q(3),Q(8))
assert q3 - 4*q2 + 4*q1 == -4
assert r3 - 4*r2 + 4*r1 == -4

print('degree_cap=60 states=',len(vals),
      'exact_degree_60_states=',max_degree_states)
print('negative_values=',len(negatives),
      'negative_outside_constraints=',negative_outside)
print('inside_constraints=',inside,'outside_constraints=',outside)
print('zero_cases=',sorted(zeros))
print('minimum_normalized=',min_ratio,'ties=',min_ties)
print('smallest_distinct_normalized_values=')
for ratio in sorted(smallest):
    print(' ',ratio,'count=',smallest[ratio][0],
          'examples=',smallest[ratio][1])
print('alpha_plus_gamma_2_count=',slice2_count,
      'zeros=',slice2_zero,'negatives=',slice2_negative,
      'minimum_normalized=',slice2_min,'ties=',slice2_ties)
print('closure moments Q=',(q1,q2,q3),
      'R=',(r1,r2,r3),'multiplication margins=',(-4,-4))
print('elapsed_seconds=',round(monotonic()-started,2))
