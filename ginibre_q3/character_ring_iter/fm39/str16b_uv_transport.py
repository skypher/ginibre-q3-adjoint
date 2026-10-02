import argparse
import random
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial

argparse.ArgumentParser(description='Exact FM-STR16 transport verifier').parse_args()
def stamp(message):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC'), message, flush=True)

def pmul(P, Q):
    R = defaultdict(Fr)
    for a, x in P.items():
        for b, y in Q.items():
            R[tuple(u + v for u, v in zip(a, b))] += x * y
    return dict(R)

@lru_cache(None)
def diagonal_matrix_coefficient(n, p):
    # D^n_pp in Sym^n(C^2), as a polynomial in (a, abar, b, bbar).
    R = defaultdict(Fr)
    for i in range(p + 1):
        j = p - i
        if 0 <= j <= n - p:
            R[(i, n - p - j, j, p - i)] += comb(p, i) * comb(n - p, j) * (-1) ** (p - i)
    return tuple(sorted(dict(R).items()))

def haar_integral(P):
    # Exact SU(2) Haar integral: int |a|^(2r)|b|^(2s) = r!s!/(r+s+1)!.
    return sum(v * Fr(factorial(a) * factorial(b), factorial(a + b + 1))
               for (a, ab, b, bb), v in P.items() if a == ab and b == bb)

@lru_cache(None)
def singlet_diagonal(items):
    # items is a sorted tuple of (highest weight n, weight k).
    if sum(k for n, k in items) != 0:
        return Fr(0)
    P = {(0, 0, 0, 0): Fr(1)}
    for n, k in items:
        P = pmul(P, dict(diagonal_matrix_coefficient(n, (n + k) // 2)))
    return haar_integral(P)

def q(items):
    return singlet_diagonal(tuple(sorted(items)))

def original_G(k, eps):
    d = {0: 1}
    for x, e in zip(k, eps):
        z = d.copy()
        for s, v in d.items():
            z[s + x] = z.get(s + x, 0) + e * v
        d = z
    return Fr(d.get(0, 0)) - Fr(d.get(2, 0) + d.get(-2, 0), 2)

def transported_G(ns, eps, k):
    # Sum_F eps_F q_F(k_F)q_Fc(k_Fc) / a_k, grouped by identical triples.
    a = q(tuple(zip(ns, k)))
    if a == 0:
        return Fr(0)
    classes = [(key, sum(v == key for v in zip(ns, k, eps)))
               for key in dict.fromkeys(zip(ns, k, eps))]
    keys, mults = [x[0] for x in classes], [x[1] for x in classes]
    total = Fr(0)
    for counts in product(*[range(m + 1) for m in mults]):
        left, right = [], []
        sign, multiplicity, left_weight = 1, 1, 0
        for (n, x, e), m, r in zip(keys, mults, counts):
            left.extend([(n, x)] * r)
            right.extend([(n, x)] * (m - r))
            sign *= e ** r
            multiplicity *= comb(m, r)
            left_weight += r * x
        if left_weight == 0:
            total += sign * multiplicity * q(left) * q(right) / a
    return total

def phi_by_character_fusion(ns, eps):
    d = {(0, 0): 1}
    for n, e in zip(ns, eps):
        z = defaultdict(int)
        for (a, b), v in d.items():
            for r in range(abs(a - n), a + n + 1, 2):
                z[r, b] += v
            for s in range(abs(b - n), b + n + 1, 2):
                z[a, s] += e * v
        d = {key: v for key, v in z.items() if v}
    return d.get((0, 0), 0)

def invariant_multiplicity(labels):
    d = {0: 1}
    for n in labels:
        z = {}
        for j, v in d.items():
            for r in range(abs(j - n), j + n + 1, 2):
                z[r] = z.get(r, 0) + v
        d = z
    return d.get(0, 0)

def phi_by_subset_multiplicities(classes):
    total = 0
    for counts in product(*[range(m + 1) for n, e, m in classes]):
        left, right, sign, multiplicity = [], [], 1, 1
        for (n, e, m), r in zip(classes, counts):
            left.extend([n] * r)
            right.extend([n] * (m - r))
            sign *= e ** r
            multiplicity *= comb(m, r)
        total += multiplicity * sign * invariant_multiplicity(left) * invariant_multiplicity(right)
    return total

def zero_sum_weights(ns):
    return (k for k in product(*[range(-n, n + 1, 2) for n in ns]) if sum(k) == 0)

def zero_weight_count(ns):
    d = {0: 1}
    for n in ns:
        z = {}
        for s, v in d.items():
            for k in range(-n, n + 1, 2):
                z[s + k] = z.get(s + k, 0) + v
        d = z
    return d.get(0, 0)

stamp('small exact transport sweep begins')
list_count = weight_count = 0
identity_mismatches = []
minimum = None
for N in range(2, 6):
    stamp(f'N={N} begins')
    level_lists = level_weights = 0
    for ns in combinations_with_replacement(range(1, 4), N):
        if sum(ns) % 2:
            continue
        for eps in product((1, -1), repeat=N):
            if eps.count(-1) % 2:
                continue
            lhs = rhs = Fr(0)
            for k in zero_sum_weights(ns):
                a = q([(n, x) for n, x in zip(ns, k)])
                g = original_G(k, eps)
                gt = transported_G(ns, eps, k)
                lhs += a * g
                rhs += a * gt
                weight_count += 1
                level_weights += 1
                if minimum is None or gt < minimum[0]:
                    minimum = (gt, (ns, eps, k, g, a))
            direct = phi_by_character_fusion(ns, eps)
            if lhs != rhs or lhs != direct:
                identity_mismatches.append((ns, eps, lhs, rhs, direct))
            list_count += 1
            level_lists += 1
    stamp(f'N={N} complete lists={level_lists} zero-weight-vectors={level_weights}')
assert list_count == 240 and weight_count == 8922 and not identity_mismatches
stamp(f'small checks pass: lists={list_count}, weight-vectors={weight_count}, mismatches=0, min Gtilde={minimum}')

# Exact counterexample to pointwise positivity of the intersection-projector transport.
ns = (1, 1, 1, 1)
eps = (1, 1, -1, -1)
k = (-1, -1, 1, 1)
a = q(list(zip(ns, k)))
g = original_G(k, eps)
gt = transported_G(ns, eps, k)
phi = phi_by_character_fusion(ns, eps)
assert (a, g, gt, phi) == (Fr(1, 3), Fr(-3), Fr(-1), 2)
stamp(f'small failure: a={a}, G={g}, Gtilde={gt}, Phi={phi}, k={k}')

# Exact check of the Gram-squared Markov transport on V_1^tensor4.
pairs = list(combinations(range(4), 2))
triples = list(combinations(range(4), 3))
A = [[int(set(S).issubset(T)) for S in pairs] for T in triples]
Pmat = []
for i, S in enumerate(pairs):
    row = []
    for j, T in enumerate(pairs):
        value = Fr(i == j)
        for r in range(4):
            for s in range(4):
                inv_entry = (Fr(1, 2) if r == s else Fr(0)) - Fr(1, 12)
                value -= A[r][i] * inv_entry * A[s][j]
        row.append(value)
    Pmat.append(row)
eps4 = (1, 1, -1, -1)
weights4 = [tuple(1 if i in S else -1 for i in range(4)) for S in pairs]
j = pairs.index((2, 3))
a4 = Pmat[j][j]
assert sum(Pmat[l][j] ** 2 for l in range(6)) == a4
markov = sum(Pmat[l][j] ** 2 * original_G(weights4[l], eps4) for l in range(6)) / a4
assert markov == -1
stamp(f'Gram-squared transport check: row sum=1, transported coefficient={markov} at k={weights4[j]}')

# The negative-G vector specified in the FM-STR16 data.
ns6, eps6, k6 = (1, 1, 1, 1, 2, 2), (1, 1, 1, 1, -1, -1), (-1, -1, -1, -1, 2, 2)
a6, g6 = q(list(zip(ns6, k6))), original_G(k6, eps6)
gt6 = transported_G(ns6, eps6, k6)
phi6 = phi_by_character_fusion(ns6, eps6)
assert (a6, g6, gt6, phi6) == (Fr(1, 5), Fr(-14), Fr(-14, 3), 28)
stamp(f'FM-STR16 negative-G witness: a={a6}, G={g6}, Gtilde={gt6}, Phi={phi6}')

# Exhaust representative weight sets in the two proved coefficient-positive classes.
def check_nonnegative_sample(ns0, eps0):
    values = [transported_G(ns0, eps0, k) for k in zero_sum_weights(ns0)]
    return len(values), min(values) if values else None
for name, ns0, eps0 in [
    ('all-plus', (1, 1, 1, 3), (1, 1, 1, 1)),
    ('eps=(-1)^n', (1, 3, 2, 4), (-1, -1, 1, 1)),
    ('all-minus odd labels', (1, 3, 5, 7), (-1, -1, -1, -1)),
]:
    count, minimum_class = check_nonnegative_sample(ns0, eps0)
    assert minimum_class is None or minimum_class >= 0
    stamp(f'class sample {name}: vectors={count}, minimum Gtilde={minimum_class}')

hard = [
    ('(+1)^2(+2)^3(-3)^14', [1] * 2 + [2] * 3 + [3] * 14, [1] * 5 + [-1] * 14,
     (1, 1, 2, -2, -2, 3, 3, 3, 3, 3, 3, 3, -3, -3, -3, -3, -3, -3, -3)),
    ('(-1,...,-8)', list(range(1, 9)), [-1] * 8, (1, 2, 3, 4, -5, -6, -7, 8)),
    ('(-1)^6(-2)^3(-3)^3(-4)(-7)', [1] * 6 + [2] * 3 + [3] * 3 + [4, 7], [-1] * 14,
     (1, 1, 1, 1, 1, 1, 2, 2, 2, -3, -3, -3, 4, -7)),
]
expected_phi = [1523121728, 980, 64930]
classes = [
    [(1, 1, 2), (2, 1, 3), (3, -1, 14)],
    [(n, -1, 1) for n in range(1, 9)],
    [(1, -1, 6), (2, -1, 3), (3, -1, 3), (4, -1, 1), (7, -1, 1)],
]
for j, (name, ns0, eps0, k0) in enumerate(hard):
    assert sum(k0) == 0 and all(k in range(-n, n + 1, 2) for n, k in zip(ns0, k0))
    direct = phi_by_character_fusion(ns0, eps0)
    subset = phi_by_subset_multiplicities(classes[j])
    g0 = original_G(k0, eps0)
    gt0 = transported_G(ns0, eps0, k0)
    assert direct == subset == expected_phi[j]
    stamp(f'hard {name}: Phi={direct} (fusion=subset); sample G={g0}; sample Gtilde={gt0}; k={k0}')

# Exact negative transported coefficients on two of the hard lists.
negative_hard = [
    ('H1', hard[0][1], hard[0][2],
     (1, 1, 0, 0, -2, -3, 3, 3, -1, -3, -1, 3, 3, -1, 1, -3, -3, -1, 3),
     Fr(-432591655232, 4465304337), Fr(6734999, 262462200), Fr(-1512)),
    ('H3', hard[2][1], hard[2][2],
     (1, 1, -1, 1, -1, 1, -2, -2, -2, 3, -3, 1, 0, 3),
     Fr(-408288, 53185), Fr(967, 37128), Fr(0)),
]
for name, ns0, eps0, k0, expect_gt, expect_a, expect_g in negative_hard:
    a0 = q(list(zip(ns0, k0)))
    gt0 = transported_G(ns0, eps0, k0)
    g0 = original_G(k0, eps0)
    assert sum(k0) == 0 and a0 == expect_a and gt0 == expect_gt and g0 == expect_g
    stamp(f'{name} negative coefficient: a={a0}, G={g0}, Gtilde={gt0}, k={k0}')

# Deterministic exact sample in H2; this is a screen, not a full census.
rng = random.Random(1610)
ns2, eps2 = hard[1][1], hard[1][2]
valid = 0
negative = []
for _ in range(500):
    sample = [rng.choice(range(-n, n + 1, 2)) for n in ns2[:-1]]
    last = -sum(sample)
    if last not in range(-ns2[-1], ns2[-1] + 1, 2):
        continue
    sample.append(last)
    valid += 1
    z = transported_G(ns2, eps2, tuple(sample))
    if z < 0:
        negative.append((tuple(sample), z))
total_h2 = zero_weight_count(ns2)
assert valid == 363 and total_h2 == 29228 and not negative
stamp(f'H2 exact screen: balanced vectors={valid}/{total_h2}; negative Gtilde=0; seed=1610')
stamp('all verifier assertions passed')
