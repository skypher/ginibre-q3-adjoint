"""FM-SEC77 (luna_max_pluto): q-Hermite (q-Gaussian) deformation of the EVEN form; exact q-positivity screens."""
import argparse
from functools import lru_cache
from itertools import combinations_with_replacement
from fractions import Fraction
import random

parser = argparse.ArgumentParser(
    description="Exact q-Hermite even-sign profile screens")
parser.parse_args()

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(a, b):
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0)
        + (b[i] if i < len(b) else 0)
        for i in range(n)
    ])

def sub(a, b):
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0)
        - (b[i] if i < len(b) else 0)
        for i in range(n)
    ])

def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)

def shift(p, k):
    return (0,) * k + tuple(p)

def eval_q(p, q):
    return sum(Fraction(c) * q**i for i, c in enumerate(p))

@lru_cache(None)
def qbinom(n, k):
    if k < 0 or k > n:
        return (0,)
    if k == 0 or k == n:
        return (1,)
    return add(qbinom(n - 1, k),
               shift(qbinom(n - 1, k - 1), n - k))

@lru_cache(None)
def qfactorial(n):
    p = (1,)
    for j in range(1, n + 1):
        p = mul(p, (1,) * j)
    return p

@lru_cache(None)
def linearization(a, b, k):
    return mul(mul(qbinom(a, k), qbinom(b, k)),
               qfactorial(k))

@lru_cache(None)
def moment_q(labels):
    # Coefficient of H_0 in the product of the listed q-Hermites.
    current = {0: (1,)}
    for n in labels:
        nxt = {}
        for d, p in current.items():
            for k in range(min(d, n) + 1):
                degree = d + n - 2 * k
                term = mul(p, linearization(d, n, k))
                nxt[degree] = add(nxt.get(degree, (0,)), term)
        current = nxt
    return current.get(0, (0,))

def all_profiles(labels):
    L = len(labels)
    values = []
    for S in range(1 << L):
        left = tuple(sorted(labels[i] for i in range(L) if S >> i & 1))
        right = tuple(sorted(labels[i] for i in range(L)
                             if not (S >> i & 1)))
        values.append(mul(moment_q(left), moment_q(right)))

    # Polynomial-valued Walsh–Hadamard transform.
    h = 1
    while h < len(values):
        for i in range(0, len(values), 2 * h):
            for j in range(i, i + h):
                a, b = values[j], values[j + h]
                values[j], values[j + h] = add(a, b), sub(a, b)
        h *= 2
    return values

def screen(labels):
    F = all_profiles(labels)
    even_T = [T for T in range(1 << len(labels))
              if T.bit_count() % 2 == 0]
    bad = [T for T in even_T if any(c < 0 for c in F[T])]
    return len(even_T), bad, F

def batch(name, words):
    profiles = bad_count = 0
    for index, labels in enumerate(words):
        count, bad, _ = screen(labels)
        profiles += count
        bad_count += len(bad)
        print(name, index, labels, "patterns", count,
              "negative coefficient patterns", len(bad), flush=True)
    print(name, "TOTAL", len(words), "lists", profiles,
          "patterns", bad_count, "negative profiles", flush=True)
    return profiles, bad_count

# Exhaustive small labels.
small_lists = small_profiles = small_bad = 0
for L in range(2, 7):
    level_lists = 0
    for labels in combinations_with_replacement(range(1, 5), L):
        if sum(labels) % 2 or sum(labels) > 16:
            continue
        count, bad, _ = screen(labels)
        small_lists += 1
        level_lists += 1
        small_profiles += count
        small_bad += len(bad)
    print("small length", L, "lists", level_lists,
          "profiles so far", small_profiles,
          "negative profiles", small_bad, flush=True)
assert (small_lists, small_profiles, small_bad) == (90, 1564, 0)

# Targeted lists and controls.
targeted = [
    (1,) * 8 + (6,),
    (1, 5, 2, 2),
    (1,) * 16,
    (1,) * 12 + (2,) * 4,
    (1,) * 8 + (2,) * 8,
    (1,) * 8 + (3, 5),
    (1, 1, 2, 2, 3, 3, 4, 4),
    (1, 1, 2, 2, 3, 3, 4, 4, 5, 7),
    (1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
]
target_profiles, target_bad = batch("targeted", targeted)
assert (target_profiles, target_bad) == (101768, 0)

# Structured controls plus 32 seeded mixed lists.
family = []
for a in (2, 4, 6, 8):
    family.append((1,) * 8 + (a,))
for n in (1, 3, 5, 7):
    for L in (4, 6, 8):
        family.append((n,) * L)
for n in (2, 4, 6, 8):
    for L in (3, 5, 7):
        family.append((n,) * L)

rng = random.Random(77477)
for _ in range(32):
    L = rng.randint(6, 10)
    labels = [rng.randint(1, 8) for _ in range(L)]
    if sum(labels) % 2:
        j = rng.randrange(L)
        labels[j] = labels[j] + 1 if labels[j] < 8 else labels[j] - 1
    family.append(tuple(labels))

family_profiles, family_bad = batch("structured-plus-random", family)
assert (len(family), family_profiles, family_bad) == (60, 7696, 0)

# Long one-high-label families.
long_high = [
    (1,) * 14 + (2,),
    (1,) * 12 + (4,),
    (1,) * 10 + (6,),
    (1,) * 8 + (8,),
]
high_profiles, high_bad = batch("long-one-high", long_high)
assert (high_profiles, high_bad) == (21760, 0)

# Exact representative polynomials and endpoint minima.
B = (1,) * 8 + (6,)
FB = all_profiles(B)
tight_labels = (1, 5, 2, 2)
FT = all_profiles(tight_labels)
assert FB[15] == (
    6, 42, 154, 400, 828, 1450, 2220, 3028, 3724, 4162, 4246,
    3962, 3382, 2636, 1868, 1196, 686, 348, 152, 54, 14, 2
)
assert FT[3] == (2, 8, 18, 30, 40, 44, 40, 30, 18, 8, 2)
even_B = [T for T in range(1 << len(B)) if T.bit_count() % 2 == 0]
print("B list minimum F(0), F(1):",
      min(FB[T][0] for T in even_B),
      min(sum(FB[T]) for T in even_B), flush=True)
print("B list T=15 polynomial:", FB[15], flush=True)
print("tight word F polynomial:", FT[3], flush=True)
print("tight word half-polynomial:", tuple(c // 2 for c in FT[3]),
      flush=True)


# ---- block ----
from fractions import Fraction
from itertools import product
from math import comb, factorial

def gaussian_moment(blocks):
    current = {0: 1}
    for n in blocks:
        nxt = {}
        for d, value in current.items():
            for k in range(min(d, n) + 1):
                degree = d + n - 2 * k
                c = comb(d, k) * comb(n, k) * factorial(k)
                nxt[degree] = nxt.get(degree, 0) + value * c
        current = nxt
    return current.get(0, 0)

def rotation_value(labels, T):
    D = sum(labels)
    if D % 2:
        return Fraction(0)
    signs = [-1 if T >> i & 1 else 1 for i in range(len(labels))]
    scale = Fraction(2) ** (len(labels) - D // 2)
    total = Fraction(0)
    for ks in product(*(range(n + 1) for n in labels)):
        if any(signs[i] != (-1) ** (labels[i] - ks[i])
               for i in range(len(labels))):
            continue
        coefficient = scale
        for n, k in zip(labels, ks):
            coefficient *= comb(n, k)
        total += coefficient * gaussian_moment(ks) * gaussian_moment(
            tuple(n - k for n, k in zip(labels, ks)))
    return total

checks = 0
for labels in ((1,) * 8 + (6,), (1, 5, 2, 2)):
    F = all_profiles(labels)  # defined in the preceding block
    for T in range(1 << len(labels)):
        if T.bit_count() % 2:
            continue
        assert rotation_value(labels, T) == sum(F[T])
        checks += 1
print("exact q=1 endpoint equalities:", checks, flush=True)
