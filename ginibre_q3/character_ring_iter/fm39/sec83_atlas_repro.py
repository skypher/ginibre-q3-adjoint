"""FM-SEC83 (luna_max_jupiter): noncrossing-matching model check, tightness atlas (L <= 10, labels <= 5), nested-chord switch obstruction."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb

@lru_cache(None)
def fusion_invariant(labels):
    state = {0: 1}
    for n in labels:
        nxt = {}
        for j, mult in state.items():
            for q in range(abs(j - n), j + n + 1, 2):
                nxt[q] = nxt.get(q, 0) + mult
        state = nxt
    return state.get(0, 0)

@lru_cache(None)
def matching_count(labels):
    # Points of each cluster are consecutive in cyclic order.
    owner = []
    for i, n in enumerate(labels):
        owner.extend([i] * n)
    size = len(owner)
    if size == 0:
        return 1
    if size % 2:
        return 0

    dp = [[0] * size for _ in range(size)]
    for length in range(2, size + 1, 2):
        for left in range(size - length + 1):
            right = left + length - 1
            total = 0
            for partner in range(left + 1, right + 1, 2):
                if owner[left] == owner[partner]:
                    continue
                inner = 1 if partner == left + 1 else dp[left + 1][partner - 1]
                outer = 1 if partner == right else dp[partner + 1][right]
                total += inner * outer
            dp[left][right] = total
    return dp[0][size - 1]

# Check the noncrossing-matching model against fusion on every label
# multiset of length at most 10, including the empty list.
matching_checks = 0
for length in range(11):
    for labels in combinations_with_replacement(range(1, 6), length):
        assert matching_count(labels) == fusion_invariant(labels), labels
        matching_checks += 1
print("matching/fusion comparisons:", matching_checks)

def fwht(values):
    out = values[:]
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2 * step):
            for j in range(start, start + step):
                x, y = out[j], out[j + step]
                out[j], out[j + step] = x + y, x - y
        step *= 2
    return out

def equal_label_groups(labels):
    groups = []
    i = 0
    while i < len(labels):
        j = i + 1
        while j < len(labels) and labels[j] == labels[i]:
            j += 1
        groups.append((i, j - i, labels[i]))
        i = j
    return groups

def canonical_even_T(labels):
    groups = equal_label_groups(labels)
    for counts in product(*(range(mult + 1) for _, mult, _ in groups)):
        if sum(counts) % 2:
            continue
        mask = 0
        t_labels = []
        orbit = 1
        for (start, mult, label), count in zip(groups, counts):
            for j in range(count):
                mask |= 1 << (start + j)
            t_labels.extend([label] * count)
            orbit *= comb(mult, count)
        yield mask, tuple(t_labels), orbit

stats = Counter()
zero_cases = defaultdict(list)
half_unit_cases = defaultdict(list)
raw_even_cases = 0
raw_parity_zeros = 0
raw_nonparity_zeros = 0
raw_half_units = 0
negative_cases = []

for length in range(1, 11):
    for labels in combinations_with_replacement(range(1, 6), length):
        full = (1 << length) - 1
        m = [
            fusion_invariant(tuple(labels[i] for i in range(length)
                                   if (S >> i) & 1))
            for S in range(1 << length)
        ]
        F = fwht([m[S] * m[full ^ S] for S in range(1 << length)])
        odd_total = sum(labels) % 2
        stats["sorted_lists"] += 1

        for T, t_labels, orbit in canonical_even_T(labels):
            value = F[T]
            stats["T_orbits"] += 1
            raw_even_cases += orbit
            assert value % 2 == 0
            if odd_total:
                assert value == 0
                stats["parity_zero_orbits"] += 1
                raw_parity_zeros += orbit
            elif value == 0:
                zero_cases[labels].append(t_labels)
                raw_nonparity_zeros += orbit
            elif value == 1:
                stats["nonempty_F1_orbits"] += 1
            elif value == 2:
                half_unit_cases[labels].append(t_labels)
                raw_half_units += orbit
            elif value < 0:
                negative_cases.append((labels, t_labels, value))

assert not negative_cases
assert stats["nonempty_F1_orbits"] == 0
print("atlas summary:", dict(stats))
print("raw even-T cases:", raw_even_cases)
print("raw parity-zero cases:", raw_parity_zeros)
print("non-parity zero cases:", raw_nonparity_zeros)
print("F=2 cases (F/2=1):", raw_half_units)
print("nonempty F=1 cases:", stats["nonempty_F1_orbits"])
print("empty-list F:", fusion_invariant(()))
print("non-parity zeros:", dict(sorted(zero_cases.items())))
print("F=2 orbit representatives:", dict(sorted(half_unit_cases.items())))

# Exact witness for the first nested-chord obstruction.
labels = (1, 1, 1, 1)
full = 0b1111
T = 0b0011          # clusters 0, 1
S = 0b0110          # clusters 1, 2

def inv_mask(mask):
    return fusion_invariant(tuple(labels[i] for i in range(4)
                                  if (mask >> i) & 1))

def term(mask):
    return (-1) ** ((mask & T).bit_count()) * inv_mask(mask) * inv_mask(full ^ mask)

F_witness = sum(term(mask) for mask in range(1 << 4))
assert inv_mask(S) == inv_mask(full ^ S) == 1
assert term(S) == -1
assert F_witness == 2

# One-cluster T flips leave the matching support.
for i in (0, 1):
    flipped = S ^ (1 << i)
    assert inv_mask(flipped) * inv_mask(full ^ flipped) == 0

def crosses(edge1, edge2):
    a, b = sorted(edge1)
    c, d = sorted(edge2)
    return (a < c < b < d) or (c < a < d < b)

D1 = (1, 2)
D2 = (0, 3)
assert not crosses(D1, D2)
print("switch witness: F =", F_witness,
      "negative term =", term(S),
      "D1 =", D1, "D2 =", D2,
      "crossing =", crosses(D1, D2))

# No supported negative summand exists for lengths 1, 2, or 3
# anywhere in the requested label range.
for length in range(1, 4):
    for small_labels in combinations_with_replacement(range(1, 6), length):
        full = (1 << length) - 1
        m = [
            fusion_invariant(tuple(small_labels[i] for i in range(length)
                                   if (S0 >> i) & 1))
            for S0 in range(1 << length)
        ]
        for T0 in range(1 << length):
            if T0.bit_count() % 2:
                continue
            for S0 in range(1 << length):
                if (S0 & T0).bit_count() % 2:
                    assert m[S0] * m[full ^ S0] == 0
print("shorter-list supported negative terms: none")
