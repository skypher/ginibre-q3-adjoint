"""FM-SEC108 (luna_max_uranus): canonical fractional straightening charges (three normalizations); fan-in overcharge at (2^5), T={0,1}; (2^6) and (1^8,6) checks."""
import runpy, io, contextlib
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from collections import Counter
from fractions import Fraction

with contextlib.redirect_stdout(io.StringIO()):
    base = runpy.run_path(
        "ginibre_q3/character_ring_iter/fm39/sec94_charging_repro.py"
    )
cross = base["cross"]
smooth = base["smooth"]
config_weight = base["config_weight"]
matchings = base["matchings"]

def crossing_pairs(M):
    return [(e, f) for e, f in combinations(M, 2) if cross(e, f)]

@lru_cache(None)
def normal_form(M, reverse=False):
    M = tuple(sorted(tuple(sorted(e)) for e in M))
    pairs = crossing_pairs(M)
    if not pairs:
        return Counter({M: 1})
    e, f = pairs[-1] if reverse else pairs[0]
    out = Counter()
    for child in smooth(M, e, f):
        out.update(normal_form(tuple(sorted(child)), reverse))
    return out

@lru_cache(None)
def path_counts(M):
    # Count paths over all choices of crossing pair at each step.
    M = tuple(sorted(tuple(sorted(e)) for e in M))
    pairs = crossing_pairs(M)
    if not pairs:
        return Counter({M: 1})
    out = Counter()
    for e, f in pairs:
        for child in smooth(M, e, f):
            out.update(path_counts(tuple(sorted(child))))
    return out

def selected_cluster(z):
    choices = []
    for side0, side1 in z[1]:
        roots = side0 | side1
        choices.append(tuple(sorted(
            edge for root in roots for edge in z[2][root]
        )))
    return min(choices)

def terminal_terms(labels, T, M, use_paths=False):
    z = config_weight(labels, T, M)
    cluster = selected_cluster(z)
    external = tuple(edge for edge in M if edge not in set(cluster))
    terms = path_counts(cluster) if use_paths else normal_form(cluster)
    return Counter({
        tuple(sorted(external + leaf)): coefficient
        for leaf, coefficient in terms.items()
    })

def profile_data(labels, T):
    valid = {}
    for M in matchings(labels):
        z = config_weight(labels, T, M)
        if z is not None:
            valid[M] = z[0]
    negative = [M for M, w in valid.items() if w < 0]
    positive = {M: w for M, w in valid.items() if w > 0}
    return valid, negative, positive

def receipts(labels, T, variant, raw_tau=False):
    valid, negative, positive = profile_data(labels, T)
    received = Counter()
    for M in negative:
        demand = -valid[M]
        terms = terminal_terms(labels, T, M, use_paths=(variant == 3))
        positive_terms = {
            N: c for N, c in terms.items() if N in positive
        }
        if variant == 2:
            if raw_tau:
                denominator = sum(terms.values())
            else:
                denominator = sum(
                    c for N, c in terms.items() if N in valid
                )
        else:
            denominator = sum(positive_terms.values())
        if denominator:
            for N, c in positive_terms.items():
                received[N] += Fraction(demand * c, denominator)
    return valid, negative, positive, received

# First failure in length order, sorted label multisets, and all even T.
first = [None, None, None]
for length in range(2, 6):
    for labels in combinations_with_replacement(range(1, 4), length):
        if sum(labels) % 2:
            continue
        for T in range(1 << length):
            if T.bit_count() % 2:
                continue
            for variant in (1, 2, 3):
                if first[variant - 1] is not None:
                    continue
                valid, negative, positive, received = receipts(
                    labels, T, variant
                )
                over = sorted(
                    (received[N] - positive[N], N,
                     received[N], positive[N])
                    for N in received if received[N] > positive[N]
                )
                if over:
                    first[variant - 1] = (
                        labels, T, len(matchings(labels)), len(valid),
                        len(negative),
                        sum(-valid[M] for M in negative),
                        sum(positive.values()), over[0]
                    )
        if all(item is not None for item in first):
            break
    if all(item is not None for item in first):
        break
for variant, result in enumerate(first, 1):
    print("first", variant, result)

# Per-source coefficients, denominators, and charges at (2^5), T={0,1}.
labels, T = (2, 2, 2, 2, 2), 3
valid, negative, positive = profile_data(labels, T)
target = ((0, 9), (1, 2), (3, 4), (5, 6), (7, 8))
totals = [Fraction(0), Fraction(0), Fraction(0)]
for M in sorted(negative):
    z = config_weight(labels, T, M)
    cluster = selected_cluster(z)
    assert normal_form(cluster) == normal_form(cluster, True)

    terms = terminal_terms(labels, T, M)
    paths = terminal_terms(labels, T, M, use_paths=True)
    positive_den = sum(c for N, c in terms.items() if N in positive)
    admissible_den = sum(c for N, c in terms.items() if N in valid)
    path_den = sum(c for N, c in paths.items() if N in positive)
    demand = -valid[M]

    charges = (
        Fraction(demand * terms[target], positive_den),
        Fraction(demand * terms[target], admissible_den),
        Fraction(demand * paths[target], path_den),
    )
    totals = [x + y for x, y in zip(totals, charges)]
    print("source", M, "weight", valid[M],
          "c_target", terms[target], "positive_den", positive_den,
          "admissible_den", admissible_den,
          "path_target", paths[target], "path_den", path_den,
          "charges", charges,
          "zero_terms", [(N, c) for N, c in terms.items()
                         if N in valid and valid[N] == 0])
print("target weight and totals", positive[target], totals)

# Requested (2^6) check and raw-formal variant-(ii) denominator.
for variant in (1, 2, 3):
    valid, negative, positive, received = receipts(
        (2, 2, 2, 2, 2, 2), (1 << 4) | (1 << 5), variant
    )
    over = sorted(
        (received[N] - positive[N], N, received[N], positive[N])
        for N in received if received[N] > positive[N]
    )
    print("2^6", variant, over[0] if over else None)
valid, negative, positive, received = receipts(
    (2, 2, 2, 2, 2, 2), (1 << 4) | (1 << 5),
    variant=2, raw_tau=True
)
over = sorted(
    (received[N] - positive[N], N, received[N], positive[N])
    for N in received if received[N] > positive[N]
)
print("2^6 raw tau", over[0] if over else None)

# All 36 level-1 pairs for (1^8,6).
labels = (1, 1, 1, 1, 1, 1, 1, 1, 6)
over_counts = [0, 0, 0]
sample = {}
for pair in combinations(range(9), 2):
    T = sum(1 << i for i in pair)
    for variant in (1, 2, 3):
        valid, negative, positive, received = receipts(labels, T, variant)
        over = sorted(
            (received[N] - positive[N], N, received[N], positive[N])
            for N in received if received[N] > positive[N]
        )
        over_counts[variant - 1] += bool(over)
        if pair == (0, 1) and over:
            sample[variant] = over[0]
print("(1^8,6) pair overcounts", over_counts)
for variant in (1, 2, 3):
    print("(1^8,6) T={0,1}", variant, sample.get(variant))
