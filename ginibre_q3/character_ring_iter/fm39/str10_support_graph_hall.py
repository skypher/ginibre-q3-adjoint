import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, prod
import random

argparse.ArgumentParser(
    description="Exact support-graph, rank, matching, and channel-table checks."
).parse_args()

def progress(message):
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(now, message, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

@lru_cache(None)
def paths(ns):
    if not ns:
        return ((0, ()),)
    return tuple((c, path+(c,))
                 for a, path in paths(ns[:-1])
                 for c in cg(a, ns[-1]))

@lru_cache(None)
def scalar_copies(word):
    ns = tuple(map(abs, word))
    out = defaultdict(lambda: [0, 0])
    for mask in range(1 << len(word)):
        ix = [i for i in range(len(word)) if mask >> i & 1]
        iy = [i for i in range(len(word)) if not (mask >> i & 1)]
        parity = sum(word[i] < 0 for i in iy) & 1
        for a, _ in paths(tuple(ns[i] for i in ix)):
            for b, _ in paths(tuple(ns[i] for i in iy)):
                out[(a, b)][parity] += 1
    return {key: tuple(value) for key, value in out.items()}

def homology_dimensions(word):
    A = scalar_copies(tuple(word[::2]))
    B = scalar_copies(tuple(word[1::2]))
    even = odd = 0
    for key in A.keys() | B.keys():
        ae, ao = A.get(key, (0, 0))
        be, bo = B.get(key, (0, 0))
        value = (ae-ao) * (be-bo)
        even += max(value, 0)
        odd += max(-value, 0)
    return even, odd

def clean(poly):
    return {e: c for e, c in poly.items() if c}

def lower(poly, ns):
    out = defaultdict(Q)
    for exponents, coefficient in poly.items():
        for i, n in enumerate(ns):
            if exponents[i] < n:
                raised = list(exponents)
                raised[i] += 1
                out[tuple(raised)] += coefficient * (n-exponents[i])
    return clean(out)

def raise_op(poly):
    out = defaultdict(Q)
    for exponents, coefficient in poly.items():
        for i, h in enumerate(exponents):
            if h:
                lowered = list(exponents)
                lowered[i] -= 1
                out[tuple(lowered)] += coefficient * h
    return clean(out)

def inner(p, q, ns):
    if len(p) > len(q):
        p, q = q, p
    return sum(
        c * q.get(e, 0) /
        prod(comb(n, h) for n, h in zip(ns, e))
        for e, c in p.items()
    )

@lru_cache(None)
def cg_basis(ns):
    if not ns:
        return ((0, (), ({(): Q(1)},)),)

    out = []
    n = ns[-1]
    for a, path, old_states in cg_basis(ns[:-1]):
        for c in cg(a, n):
            j = (a+n-c)//2
            top = defaultdict(Q)
            for h in range(j+1):
                for exponents, coefficient in old_states[h].items():
                    top[exponents+(j-h,)] += (
                        (-1)**h * comb(j, h) * coefficient
                    )
            top = clean(top)
            assert not raise_op(top)

            states = [top]
            for _ in range(c):
                next_state = lower(states[-1], ns)
                divisor = c-len(states)+1
                states.append({
                    exponents: coefficient/Q(divisor)
                    for exponents, coefficient in next_state.items()
                })
            assert not lower(states[-1], ns)
            out.append((c, path+(c,), tuple(states)))

    assert sum(c+1 for c, _, _ in out) == prod(n+1 for n in ns)
    return tuple(out)

@lru_cache(None)
def half_copies(word):
    ns = tuple(map(abs, word))
    out = defaultdict(lambda: [[], []])
    for mask in range(1 << len(word)):
        ix = [i for i in range(len(word)) if mask >> i & 1]
        iy = [i for i in range(len(word)) if not (mask >> i & 1)]
        parity = sum(word[i] < 0 for i in iy) & 1
        for a, path_x, _ in cg_basis(tuple(ns[i] for i in ix)):
            for b, path_y, _ in cg_basis(tuple(ns[i] for i in iy)):
                out[(a, b)][parity].append((mask, path_x, path_y))
    for even, odd in out.values():
        even.sort()
        odd.sort()
    return {key: (tuple(e), tuple(o)) for key, (e, o) in out.items()}

def unmatched(copies):
    out = {}
    for channel, (even, odd) in copies.items():
        matched = min(len(even), len(odd))
        by_parity = {}
        if len(even) > matched:
            by_parity[0] = even[matched:]
        if len(odd) > matched:
            by_parity[1] = odd[matched:]
        if by_parity:
            out[channel] = by_parity
    return out

def tensor(p, q):
    out = defaultdict(Q)
    for x, cx in p.items():
        for y, cy in q.items():
            out[tuple(a+b for a, b in zip(x, y))] += cx*cy
    return clean(out)

def embed(poly, indices, total):
    out = {}
    for exponents, coefficient in poly.items():
        vector = [0]*total
        for i, h in zip(indices, exponents):
            vector[i] = h
        out[tuple(vector)] = coefficient
    return out

def states_for(ns, path):
    return next(states for _, p, states in cg_basis(ns) if p == path)

def pair_invariant(ns, left, right, left_path, right_path, spin):
    total = len(ns)
    L = states_for(tuple(ns[i] for i in left), left_path)
    R = states_for(tuple(ns[i] for i in right), right_path)
    out = defaultdict(Q)
    for h in range(spin+1):
        p = tensor(embed(L[h], left, total),
                   embed(R[spin-h], right, total))
        factor = (-1)**h * comb(spin, h)
        for exponents, coefficient in p.items():
            out[exponents] += factor*coefficient
    return clean(out)

def harmonic_vector(word, A, B, channel, copy_a, copy_b):
    ns = tuple(map(abs, word))
    total = len(ns)
    mask_a, path_ax, path_ay = copy_a
    mask_b, path_bx, path_by = copy_b

    ax = tuple(A[j] for j in range(len(A)) if mask_a >> j & 1)
    ay = tuple(A[j] for j in range(len(A)) if not (mask_a >> j & 1))
    bx = tuple(B[j] for j in range(len(B)) if mask_b >> j & 1)
    by = tuple(B[j] for j in range(len(B)) if not (mask_b >> j & 1))

    x_part = pair_invariant(ns, ax, bx, path_ax, path_bx, channel[0])
    y_part = pair_invariant(ns, ay, by, path_ay, path_by, channel[1])
    vector = tensor(x_part, y_part)
    assert vector and not raise_op(vector) and not lower(vector, ns)
    y_mask = sum(1 << i for i in ay+by)
    return y_mask, vector

def harmonic_basis(word, parity):
    A = tuple(range(0, len(word), 2))
    B = tuple(range(1, len(word), 2))
    ca = unmatched(half_copies(tuple(word[i] for i in A)))
    cb = unmatched(half_copies(tuple(word[i] for i in B)))
    out = []
    for channel in sorted(ca.keys() & cb.keys()):
        for pa, copies_a in ca[channel].items():
            for pb, copies_b in cb[channel].items():
                if pa ^ pb == parity:
                    for copy_a in copies_a:
                        for copy_b in copies_b:
                            vector = harmonic_vector(
                                word, A, B, channel, copy_a, copy_b
                            )
                            out.append((channel, vector))
    return out

def cut_weight(S, T):
    common = S & T
    value = 1
    i = 0
    while common:
        if common & 1:
            value *= i+2
        common >>= 1
        i += 1
    return value

def rational_rank(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    if not a:
        return 0
    r = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        value = a[r][column]
        a[r] = [x/value for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][column]:
                value = a[i][column]
                a[i] = [x-value*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def weighted_matching(source_caps, target_caps, graph):
    left = [channel for channel, cap in source_caps.items()
            for _ in range(cap)]
    right = [channel for channel, cap in target_caps.items()
             for _ in range(cap)]
    matched = {}

    def augment(i, seen):
        for target in graph[left[i]]:
            for j, channel in enumerate(right):
                if channel != target or j in seen:
                    continue
                seen.add(j)
                if j not in matched or augment(matched[j], seen):
                    matched[j] = i
                    return True
        return False

    return sum(augment(i, set()) for i in range(len(left)))

def graph_data(word):
    even_basis = harmonic_basis(word, 0)
    odd_basis = harmonic_basis(word, 1)
    ns = tuple(map(abs, word))

    matrix = [
        [cut_weight(S, T)*inner(p, q, ns)
         for _, (T, q) in odd_basis]
        for _, (S, p) in even_basis
    ]

    source = Counter(channel for channel, _ in odd_basis)
    target = Counter(channel for channel, _ in even_basis)
    graph = defaultdict(set)
    for i, (beta, _) in enumerate(even_basis):
        for j, (alpha, _) in enumerate(odd_basis):
            if matrix[i][j]:
                graph[alpha].add(beta)

    radius = max(
        (max(abs(alpha[0]-beta[0]), abs(alpha[1]-beta[1]))
         for alpha in source for beta in graph[alpha]),
        default=0
    )
    return (even_basis, odd_basis, matrix, dict(sorted(source.items())),
            dict(sorted(target.items())),
            {a: tuple(sorted(graph[a])) for a in sorted(source)},
            rational_rank(matrix),
            weighted_matching(source, target, graph), radius)

representatives = [
    (-1,-1,-1,-2,-2,-2,3),
    (-1,-2,-2,-2,3,3,3),
    (1,1,1,1,1,-2,-2,3),
    (1,1,1,-2,-2,3,3,3),
    (1,-2,-2,3,3,3,3,3),
]
partners = [
    (1,1,1,-2,-2,-2,-3),
    (1,-2,-2,-2,-3,-3,-3),
    (-1,-1,-1,-1,-1,-2,-2,-3),
    (-1,-1,-1,-2,-2,-3,-3,-3),
    (-1,-2,-2,-3,-3,-3,-3,-3),
]
base_words = [w for pair in zip(representatives, partners) for w in pair]
extra_words = [
    (1,-2,-2,3,-4,-4),
    (1,2,2,3,-4,-4),
    (-1,2,2,-4,-4,-5),
    (1,-2,-4,-4,-4,5),
    (2,-3,-3,4,5,5),
    (-2,-3,4,5,5,5),
]

profiles = []
for counts in product(range(9), repeat=3):
    if sum(counts) > 8:
        continue
    present = [i+1 for i, count in enumerate(counts) if count]
    for signs in range(1 << len(present)):
        word = tuple(
            (-n if signs >> j & 1 else n)
            for j, n in enumerate(present)
            for _ in range(counts[n-1])
        )
        if sum(x < 0 for x in word) % 2 == 0:
            profiles.append(word)
open_profiles = [w for w in profiles if homology_dimensions(w)[1] > 0]
assert len(profiles) == 481
assert len(open_profiles) == 10
assert set(open_profiles) == set(base_words)
progress("profile box: 481 total, 10 with odd d0 homology")

for word in base_words + extra_words:
    progress("graph start " + str(word))
    he, ho = homology_dimensions(word)
    E, O, M, source, target, graph, rank, matching, radius = graph_data(word)
    assert (len(E), len(O)) == (he, ho)
    assert rank == matching == ho
    print("GRAPH", word, "H0", (he, ho),
          "rank", rank, "weighted matching", matching,
          "source caps", source, "target caps", target,
          "neighbors", graph, "Linf radius", radius, flush=True)

def channel_table(word):
    d = {(0, 0, 0): 1}
    for signed in word:
        n = abs(signed)
        q = defaultdict(int)
        for (a, b, parity), value in d.items():
            for c in cg(a, n):
                q[(c, b, parity)] += value
            for c in cg(b, n):
                q[(a, c, parity ^ (signed < 0))] += value
        d = q
    even = defaultdict(int)
    odd = defaultdict(int)
    for (a, b, parity), value in d.items():
        (odd if parity else even)[(a, b)] += value
    return {key: even[key]-odd[key] for key in even.keys() | odd.keys()}

def dimension_screen(word):
    A = channel_table(word[::2])
    B = channel_table(word[1::2])
    products = {
        key: A.get(key, 0)*B.get(key, 0)
        for key in A.keys() | B.keys()
    }
    plus = sum(max(value, 0) for value in products.values())
    minus = sum(max(-value, 0) for value in products.values())
    nplus = sum(value > 0 for value in products.values())
    nminus = sum(value < 0 for value in products.values())
    return plus, minus, plus-minus, nplus, nminus

for k in range(1, 21):
    final_sign = 1 if k % 2 == 0 else -1
    word = tuple([-i for i in range(1, k+1)] + [final_sign*(k+1)])
    progress("run " + str(k) + " " + str(dimension_screen(word)))

special = [
    ("sign-min", (-1,2,3,4,-5,6,7,8)),
    ("F1", (1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)),
]
for name, word in special:
    progress(name + " W=" + str(sum(map(abs, word))) +
             " " + str(dimension_screen(word)))

rng = random.Random(20261002)
samples = []
seen = set()
while len(samples) < 100:
    length = rng.randint(4, 12)
    magnitudes = [rng.randint(1, 14) for _ in range(length)]
    counts = Counter(magnitudes)
    signs = {n: (-1 if rng.randrange(2) else 1) for n in counts}
    if sum(counts[n] for n in counts if signs[n] < 0) % 2:
        continue
    word = tuple(
        signs[n]*n for n in sorted(counts) for _ in range(counts[n])
    )
    if sum(map(abs, word)) > 80 or word in seen:
        continue
    seen.add(word)
    samples.append(word)

positive_odd = []
vacuous = []
for word in samples:
    plus, minus, gap, _, _ = dimension_screen(word)
    assert gap >= 0
    if minus:
        positive_odd.append((gap, word, plus, minus))
    else:
        vacuous.append((word, plus))
progress("random W<=80: " + str(len(positive_odd)) +
         " nonzero odd cases, minimum " + str(min(positive_odd)) +
         ", " + str(len(vacuous)) + " vacuous")
print("PASS", flush=True)
