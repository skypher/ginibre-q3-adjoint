from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
import random

@lru_cache(None)
def mult(ns):
    ns = tuple(sorted(ns))
    k = len(ns)
    if k == 0:
        return 1
    total = sum(ns)
    if total & 1:
        return 0
    if k == 1:
        return 0
    if k == 2:
        return int(ns[0] == ns[1])
    half = total // 2
    ans = 0
    for mask in range(1 << k):
        shift = sum(ns[i] + 1 for i in range(k) if mask >> i & 1)
        top = half - shift + k - 2
        if top >= k - 2:
            ans += (-1 if mask.bit_count() & 1 else 1) * comb(top, k - 2)
    assert ans >= 0
    return ans

def subset(ns, mask):
    return tuple(ns[i] for i in range(8) if mask >> i & 1)

def reduced_terms(ns):
    ms = [mult(subset(ns, m)) for m in range(256)]
    f = [0] * 256
    f[0] = ms[255]
    for m in range(1, 255):
        k = m.bit_count()
        if k in (2, 3) or (k == 4 and (m & 1)):
            f[m] = ms[m] * ms[255 ^ m]
    return f, ms

def fwht(v):
    v = list(v)
    h = 1
    while h < 256:
        for lo in range(0, 256, 2 * h):
            for i in range(lo, lo + h):
                x, y = v[i], v[i + h]
                v[i], v[i + h] = x + y, x - y
        h *= 2
    return v

def direct_bivariate(ns, signs):
    # Multiply each factor into the x-coordinate or y-coordinate character row.
    tab = {(0, 0): 1}
    for n, eps in zip(ns, signs):
        out = {}
        for (x, y), c in tab.items():
            for z in range(abs(x - n), x + n + 1, 2):
                out[(z, y)] = out.get((z, y), 0) + c
            for z in range(abs(y - n), y + n + 1, 2):
                out[(x, z)] = out.get((x, z), 0) + eps * c
        tab = out
    return tab.get((0, 0), 0)

# Every label-1..8 absolute profile, with every positional signing.
profiles = signings = odd_zero = 0
minimum = None
for ns in combinations_with_replacement(range(1, 9), 8):
    profiles += 1
    vals = fwht(reduced_terms(ns)[0])
    for sg in range(256):
        signings += 1
        if sg.bit_count() & 1:
            odd_zero += 1
            continue
        phi = 2 * vals[sg]
        assert phi >= 0, (ns, sg, phi)
        minimum = phi if minimum is None else min(minimum, phi)
assert (profiles, signings, odd_zero, minimum) == (6435, 1647360, 823680, 0)
print("small exhaustive:", profiles, "profiles;", signings, "signings;",
      odd_zero, "odd-parity zeros; min Phi =", minimum)

# Compare the inclusion-exclusion evaluator with direct two-variable
# multiplication and the uncompressed subset sum.
rng = random.Random(160102)
sample = [((1, 1, 1, 2, 3, 8, 8, 8),
           (1, 1, 1, -1, 1, -1, -1, -1))]
for _ in range(120):
    ns = tuple(sorted(rng.randint(1, 40) for _ in range(8)))
    signs = tuple(rng.choice((-1, 1)) for _ in range(8))
    sample.append((ns, signs))
for ns, signs in sample:
    _, ms = reduced_terms(ns)
    direct = 0
    for m in range(256):
        neg = sum(signs[i] < 0 for i in range(8) if m >> i & 1)
        direct += (-1 if neg & 1 else 1) * ms[m] * ms[255 ^ m]
    assert direct == direct_bivariate(ns, signs)
    assert direct >= 0
print("direct bivariate cross-check:", len(sample), "words PASS")

# Independent enumeration of the finite residual cutoff.
words = assignments = 0
minimum = None
witness = None
for a0 in range(1, 3):
    for a1 in range(a0, 5):
        for a2 in range(a1, 5):
            for a3 in range(a2, 5):
                D = a0 + a1 + a2 + a3
                for a4 in range(a3, 2 * D + 1):
                    for a5 in range(a4, a4 + D + 1):
                        for a6 in range(a5, a4 + D + 1):
                            lo = max(a6, 2 * a4 - D)
                            hi = min(a4 + D, D + a4 + a5 - a6)
                            for a7 in range(lo, hi + 1):
                                ns = (a0, a1, a2, a3, a4, a5, a6, a7)
                                odd = sum(x & 1 for x in ns)
                                sector = ((a0 == 1 and odd in (4, 6))
                                          or (a0 == 2 and odd == 0))
                                total = sum(ns)
                                delta = total // 2 - a7
                                if (not sector or total & 1 or a7 < 6
                                        or delta < 8 or a6 > delta
                                        or sum(x >= 3 for x in ns[:7]) < 2):
                                    continue
                                words += 1
                                _, ms = reduced_terms(ns)
                                terms = [(m, ms[m] * ms[255 ^ m])
                                         for m in range(1, 255)
                                         if m.bit_count() in (2, 3)
                                         or (m.bit_count() == 4 and m & 1)]
                                classes = []
                                for i, x in enumerate(ns):
                                    if i == 0 or ns[i - 1] != x:
                                        classes.append(0)
                                    classes[-1] |= 1 << i
                                for choices in range(1 << len(classes)):
                                    sg = 0
                                    for j, cls in enumerate(classes):
                                        if choices >> j & 1:
                                            sg |= cls
                                    if sg.bit_count() & 1:
                                        continue
                                    assignments += 1
                                    half = ms[255]
                                    for m, coef in terms:
                                        half += (-coef if (sg & m).bit_count() & 1
                                                 else coef)
                                    phi = 2 * half
                                    assert phi >= 0, (ns, sg, phi)
                                    if minimum is None or phi < minimum:
                                        minimum = phi
                                        witness = (ns, sg)
assert (words, assignments, minimum) == (14903, 517426, 270)
print("cutoff census:", words, "words;", assignments,
      "pair-free compatible signings; min Phi =", minimum, "at", witness)

# Random exact screen beyond label 8.
for _ in range(300):
    ns = tuple(sorted(rng.randint(9, 400) for _ in range(8)))
    signs = tuple(rng.choice((-1, 1)) for _ in range(8))
    _, ms = reduced_terms(ns)
    phi = 0
    for m in range(256):
        neg = sum(signs[i] < 0 for i in range(8) if m >> i & 1)
        phi += (-1 if neg & 1 else 1) * ms[m] * ms[255 ^ m]
    assert phi >= 0, (ns, signs, phi)
print("random larger labels: 300 exact words, max label 400 PASS")
