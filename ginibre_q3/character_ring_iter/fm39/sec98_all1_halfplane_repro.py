"""FM-SEC98 (luna_max_pluto): all-1 sector on the half-plane; grid screens; tensor Bernstein certificates above the diagonal; ordered {1,2} and {1,2,3} label screens."""
import itertools
import math
from fractions import Fraction
from math import comb
import numpy as np

N = 3
D = 21

def pairings(legs):
    if not legs:
        yield ()
        return
    a = legs[0]
    for j in range(1, len(legs)):
        b = legs[j]
        rest = legs[1:j] + legs[j+1:]
        for m in pairings(rest):
            yield ((a, b),) + m

def tensor_bernstein(P):
    # Substitute s=q+(1-q)t in the exact power polynomial.
    C = {}
    for i in range(22):
        for j in range(22 - i):
            v = int(P[i, j])
            if not v:
                continue
            for ell in range(j + 1):
                base = v * comb(j, ell)
                for h in range(ell + 1):
                    key = (i + j - ell + h, ell)
                    C[key] = C.get(key, 0) + base * (-1)**h * comb(ell, h)
    C = {key: v for key, v in C.items() if v}
    if not C:
        return [Fraction(0)]
    nq = max(r for r, ell in C)
    nt = max(ell for r, ell in C)
    Q = {}
    for u in range(nq + 1):
        for ell in range(nt + 1):
            Q[u, ell] = sum(
                Fraction(C.get((r, ell), 0) * comb(u, r), comb(nq, r))
                for r in range(u + 1)
            )
    return [
        sum(
            Q[u, ell] * Fraction(comb(k, ell), comb(nt, ell))
            for ell in range(k + 1)
        )
        for u in range(nq + 1)
        for k in range(nt + 1)
    ]

def screen_length(L):
    m = L // 2
    pair_count = math.prod(range(1, L, 2))
    ncolors = 1 << m
    nassign = pair_count * ncolors
    assert nassign * N**D < 2**63

    masks = np.arange(ncolors, dtype=np.int32)
    bits = ((masks[:, None] >> np.arange(m, dtype=np.int32)) & 1).astype(np.int32)
    ids = np.empty(nassign, dtype=np.int32)
    offset = 0

    for M in pairings(tuple(range(L))):
        chord_bits = np.array([(1 << a) | (1 << b) for a, b in M],
                              dtype=np.int32)
        edges = []
        for u, v in itertools.combinations(range(m), 2):
            a, b = M[u]
            c, d = M[v]
            if (a < c < b < d) or (c < a < d < b):
                edges.append((u, v))

        diff = np.zeros(ncolors, dtype=np.int32)
        for u, v in edges:
            diff += bits[:, u] ^ bits[:, v]
        same = len(edges) - diff
        vertex_mask = bits @ chord_bits
        ids[offset:offset+ncolors] = vertex_mask * 484 + same * 22 + diff
        offset += ncolors

    assert offset == nassign
    counts = np.bincount(ids, minlength=(1 << L) * 484)
    counts = counts.reshape(1 << L, 22, 22)
    del ids

    # Walsh-Hadamard transform from colour subsets S to sign masks T.
    F = np.zeros_like(counts)
    for i in range(22):
        for j in range(22 - i):
            vals = counts[:, i, j].copy()
            h = 1
            while h < (1 << L):
                for start in range(0, 1 << L, 2*h):
                    left = vals[start:start+h].copy()
                    right = vals[start+h:start+2*h].copy()
                    vals[start:start+h] = left + right
                    vals[start+h:start+2*h] = left - right
                h *= 2
            F[:, i, j] = vals
    del counts

    Tlist = [T for T in range(1 << L) if T.bit_count() % 2 == 0]
    Fmat = F[Tlist].reshape(len(Tlist), -1)
    bad_grid = None
    grid_points = 0

    for qi in range(N + 1):
        for si in range(-qi, N + 1):
            w = np.zeros((22, 22), dtype=np.int64)
            for i in range(22):
                for j in range(22):
                    if i + j <= D:
                        w[i, j] = qi**i * si**j * N**(D-i-j)
            values = Fmat @ w.reshape(-1)
            grid_points += 1
            if np.any(values < 0) and bad_grid is None:
                ix = int(np.where(values < 0)[0][0])
                bad_grid = (Tlist[ix], qi, si, int(values[ix]))

    bad_bernstein = None
    bernstein_count = 0
    for T in Tlist:
        coeffs = tensor_bernstein(F[T])
        bernstein_count += len(coeffs)
        negatives = [x for x in coeffs if x < 0]
        if negatives and bad_bernstein is None:
            bad_bernstein = (T, min(negatives))

    named = None
    if L == 4:
        named = {}
        for T in (0b1100, 0b1010, 0b0101):
            named[T] = [
                (i, j, int(F[T, i, j]))
                for i in range(22) for j in range(22)
                if F[T, i, j] != 0
            ]
    return len(Tlist), grid_points, len(Tlist)*grid_points, bad_grid, \
           bernstein_count, bad_bernstein, named

total_words = total_grid = total_bernstein = 0
for L in range(2, 15, 2):
    result = screen_length(L)
    words, points, evaluations, bad_grid, nbern, bad_bern, named = result
    total_words += words
    total_grid += evaluations
    total_bernstein += nbern
    print("L", L, "words", words, "grid points/word", points,
          "evaluations", evaluations, "grid negative", bad_grid,
          "Bernstein coefficients", nbern, "negative coefficient", bad_bern)
    if named is not None:
        print("L=4 named polynomials:", named)
        q, s = Fraction(0), Fraction(1, 2)
        for T in (0b1100, 0b1010, 0b0101):
            Fval = sum(
                Fraction(c) * q**i * s**j
                for i, j, c in named[T]
            )
            print("T", T, "F", Fval, "phi", Fval / 4)
print("totals", total_words, total_grid, total_bernstein)


# ---- block ----
import itertools
from collections import defaultdict
from fractions import Fraction
from math import comb

def matchings(legs, owner):
    if not legs:
        yield ()
        return
    a = legs[0]
    for j in range(1, len(legs)):
        b = legs[j]
        if owner[a] == owner[b]:
            continue
        rest = legs[1:j] + legs[j+1:]
        for m in matchings(rest, owner):
            yield ((a, b),) + m

def crossing(e, f):
    a, b = e
    c, d = f
    return (a < c < b < d) or (c < a < d < b)

def Fqs(labels, T):
    owner = tuple(i for i, n in enumerate(labels) for _ in range(n))
    poly = defaultdict(int)
    for S in range(1 << len(labels)):
        sign = (-1) ** ((S & T).bit_count())
        left = tuple(k for k in range(len(owner)) if (S >> owner[k]) & 1)
        right = tuple(k for k in range(len(owner)) if not ((S >> owner[k]) & 1))
        if len(left) % 2:
            continue
        for M1 in matchings(left, owner):
            c1 = sum(crossing(e, f) for e, f in itertools.combinations(M1, 2))
            for M2 in matchings(right, owner):
                c2 = sum(crossing(e, f) for e, f in itertools.combinations(M2, 2))
                c12 = sum(crossing(e, f) for e in M1 for f in M2)
                poly[c1 + c2, c12] += sign
    return {k: v for k, v in poly.items() if v}

def tensor_bernstein(poly):
    C = defaultdict(int)
    for (i, j), v in poly.items():
        for ell in range(j + 1):
            for h in range(ell + 1):
                C[i + j - ell + h, ell] += (
                    v * comb(j, ell) * (-1)**h * comb(ell, h)
                )
    C = {k: v for k, v in C.items() if v}
    if not C:
        return [Fraction(0)]
    nq = max(r for r, ell in C)
    nt = max(ell for r, ell in C)
    Q = {}
    for u in range(nq + 1):
        for ell in range(nt + 1):
            Q[u, ell] = sum(
                Fraction(C.get((r, ell), 0) * comb(u, r), comb(nq, r))
                for r in range(u + 1)
            )
    return [
        sum(
            Q[u, ell] * Fraction(comb(k, ell), comb(nt, ell))
            for ell in range(k + 1)
        )
        for u in range(nq + 1)
        for k in range(nt + 1)
    ]

def exact_value(poly, q, s):
    return sum(Fraction(v) * q**i * s**j
               for (i, j), v in poly.items())

# Ordered label lists in {1,2}, lengths 2 through 6.
profiles = []
for L in range(2, 7):
    for labels in itertools.product((1, 2), repeat=L):
        if sum(labels) % 2:
            continue
        for T in range(1 << L):
            if T.bit_count() % 2 == 0:
                profiles.append((labels, T, Fqs(labels, T)))

N = 10
grid_points = sum(N + qi + 1 for qi in range(N + 1))
negative = None
tested = 0
bernstein_count = 0
negative_bernstein = None

for labels, T, poly in profiles:
    coeffs = tensor_bernstein(poly)
    bernstein_count += len(coeffs)
    neg_coeffs = [x for x in coeffs if x < 0]
    if neg_coeffs and negative_bernstein is None:
        negative_bernstein = (labels, T, min(neg_coeffs))
    for qi in range(N + 1):
        for si in range(-qi, N + 1):
            value = exact_value(poly, Fraction(qi, N), Fraction(si, N))
            tested += 1
            if value < 0 and negative is None:
                negative = (labels, T, qi, si, value)

print("ordered {1,2} profiles", len(profiles),
      "grid points/profile", grid_points, "evaluations", tested,
      "negative value", negative,
      "Bernstein coefficients", bernstein_count,
      "negative Bernstein coefficient", negative_bernstein)
print("(2,2), T=empty:", Fqs((2, 2), 0))

# General labels through four blocks: tensor Bernstein screen.
profiles3 = 0
coefficients3 = 0
negative3 = None
for L in range(2, 5):
    for labels in itertools.product((1, 2, 3), repeat=L):
        if sum(labels) % 2:
            continue
        for T in range(1 << L):
            if T.bit_count() % 2:
                continue
            coeffs = tensor_bernstein(Fqs(labels, T))
            profiles3 += 1
            coefficients3 += len(coeffs)
            neg = [x for x in coeffs if x < 0]
            if neg and negative3 is None:
                negative3 = (labels, T, min(neg))
print("ordered {1,2,3}, L<=4:", profiles3, "profiles,",
      coefficients3, "coefficients, first negative", negative3)
