import argparse
argparse.ArgumentParser(
    description="Independent exact check of the STR8b channel complex."
).parse_args()

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, prod
import random

def channels(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def fusion_paths(ns):
    if not ns:
        return ((0, ()),)
    return tuple(
        (c, path+(c,))
        for a, path in fusion_paths(ns[:-1])
        for c in channels(a, ns[-1])
    )

def fusion_counts(ns):
    return Counter(a for a, _ in fusion_paths(tuple(ns)))

def inv0(ns):
    d = {0: 1}
    for n in ns:
        e = defaultdict(int)
        for a, count in d.items():
            for c in channels(a, n):
                e[c] += count
        d = dict(e)
    return d.get(0, 0)

def phi(word):
    ns = tuple(abs(x) for x in word)
    total = 0
    for mask in range(1 << len(word)):
        X = [ns[i] for i in range(len(word)) if not (mask >> i) & 1]
        Y = [ns[i] for i in range(len(word)) if (mask >> i) & 1]
        sign = (-1) ** sum(
            word[i] < 0 for i in range(len(word)) if (mask >> i) & 1
        )
        total += sign * inv0(X) * inv0(Y)
    return total

def chain_dims(word):
    ns = tuple(abs(x) for x in word)
    dims = [0, 0]
    for mask in range(1 << len(word)):
        X = [ns[i] for i in range(len(word)) if not (mask >> i) & 1]
        Y = [ns[i] for i in range(len(word)) if (mask >> i) & 1]
        parity = sum(
            word[i] < 0 for i in range(len(word)) if (mask >> i) & 1
        ) & 1
        dims[parity] += inv0(X) * inv0(Y)
    return tuple(dims)

@lru_cache(None)
def scalar_copies(word):
    ns = tuple(abs(x) for x in word)
    groups = defaultdict(lambda: [0, 0])
    for mask in range(1 << len(word)):
        X = [i for i in range(len(word)) if not (mask >> i) & 1]
        Y = [i for i in range(len(word)) if (mask >> i) & 1]
        mx = fusion_counts(tuple(ns[i] for i in X))
        my = fusion_counts(tuple(ns[i] for i in Y))
        parity = sum(word[i] < 0 for i in Y) & 1
        for a, ca in mx.items():
            for b, cb in my.items():
                groups[(a, b)][parity] += ca * cb
    return {t: tuple(v) for t, v in groups.items()}

def unmatched_counts(groups):
    result = {}
    for typ, (even, odd) in groups.items():
        matched = min(even, odd)
        if even > matched:
            result[(typ, 0)] = even - matched
        if odd > matched:
            result[(typ, 1)] = odd - matched
    return result

def d0_homology(A, B):
    fA = {t: e-o for t, (e, o) in A.items()}
    fB = {t: e-o for t, (e, o) in B.items()}
    he = ho = 0
    for typ in set(fA) | set(fB):
        value = fA.get(typ, 0) * fB.get(typ, 0)
        he += max(value, 0)
        ho += max(-value, 0)
    sa, sb = unmatched_counts(A), unmatched_counts(B)
    actual = tuple(
        sum(nA * sb.get((typ, p ^ q), 0)
            for (typ, q), nA in sa.items())
        for p in (0, 1)
    )
    assert actual == (he, ho)
    return he, ho, sa, sb

def verify_d0_blocks(groups):
    D = ((0, 1), (0, 0))
    G = (1, -1)
    for even, odd in groups.values():
        for _ in range(min(even, odd)):
            assert all(
                sum(D[i][k]*D[k][j] for k in range(2)) == 0
                for i in range(2) for j in range(2)
            )
            assert all(
                D[i][j]*G[j] + G[i]*D[i][j] == 0
                for i in range(2) for j in range(2)
            )
            assert G[0]*G[0] == G[1]*G[1] == 1

def clean(poly):
    return {e: c for e, c in poly.items() if c}

def lowering(poly, ns):
    out = {}
    for e, c in poly.items():
        for i, n in enumerate(ns):
            if e[i] < n:
                q = list(e)
                q[i] += 1
                q = tuple(q)
                out[q] = out.get(q, Fraction(0)) + c*(n-e[i])
    return clean(out)

def raising(poly, ns):
    out = {}
    for e, c in poly.items():
        for i, n in enumerate(ns):
            if e[i]:
                q = list(e)
                q[i] -= 1
                q = tuple(q)
                out[q] = out.get(q, Fraction(0)) + c*e[i]
    return clean(out)

@lru_cache(None)
def cg_paths(ns):
    if not ns:
        return ((0, (), ({(): Fraction(1)},)),)
    n = ns[-1]
    result = []
    for a, path, old_states in cg_paths(ns[:-1]):
        for c in channels(a, n):
            j = (a+n-c)//2
            top = {}
            for h in range(j+1):
                for e, value in old_states[h].items():
                    q = e+(j-h,)
                    top[q] = top.get(q, Fraction(0)) + (-1)**h*comb(j, h)*value
            top = clean(top)
            assert not raising(top, ns)
            states = [top]
            for h in range(c):
                nxt = lowering(states[-1], ns)
                states.append({
                    e: value/Fraction(c-h) for e, value in nxt.items()
                })
            assert not lowering(states[-1], ns)
            result.append((c, path+(c,), tuple(states)))
    assert sum(c+1 for c, _, _ in result) == prod(n+1 for n in ns)
    return tuple(result)

@lru_cache(None)
def cg_lookup(ns):
    return {path: (spin, states) for spin, path, states in cg_paths(ns)}

@lru_cache(None)
def polynomial_copies(word):
    ns = tuple(abs(x) for x in word)
    groups = defaultdict(lambda: [[], []])
    for mask in range(1 << len(word)):
        X = [i for i in range(len(word)) if not (mask >> i) & 1]
        Y = [i for i in range(len(word)) if (mask >> i) & 1]
        nx, ny = tuple(ns[i] for i in X), tuple(ns[i] for i in Y)
        parity = sum(word[i] < 0 for i in Y) & 1
        for a, px, _ in cg_paths(nx):
            for b, py, _ in cg_paths(ny):
                groups[(a, b)][parity].append((mask, px, py))
    return {
        t: (tuple(sorted(e)), tuple(sorted(o)))
        for t, (e, o) in groups.items()
    }

def unmatched_paths(groups):
    result = {}
    for typ, (even, odd) in groups.items():
        matched = min(len(even), len(odd))
        if len(even) > matched:
            result[(typ, 0)] = even[matched:]
        if len(odd) > matched:
            result[(typ, 1)] = odd[matched:]
    return result

def embed(poly, positions, total):
    out = {}
    for e, c in poly.items():
        q = [0]*total
        for j, i in enumerate(positions):
            q[i] = e[j]
        q = tuple(q)
        out[q] = out.get(q, Fraction(0)) + c
    return clean(out)

def pair_invariant(ns, left, right, left_path, right_path, spin):
    nl = tuple(ns[i] for i in left)
    nr = tuple(ns[i] for i in right)
    sl, L = cg_lookup(nl)[left_path]
    sr, R = cg_lookup(nr)[right_path]
    assert sl == sr == spin
    total = len(ns)
    out = {}
    for h in range(spin+1):
        x = embed(L[h], left, total)
        y = embed(R[spin-h], right, total)
        factor = (-1)**h * comb(spin, h)
        for ex, cx in x.items():
            for ey, cy in y.items():
                q = tuple(i+j for i, j in zip(ex, ey))
                out[q] = out.get(q, Fraction(0)) + factor*cx*cy
    return clean(out)

def harmonic_vector(word, posA, posB, copyA, copyB, typ):
    ns = tuple(abs(x) for x in word)
    ma, pax, pay = copyA
    mb, pbx, pby = copyB
    Ax = [posA[i] for i in range(len(posA)) if not (ma >> i) & 1]
    Ay = [posA[i] for i in range(len(posA)) if (ma >> i) & 1]
    Bx = [posB[i] for i in range(len(posB)) if not (mb >> i) & 1]
    By = [posB[i] for i in range(len(posB)) if (mb >> i) & 1]
    x = pair_invariant(ns, Ax, Bx, pax, pbx, typ[0])
    y = pair_invariant(ns, Ay, By, pay, pby, typ[1])
    out = {}
    for ex, cx in x.items():
        for ey, cy in y.items():
            q = tuple(i+j for i, j in zip(ex, ey))
            out[q] = out.get(q, Fraction(0)) + cx*cy
    global_mask = sum(1 << i for i in Ay+By)
    return clean(out), global_mask

def harmonic_basis(word, posA, posB, sa, sb, parity):
    basis = []
    for (typ, pA), copiesA in sa.items():
        copiesB = sb.get((typ, parity ^ pA), ())
        for copyA in copiesA:
            for copyB in copiesB:
                vector, mask = harmonic_vector(
                    word, posA, posB, copyA, copyB, typ
                )
                basis.append((mask, vector))
    return basis

def inner(p, q, ns):
    ans = Fraction(0)
    for e, c in p.items():
        if e in q:
            norm = Fraction(1)
            for i, n in enumerate(ns):
                norm /= comb(n, e[i])
            ans += c*q[e]*norm
    return ans

def rational_rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return 0
    rows, cols, rank = len(a), len(a[0]), 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        value = a[rank][col]
        a[rank] = [x/value for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                value = a[i][col]
                a[i] = [x-value*y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank

def correction_matrix(word, even_basis, odd_basis):
    ns = tuple(abs(x) for x in word)
    matrix = []
    for maskE, vE in even_basis:
        row = []
        for maskO, vO in odd_basis:
            shared = maskE & maskO
            weight = 1
            i = 0
            while shared:
                if shared & 1:
                    weight *= i+2
                shared >>= 1
                i += 1
            row.append(weight * inner(vE, vO, ns))
        matrix.append(row)
    return matrix

def corrected_rank(word, expected_he, expected_ho):
    posA = list(range(0, len(word), 2))
    posB = list(range(1, len(word), 2))
    A = tuple(word[i] for i in posA)
    B = tuple(word[i] for i in posB)
    PA, PB = polynomial_copies(A), polynomial_copies(B)
    SA, SB = scalar_copies(A), scalar_copies(B)
    assert {t: tuple(map(len, v)) for t, v in PA.items()} == SA
    assert {t: tuple(map(len, v)) for t, v in PB.items()} == SB
    sa, sb = unmatched_paths(PA), unmatched_paths(PB)
    even = harmonic_basis(word, posA, posB, sa, sb, 0)
    odd = harmonic_basis(word, posA, posB, sa, sb, 1)
    assert (len(even), len(odd)) == (expected_he, expected_ho)
    matrix = correction_matrix(word, even, odd)
    return rational_rank(matrix), (len(matrix), len(matrix[0]) if matrix else 0)

def profiles_1_to_3(maxlen):
    for counts in product(range(maxlen+1), repeat=3):
        if sum(counts) > maxlen:
            continue
        present = [i+1 for i, count in enumerate(counts) if count]
        for signs in range(1 << len(present)):
            word = tuple(
                (-n if (signs >> j) & 1 else n)
                for j, n in enumerate(present)
                for _ in range(counts[n-1])
            )
            if sum(x < 0 for x in word) % 2 == 0:
                yield word

def display(word):
    counts = Counter(word)
    return ' '.join(
        (f'{n:+d}^{m}' if m > 1 else f'{n:+d}')
        for n, m in sorted(counts.items(), key=lambda z:(abs(z[0]), z[0]))
    ) or 'empty'

profiles = list(profiles_1_to_3(8))
assert len(profiles) == 481
open_cases = []
odd_chain_profiles = 0
for word in profiles:
    A = scalar_copies(tuple(word[i] for i in range(0, len(word), 2)))
    B = scalar_copies(tuple(word[i] for i in range(1, len(word), 2)))
    verify_d0_blocks(A)
    verify_d0_blocks(B)
    he, ho, _, _ = d0_homology(A, B)
    assert he-ho == phi(word)
    dims = chain_dims(word)
    odd_chain_profiles += dims[1] > 0
    assert dims[0]-dims[1] == phi(word)
    if ho:
        open_cases.append((word, he, ho))
assert len(open_cases) == 10
print("SCREEN profiles=481 odd-chain=88 nonzero-Hodd(d0)=10 Hodd(d)=0")
for word, he, ho in open_cases:
    rank, shape = corrected_rank(word, he, ho)
    assert rank == ho
    print("CORRECTION", display(word), "H0=", (he, ho),
          "rank=", rank, "H=", (he-rank, ho-rank), "matrix=", shape)

# Eighth distinct unconditioned random length-9 draw, seed 20261002.
candidate_rng = random.Random(20261002)
def random_word(length, rng):
    while True:
        counts = [0]*4
        for _ in range(length):
            counts[rng.randrange(4)] += 1
        present = [i+1 for i, count in enumerate(counts) if count]
        signs = {n: (-1 if rng.randrange(2) else 1) for n in present}
        word = tuple(signs[n]*n for n in present for _ in range(counts[n-1]))
        if sum(x < 0 for x in word) % 2 == 0:
            return word

seen = set()
draws = []
while len(draws) < 8:
    word = random_word(9, candidate_rng)
    if word not in seen:
        seen.add(word)
        draws.append(word)
extra = draws[-1]
assert extra == (-1, -1, 2, 2, 2, 3, 3, 4, 4)
posA = list(range(0, len(extra), 2))
posB = list(range(1, len(extra), 2))
A = scalar_copies(tuple(extra[i] for i in posA))
B = scalar_copies(tuple(extra[i] for i in posB))
he, ho, _, _ = d0_homology(A, B)
assert (he, ho) == (638, 2) and he-ho == phi(extra)
rank, shape = corrected_rank(extra, he, ho)
assert rank == ho
print("RANDOM_POSITIVE", display(extra), "H0=", (he, ho),
      "rank=", rank, "Hodd(d)=", ho-rank, "matrix=", shape)

rng = random.Random(20261002)
for length in (9, 10):
    accepted = []
    draws = 0
    while len(accepted) < 12:
        draws += 1
        word = random_word(length, rng)
        A = scalar_copies(tuple(word[i] for i in range(0, length, 2)))
        B = scalar_copies(tuple(word[i] for i in range(1, length, 2)))
        he, ho, _, _ = d0_homology(A, B)
        assert he-ho == phi(word)
        if ho == 0:
            accepted.append(word)
    print("RANDOM_ZERO_HODD0", "length=", length, "n=", len(accepted),
          "draws=", draws, "all Hodd(d)=0")
    for word in accepted:
        print("SAMPLE", display(word))
print("PASS")
