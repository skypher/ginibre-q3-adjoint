# FM-MECH30 (astra_max_ceres) printed reproducer: Hypothesis B certificates (subspace mixtures), obstructions.
import argparse
from functools import lru_cache
from itertools import combinations_with_replacement
from fractions import Fraction
from math import comb
import random
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix

parser = argparse.ArgumentParser(
    description="FM-MECH30 exact certificate verification")
parser.add_argument("--max-count", type=int, default=8)
args = parser.parse_args()

@lru_cache(None)
def fusion(labels):
    b = {0: 1}
    for n in labels:
        out = {}
        for p, v in b.items():
            for q in range(abs(p-n), p+n+1, 2):
                out[q] = out.get(q, 0) + v
        b = out
    return b

def inv(labels):
    return fusion(labels).get(0, 0)

def table(labels):
    n = len(labels)
    return [
        inv(tuple(labels[i] for i in range(n) if S >> i & 1))
        * inv(tuple(labels[i] for i in range(n)
                    if not (S >> i & 1)))
        for S in range(1 << (n-1))
    ]

def subspaces(d):
    ans = [(0,)]
    layer = {(0,)}
    for _ in range(d):
        nxt = set()
        for H in layer:
            seen = set(H)
            for v in range(1, 1 << d):
                if v in seen:
                    continue
                coset = {v ^ x for x in H}
                seen.update(coset)
                nxt.add(tuple(sorted(set(H) | coset)))
        ans.extend(sorted(nxt))
        layer = nxt
    return ans

def solve_exact(A, b):
    M = [list(map(int, row)) + [int(v)]
         for row, v in zip(A, b)]
    rows, cols, prev = len(M), len(A[0]), 1
    for k in range(cols):
        p = next(i for i in range(k, rows) if M[i][k])
        M[p], M[k] = M[k], M[p]
        pivot = M[k][k]
        for i in range(k+1, rows):
            for j in range(k+1, cols+1):
                z = pivot*M[i][j] - M[i][k]*M[k][j]
                assert z % prev == 0
                M[i][j] = z // prev
            M[i][k] = 0
        prev = pivot
    assert all(M[i][cols] == 0 for i in range(cols, rows))
    x = [Fraction(0)] * cols
    for i in range(cols-1, -1, -1):
        x[i] = (
            Fraction(M[i][cols])
            - sum((M[i][j]*x[j] for j in range(i+1, cols)),
                  Fraction(0))
        ) / M[i][i]
    return x

def certify(A, b):
    result = linprog(
        np.zeros(A.shape[1]), A_eq=A,
        b_eq=np.array(b, dtype=float)/b[0],
        bounds=(0, None), method="highs")
    assert result.success, (
        "Search failed; this is not an exact separation.")
    active = [i for i, v in enumerate(result.x) if v > 1e-8]
    M = A[:, active].toarray().astype(int).tolist()
    coeff = solve_exact(M, b)
    assert all(v >= 0 for v in coeff)
    assert all(
        sum((Fraction(M[i][j])*coeff[j]
             for j in range(len(active))), Fraction(0)) == b[i]
        for i in range(len(b)))
    return len(active)

counts = [0, 0]
for n in range(2, args.max_count+1):
    codes = subspaces(n-1)
    labels_max = 4 if n <= 7 else 3
    for labels in combinations_with_replacement(
            range(1, labels_max+1), n):
        f = table(labels)
        if not f[0]:
            assert not any(f)
            counts[1] += 1
            continue
        allowed = [H for H in codes if all(f[x] for x in H)]
        rows, cols = [], []
        for j, H in enumerate(allowed):
            rows.extend(H)
            cols.extend([j]*len(H))
        A = csc_matrix(
            (np.ones(len(rows)), (rows, cols)),
            shape=(len(f), len(allowed)))
        certify(A, f)
        counts[0] += 1
print("general label lists: nonzero certificates, zero tables =",
      counts, flush=True)

rng = random.Random(303012)

def span(basis):
    H = {0}
    for v in basis:
        H |= {x ^ v for x in H}
    return H

repeated = 0
for k in range(1, 13):
    R = [inv((2,)*j) for j in range(k+2)]
    target = [
        comb(k, j)*R[j]*(R[k-j]+R[k-j+1])
        for j in range(k+1)
    ]
    K = [
        [
            sum((-1)**b*comb(w, b)*comb(k-w, j-b)
                for b in range(max(0, j-k+w), min(j, w)+1))
            for j in range(k+1)
        ]
        for w in range(k+1)
    ]
    enums = {(1,)+(0,)*k}

    def put_generator(basis):
        e = [0]*(k+1)
        for x in span(basis):
            e[x.bit_count()] += 1
        if not e[1]:
            enums.add(tuple(e))

    for size in range(2, k+1):
        put_generator([(1 << i)^1 for i in range(1, size)])

    for batch in range(10):
        for _ in range(3000):
            if rng.randrange(2):
                rank = rng.randrange(1, min(k, 7)+1)
                put_generator([
                    rng.randrange(1, 1 << k) for _ in range(rank)])
            else:
                h = rng.randrange(1, min(k, 7)+1)
                columns = [
                    rng.randrange(1, 1 << h) for _ in range(k)]
                weights = [
                    sum((u&t).bit_count() % 2 for t in columns)
                    for u in range(1 << h)
                ]
                nums = [
                    sum(K[w][j] for w in weights)
                    for j in range(k+1)
                ]
                assert all(v % (1 << h) == 0 for v in nums)
                e = tuple(v // (1 << h) for v in nums)
                assert e[1] == 0
                enums.add(e)

        E = sorted(enums)
        A = csc_matrix(np.array(E, dtype=int).T)
        result = linprog(
            np.zeros(len(E)), A_eq=A,
            b_eq=np.array(target, dtype=float)/target[0],
            bounds=(0, None), method="highs")
        if result.success:
            certify(A, target)
            repeated += 1
            break
    else:
        raise AssertionError(("dictionary insufficient", k))

print("lists (1,1,2^k), k=1..12: exact certificates =",
      repeated, flush=True)

B = {1, 2, 4, 8, 15}
g = [3 if x == 0 else int(x in B) for x in range(16)]
F = [
    sum((-1)**((u&x).bit_count())*g[x] for x in range(16))
    for u in range(16)
]
assert min(F) == 0
assert all((x^y) not in B for x in B for y in B if x != y)
assert all(
    sum((u&x).bit_count() % 2 == 0 for x in B) <= 3
    for u in range(1, 16))
print("gluing obstruction: Fourier min =", min(F),
      "; origin needs 5, has 3")

labels = (1, 1, 2)
v = []
for S in range(4):
    a = fusion(tuple(labels[i] for i in range(3) if S >> i & 1))
    b = fusion(tuple(labels[i] for i in range(3)
                     if not (S >> i & 1)))
    v.append(sum(m*b.get(j, 0) for j, m in a.items() if j <= 1))
assert v == [1, 1, 1, 0]
assert sum((-1)**S.bit_count()*v[S] for S in range(4)) == -1
print("spin cutoff J=1: values", v, "; Fourier value = -1")

def mul(f, g):
    out = {}
    for (i, j), v in f.items():
        for (k, l), w in g.items():
            out[i+k, j+l] = out.get((i+k, j+l), 0) + v*w
    return {ij: v for ij, v in out.items() if v}

def power(f, n):
    out = {(0, 0): 1}
    for _ in range(n):
        out = mul(out, f)
    return out

def h(k):
    return {
        (k-2*b-j, j): (-1)**b*comb(k+1-b, b)
        for b in range(k//2+1)
        for j in range(k+1-2*b)
    }

def hat(k):
    out = {}
    for b in range(k//2+1):
        j, v = k-2*b, (-1)**b*comb(k-b, b)
        out[j, 0] = out.get((j, 0), 0) + v
        out[0, j] = out.get((0, j), 0) + v
    return out

def cat(k):
    return 0 if k % 2 else comb(k, k//2)//(k//2+1)

def phi(w, r):
    return Fraction(
        sum(v*(-1)**b*comb(2*r, b)*cat(i+2*r-b)*cat(j+b)
            for (i, j), v in w.items()
            for b in range(2*r+1)), 2)

print("phi_1(S2^k), k=1..12:",
      [int(phi(power(hat(2), k), 1)) for k in range(1, 13)])
print("phi_1(h2^4), phi_2(h2^6):",
      [int(phi(power(h(2), 2*r+2), r)) for r in (1, 2)])
w = mul(mul(h(3), power(hat(3), 2)), power(h(1), 3))
print("phi_r(h3 S3^2 h1^3), r=0..4:",
      [int(phi(w, r)) for r in range(5)])
