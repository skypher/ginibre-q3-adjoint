"""FM-SEC121 (luna_max_venus): B census for labels <= 4 at L = 7, 8, 9 (subgroup-orbit certificates), Fourier screens, L = 10 samples."""
from functools import lru_cache
from itertools import combinations_with_replacement

@lru_cache(None)
def inv_counts(counts, labs):
    row = {0: 1}
    for n, c in zip(labs, counts):
        for _ in range(c):
            out = {}
            for p, v in row.items():
                for q in range(abs(p-n), p+n+1, 2):
                    out[q] = out.get(q, 0) + v
            row = out
    return row.get(0, 0)

def wht(a):
    a = list(a)
    h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            for j in range(i, i+h):
                x, y = a[j], a[j+h]
                a[j], a[j+h] = x+y, x-y
        h *= 2
    return a

def fourier_screen(L, maxlab):
    lists = negative = 0
    for labels in combinations_with_replacement(range(1, maxlab+1), L):
        if sum(labels) % 2:
            continue
        labs = tuple(sorted(set(labels)))
        cache = {}
        def moment(mask):
            counts = tuple(sum(1 for i, n in enumerate(labels)
                               if n == lab and (mask >> i) & 1)
                           for lab in labs)
            if counts not in cache:
                cache[counts] = inv_counts(counts, labs)
            return cache[counts]
        full = (1 << L) - 1
        f = [moment(S) * moment(full ^ S) for S in range(1 << L)]
        F = wht(f)
        assert all(z == 0 for T, z in enumerate(F) if T.bit_count() % 2)
        negative += sum(z < 0 for T, z in enumerate(F)
                        if T.bit_count() % 2 == 0)
        lists += 1
    print("L", L, "label bound", maxlab,
          "even-total lists", lists, "negative Fourier entries", negative)

for L, bound in ((9, 3), (10, 3), (7, 4), (8, 4), (9, 4)):
    fourier_screen(L, bound)


# ---- block ----
from itertools import combinations, combinations_with_replacement
from functools import lru_cache
from fractions import Fraction as Q
from random import Random
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

@lru_cache(None)
def inv(labels):
    row = {0: 1}
    for n in labels:
        out = {}
        for p, v in row.items():
            for q in range(abs(p-n), p+n+1, 2):
                out[q] = out.get(q, 0) + v
        row = out
    return row.get(0, 0)

def subspaces(d):
    for k in range(d+1):
        for piv in combinations(range(d), k):
            free = [(i,j) for i,p in enumerate(piv)
                    for j in range(p+1,d) if j not in piv]
            for mask in range(1 << len(free)):
                basis = [1 << p for p in piv]
                for h,(i,j) in enumerate(free):
                    if (mask >> h) & 1:
                        basis[i] |= 1 << j
                members = [0]
                for b in basis:
                    members += [x ^ b for x in members]
                yield members

def profile_data(labels):
    L = len(labels)
    G = 1 << (L-1)
    labs = tuple(sorted(set(labels)))
    totals = tuple(labels.count(n) for n in labs)

    def canonical(x):
        c = tuple(sum(1 for i,n in enumerate(labels)
                      if n == lab and ((x >> i) & 1))
                  for lab in labs)
        return min(c, tuple(totals[j]-c[j] for j in range(len(labs))))

    p_of_x = [canonical(x) for x in range(G)]
    profiles = sorted(set(p_of_x))
    ix = {p:i for i,p in enumerate(profiles)}
    qprof = [ix[p] for p in p_of_x]
    sizes = [0]*len(profiles)
    for q in qprof:
        sizes[q] += 1

    target = []
    for p in profiles:
        left = tuple(labs[j] for j,n in enumerate(p) for _ in range(n))
        right = tuple(labs[j] for j,n in enumerate(p)
                      for _ in range(totals[j]-n))
        target.append(inv(left)*inv(right)*sizes[len(target)])
    return qprof, profiles, sizes, target

def exact_certificate(cols, target, seed):
    if not any(target):
        return (), Q(0)
    A = np.asarray(cols, dtype=float).T
    b = np.asarray(target, dtype=float)
    rng = Random(seed)
    for _ in range(12):
        cost = np.asarray([1.0 + rng.randrange(1_000_000)/1_000_000
                           for _ in range(len(cols))])
        lp = linprog(cost, A_eq=A, b_eq=b, bounds=(0,None), method="highs")
        if not lp.success:
            return None, ("LP infeasible", lp.message)
        active = [j for j,x in enumerate(lp.x) if x > 1e-8]
        mat = Matrix([[cols[j][i] for j in active]
                      for i in range(len(target))])
        try:
            sol, params = mat.gauss_jordan_solve(Matrix(target))
        except ValueError:
            continue
        if params.rows:
            continue
        coeff = [Q(int(v.p), int(v.q)) for v in sol]
        if any(v < 0 for v in coeff):
            continue
        got = [sum(coeff[h]*cols[j][i] for h,j in enumerate(active))
               for i in range(len(target))]
        if got == [Q(x) for x in target]:
            return tuple((j,v) for j,v in zip(active,coeff) if v), min(coeff)
    return None, "no exact basic support after 12 cost perturbations"

def certify(labels):
    L = len(labels)
    qprof, profiles, sizes, target = profile_data(labels)
    seen = set()
    nsub = 0
    for members in subspaces(L-1):
        nsub += 1
        hist = [0]*len(profiles)
        for x in members:
            hist[qprof[x]] += 1
        seen.add(tuple(hist))
    cols = list(seen)
    support, minc = exact_certificate(
        cols, target, hash(labels) & 0x7fffffff)
    if support is None:
        return {"ok": False, "labels": labels, "subspaces": nsub,
                "profiles": len(profiles), "columns": len(cols),
                "reason": minc}
    return {"ok": True, "labels": labels, "subspaces": nsub,
            "profiles": len(profiles), "columns": len(cols),
            "terms": len(support), "min_coefficient": minc,
            "zero": not any(target)}

for L in (7,8,9):
    words = [w for w in combinations_with_replacement(range(1,5), L)
             if sum(w) % 2 == 0]
    max_terms = max_cols = zero = 0
    min_c = None
    for k,w in enumerate(words, 1):
        rec = certify(w)
        assert rec["ok"], rec
        max_terms = max(max_terms, rec["terms"])
        max_cols = max(max_cols, rec["columns"])
        zero += int(rec["terms"] == 0)
        if rec["terms"] and (min_c is None or rec["min_coefficient"] < min_c):
            min_c = rec["min_coefficient"]
        if k % 10 == 0 or k == len(words):
            print("PROGRESS L", L, k, "/", len(words),
                  "max_terms", max_terms, "max_columns", max_cols, flush=True)
    print("SUMMARY L", L, "lists", len(words), "zero_tables", zero,
          "max_terms", max_terms, "max_columns", max_cols,
          "minimum_coefficient", min_c, flush=True)

samples = [
    (1,)*8 + (3,3),
    (1,)*7 + (2,2,3),
    (1,)*6 + (2,2,3,3),
]
for w in samples:
    print("SAMPLE L10", certify(w), flush=True)
