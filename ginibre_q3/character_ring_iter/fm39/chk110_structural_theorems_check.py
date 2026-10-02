#!/usr/bin/env python3
"""Fresh exact checks of FM-STR1b Theorem 3 and FM-STR7 Proposition 2."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product

ap = argparse.ArgumentParser(description="Exact structural theorem checks for FM3")
ap.add_argument("--help-only", action="store_true", help="print help and exit")
args = ap.parse_args()
if args.help_only:
    raise SystemExit(0)

@lru_cache(None)
def cg(a, n):
    return tuple(range(abs(a-n), a+n+1, 2))

@lru_cache(None)
def fusion(ns):
    d = {0: 1}
    for n in ns:
        e = defaultdict(int)
        for a, multiplicity in d.items():
            for c in cg(a, n):
                e[c] += multiplicity
        d = dict(e)
    return d

def mult(ns, j):
    return fusion(tuple(sorted(ns))).get(j, 0)

@lru_cache(None)
def path_basis(ns, target):
    """Multiplicity-free fusion-tree paths from V_0 to V_target."""
    ns = tuple(ns)
    if not ns:
        return ((),) if target == 0 else ()
    suffix = [fusion(ns[k:]) for k in range(len(ns)+1)]
    out = []
    def rec(k, a, path):
        if k == len(ns):
            if a == target:
                out.append(tuple(path))
            return
        for c in cg(a, ns[k]):
            if any(target in cg(c, b) for b in suffix[k+1]):
                rec(k+1, c, path+[c])
    rec(0, 0, [])
    return tuple(out)

def signed_dims(minus_labels, plus_labels):
    labels = tuple([-n for n in minus_labels] + list(plus_labels))
    positive = negative = 0
    for mask in range(1 << len(labels)):
        x = tuple(sorted(abs(labels[i]) for i in range(len(labels))
                         if mask >> i & 1))
        y = tuple(sorted(abs(labels[i]) for i in range(len(labels))
                         if not (mask >> i & 1)))
        term = mult(x, 0) * mult(y, 0)
        minus_parity = sum(labels[i] < 0 for i in range(len(labels))
                           if not (mask >> i & 1)) & 1
        if minus_parity:
            negative += term
        else:
            positive += term
    return positive, negative

def pb_term(labels, mask):
    S = tuple(labels[i] for i in range(len(labels)) if mask >> i & 1)
    T = tuple(labels[i] for i in range(len(labels)) if not (mask >> i & 1))
    complement_factor = (mult(T, 0) + 3*mult(T, 2)
                         + 2*mult(T, 4) + mult(T, 6))
    return mult(S, 2) * complement_factor

def check_pb(labels):
    L = len(labels)
    for mask in range(1, (1 << L)-1):
        assert pb_term(labels, mask) == 0, (labels, mask, pb_term(labels, mask))
    return (1 << L)-2

def rank_q(matrix):
    A = [[Fraction(x) for x in row] for row in matrix]
    if not A:
        return 0
    rank = 0
    for col in range(len(A[0])):
        pivot = next((r for r in range(rank, len(A)) if A[r][col]), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        z = A[rank][col]
        A[rank] = [v/z for v in A[rank]]
        for r in range(len(A)):
            if r != rank and A[r][col]:
                z = A[r][col]
                A[r] = [u-z*v for u,v in zip(A[r], A[rank])]
        rank += 1
        if rank == len(A):
            break
    return rank

def det_q(matrix):
    A = [[Fraction(x) for x in row] for row in matrix]
    d = Fraction(1)
    for c in range(len(A)):
        pivot = next((r for r in range(c, len(A)) if A[r][c]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != c:
            A[c], A[pivot] = A[pivot], A[c]
            d = -d
        z = A[c][c]
        d *= z
        A[c] = [v/z for v in A[c]]
        for r in range(c+1, len(A)):
            z = A[r][c]
            A[r] = [u-z*v for u,v in zip(A[r], A[c])]
    return d

# Operator basis: source (orientation, distinguished minus, beta fusion path).
pair_list = list(combinations(range(4), 2))
base_rows = []
for i,j in pair_list:
    row = [0]*4
    row[i], row[j] = 1, -2
    base_rows.append(row)
assert rank_q(base_rows) == 4
assert det_q([base_rows[k] for k in (0,1,3,2)]) == 4

def operator_matrix(b):
    source = [(o,i,k) for o in range(2) for i in range(4) for k in range(b)]
    target = [(o,p,k) for o in range(2) for p in range(6) for k in range(b)]
    M = [[0]*len(source) for _ in target]
    for r,(o,p,k) in enumerate(target):
        i,j = pair_list[p]
        for col,key in enumerate(source):
            oo,s,kb = key
            if oo == o and kb == k:
                if s == i:
                    M[r][col] += 1
                if s == j:
                    M[r][col] -= 2
    return source,target,M

# SO(3) coordinate tensors:
# h_i=(-1)^i <beta,x_i> epsilon(x_0,...,omit x_i,...,x_3).
def levi(v):
    if len(set(v)) < 3:
        return 0
    return (-1)**sum(v[i] > v[j] for i in range(3) for j in range(i+1,3))

def make_h(i):
    H = {}
    for z in product(range(3), repeat=5):
        beta, xs = z[0], z[1:]
        omitted = tuple(xs[j] for j in range(4) if j != i)
        c = (-1)**i * (beta == xs[i]) * levi(omitted)
        if c:
            H[z] = c
    return H

def contract_delta(H, i, j):
    out = defaultdict(int)
    for z,c in H.items():
        if z[i+1] == z[j+1]:
            out[(z[0],)+tuple(z[k+1] for k in range(4) if k not in (i,j))] += c
    return {z:c for z,c in out.items() if c}

H = [make_h(i) for i in range(4)]
for z in product(range(3), repeat=5):
    assert sum(P.get(z,0) for P in H) == 0
for i,j in pair_list:
    C = [contract_delta(P,i,j) for P in H]
    assert all(not C[k] for k in range(4) if k not in (i,j))
    assert C[i] == {z:-v for z,v in C[j].items()}
    Mpsi = [[C[i].get((beta,u,v),0) for beta in range(3)]
            for u,v in product(range(3), repeat=2)]
    assert rank_q(Mpsi) == 3
    gram = [[sum(Mpsi[r][a]*Mpsi[r][b] for r in range(9))
             for b in range(3)] for a in range(3)]
    assert gram == [[2,0,0],[0,2,0],[0,0,2]]

# PB and exact dimensions for q,t <= 12, q,t != 2.
qt_cases = qt_nonzero_b = qt_sum_phi = 0
for q in range(1,13):
    for t in range(q,13):
        if q == 2 or t == 2:
            continue
        labels = (q,t)
        assert check_pb(labels) == 2
        beta_paths = path_basis(labels, 2)
        b = mult(labels, 2)
        assert len(beta_paths) == b
        even,odd = signed_dims((2,2,2,2), labels)
        assert odd == 8*b and even >= odd
        source,target,M = operator_matrix(b)
        assert len(source) == 8*b and len(target) == 12*b
        assert rank_q(M) == 8*b
        qt_cases += 1
        qt_nonzero_b += bool(b)
        qt_sum_phi += even-odd

# PB family A_L=(4,8,...,2^L,2^(L+1)-2), checked through L=8.
family = []
proper_subsets = 0
for L in range(2,9):
    labels = tuple(4*(1 << j) for j in range(L-1))
    labels += (sum(labels)+2,)
    assert labels[-1] - sum(labels[:-1]) == 2
    for mask in range(1, (1 << L)-1):
        S = tuple(labels[i] for i in range(L) if mask >> i & 1)
        # For even labels, a gap greater than 2 excludes V_2.
        assert 2*max(S)-sum(S) > 2
        assert mult(S,2) == 0
        assert pb_term(labels,mask) == 0
        proper_subsets += 1
    assert check_pb(labels) == (1 << L)-2
    beta_paths = path_basis(labels,2)
    b = mult(labels,2)
    assert b == 1 and len(beta_paths) == 1
    even,odd = signed_dims((2,2,2,2), labels)
    assert odd == 8*b and even >= odd
    source,target,M = operator_matrix(b)
    assert rank_q(M) == 8*b
    family.append((L,labels,b,beta_paths[0],even,odd,even-odd))

# Folded path pairing for all words with each signed label occurring evenly.
def signed_label(code):
    n = code//2+1
    return -n if code & 1 else n

def half_paths(Hword):
    fibers = defaultdict(list)
    def rec(k,a,b,sign,key):
        if k == len(Hword):
            fibers[(a,b)].append((key,sign))
            return
        n = abs(Hword[k])
        eps = -1 if Hword[k] < 0 else 1
        for c in cg(a,n):
            rec(k+1,c,b,sign,key+((0,c),))
        for c in cg(b,n):
            rec(k+1,a,c,sign*eps,key+((1,c),))
    rec(0,0,0,1,())
    for v in fibers:
        fibers[v].sort(key=lambda p:p[0])
    return fibers

def direct_phi(word):
    d = {(0,0):1}
    for lab in word:
        n = abs(lab)
        eps = -1 if lab < 0 else 1
        e = defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n):
                e[(c,b)] += v
            for c in cg(b,n):
                e[(a,c)] += eps*v
        d = dict(e)
    return d.get((0,0),0)

def unsigned_path_count(word):
    d = {(0,0):1}
    for lab in word:
        e = defaultdict(int)
        n = abs(lab)
        for (a,b),v in d.items():
            for c in cg(a,n):
                e[(c,b)] += v
            for c in cg(b,n):
                e[(a,c)] += v
        d = dict(e)
    return d.get((0,0),0)

fold_words = fold_atoms = fold_fixed = 0
for hlen in range(5):
    for half_codes in combinations_with_replacement(range(8), hlen):
        Hword = tuple(signed_label(c) for c in half_codes)
        full = Hword + tuple(reversed(Hword))
        fibers = half_paths(Hword)
        fixed_count = atom_pairs = 0
        for endpoint,paths in fibers.items():
            plus = [p for p in paths if p[1] == 1]
            minus = [p for p in paths if p[1] == -1]
            mates = {}
            for P,Q in zip(plus,minus):
                mates[P[0]] = Q[0]
                mates[Q[0]] = P[0]
            assert all(mates[mates[k]] == k for k in mates)
            sign_by_key = {key:sgn for key,sgn in paths}
            assert all(sign_by_key[k] == -sign_by_key[v]
                       for k,v in mates.items())
            unpaired = [(key,sgn) for key,sgn in paths if key not in mates]
            assert len(unpaired) == abs(len(plus)-len(minus))
            if unpaired:
                assert len({sgn for _,sgn in unpaired}) == 1
                assert all(sgn*sgn == 1 for _,sgn in unpaired)
            fixed_count += len(unpaired)**2
            atom_pairs += len(paths)**2
        assert atom_pairs == unsigned_path_count(full)
        assert fixed_count == direct_phi(full)
        fold_words += 1
        fold_atoms += atom_pairs
        fold_fixed += fixed_count

assert qt_cases == 66 and qt_nonzero_b == 20
assert proper_subsets == 494
assert fold_words == 495 and fold_atoms == 515337 and fold_fixed == 257473
print("STR1b SO(3): Schouten, six delta contractions, Psi Gram=2I: PASS")
print("STR1b operator: rank=8b; 4x4 minor determinant=4")
print("STR1b q,t cases:",qt_cases,"nonzero b:",qt_nonzero_b,"sum Phi:",qt_sum_phi)
print("STR1b A_L rows (L, labels, b, beta path, I+, I-, Phi):")
for row in family:
    print(" ",row)
print("STR1b PB proper-subset checks:",proper_subsets,"PASS")
print("STR7 folded words:",fold_words,"path-pairs:",fold_atoms,
      "positive fixed points:",fold_fixed)
print("STR7 path-pair and fixed-point counts match signed fusion: PASS")
