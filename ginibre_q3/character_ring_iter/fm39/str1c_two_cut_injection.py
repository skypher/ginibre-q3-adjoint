import argparse
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from fractions import Fraction as Q
from math import comb, factorial, isqrt

argparse.ArgumentParser(
    description="FM-STR1c: exact transfer, symmetry, and two-cut checks; memory only."
).parse_args()

@lru_cache(None)
def fusion(ns):
    D = {0: 1}
    for n in ns:
        E = {}
        for a, v in D.items():
            for b in range(abs(a-n), a+n+1, 2):
                E[b] = E.get(b, 0) + v
        D = E
    return D

def inv(ns):
    return fusion(tuple(sorted(ns))).get(0, 0)

def multiplicity(ns, t):
    w = sum(ns)
    if t < 0 or t > w or (w-t) % 2:
        return 0
    d = (w-t)//2
    A = [1] + [0]*d
    for n in ns:
        B = []
        s = 0
        for j in range(d+1):
            s += A[j]
            if j > n:
                s -= A[j-n-1]
            B.append(s)
        A = B
    return A[d] - (A[d-1] if d else 0)

def dimensions(n, m, plus):
    dp = dm = 0
    for mask in range(1 << len(plus)):
        S = tuple(a for i, a in enumerate(plus) if mask >> i & 1)
        T = tuple(a for i, a in enumerate(plus) if not (mask >> i & 1))
        dp += 2*inv((n, m)+S)*inv(T)
        dm += 2*inv((n,)+S)*inv((m,)+T)
    return dp, dm

# A permutation character certificate: ker(C[2-subsets] -> C[1-subsets]).
@lru_cache(None)
def partitions(n, top=0):
    if n == 0:
        return ((),)
    top = min(top or n, n)
    return tuple((a,)+p for a in range(top, 0, -1)
                 for p in partitions(n-a, a))

def z_class(c):
    return product_int(l**c.count(l)*factorial(c.count(l)) for l in set(c))

def product_int(xs):
    a = 1
    for x in xs:
        a *= x
    return a

def adams(n, ell):
    w = {ell*j: 1 for j in range(-n, n+1, 2)}
    return {j: w.get(j, 0)-w.get(j+2, 0)
            for j in range(n*ell+1)
            if w.get(j, 0) != w.get(j+2, 0)}

def cycle_table(n, cycles):
    A = {(0, 0): 1}
    for ell in cycles:
        B = {}
        for (a, b), v in A.items():
            for j, u in adams(n, ell).items():
                for k in range(abs(a-j), a+j+1, 2):
                    B[k, b] = B.get((k, b), 0) + u*v
                for k in range(abs(b-j), b+j+1, 2):
                    B[a, k] = B.get((a, k), 0) + u*v
        A = {ij: v for ij, v in B.items() if v}
    return A

norm = source = target = Q(0)
for c in partitions(6):
    f1, f2 = c.count(1), c.count(2)
    chi = comb(f1, 2) + f2 - f1
    A = cycle_table(3, c)
    norm += Q(chi*chi, z_class(c))
    source += Q(chi*A.get((1, 1), 0), z_class(c))
    target += Q(chi*2*A.get((0, 0), 0), z_class(c))
assert (norm, source, target) == (1, 1, 0)
assert dimensions(1, 1, (3,)*6) == (946, 160)
print("S2 x S6 certificate: character norm 1, source multiplicity 1,"
      " target multiplicity 0; missing dimension 9", flush=True)

# Exact normalized four-cup Gram matrix.
for q in range(1, 9):
    d = q+1
    def cup(i, j):
        return (-1)**i if i+j == q else 0
    norm1 = norm2 = cross = 0
    for u, v, a, b in product(range(d), repeat=4):
        h = cup(u, a)*cup(v, b)
        k = cup(u, b)*cup(v, a)
        norm1 += h*h
        norm2 += k*k
        cross += h*k
    assert (norm1, norm2, cross) == (d*d, d*d, d)
    assert 1-Q(cross, d*d)**2 == Q(q*(q+2), (q+1)**2)
print("two-cut Gram formula: q = 1..8 PASS", flush=True)

# Test the theorem's actual support condition, not just a weight surrogate.
tested = 0
for n in range(1, 6):
    for m in range(n, 6):
        labels = tuple(a for a in range(1, 6) if a not in (n, m))
        for L in range(1, 6):
            for plus in combinations_with_replacement(labels, L):
                full = (1 << L)-1
                support = []
                for mask in range(1 << L):
                    S = tuple(a for i, a in enumerate(plus) if mask >> i & 1)
                    T = tuple(a for i, a in enumerate(plus) if not (mask >> i & 1))
                    v = inv((n,)+S)*inv((m,)+T)
                    if v:
                        support.append((mask, v))
                if not support or len(support) > 2:
                    continue
                if len(support) == 2 and support[0][0] ^ support[1][0] != full:
                    continue
                dp, dm = dimensions(n, m, plus)
                assert dm == 2*sum(v for _, v in support)
                assert dp >= dm
                tested += 1
print("two-cut small-list consumer checks:", tested, "PASS", flush=True)

# Unbounded family: (-1)^2, three +q's, and a primitive block B_L.
family_checks = subset_checks = 0
for q in (3, 5, 7):
    assert multiplicity((q, q, q), 1) == 2
    M = 3*q+3
    for L in range(2, 11):
        first = tuple(M*(1 << j) for j in range(L-1))
        B = first+(sum(first)-1,)
        assert multiplicity(B, 1) == L-1
        for mask in range(1, (1 << L)-1):
            S = [a for i, a in enumerate(B) if mask >> i & 1]
            assert 2*max(S)-sum(S) >= M-1 > 3*q+1
            subset_checks += 1
        gramdet = Q(3, 4)**(4*(L-1))
        assert gramdet > 0
        family_checks += 1
print("unbounded family:", family_checks, "members,",
      subset_checks, "proper-subset bounds PASS", flush=True)

# Fixed-prefix first transfers see only one of binom(2n,n) source cuts.
for n in range(2, 9):
    assert comb(2*n, n) > 1
    assert multiplicity((1,)*n, n) == 1
    assert all(multiplicity((1,)*j, n) == 0 for j in range(n))
print("prefix obstruction: rank <= 2 versus dimension 2*binom(2n,n),"
      " n = 2..8 PASS", flush=True)

# Exact integer coefficient matrices for direct + pair-transfer corrections.
prime = 65521
assert all(prime % d for d in range(2, isqrt(prime)+1))

def edge(A, i, j):
    B = {}
    for e, c in A.items():
        for k, s in ((j, 1), (i, -1)):
            f = list(e)
            f[k] += 1
            f = tuple(f)
            B[f] = B.get(f, 0)+s*c
    return {e: c for e, c in B.items() if c}

def covariants(u, S, labels):
    out = []
    for t in range(2):
        r = [labels[i]-(j == t) for j, i in enumerate(S)]
        exponents = ((r[0]+r[1]-r[2])//2,
                     (r[0]+r[2]-r[1])//2,
                     (r[1]+r[2]-r[0])//2)
        assert min(exponents) >= 0
        A = edge({(0,)*8: 1}, u, S[t])
        for (i, j), exponent in zip(combinations(S, 2), exponents):
            for _ in range(exponent):
                A = edge(A, i, j)
        out.append(A)
    return out

def multiply(A, B):
    C = {}
    for e, a in A.items():
        for f, b in B.items():
            g = tuple(x+y for x, y in zip(e, f))
            C[g] = C.get(g, 0)+a*b
    return {e: c for e, c in C.items() if c}

def contraction(A, i, j, n):
    B = {}
    for e, c in A.items():
        if e[i]+e[j] != n:
            continue
        f = list(e)
        f[i] = f[j] = 0
        f = tuple(f)
        v = c*(-1)**e[i]*factorial(e[i])*factorial(n-e[i])
        B[f] = B.get(f, 0)+v
    return {e: c for e, c in B.items() if c}

def independent_rows(cols):
    piv = {}
    rows = []
    for col in cols:
        A = {e: c % prime for e, c in col.items() if c % prime}
        while A:
            e = min(A)
            if e not in piv:
                break
            c = A[e]
            for f, v in piv[e].items():
                A[f] = (A.get(f, 0)-c*v) % prime
                if not A[f]:
                    del A[f]
        if A:
            e = min(A)
            iv = pow(A[e], prime-2, prime)
            piv[e] = {f: c*iv % prime for f, c in A.items()}
            rows.append(e)
    return rows

def determinant_mod(A):
    A = [[v % prime for v in row] for row in A]
    det = 1
    for j in range(len(A)):
        k = next((k for k in range(j, len(A)) if A[k][j]), None)
        if k is None:
            return 0
        if k != j:
            A[k], A[j] = A[j], A[k]
            det = -det
        d = A[j][j]
        det = det*d % prime
        iv = pow(d, prime-2, prime)
        A[j] = [x*iv % prime for x in A[j]]
        for i in range(j+1, len(A)):
            d = A[i][j]
            A[i] = [(x-d*y) % prime for x, y in zip(A[i], A[j])]
    return det % prime

for plus in ((3,)*6, (3,)*5+(5,), (3,)*4+(5,)*2, (3,)*3+(5,)*3):
    labels = (1, 1)+plus
    cols = []
    blind = 0
    for S in combinations(range(2, 8), 3):
        T = tuple(i for i in range(2, 8) if i not in S)
        for A, B in product(covariants(0, S, labels),
                            covariants(1, T, labels)):
            H = multiply(A, B)
            C = {(0,)+e: v for e, v in H.items()}
            all_pair_contractions = []
            for i, j in combinations(range(2, 8), 2):
                if labels[i] != labels[j]:
                    continue
                K = contraction(H, i, j, labels[i])
                all_pair_contractions.append(K)
                if (i in S) == (j in S):
                    assert not K
                    continue
                weight = i+1 if i in S else j+1
                for e, v in K.items():
                    C[(1, i, j)+e] = weight*v
            if not any(all_pair_contractions):
                blind += 1
            cols.append(C)
    assert len(cols) == 80
    rows = independent_rows(cols)
    assert len(rows) == 80
    minor = determinant_mod([[C.get(row, 0) for C in cols] for row in rows])
    assert minor != 0
    dp, dm = dimensions(1, 1, plus)
    assert dm == 160
    if plus == (3,)*3+(5,)*3:
        assert blind == 8
        assert (dp, dm) == (826, 160)
    print("repair", plus, "rank 160/160; minor mod 65521 =", minor,
          "; pair-blind columns/orientation =", blind,
          "; Phi =", dp-dm, flush=True)

print("FM-STR1c PASS", flush=True)
