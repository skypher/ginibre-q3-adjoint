import argparse
from collections import defaultdict
from itertools import combinations_with_replacement
from math import comb

ap = argparse.ArgumentParser(
    description="FM-MECH107: exact cone certificates and insertion controls."
)
ap.add_argument("--max-core", type=int, default=10)
ap.add_argument("--box", type=int, default=3)
ap.add_argument("--basis", type=int, default=18)
args = ap.parse_args()

def clean(F):
    return {q: v for q, v in F.items() if v}

def add(F, G, c=1):
    out = dict(F)
    for q, v in G.items():
        out[q] = out.get(q, 0) + c*v
    return clean(out)

def step(F, n, eps=1):
    out = defaultdict(int)
    for (i, j), v in F.items():
        for r in range(abs(i-n), i+n+1, 2):
            out[r, j] += v
        for r in range(abs(j-n), j+n+1, 2):
            out[i, r] += eps*v
    return clean(out)

def word(B):
    F = {(0, 0): 1}
    for z in B:
        F = step(F, abs(z), 1 if z > 0 else -1)
    return F

def anti(i, j):
    if i < 0 or j < 0 or i == j:
        return {}
    return {(i, j): 1, (j, i): -1}

def acone(F):
    return (
        all(F.get((j, i), 0) == -v and (i < j or v >= 0)
            for (i, j), v in F.items())
        and all(i != j for i, j in F)
    )

def P(F):
    return add(step(F, 2), F, -1)

def multi_step(F, labels):
    for n in labels:
        F = step(F, n)
    return F

def columns_nonnegative(F):
    return all(v >= 0 for (i, j), v in F.items() if j == 0)

basis_checks = wide_checks = 0
for i in range(1, args.basis+1):
    for j in range(i):
        A = anti(i, j)
        rhs = {}
        for r, s in ((i+1, j), (i-1, j), (i, j+1), (i, j-1)):
            rhs = add(rhs, anti(r, s))
        assert step(A, 1) == rhs and acone(rhs)
        basis_checks += 1

        if (i-j) % 2 == 0:
            rhs = {}
            for r, s in ((i+2, j), (i-2, j), (i, j+2)):
                rhs = add(rhs, anti(r, s))
            if j >= 2:
                rhs = add(rhs, anti(i, j-2))
            if j:
                rhs = add(rhs, A)
            assert P(A) == rhs and acone(rhs)
            basis_checks += 1

        for ell in (2, 4, 6):
            if i-j >= ell:
                assert acone(add(step(A, ell), A, -1))
                wide_checks += 1
        for ell, m in ((1, 1), (1, 3), (3, 3), (3, 5)):
            if i-j >= ell+m:
                assert acone(add(multi_step(A, (ell, m)), A, -1))
                wide_checks += 1

# Bounded bridge tests; the parameter-uniform proof is algebraic.
plus_lists = [()]
for L in range(1, 4):
    plus_lists += [
        ns for ns in combinations_with_replacement(range(3, 7), L)
        if sum(ns) <= args.max_core
    ]

backgrounds = column_checks = general_insertions = 0
for ns in plus_lists:
    M = sum(ns)
    for n in range(max(3, M), args.max_core+1):
        if n in ns:
            continue
        core = word((-n,) + ns)
        assert acone(core)
        eta = (n+M) % 2

        for a in range(eta, args.box+1):
            F = multi_step(core, (1,)*a)
            for b in range(args.box+1):
                assert columns_nonnegative(F)
                assert columns_nonnegative(P(F))

                B = (
                    (-1,)*a + (2,)*b
                    + tuple(z*((-1)**abs(z)) for z in (-n,) + ns)
                )
                reflected = word(B)
                assert reflected == {
                    (i, j): v*((-1)**j) for (i, j), v in F.items()
                }
                assert columns_nonnegative(P(reflected))
                backgrounds += 2
                column_checks += max(i+j for i, j in F) + 3

                if b <= 1:
                    for ins in ((4,), (6,), (1, 1), (1, 3), (3, 3)):
                        if n >= M+sum(ins) and n not in ins:
                            gain = add(multi_step(F, ins), F, -1)
                            assert columns_nonnegative(gain)
                            general_insertions += 1
                F = step(F, 2)

# All-plus backgrounds: even insertions and odd pairs.
all_plus_checks = 0
for L in range(4):
    for ns in combinations_with_replacement(range(1, 6), L):
        F = word(ns)
        for ins in ((2,), (4,), (6,), (1, 1), (1, 3), (3, 3), (3, 5)):
            gain = add(multi_step(F, ins), F, -1)
            for p in range(max((3,) + ns + ins),
                           sum(ns)+sum(ins)+1):
                assert gain.get((p, 0), 0) >= 0
                all_plus_checks += 1

# Complete one-core corollary, including odd n with a=0.
single_core_checks = 0
for n in range(3, args.max_core+1):
    for a in range(args.box+1):
        for b in range(args.box+1):
            for sig in (-1, 1):
                for eps in (-1, 1):
                    F = word((sig,)*a + (2,)*b + (eps*n,))
                    gain = P(F)
                    for p in range(n, n+a+2*b+3):
                        assert gain.get((p, 0), 0) >= 0
                        single_core_checks += 1

identity_checks = 0
for B in ((-1, 3, 3), (-1, 2, 2, 2, 2), (-4, 6),
          (-3, -3, 4), (1, -2, 3, -4)):
    F = word(B)
    Q = P(F)
    for p in range(2, sum(map(abs, B))+3):
        rhs = (F.get((p-2, 0), 0) + F.get((p+2, 0), 0)
               + F.get((p, 2), 0))
        assert Q.get((p, 0), 0) == rhs
        identity_checks += 1

# The unused two-sided strengthening fails.
F = word((-1, 2, 2, 2, 2))
p = 3
A = F.get((p-2, 0), 0) + F.get((p+2, 0), 0)
C = F.get((p, 2), 0)
assert (A, C, A+C) == (43, 44, 87)

# Actual two-core obstruction to coordinatewise cone positivity.
F = word((-4, 6))
rhs = {}
for r in (2, 4, 6, 8, 10):
    rhs = add(rhs, anti(r, 0))
rhs = add(rhs, anti(6, 4), -1)
assert F == rhs and not acone(F)
two_core_columns = [
    (p, P(F).get((p, 0), 0)) for p in (6, 8, 10, 12)
]
assert two_core_columns == [(6, 2), (8, 2), (10, 1), (12, 1)]
assert P(step(anti(2, 1), 2)).get((3, 0), 0) == -2

def fusion(ns):
    row = {0: 1}
    for n in ns:
        nxt = defaultdict(int)
        for j, v in row.items():
            for k in range(abs(j-n), j+n+1, 2):
                nxt[k] += v
        row = dict(nxt)
    return row

# Defect-table identity and subgroup-mixture obstruction.
table_entries = 0
for N in range(5, 13):
    p = N-2
    for w in range(N+1):
        L = fusion((1,)*w)
        R = fusion((1,)*(N-w))
        h = (
            L.get(2, 0)*R.get(p, 0)
            + L.get(0, 0)*(R.get(p-2, 0)+R.get(p+2, 0))
        )
        expected = {
            0: (N-1)*(N-2)//2,
            2: N-2,
            4: 2
        }.get(w, 0)
        assert h == expected
        table_entries += 1

    h0 = (N-1)*(N-2)//2
    total = h0 + comb(N, 2)*(N-2) + 2*comb(N, 4)
    assert P(word((1,)*N)).get((p, 0), 0) == total
    if N >= 9:
        assert total > 16*h0
        assert 12*(total-16*h0) == (
            (N-1)*(N-2)*(N*N+3*N-90)
        )

print("cone identities:", basis_checks)
print("wide-gap insertion certificates:", wide_checks)
print("dominant-core backgrounds:", backgrounds)
print("column checks:", column_checks)
print("general even/odd-pair insertion bridges:", general_insertions)
print("all-plus insertion checks:", all_plus_checks)
print("complete one-core checks:", single_core_checks)
print("arbitrary-background identities:", identity_checks)
print("two-sided control:", (A, C, A+C))
print("two-core columns:", two_core_columns)
print("defect-table entries:", table_entries)
print("subgroup obstruction N=9,p=7: origin=28, mass=532, capacity=448")
print("PASS")
