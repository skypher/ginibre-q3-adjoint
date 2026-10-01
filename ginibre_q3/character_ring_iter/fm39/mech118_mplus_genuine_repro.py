import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product

ap = argparse.ArgumentParser(
    description="Verify even insertions and the two-core ordering screen"
)
ap.add_argument("--a-max", type=int, default=4)
ap.add_argument("--b-max", type=int, default=3)
ap.add_argument("--core-max", type=int, default=6)
args = ap.parse_args()

def fuse(i, j):
    return range(abs(i-j), i+j+1, 2)

def mul(A, B):
    out = defaultdict(int)
    for (i, j), u in A.items():
        for (k, ell), v in B.items():
            for r in fuse(i, k):
                for s in fuse(j, ell):
                    out[r, s] += u*v
    return {key: value for key, value in out.items() if value}

def add(A, B, scale=1):
    out = defaultdict(int, A)
    for key, value in B.items():
        out[key] += scale*value
    return {key: value for key, value in out.items() if value}

def char_factor(n, eps):
    return {(n, 0): 1, (0, n): eps}

def character(factors):
    out = {(0, 0): 1}
    for n, eps in factors:
        out = mul(out, char_factor(n, eps))
    return out

def insertion_formula(C, p, ell):
    assert ell >= 2 and ell % 2 == 0 and p >= ell
    q = ell // 2
    return sum(
        C.get((p-2*j, 0), 0) + C.get((p+2*j, 0), 0)
        for j in range(1, q+1)
    ) + C.get((p, ell), 0)

# Check (1) on pair-free signed backgrounds.
insert_checks = 0
for a in range(args.a_max + 1):
    for eps1 in (-1, 1):
        for b in range(args.b_max + 1):
            for count in range(3):
                for labels in combinations(range(3, args.core_max + 1), count):
                    for signs in product((-1, 1), repeat=count):
                        cores = tuple(zip(labels, signs))
                        factors = [(1, eps1)] * a + [(2, 1)] * b + list(cores)
                        C = character(factors)
                        weight = a + 2*b + sum(labels)
                        max_label = max([2] + list(labels))
                        for ell in (2, 4, 6):
                            if any(n == ell and eps == -1 for n, eps in cores):
                                continue
                            for p in range(max(3, max_label, ell), weight+ell+1):
                                if any(n == p and eps == -1 for n, eps in cores):
                                    continue
                                after = mul(C, char_factor(ell, 1))
                                actual = after.get((p, 0), 0) - C.get((p, 0), 0)
                                assert actual == insertion_formula(C, p, ell)
                                insert_checks += 1

# K_{a,b}(p;i,j) = E[s^a Z^b (Z-1) A_ij U_p(x)].
def plus(n):
    return {(n, 0): 1, (0, n): 1}

@lru_cache(None)
def qbase(a, b):
    C = {(0, 0): 1}
    for _ in range(a):
        C = mul(C, plus(1))
    for _ in range(b):
        C = mul(C, plus(2))
    out = mul(C, add(plus(2), {(0, 0): 1}, -1))
    return tuple(sorted(out.items()))

@lru_cache(None)
def qtimesp(a, b, p):
    out = mul(dict(qbase(a, b)), {(p, 0): 1})
    return tuple(sorted(out.items()))

def K(a, b, p, i, j):
    P = dict(qtimesp(a, b, p))
    return P.get((i, j), 0) - P.get((j, i), 0)

ordering_checks = 0
zero_slacks = 0
min_slack = None
min_record = None
min_individual_K = None

for m in range(3, 9):
    for n in range(m+1, 10):
        for a in range((m+n) % 2, 8):
            for b in range(6):
                degree = a + 2*b + m + n
                for p in range(n+1, degree+1):
                    if (p-a-m-n) % 2:
                        continue
                    terms = [K(a, b, p, r, 0)
                             for r in range(n-m, n+m+1, 2)]
                    rhs = K(a, b, p, n, m)
                    assert min(terms + [rhs]) >= 0
                    current_min = min(terms + [rhs])
                    min_individual_K = (
                        current_min if min_individual_K is None
                        else min(min_individual_K, current_min)
                    )
                    slack = sum(terms) - rhs
                    assert slack >= 0, (a, b, m, n, p, terms, rhs)
                    ordering_checks += 1
                    if slack == 0:
                        zero_slacks += 1
                    if min_slack is None or slack < min_slack:
                        min_slack = slack
                        min_record = (
                            a, b, m, n, p, sum(terms), rhs, slack
                        )

print("even insertion identity checks:", insert_checks)
print("two-core ordering checks:", ordering_checks)
print("zero slacks:", zero_slacks)
print("minimum individual K:", min_individual_K)
print("minimum ordering slack record:", min_record)
print("PASS")
