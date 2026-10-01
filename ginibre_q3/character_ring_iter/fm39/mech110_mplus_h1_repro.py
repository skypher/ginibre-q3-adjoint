import argparse
from collections import defaultdict
from itertools import combinations, product

parser = argparse.ArgumentParser(description="Verify pair-free H=1 insertion identities")
parser.add_argument("--max-a", type=int, default=8)
parser.add_argument("--max-b", type=int, default=6)
args = parser.parse_args()

def fuse(m, n):
    return range(abs(m - n), m + n + 1, 2)

def multiply(A, B):
    C = defaultdict(int)
    for (i, j), x in A.items():
        for (k, ell), y in B.items():
            for r in fuse(i, k):
                for s in fuse(j, ell):
                    C[r, s] += x * y
    return {key: value for key, value in C.items() if value}

def factor(n, eps):
    return {(n, 0): 1, (0, n): eps}

def character(factors):
    C = {(0, 0): 1}
    for n, eps in factors:
        C = multiply(C, factor(n, eps))
    return C

def reflect_y(C):
    return {(i, j): ((-1) ** j) * c for (i, j), c in C.items() if c}

def channels(C, p):
    return C.get((p - 2, 0), 0), C.get((p + 2, 0), 0), C.get((p, 2), 0)

checks = 0
for a in range(args.max_a + 1):
    for b in range(args.max_b + 1):
        plus = character([(1, 1)] * a + [(2, 1)] * b)
        minus = character([(1, -1)] * a + [(2, 1)] * b)
        assert reflect_y(minus) == plus
        assert all(c >= 0 for c in plus.values())
        plus2 = multiply(plus, factor(2, 1))
        for p in range(3, a + 2 * b + 5):
            inc = plus2.get((p, 0), 0) - plus.get((p, 0), 0)
            assert inc == sum(channels(plus, p)) >= 0
            checks += 1

core_checks = 0
core_lists = [(), (3,), (4,), (3, 4), (3, 5, 6)]
for cores in core_lists:
    for a in range(4):
        for b in range(3):
            plus = character(
                [(1, 1)] * a + [(2, 1)] * b + [(n, 1) for n in cores]
            )
            reflected = character(
                [(1, -1)] * a
                + [(2, 1)] * b
                + [(n, (-1) ** n) for n in cores]
            )
            assert reflect_y(reflected) == plus
            assert all(c >= 0 for c in plus.values())
            plus2 = multiply(plus, factor(2, 1))
            top = a + 2 * b + sum(cores) + 4
            for p in range(3, top + 1):
                inc = plus2.get((p, 0), 0) - plus.get((p, 0), 0)
                assert inc == sum(channels(plus, p)) >= 0
                core_checks += 1

witness = character([(1, 1), (1, 1), (3, -1), (4, -1)])
inserted = multiply(witness, factor(2, 1))
witness_tri = channels(witness, 5)
assert witness_tri == (5, 4, -1)
assert witness.get((5, 0), 0) == 5 and inserted.get((5, 0), 0) == 13

first_bad_channel = None
first_bad_increment = None
census = 0
for weight in range(10):
    for a in range(weight + 1):
        for b in range(weight // 2 + 1):
            for q in range(weight // 3 + 1):
                for ns in combinations(range(3, weight + 1), q):
                    if a + 2 * b + sum(ns) != weight:
                        continue
                    for signs in product((-1, 1), repeat=q):
                        cores = tuple(zip(ns, signs))
                        for eps1 in (-1, 1):
                            nminus = (a if eps1 == -1 else 0) + sum(
                                s == -1 for s in signs
                            )
                            if nminus % 2:
                                continue
                            C = character(
                                [(1, eps1)] * a + list(cores) + [(2, 1)] * b
                            )
                            max_label = max([2] + list(ns))
                            for p in range(max(3, max_label), weight + 3):
                                if any(n == p and s == -1 for n, s in cores):
                                    continue
                                if (p - weight) % 2:
                                    continue
                                tri = channels(C, p)
                                inc = sum(tri)
                                census += 1
                                rec = (weight, a, eps1, cores, b, p, tri, inc)
                                if min(tri) < 0 and first_bad_channel is None:
                                    first_bad_channel = rec
                                if inc < 0 and first_bad_increment is None:
                                    first_bad_increment = rec

expected = (9, 2, 1, ((3, -1), (4, -1)), 0, 5, (5, 4, -1), 8)
assert first_bad_channel == expected
assert first_bad_increment is None

print("H=1 insertion/reflection checks:", checks)
print("character-valued core extension checks:", core_checks)
print("pair-free even-minus census entries through weight 9:", census)
print("signed-core witness (g_before, channels, g_after):",
      witness.get((5, 0), 0), witness_tri, inserted.get((5, 0), 0))
print("first negative channel in census:", first_bad_channel)
print("negative insertion increments in this census:", first_bad_increment)
print("PASS")
