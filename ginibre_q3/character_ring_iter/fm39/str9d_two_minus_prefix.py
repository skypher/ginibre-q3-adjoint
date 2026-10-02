import sys
if any(arg in ("-h", "--help") for arg in sys.argv[1:]):
    print("FM-STR9d verifier: exact SU(2)xSU(2) character tables and two-minus prefixes.")
    raise SystemExit(0)

from datetime import datetime, timezone
from itertools import combinations_with_replacement

def log(*items):
    print(datetime.now(timezone.utc).isoformat(timespec="seconds"), *items, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def multiply_factor(table, n, eps=1):
    out = {}
    for (r, s), value in table.items():
        for t in cg(r, n):
            out[t, s] = out.get((t, s), 0) + value
        for t in cg(s, n):
            out[r, t] = out.get((r, t), 0) + eps*value
    return {key: value for key, value in out.items() if value}

def character(factors):
    table = {(0, 0): 1}
    for n, eps in factors:
        table = multiply_factor(table, n, eps)
    return table

def interior_split(plus_labels):
    ordered = sorted(enumerate(plus_labels), key=lambda z: (-z[1], z[0]))
    A, B = [], []
    weight_A = weight_B = 0
    for _, n in ordered:
        if weight_A <= weight_B:
            A.append((n, 1))
            weight_A += n
        else:
            B.append((n, 1))
            weight_B += n
    if weight_A > weight_B:
        A, B = B, A
    return A, B

def pair_prefix(plus_labels, n, m):
    A, B = interior_split(plus_labels)
    Bprime = B + [(n, -1), (m, -1)]
    fa, fb = character(A), character(Bprime)
    degree_A = max((r+s for r, s in fa), default=0)
    levels = []
    for height in range(degree_A + 1):
        levels.append(sum(value*fb.get((r, s), 0)
                          for (r, s), value in fa.items()
                          if r+s == height))
    prefixes, total = [], 0
    for value in levels:
        total += value
        prefixes.append(total)
    full = character([(n, -1), (m, -1)] + [(q, 1) for q in plus_labels])
    phi = full.get((0, 0), 0)
    assert prefixes[-1] == phi
    return prefixes, phi, A, B

def add_tables(*terms):
    out = {}
    for coefficient, table in terms:
        for key, value in table.items():
            out[key] = out.get(key, 0) + coefficient*value
    return {key: value for key, value in out.items() if value}

log("checking Sp(4) factor identities")
one = {(0, 0): 1}
h2 = {(2, 0): 1, (1, 1): 1, (0, 2): 1}
chi11 = {(0, 0): 1, (1, 1): 1}
s2 = character([(2, 1)])
assert add_tables((1, h2), (-1, chi11), (1, one)) == s2
d_s2 = multiply_factor(s2, 1, -1)
assert d_s2 == {(3, 0): 1, (0, 3): -1,
                (1, 0): 1, (0, 1): -1,
                (2, 1): -1, (1, 2): 1}
for n in range(1, 31):
    K_n = {(n-1-j, j): 1 for j in range(n)}
    assert multiply_factor(K_n, 1, -1) == character([(n, -1)])
log("PASS: S2=chi(2,0)-chi(1,1)+chi(0,0); d K_n=U_n(x)-U_n(y), n<=30")

log("checking the uniform (+1)^a two-minus class")
class_checks = 0
for a in range(17):
    for n in range(1, 13):
        for m in range(n, 13):
            if a and (n == 1 or m == 1):
                continue
            prefixes, phi, A, B = pair_prefix([1]*a, n, m)
            assert all(value >= 0 for value in prefixes)
            class_checks += 1
    log("fundamental suffix length", a, "cases", class_checks)
assert class_checks == 1134
log("PASS fundamental-background screen", class_checks)

log("checking arbitrary plus labels on a bounded box")
screen_checks = 0
for length in range(7):
    for plus_labels in combinations_with_replacement(range(1, 8), length):
        for n in range(1, 8):
            for m in range(n, 8):
                if n in plus_labels or m in plus_labels:
                    continue
                prefixes, phi, A, B = pair_prefix(list(plus_labels), n, m)
                assert all(value >= 0 for value in prefixes)
                screen_checks += 1
    log("arbitrary-plus length", length, "cumulative cases", screen_checks)
assert screen_checks == 16170
log("PASS bounded arbitrary-plus screen", screen_checks)
log("PASS FM-STR9d exact verifier")
