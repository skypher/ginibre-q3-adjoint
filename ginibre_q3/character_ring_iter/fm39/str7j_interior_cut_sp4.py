from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import time

def stamp():
    return time.strftime("%H:%M:%S", time.gmtime())

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def mul_factor(tab, z):
    n = abs(z)
    eps = 1 if z > 0 else -1
    out = defaultdict(int)
    for (r, s), value in tab.items():
        for t in cg(r, n):
            out[t, s] += value
        for t in cg(s, n):
            out[r, t] += eps * value
    return {k: v for k, v in out.items() if v}

@lru_cache(None)
def table_cached(word):
    tab = {(0, 0): 1}
    for z in word:
        tab = mul_factor(tab, z)
    return tab

def table(word):
    return table_cached(tuple(sorted(word, key=lambda z: (-abs(z), z))))

def interior_cut(word):
    A, B = [], []
    wa = wb = 0
    for z in sorted(word, key=lambda z: (-abs(z), z)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    return A, B

def pair_prefix(L, i, j):
    C = [z for k, z in enumerate(L) if k not in (i, j)]
    A, B = interior_cut(C)
    Bp = B + [L[i], L[j]]
    FA, FB = table(A), table(Bp)
    layer = defaultdict(int)
    for (r, s), value in FA.items():
        layer[r+s] += value * FB.get((r, s), 0)
    running = 0
    minimum = None
    argmin = None
    for height in sorted(layer):
        running += layer[height]
        if minimum is None or running < minimum:
            minimum, argmin = running, height
        assert running >= 0, (L, i, j, height, running)
    imbalance = abs(sum(map(abs, A)) - sum(map(abs, Bp)))
    return minimum, argmin, running, A, Bp, imbalance

def audit(name, L):
    start = time.time()
    count = 0
    least = None
    where = None
    max_imbalance = 0
    for i, j in combinations(range(len(L)), 2):
        minimum, height, phi, A, Bp, imbalance = pair_prefix(L, i, j)
        count += 1
        max_imbalance = max(max_imbalance, imbalance)
        if least is None or minimum < least:
            least = minimum
            where = (L[i], L[j], height, tuple(A), tuple(Bp), phi)
    print(stamp(), name, "pairs", count, "minimum_prefix", least,
          "at", where, "max_weight_imbalance", max_imbalance,
          "seconds", round(time.time()-start, 2), flush=True)
    return count

def multiply_tables(A, B):
    out = defaultdict(int)
    for (r, s), va in A.items():
        for (u, v), vb in B.items():
            for x in cg(r, u):
                for y in cg(s, v):
                    out[x, y] += va * vb
    return {k: value for k, value in out.items() if value}

def multiply_d(tab):
    out = defaultdict(int)
    for (r, s), value in tab.items():
        for x in cg(r, 1):
            out[x, s] += value
        for y in cg(s, 1):
            out[r, y] -= value
    return {k: value for k, value in out.items() if value}

def quotient_table(word):
    # K_n = sum_{j=0}^{n-1} U_j(x) U_{n-1-j}(y);
    # S_n = U_n(x) + U_n(y).
    q = {(0, 0): 1}
    for z in word:
        n = abs(z)
        factor = ({(j, n-1-j): 1 for j in range(n)} if z < 0
                  else {(n, 0): 1, (0, n): 1})
        q = multiply_tables(q, factor)
    return q

def character_coefficients(symmetric_tab):
    # d*chi_(alpha,beta) = U_(alpha+1)(x)U_beta(y)
    #                         - U_beta(x)U_(alpha+1)(y).
    alternating = multiply_d(symmetric_tab)
    return {(i-1, j): value for (i, j), value in alternating.items()
            if i > j}

def top_pair(B):
    best = None
    for i, j in combinations(range(len(B)), 2):
        if (abs(B[i]) - abs(B[j])) % 2:
            continue
        key = (abs(B[i]) + abs(B[j]), max(abs(B[i]), abs(B[j])))
        if best is None or key > best[0]:
            best = (key, i, j)
    assert best is not None
    return best[1], best[2]

# Exact segregated HPP failure from FM-SEC181.
bad_A = [1, -2, -4, -4, -4]
bad_B = [1] + [3]*6 + [5]*4 + [6]
TA, TB = table(bad_A), table(bad_B)
bad_layers = defaultdict(int)
for (r, s), value in TA.items():
    bad_layers[r+s] += value * TB.get((r, s), 0)
pi11 = sum(value for h, value in bad_layers.items() if h <= 11)
phi_bad = sum(bad_layers.values())
assert (pi11, phi_bad) == (-24695910, 120541550)
print(stamp(), "segregated Pi_11", pi11, "Phi", phi_bad, flush=True)

words = [
    ("SEC181-W1", [1,1,-2] + [3]*6 + [-4]*3 + [5]*4 + [6]),
    ("SEC181-W2", [1,1,-2] + [3]*5 + [-4]*3 + [5]*5 + [6]),
    ("SEC181-W3", [1,1,1,-2,2] + [3]*5 + [-4]*3 + [4]*2 + [5]*2 + [8]),
]
checks = 0
for name, word in words:
    checks += audit(name, word)

F1 = [1,1,-2] + [3]*4 + [-4]*3 + [5]*4 + [8]
checks += audit("F1", F1)
tight = [-1,-3,-3,-4,-5,-5,-6,-7]
checks += audit("census-tight", tight)

L = words[0][1]
i, j = top_pair(L[:-1])
A, B = interior_cut([z for k, z in enumerate(L) if k not in (i, j)])
Bp = B + [L[i], L[j]]
a, b = sum(z < 0 for z in A), sum(z < 0 for z in Bp)
assert (L[i], L[j], tuple(A), tuple(Bp), a, b) == (
    5,5,(6,-4,-4,3,3,3,-2),(5,5,-4,3,3,3,1,1,5,5),3,1)
topdata = pair_prefix(L, i, j)
assert (topdata[0], topdata[1], topdata[2], topdata[5]) == (
    3938032, 1, 120541550, 10)
GA = multiply_d(multiply_d(quotient_table(A)))
GB = quotient_table(Bp)
ca, cb = character_coefficients(GA), character_coefficients(GB)
assert (ca[(5,5)], cb[(16,16)]) == (-433, -12)
print(stamp(), "W1 TopPair prefix min", topdata[0], "at height", topdata[1],
      "weights", sum(map(abs, A)), sum(map(abs, Bp)),
      "reserve r=(1,1); chi_(5,5)=-433; chi_(16,16)=-12", flush=True)

run_cases = 0
for k in range(2, 11):
    sign = (-1)**k
    for p0 in range(k+1, k+4):
        if (k*(k+1)//2 + p0) % 2:
            continue
        run = [-n for n in range(1, k+1)] + [sign*p0]
        checks += audit(f"run-k{k}-p{sign*p0}", run)
        run_cases += 1

print(stamp(), "exact pair-prefix checks", checks, "run cases", run_cases,
      "table cache entries", table_cached.cache_info().currsize, "PASS",
      flush=True)
