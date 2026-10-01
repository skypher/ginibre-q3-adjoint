from collections import Counter, defaultdict
from itertools import product
from math import comb
import random

def U(n):
    return {n - 2*j: (-1)**j * comb(n-j, j)
            for j in range(n//2 + 1)}

# Character-basis walk recurrence used for A_ab.
def entry(L, A, B):
    L = sorted(L, key=lambda z: -abs(z))
    rem = sum(map(abs, L))
    cur = {(0, 0): 1}
    for x in L:
        n = abs(x)
        positive = x > 0
        rem -= n
        nxt = {}
        for (s, t), v in cur.items():
            if positive:
                if abs(t-B) <= rem:
                    for z in range(abs(s-n), s+n+1, 2):
                        if abs(z-A) <= rem:
                            nxt[(z, t)] = nxt.get((z, t), 0) + v
            elif abs(s-A) <= rem:
                for z in range(abs(t-n), t+n+1, 2):
                    if abs(z-B) <= rem:
                        nxt[(s, z)] = nxt.get((s, z), 0) - v
        cur = nxt
    return cur.get((A, B), 0)

def cat(k):
    return comb(2*k, k) // (k+1)

# Independent monomial/Catalan evaluation of A_ab.
def moment_A(C, a, b):
    P = {(0, 0): 1}
    for x in C:
        Q = defaultdict(int)
        for (i, j), v in P.items():
            for d, c in U(abs(x)).items():
                key = (i+d, j) if x > 0 else (i, j+d)
                Q[key] += v*c*(1 if x > 0 else -1)
        P = {k: v for k, v in Q.items() if v}
    for axis, n in ((0, a), (1, b)):
        Q = defaultdict(int)
        for (i, j), v in P.items():
            for d, c in U(n).items():
                key = (i+d, j) if axis == 0 else (i, j+d)
                Q[key] += v*c
        P = {k: v for k, v in Q.items() if v}
    return sum(v*cat(i//2)*cat(j//2)
               for (i, j), v in P.items()
               if i % 2 == 0 and j % 2 == 0)

def scan(L, mode):
    vals = list(dict.fromkeys(L))
    for i, a in enumerate(vals):
        for j in range(i, len(vals)):
            b = vals[j]
            if i == j and L.count(a) < 2:
                continue
            ca, cb = L.count(a), L.count(b)
            preserves = (ca == 2 if a == b else ca == 1 and cb == 1)
            if mode == "keep" and not preserves:
                continue
            if mode == "create" and preserves:
                continue
            C = L.copy()
            C.remove(a)
            C.remove(b)
            v = entry(C, abs(a), abs(b))
            if b < 0:
                v = -v
            if ((mode == "keep" and v > 0)
                or (mode == "create" and v >= 0)
                or (mode == "all" and v >= 0)):
                return a, b, v, preserves
    return None

def residual(B, p):
    W = sum(map(abs, B))
    mx = max(map(abs, B))
    delta = (W-p)//2
    sig = (-1)**sum(x < 0 for x in B)
    return (p >= max(6, mx) and W >= p and (W-p) % 2 == 0
            and delta >= 8 and mx <= delta
            and sum(abs(x) >= 3 for x in B) >= 2
            and not any(abs(x) == p and x == -sig*p for x in B))

profiles = [
    ("B56", [2]*8 + [3]*8 + [4]*4),
    ("B102", [2]*9 + [3]*10 + [5]*6 + [6]*4),
    ("B182", [2]*8 + [3]*10 + [5]*8 + [8]*12),
    ("B292", [2]*8 + [3]*8 + [5]*8 + [8]*10 + [12]*11),
    ("B230_N80", [2]*30 + [3]*30 + [4]*20),
    ("B280_N80", [2]*20 + [3]*20 + [4]*20 + [5]*20),
    ("B360_N80", [3]*20 + [4]*20 + [5]*20 + [6]*20),
]

summaries = []
positive = None
for name, mag in profiles:
    W = sum(mag)
    delta = max(8, max(mag))
    p = W - 2*delta
    labels = sorted(set(mag))
    counts, first_edges = Counter(), Counter()
    for signs in product((1, -1), repeat=len(labels)):
        sign = dict(zip(labels, signs))
        B = [sign[n]*n for n in mag]
        sigma = (-1)**sum(x < 0 for x in B)
        L = B + [sigma*p]
        assert residual(B, p)
        assert scan(L, "keep") is None
        e = scan(L, "all")
        if e is None:
            counts["no_nonnegative_flip"] += 1
            continue
        a, b, v, preserves = e
        assert not preserves
        counts["paircreating_nonnegative"] += 1
        counts["zero" if v == 0 else "positive"] += 1
        first_edges[(a, b, v)] += 1
        if name == "B292" and signs == (1, -1, -1, -1, -1):
            C = L.copy()
            C.remove(a)
            C.remove(b)
            assert (a, b, v) == (2, 2, 26942364)
            assert moment_A(C, abs(a), abs(b)) == v
            positive = (B, sigma*p, a, b, v)
    summaries.append((name, len(mag), W, p, sum(first_edges.values()),
                      dict(counts), first_edges.most_common(5)))

# Random local starts on three consecutive/mixed profiles.
def flip(L, a, b):
    Q = L.copy()
    Q.remove(a)
    Q.remove(b)
    return Q + [-a, -b]

def local_trials(mag, p, rng):
    stats = Counter()
    for trial in range(4):
        labels = sorted(set(map(abs, mag)))
        signs = {n: (1 if rng.randrange(2) else -1) for n in labels}
        B = [signs[abs(n)]*abs(n) for n in mag]
        sigma = (-1)**sum(x < 0 for x in B)
        L = B + [sigma*p]
        moves = 0
        while moves < 30:
            e = scan(L, "keep")
            if e is None:
                break
            L = flip(L, e[0], e[1])
            moves += 1
        ec, ea = scan(L, "create"), scan(L, "all")
        stats["starts"] += 1
        stats["positive_preserving_moves"] += moves
        stats["paircreating_nonnegative"] += ec is not None
        stats["all_pair_nonnegative"] += ea is not None
        stats["preserving_nonnegative"] += ea is not None and ea[3]
        assert ea is not None
        if trial == 0:
            a, b, v, _ = ea
            C = L.copy()
            C.remove(a)
            C.remove(b)
            assert moment_A(C, abs(a), abs(b)) == (v if b > 0 else -v)
    return dict(stats)

rng = random.Random(16969)
runs = [
    ("run14", list(range(2, 15)) + [2]*5, 86),
    ("run20", list(range(2, 21)), 169),
    ("run24", list(range(2, 25)), 251),
]
run_results = []
for name, mag, p in runs:
    run_results.append((name, len(mag), sum(mag), p,
                        local_trials(mag, p, rng)))

# Six starts with two copies of each label 2,...,19.
mag = [n for n in range(2, 20) for _ in range(2)]
p = sum(mag) - 38
rng = random.Random(169204)
twocopy = Counter()
for _ in range(6):
    signs = {n: (1 if rng.randrange(2) else -1) for n in range(2, 20)}
    B = [signs[abs(n)]*abs(n) for n in mag]
    sigma = (-1)**sum(x < 0 for x in B)
    L = B + [sigma*p]
    e = scan(L, "create")
    ea = scan(L, "all")
    assert e is not None and ea is not None
    twocopy["starts"] += 1
    twocopy["paircreating_nonnegative"] += 1
    if twocopy["starts"] == 1:
        a, b, v, _ = e
        C = L.copy()
        C.remove(a)
        C.remove(b)
        assert moment_A(C, abs(a), abs(b)) == (v if b > 0 else -v)

# Independent exact checks: an all-positive run and an N=80 profile.
B20 = list(range(2, 21))
p20 = 169
L20 = B20 + [p20]
e20 = scan(L20, "all")
assert e20 == (2, 3, 0, True)
C20 = L20.copy()
C20.remove(2)
C20.remove(3)
assert entry(B20, p20, 0) == 11612114573
assert moment_A(B20, p20, 0) == 11612114573
assert moment_A(C20, 2, 3) == 0

B80 = [2]*30 + [3]*30 + [4]*20
L80 = B80 + [214]
e80 = scan(L80, "all")
assert e80 == (2, 2, 0, False)
C80 = L80.copy()
C80.remove(2)
C80.remove(2)
assert moment_A(C80, 2, 2) == 0

print("BLOCK_PROFILES")
for row in summaries:
    print(row)
print("RANDOM_LOCAL", run_results)
print("TWO_COPY_RUN", len(mag), sum(mag), p, dict(twocopy))
print("POSITIVE_EDGE", positive)
print("ALLPLUS_ZERO", len(B20), sum(B20), p20, e20, "g", 11612114573)
print("N80_ZERO", len(B80), sum(B80), 214, e80)
