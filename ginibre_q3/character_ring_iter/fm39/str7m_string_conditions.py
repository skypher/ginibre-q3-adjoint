from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations
from datetime import datetime, timezone

def stamp(msg):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
          msg, flush=True)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def mul(P, Q):
    out = defaultdict(int)
    for (a, b), x in P.items():
        for (c, d), y in Q.items():
            for i in cg(a, c):
                for j in cg(b, d):
                    out[i, j] += x*y
    return {k: v for k, v in out.items() if v}

def dmul(P):
    out = defaultdict(int)
    for (i, j), v in P.items():
        out[i+1, j] += v
        if i:
            out[i-1, j] += v
        out[i, j+1] -= v
        if j:
            out[i, j-1] -= v
    return {k: v for k, v in out.items() if v}

@lru_cache(None)
def quotient(word):
    P = {(0, 0): 1}
    for z in word:
        n = abs(z)
        factor = ({(j, n-1-j): 1 for j in range(n)} if z < 0
                  else {(n, 0): 1, (0, n): 1})
        P = mul(P, factor)
    return P

def spvec(word, extra_d):
    P = quotient(tuple(word))
    for _ in range(extra_d):
        P = dmul(P)
    D = dmul(P)
    assert all(D.get((b, a), 0) == -v for (a, b), v in D.items())
    return {(a-1, b): v for (a, b), v in D.items() if a > b and v}

def cut(C):
    A, B = [], []
    wa = wb = 0
    for z in sorted(C, key=lambda z: (-abs(z), z)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    key = lambda w: tuple(sorted(w, key=lambda z: (-abs(z), z)))
    return key(A), key(B)

def string_flags(A, B):
    ga, gb = defaultdict(dict), defaultdict(dict)
    for (a, b), v in A.items():
        ga[a-b][b] = v
    for (a, b), v in B.items():
        gb[a-b][b] = v

    bad1, bad2 = [], []
    for ell, row in ga.items():
        prev = 0
        for j in range(max(row, default=-1) + 1):
            cur = row.get(j, 0)
            if cur < prev:
                bad1.append((ell, j, prev-cur))
            prev = cur

    for ell, row in gb.items():
        tail = sum(row.values())
        if tail < 0:
            bad2.append((ell, 0, -tail))
        for J in range(1, max(row, default=-1) + 1):
            tail -= row.get(J-1, 0)
            if tail < 0:
                bad2.append((ell, J, -tail))
    return bad1, bad2

def layer20(A, B):
    ga, gb = defaultdict(dict), defaultdict(dict)
    for (a, b), v in A.items():
        ga[a-b][b] = v
    for (a, b), v in B.items():
        gb[a-b][b] = v

    out = defaultdict(int)
    for ell in set(ga) | set(gb):
        ar, br = ga.get(ell, {}), gb.get(ell, {})
        prev, tail = 0, sum(br.values())
        for J in range(max(max(ar, default=-1), max(br, default=-1)) + 1):
            cur = ar.get(J, 0)
            out[ell+2*J] += 2*(cur-prev)*tail
            prev = cur
            tail -= br.get(J, 0)
    return {t: v for t, v in out.items() if v}

def char_table(word):
    P = {(0, 0): 1}
    for z in word:
        n = abs(z)
        eps = 1 if z > 0 else -1
        out = defaultdict(int)
        for (a, b), v in P.items():
            for c in cg(a, n):
                out[c, b] += v
            for c in cg(b, n):
                out[a, c] += eps*v
        P = {k: v for k, v in out.items() if v}
    return P

stamp("exact all-minus witness begins")
L = (-1, -2, -3, -3, -3, -4, -4, -5, -5, -6)
seen, modes, positive = set(), Counter(), Counter()

for i, j in combinations(range(len(L)), 2):
    sig = tuple(sorted((L[i], L[j])))
    if sig in seen:
        continue
    seen.add(sig)
    C = list(L)
    C.remove(L[i])
    C.remove(L[j])
    A, B = cut(C)
    Bp = tuple(sorted(B+sig, key=lambda z: (-abs(z), z)))
    ma = sum(z < 0 for z in A)
    mb = sum(z < 0 for z in Bp)
    assert ma % 2 == mb % 2 == 0 and ma >= 2 and mb >= 2

    X, Y = spvec(A, ma-2), spvec(Bp, mb)
    b1, b2 = string_flags(X, Y)
    modes["20"] += 1
    modes["20_strict"] += not b1 and not b2
    positive["20"] += all(v >= 0 for v in layer20(X, Y).values())

    X2, Y2 = spvec(Bp, mb-2), spvec(A, ma)
    b1, b2 = string_flags(X2, Y2)
    modes["02"] += 1
    modes["02_strict"] += not b1 and not b2
    positive["02"] += all(v >= 0 for v in layer20(X2, Y2).values())

assert len(seen) == 18
assert modes == Counter({"20": 18, "02": 18,
                         "20_strict": 0, "02_strict": 0})
assert positive == Counter({"20": 17, "02": 17})

sig = (-2, -1)
C = list(L)
C.remove(-2)
C.remove(-1)
A, B = cut(C)
Bp = tuple(sorted(B+sig, key=lambda z: (-abs(z), z)))
assert A == (-5, -5, -3, -3)
assert Bp == (-6, -4, -4, -3, -2, -1)

X, Y = spvec(A, 2), spvec(Bp, 6)
rowA = [X.get((2+j, j), 0) for j in range(7)]
rowB = [Y.get((6+j, j), 0) for j in range(8)]
tails = [sum(rowB[J:]) for J in range(8)]
assert rowA == [12, 11, 9, -2, -4, -7, -2]
assert rowB == [197, -135, 17, 17, -33, 31, -11, 1]
assert tails == [84, -113, 22, 5, -12, 21, -10, 1]

ly = layer20(X, Y)
assert sorted(ly.items()) == [
    (0, 360), (2, 1200), (4, 2408), (6, 2712), (8, 1504),
    (10, 562), (12, 242), (14, 136), (16, 30)
]
TA, TB = char_table(A), char_table(Bp)
direct = defaultdict(int)
for (r, s), v in TA.items():
    direct[r+s] += v*TB.get((r, s), 0)
assert {t: v for t, v in direct.items() if v} == ly

stamp("18 pair types: S1+S2 pass 0/18 per allocation; "
      "layer-positive 17/18 per allocation")
stamp(f"pair {sig}: A ell=2 {rowA}; B ell=6 {rowB}; tails {tails}")
stamp(f"pair {sig}: direct height layers {sorted(ly.items())}; PASS")

# Exact (1,1) slice failure for a separate census cut.
L1 = (-1, -2, -3, -3, -3, -4, -4, -5, -5, -5, -5, -20)
sig1 = (-5, -5)
C = list(L1)
C.remove(-5)
C.remove(-5)
A1, B1 = cut(C)
B1p = tuple(sorted(B1+sig1, key=lambda z: (-abs(z), z)))
assert A1 == (-20, -3, -2)
assert B1p == (-5, -5, -5, -5, -4, -4, -3, -3, -1)

X1, Y1 = spvec(A1, 2), spvec(B1p, 8)
dots = defaultdict(int)
for lam, v in X1.items():
    dots[sum(lam)] += v*Y1.get(lam, 0)
assert dots[24] == -34
assert min(X1.values()) == -1 and min(Y1.values()) == -5538
stamp(f"(1,1) pair {sig1}: q=24 slice dot={dots[24]}, "
      f"minimum coefficients=({min(X1.values())}, {min(Y1.values())}); PASS")
