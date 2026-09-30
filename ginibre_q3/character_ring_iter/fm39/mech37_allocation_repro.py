"""FM-MECH37 (astra_max_ceres): channel allocations for the H_AC insertion step; channel-preservation obstruction at (1,1,1;1,2);
allocation polytope; harmonic-sum obstruction; q-boundary certificate."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations

argparse.ArgumentParser(
    description="FM-MECH37: matching-channel obstructions and q-certificate"
).parse_args()

def fusion(mu):
    rows = [{0: 1}]
    for n in mu:
        for row in rows[:]:
            nxt = {}
            for j, a in row.items():
                for k in range(abs(j-n), j+n+1, 2):
                    nxt[k] = nxt.get(k, 0) + a
            rows.append(nxt)
    return rows

def walsh(a):
    a = list(a)
    h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            for j in range(i, i+h):
                a[j], a[j+h] = a[j]+a[j+h], a[j]-a[j+h]
        h *= 2
    return a

def allocation(mu, h, n):
    L = len(mu)
    owner = tuple(
        i for i, k in enumerate(mu+(h, n)) for _ in range(k)
    )

    @lru_cache(None)
    def nc(legs):
        if not legs:
            return ((),)
        out = []
        a = legs[0]
        for z in range(1, len(legs), 2):
            b = legs[z]
            if owner[a] == owner[b]:
                continue
            for p in nc(legs[1:z]):
                for q in nc(legs[z+1:]):
                    out.append(tuple(sorted(((a, b),)+p+q)))
        return tuple(out)

    @lru_cache(None)
    def resolve(M):
        if any(owner[a] == owner[b] for a, b in M):
            return ()
        for i, (a, b) in enumerate(M):
            for j in range(i+1, len(M)):
                c, d = M[j]
                if a < c < b < d:
                    rest = M[:i]+M[i+1:j]+M[j+1:]
                    out = Counter()
                    for edges in (
                        ((a, c), (b, d)),
                        ((a, d), (c, b))
                    ):
                        for k, v in resolve(tuple(sorted(rest+edges))):
                            out[k] += v
                    return tuple(sorted(out.items()))
        caps = sum(
            {owner[a], owner[b]} == {L, L+1} for a, b in M
        )
        return ((h+n-2*caps, 1),)

    b = fusion(mu)
    last = (1 << L)-1
    js = tuple(range(abs(h-n), h+n+1, 2))
    F = {
        j: [
            b[S].get(0, 0)*b[last ^ S].get(j, 0)
            for S in range(last+1)
        ] for j in js
    }
    C = [
        b[S].get(h, 0)*b[last ^ S].get(n, 0)
        for S in range(last+1)
    ]
    CJ = {j: [Q(0)]*(last+1) for j in js}

    for S in range(last+1):
        left = tuple(
            i for i, v in enumerate(owner)
            if v == L or (v < L and S >> v & 1)
        )
        right = tuple(
            i for i, v in enumerate(owner)
            if v == L+1 or (v < L and not (S >> v & 1))
        )
        count = 0
        for A in nc(left):
            for B in nc(right):
                hist = resolve(tuple(sorted(A+B)))
                z = sum(v for j, v in hist)
                assert z > 0
                for j, v in hist:
                    CJ[j][S] += Q(v, z)
                count += 1
        assert count == C[S]
        assert sum(CJ[j][S] for j in js) == C[S]

    margin = min(
        (u-abs(v), j, T)
        for j in js
        for T, (u, v) in enumerate(zip(walsh(F[j]), walsh(CJ[j])))
    )
    return F, C, CJ, margin

cases = [
    ((1, 1, 1), 1, 2, Q(-5, 6)),
    ((1,)*5+(2,), 1, 2, Q(-14677, 798)),
    ((1,)*7, 1, 6, Q(-223, 140)),
    ((2,)*3, 2, 2, Q(-6, 5)),
    ((1,)*6+(2,), 2, 4, Q(-53482133, 3063060))
]
for mu, h, n, expected in cases:
    F, C, CJ, margin = allocation(mu, h, n)
    assert margin[0] == expected
    print(
        mu, ";", h, n, ": objects", sum(C),
        "minimum margin", str(margin[0]),
        "channel", margin[1], "Walsh index", margin[2]
    )

for L in range(2, 10):
    F, C, CJ, margin = allocation((1,)*L, 1, L-1)
    assert F[L] == [1]+[0]*((1 << L)-1)
    assert C == [
        int(S.bit_count() == 1) for S in range(1 << L)
    ]
    assert all(CJ[L][1 << i] == Q(1, L-i) for i in range(L))
    H = sum((Q(1, k) for k in range(1, L+1)), Q(0))
    assert sum(CJ[L]) == H > 1
print("harmonic family L=2..9: 8 exact checks")

F, C, CJ, _ = allocation((1, 1, 1), 1, 2)
forced = [0]*8
forced[4] = 1
low = [c-d for c, d in zip(C, forced)]
assert sum(forced) == sum(F[3]) == 1
assert walsh(F[1])[4]-walsh(low)[4] == -1
print("planar-channel fidelity: forced Fourier margin -1")

# Direct q-weighted matchings for lambda=(1,1,1,1,2).
owner = (0, 1, 2, 3, 4, 4)

@lru_cache(None)
def pairs(legs):
    if not legs:
        return ((),)
    a = legs[0]
    out = []
    for z in range(1, len(legs)):
        b = legs[z]
        if owner[a] == owner[b]:
            continue
        for M in pairs(legs[1:z]+legs[z+1:]):
            out.append(tuple(sorted(((a, b),)+M)))
    return tuple(out)

def moment(legs):
    out = Counter()
    for M in pairs(legs):
        crossings = sum(
            a < c < b < d
            for i, (a, b) in enumerate(M)
            for c, d in M[i+1:]
        )
        out[crossings] += 1
    return out

def square(points, N):
    out = [Q(0)]*N
    for x, a in points.items():
        for y, b in points.items():
            out[x ^ y] += a*b
    return out

sq = square({1: 1, 2: 1, 4: 1, 8: 1}, 16)
checks = 0
for S in range(16):
    left = tuple(
        i for i, v in enumerate(owner) if v < 4 and S >> v & 1
    )
    right = tuple(i for i in range(6) if i not in left)
    direct = Counter()
    for i, a in moment(left).items():
        for j, b in moment(right).items():
            direct[i+j] += a*b

    lowpoly = (
        (2, 3, 1, 0) if S == 0 else
        ((1, 1, 0, 0) if S.bit_count() == 2 else (0, 0, 0, 0))
    )
    highpoly = (1, 2, 2, 1) if S == 0 else (0, 0, 0, 0)
    for d in range(4):
        cert1 = (
            (sq[S]/2 if d in (0, 1) else 0)
            + int(S == 0 and d in (1, 2))
        )
        cert3 = highpoly[d]
        assert cert1 == lowpoly[d]
        assert direct[d] == cert1+cert3
        checks += 2
    assert all(d < 4 for d in direct)
print("q boundary: 16 direct table entries; CP coefficient checks", checks)

# Exact vertices of the allocation polytope.
def indicator(basis):
    W = [0]
    for b in basis:
        W += [x ^ b for x in W]
    assert len(W) == len(set(W))
    return [int(x in W) for x in range(16)]

def halfsum(rows):
    return [
        sum((Q(row[x], 2) for row in rows), Q(0))
        for x in range(16)
    ]

F1 = [2, 0, 0, 1, 0, 1, 1, 0]
F3 = [1, 0, 0, 0, 0, 0, 0, 0]
C = [0, 1, 1, 0, 1, 0, 0, 0]
for i, j in combinations((1, 2, 4), 2):
    k = 7 ^ i ^ j
    D = [Q(0)]*8
    D[i] = D[j] = Q(1, 2)
    G1 = halfsum([
        indicator((i ^ j, i ^ k)),
        indicator((i ^ k, i ^ 8)),
        indicator((j ^ k, j ^ 8)),
        indicator((i ^ j,))
    ])
    G3 = halfsum([
        indicator((i ^ 8,)), indicator((j ^ 8,))
    ])
    assert G1 == F1+[c-d for c, d in zip(C, D)]
    assert G3 == F3+D
print("allocation polytope: 3 nonzero vertices, 96 exact entries")
print("all assertions passed")
