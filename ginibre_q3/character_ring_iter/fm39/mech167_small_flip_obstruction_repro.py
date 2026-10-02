import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from math import comb
from fractions import Fraction as Q

ap = argparse.ArgumentParser(description="FM-MECH167 exact, memory-only verifier.")
ap.parse_args()
def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def polynomial(word):
    f = {(0, 0): 1}
    for z in word:
        n, e = abs(z), 1 if z > 0 else -1
        g = defaultdict(int)
        for (i, j), v in f.items():
            for k in cg(i, n): g[k, j] += v
            for k in cg(j, n): g[i, k] += e*v
        f = {ij: v for ij, v in g.items() if v}
    return f

def check(m, A, fixed, p, q, cutoff, expected_W):
    powers = [[1]]
    for t in range(A):
        row = [0]*(len(powers[-1])+m)
        for i, v in enumerate(powers[-1]):
            if v:
                for j in cg(i, m): row[j] += v
        powers.append(row)

    @lru_cache(None)
    def moment(t, i, j):
        if (i | j) & 1: return 0
        return sum(comb(t, u) *
                   (powers[u][i] if i < len(powers[u]) else 0) *
                   (powers[t-u][j] if j < len(powers[t-u]) else 0)
                   for u in range(t+1))

    @lru_cache(None)
    def entry(t, word, a=0, b=0):
        return sum(v*sum(moment(t, x, y)
                         for x in cg(i, a) for y in cg(j, b))
                   for (i, j), v in polynomial(word).items())

    def ev(t, word, a=0, b=0):
        return entry(t, tuple(sorted(word)), a, b)

    def remove(t, word, pair):
        word = list(word)
        for z in pair:
            if z == m: t -= 1
            else: word.remove(z)
        assert t >= 0
        return t, word

    full = [m]*A + list(fixed)
    B = full[:]; B.remove(p)
    W = sum(map(abs, B))
    assert W == expected_W and (W-p) % 2 == 0
    assert p >= max(map(abs, B)) and p >= 6
    assert (W-p)//2 >= max(8, max(map(abs, B)))
    assert all(-z not in full for z in full)
    assert sum(z < 0 for z in full) % 2 == 0
    assert max((abs(a)+abs(b), max(abs(a), abs(b)))
               for i, a in enumerate(B) for b in B[i+1:]
               if (a-b) % 2 == 0) == (2*q, q)
    phi = ev(A, fixed)
    t, child = remove(A, fixed, (-q, -q))
    phi_child = ev(t, child)
    assert phi > 0 and phi % 2 == phi_child % 2 == 0
    assert 3*phi < 2*phi_child and phi_child < 2*phi

    counts = Counter(full)
    classes = sorted(counts, key=abs)
    ds = {}
    for i, a in enumerate(classes):
        for b in classes[i:]:
            if a == b and counts[a] < 2: continue
            t, rest = remove(A, fixed, (a, b))
            flipped = rest + [-a, -b]
            difference = phi-ev(t, flipped)
            assert difference % 4 == 0
            d = difference//4
            assert d == (1 if b > 0 else -1)*ev(t, rest, abs(a), abs(b))
            ds[a, b] = d
    small = {ab: d for ab, d in ds.items()
             if min(map(abs, ab)) <= cutoff}
    assert small and max(small.values()) < 0
    assert ds[-q, -q] < 0 and ds[-q, p] < 0
    assert 0 < 4*ds[m, m] < phi

    r = full.count(1)
    if r == 2:
        rest = list(fixed); rest.remove(1); rest.remove(1)
        rhs = 0
        rc = Counter(rest); rc[m] += A
        for z, mult in rc.items():
            t, sub = remove(A, rest, (z,))
            rhs += 4*mult*(abs(z)+1)*(1 if z > 0 else -1)*ev(t, sub, 1, abs(z)-1)
        assert (W+p+4)*ds[1, 1] == rhs
    else:
        assert r == 1
        rest = list(fixed); rest.remove(1)
        rhs = 0
        rc = Counter(rest); rc[m] += A
        for z, mult in rc.items():
            n, e = abs(z), 1 if z > 0 else -1
            t, sub = remove(A, rest, (z,))
            lower = ev(t, sub+[e*(n-1)])
            rhs += 2*mult*((n+1)*lower+n*ds[1, z])
        assert (W+p+2)*phi == rhs

    # Independently compute the both-in-x term of the triple-flip identity.
    child.remove(p)
    T1 = sum(ev(A, child, s, 0)
             for c in cg(q, q) for s in cg(p, c))
    assert phi == 2*T1 + ds[-q, -q] + 2*ds[-q, p]
    assert T1 >= phi_child//2
    K = Q((m*(m+2))**2,
          4*(q+1)**2*(q*(q+2))**2*((4 if m == 2 else 6)**2-1))
    assert K == (Q(1, 13500) if m == 2 else Q(1, 27440))
    print("W =", W, "delta =", (W-p)//2, "L =", len(full),
          "small flips =", len(small), "all negative;",
          "3/2 < child/parent < 2; carrier flip positive; PASS")

check(2, 128, (1,-3,-3,-3,-3,-4,-4,7), 7, 4, 1, 277)
check(2, 128, (1,1,-3,-3,-3,-3,-4,-4,6), 6, 4, 1, 278)
check(4, 192, (1,-2,-2,-3,-3,-3,-3,-6,-6,9), 9, 6, 2, 797)
check(4, 192, (1,1,-2,-2,-3,-3,-3,-3,-6,-6,8), 8, 6, 2, 798)
print("ALL EXACT CHECKS PASS")
