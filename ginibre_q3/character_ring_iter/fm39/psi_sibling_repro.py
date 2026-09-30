
from math import comb
from functools import lru_cache

@lru_cache(None)
def row(a, e):
    c = [0] * (a + e + 1)
    for r in range(e + 1):
        for s in range(a + 1):
            c[r + s] += (-1) ** r * comb(e, r) * comb(a, s)
    return c

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def values(a, e, j, i):
    c = row(a, e)
    C = lambda k: at(c, k)
    B = lambda k: C(k - 1) + C(k + 1)
    D = lambda k: C(k)**2 - C(k - 1) * C(k + 1)
    delta = lambda k: D(k) - D(k + 1)

    W = B(i) * C(j) - C(i) * B(j)
    s = i - j
    return {
        "coefficients": c,
        "W": W,
        "delta_j": delta(j),
        "delta_i_minus_1": delta(i - 1),
        "square_left": W * W,
        "square_right": s * s * delta(j) * delta(i - 1),
        "area_slack": sum(delta(k) for k in range(j, i)) - abs(W),
    }

# Verify the sibling-row coefficient identities on a finite exact grid.
sibling_checks = 0
for a in range(13):
    for e in range(13):
        c = row(a, e)
        F = row(a + 2, e)
        G = row(a, e + 2)
        C = lambda k: at(c, k)
        B = lambda k: C(k - 1) + C(k + 1)

        for k in range(-1, a + e + 2):
            fk = at(F, k + 1)
            gk = at(G, k + 1)
            assert fk - gk == 4 * C(k)
            assert fk + gk == 2 * B(k)
            sibling_checks += 2
print("sibling coefficient identities:", sibling_checks, "PASS")

for args in ((3, 5, 7, 9), (3, 5, 5, 8), (0, 12, 6, 8)):
    v = values(*args)
    print(args, v)
    if args == (0, 12, 6, 8):
        assert v["square_left"] - v["square_right"] == 356303376

# Root crossing for N=8, n=N+2=10.
def K(r, X, n):
    if r == 0:
        return 1
    previous, current = 1, X
    for m in range(1, r):
        previous, current = current, X * current - m * (n - m + 1) * previous
    return current

N = 8
X = lambda k: 2 * k - N
assert (X(5), X(8)) == (2, 8)
assert K(5, 4, N + 2) == 0
print("Krawtchouk root crossing:", (X(5), 4, X(8)), "K5(4;10)=0")
