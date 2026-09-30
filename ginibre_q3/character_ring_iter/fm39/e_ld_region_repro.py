
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

def ld_region(L, s, x, y):
    # Exact test for L >= s*(sqrt(x)+sqrt(y)), with nonnegative integers.
    if L < 0:
        return False
    z = L * L - s * s * (x + y)
    return z >= 0 and z * z >= 4 * s**4 * x * y

def data(a, e, j, i):
    c = row(a, e)
    C = lambda k: at(c, k)
    D = lambda k: C(k)**2 - C(k - 1) * C(k + 1)
    B = lambda k: C(k - 1) + C(k + 1)
    S = lambda p, q: C(p) * C(q - 1) - C(p - 1) * C(q)

    s = i - j
    L = D(j) - D(i)
    W = B(i) * C(j) - C(i) * B(j)
    x = D(j) * D(i - 1)
    y = D(j + 1) * D(i)
    return (
        L, W, L - abs(W), s, x, y,
        S(j, i), S(j + 1, i + 1),
        D(j), D(j + 1), D(i - 1), D(i)
    )

all_pairs = ld_pairs = region_pairs = failures = region_failures = 0
for a in range(41):
    for e in range(41):
        N = a + e
        for j in range((N + 1) // 2, N + 1):
            for i in range(j + 1, N + 2):
                L, W, slack, s, x, y, S0, S1, A, B1, C1, Di = data(a, e, j, i)

                assert A >= B1 >= 0
                assert C1 >= Di >= 0 and A >= Di >= 0

                # Theorem LD endpoint bounds, checked exactly after squaring.
                assert S0 * S0 <= s * s * x
                assert S1 * S1 <= s * s * y
                assert W == S0 - S1
                ld_pairs += 2

                all_pairs += 1
                failures += slack < 0
                if ld_region(L, s, x, y):
                    region_pairs += 1
                    region_failures += slack < 0

assert failures == 0 and region_failures == 0
print("consumer pairs a,e<=40:", all_pairs, "screen failures:", failures)
print("Theorem LD endpoint checks:", ld_pairs, "PASS")
print("pairs in the proved LD-energy region:", region_pairs,
      "failures:", region_failures)

# q=0: W_(j+1,j) = D_j - D_(j+1).
q0_checks = 0
for a in range(41):
    for e in range(41):
        N = a + e
        c = row(a, e)
        C = lambda k: at(c, k)
        D = lambda k: C(k)**2 - C(k - 1) * C(k + 1)
        B = lambda k: C(k - 1) + C(k + 1)
        for j in range((N + 1) // 2, N + 1):
            W = B(j + 1) * C(j) - C(j + 1) * B(j)
            assert W == D(j) - D(j + 1)
            q0_checks += 1
print("q=0 identity checks:", q0_checks, "PASS")

for args in ((3, 5, 7, 9), (3, 5, 5, 8)):
    L, W, slack, s, x, y, S0, S1, A, B1, C1, Di = data(*args)
    print("boundary", args, "L,W,slack=", (L, W, slack),
          "LD products=", (x, y), "in LD region=", ld_region(L, s, x, y),
          "chords=", (S0, S1), "D values=", (A, B1, C1, Di))

# Krawtchouk root crossing for N=8, n=N+2=10.
def K(r, X, n):
    if r == 0:
        return 1
    prev, cur = 1, X
    for m in range(1, r):
        prev, cur = cur, X * cur - m * (n - m + 1) * prev
    return cur

N = 8
X = lambda k: 2 * k - N
assert (X(5), X(8)) == (2, 8)
assert K(5, 4, N + 2) == 0
print("root crossing:", (X(5), 4, X(8)), "K5(4;10)=0")
