from argparse import ArgumentParser
from itertools import permutations, product
from math import comb

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def lam(p, q, a):
    total = 0
    for eta, target in zip((5, -3, 1), ((2, 1), (3, 2), (4, 1))):
        for perm in permutations((0, 1)):
            ps = 1 if perm == (0, 1) else -1
            for e1, e2 in product((-1, 1), repeat=2):
                u = e1 * target[perm[0]]
                v = e2 * target[perm[1]]
                X, Y = p + 2 - u, q + 1 - v
                h, ell = a + X + Y, a - X + Y
                if h % 2 == 0 and ell % 2 == 0:
                    total += eta * ps * e1 * e2 * C(a, h // 2) * C(a, ell // 2)
    return total

def tail(d, s, a):
    if (a - d) % 2 or d > a:
        return 0
    n, k = (a - d) // 2, (a + d) // 2
    A = -3 * (C(a, n) - C(a, n - 1)) + C(a, n + 1) - C(a, n - 2)
    B = 5 * (C(a, n) - C(a, n - 1)) + C(a, n - 3) - C(a, n + 2)
    D = 5 * (C(a, n - 2) - C(a, n + 1)) - 3 * (C(a, n - 3) - C(a, n + 2))
    return (A * sum(C(a, k - 1 + s + i) for i in range(5))
            + B * sum(C(a, k + s + i) for i in range(3))
            + D * C(a, k + 1 + s))

def main():
    parser = ArgumentParser(description="Exact finite Weyl/tail check for TR.")
    parser.add_argument("--max-a", type=int, default=48)
    limit = parser.parse_args().max_a
    checked = 0
    for a in range(limit + 1):
        top = a + 2
        for d in range(top + 1):
            vals = [lam(q + d, q, a) for q in range(top + 1)]
            pref = [0]
            for z in vals:
                pref.append(pref[-1] + z)
            for t in range(top + 1):
                value = pref[t + 1]
                if (a - d) % 2:
                    assert value == 0, (a, d, t, value)
                elif d <= a:
                    assert value == tail(d, 0, a) - tail(d, t + 1, a)
                elif d == a + 2:
                    assert value == 1, (a, d, t, value)
                else:
                    assert value == 0, (a, d, t, value)
                assert value >= 0, (a, d, t, value)
                checked += 1
        if a % 8 == 0 or a == limit:
            print(f"progress a={a}/{limit} checked_triples={checked}", flush=True)
    print(f"PASS checked_triples={checked} max_a={limit}", flush=True)

if __name__ == "__main__":
    main()
