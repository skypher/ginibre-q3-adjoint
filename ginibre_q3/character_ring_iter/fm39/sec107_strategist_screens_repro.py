"""FM-SEC107 (luna_max_vesta): character-defect matrices (U-FM3 with one extra pair), Lee-Yang root condition and split-count descent screens."""
from functools import lru_cache
from itertools import combinations_with_replacement
from sympy import Poly, symbols

@lru_cache(None)
def fusion_count(ns):
    state = {0: 1}
    for n in ns:
        nxt = {}
        for a, v in state.items():
            for c in range(abs(a - n), a + n + 1, 2):
                nxt[c] = nxt.get(c, 0) + v
        state = nxt
    return state.get(0, 0)

def profile(ns, T):
    L = len(ns)
    f = []
    for S in range(1 << L):
        left = tuple(ns[i] for i in range(L) if (S >> i) & 1)
        right = tuple(ns[i] for i in range(L) if not ((S >> i) & 1))
        f.append(fusion_count(left) * fusion_count(right))
    t = T.bit_count()
    a = [
        sum(f[S] for S in range(1 << L)
            if ((S & T).bit_count() == j))
        for j in range(t + 1)
    ]
    value = sum((-1) ** j * a[j] for j in range(t + 1))
    return a, value

def circle_test(a):
    if not any(a):
        return True, None, 0, 0
    if a[0] == 0:
        return False, "zero at origin", 0, 0

    m = (len(a) - 1) // 2
    u = symbols("u")
    if m == 0:
        return True, Poly(a[0], u), 1, 1

    c0 = Poly(2, u)
    c1 = Poly(u, u)
    q = Poly(a[m], u) + a[m - 1] * c1
    prev, cur = c0, c1
    for k in range(2, m + 1):
        nxt = Poly(u, u) * cur - prev
        q += a[m - k] * nxt
        prev, cur = cur, nxt

    roots_in_range = q.count_roots(-2, 2)
    return roots_in_range == q.degree(), q, roots_in_range, q.degree()

profiles = [
    ns
    for L in range(1, 8)
    for ns in combinations_with_replacement((1, 2, 3), L)
]

total = radial_bad = circle_bad = value_negative = 0
first_radial = first_circle = None

for ns in profiles:
    for T in range(1 << len(ns)):
        if T.bit_count() % 2:
            continue
        a, value = profile(ns, T)
        total += 1

        if any(a[j] < a[j + 1] for j in range(len(a) // 2)):
            radial_bad += 1
            if first_radial is None:
                first_radial = (ns, T, a, value)

        if value < 0:
            value_negative += 1

        ok, q, roots, degree = circle_test(a)
        if not ok:
            circle_bad += 1
            if first_circle is None:
                qexpr = q.as_expr() if isinstance(q, Poly) else q
                first_circle = (ns, T, a, value, qexpr, roots, degree)

print("profile cases", total)
print("descent failures", radial_bad, "first", first_radial)
print("negative EVEN values", value_negative)
print("circle failures", circle_bad, "first", first_circle)

for ns in [(1, 5, 2, 2), (1,) * 8 + (6,)]:
    count = descent_pass = circle_pass = 0
    values = []
    first_descent = first_circle = None
    for T in range(1 << len(ns)):
        if T.bit_count() % 2:
            continue
        a, value = profile(ns, T)
        count += 1
        values.append(value)

        if not any(a[j] < a[j + 1] for j in range(len(a) // 2)):
            descent_pass += 1
        elif first_descent is None:
            first_descent = (T, a, value)

        ok, q, roots, degree = circle_test(a)
        if ok:
            circle_pass += 1
        elif first_circle is None:
            qexpr = q.as_expr() if isinstance(q, Poly) else q
            first_circle = (T, a, value, qexpr, roots, degree)

    print("tight", ns, "T cases", count,
          "descent pass", descent_pass,
          "circle pass", circle_pass,
          "minimum EVEN value", min(values),
          "first descent failure", first_descent,
          "first circle failure", first_circle)

@lru_cache(None)
def cg(a, b):
    return tuple(range(abs(a - b), a + b + 1, 2))

def multiply_factor(R, n, sign):
    out = {}
    for (a, b), v in R.items():
        for c in cg(a, n):
            out[c, b] = out.get((c, b), 0) + v
        for c in cg(b, n):
            out[a, c] = out.get((a, c), 0) + sign * v
    return out

def multiply_B(R, n):
    out = {}
    for (a, b), v in R.items():
        for i in range(n):
            j = n - 1 - i
            for c in cg(a, i):
                for d in cg(b, j):
                    out[c, d] = out.get((c, d), 0) + v
    return out

def multiply_S(R, n):
    out = {}
    for (a, b), v in R.items():
        for c in cg(a, n):
            out[c, b] = out.get((c, b), 0) + v
        for d in cg(b, n):
            out[a, d] = out.get((a, d), 0) + v
    return out

def defect_check(Q):
    support = [(a, b, v) for (a, b), v in Q.items() if v]
    if not support:
        return (0, None), None

    D = max(max(a, b) for a, b, _ in support)
    rho = [Q.get((c, 0), 0) for c in range(D + 1)]

    # Prefix sums for each parity, with indices kept in character labels.
    pref = [[0] * (D + 2) for _ in range(2)]
    for parity in range(2):
        for c in range(D + 1):
            pref[parity][c + 1] = (
                pref[parity][c] +
                (rho[c] if c % 2 == parity else 0)
            )

    def interval_sum(d, u):
        if d > u:
            return 0
        parity = d % 2
        return pref[parity][u + 1] - pref[parity][d]

    best = None
    first_bad = None

    # Finite interior.
    for a in range(1, D + 1):
        for b in range(1, D + 1):
            gap = (
                interval_sum(abs(a - b), min(a + b, D))
                - Q.get((a, b), 0)
            )
            item = (gap, a, b)
            if best is None or gap < best[0]:
                best = item
            if gap < 0 and first_bad is None:
                first_bad = item

    # If max(a,b)>D, the gap is a tail sum indexed by |a-b|.
    for d in range(D + 1):
        gap = interval_sum(d, D)
        item = (gap, "tail", d)
        if best is None or gap < best[0]:
            best = item
        if gap < 0 and first_bad is None:
            first_bad = item

    return best, first_bad

def factored_core(ns, T):
    F = {(0, 0): 1}
    t = T.bit_count()
    for i, n in enumerate(ns):
        F = multiply_B(F, n) if ((T >> i) & 1) else multiply_S(F, n)
    for _ in range(t - 2):
        F = multiply_factor(F, 1, -1)
    return F

checks = failures = 0
first = None
minimum = None

for ns in profiles:
    N = sum(ns)
    for T in range(1 << len(ns)):
        t = T.bit_count()
        if t % 2 or t == 0:
            continue
        core = factored_core(ns, T)
        for k in range(N + 1):
            best, bad = defect_check(multiply_S(core, k))
            checks += 1
            if bad is not None:
                failures += 1
                if first is None:
                    first = (ns, T, k, bad)
            if minimum is None or best[0] < minimum[0]:
                minimum = (best[0], ns, T, k, best[1:])

print("defect matrices", checks, "with negative entries", failures)
print("first", first, "minimum gap", minimum)

for ns in [(1, 5, 2, 2), (1,) * 8 + (6,)]:
    checks = failures = 0
    first = None
    minimum = None
    N = sum(ns)
    for T in range(1 << len(ns)):
        t = T.bit_count()
        if t % 2 or t == 0:
            continue
        core = factored_core(ns, T)
        for k in range(N + 1):
            best, bad = defect_check(multiply_S(core, k))
            checks += 1
            if bad is not None:
                failures += 1
                if first is None:
                    first = (T, k, bad)
            if minimum is None or best[0] < minimum[0]:
                minimum = (best[0], T, k, best[1:])
    print("tight", ns, "defect matrices", checks,
          "with negative entries", failures,
          "first", first, "minimum gap", minimum)
