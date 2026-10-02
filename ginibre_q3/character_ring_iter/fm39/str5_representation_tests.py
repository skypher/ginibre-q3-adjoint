from functools import lru_cache
from itertools import combinations_with_replacement
from collections import defaultdict
from fractions import Fraction
import sympy as sp

@lru_cache(None)
def pc(ns):
    d = {0: 1}
    for n in ns:
        q = {}
        for a, v in d.items():
            for c in range(abs(a - n), a + n + 1, 2):
                q[c] = q.get(c, 0) + v
        d = q
    return d.get(0, 0)

@lru_cache(None)
def paths(ns):
    d = {0: [()]}
    for n in ns:
        q = {}
        for a, ps in d.items():
            for c in range(abs(a - n), a + n + 1, 2):
                q.setdefault(c, []).extend(p + (c,) for p in ps)
        d = q
    return tuple(d.get(0, ()))

def dim(ns):
    L = len(ns)
    return sum(
        pc(tuple(ns[i] for i in range(L) if s >> i & 1)) *
        pc(tuple(ns[i] for i in range(L) if not (s >> i & 1)))
        for s in range(1 << L)
    )

def basis(ns):
    L = len(ns)
    out = []
    for s in range(1 << L):
        a = tuple(ns[i] for i in range(L) if s >> i & 1)
        b = tuple(ns[i] for i in range(L) if not (s >> i & 1))
        out.extend((s, p, q) for p in paths(a) for q in paths(b))
    return out

print("CENSUS L, profiles, sum dimensions, maximum:")
for L in range(1, 9):
    rows = [(dim(ns), ns)
            for ns in combinations_with_replacement(range(1, 5), L)]
    assert all(d % 2 == 0 for d, _ in rows)
    print(L, len(rows), sum(x for x, _ in rows), max(rows))
print("BASIS (1,1,1,1):", basis((1, 1, 1, 1)))

def side(ns, s):
    L = len(ns)
    return (
        paths(tuple(ns[i] for i in range(L) if s >> i & 1)),
        paths(tuple(ns[i] for i in range(L) if not (s >> i & 1)))
    )

def canon(s, L):
    return s if not (s >> (L - 1) & 1) else s ^ ((1 << L) - 1)

def slots(ns, T):
    L = len(ns)
    reps = range(1 << (L - 1))
    d = {s: len(side(ns, s)[0]) * len(side(ns, s)[1]) for s in reps}
    pos = [s for s in reps if not ((s & T).bit_count() % 2) and d[s]]
    neg = [s for s in reps if (s & T).bit_count() % 2 and d[s]]
    return pos, neg, d

def rank(M):
    A = [[Fraction(x) for x in row] for row in M]
    if not A:
        return 0
    m, n = len(A), len(A[0])
    i = 0
    for j in range(n):
        q = next((k for k in range(i, m) if A[k][j]), None)
        if q is None:
            continue
        A[i], A[q] = A[q], A[i]
        v = A[i][j]
        A[i] = [x / v for x in A[i]]
        for k in range(m):
            if k != i and A[k][j]:
                v = A[k][j]
                A[k] = [x - v * y for x, y in zip(A[k], A[i])]
        i += 1
        if i == m:
            break
    return i

def channel(ns, T, rule):
    L = len(ns)
    pos, neg, d = slots(ns, T)
    ro, co = {}, {}
    q = 0
    for s in pos:
        ro[s] = q
        q += d[s]
    nr = q
    q = 0
    for s in neg:
        co[s] = q
        q += d[s]
    M = [[0] * q for _ in range(nr)]
    for s in neg:
        es = []
        for i in range(L):
            if not (T >> i & 1):
                continue
            for j in range(L):
                if T >> j & 1 or ns[i] != ns[j]:
                    continue
                if bool(s >> i & 1) == bool(s >> j & 1):
                    continue
                t = canon(s ^ (1 << i) ^ (1 << j), L)
                if t in ro:
                    es.append((i, j, t))
        if rule == "first" and es:
            es = [min(es, key=lambda x: (x[0], x[1]))]
        A, B = side(ns, s)
        src = [(x, y) for x in A for y in B]
        for i, j, t in es:
            C, D = side(ns, t)
            dst = [(x, y) for x in C for y in D]
            flip = (s ^ (1 << i) ^ (1 << j)) != t
            coef = (-1) ** ns[i] if rule == "FS" else 1
            for k, (x, y) in enumerate(src):
                target = (y, x) if flip else (x, y)
                M[ro[t] + dst.index(target)][co[s] + k] += coef
    return M

def map_report(ns, T, show):
    pos, neg, d = slots(ns, T)
    dp = sum(d[s] for s in pos)
    dm = sum(d[s] for s in neg)
    print("MAP", ns, "T=", T, "canonical E+/E-=", dp, dm,
          "full parity=", 2 * dp, 2 * dm,
          "tau eigendims=", dim(ns) // 2, dim(ns) // 2)
    for rule in ("sum", "FS", "first"):
        M = channel(ns, T, rule)
        r = rank(M)
        print(" ", rule, "rank=", r,
              ("matrix=" + str(M) if show
               else "matrix_zero=" + str(not any(map(any, M)))))
        assert r == (1 if show else 0)

map_report((1, 1, 1, 1), 3, True)
map_report((1, 1, 1, 1, 2, 2), (1 << 4) | (1 << 5), False)

z, w = sp.symbols("z w")

def U(n, q):
    return sum(q ** (n - 2 * k) for k in range(n + 1))

def alt(a, b):
    return ((z**a - z**(-a)) * (w**b - w**(-b)) -
            (z**b - z**(-b)) * (w**a - w**(-a)))

def chi(a, b):
    return sp.cancel(alt(a + 2, b + 1) / alt(2, 1))

def h(n):
    return sum(U(a, z) * U(n - 1 - a, w) for a in range(n))

def S(n):
    return U(n, z) + U(n, w)

for n in range(1, 5):
    assert sp.cancel(h(n) - chi(n - 1, 0)) == 0
for n in range(1, 5):
    rhs = chi(1, 0) if n == 1 else chi(n, 0) - chi(n - 1, 1) + chi(n - 2, 0)
    assert sp.cancel(S(n) - rhs) == 0
print("SP4 character formulas exact for h_1..h_4 and S_1..S_4")

A, B = sp.symbols("A B", real=True)
so5 = (sp.sin(A / 2)**2 * sp.sin(B / 2)**2 *
       sp.sin((A + B) / 2)**2 * sp.sin((A - B) / 2)**2)
assert sp.trigsimp(
    (sp.cos(A) - sp.cos(B))**2 -
    4 * sp.sin((A + B) / 2)**2 * sp.sin((A - B) / 2)**2
) == 0
assert sp.trigsimp(
    sp.sin(A / 2)**2 * sp.sin(B / 2)**2 *
    (sp.cos(A) - sp.cos(B))**2 - 4 * so5
) == 0

su3 = (sp.sin((A - B) / 2)**2 *
       sp.sin((2 * A + B) / 2)**2 *
       sp.sin((A + 2 * B) / 2)**2)
g2 = sp.prod(sp.sin(q / 2)**2
             for q in (A, B, A + B, 2 * A + B, 3 * A + B, 3 * A + 2 * B))
print("DENSITY r=0 SO4; r=1 SO5/Sp4 ratio 4")
print("SU3 mismatch witness density:",
      sp.simplify(su3.subs({A: 0, B: sp.pi / 3})))
print("G2 mismatch witness density:",
      sp.simplify(g2.subs({A: sp.pi / 3, B: sp.pi / 3})))
assert sp.trigsimp(
    sp.sin(A / 2)**4 * sp.sin(B / 2)**4 *
    (sp.cos(A) - sp.cos(B))**2 -
    4 * so5 * sp.sin(A / 2)**2 * sp.sin(B / 2)**2
) == 0
print("DENSITY r=2 / SO5 = 4 sin(A/2)^2 sin(B/2)^2")

supp = {(1, 0), (-1, 0), (0, 1), (0, -1)}
su3ref = {(-a, b - a) for a, b in supp}
g2ref = {(-a + 3 * b, b) for a, b in supp}
assert su3ref != supp and g2ref != supp
print("S1 support not Weyl invariant:", sorted(su3ref), sorted(g2ref))

def fourier(n, eps):
    c = defaultdict(int)
    for m in range(-n, n + 1, 2):
        c[(m, m)] += 1
        c[(m, -m)] += eps
    return {k: v for k, v in c.items() if v}

def conv(a, b):
    q = defaultdict(int)
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            q[(i + k, j + l)] += x * y
    return dict(q)

for n in range(1, 9):
    assert min(fourier(n, 1).values()) >= 0
    assert fourier(n, -1)[(n, -n)] == -1
assert conv(fourier(1, -1), fourier(1, -1))[(2, 0)] == -2
print("PD: plus nonnegative; minus coefficient -1 at (n,-n); "
      "W_(1,-)^2 coefficient (2,0)=-2")

f0, f1, f2 = Fraction(0), Fraction(1, 2), Fraction(0)
qform = f0 * f0 + f0 * f2 - 2 * f1 * f1
assert qform == Fraction(-1, 2)
print("Ginibre Q(cos A) =", qform)
