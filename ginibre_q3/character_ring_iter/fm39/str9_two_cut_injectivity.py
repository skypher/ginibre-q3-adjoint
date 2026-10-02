import sys
if any(a in ("-h", "--help") for a in sys.argv[1:]):
    print("Usage: python3 -u -B verifier.py; exact FM-STR9 certificates, no files written.")
    raise SystemExit
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from collections import Counter
from math import comb, prod

@lru_cache(None)
def fusion(ns):
    d = {0: 1}
    for n in ns:
        z = {}
        for a, c in d.items():
            for b in range(abs(a-n), a+n+1, 2):
                z[b] = z.get(b, 0) + c
        d = z
    return d

@lru_cache(None)
def table(word):
    d = {(0, 0): 1}
    for v in word:
        n, s = abs(v), (1 if v > 0 else -1)
        z = {}
        for (a, b), c in d.items():
            for t in range(abs(a-n), a+n+1, 2):
                z[t, b] = z.get((t, b), 0) + c
            for t in range(abs(b-n), b+n+1, 2):
                z[a, t] = z.get((a, t), 0) + s*c
        d = {k: c for k, c in z.items() if c}
    return d

def hdim(word):
    a, b = table(word[::2]), table(word[1::2])
    vals = [c*b.get(k, 0) for k, c in a.items()]
    return sum(max(c, 0) for c in vals), sum(max(-c, 0) for c in vals)

@lru_cache(None)
def no_singlet(ns):
    return all(fusion(tuple(n for i, n in enumerate(ns) if s >> i & 1)).get(0, 0) == 0
               for s in range(1, 1 << len(ns)))

def hypotheses(word):
    ns = tuple(map(abs, word))
    minus = [abs(v) for v in word if v < 0]
    assert len(minus) == 2
    n, m = minus
    plus = tuple(v for v in word if v > 0)
    full = (1 << len(plus))-1
    cuts = []
    for s in range(full+1):
        a = tuple(v for i, v in enumerate(plus) if s >> i & 1)
        b = tuple(v for i, v in enumerate(plus) if not s >> i & 1)
        if fusion(a).get(n, 0)*fusion(b).get(m, 0):
            cuts.append(s)
    if not cuts or len(cuts) > 2:
        return False
    if len(cuts) == 2 and cuts[0] ^ cuts[1] != full:
        return False
    return no_singlet(ns[::2]) and no_singlet(ns[1::2])

def det(a):
    a = [list(map(Q, row)) for row in a]
    n = len(a)
    if not n:
        return Q(1)
    ans = Q(1)
    for k in range(n):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return Q(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            ans = -ans
        t = a[k][k]
        ans *= t
        for i in range(k+1, n):
            s = a[i][k]/t
            for j in range(k+1, n):
                a[i][j] -= s*a[k][j]
    return ans

def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]

def transpose(a):
    return list(map(list, zip(*a)))

# Unit-vector cup contractions and the exact four-block determinant.
energy_checks = 0
for q in range(1, 9):
    ra, rb = {}, {}
    for i in range(q+1):
        for j in range(q+1):
            c = Q((-1)**(i+j), q+1)
            ra[i, j, q-i, q-j] = c
            rb[i, j, q-j, q-i] = c
    assert sum(c*c for c in ra.values()) == 1
    assert sum(c*c for c in rb.values()) == 1
    rho = sum(c*rb.get(k, 0) for k, c in ra.items())
    assert rho == Q(1, q+1)
    for ell in range(2, 10):
        k = prod(range(7, ell+7))
        a, b, c, d = 240, 3*k, 2*k, 360
        z = [[1, 1, 0, 0], [0, 0, 1, 1],
             [a, b, 0, 0], [0, 0, c, d]]
        h = [[1, rho, 0, 0], [rho, 1, 0, 0],
             [0, 0, 1, rho], [0, 0, rho, 1]]
        gram = mm(mm(transpose(z), h), z)
        assert det(gram) == (1-rho*rho)**2*((b-a)*(d-c))**2
        bound = (1-rho)*min(Q((b-a)**2, 2+a*a+b*b),
                            Q((d-c)**2, 2+c*c+d*d))
        assert bound > 0
        rem = [[gram[i][j]-(bound if i == j else 0) for j in range(4)]
               for i in range(4)]
        for size in range(1, 5):
            for ix in combinations(range(4), size):
                assert det([[rem[i][j] for j in ix] for i in ix]) >= 0
        energy_checks += 1
print("cup overlap and Gram/energy certificates:", energy_checks, "PASS", flush=True)

# The unbounded-factor family, including its arithmetic separation conditions.
family_checks = 0
for low in ((3, 5, 7), (3, 5, 9), (3, 7, 9), (5, 7, 9)):
    r = fusion(low).get(1, 0)
    assert r in (1, 2)
    for ell in range(2, 11):
        m = sum(low)+3
        high = tuple(m*2**i for i in range(ell-1)) + (m*(2**(ell-1)-1)-1,)
        for s in range(1, (1 << ell)-1):
            z = [v for i, v in enumerate(high) if s >> i & 1]
            assert 2*max(z)-sum(z) >= m-1 > sum(low)+1
        for s in range((1 << 3)-1):
            z = tuple(v for i, v in enumerate(low) if s >> i & 1)
            assert fusion(z).get(1, 0) == 0
        word = (-1, -1)+low+tuple(sorted(high))
        for half in (word[::2], word[1::2]):
            small = tuple(abs(v) for v in half if abs(v) < m-1)
            assert no_singlet(small)
            assert 0 < sum(abs(v) >= m-1 for v in half) < ell
        k = prod(range(7, ell+7))
        assert 3*k != 240 and 2*k != 360
        if ell <= 5:
            assert fusion(tuple(sorted(high))).get(1, 0) == ell-1
        family_checks += 1
print("unbounded family:", family_checks, "parameter instances PASS", flush=True)
for ell in (2, 3, 4):
    high = tuple(18*2**i for i in range(ell-1)) + (18*(2**(ell-1)-1)-1,)
    word = (-1, -1, 3, 5, 7)+tuple(sorted(high))
    he, ho = hdim(word)
    assert ho >= 4*(ell-1)
    assert he-ho == table(word).get((0, 0), 0)
    print("family ell =", ell, "H0 =", (he, ho), flush=True)

# Rational CG realization: STR8b's path ordering and matching.
def clean(p):
    return {v: c for v, c in p.items() if c}

def lower(p, ns):
    z = {}
    for v, c in p.items():
        for i, n in enumerate(ns):
            if v[i] < n:
                w = list(v)
                w[i] += 1
                w = tuple(w)
                z[w] = z.get(w, 0)+c*(n-v[i])
    return clean(z)

def raising(p):
    z = {}
    for v, c in p.items():
        for i, h in enumerate(v):
            if h:
                w = list(v)
                w[i] -= 1
                w = tuple(w)
                z[w] = z.get(w, 0)+c*h
    return clean(z)

def inner(p, q, ns):
    if len(p) > len(q):
        p, q = q, p
    return sum((c*q.get(v, 0)/prod(comb(n, h) for n, h in zip(ns, v))
                for v, c in p.items()), Q(0))

@lru_cache(None)
def cg(ns):
    if not ns:
        return ((0, (), ({(): Q(1)},)),)
    n = ns[-1]
    out = []
    for a, path, states in cg(ns[:-1]):
        for c in range(abs(a-n), a+n+1, 2):
            j = (a+n-c)//2
            hv = {}
            for h in range(j+1):
                for v, z in states[h].items():
                    w = v+(j-h,)
                    hv[w] = hv.get(w, 0)+(-1)**h*comb(j, h)*z
            hv = clean(hv)
            assert not raising(hv)
            st = [hv]
            for h in range(c):
                st.append({v: z/Q(c-h) for v, z in lower(st[-1], ns).items()})
            assert not lower(st[-1], ns)
            out.append((c, path+(c,), tuple(st)))
    assert sum(c+1 for c, _, _ in out) == prod(n+1 for n in ns)
    return tuple(out)

@lru_cache(None)
def copies(word):
    ns = tuple(map(abs, word))
    d = {}
    for mask in range(1 << len(ns)):
        ix = tuple(i for i in range(len(ns)) if mask >> i & 1)
        iy = tuple(i for i in range(len(ns)) if not mask >> i & 1)
        parity = sum(word[i] < 0 for i in iy) % 2
        for a, pa, _ in cg(tuple(ns[i] for i in ix)):
            for b, pb, _ in cg(tuple(ns[i] for i in iy)):
                d.setdefault((a, b), [[], []])[parity].append((mask, pa, pb))
    for z in d.values():
        for v in z:
            v.sort()
    return d

def surplus(word):
    d = {}
    for ab, (ev, od) in copies(word).items():
        r = min(len(ev), len(od))
        if len(ev) > r:
            d[ab] = (0, ev[r:])
        if len(od) > r:
            d[ab] = (1, od[r:])
    return d

def embed(p, ix, L):
    z = {}
    for v, c in p.items():
        w = [0]*L
        for i, h in zip(ix, v):
            w[i] = h
        z[tuple(w)] = c
    return z

def multiply(p, q):
    z = {}
    for a, c in p.items():
        for b, d in q.items():
            v = tuple(x+y for x, y in zip(a, b))
            z[v] = z.get(v, 0)+c*d
    return clean(z)

def states(ns, path):
    return next(st for _, pa, st in cg(ns) if pa == path)

def cup(ns, left, right, pl, pr, spin):
    ls = states(tuple(ns[i] for i in left), pl)
    rs = states(tuple(ns[i] for i in right), pr)
    z = {}
    for h in range(spin+1):
        p = multiply(embed(ls[h], left, len(ns)), embed(rs[spin-h], right, len(ns)))
        for v, c in p.items():
            z[v] = z.get(v, 0)+(-1)**h*comb(spin, h)*c
    return clean(z)

def odd_vectors(word):
    ns = tuple(map(abs, word))
    A, B = tuple(range(0, len(ns), 2)), tuple(range(1, len(ns), 2))
    aa = surplus(tuple(word[i] for i in A))
    bb = surplus(tuple(word[i] for i in B))
    out = []
    for ab in sorted(aa.keys() & bb.keys()):
        ep, al = aa[ab]
        eq, bl = bb[ab]
        if ep ^ eq != 1:
            continue
        for sx, pa, qa in al:
            for tx, pb, qb in bl:
                ax = tuple(i for j, i in enumerate(A) if sx >> j & 1)
                ay = tuple(i for j, i in enumerate(A) if not sx >> j & 1)
                bx = tuple(i for j, i in enumerate(B) if tx >> j & 1)
                by = tuple(i for j, i in enumerate(B) if not tx >> j & 1)
                p = multiply(cup(ns, ax, bx, pa, pb, ab[0]),
                             cup(ns, ay, by, qa, qb, ab[1]))
                assert p and not raising(p) and not lower(p, ns)
                out.append((sum(1 << i for i in ay+by), p))
    return out

def pure_minor(word):
    ns = tuple(map(abs, word))
    A, B = tuple(range(0, len(ns), 2)), tuple(range(1, len(ns), 2))
    od = odd_vectors(word)
    assert len(od) == hdim(word)[1]
    pivots, chosen = {}, []
    full = (1 << len(ns))-1
    for a, pa, _ in cg(tuple(ns[i] for i in A)):
        for b, pb, _ in cg(tuple(ns[i] for i in B)):
            if a != b:
                continue
            p = cup(ns, A, B, pa, pb, a)
            for T in (0, full):
                row = [inner(p, q, ns)*prod(i+2 for i in range(len(ns)) if (S & T) >> i & 1)
                       for S, q in od]
                z = row[:]
                for j, v in sorted(pivots.items()):
                    if z[j]:
                        t = z[j]
                        z = [x-t*y for x, y in zip(z, v)]
                j = next((i for i, x in enumerate(z) if x), None)
                if j is not None:
                    t = z[j]
                    pivots[j] = [x/t for x in z]
                    chosen.append(row)
                if len(chosen) == len(od):
                    answer = det(chosen)
                    assert answer
                    return answer
    raise AssertionError(("missing minor", word, len(chosen), len(od)))

# Every pair-free two-minus list in this stated small box.
eligible = tested = 0
for L in range(2, 8):
    for ns in combinations_with_replacement(range(1, 8), L):
        if sum(ns) > 24:
            continue
        count = Counter(ns)
        negsets = [(v,) for v, c in count.items() if c == 2]
        negsets += list(combinations([v for v, c in count.items() if c == 1], 2))
        for neg in negsets:
            word = tuple(-v if v in neg else v for v in ns)
            if not hypotheses(word):
                continue
            eligible += 1
            he, ho = hdim(word)
            assert he-ho == table(word).get((0, 0), 0) >= 0
            if ho:
                d = pure_minor(word)
                tested += 1
                if tested % 5 == 0:
                    print("exact projected minors:", tested, flush=True)
assert (eligible, tested) == (94, 29)
print("small box: 94 qualifying lists; 29 nonzero harmonic sources; all minors nonzero", flush=True)

# Actual two-minus obstruction to retaining only pure-colour targets.
w = (2, 2, 2, 3, 3, -4, 5, 5, 5, 5, 6, 6, 7, -9)
he, ho = hdim(w)
capacity = 2*fusion(tuple(map(abs, w))).get(0, 0)
assert (he, ho, capacity) == (56676672, 17240634, 15925636)
assert ho > capacity
assert he-ho == 39436038
print("pure-target obstruction:", (he, ho, capacity), "deficiency", ho-capacity, flush=True)
print("FM-STR9 PASS", flush=True)
