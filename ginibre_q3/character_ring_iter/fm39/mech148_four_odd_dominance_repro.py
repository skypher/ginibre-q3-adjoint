import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

argparse.ArgumentParser(
    description="FM-MECH148: four-odd minima, positive descent, and obstruction."
).parse_args()

@lru_cache(None)
def cg(a, b):
    return tuple(range(abs(a-b), a+b+1, 2))

def add(F, G, c=1):
    R = dict(F)
    for ij, v in G.items():
        R[ij] = R.get(ij, 0) + c*v
    return {ij: v for ij, v in R.items() if v}

def step(F, n, eps=1):
    R = defaultdict(int)
    for (i, j), v in F.items():
        for k in cg(i, n):
            R[k, j] += v
        for k in cg(j, n):
            R[i, k] += eps*v
    return {ij: v for ij, v in R.items() if v}

@lru_cache(None)
def word(B):
    F = {(0, 0): 1}
    for z in B:
        F = step(F, abs(z), 1 if z > 0 else -1)
    return F

def phi(B):
    return word(tuple(B)).get((0, 0), 0)

def anti(i, j):
    if min(i, j) < 0 or i == j:
        return {}
    return {(i, j): 1, (j, i): -1}

def cone(F, equal_parity=True):
    return all(i != j and F.get((j, i), 0) == -v
               and (i < j or v >= 0)
               and (not equal_parity or (i-j) % 2 == 0)
               for (i, j), v in F.items())

def dot(F, G):
    return sum(v*G.get(ij, 0) for ij, v in F.items())

def reflect(B):
    return tuple(-z if abs(z) % 2 else z for z in B)

# Independent monomial/Catalan routines.
@lru_cache(None)
def U(n):
    return tuple((n-2*j, (-1)**j*comb(n-j, j))
                 for j in range(n//2+1))

def moment(n):
    return 0 if n % 2 else comb(n, n//2)//(n//2+1)

def pmul(F, G):
    R = defaultdict(int)
    for (i, j), v in F.items():
        for (k, l), w in G.items():
            R[i+k, j+l] += v*w
    return {ij: v for ij, v in R.items() if v}

def pword(B, p=None):
    F = {(0, 0): 1} if p is None else {(i, 0): v for i, v in U(p)}
    for z in B:
        G = defaultdict(int)
        for i, v in U(abs(z)):
            G[i, 0] += v
            G[0, i] += (1 if z > 0 else -1)*v
        F = pmul(F, G)
    return F

def direct(B):
    return sum(v*moment(i)*moment(j)
               for (i, j), v in pword(tuple(B)).items())

# The exact cone identities used in the uniform proof.
basis_checks = 0
for i in range(2, 21):
    for j in range(i):
        if (i-j) % 2:
            continue
        F = anti(i, j)
        rhs = {}
        for u, v in ((i+2, j), (i-2, j), (i, j+2)):
            rhs = add(rhs, anti(u, v))
        if j >= 2:
            rhs = add(rhs, anti(i, j-2))
        if j >= 1:
            rhs = add(rhs, F)
        assert add(step(F, 2), F, -1) == rhs and cone(rhs)
        basis_checks += 1

pair_checks = 0
for a in (1, 3, 5, 7):
    for b in (1, 3, 5, 7):
        for j in range(4):
            for gap in (a+b, a+b+2, a+b+4):
                i = j+gap
                F = anti(i, j)
                R = add(step(step(F, a), b), F, -min(a, b)-1)
                assert cone(R)
                pair_checks += 1

profiles = direct_checks = 0
for odds in combinations_with_replacement((1, 3, 5), 4):
    A, B = odds[:2], odds[2:]
    for off in (0, 2):
        p, q = sum(A)+off, sum(B)+2-off
        F, H = word((-p,)+A), word((-q,)+B)
        F0 = anti(p, 0)
        m = min(A)+1
        R = add(F, F0, -m)
        assert cone(F) and cone(H) and cone(R)
        previous = None
        for t in range(5):
            L = (-p,)+A+(-q,)+B+(2,)*t
            C = (-p, -q)+B+(2,)*t
            value, child = dot(F, H), dot(F0, H)
            assert value == phi(L) == phi(reflect(L))
            assert child == phi(C)
            assert value == m*child+dot(R, H)
            assert value >= m*child >= 0
            if previous is not None:
                assert value >= previous
            previous = value
            profiles += 1
            if profiles % 31 == 0:
                assert direct(L) == value
                direct_checks += 1
            H = step(H, 2)

# Both block parities, and an extra even plus label inside a gap budget.
for p, A, q, B in (
        (5, (3,), 8, (1, 7)),
        (7, (1, 4), 9, (5,)),
        (12, (1, 7, 4), 10, (3, 5)),
        (11, (1, 3, 5), 6, ())):
    assert p >= sum(A) and q >= sum(B)
    assert (p-sum(A)) % 2 == (q-sum(B)) % 2 == 0
    F, H = word((-p,)+A), word((-q,)+B)
    assert cone(F) and cone(H)
    previous = None
    for t in range(5):
        L = (-p,)+A+(-q,)+B+(2,)*t
        assert sum(abs(z) % 2 for z in L) == 4
        value = dot(F, H)
        assert value == phi(L) == phi(reflect(L)) and value >= 0
        if previous is not None:
            assert value >= previous
        previous = value
        profiles += 1
        H = step(H, 2)
print("cone and pair identities:", basis_checks, pair_checks)
print("four-odd sector profiles:", profiles, "Catalan bridges:", direct_checks)

example = []
for t in range(7):
    L = (-8, -8, -1, -3, -5, -7)+(2,)*t
    C = (-8, -8, -3, -5)+(2,)*t
    value, child = phi(L), phi(C)
    fused = sum(phi((-8, -8, j)+(2,)*t) for j in cg(3, 5))
    assert child == fused and value >= 2*child >= 0
    example.append((t, value, child, value-2*child))
print("six-minus example (t,parent,child,remainder):", example)

# Exact Walsh transform and all pair-free even-sign patterns.
@lru_cache(None)
def fusion(ns):
    F = {0: 1}
    for n in ns:
        R = defaultdict(int)
        for i, v in F.items():
            for j in cg(i, n):
                R[j] += v
        F = dict(R)
    return F

def spectrum(ns):
    full = (1 << len(ns))-1
    m = [fusion(tuple(n for i, n in enumerate(ns) if s >> i & 1)).get(0, 0)
         for s in range(full+1)]
    f = [m[s]*m[full ^ s] for s in range(full+1)]
    for i in range(len(ns)):
        for s in range(full+1):
            if not (s >> i & 1):
                t = s | (1 << i)
                u, v = f[s], f[t]
                f[s], f[t] = u+v, u-v
    return f

def masks(ns):
    classes = {n: i for i, n in enumerate(sorted(set(ns)))}
    for bits in range(1 << len(classes)):
        s = sum(1 << i for i, n in enumerate(ns) if bits >> classes[n] & 1)
        if s.bit_count() % 2 == 0:
            yield s

tables = failures = 0
for r in range(6):
    for ev in combinations_with_replacement((2, 4, 6, 8), r):
        ns = (1, 3, 5, 7)+ev
        p = max(ns)
        rest = list(ns)
        rest.remove(p)
        delta = (sum(rest)-p)//2
        if delta < 8 or max(rest) > delta:
            continue
        f = spectrum(ns)
        allowed = list(masks(ns))
        z = min(f[s] for s in allowed)
        minimizers = [s for s in allowed if f[s] == z]
        tables += 1
        if not any(all(not (s >> i & 1) for i in range(4, len(ns)))
                   for s in minimizers):
            failures += 1
assert (tables, failures) == (97, 15)
print("global residual Walsh census:", tables,
      "without an all-even-plus minimum:", failures)

for ev, value, signs in (
        ((4, 6, 6), 342, {3, 12}),
        ((2, 6, 8), 224, {5, 10}),
        ((4, 4, 6), 276, {6, 9}),
        ((4, 6, 8), 350, {23, 24, 48, 63})):
    ns = (1, 3, 5, 7)+ev
    f = spectrum(ns)
    admissible = list(masks(ns))
    minimum = min(f[s] for s in admissible)
    assert minimum == value
    assert {s for s in admissible if f[s] == minimum} == signs
print("global minimum witnesses: 4")

conditional = [
    ((2, -6, -6), (230, 302, 238, 246)),
    ((4, 6, 6), (470, 342, 366, 350)),
    ((2, 6, 8), (304, 236, 224, 244)),
    ((2, -6, -8), (240, 292, 240, 236)),
    ((2, 6, -8), (232, 236, 240, 268)),
    ((2, -6, 8), (248, 228, 256, 244)),
    ((-4, -8, -8), (424, 432, 416, 424)),
    ((-2, 6, 6), (246, 270, 238, 230)),
]
odd = (1, 3, 5, 7)
formula_checks = 0
for C, expected in conditional:
    parity = sum(z < 0 for z in C) % 2
    representatives = (0, 3, 5, 6) if parity == 0 else (1, 2, 4, 7)
    G = word(C)
    h = sum(v*G.get((j, 0), 0) for j, v in fusion(odd).items())
    js = []
    for pair, other in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        js.append(sum(G.get((u, v), 0)
                      for u in cg(odd[pair[0]], odd[pair[1]])
                      for v in cg(odd[other[0]], odd[other[1]])))
    got = []
    for s in range(16):
        if s.bit_count() % 2 != parity:
            continue
        es = tuple(-1 if s >> i & 1 else 1 for i in range(4))
        L = C+tuple(e*n for e, n in zip(es, odd))
        value = 2*(h+es[2]*es[3]*js[0]+es[1]*es[3]*js[1]+es[1]*es[2]*js[2])
        assert value == phi(L) == direct(L)
        if s in representatives:
            got.append((s, value))
        formula_checks += 1
    assert tuple(v for s, v in sorted(got)) == expected
print("conditional four-value identities:", formula_checks)

# Large exact obstruction, evaluated two different ways.
bmax, cutoff = 69, 15
walk = [[1]+[0]*cutoff]
row = [1]
for h in range(bmax):
    new = [0]*(len(row)+1)
    for j, v in enumerate(row):
        new[j+1] += v
        if j:
            new[j] += v
            new[j-1] += v
    row = new
    walk.append((row+[0]*(cutoff+1))[:cutoff+1])

@lru_cache(None)
def background(b, i, j):
    if i % 2 or j % 2:
        return 0
    i, j = i//2, j//2
    assert max(i, j) <= cutoff
    return sum(comb(b, h)*walk[h][i]*walk[b-h][j] for h in range(b+1))

def gp(C, p, b):
    F = {(p, 0): 1}
    for z in C:
        F = step(F, abs(z), 1 if z > 0 else -1)
    return sum(v*background(b, i, j) for (i, j), v in F.items())

@lru_cache(None)
def eta(h, i):
    if i % 2:
        return 0
    return sum((-1)**(h-s)*comb(h, s)*moment(i+2*s) for s in range(h+1))

@lru_cache(None)
def monomial_background(b, i, j):
    return sum(comb(b, h)*eta(h, i)*eta(b-h, j) for h in range(b+1))

def direct_gp(C, p, b):
    return sum(v*monomial_background(b, i, j)
               for (i, j), v in pword(C, p).items())

C, p, b = (-1, -1, 3, 3), 18, 69
cores = (C, (3, 3), (-1, -1), (-1, 3))
values = [gp(T, p, b) for T in cores]
assert values == [
    3985557531774264692874495926119213973103385731444,
    130001226178799539823573056589784101075638988851760,
    34023877698722597580343855134576621372792413564864,
    3995480696683637832196445307409611263264547205760,
]
for T, value in zip(cores, values):
    assert direct_gp(T, p, b) == value
assert 0 < values[0] < min(values[1:])
working = gp(C, p, b-1)
assert working == direct_gp(C, p, b-1)
assert working == 705441163727343646593009249729338359858046924208
assert working < values[0]
print("odd-removal witness parent:", values[0])
print("three child values:", values[1:])
print("working +2 removal:", working)
print("PASS")
