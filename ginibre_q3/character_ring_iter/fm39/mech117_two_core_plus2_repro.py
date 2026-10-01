import argparse
from collections import defaultdict
from functools import lru_cache
from math import comb

ap = argparse.ArgumentParser(
    description="FM-MECH117: reanchoring and even-plus insertion.")
ap.add_argument("--max-core", type=int, default=10)
ap.add_argument("--max-a", type=int, default=3)
ap.add_argument("--max-b", type=int, default=2)
args = ap.parse_args()

def clean(F):
    return {q: v for q, v in F.items() if v}

def add(F, G, c=1):
    R = dict(F)
    for q, v in G.items():
        R[q] = R.get(q, 0) + c*v
    return clean(R)

def step(F, n, eps=1):
    R = defaultdict(int)
    for (i, j), v in F.items():
        for r in range(abs(i-n), i+n+1, 2):
            R[r, j] += v
        for r in range(abs(j-n), j+n+1, 2):
            R[i, r] += eps*v
    return clean(R)

@lru_cache(None)
def word(B):
    F = {(0, 0): 1}
    for n in B:
        F = step(F, abs(n), 1 if n > 0 else -1)
    return F

def gain(F, ell=2):
    return add(step(F, ell), F, -1)

def anti(i, j):
    if min(i, j) < 0 or i == j:
        return {}
    return {(i, j): 1, (j, i): -1}

def cone(F):
    return all(i != j and F.get((j, i), 0) == -v
               and (i < j or v >= 0) for (i, j), v in F.items())

def refl(F):
    return {(i, j): (-1)**j*v for (i, j), v in F.items()}

# Algebraic closure identities.
closure = 0
for i in range(1, 19):
    for j in range(i):
        A = anti(i, j)
        rhs = {}
        for r, s in ((i+1, j), (i-1, j), (i, j+1), (i, j-1)):
            rhs = add(rhs, anti(r, s))
        assert step(A, 1) == rhs and cone(rhs)
        closure += 1

        if (i-j) % 2 == 0:
            rhs = {}
            for r, s in ((i+2, j), (i-2, j), (i, j+2)):
                rhs = add(rhs, anti(r, s))
            if j >= 2:
                rhs = add(rhs, anti(i, j-2))
            if j >= 1:
                rhs = add(rhs, A)
            assert gain(A) == rhs and cone(rhs)
            closure += 1

        for ell in (2, 4, 6):
            if i-j >= ell:
                assert cone(gain(A, ell))
                closure += 1

# Independent monomial/Catalan evaluator.
@lru_cache(None)
def U(n):
    return tuple((n-2*j, (-1)**j*comb(n-j, j))
                 for j in range(n//2+1))

def pmul(F, G):
    R = defaultdict(int)
    for (i, j), v in F.items():
        for (k, l), w in G.items():
            R[i+k, j+l] += v*w
    return clean(R)

@lru_cache(None)
def pword(B):
    F = {(0, 0): 1}
    for n in B:
        G = {}
        for i, v in U(abs(n)):
            G[i, 0] = G.get((i, 0), 0) + v
            G[0, i] = G.get((0, i), 0) + (1 if n > 0 else -1)*v
        F = pmul(F, clean(G))
    return F

@lru_cache(None)
def moment(n):
    return 0 if n % 2 else comb(n, n//2)//(n//2+1)

def direct(B, p):
    value = delta = 0
    for (i, j), v in pword(B).items():
        for h, w in U(p):
            x = moment(i+h)
            value += v*w*x*moment(j)
            delta += v*w*(moment(i+h+2)*moment(j)
                          + x*moment(j+2)-3*x*moment(j))
    return value, delta

profiles = bridges = direct_checks = 0
for m in range(3, args.max_core):
    for n in range(m+1, args.max_core+1):
        core = word((-m, n))
        rhs = {q: -v for q, v in anti(n, m).items()}
        for r in range(n-m, n+m+1, 2):
            rhs = add(rhs, anti(r, 0))
        assert core == rhs

        for a in range(args.max_a+1):
            for b in range(args.max_b+1):
                B = (1,)*a + (2,)*b + (-m, n)
                F = word(B)
                Q = gain(F)
                for p in range(n+1, n+6):
                    H = word((1,)*a + (2,)*b + (-p, n))
                    value = F.get((p, 0), 0)
                    delta = Q.get((p, 0), 0)
                    assert value == H.get((m, 0), 0)
                    assert delta == gain(H).get((m, 0), 0)
                    assert value >= 0 and delta >= 0

                    if (p+a+m+n) % 2:
                        assert value == delta == 0
                    elif a or m % 2 == 0:
                        assert a >= (p+n) % 2
                        assert cone(H) and cone(gain(H))
                    elif n % 2:
                        J = word((2,)*b + (-n, m))
                        assert refl(F) == J
                        assert cone(J) and cone(gain(J))
                    else:
                        J = word((2,)*b + (m, n))
                        assert refl(F) == J
                        assert min(J.values()) >= 0
                        assert min(gain(J).values()) >= 0

                    profiles += 1
                    bridges += 2
                    if profiles % 23 == 0:
                        assert direct(B, p) == (value, delta)
                        direct_checks += 1

# Bounded checks of the uniform multiple-plus-core lemma.
dominant = even_insertions = 0
for ns in ((), (3,), (4,), (3, 3), (3, 4), (4, 4), (3, 3, 3)):
    M = sum(ns)
    for off in (0, 1, 4, 5, 6, 7):
        p = M+off
        if p < 3:
            continue
        eta = (p+M) % 2
        for a in (eta, eta+1):
            for b in range(3):
                H = word((1,)*a + (2,)*b + (-p,) + ns)
                assert cone(H) and cone(gain(H))
                dominant += 1
                for ell in (4, 6):
                    if p >= M+ell:
                        assert cone(gain(H, ell))
                        even_insertions += 1

# Cone obstructions are not negative consumer values.
assert gain(anti(2, 0), 4) == add(
    add(anti(6, 0), anti(4, 0)), anti(4, 2), -1)

Q4 = gain(word((1, -8, 6)), 4)
assert Q4.get((3, 2)) == -1
assert gain(word((1, -3, 6)), 4).get((8, 0)) == 5

B = (1, 2, 2, -8, 9, 10)
H = (1, 2, 2, -12, 9, 10)
assert gain(word(B)).get((12, 0)) == 403
assert gain(word(H)).get((8, 0)) == 403
assert gain(word(H)).get((8, 6)) == -82
assert (sum(map(abs, B))-12)//2 == 10

print("cone identities/closures:", closure)
print("two-core profiles:", profiles)
print("reanchoring equalities:", bridges)
print("direct Catalan agreements:", direct_checks)
print("multiple-plus-core certificates:", dominant)
print("even-insertion certificates:", even_insertions)
print("+4 control: cone coefficient -1, consumer gain 5")
print("three-core control: cone coefficient -82, consumer gain 403, distance 10")
print("PASS")