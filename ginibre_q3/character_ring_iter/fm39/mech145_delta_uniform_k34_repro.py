import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial
from random import Random

ap = argparse.ArgumentParser(description="FM-MECH145 exact controls; memory only")
ap.add_argument("--samples", type=int, default=80)
args = ap.parse_args()

def choose(n, j):
    return comb(n, j) if 0 <= j <= n else 0

@lru_cache(None)
def mult(ns, p):
    if p < 0 or (sum(ns)-p) % 2:
        return 0
    if not ns:
        return int(p == 0)
    if len(ns) == 1:
        return int(ns[0] == p)
    d, k = (sum(ns)-p)//2, len(ns)
    return sum((-1)**mask.bit_count() *
               choose(d-sum(ns[i]+1 for i in range(k) if mask>>i&1)+k-2, k-2)
               for mask in range(1<<k))

@lru_cache(None)
def row(a, b):
    t, c = a+2*b, [1]
    for j in range(t):
        v = a*c[j]-(t-j+1)*(c[j-1] if j else 0)
        assert v % (j+1) == 0
        c.append(v//(j+1))
    return tuple(c)

@lru_cache(None)
def base(a, b):
    t, c = a+2*b, row(a, b)
    C = lambda j: c[j] if 0 <= j <= t else 0
    B = lambda j: C(j-1)+C(j+1)
    out = {}
    for h in range(t+1):
        for l in range((t-h)%2, t-h+1, 2):
            i, j = (t+h+l)//2+1, (t+h-l)//2
            v = B(i)*C(j)-C(i)*B(j)
            if v:
                out[h, l] = v
    return out

def edge(a, b):
    return factorial(2*a+2*b+2)*factorial(2*b+2)//(
        2*factorial(a+b+1)*factorial(b+1)*factorial(a+2*b+2))

def direct(word):
    out = {(0, 0): 1}
    for n in word:
        nxt = defaultdict(int)
        for (h, l), v in out.items():
            for r in range(abs(h-abs(n)), h+abs(n)+1, 2):
                nxt[r, l] += v
            for r in range(abs(l-abs(n)), l+abs(n)+1, 2):
                nxt[h, r] += (1 if n > 0 else -1)*v
        out = {i: v for i, v in nxt.items() if v}
    return out

def gp(a, b, cores, p):
    ans, k = 0, len(cores)
    for mask in range(1<<k):
        sel = tuple(abs(cores[i]) for i in range(k) if mask>>i&1)
        rem = tuple(abs(cores[i]) for i in range(k) if not mask>>i&1)
        sign = (-1)**sum(cores[i] < 0 for i in range(k) if mask>>i&1)
        for (h, l), v in base(a, b).items():
            ans += sign*v*mult(sel, l)*mult(rem+(h,), p)
    return ans

def long_ratio(J, t, h):
    v = Q(1)
    for r in range(h):
        v *= Q(J+1-r, t-J+2+r)
    return v*Q(t-J+h+1, t-2*J+2*h)

def newton_xi(d, ns):
    s, k = Q(d+1, d+5), len(ns)
    r, kap = s*s, Q(1)
    for j in range((min(ns)+1)//2):
        kap *= Q((d+1-j)*(d+5), (d+5+j)*(d+1))
    v = Q(1)
    for n in ns:
        v *= 1+r**(n+1)+(n+1)*(1+r)*s**n
    return kap*(v-1)/(1-r)**(k-1)

def gaussian_beta2(d, t, ns):
    J = min(d, t//2)
    h = (min(ns)+1)//2-(d-J)
    if not 1 <= h <= J:
        return None
    A = choose(d+len(ns)-1, len(ns)-1)
    v = 1
    for n in ns:
        v *= 2*n+4
    return long_ratio(J, t, h)*(A*(v-1))**2

for a in range(11):
    for b in range(7):
        t, q = a+2*b, a+2*b+1
        c, L = row(a, b), edge(a, b)
        C = lambda j: c[j] if 0 <= j <= t else 0
        D, S = q*q-a*a, sum(v*v for v in c)
        rho2 = Q(q+1, D)
        lam = Q(2*a, q)
        var = sum((C(j-1)+C(j+1)-lam*C(j))**2 for j in range(-1, t+2))
        assert Q(L, S) == Q(2*D, q*(q+1))
        assert var == Q(4*D, q*q*(q+1))*S
        st = base(a, b)
        assert all(v >= 0 for (h, l), v in st.items() if l == 0)
        assert sum((h+1)*v for (h, l), v in st.items() if l == 0) == L
        for axis in (0, 1):
            for j in range(t+1):
                s = sum(abs(v) for hl, v in st.items() if hl[axis] == j)
                assert s*s <= 4*rho2*L*L
        der = (a*edge(a-1, b) if a else 0)+(4*b*edge(a, b-1) if b else 0)
        assert Q(der, L) == Q((q+1)*(D-q), 2*D)
        P = [C(j)**2-C(j-1)*C(j+1) for j in range(t+1)]
        assert P == P[::-1] and min(P) >= 0
        for J in range(1, t//2+1):
            assert P[J-1] <= P[J]
            for h in range(1, J+1):
                assert P[J-h] <= long_ratio(J, t, h)*P[J]
        if a <= 5 and b <= 3:
            assert st == direct((1,)*a+(-2,)*b)

rng = Random(145)
for k in (3, 4):
    for trial in range(args.samples):
        b = rng.randrange(1, 6)
        a = rng.randrange(4*b+1)
        ns = tuple(sorted(rng.randrange(3, 13) for _ in range(k)))
        signs = {n: rng.choice((-1, 1)) for n in ns}
        cores = tuple(n*signs[n] for n in ns)
        t, q = a+2*b, a+2*b+1
        W = t+sum(ns)
        pmin = max(ns)+(max(ns)-W)%2
        p = rng.randrange(pmin, W+1, 2)
        M = mult(ns, p+t%2)
        F, L = gp(a, b, cores, p), edge(a, b)
        C = 21 if k == 3 else 16*(ns[1]+1)+14*(q+1)
        assert (F-M*L)**2 <= C*C*q*L*L
        if trial < 8:
            assert F == direct((1,)*a+(-2,)*b+cores).get((p, 0), 0)

for k in (3, 4):
    for ns in combinations_with_replacement(range(3, 8), k):
        def ext(p):
            return mult(ns, p) if p >= 0 else -mult(ns, -p-2)
        lip = 1 if k == 3 else ns[0]+1
        for p in range(-sum(ns)-4, sum(ns)+4):
            if (p-sum(ns))%2 == 0:
                assert abs(ext(p+2)-ext(p)) <= 2*lip

def schur(A, l, n):
    return Q((n+1)*choose(A+1, l)*choose(A+1, l+n+1), A+1)

def monomial(d, e, h, degree=0, spin=0):
    A, ans = 2*d+4, Q(0)
    for n in range(h%2, h+1, 2):
        w = choose(h, (h-n)//2)-choose(h, (h-n)//2-1)
        for u in range(abs(n-spin), n+spin+1, 2):
            ans += w*schur(A, d-(e+degree+u)//2, u)
    return ans

for trial in range(400):
    d = rng.randrange(2, 24)
    h = rng.randrange(2*d+1)
    e = h+2*rng.randrange((2*d-h)//2+1)
    degree = rng.randrange(1, 2*d+1)
    spin = rng.randrange(degree%2, degree+1, 2)
    s, kap = Q(d+1, d+5), Q(1)
    for j in range((degree+1)//2):
        kap *= Q((d+1-j)*(d+5), (d+5+j)*(d+1))
    assert monomial(d, e, h, degree, spin) <= (
        (spin+1)*kap*s**degree*monomial(d, e, h))

for d in (512, 1024, 2048):
    for k in (3, 4):
        X = newton_xi(d, ((d+1)//2,)*k)
        assert X < 1
        if d >= 1024 or k == 3:
            assert 128*X < 1
for k, n in ((3, 6), (4, 7)):
    r = Q(1, 4)
    E = ((1+r**(n+1)+(n+1)*(1+r)/2**n)**k-1)/(1-r)**(k-1)
    assert E < Q(7, 8)


def plus_coeff(a, b, j):
    ans = 0
    for v in range(0, a+1, 2):
        for u in range(a-v+1):
            w = j-u-v//2
            for s in range(max(0, w//2+1)):
                r = w-2*s
                if r+s <= b:
                    m = r+v//2
                    ans += (choose(a, v)*choose(a-v, u)*choose(b, s)*
                            choose(b-s, r)*choose(2*m, m)//(m+1))
    return ans

for d in range(1, 13):
    for a in (0, 3, 8):
        Pp = [plus_coeff(a, 2*d+4, j) for j in range(d+1)]
        assert all(4*Pp[j-1] <= Pp[j] for j in range(1, d+1))

d, A = 2048, 4100
P = [1]
for j in range(1, d+1):
    top = P[-1]*(A+2-j)*(A+1-j)
    assert top % (j*(j+1)) == 0
    P.append(top//(j*(j+1)))
for k in (3, 4):
    F = sum(choose(d-j+k-2, k-2)*P[j] for j in range(d+1))
    F += -k*choose(A, d)+choose(k, 2)
    assert 128*F >= 127*P[-1]


c = row(0, d)
Pm = [c[j]**2-(c[j-1] if j else 0)*c[j+1] for j in range(d+1)]
for k in (3, 4):
    assert gaussian_beta2(d, 2*d, (d,)*k) < Q(1, 128**2)
    F = sum(choose(d-j+k-2, k-2)*Pm[j] for j in range(d+1))
    F += -k*c[d]+choose(k, 2)
    assert 128*F >= 127*Pm[-1]

for ns, p in (((165, 166, 167), 168), ((250, 251, 252, 253), 254)):
    a, b, q, k = 40, 10, 61, len(ns)
    M = mult(ns, p)
    C = 21 if k == 3 else 16*(ns[1]+1)+14*(q+1)
    assert M*M >= C*C*q
    F = gp(a, b, tuple(-n for n in ns), p)
    assert F > 0
    print("exact consumer:", (a, b, ns, p), "M =", M, "g =", F)

for k in (3, 4):
    d = 8192
    assert gaussian_beta2(d, 2*d, (d//2,)*k) < Q(1, 128**2)
d = 65536
assert (d//8)**2 >= 21**2*(2*d)
assert 23328*2**(15*16) < 2**(508-14)
assert 2*choose(14, 7)//8-choose(16, 8)//9 == -572
assert 2204**2-2203*2204-881 == 1323
print("PASS: identities, both uniform certificates, interior bounds, consumers")
