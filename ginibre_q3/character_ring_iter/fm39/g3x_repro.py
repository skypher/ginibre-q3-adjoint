
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

p = argparse.ArgumentParser(description="FM-MECH22 exact verification")
p.parse_args()

@lru_cache(None)
def row(e, a):
    N = e+a
    if N == 0:
        return (1,)
    c = [1, a-e]
    for k in range(1, N):
        v = (a-e)*c[k]-(N-k+1)*c[k-1]
        assert v % (k+1) == 0
        c.append(v//(k+1))
    return tuple(c)

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def data(e, a, parts):
    cs = row(e, a)
    c = lambda k: at(cs, k)
    A, B, C = sorted((k+1 for k in parts), reverse=True)
    N = e+a
    assert (N+A-B-C) % 2 == 0
    l = (N+A-B-C)//2
    t = l+B+1

    g = lambda k: c(k+1)-c(k-1)
    f = lambda k: c(k)-c(k-1)
    D = lambda k: c(k)**2-c(k-1)*c(k+1)
    AA = lambda x: c(x)-c(x+C)
    BB = lambda x: c(x-1)-c(x+C+1)
    P = lambda x: (
        sum(D(k) for k in range(x, x+C+1))
        -c(x)*c(x+C)+c(x-1)*c(x+C+1)
    )
    mu = lambda x, y: g(x)*g(y)-g(x-1)*g(y+1)
    T = lambda x, L: sum(
        mu(x+2*i-1, x+2*j-1)
        for i in range(1, L//2+1)
        for j in range(i, L//2+1)
    )

    def fterms(x):
        return [
            f(x+s//2)*f(x+(s+1)//2)
            -f(x+max(0, s-C-1))*f(x+min(s, C+1))
            for s in range(2, 2*C+1)
        ]

    for x in (l, t):
        assert P(x) == sum(fterms(x))
        if C % 2 == 0:
            assert P(x) == T(x, C)
        else:
            assert P(x) == (
                T(x, C-1)+D(x+C)
                -c(x)*f(x+C)+c(x-1)*f(x+C+1)
            )

    J = AA(l)*BB(t)-BB(l)*AA(t)
    value = P(l)-P(t)-J
    al, be, ga, de = l, l+C, t-1, t+C-1
    M = sum(D(k) for k in range(al, be+1))
    M -= sum(D(k) for k in range(ga+1, de+2))
    R = (
        c(al)*c(be)-c(al-1)*c(be+1)
        -c(ga+1)*c(de+1)+c(ga)*c(de+2),
        c(al)*c(ga)-c(al-1)*c(ga+1)
        -c(be+1)*c(de+1)+c(be)*c(de+2),
        -c(be)*c(ga)+c(be+1)*c(ga+1)
        +c(al-1)*c(de+1)-c(al)*c(de+2)
    )
    assert value == M-sum(R)

    if C % 2 == 0:
        h = C//2
        rect = sum(
            mu(l+2*i-1, t+2*j)
            for i in range(1, h+1) for j in range(h)
        )
        tri = sum(
            mu(l+2*i-1+s, t+C-s)
            for i in range(1, h+1)
            for s in range(1, C-2*i+2)
        )
        assert J == rect-tri
        assert value == T(l, C)-T(t, C)-rect+tri
    return value, (al, be, ga, de), R, g, fterms

@lru_cache(None)
def hp(k):
    return {
        (k-2*d-j, j): (-1)**d*comb(k+1-d, d)
        for d in range(k//2+1)
        for j in range(k+1-2*d)
    }

def mul(f, g):
    out = defaultdict(int)
    for (i, j), x in f.items():
        for (k, l), y in g.items():
            out[i+k, j+l] += x*y
    return {q: v for q, v in out.items() if v}

@lru_cache(None)
def cat(k):
    return 0 if k % 2 else comb(k, k//2)//(k//2+1)

@lru_cache(None)
def moment(r, i, j):
    return sum(
        (-1)**d*comb(2*r, d)*cat(i+2*r-d)*cat(j+d)
        for d in range(2*r+1)
    )

def direct(r, a, parts):
    f = {(a-j, j): comb(a, j) for j in range(a+1)}
    for k in parts:
        f = mul(f, hp(k))
    return Fraction(
        sum(v*moment(r, i, j) for (i, j), v in f.items()), 2
    )

bridges = far = odd = 0
for r in range(2, 5):
    for a in range(5):
        for parts in combinations_with_replacement(range(1, 6), 3):
            if (sum(parts)+a) % 2:
                continue
            val, ends, R, g, ft = data(2*r-3, a, parts)
            assert val == direct(r, a, parts)
            bridges += 1
            far += ends[2] <= a+2*r-3 and all(R)
            odd += (min(parts)+1) % 2 == 1
print("Catalan bridges:", bridges,
      "all three corrections nonzero:", far, "C odd:", odd)

cross = minor_count = 0
for m in (4, 6, 8):
    for C in (4, 6, 8):
        for gap in (1, 3):
            Y = C+gap
            for extra in (0, 2, 10):
                e = 3*m*Y*Y-3+extra
                a = m
                N = e+a
                B = (N+C+2-gap)//2
                A = B+gap
                val, ends, R, g, ft = data(
                    e, a, (A-1, B-1, C-1)
                )
                al, be, ga, de = ends
                assert ga == N+1 and 0 <= al < be <= N
                assert gap-C < 0 and 3*m*Y*Y <= N-m+4
                folded = lambda k: (-1)**(k+1)*g(k)
                assert folded(al)*folded(be) < 0
                for i in range(1, C//2+1):
                    for j in range(i, C//2+1):
                        p0, q0 = al+2*i-1, al+2*j-1
                        assert (
                            g(p0)*g(q0)-g(p0-1)*g(q0+1) > 0
                        )
                        minor_count += 1
                assert val > 0
                cross += 1
print("Central-root instances:", cross,
      "positive g minors:", minor_count)

odd_cases = balanced_terms = 0
for e in (3, 5, 7, 9):
    for C in (3, 5, 7):
        for gap in (0, 2, 4):
            Y = C+gap
            a = (e+1)*(Y+1)**2//2-2
            N = e+a
            B = (N+C+2-gap)//2
            A = B+gap
            val, ends, R, g, ft = data(
                e, a, (A-1, B-1, C-1)
            )
            al, be, ga, de = ends
            assert ga == N+1 and 0 <= al < be <= N
            assert (e+1)*(Y+1)**2 <= 2*(a+2)
            terms = ft(al)
            assert all(v > 0 for v in terms) and val > 0
            odd_cases += 1
            balanced_terms += len(terms)
print("C-odd central instances:", odd_cases,
      "positive balanced f terms:", balanced_terms)

def kraw(degree, n):
    prev = [1]
    if degree == 0:
        return prev
    cur = [0, 1]
    for j in range(1, degree):
        nxt = [0]+cur
        for k, b in enumerate(prev):
            nxt[k] -= j*(n-j+1)*b
        prev, cur = cur, nxt
    return cur

root_bounds = curvature_bounds = 0
for m in range(2, 22, 2):
    for n in range(2*m+2, 2*m+23):
        K = kraw(m, n)
        H = n-m+2
        reciprocal = Fraction(-K[2], K[0])
        assert 0 < reciprocal <= Fraction(m, 2*H)
        root_bounds += 1
        z = Fraction(H, 3*m)
        assert z <= Fraction(n, 6)
        bnd = Fraction(n)*z/(n*n-z)
        rbnd = Fraction(42, 25)*m*z/H
        assert bnd <= Fraction(6, 35)
        assert rbnd <= Fraction(14, 25)
        assert bnd+rbnd <= Fraction(128, 175) < 1
        curvature_bounds += 1
print("Reciprocal-root bounds:", root_bounds,
      "curvature bounds:", curvature_bounds)

odd_face = cross_face = 0
for r in (4, 6, 8):
    e = 2*r-3
    for C in (3, 5, 7):
        for gap in (0, 2, 4):
            Y = C+gap
            a = (e+1)*(Y+1)**2//2-2
            N = e+a
            B = (N+C-gap)//2
            A = B+gap
            val, ends, R, g, ft = data(
                e, a, (A-1, B-1, C-1)
            )
            al, be, ga, de = ends
            c = row(e, a)
            assert ga == N and at(c, al)-at(c, be) > 0
            assert val == sum(ft(al))+at(c, al)-at(c, be) > 0
            assert all(R)
            odd_face += 1

for m in (4, 6, 8):
    for C in (4, 6, 8):
        for gap in (1, 3):
            Y = C+gap
            for extra in (0, 2, 4, 6):
                e = 3*m*Y*Y-3+extra
                a = m
                N = e+a
                B = (N+C-gap)//2
                A = B+gap
                val, ends, R, g, ft = data(
                    e, a, (A-1, B-1, C-1)
                )
                al, be, ga, de = ends
                if (al+m//2) % 2:
                    continue
                c = row(e, a)
                assert ga == N and at(c, al)-at(c, be) > 0
                assert val > 0 and all(R)
                cross_face += 1
print("C-odd gamma=N instances:", odd_face)
print("Central-root gamma=N instances:", cross_face)

for r, a, parts in (
    (3, 30, (18, 18, 2)),
    (4, 46, (26, 26, 2)),
):
    val = direct(r, a, parts)
    assert val == data(2*r-3, a, parts)[0]
    print((r, a, parts), val)

v = data(297, 4, (153, 152, 3))[0]
print("Root-crossing example positive:", v > 0,
      "integer digits:", len(str(v)))
