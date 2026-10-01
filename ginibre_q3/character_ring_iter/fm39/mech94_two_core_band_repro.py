import argparse
from math import comb
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from random import Random
import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument("--max-N", type=int, default=30)
ap.add_argument("--max-b", type=int, default=10)
args = ap.parse_args()

def at(c, j):
    return c[j] if 0 <= j < len(c) else 0

@lru_cache(None)
def row(a, e):
    N = a + e
    c = [1]
    for j in range(N):
        z, r = divmod(
            (a-e)*c[-1] - (N-j+1)*(c[-2] if j else 0),
            j+1
        )
        assert r == 0
        c.append(z)
    return tuple(c)

@lru_cache(None)
def poly(a, e, b):
    c = list(row(a, e))
    for _ in range(b):
        c = [at(c,j)+at(c,j-2) for j in range(len(c)+2)]
    return tuple(c)

def Det(c, j):
    return at(c,j)**2 - at(c,j-1)*at(c,j+1)

def K(c, j):
    return at(c,j)**2 - at(c,j-2)*at(c,j+2)

def T(a, e, b, d):
    return sum(
        2**h*comb(b,h)*Det(poly(a,e,b-h),d-h)
        for h in range(b+1)
    )

def C(a, e, b, d):
    c = poly(a,e,b)
    return ((a-e)*at(c,d)-at(c,d-1)-at(c,d+1)
            + (2*b*at(poly(a,e,b-1),d-1) if b else 0))

def S(a, e, b, d):
    return sum(
        2**h*comb(b,h)*K(poly(a,e,b-h),d-h)
        for h in range(b+1)
    )

def W(a, e, b, d, r):
    out = 0
    for h in range(b+1):
        c = poly(a,e,b-h)
        j, i = d-h, r-h
        out += 2**h*comb(b,h)*(
            (at(c,i-1)+at(c,i+1))*at(c,j)
            - at(c,i)*(at(c,j-1)+at(c,j+1))
        )
    return out

# Universal identity used in the quantitative Newton proof.
def elem(xs, j):
    if j < 0 or j > len(xs):
        return sp.Integer(0)
    return sum(
        (sp.prod(v) for v in combinations(xs,j)),
        sp.Integer(0)
    )

identities = 0
for M in range(1,7):
    xs = sp.symbols("x0:"+str(M))
    for j in range(1,min(4,M)+1):
        dj = elem(xs,j)**2-elem(xs,j-1)*elem(xs,j+1)
        rhs = elem(xs,j)**2
        for i in range(M):
            yy = xs[:i]+xs[i+1:]
            rhs += xs[i]**2*(
                elem(yy,j-1)**2-elem(yy,j-2)*elem(yy,j)
            )
        assert sp.expand((j+1)*dj-rhs) == 0
        identities += 1
assert identities == 18
for j in range(1,101):
    assert comb(2*j,j) >= j*(j+1)
print("18 universal square identities and binomial bound PASS",
      flush=True)

count = 0
for N in range(args.max_N+1):
    for e in range(N+1):
        a = N-e
        for b in range(args.max_b+1):
            D = N+2*b
            cp = poly(a,e,b)
            for d in range(D+1):
                assert K(cp,d) >= 2*abs(cp[d])-1
                assert T(a,e,b,d)-1 >= abs(C(a,e,b,d))
                if b and 2 <= d <= D-2:
                    ss = S(a,e,b-1,d-1)
                    aa = at(poly(a,e,b-1),d-1)
                    assert T(a,e,b,d) == (
                        T(a,e,b-1,d)+T(a,e,b-1,d-2)+ss
                    )
                    assert C(a,e,b,d) == (
                        C(a,e,b-1,d)+C(a,e,b-1,d-2)+2*aa
                    )
                    assert ss+1 >= 2*abs(aa)
                count += 1
    print("completed N",N,"coefficient checks",count,flush=True)
if (args.max_N,args.max_b) == (30,10):
    assert count == 169136

# General-deficit transfer.
rng = Random(9421)
for case in range(200):
    a = rng.randrange(9)
    e = rng.randrange(9)
    b = rng.randrange(1,6)
    D = a+e+2*b
    d = rng.randrange(D+1)
    r = rng.randrange(-1,d+1)
    assert W(a,e,b,d,0) == C(a,e,b,d)
    rhs = (
        W(a,e,b-1,d,r)+W(a,e,b-1,d-2,r)
        + W(a,e,b-1,d,r-2)+W(a,e,b-1,d-2,r-2)
        + 2*W(a,e,b-1,d-1,r-1)
    )
    assert W(a,e,b,d,r) == rhs
    rem = (
        S(a,e,b-1,d-1)+T(a,e,b-1,r)
        - T(a,e,b-1,r-2)-S(a,e,b-1,r-1)
    )
    assert T(a,e,b,d)-T(a,e,b,r) == (
        T(a,e,b-1,d)-T(a,e,b-1,r)
        + T(a,e,b-1,d-2)-T(a,e,b-1,r)+rem
    )
print("200 general-deficit transfer identities PASS",flush=True)

# Independent character-ring calculation.
def cgbuild(a, e, b, cap):
    out = {(0,0):1}
    for n,eps,ct in ((1,1,a),(1,-1,e),(2,1,b)):
        for _ in range(ct):
            new = {}
            for (z,j),v in out.items():
                for h in range(n+1):
                    if z+2*h <= cap:
                        key = (z+2*h,j)
                        new[key] = new.get(key,0)+v
                if z+n <= cap:
                    for ell in range(abs(j-n),j+n+1,2):
                        key = (z+n,ell)
                        new[key] = new.get(key,0)+eps*v
            out = {ij:v for ij,v in new.items() if v}
    return out

rng = Random(9420)
bridges = 0
for case in range(120):
    a = rng.randrange(10)
    e = rng.randrange(10)
    b = rng.randrange(7)
    D = a+e+2*b
    if D < 8:
        continue
    d = rng.randrange(4,D//2+1)
    n, p = d-1, D-d-1
    tab = cgbuild(a,e,b,2*d)
    assert T(a,e,b,d) == tab.get((2*d,0),0)
    assert C(a,e,b,d) == (
        tab.get((d+1,d-1),0)-tab.get((d-1,d-1),0)
    )
    for eps in (-1,1):
        new = {}
        for (z,j),v in tab.items():
            for h in range(n+1):
                if z+2*h <= 2*d:
                    key = (z+2*h,j)
                    new[key] = new.get(key,0)+v
            if z+n <= 2*d:
                for ell in range(abs(j-n),j+n+1,2):
                    key = (z+n,ell)
                    new[key] = new.get(key,0)+eps*v
        f = new.get((2*d,0),0)-new.get((2*d-2,0),0)
        assert f == T(a,e,b,d)-1+eps*C(a,e,b,d) >= 0
    bridges += 1
assert bridges == 99
print("99 independent character-ring bridges, both signs PASS",
      flush=True)

@lru_cache(None)
def inv(ns):
    if not ns:
        return 1
    if sum(ns)%2 or 2*max(ns) > sum(ns):
        return 0
    rr = {0:1}
    for n in ns:
        new = {}
        for j,v in rr.items():
            for ell in range(abs(j-n),j+n+1,2):
                new[ell] = new.get(ell,0)+v
        rr = new
    return rr.get(0,0)

def even(word):
    ans = 0
    for mask in range(1 << len(word)):
        left, right = [], []
        sgn = 1
        for i,(n,s) in enumerate(word):
            if mask >> i & 1:
                left.append(n)
                sgn *= s
            else:
                right.append(n)
        ans += (sgn*inv(tuple(sorted(left)))
                * inv(tuple(sorted(right))))
    return ans

for a,e,b,d in (
    (4,2,2,4), (3,3,2,5), (5,1,3,5),
    (4,4,1,4), (2,2,4,6)
):
    D = a+e+2*b
    n,p = d-1,D-d-1
    for eps in (-1,1):
        word = ([(1,1)]*a+[(1,-1)]*e+[(2,1)]*b
                +[(n,eps),(p,eps*(-1)**e)])
        assert even(word) == 2*(
            T(a,e,b,d)-1+eps*C(a,e,b,d)
        )
print("10 direct EVEN split checks PASS",flush=True)

for r in (13,20,31):
    a,e,b,d = 2*r,r,r,2*r
    rr,ww = T(a,e,b,d)-1,C(a,e,b,d)
    print("unbounded family r",r,"R",rr,"C",ww,
          "minimum F",rr-abs(ww),flush=True)
    assert rr >= abs(ww)

# Endpoint normalization is necessary.
q = [Q(comb(4,j),6) for j in range(5)]
assert Det(q,2) == Q(5,9) < 2*abs(q[2])-1

# Separate next-band remainder positivity fails.
for s in range(2,31):
    m,j = 3*s,2*s-1
    lhs = S(m,m,0,j)+T(m,m,0,1)-1
    cp = poly(m,m,0)
    ww = at(cp,j-1)+at(cp,j+1)+2*C(m,m,0,j)
    assert lhs == 3*s-1
    assert abs(ww) == comb(3*s,s)-comb(3*s,s-1) > lhs
    childR = T(m,m,1,2*s)-T(m,m,1,1)
    childW = W(m,m,1,2*s,1)
    assert childR >= abs(childW)
    cs,ct = comb(m,s),comb(m,s-1)
    assert childR == cs*cs+ct*ct-m-1
    assert abs(childW) == m*(cs-ct)
print("29 next-band remainder witnesses; full values positive PASS",
      flush=True)
print("ALL CHECKS PASS",flush=True)