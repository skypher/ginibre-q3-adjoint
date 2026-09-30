
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
from fractions import Fraction

argparse.ArgumentParser(
    description="Bounded exact checks of the three-factor construction."
).parse_args()

def mul(f, g):
    out = defaultdict(int)
    for (i,j), c in f.items():
        for (k,l), d in g.items():
            out[i+k,j+l] += c*d
    return {ij:c for ij,c in out.items() if c}

@lru_cache(None)
def up(k):
    return tuple((k-2*j, (-1)**j*comb(k-j,j))
                 for j in range(k//2+1))

@lru_cache(None)
def h(k):
    out = defaultdict(int)
    for d,c in up(k+1):
        for j in range(d):
            out[d-1-j,j] += c
    return dict(out)

@lru_cache(None)
def hat(k):
    out = defaultdict(int)
    for d,c in up(k):
        out[d,0] += c
        out[0,d] += c
    return dict(out)

def word(ks,a=0):
    f = {(a-j,j):comb(a,j) for j in range(a+1)}
    for k in ks:
        f = mul(f,h(k))
    return f

@lru_cache(None)
def moment(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def kernel(r,i,j):
    return sum((-1)**b*comb(2*r,b)
               *moment(i+2*r-b)*moment(j+b)
               for b in range(2*r+1))

def phi(r,f):
    return Fraction(sum(c*kernel(r,i,j)
                        for (i,j),c in f.items()),2)

def branch(lam,p,q):
    if p < q or q < 0 or (sum(lam)-p-q)%2:
        return 0
    al,be,ga = lam[0]-lam[1],lam[1]-lam[2],lam[2]-lam[3]
    return int(2*be+al+ga >= p+q >= al+ga
               >= p-q >= abs(ga-al))

def B(lam,p):
    return branch(lam,p,0)-branch(lam,p-1,1)+branch(lam,p-2,0)

def epsilon(lam,p):
    al,be,ga = lam[0]-lam[1],lam[1]-lam[2],lam[2]-lam[3]
    return int((al+ga==p and al*ga==0)
               or (al+ga==p-2 and be==0))

def C(lam,p,q):
    return (sum(B(lam,k) for k in range(p-q,p+q+1,2))
            + branch(lam,p,q)+branch(lam,p-2,q)
            - branch(lam,p-1,q+1)-branch(lam,p-1,q-1))

def add_fundamentals(data,n):
    for _ in range(n):
        nxt = defaultdict(int)
        for lam,c in data.items():
            for j in range(4):
                if j == 0 or lam[j] < lam[j-1]:
                    mu = list(lam)
                    mu[j] += 1
                    nxt[tuple(mu)] += c
        data = dict(nxt)
    return data

slacks = 0
for ascending in combinations_with_replacement(range(9),4):
    lam = tuple(reversed(ascending))
    for p in range(2,11):
        assert B(lam,p) == epsilon(lam,p)
        slacks += 1

counts = 0
for r in range(2,7):
    for b,d in combinations_with_replacement(range(1,8,2),2):
        initial = {(b+d-j,j,0,0):1 for j in range(min(b,d)+1)}
        data = add_fundamentals(initial,2*r-3)
        for c in range(2,9,2):
            value = sum(mult*epsilon(lam,c+1)
                        for lam,mult in data.items())
            assert value == phi(r,word((b,c,d))) >= 0
            counts += 1

negative_family = 0
for p in range(5,12,2):
    for q in range(3,p,2):
        assert C((p,q-1,1,0),p,q) == -1
        negative_family += 1

# Jacobi-Trudi: s_521 = h5*h2*h1 - h5*h3 - h6*h1^2 + h7*h1.
s521 = defaultdict(int)
for sign,parts in [(1,(5,2,1)),(-1,(5,3)),
                   (-1,(6,1,1)),(1,(7,1))]:
    for ij,c in word(parts).items():
        s521[ij] += sign*c
direct_negative = phi(1,mul(mul(dict(s521),hat(5)),hat(3)))
assert direct_negative == -1

data = add_fundamentals({(3,0,0,0):1},5)
positive = sum(c*max(0,C(lam,5,3)) for lam,c in data.items())
negative = sum(c*min(0,C(lam,5,3)) for lam,c in data.items())
direct = phi(1,mul(mul(word((3,),5),hat(5)),hat(3)))
assert positive+negative == direct == phi(4,word((2,3,4),1))

print("branching-slack indicator checks:",slacks)
print("positive Schur-count vs Catalan checks:",counts)
print("negative two-hat constituent family checks:",negative_family)
print("direct phi_1(hatS5*hatS3*s_521):",direct_negative)
print("[s_521](h3*h1^5):",data[(5,2,1,0)])
print("phi_4(h2*h3*h4*h1):",positive,"+",negative,"=",direct)
