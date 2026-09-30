
import argparse
from math import comb
from fractions import Fraction
from functools import lru_cache
from collections import defaultdict
from itertools import combinations_with_replacement as triples

argparse.ArgumentParser(
    description="Exact outward Pascal induction and Catalan bridges."
).parse_args()

@lru_cache(None)
def seq(e,a):
    return tuple(
        sum((-1)**i*comb(e,i)*comb(a,k-i)
            for i in range(max(0,k-a),min(e,k)+1))
        for k in range(e+a+1))

def get(cs,k):
    return cs[k] if 0 <= k < len(cs) else 0

def D(c,k):
    return c(k)**2-c(k-1)*c(k+1)

def data(c,N,parts):
    A,B,C = sorted((v+1 for v in parts),reverse=True)
    if (N+A-B-C)%2:
        return None
    al = (N+A-B-C)//2
    be,ga,de = al+C,al+B,al+B+C
    M = sum(D(c,k) for k in range(al,be+1))
    M -= sum(D(c,k) for k in range(ga+1,de+2))
    RC = (c(al)*c(be)-c(al-1)*c(be+1)
          -c(ga+1)*c(de+1)+c(ga)*c(de+2))
    RB = (c(al)*c(ga)-c(al-1)*c(ga+1)
          -c(be+1)*c(de+1)+c(be)*c(de+2))
    RA = (-c(be)*c(ga)+c(be+1)*c(ga+1)
          +c(al-1)*c(de+1)-c(al)*c(de+2))
    return M,(RC,RB,RA),(al,be,ga,de),(A,B,C)

def margin(M,R,e,a,ABC):
    _,B,C = ABC
    terms = []
    if e%2:
        terms.append(sum(R))
    if a%2:
        terms.append((-1)**C*R[0]+(-1)**B*R[1]
                     +(-1)**(B+C)*R[2])
    return M-max(terms)

identities = points = 0
failures = {}
for e in range(17):
    for a in range(e,e+18):
        if e%2 == a%2 == 0:
            continue
        N = e+a
        old,new = seq(e,a),seq(e,a+2)
        c = lambda k: get(old,k)
        q = lambda k: c(k-1)+2*c(k)+c(k+1)
        E = lambda k: c(k)*c(k-1)-c(k-2)*c(k+1)
        H = lambda k: c(k)**2-c(k-2)*c(k+2)
        L = lambda k: (
            D(c,k-1)+D(c,k)+D(c,k+1)
            +2*E(k)+2*E(k+1)+H(k))

        for k in range(-2,N+3):
            assert q(k) == get(new,k+1)
            assert E(k) >= 0 and H(k) >= 0
            assert D(q,k)-D(c,k) == L(k) >= 0
            identities += 1
            if 2*k >= N and L(k) < L(k+1):
                failures.setdefault(
                    "L right monotonicity",
                    (e,a,k,L(k),L(k+1)))

        for parts in triples(range(1,13),3):
            original = data(c,N,parts)
            if original is None:
                continue
            M,R,ends,ABC = original
            M1,R1,_,_ = data(q,N,parts)
            al,be,ga,de = ends
            inc = sum(L(k) for k in range(al,be+1))
            inc -= sum(L(k) for k in range(ga+1,de+2))
            assert M1-M == inc
            g = margin(M,R,e,a,ABC)
            g1 = margin(M1,R1,e,a+2,ABC)
            if g1 < g:
                failures.setdefault(
                    "consumer margin",(e,a,parts,g,g1))
            if inc < 0:
                failures.setdefault(
                    "main-term increment",(e,a,parts,inc))
            points += 1

print("Centered Pascal identities:",identities)
print("Exact consumer-margin steps:",points)
print("First failures:",failures)

# Independent polynomial expansion and Catalan moments.
def mul(f,g):
    out = defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            out[i+k,j+l] += x*y
    return {p:v for p,v in out.items() if v}

@lru_cache(None)
def hp(k):
    return {
        (k-2*d-j,j):(-1)**d*comb(k+1-d,d)
        for d in range(k//2+1)
        for j in range(k+1-2*d)}

@lru_cache(None)
def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r,i,j):
    return sum(
        (-1)**d*comb(2*r,d)*cat(i+2*r-d)*cat(j+d)
        for d in range(2*r+1))

@lru_cache(None)
def direct(r,a,parts):
    f = {(a-j,j):comb(a,j) for j in range(a+1)}
    for k in parts:
        f = mul(f,hp(k))
    return Fraction(
        sum(v*moment(r,i,j) for (i,j),v in f.items()),2)

bridges = steps = 0
for e in range(6):
    for a in range(e,e+7):
        if e%2 == a%2 == 0:
            continue
        for parts in triples(range(1,6),3):
            if (e+a+sum(parts)+3)%2:
                continue
            values = []
            for aa in (a,a+2):
                cs = seq(e,aa)
                c = lambda k: get(cs,k)
                M,R,_,ABC = data(c,e+aa,parts)
                _,B,C = ABC
                exact = []
                if e%2:
                    val = direct((e+3)//2,aa,parts)
                    assert val == M-sum(R)
                    exact.append(val)
                    bridges += 1
                if aa%2:
                    val = direct((aa+3)//2,e,parts)
                    folded = (M-(-1)**C*R[0]
                              -(-1)**B*R[1]
                              -(-1)**(B+C)*R[2])
                    assert val == folded
                    exact.append(val)
                    bridges += 1
                assert min(exact) == margin(M,R,e,aa,ABC)
                values.append(min(exact))
            assert values[1] >= values[0]
            steps += 1

print("Independent orientation bridges:",bridges)
print("Independent consumer-margin steps:",steps)
