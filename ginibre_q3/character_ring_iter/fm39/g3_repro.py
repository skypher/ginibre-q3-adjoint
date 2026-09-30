
import argparse
from math import comb
from functools import lru_cache
from collections import defaultdict
from fractions import Fraction

argparse.ArgumentParser(
    description="FM-MECH21 independent caller checks"
).parse_args()

def row(e,a):
    N=e+a
    if N==0: return (1,)
    out=[1,a-e]
    for k in range(1,N):
        num=(a-e)*out[k]-(N-k+1)*out[k-1]
        assert num%(k+1)==0
        out.append(num//(k+1))
    return tuple(out)

def at(cs,k):
    return cs[k] if 0<=k<len(cs) else 0

def product(f,g):
    out=defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            out[i+k,j+l]+=x*y
    return {p:v for p,v in out.items() if v}

def hp(k):
    return {(k-2*d-j,j):(-1)**d*comb(k+1-d,d)
            for d in range(k//2+1)
            for j in range(k+1-2*d)}

@lru_cache(None)
def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r,i,j):
    return sum((-1)**d*comb(2*r,d)
               *cat(i+2*r-d)*cat(j+d)
               for d in range(2*r+1))

def direct(r,a,parts):
    f={(a-j,j):comb(a,j) for j in range(a+1)}
    for k in parts:
        f=product(f,hp(k))
    return Fraction(
        sum(v*moment(r,i,j) for (i,j),v in f.items()),2
    )

def check(r,a,parts):
    e=2*r-3
    N=e+a
    m=min(e,a)
    A,B,C=sorted((k+1 for k in parts),reverse=True)
    d=A-B-C
    Y=A-B+C
    assert (N+d)%2==0 and C%2==0

    al=(N+d)//2
    be,ga=al+C,al+B
    assert 0<=al<be<=N and ga>=N+1

    central=m%2==1 and (m+1)*Y*Y<=2*(N-m+3)
    outer=d>=0 and d*d>=4*m*(N-m+3)
    assert central or outer

    cs=row(e,a)
    c=lambda k:at(cs,k)
    g=lambda k:c(k+1)-c(k-1)
    D=lambda k:c(k)**2-c(k-1)*c(k+1)

    F=sum(D(k) for k in range(al,be+1))
    F-=c(al)*c(be)-c(al-1)*c(be+1)

    terms=[
        g(al+2*i-1)*g(al+2*j-1)
        -g(al+2*i-2)*g(al+2*j)
        for i in range(1,C//2+1)
        for j in range(i,C//2+1)
    ]

    folded=row(m,N-m)
    gs=[at(folded,k+1)-at(folded,k-1)
        for k in range(al,be+1)]
    assert min(gs)>0 or max(gs)<0
    assert F==sum(terms)==direct(r,a,parts)
    assert all(t>0 for t in terms)
    return F

for args in (
    (3,29,(18,18,3)),
    (20,4,(38,7,3)),
    (3,42,(39,10,3)),
):
    print(args,check(*args))
