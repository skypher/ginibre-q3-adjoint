
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement as triples
from math import comb
argparse.ArgumentParser(
    description="FM-MECH16 exact endpoint and sector verifier."
).parse_args()

def mul(f,g):
    out=defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            out[i+k,j+l]+=x*y
    return {ij:x for ij,x in out.items() if x}

@lru_cache(None)
def H(k):
    return {(k-2*b-j,j):(-1)**b*comb(k+1-b,b)
            for b in range(k//2+1)
            for j in range(k+1-2*b)}

@lru_cache(None)
def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r,i,j):
    return sum((-1)**b*comb(2*r,b)
               *cat(i+2*r-b)*cat(j+b)
               for b in range(2*r+1))

def direct(r,a,parts):
    f={(a-j,j):comb(a,j) for j in range(a+1)}
    for k in parts:
        f=mul(f,H(k))
    return Fraction(
        sum(x*moment(r,i,j) for (i,j),x in f.items()),2)

@lru_cache(None)
def coeffs(e,a):
    return tuple(
        sum((-1)**i*comb(e,i)*comb(a,k-i)
            for i in range(max(0,k-a),min(e,k)+1))
        for k in range(a+e+1))

def table(cs):
    N=len(cs)-1
    def c(k):
        return cs[k] if 0<=k<=N else 0
    def D(k):
        return c(k)**2-c(k-1)*c(k+1)
    def Wij(i,j):
        return ((c(i-1)+c(i+1))*c(j)
                -c(i)*(c(j-1)+c(j+1)))
    def W(p,q):
        return (0 if (N+p+q)%2 else
                Wij((N+p+q)//2+1,(N+p-q)//2))
    return c,D,W,Wij

def fusion(p,q):
    return range(abs(p-q),p+q+1,2)

def Q4(x,y,z,t):
    return (x+t)*(y+z)-x*t-y*z

def endpoint(cs,parts):
    N=len(cs)-1
    A,B,C=sorted((k+1 for k in parts),reverse=True)
    if (N+A+B+C)%2:
        return 0,0
    c,D,_,_=table(cs)
    al=(N+A-B-C)//2
    be=al+C; ga=al+B; de=al+B+C
    M=(sum(D(k) for k in range(al,be+1))
       -sum(D(k) for k in range(ga+1,de+2)))
    R=(Q4(c(al),c(be),c(ga),c(de+2))
       -Q4(c(al-1),c(be+1),c(ga+1),c(de+1)))
    return M,R

def split(cs,parts):
    A,B,C=sorted((k+1 for k in parts),reverse=True)
    _,_,W,_=table(cs)
    M=sum(W(l,0)
          for d in fusion(A,B) for l in fusion(d,C))
    R=sum(W(d,C) for d in fusion(A,B))
    R+=sum(W(d,B) for d in fusion(A,C))
    R+=sum(W(d,A) for d in fusion(B,C))
    return M,R

def mixed(cs,p,q):
    N=len(cs)-1
    _,D,W,_=table(cs)
    if (N+p+q)%2:
        return 0
    return (D((N+abs(p-q))//2)
            -D((N+p+q)//2+1)-W(p,q))

nd=ns=0
for r in range(2,5):
    for a in range(7):
        cs=coeffs(2*r-3,a)
        for parts in triples(range(1,7),3):
            M,R=endpoint(cs,parts)
            assert M-R==direct(r,a,parts)
            nd+=1
for r in range(2,8):
    for a in range(11):
        cs=coeffs(2*r-3,a)
        for parts in triples(range(1,10),3):
            assert endpoint(cs,parts)==split(cs,parts)
            ns+=1
print("Direct Catalan identities:",nd)
print("Separate main/cross identities:",ns)

np=0
for r in range(2,18):
    e=2*r-3; cs=coeffs(e,e)
    for parts in triples(range(1,20),3):
        if sum(k%2 for k in parts)!=1:
            continue
        C=next(k+1 for k in parts if k%2)
        A,B=[k+1 for k in parts if k%2==0]
        _,_,W,_=table(cs)
        assert sum(W(d,A) for d in fusion(B,C))==0
        assert sum(W(d,B) for d in fusion(A,C))==0
        rows=[mixed(cs,d,C) for d in fusion(A,B)]
        M,R=endpoint(cs,parts)
        assert min(rows)>=0 and sum(rows)==M-R
        np+=1
print("Parity-cancellation identities:",np)

no=0
for r in range(2,21):
    e=2*r-3; cs=coeffs(e,e); c,_,_,_=table(cs)
    for parts in triples(range(1,24,2),3):
        u,v,w=sorted(parts,reverse=True)
        if v%4!=1 or w%4!=1:
            continue
        A,B,C=u+1,v+1,w+1
        al=(2*e+A-B-C)//2
        be=al+C; ga=al+B; de=al+B+C
        if al%2==0:
            x,y,z,t=map(abs,
                (c(al),c(be),c(ga),c(de+2)))
            expected=-(x*y-z*t)-(x+y)*(z-t)
        else:
            x,y,z,t=map(abs,
                (c(al-1),c(be+1),c(ga+1),c(de+1)))
            expected=-(x*y-z*t)-(x-y)*(z+t)
        assert x>=y>=z>=t>=0
        assert endpoint(cs,parts)[1]==expected<=0
        no+=1
print("Ordered-binomial sign identities:",no)

for mode in ("strip","small"):
    count=0
    for r in range(2,18 if mode=="strip" else 11):
        e=2*r-3
        aset=(e-1,e,e+1) if mode=="strip" else range(19)
        labels=range(1,20 if mode=="strip" else 15)
        for a in aset:
            if mode=="small" and min(e,a)>2:
                continue
            N=e+a; cs=coeffs(e,a)
            for parts in triples(labels,3):
                if (a+sum(parts))%2:
                    continue
                A,B,C=sorted(
                    (k+1 for k in parts),reverse=True)
                if A+B-C<=N:
                    continue
                rows=[mixed(cs,d,C) for d in fusion(A,B)]
                M,R=endpoint(cs,parts)
                assert min(rows)>=0 and sum(rows)==M-R
                count+=1
    print("Upper-fusion",mode,"identities:",count)

fake=tuple(
    (-1)**(k//2) if k%2==0 else 0 for k in range(11))
c,D,_,Wij=table(fake)
assert all(c(10-k)==-c(k) for k in range(11))
assert all(D(k)>=D(k+1)>=0 for k in range(5,12))
slacks=[D(j)-D(i)-abs(Wij(i,j))
        for j in range(5,12) for i in range(j+1,13)]
assert min(slacks)>=0
M,R=endpoint(fake,(3,3,3))
assert (M,R,M-R)==(2,3,-1)
print("Artificial E checks:",len(slacks),
      "; minimum slack:",min(slacks))
print("Artificial sequence M,R,M-R:",M,R,M-R)

for r,a,parts in [
    (4,5,(3,4,4)), (4,5,(5,5,5)), (4,5,(3,3,3))
]:
    M,R=endpoint(coeffs(2*r-3,a),parts)
    assert M-R==direct(r,a,parts)
    print((r,a,parts),"M,R,phi:",M,R,M-R)
