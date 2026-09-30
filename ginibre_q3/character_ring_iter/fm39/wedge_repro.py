
import argparse
from math import comb
from fractions import Fraction as Q
from collections import defaultdict
from functools import lru_cache
import sympy as s

argparse.ArgumentParser(
    description="Exact verification of the balanced-family proof."
).parse_args()

def mul(f,g):
    h=defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            h[i+k,j+l]+=x*y
    return {ij:c for ij,c in h.items() if c}

@lru_cache(None)
def H(k):
    h={}
    for b in range((k+1)//2+1):
        d=k+1-2*b
        for j in range(d):
            h[d-1-j,j]=(-1)**b*comb(k+1-b,b)
    return h

def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

def direct(r,m,odd_count):
    f={(2*r-j,j):comb(2*r,j) for j in range(2*r+1)}
    labels=[2*m+1]*odd_count+[2*m]*(3-odd_count)
    for k in labels:
        f=mul(f,H(k))
    return Q(sum(
        v*sum((-1)**b*comb(2*r,b)*cat(i+2*r-b)*cat(j+b)
              for b in range(2*r+1))
        for (i,j),v in f.items()),2)

def cg(p,q):
    return range(abs(p-q),p+q+1,2)

def evaluate(r,m,odd_count=2):
    e=2*r-3
    N=4*r-3
    c=[0]*(N+1)
    for j in range(e+1):
        for k in range(4):
            c[2*j+k]+=(-1)**j*comb(e,j)*comb(3,k)

    def C(k):
        return c[k] if 0<=k<=N else 0

    def W(p,q):
        if p+q>N or (p+q-N)%2:
            return 0
        if p<q:
            return -W(q,p)
        i=(N+p+q)//2+1
        j=(N+p-q)//2
        return ((C(i-1)+C(i+1))*C(j)
                -C(i)*(C(j-1)+C(j+1)))

    K=2*r+m
    if odd_count==2:
        A=2*m+2
        B=2*m+1
        total=sum(W(l,0) for d in cg(A,A) for l in cg(d,B))
        c1=sum(W(d,B) for d in cg(A,A))
        c2=sum(W(d,A) for d in cg(A,B))
        correction=(C(K)**2-C(K-1)**2
                    +2*C(K)*(C(K+1)-C(K-1)))
        assert c1==C(K)**2-C(K-1)**2
        assert c2==C(K)*(C(K+1)-C(K-1))
        assert correction==c1+2*c2
        cap=24
    else:
        assert odd_count==0
        A=2*m+1
        total=sum(W(l,0) for d in cg(A,A) for l in cg(d,A))
        cross=sum(W(d,A) for d in cg(A,A))
        assert cross==C(K)**2-C(K-1)**2
        correction=3*cross
        cap=27

    central=comb(e,r-1)
    T=comb(e,(3*r+1)//2-1)
    tau=Q(4*(2*r-1)*central**2,r*r)
    assert tau==W(1,0)
    assert total>=2*tau
    if T:
        assert Q(central,T)>Q(3*r,4)
    else:
        assert correction<=0
    assert 2*tau>cap*T*T
    assert correction<=cap*T*T
    if K%2:
        assert correction<=0
    return total-correction,2*tau-cap*T*T

identities=bridges=0
for r in range(2,15):
    for m in range(r,2*r+1):
        for odd_count in (0,2):
            value,lower=evaluate(r,m,odd_count)
            assert value>=lower>0
            identities+=1

for r in range(2,9):
    for m in sorted({r,r+1,2*r-1}):
        for odd_count in (0,2):
            value,lower=evaluate(r,m,odd_count)
            assert value==direct(r,m,odd_count)
            bridges+=1

for r in range(4,41):
    e=2*r-3
    b=comb(e,r-1)
    T=comb(e,(3*r+1)//2-1)
    assert Q(b,T)>Q(3*r,4)

a,v=s.symbols("a v")
b=3-a
z=1-3*a
even=(a-3)**2-v**2+2*(a-3)*((3*a-1)-v)
completed=2*b*b+2*b*z-(v-b)**2
assert s.expand(even-completed)==0

A,B,D=s.symbols("A B D")
odd=((3*A-B)**2-(A-3*B)**2
     +2*(3*A-B)*((3*A-D)-(A-3*B)))
assert s.expand(
    odd-(8*(A*A-B*B)+2*(3*A-B)*(2*A+3*B-D))
)==0

print("General split/tail certificate checks:",identities)
print("Independent direct Catalan bridges:",bridges)
print("Central/tail ratio checks:",37)
print("Even completed-square and odd sign identities: passed")
print("Both adjacent-label parity branches passed.")
for n in range(4,9):
    value,lower=evaluate(n,n)
    print("n =",n,"phi =",value,"uniform lower bound =",lower)
