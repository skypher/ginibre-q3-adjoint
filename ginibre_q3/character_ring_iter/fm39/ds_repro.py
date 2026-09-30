
import argparse
from math import comb
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
from itertools import combinations_with_replacement as triples
argparse.ArgumentParser(
    description="FM-MECH17: diagonal-strip proof and consumer checks."
).parse_args()

def b(e,k):
    return comb(e,k) if 0<=k<=e else 0

@lru_cache(None)
def row(e,P,Q,S):
    al=e+P-Q-S
    be=al+2*S; ga=al+2*Q; de=al+2*Q+2*S
    def c(k):
        return ((-1)**(k//2)*b(e,k//2)
                if k%2==0 and 0<=k<=2*e else 0)
    def D(k):
        return c(k)**2-c(k-1)*c(k+1)
    def q(x,y,z,t):
        return (x+t)*(y+z)-x*t-y*z
    M=(sum(D(k) for k in range(al,be+1))
       -sum(D(k) for k in range(ga+1,de+2)))
    R=(q(c(al),c(be),c(ga),c(de+2))
       -q(c(al-1),c(be+1),c(ga+1),c(de+1)))
    h=al//2
    if al%2==0:
        x,y,z,t=b(e,h),b(e,h+S),b(e,h+Q),b(e,h+Q+S+1)
        upper=x*y+x*z-z*t-y*t-y*z+x*t
    else:
        x,y,z,t=b(e,h),b(e,h+S+1),b(e,h+Q+1),b(e,h+Q+S+1)
        upper=x*y+x*z-z*t-y*t+y*z-x*t
    return M,R,upper,al,h,(x,y,z,t)

total=base=joint=scalar=small=large=0
for e in range(1,28,2):
    for P in range(1,19):
        for Q in range(1,P+1):
            for S in range(1,Q+1):
                M,R,T,al,h,(x,y,z,t)=row(e,P,Q,S)
                assert all(type(v) is int
                           for v in (M,R,T,al,h,x,y,z,t))
                G=M-T
                assert x>=y>=z>=t>=0 and R<=T<=M
                total+=1
                if S==1:
                    if al%2==0:
                        v=b(e,h+Q+1); d=x-y
                        f=(d*d+d*(2*y-z-t)+2*(y*y-z*z)
                           +(z-v)*(2*z+v+t))
                        assert x>=y>=z>=v>=t>=0 and f==G
                    else:
                        v=b(e,h+1)
                        assert x>=v>=y>=z>=t>=0
                        assert v*v>=x*y
                        assert G==(v*v+v*(x+y)-x*y-z*z-t*t
                                   -(x+y)*(z-t))
                        if (P,Q)==(2,1):
                            n=(e-1)//2
                            f=F(8*(n+1)*(n*n+3*n+9),
                                (n+2)**2*(n+3)**2)
                            assert G==b(e,n)**2*f
                        else:
                            assert x+t>=2*z
                            rhs=z*z+t*t+(x+y)*(z-t)
                            assert (x+y)**2*x*y>=rhs**2
                            scalar+=1
                    base+=1
                if Q==P:
                    continue
                xm=b(e,h-1); tp=b(e,h+Q+S+2)
                assert min(x,xm)>=y>=z>=t>=tp>=0
                H=xm-x+t-tp
                if al%2==0:
                    inc=(xm+t)*(xm+x-t-tp)-H*(y+z)
                else:
                    inc=(x+tp)*(xm+x-t-tp)-H*(y+z)
                    assert x+tp>=2*t
                    if x and H>0:
                        u,c,d=F(xm,x),F(t,x),F(tp,x)
                        if u<=F(3,2):
                            low=((1+d)*(1+u-c-d)
                                 -2*min(1,u)*(u-1+c-d))
                            assert F(inc,x*x)>=low>=0
                            small+=1
                        else:
                            assert F(y,x)<=1/u**2
                            assert F(z,x)<=1/u**2
                            assert c<=1/u**3
                            low=u+1-2/u**3-2/u
                            assert F(inc,x*x)>=low>=F(31,54)
                            large+=1
                M1,_,T1,_,_,_=row(e,P,Q+1,S+1)
                assert inc==M1-T1-G and inc>=0
                joint+=1
print("R <= upper <= M:",total)
print("Base identities:",base)
print("Joint-increment identities:",joint)
print("Scalar comparisons:",scalar)
print("Positive-H small/large-ratio branches:",small,large)

def mul(f,g):
    out=defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            out[i+k,j+l]+=x*y
    return {p:v for p,v in out.items() if v}

@lru_cache(None)
def Hpoly(k):
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

@lru_cache(None)
def phi(r,a,parts):
    f={(a-j,j):comb(a,j) for j in range(a+1)}
    for k in parts:
        f=mul(f,Hpoly(k))
    return F(sum(v*moment(r,i,j)
                 for (i,j),v in f.items()),2)

def closed(r,parts):
    e=2*r-3
    if (e+sum(parts))%2:
        return 0
    A,B,C=sorted((k+1 for k in parts),reverse=True)
    al=(2*e+A-B-C)//2
    be=al+C; ga=al+B; de=al+B+C
    def c(k):
        return ((-1)**(k//2)*b(e,k//2)
                if k%2==0 and 0<=k<=2*e else 0)
    def D(k):
        return c(k)**2-c(k-1)*c(k+1)
    def q(x,y,z,t):
        return (x+t)*(y+z)-x*t-y*z
    M=(sum(D(k) for k in range(al,be+1))
       -sum(D(k) for k in range(ga+1,de+2)))
    R=(q(c(al),c(be),c(ga),c(de+2))
       -q(c(al-1),c(be+1),c(ga+1),c(de+1)))
    return M-R

nd=na=0
for r in range(2,7):
    e=2*r-3
    for parts in triples(range(1,8),3):
        val=phi(r,e,parts)
        assert val==closed(r,parts) and val>=0
        nd+=1
    for u,v in triples(range(1,8,2),2):
        for w in range(2,9,2):
            lhs=phi(r,e+1,(u,v,w))
            rhs=(phi(r,e,(u,v,w-1))
                 +phi(r,e,(u,v,w+1)))
            assert lhs==rhs and lhs>=0
            na+=1
print("Direct Catalan-moment bridge:",nd)
print("Adjacent-strip transfers:",na)

for r,parts in [
    (4,(3,3,3)), (5,(3,7,11)),
    (6,(7,7,7)), (7,(7,11,15))
]:
    value=phi(r,2*r-3,parts)
    assert value==closed(r,parts)
    print((r,2*r-3,parts),value)
