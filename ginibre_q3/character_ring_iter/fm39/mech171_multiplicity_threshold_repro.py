import argparse, re
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb, prod
from pathlib import Path

ap=argparse.ArgumentParser(description="FM-MECH171 exact verifier; memory only.")
ap.add_argument("--no-census",action="store_true")
args=ap.parse_args()
def cg(a,b): return range(abs(a-b),a+b+1,2)
def ceil2(n): return (n-1).bit_length()
def U(n,x):
    a,b=F(1),x
    if n==0:return a
    for _ in range(1,n):a,b=b,x*b-a
    return b
for n in range(1,17):
    for j in range(-32,33):
        x=F(j,16)
        assert abs(U(n,x))/(n+1)<=1-(2-abs(x))/3

constant_checks=0
for q in range(2,25):
    ell=ceil2(q)
    for n in range(1,q):
        h=n+1
        for t in range(9):
            for minus in (False,True):
                if minus and n%2:continue
                eta=F(1,q*q*(h*h if minus else 1))
                den=8*q*q*(h**4 if minus else h*h)
                K=(16+3*ceil2(den)+t*ceil2(2*h*h) if minus
                   else 16+4*t+(t+3)*ceil2(den))
                r0=2+128*q*q*(h*h if minus else 1)*(t+3)*ell
                block=8/eta
                assert block.denominator==1
                ratio=(1-eta/24)/(1-eta/6)
                assert block*(ratio-1)>=1
                assert r0-2>=block*K
                constant_checks+=1
print("character-gap and threshold controls:",constant_checks,"PASS")

@lru_cache(None)
def poly(C):
    out={(0,0):1}
    for z in C:
        nxt=defaultdict(int);n=abs(z);e=1 if z>0 else -1
        for (i,j),v in out.items():
            for k in cg(i,n):nxt[k,j]+=v
            for k in cg(j,n):nxt[i,k]+=e*v
        out={k:v for k,v in nxt.items() if v}
    return out
@lru_cache(None)
def choose(n):return tuple(comb(n,k) for k in range(n+1))
def mu1(t,i):
    if i>t or (t-i)%2:return 0
    k=(t-i)//2
    return comb(t,k)-(comb(t,k-1) if k else 0)
@lru_cache(None)
def moment1(t,i,j):
    return sum(v*mu1(k,i)*mu1(t-k,j) for k,v in enumerate(choose(t)))
def entry(moment,t,C,a=0,b=0):
    return sum(v*sum(moment(t,x,y) for x in cg(i,a) for y in cg(j,b))
               for (i,j),v in poly(tuple(sorted(C))).items())
def J(m,h):
    v=F(comb(2*m,m)*comb(2*h,h),comb(m+h,m))
    z=F(2*(2*m+1)*(2*h+1),(m+h+1)**2*(m+h+2))*v*v
    assert z.denominator==1
    return z.numerator
for C,r in (((-2,-2),24),((2,2),34),((-3,-3),52),
            ((-2,3,-5),128),((-3,4,-7,8),496)):
    h=sum(z<0 for z in C)//2;b=2*h+3
    K=sum(5*z*(z+2) if z>0 else 3*(-z-1)*(-z+3) for z in C)
    score=15*(r*r-(b-1)**2)-K*b*(2*r+b)
    threshold=((2*K+15)*b-15+14)//15
    assert r>=threshold and score>=0
    c=prod(2*(z+1) if z>0 else (-z)*(-z+1)*(-z+2)//6 for z in C)
    d=entry(moment1,r-2,C,1,1)
    lower=F(4*c*J((r-2)//2,h)*score,15*((r+b)**2-1))
    assert d>=lower>=0 and entry(moment1,r,C)>0
    rhs=4*(r-2)*entry(moment1,r-2,C)
    for z,mult in Counter(C).items():
        R=list(C);R.remove(z)
        rhs+=4*mult*(abs(z)+1)*(1 if z>0 else -1)*entry(
            moment1,r-2,R,1,abs(z)-1)
    assert (r+sum(map(abs,C))+4)*d==rhs
    print("fundamental threshold:",C,r,"PASS")

# Multiplicity-five obstruction to all flips incident to that class.
A=512
C=(1,)*5+(-3,)*8+(-4,-4,7)
powers=[[1]]
for t in range(A):
    nxt=[0]*(len(powers[-1])+1)
    for j,v in enumerate(powers[-1]):
        if j==0:nxt[1]+=v
        else:
            nxt[j-1]+=v;nxt[j]+=v;nxt[j+1]+=v
    powers.append(nxt)
@lru_cache(None)
def moment2(t,i,j):
    if (i|j)&1:return 0
    i//=2;j//=2
    return sum(v*(powers[k][i] if i<len(powers[k]) else 0)*
               (powers[t-k][j] if j<len(powers[t-k]) else 0)
               for k,v in enumerate(choose(t)))
phi=entry(moment2,A,C)
assert phi>0
ds={}
for z in (1,2,-3,-4,7):
    R=list(C);R.remove(1);a=A
    if z==2:a-=1
    else:R.remove(z)
    d=(1 if z>0 else -1)*entry(moment2,a,R,1,abs(z))
    difference=phi-entry(moment2,a,R+[-1,-z])
    assert difference==4*d and 16*d < -phi
    ds[z]=d
d22=entry(moment2,A-2,C,2,2)
assert 0<4*d22<phi
assert phi-entry(moment2,A-2,C+(-2,-2))==4*d22
assert 5+2*A+3*8+4*2==1061
print("r=5, W=1061, p=7: all five incident flips < -Phi/16; D22>0 PASS")

if not args.no_census:
    path=Path("/tmp/claude-1006/-home-yang-q3adjoint/"
              "2d613d32-0be2-46ab-91e8-7dc4d59487ae/"
              "scratchpad/flipdesc/fx6_52_noflip.log")
    counts=Counter();N=0
    for line in path.read_text().splitlines():
        m=re.search(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?) phi=(\d+)",line)
        if not m:continue
        B=tuple(map(int,m[3].split()));L=B+(int(m[2]),)
        assert sum(map(abs,B))==int(m[1])<=52
        assert all(-z not in L for z in L)
        mult=max(Counter(L).values())
        assert mult<=4
        counts[mult]+=1;N+=1
    assert N==75532
    assert dict(counts)=={1:53864,2:20208,3:1406,4:54}
    print("census:",N,dict(sorted(counts.items())),"PASS")
print("ALL EXACT CHECKS PASS")
