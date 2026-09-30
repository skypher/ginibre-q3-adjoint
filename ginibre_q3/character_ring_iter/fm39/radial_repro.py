
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb
from fractions import Fraction

argparse.ArgumentParser(
    description="Exact bounded checks of the whole-word radial bound."
).parse_args()

def mul(f,g):
    out=defaultdict(int)
    for (i,j),c in f.items():
        for (k,l),d in g.items():
            out[i+k,j+l]+=c*d
    return {ij:c for ij,c in out.items() if c}

@lru_cache(None)
def h(k):
    out=defaultdict(int)
    for b in range((k+1)//2+1):
        degree=k+1-2*b
        c=(-1)**b*comb(k+1-b,b)
        for j in range(degree):
            out[degree-1-j,j]+=c
    return dict(out)

@lru_cache(None)
def fund(n):
    return {(n-j,j):comb(n,j) for j in range(n+1)}

@lru_cache(None)
def mom(n):
    return 0 if n%2 else comb(n,n//2)//(n//2+1)

@lru_cache(None)
def ker(r,i,j):
    return sum((-1)**b*comb(2*r,b)
               *mom(i+2*r-b)*mom(j+b)
               for b in range(2*r+1))

def phi(r,f):
    return Fraction(sum(c*ker(r,i,j)
                        for (i,j),c in f.items()),2)

def quotient_sd(core,t):
    # x=(s+d)/2, y=(s-d)/2; remove the polynomial factor s^t.
    out=defaultdict(Fraction)
    for (i,j),c in core.items():
        for b in range(i+1):
            for e in range(j+1):
                out[i+j-b-e,b+e]+=Fraction(
                    c*comb(i,b)*comb(j,e)*(-1)**e,2**(i+j))
    ans={}
    for (sdeg,ddeg),c in out.items():
        if not c:
            continue
        assert ddeg%2==0 and sdeg>=t and (sdeg-t)%2==0
        ans[(sdeg-t)//2,ddeg//2]=c
    return ans

def moment_ratio(r,m,dr,dm):
    # M(r,m)=phi_r(h1^(2m)); return M(r+dr,m+dm)/M(r,m).
    out=Fraction(1)
    for _ in range(dr):
        out*=Fraction(4*(2*r+1)*(2*r+3),
                      (r+m+2)*(r+m+3))
        r+=1
    for _ in range(dm):
        out*=Fraction(4*(2*m+1)*(2*m+3),
                      (r+m+2)*(r+m+3))
        m+=1
    return out

def normalized_value(r,a,t,F):
    m=(a+t)//2
    return sum(c*moment_ratio(r,m,j,i)
               for (i,j),c in F.items())

moment_checks=0
for r in range(7):
    for m in range(7):
        assert phi(r,fund(2*m)) == (
            Fraction(1,2)*moment_ratio(0,0,r,m))
        moment_checks+=1

bridge_checks=lower_checks=threshold_checks=0
for ks in combinations_with_replacement(range(1,6),3):
    t=sum(k%2 for k in ks)
    d=(sum(ks)-t)//2
    B=Fraction(1)
    core={(0,0):1}
    for k in ks:
        core=mul(core,h(k))
        B*=Fraction(comb(k+3,3),4 if k%2 else 1)
    assert B.denominator==1
    F=quotient_sd(core,t)
    K=16*B*sum(j*j for j in range(1,d+1))

    for r in (2,3,4):
        for a in range(6):
            if (a+t)%2:
                continue
            ratio=normalized_value(r,a,t,F)
            direct=phi(r,mul(core,fund(a)))
            base=phi(r,fund(a+t))
            assert ratio==direct/base
            assert ratio>=1-Fraction(K,2*r+a+t+6)
            bridge_checks+=1
            lower_checks+=1

    a=1 if t%2 else 2
    r=max(2,(int(K)-a-t-6+1)//2)
    ratio=normalized_value(r,a,t,F)
    assert 2*r+a+t+6>=K
    assert ratio>=1-Fraction(K,2*r+a+t+6)>=0
    threshold_checks+=1

print("fundamental product-formula checks:",moment_checks)
print("whole-word quotient/moment bridges:",bridge_checks)
print("radial lower-bound checks:",lower_checks)
print("bounded-label exact threshold checks:",threshold_checks)
