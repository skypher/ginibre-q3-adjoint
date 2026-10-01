import argparse
argparse.ArgumentParser(
    description="FM-MECH121 uniform bounds; exact arithmetic, memory only"
).parse_args()
from fractions import Fraction as Q
from math import comb, factorial
from functools import lru_cache
from random import Random

def fusion(ns):
    a={0:1}
    for n in ns:
        b={}
        for j,v in a.items():
            for h in range(abs(j-n),j+n+1,2):
                b[h]=b.get(h,0)+v
        a=b
    return a

@lru_cache(None)
def mom(r,n=0):
    return fusion((1,)*r).get(n,0)

def crow(t,x,d):
    c=[Q(1),Q(x)]
    for j in range(1,d):
        c.append((x*c[j]-(t-j+1)*c[j-1])/(j+1))
    return c

def alpha(l,j,n,t):
    if j==l:
        return Q(n+1,l+n+1)
    v=Q((n+1)*(t-2*j-n)*factorial(j+n),
        factorial(l+n+1))
    for h in range(j,l-1):
        v*=t-h
    return v

checks=0
for t in range(1,19):
    for x in range(-t,t+1,2):
        c=crow(t,x,14)
        for l in range(6):
            for n in range(6):
                W=c[l]*c[l+n]-(c[l-1]*c[l+n+1] if l else 0)
                rhs=sum(alpha(l,j,n,t)*c[j]*c[j+n]
                        for j in range(l+1))
                assert W==rhs
                checks+=1
print("same-gap identities",checks,"PASS",flush=True)

checks=0
for d in range(1,17):
    for t in (2*d,3*d+7,11*d+7):
        for x in (-t,-t//2,0,t//3,t):
            c=crow(t,x,d+1)
            P=[c[j]**2-(c[j-1]*c[j+1] if j else 0)
               for j in range(d+1)]
            assert all(v>=0 for v in P)
            for j in range(1,d+1):
                assert (t-2*j+2)*P[j-1] <= (j+1)*P[j]
                checks+=1

for m in range(25):
    for n in range(25):
        if m>=n and (m-n)%2==0:
            assert mom(m,n)==Q(
                2*(n+1)*comb(m,(m-n)//2),m+n+2)
        else:
            assert mom(m,n)==0
print("Turan ratios",checks,
      "and kernel moment identities PASS",flush=True)

def cost(k,s):
    r=s*s
    return ((1+r**4+4*s**3*(1+r))**k-1)/(1-r)**(k-1)

assert cost(3,Q(1,3))==Q(163864516807,223154201664)<Q(3,4)
assert cost(4,Q(1,4))==Q(
    72151010412074339,202661983231672320)<Q(3,8)
assert Q(1081,6561)<Q(1,6)
assert Q(9,8)*Q(1,5)==Q(9,40)
print("uniform core-budget constants PASS",flush=True)

def mul(a,b,d):
    out={}
    for (i,j),v in a.items():
        for (h,l),w in b.items():
            if i+h+2*(j+l)<=2*d:
                m=(i+h,j+l)
                out[m]=out.get(m,0)+v*w
    return {m:v for m,v in out.items() if v}

def power(a,n,d):
    b={(0,0):1}
    while n:
        if n&1:
            b=mul(b,a,d)
        n//=2
        if n:
            a=mul(a,a,d)
    return b

def block(n,e,d):
    a={(0,j):1 for j in range(min(n,d)+1)}
    if n<=2*d:
        for j in range(n//2+1):
            m=(n-2*j,j)
            a[m]=a.get(m,0)+e*(-1)**j*comb(n-j,j)
    return {m:v for m,v in a.items() if v}

def moments(poly,d):
    P=[0]*(d+1)
    for (h,l),v in poly.items():
        if h%2==0:
            P[h//2+l]+=v*(comb(h,h//2)//(h//2+1))
    return P

def small(A,E,b,d):
    p=power(block(1,1,d),A,d)
    p=mul(p,power(block(1,-1,d),E,d),d)
    return mul(p,power(block(2,1,d),b,d),d)

def consumer(base,cores,d):
    p=base
    for n,e in cores:
        p=mul(p,block(n,e,d),d)
    P=moments(p,d)
    return P[d]-P[d-1]

rng=Random(121)
checks=0
for d in range(1,13):
    for _ in range(5):
        a=rng.randrange(2*d+8)
        b=rng.randrange(2*d+8)
        P=moments(small(a,0,b,d),d)
        for j in range(1,d+1):
            if a+b-2*j+2>=0:
                assert (a+b-2*j+2)*P[j-1]<=j*P[j]
                checks+=1
print("plus-background insertion ratios",checks,"PASS",flush=True)

kernel_checks=0
for _ in range(12):
    d=8
    base=small(rng.randrange(7),rng.randrange(7),
               rng.randrange(7),d)
    P=moments(base,d)
    assert all(v>=0 for v in P)
    for n in range(d+1):
        for r in range(d-n+1):
            T=sum(v*mom(h,n) for (h,j),v in base.items()
                  if h+2*j==n+2*r)
            assert T*T<=(n+1)**2*P[r]*P[r+n]
            kernel_checks+=1
print("integrated kernel controls",kernel_checks,"PASS",flush=True)

for k,d in ((3,5),(3,8),(3,12),(4,5),(4,8),(4,12)):
    ns=[rng.randrange(3,d+4) for _ in range(k)]
    es=[rng.choice((-1,1)) for _ in range(k)]
    cores=list(zip(ns,es))
    if k==3:
        t=11*d+7
        ab=11*d-2
        aa=4*d+2
        margin=Q(1,4)
    else:
        t=18*d+14
        ab=18*d-2
        aa=5*d+3
        margin=Q(5,8)
    for tag,A,E,b in (
        ("minus",t//2,t-t//2,0),
        ("plus_count",ab//3,0,ab-ab//3),
        ("plus_Newton",aa,0,2)):
        base=small(A,E,b,d)
        P=moments(base,d)
        F=consumer(base,cores,d)
        assert F>=margin*P[d]>0,(k,d,tag,F,P[d])
    print("consumer k,d",k,d,
          "three uniform sectors PASS",flush=True)

# Mixed forward-difference multipliers.
for k,d in ((3,6),(3,10),(3,16),(4,6),(4,10),(4,16)):
    aa=(4*d+2) if k==3 else (5*d+3)
    margin=Q(1,4) if k==3 else Q(5,8)
    cores=[(rng.randrange(3,d+3),rng.choice((-1,1)))
           for _ in range(k)]
    for _ in range(5):
        j=rng.randrange(d+1)
        i=rng.randrange(2*d-2*j+1)
        H=power({(1,0):1,(0,1):1},i,d)
        H=mul(H,power({(2,0):1,(0,2):1},j,d),d)
        base=mul(power(block(1,1,d),aa,d),H,d)
        P=moments(base,d)
        F=consumer(base,cores,d)
        assert F>=margin*P[d]>=0
print("30 mixed Newton coefficients at uniform shifts PASS",flush=True)

for kk in range(1,61):
    r=Q(1,9*kk)
    upper_u=r**4+4*r*Q(1,3)*(1+r)
    assert kk*upper_u<Q(1,6)
    assert (1-r)**(kk-1)>=Q(8,9)
    for n in range(3,31):
        assert (n+1)**2*r**n<=16*r**3
print("general-k rational majorants PASS",flush=True)

def W(a,l,n):
    if l<0 or l+n>a:
        return 0
    return Q((n+1)*comb(a+1,l)*comb(a+1,l+n+1),a+1)

for a in range(1,21):
    c=crow(a,a,a+1)
    for l in range(a+1):
        for n in range(a-l+1):
            assert W(a,l,n)==(
                c[l]*c[l+n]-(c[l-1]*c[l+n+1] if l else 0))
print("binomial Schur formula PASS",flush=True)

# Fixed-reference obstruction at the requested shift.
C={(0,0):1,(0,1):-1}
for _ in range(3):
    C=mul(C,block(3,1,12),12)

def lower_ratio(N,K,h):
    v=Q(1)
    for j in range(h):
        v*=Q(K-j,N-K+j+1)
    return v

def wratio(d,l,n):
    if l<0:
        return Q(0)
    return ((n+1)*lower_ratio(2*d+5,d,d-l)
            *lower_ratio(2*d+5,d+1,d-l-n))

def endpoint(d):
    return sum(
        v*sum(m*wratio(d,d-j-(h+n)//2,n)
              for n,m in fusion((1,)*h).items())
        for (h,j),v in C.items())

for d in (16,64,256):
    a=2*d+4
    F=sum(
        v*sum(m*W(a,d-j-(h+n)//2,n)
              for n,m in fusion((1,)*h).items())
        for (h,j),v in C.items())
    assert endpoint(d)==F/W(a,d,0)

ratios=[endpoint(d) for d in (256,4096,2**20,2**24)]
assert all(r>0 for r in ratios)
assert all(ratios[j]>ratios[j+1] for j in range(3))
assert ratios[-1]<Q(1,1000)
print("endpoint d=2^24: 0 < F/P_d < 1/1000 PASS",flush=True)
print("ALL CHECKS PASS",flush=True)