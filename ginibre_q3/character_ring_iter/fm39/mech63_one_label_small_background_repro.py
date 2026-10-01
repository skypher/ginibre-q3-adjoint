import argparse, itertools, random
from fractions import Fraction as Q
from math import comb
import sympy as sp

argparse.ArgumentParser(
    description='FM-MECH63 exact verifier; no file output'
).parse_args()

s,t,z=sp.symbols('s t z')
assert sp.expand(
    1-(1-s)*(1-3*s)**2-sp.Rational(3,4)*s
    -s*(3*s-sp.Rational(5,2))**2
)==0
assert sp.expand(
    1-(1-s)*(1-sp.Rational(6,5)*s)**2-sp.Rational(3,4)*s
    -s*(sp.Rational(1,4)+sp.Rational(24,25)*(1-s)
        +sp.Rational(36,25)*(1-s)**2)
)==0
assert sp.expand(
    (1-s)*(1-4*s)+1-s/3
    -4*(s-sp.Rational(2,3))**2-sp.Rational(2,9)
)==0
D=6-12*t+16*t*t
P=16*(1+t*t)*z*z-12*(1+t)*z+2
assert sp.expand(
    D-16*(t-sp.Rational(3,8))**2-sp.Rational(15,4)
)==0
assert sp.factor(
    P+sp.Rational(5,2)
    -16*(1+t*t)*(z-3*(1+t)/(8*(1+t*t)))**2
    -9*(t-1)**2/(4*(1+t*t))
)==0

bounds=(
    Q(663552,2**32),
    Q(55296*256**2,2**42),
    Q(165888*257**2,2**42)
)
assert all(v<Q(1,4) for v in bounds)
assert 16*90**2<2**17-256
assert 380*36>17*256
assert 4096*36>=2**17
print('Polynomial and cutoff certificates: PASS',bounds)

def U(n,x):
    a,b=Q(1),x
    if n==0:
        return a
    for _ in range(2,n+1):
        a,b=b,x*b-a
    return b

def dP(n,t):
    a,b,da,db=Q(1),2*t,Q(0),Q(2)
    for _ in range(2,n+1):
        a,b,da,db=b,2*t*b-a,db,2*b+2*t*db-da
    return db

def R(n,e,t,q):
    den=n+1+e*U(n,2*q)
    if den:
        return (U(n,2*t)+e*U(n,2*q*t))/den
    return t*dP(n,t)/dP(n,Q(1))

small=odd=0
for n in (3,4):
    for q,e,t in itertools.product(
        [Q(j,8) for j in range(-8,9)],
        (-1,1), [Q(j,16) for j in range(17)]
    ):
        assert abs(R(n,e,t,q))<=1-(1-t*t)/3
        small+=1

for n in range(5,66,2):
    for q,e,t in itertools.product(
        [Q(j,8) for j in range(-8,9)], (-1,1),
        (Q(1,1000),Q(1,n+1),Q(1,8),Q(1,4),
         Q(1,2),Q(3,4),Q(1))
    ):
        assert abs(R(n,e,t,q))<=4*t
        odd+=1
print('Small profiles:',small,'odd exceptional profiles:',odd)

resources=0
for p1,m1,m2,p3,m3,m4 in itertools.product(range(3),repeat=6):
    V=Q(p1+m1,2)+m2
    for n,e in itertools.product((5,6),(-1,1)):
        E=m1+m2+m3+m4+(e==-1)
        A=p1+m2+p3+m4+int(
            (n%2 and e==1) or (n%2==0 and e==-1)
        )
        if E==0 or E%2 or A%2:
            continue
        N=Q(A+E,2)
        nu=Q(1,2) if n%2 else Q(e==-1)
        assert N==V+Q(p3+m3,2)+m4+nu

        if V<1:
            need=1-V
            raw=1
            removed=0
            if m4:
                need-=1
                raw*=3
                removed+=1
            else:
                for _ in range(min(p3+m3,2)):
                    if need<=0:
                        break
                    need-=Q(1,2)
                    raw*=2
                    removed+=1
            if need>0:
                assert n%2==1 and e==-1
                need-=Q(1,2)
                raw*=4
            assert need<=0
            assert raw**6*3**int(2*V+2*removed)<=12**6
        resources+=1
print('Consumer resource cases:',resources)

def mono(n):
    return {
        n-2*j:(-1)**j*comb(n-j,j)
        for j in range(n//2+1)
    }

def mul(a,b,cap=None):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            if cap is None or i+k<=cap:
                c[i+k,j+l]=c.get((i+k,j+l),0)+v*w
    return {ij:v for ij,v in c.items() if v}

def pair(n,e):
    p={}
    for j,c in mono(n).items():
        p[j,0]=p.get((j,0),0)+c
        p[0,j]=p.get((0,j),0)+e*c
    return p

def moment(j):
    return 0 if j%2 else comb(j,j//2)//(j//2+1)

def direct(word):
    p={(0,0):1}
    for n,e in word:
        p=mul(p,pair(n,e))
    return sum(
        c*moment(i)*moment(j)
        for (i,j),c in p.items()
    )

def coefficient(background,n):
    cap=sum(j for j,e in background)-n
    if cap<0 or cap%2:
        return 0
    p={(0,0):1}
    for j,e in background:
        a={(2*k,0):1 for k in range(j+1)}
        for k,c in mono(j).items():
            a[j,k]=a.get((j,k),0)+e*c
        p=mul(p,a,cap)
    return 2*sum(
        c*moment(j)*((i==cap)-(i==cap-2))
        for (i,j),c in p.items()
    )

rng=random.Random(63)
least=None
for _ in range(180):
    bg=[
        (rng.randrange(1,5),rng.choice((-1,1)))
        for _ in range(rng.randrange(1,10))
    ]
    n=rng.randrange(5,15)
    e=(-1)**sum(sg==-1 for j,sg in bg)
    value=direct(bg+[(n,e)])
    assert value==coefficient(bg,n) and value>=0
    if value>0:
        least=value if least is None else min(least,value)
print('Distance/Catalan bridges: 180; least positive:',least)

h,n=4,20
p={
    h+1+j:(-1)**j*comb(2*h,j)*2**j
    for j in range(2*h+1)
}
u={j//2:c*2**j for j,c in mono(n).items()}
mu=lambda j:Q(1,(j+1)*(j+2))
ray=(
    sum(a*b*mu(i+j) for i,a in p.items() for j,b in u.items())
    /sum(a*mu(i) for i,a in p.items())
)
assert ray==-Q(3277672,5980345)
assert direct([(1,-1)]*2+[(3,1)]*8+[(20,1)])==406
print('Ray/full witness:',ray,406)

term=upper=Q(1)
for j in range(40):
    term*=-Q(20*(j+2),(2*j+3)*(2*j+2))
    upper+=term
assert upper<-Q(1,50)

for n in (5,10,18,20,100):
    print('n',n,'B cutoff',max(384*(n+1)**2,2**17))
print('PASS')
