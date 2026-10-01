import argparse,itertools,random
from fractions import Fraction as Q
from math import comb

argparse.ArgumentParser(
    description='FM-MECH69 exact verifier; memory only'
).parse_args()

# Constants in the even-minus derivative and decay argument.
assert 4*Q(3,4)**3>=1
assert Q(3,2)/Q(3,4)**2<=3
assert Q(3,2)**2<=16*Q(3,4)**5
assert 2*3**2-3*3-4==5
assert 64<75 and 30<32

def U(n,x):
    a,b=Q(1),x
    if n==0:return a
    for _ in range(2,n+1):
        a,b=b,x*b-a
    return b

def der(n,t):
    a,b,da,db=Q(1),2*t,Q(0),Q(2)
    for _ in range(2,n+1):
        a,b,da,db=b,2*t*b-a,db,2*b+2*t*db-da
    return db

def R(n,e,t,q):
    den=n+1+e*U(n,2*q)
    if den:
        return (U(n,2*t)+e*U(n,2*q*t))/den
    return t*der(n,t)/der(n,Q(1))

checks=0
for n in range(6,81,2):
    M=n+1
    for q in [Q(j,12) for j in range(-12,13)]:
        assert M-U(n,2*q)>=Q(M,3)*(1-q*q)
        for t in (
            Q(0),Q(1,M*M),Q(1,M),Q(1,16),
            Q(1,4),Q(1,2),Q(3,4),Q(1)
        ):
            r=R(n,-1,t,q)
            assert r*r<=32*t*t
            if t<=Q(1,2):
                assert abs(r)<=6*M*t*t
                assert r*r*M*M*(1-t*t)<=16
            checks+=1
print('Even-minus profile checks:',checks)

cases=0
for p1,m1,m2,m3,m4,lm in itertools.product(range(3),repeat=6):
    V=Q(p1+m1,2)+m2
    E=m1+m2+m3+m4+lm
    if E<2 or E%2:
        continue
    if V<1:
        if m4:
            raw2=9
            removed=1
        else:
            need=int(2*(1-V))
            taken=min(need,m3)
            raw2=4**taken
            removed=taken
            need-=taken
            assert need<=lm
            raw2*=32**need
        assert raw2**3*3**int(2*V+2*removed)<=32**6
    cases+=1
print('Consumer resource configurations:',cases)

near=Q(663552,2**32)
far1=Q(4608*32*256**2,2**42)
far2=Q(13824*32*257**2,2**42)
assert all(v<Q(1,4) for v in (near,far1,far2))
assert 360**2<2**17-256
assert near+Q(2,3)<1
assert 380*36>17*256
assert 4096*36>=2**17
assert 384*341<2**17<=384*342
print('Tail certificates:',near,far1,far2)
print('Cutoff: max(384*S, 131072); homogeneous: 4096*S')

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

def coefficient(bg,n):
    cap=sum(j for j,e in bg)-n
    if cap<0 or cap%2:
        return 0
    p={(0,0):1}
    for j,e in bg:
        a={(2*h,0):1 for h in range(j+1)}
        for h,c in mono(j).items():
            a[j,h]=a.get((j,h),0)+e*c
        p=mul(p,a,cap)
    return 2*sum(
        c*moment(j)*((i==cap)-(i==cap-2))
        for (i,j),c in p.items()
    )

rng=random.Random(69)
least=None
positive=0
for _ in range(80):
    ns=sorted(
        rng.randrange(5,13)
        for _ in range(rng.randrange(2,5))
    )
    N=ns.pop()
    bg=[(j,rng.choice((-1,1))) for j in ns]
    bg += [
        (rng.randrange(1,5),rng.choice((-1,1)))
        for _ in range(rng.randrange(7))
    ]
    e=(-1)**sum(s==-1 for j,s in bg)
    value=direct(bg+[(N,e)])
    assert value==coefficient(bg,N) and value>=0
    if value:
        positive+=1
        least=value if least is None else min(least,value)
print('Multi-label bridges:',80,'positive:',positive,'least:',least)

for ell in range(2,17):
    n=4*ell
    assert Q(4*mono(n)[2],n)==-(2*ell+1)
print('Unbounded extracted even-minus profiles: 15 exact checks')

for Ms,u in (
    ((4,5,6,23),Q(1,7)),
    ((6,6,6,30),Q(3,2)),
    ((6,7,11),Q(1,5))
):
    top=max(Ms)
    lhs=Q(1)
    for M in Ms:
        lhs*=(1+Q(M*M,4096)*u)**(top*top)
    rhs=(1+Q(top*top,4096)*u)**sum(M*M for M in Ms)
    assert lhs>=rhs
print('Concavity product certificates: 3')

for m in (5,6,10,100,10000):
    S=(m+1)**2+36*m
    assert Q(S,(m+1)**2)<=6<2**21
    assert 5*(m+1)**2-36*m==(m-5)*(5*m-1)
print('Unbounded ordinary-count family: PASS')
print('PASS')
