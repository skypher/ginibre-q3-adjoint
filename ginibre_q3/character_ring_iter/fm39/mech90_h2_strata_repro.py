import numpy as np
from fractions import Fraction as Q
from math import comb
from collections import defaultdict
from itertools import product
from random import Random

# Finite remainder: balanced backgrounds with a+r < 44.
S=91
cat=[comb(2*j,j)//(j+1) for j in range(134)]
ct=np.array(cat[:S],dtype=object)
ce=np.array([cat[j+2]-3*cat[j+1]+cat[j]
             for j in range(S)],dtype=object)
one=[(0,0,1),(1,0,2),(2,0,1),(1,1,-1)]
three=[(j,0,j+1 if j<=3 else 7-j) for j in range(7)]
three += [(3,1,-4),(3,2,4),(3,3,-1)]

def times(A,terms):
    B=np.zeros((S,S),dtype=object)
    for i,j,c in terms:
        B[i:,j:]+=c*A[:S-i,:S-j]
    return B

A=np.zeros((S,S),dtype=object);A[0,0]=1
Va=np.array([1],dtype=object)
profiles=checks=0
for a in range(44):
    P=A.copy();V=Va.copy()
    for r in range(44-a):
        profiles+=1
        C=P@ct;E=P@ce
        p1=sum(v*cat[j] for j,v in enumerate(V))
        e1=sum(v*(cat[j+2]-3*cat[j+1]+cat[j])
               for j,v in enumerate(V))
        delta=(3*a+1)//2+4*r
        for sg in (-1,0,1):
            D=2*a+6*r+(4 if sg else 0)
            total=5*p1+sg*e1 if sg else p1
            prefix=0
            for g in range(D-delta+1):
                if g:
                    j=g-1
                    prefix += (sum(C[j-k] for k in range(min(4,j)+1))
                               if sg else C[j])
                    if sg and j>=2:
                        prefix+=sg*E[j-2]
                assert total>=prefix,(a,r,sg,D-g)
                checks+=1
        if r+1<44-a:
            P=times(P,three)
            V=np.convolve(V,np.array([16,-4,4,-1],dtype=object))
    if a<43:
        A=times(A,one)
        Va=np.convolve(Va,np.array([4,-1],dtype=object))
expected=sum(3*(a//2+2*r+1)+8
             for a in range(44) for r in range(44-a))
assert (profiles,checks)==(990,116589)==(990,expected)
print("Finite remainder:",profiles,checks,"PASS")

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def epoly(word,cap=None):
    D=sum(n for n,e in word)
    if cap is None:cap=2*D
    T={(0,0):1}
    for n,e in word:
        N=defaultdict(int)
        for (a,j),v in T.items():
            for k in range(n+1):
                if a+2*k<=cap:N[a+2*k,j]+=v
            if a+n<=cap:
                for k in range(abs(j-n),j+n+1,2):
                    N[a+n,k]+=e*v
        T={ij:v for ij,v in N.items() if v}
    return [T.get((2*j,0),0) for j in range(cap//2+1)]

def gram(word):
    D=sum(n for n,e in word);out=[0]*(D+1)
    choices=[list(range(0 if n%2==0 and e>0
                        else (2 if n%2==0 else 1),n+1,2))
             for n,e in word]
    for ks in product(*choices):
        a=[1];weight=1
        for (n,e),k in zip(word,ks):
            f=[0]*(n+1);sh=(n-k)//2;f[sh]=1
            if k:f[sh+k]=e
            else:weight*=2
            a=mul(a,f)
        for j in range(D+1):
            out[j]+=weight*(a[j]**2
                -(a[j-1] if j else 0)*(a[j+1] if j<D else 0))
    return out

rng=Random(9003)
for _ in range(100):
    word=[(rng.randrange(1,7),rng.choice((-1,1)))
          for _ in range(rng.randrange(1,7))]
    assert epoly(word)==gram(word)
print("100 half-angle identities: PASS")

# Exact identities of quadratic polynomials.
def add(q,i,j,v):
    if v:q[tuple(sorted((i,j)))]+=v
def clean(q):
    return {ij:v for ij,v in q.items() if v}
def energy_form(D,d,sg):
    out=defaultdict(int)
    def fold(j):
        if j<0 or j>D or (2*j==D and sg<0):return None,0
        return (j,1) if 2*j<=D else (D-j,sg)
    for j in range(d+1):
        i,s=fold(j)
        if s:add(out,i,i,2)
        i,s=fold(j-1);k,t=fold(j+1)
        if s*t:add(out,i,k,-2*s*t)
    return clean(out)
def sos(terms):
    out=defaultdict(int)
    for weight,f in terms:
        ids=sorted(f)
        for ii,i in enumerate(ids):
            for j in ids[ii:]:
                add(out,i,j,weight*f[i]*f[j]*(1 if i==j else 2))
    return clean(out)
def unit(i):return {i:1}

cnt=centers=0
for D in range(8,241):
    m=D//2
    if D%2:
        for sg in (-1,1):
            terms=[(1,unit(0)),(1,unit(1))]
            terms += [(1,{i:1,i+2:-1}) for i in range(m-1)]
            terms += [(1,{m:1,m-1:-sg})]
            assert energy_form(D,m,sg)==sos(terms)
            cnt+=1
    else:
        b,c,d,e=m-2,m-3,m-4,m-1
        base=[(1,unit(i)) for i in (0,1,b,e)]
        base += [(1,{i:1,i+2:-1}) for i in range(m-2)]
        for off in (-1,0,1,2,3):
            terms=base[:]
            if off==0:terms.append((2,unit(e)))
            if off==1:terms.append((4,unit(e)))
            if off==2:
                terms.remove((1,unit(e)))
                terms.remove((1,{c:1,e:-1}))
                terms += [(1,{c:1,e:-2}),(2,unit(e)),(2,unit(b))]
            if off==3:
                for term in ((1,unit(e)),(1,{c:1,e:-1}),
                             (1,unit(b)),(1,{d:1,b:-1})):
                    terms.remove(term)
                terms += [(1,{c:1,e:-2}),(2,unit(e)),
                          (2,unit(c)),(1,{d:1,b:-2})]
            assert energy_form(D,m+off,-1)==sos(terms)
            cnt+=1
        for j in (m-1,m):
            Z=defaultdict(int,energy_form(D,j,-1))
            for ij,v in energy_form(D,j-1,-1).items():Z[ij]-=v
            assert clean(Z)=={(e,e):2}
            centers+=1
assert (cnt,centers)==(817,234)
print("817 SOS and 234 central identities: PASS")

# Rational complex checks of the ellipse lemma.
def cmul(z,w):
    return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def Uc(n,x):
    a=(Q(1),Q(0));b=x
    if not n:return a
    for _ in range(n-1):
        c=cmul(x,b);a,b=b,(c[0]-a[0],c[1]-a[1])
    return b
def Ur(n,x):
    a,b=Q(1),x
    if not n:return a
    for _ in range(n-1):a,b=b,x*b-a
    return b

ctest=0
for s in (Q(1,4),Q(1,2),Q(3,4),Q(1)):
    for t in (Q(0),Q(1,4),Q(1,2),Q(1),Q(2),Q(4)):
        co=(1-t*t)/(1+t*t);si=2*t/(1+t*t)
        x=((s+1/s)*co,(s-1/s)*si)
        for eta in (Q(j,2) for j in range(-4,5)):
            for n in range(1,5):
                z=Uc(n,x);z=cmul(z,z)
                z=(z[0]-Ur(n,eta)**2,z[1])
                assert z[0]**2+z[1]**2<=Ur(n,s+1/s)**4
                ctest+=1
assert ctest==864
print("864 rational ellipse checks: PASS")

q=Q(16,25);s=Q(4,5)
G=lambda n:sum(q**j for j in range(n+1))
gamma=G(1)**2/(4*s)
gamma3=G(3)**2/(16*q*q)
C=Q(2,3)*(G(4)+5*q*q)/((1-q)*q**3)
assert gamma==Q(1681,2000) and gamma3<gamma
assert C==Q(1768561,55296)
assert C*45*gamma**44<1
assert Q(46,45)*gamma<1
assert Q(8,3)*q/(1-q)<C
print("H>=44 induction: PASS")

def p1(r):
    a=[1]
    for _ in range(r):a=mul(a,[16,-4,4,-1])
    a=mul(a,[6,-3,1])
    return sum(v*comb(2*j,j)//(j+1) for j,v in enumerate(a))

r=100;D=604;delta=400;P1=p1(r)
err=G(3)**(2*r)*(G(4)+5*q*q)/((1-q)*q**(D-delta-1)*P1)
assert 0<err<Q(1,10**6)
old=Q(10*(D+1)**2,3*(2*D+2-2*r))*16**(D-delta-1)
old*=Q(4625,16384)**(2*r)*Q(71185,262144)
assert old>1
print("D=604, delta=400: new error < 10^-6; old error > 1")

# Actual negative channels in the unresolved family.
S=[1,-1,0,1,-2,1,0,-1,1]
for m in range(1,11):
    r=6*m+3;R=r-2;a=[0]*(r+2)
    for t in range(m+1):
        for j,v in enumerate(S):
            if 4+6*t+j<len(a):
                a[4+6*t+j]+=(-1)**t*comb(R,t)*v
    value=sum(a[j]**2-(a[j-1] if j else 0)*a[j+1]
              for j in range(r+1))
    expected=sum(11*comb(R,t)**2-7*comb(R,t)*comb(R,t+1)
                 for t in range(m-1))
    expected+=8*comb(R,m-1)**2-2*comb(R,m-1)*comb(R,m)
    assert value==expected<0
    if m==2:assert value==-756

word=[(3,1)]*15+[(3,-1)]*15+[(4,1)]
assert sum(epoly(word,30))==309937000389
central=epoly(word,94)
assert central[46]==central[47]==198943471105586480
print("Negative-channel control: full V =",2*309937000389)
print("h=1 central plateau: PASS")