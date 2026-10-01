from fractions import Fraction as Q
from math import comb, prod
from functools import lru_cache
from collections import defaultdict
from random import Random

rho=(Q(25,32),Q(321,512),Q(4625,16384),Q(71185,262144))
for n,r in enumerate(rho,1):
    assert r<=Q(4,5)**n
assert Q(19,10)**2<Q(15,4)
assert Q(19,10)*Q(7,44)>Q(3,10)
E=lambda D,g:Q(10*(D+1),3)*16**(g-1)*Q(4,5)**D
assert E(61,3)<1 and E(64,4)<Q(3,5)
assert Q(81,65)*16*Q(4,5)**16<1

def qmult(ns,j):
    a,b,c,d=ns
    il,ih=abs(a-b),a+b
    jl,jh=abs(c-d),c+d
    if (sum(ns)-j)%2:return 0
    ans=0
    for off in range(-j,j+1,2):
        lo=max(il,jl-off,(j-off+1)//2)
        hi=min(ih,jh-off)
        lo+=(il-lo)%2
        if lo<=hi:ans+=(hi-lo)//2+1
    return ans

rng=Random(8704)
checks=0
for _ in range(2000):
    ns=sorted(rng.randrange(1,81) for _ in range(4))
    sigma=sum(ns)%2
    A=max(Q(0),1+min(Q(ns[0]),Q(sum(ns)-2*ns[-1],2)))
    assert A==Q(qmult(ns,sigma),sigma+1)
    for j in range(sigma,22,2):
        v=j//2
        assert abs(qmult(ns,j)-(j+1)*A)<=Q(3,2)*v*(v+1)
        checks+=1
assert checks==22000

def bgcoef(word):
    a={(0,0):1}
    for n,e in word:
        b=defaultdict(int)
        for (j,k),v in a.items():
            for t in range(abs(j-n),j+n+1,2):b[t,k]+=v
            for t in range(abs(k-n),k+n+1,2):b[j,t]+=e*v
        a={ij:v for ij,v in b.items() if v}
    return a

bg=[(1,1)]*2+[(1,-1)]*2+[(3,1),(4,-1)]
B=bgcoef(bg);D=11
P1=sum((j+1)*v for (j,k),v in B.items() if k==0)
C0=prod(2*(n+1) for n,e in bg)
H=sum(Q((j-1)*(j+1),4)*v
      for (j,k),v in B.items() if k==0 and j%2)
K=sum(v for (j,k),v in B.items() if j%2 and k%2==0)
assert (P1,C0,H,K)==(196,20480,244,22)

@lru_cache(None)
def inv(ns):
    if not ns:return 1
    if sum(ns)%2 or 2*max(ns)>sum(ns):return 0
    if len(ns)==1:return 0
    L=len(ns);d=sum(ns)//2;out=0
    for mask in range(1<<L):
        k=d-sum(ns[j]+1 for j in range(L) if mask>>j&1)
        if k>=0:
            out+=(-1)**mask.bit_count()*comb(k+L-2,L-2)
    return out

def even(word):
    ns=tuple(n for n,e in word);L=len(ns);full=(1<<L)-1
    m=[inv(tuple(sorted(ns[j] for j in range(L) if S>>j&1)))
       for S in range(1<<L)]
    neg=sum(1<<j for j,(n,e) in enumerate(word) if e<0)
    return sum((-1)**((S&neg).bit_count())*m[S]*m[full^S]
               for S in range(1<<L))

for family,n in ((0,784),(1,575),(0,12),(1,12),(0,13),(1,13)):
    large=[n,n,n,n+1] if family==0 else [n,n,2*n,2*n+1]
    A=Q(2*n+1,2);values=[]
    for mask in range(16):
        if mask.bit_count()%2==0:continue
        es=[-1 if mask>>j&1 else 1 for j in range(4)]
        V=even(bg+list(zip(large,es)))
        if family==0:
            T=es[0]*es[1]+es[0]*es[2]+es[1]*es[2]
            expected=2*(A*P1-Q(3,2)*H+T*K)
            r=3
        else:
            expected=2*(A*P1-H+es[0]*es[1]*K)
            r=1
        assert V==expected>0
        if n in (784,575):
            assert A*P1>=(Q(3*(D+1),8)+r)*C0
        values.append(V)
    print("quartet",family,n,sorted(set(values)))

for n in range(5,12):
    for family in (0,1):
        large=[n,n,n,n+1] if family==0 else [n,n,2*n,2*n+1]
        for mask in range(16):
            if mask.bit_count()%2:
                word=bg+[(m,-1 if mask>>j&1 else 1)
                         for j,m in enumerate(large)]
                assert even(word)>0

a=8
P={0:1,a:-3,2*a:5,3*a:-3,4*a:1}
assert sum(P.values())==1
assert sum(v for k,v in P.items() if k<=a)==-2
for j in range(1,41):
    assert sum(v*(k-2*a)**(2*j) for k,v in P.items()) \
           ==2*a**(2*j)*(2**(2*j)-3)>0

def build(word,maxdeg):
    T={(0,0):1}
    for n,e in word:
        N=defaultdict(int)
        for (a,j),v in T.items():
            for k in range(n+1):
                if a+2*k<=maxdeg:N[a+2*k,j]+=v
            if a+n<=maxdeg:
                for k in range(abs(j-n),j+n+1,2):N[a+n,k]+=e*v
        T={ij:v for ij,v in N.items() if v}
    return T

for r,expected in ((9,10658882),(11,407294418)):
    background=[(3,1)]*r+[(3,-1)]*r+[(4,1)]
    T=build(background,2*r)
    V=2*sum(T.get((2*d,0),0) for d in range(r+1))
    assert V==expected
    print("saturation witness",r,V)
print("PASS: analytic constants, 22000 channels, families, obstruction")