import argparse
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from math import comb
from random import Random

argparse.ArgumentParser(
    description="FM-MECH112: affine quartet descent and channel obstruction."
).parse_args()

def step(A,n,e):
    B=defaultdict(int)
    for (a,b),v in A.items():
        for j in range(abs(a-n),a+n+1,2): B[j,b]+=v
        for j in range(abs(b-n),b+n+1,2): B[a,j]+=e*v
    return {k:v for k,v in B.items() if v}

def character(word):
    A={(0,0):1}
    for z in word: A=step(A,abs(z),1 if z>0 else -1)
    return A

@lru_cache(None)
def fusion(ns):
    A={0:1}
    for n in ns:
        B=defaultdict(int)
        for a,v in A.items():
            for b in range(abs(a-n),a+n+1,2): B[b]+=v
        A=dict(B)
    return A

def beta(a,t):
    assert (a-sum(t))%2==0
    lo1,lo2=abs(t[0]-t[1]),abs(t[2]-t[3])
    hi1,hi2=t[0]+t[1],t[2]+t[3]
    ans=0
    for d in range(-a,a+1,2):
        lo=max(lo2,lo1-d,(a-d)//2)
        lo+=(hi2-lo)%2
        hi=min(hi2,hi1-d)
        assert (hi-lo)%2==0
        ans+=1+(hi-lo)//2
    return ans

def quartet(B,t,e):
    D=sum(map(abs,B)); M=max(t); q0=D+M+1
    assert (D+sum(t))%2==0 and len(t)==len(e)==4
    eta=(-1)**sum(z<0 for z in B)
    assert e[0]*e[1]*e[2]*e[3]==eta
    A=character(B)
    G=sum((a+1)*v for (a,b),v in A.items() if b==0)
    J=sum(beta(a,t)*v for (a,b),v in A.items() if b==0)
    def low(a,i,j):
        return a>=abs(t[i]-t[j]) and (a-t[i]-t[j])%2==0
    K=0
    for S in combinations(range(4),2):
        C=tuple(i for i in range(4) if i not in S)
        K+=e[S[0]]*e[S[1]]*sum(
            v for (a,b),v in A.items() if low(a,*C) and low(b,*S))
    assert K%2==0 and G>0
    C=J+K//2
    Q=max(q0,(-C+G-1)//G)
    return A,G,C,q0,Q

def assigned_integral(A,t,e,q):
    ans=0
    for mask in range(16):
        X=tuple(q+t[i] for i in range(4) if not(mask>>i&1))
        Y=tuple(q+t[i] for i in range(4) if mask>>i&1)
        sign=1
        for i in range(4):
            if mask>>i&1: sign*=e[i]
        fx,fy=fusion(X),fusion(Y)
        ans+=sign*sum(v*fx.get(a,0)*fy.get(b,0)
                      for (a,b),v in A.items())
    return ans

rng=Random(112)
checks=0
for _ in range(100):
    B=tuple(rng.choice((-1,1))*rng.randint(1,5)
            for _ in range(rng.randrange(7)))
    D=sum(map(abs,B))
    t=[0]+sorted(rng.randrange(6) for _ in range(3))
    if (D+sum(t))%2: t[-1]+=1
    t=tuple(t)
    e=[rng.choice((-1,1)) for _ in range(3)]
    e.append((-1)**sum(z<0 for z in B)*e[0]*e[1]*e[2])
    A,G,C,q0,Q=quartet(B,t,e)
    for q in (q0,q0+1):
        row=fusion(tuple(q+s for s in t))
        for a in range(D+1):
            predicted=((a+1)*q+beta(a,t)) if (a-sum(t))%2==0 else 0
            assert row.get(a,0)==predicted
            checks+=1
        assert assigned_integral(A,t,e,q)==2*(q*G+C)
    assert Q*G+C>=0
print("100 SIGNED QUARTETS PASS;",checks,"LOW-SPIN IDENTITIES",flush=True)

examples=[
    ((-1,-1,-1,-2,-3,-4),(0,0,0,0),(-1,)*4,(670,-518,13)),
    ((-1,-2,-3,-4),(0,0,0,2),(-1,)*4,(114,-126,13)),
    ((-1,-2,-3),(0,0,0,2),(-1,-1,-1,1),(23,0,9)),
    ((-1,)*5+(-3,)*6+(-5,),(0,0,0,2),(-1,)*4,
     (23865996,73098810,31))
]
for B,t,e,expected in examples:
    A,G,C,q0,Q=quartet(B,t,e)
    assert (G,C,Q)==expected
    for q in (Q,Q+1,Q+5):
        assert assigned_integral(A,t,e,q)==2*(q*G+C)>0
    print("AFFINE",B,t,"g(q) =",G,"* q +",C,"for q >=",Q,flush=True)

# Independent monomial/Catalan integration.
def poly_mul(A,B):
    C=defaultdict(int)
    for (i,j),v in A.items():
        for (k,l),w in B.items(): C[i+k,j+l]+=v*w
    return {k:v for k,v in C.items() if v}

def un(n):
    a,b={0:1},{1:1}
    if n==0:return a
    for _ in range(1,n):
        c=defaultdict(int)
        for k,v in b.items(): c[k+1]+=v
        for k,v in a.items(): c[k]-=v
        a,b=b,dict(c)
    return b

def catalan_moment(word):
    A={(0,0):1}
    for z in word:
        B=defaultdict(int);e=1 if z>0 else -1
        for k,v in un(abs(z)).items(): B[k,0]+=v;B[0,k]+=e*v
        A=poly_mul(A,B)
    def cat(n):return comb(2*n,n)//(n+1)
    return sum(v*cat(a//2)*cat(b//2) for (a,b),v in A.items()
               if a%2==b%2==0)

base=(-1,-1,-1,-2,-3,-4)
for q in (13,14):
    assert catalan_moment(base+(-q,)*4)==2*(670*q-518)
print("INDEPENDENT CATALAN CHECKS PASS",flush=True)

# Exact invariant fusion path in polynomial coordinates.
@lru_cache(None)
def embed(a,b,j,t):
    k=(a+b-j)//2
    A={(r,k-r):F((-1)**r*comb(k,r)) for r in range(k+1)}
    for l in range(t):
        B=defaultdict(F)
        for (u,v),c in A.items():
            if u<a:B[u+1,v]+=c*(a-u)/(j-l)
            if v<b:B[u,v+1]+=c*(b-v)/(j-l)
        A=dict(B)
    return A

@lru_cache(None)
def amplitudes(ns,path,t):
    if len(ns)==1:return {(t,):F(1)}
    out={}
    for (u,v),c in embed(path[-2],ns[-1],path[-1],t).items():
        for state,d in amplitudes(ns[:-1],path[:-1],u).items():
            out[state+(v,)]=c*d
    return out

def laurent_mul(A,B):
    C=defaultdict(F)
    for i,a in A.items():
        for j,b in B.items():C[i+j]+=a*b
    return {j:v for j,v in C.items() if v}

ns=(1,1,1,2,3,4); path=(1,2,3,1,4,0)
V=amplitudes(ns,path,0)
raised=defaultdict(F)
P=defaultdict(F); norm=F(0)
for state,c in V.items():
    assert sum(n-2*t for n,t in zip(ns,state))==0
    for i,t in enumerate(state):
        if t:
            s=list(state);s[i]-=1;raised[tuple(s)]+=t*c
    wt=1;pol={0:F(1)}
    for n,t in zip(ns,state):
        wt*=comb(n,t);m=n-2*t
        term=defaultdict(F);term[m]+=1;term[-m]-=1
        pol=laurent_mul(pol,term)
    w=c*c/wt;norm+=w
    for j,v in pol.items():P[j]+=w*v
assert all(v==0 for v in raised.values())
P={j:v/norm for j,v in P.items() if v}
assert P[0]==-F(8,5) and P[4]==F(2,5)
for q in (8,14,20):
    S={0:F(2)}
    for m in range(-q,q+1,2):S[2*m]=S.get(2*m,0)-F(2,q+1)
    T=laurent_mul(P,laurent_mul(S,S))
    channel=T.get(0,0)-T.get(4,0)
    assert channel==-8-F(32,5*(q+1)**2)<0
print("NEGATIVE COMPLETE FUSION CHANNEL: -8-32/[5(q+1)^2]",flush=True)
print("FM-MECH112 EXACT VERIFIER PASS",flush=True)