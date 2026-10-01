import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, prod
from fractions import Fraction
from random import Random

argparse.ArgumentParser(description="FM-MECH142 exact verifier; no file writes.").parse_args()

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def fusion(ns):
    A={0:1}
    for n in ns:
        T=defaultdict(int)
        for a,v in A.items():
            for c in cg(a,n): T[c]+=v
        A=dict(T)
    return A

def character(word):
    A={(0,0):1}
    for z in word:
        n=abs(z); e=1 if z>0 else -1
        T=defaultdict(int)
        for (a,b),v in A.items():
            for c in cg(a,n): T[c,b]+=v
            for c in cg(b,n): T[a,c]+=e*v
        A={k:v for k,v in T.items() if v}
    return A

def phi(word):
    return character(word).get((0,0),0)

PARTS=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))

def descent(ns,D,eta):
    row=fusion(ns); N=[0]*(D+1); Z=[0]*(D+1); children=[]
    for I,J in PARTS:
        ex=prod(eta[i] for i in I); ey=prod(eta[i] for i in J)
        for a in fusion(tuple(ns[i] for i in I)):
            for b in fusion(tuple(ns[i] for i in J)):
                if a+b>D: continue
                if a==0 or b==0:
                    assert (a!=0 or ex==1) and (b!=0 or ey==1)
                    Z[a+b]+=1
                else:
                    children.append((ex*a,ey*b))
                    for c in cg(a,b): N[c]+=1
    assert all(row.get(c,0)>=N[c] for c in range(D+1))
    return [row.get(c,0)+Z[c]-N[c] for c in range(D+1)],children

rng=Random(142); checked=0
for _ in range(120):
    C=tuple(rng.choice((-1,1))*rng.randint(1,3)
            for _ in range(rng.randrange(5)))
    D=sum(map(abs,C)); t=[0]+sorted(rng.randrange(5) for _ in range(3))
    if (D+sum(t))%2: t[-1]+=1
    q=D+max(t)+1; eta=[rng.choice((-1,1)) for _ in range(4)]
    for i in range(4):
        for j in range(i):
            if t[i]==t[j]: eta[i]=eta[j]
    sigma=(-1)**sum(z<0 for z in C)
    if prod(eta)!=sigma: continue
    ns=tuple(q+x for x in t); H,children=descent(ns,D,eta)
    rhs=sum(phi(C+E) for E in children)+2*H[0]*phi(C)
    rhs+=sum(H[c]*phi(C+(sigma*c,)) for c in range(1,D+1))
    assert phi(C+tuple(e*n for e,n in zip(eta,ns)))==rhs
    checked+=1
print("positive descent identities:",checked)

def mono_phi(ns):
    A={(0,0):1}
    for n in ns:
        T=defaultdict(int)
        for j in range(n//2+1):
            k=n-2*j; v=(-1)**j*comb(n-j,j)
            for (a,b),w in A.items():
                T[a+k,b]+=v*w
                T[a,b+k]-=v*w
        A={k:v for k,v in T.items() if v}
    cat=lambda n:comb(2*n,n)//(n+1)
    return sum(v*cat(a//2)*cat(b//2) for (a,b),v in A.items()
               if a%2==b%2==0)

tests=[
 ((1,)*12+(6,6,6,8),(1,)*12+(5,5,5,7),11273148,12382358),
 ((1,)*26+(2,3,3,4)+(7,)*4,
  (1,)*26+(2,3,3,4)+(5,)*4,35495241211284322,37223085948167400)
]
for A,B,x,y in tests:
    assert phi(tuple(-n for n in A))==mono_phi(A)==x
    assert phi(tuple(-n for n in B))==mono_phi(B)==y
    assert x<y
    print("shift obstruction:",x,y)

D=18
dual={
 (0,0):2895,(0,2):4788,(0,4):406,(0,8):256,
 (0,14):524,(0,16):104,(0,18):80,
 (1,1):2987,(1,7):64,(1,17):92,
 (2,2):2078,(2,8):244,(2,10):256,
 (3,9):76,(3,11):-152,(3,13):-628,
 (4,6):-2230,(4,10):640,(4,14):-92,
 (5,5):-523,(5,11):-76,(5,13):-964,
 (6,6):-229,(6,8):510,(6,10):1290,(6,12):-964,
 (7,7):5614,(7,9):352,(8,8):2953,(8,10):5962,(9,9):5968
}
coords=[(a,b) for a in range(D+1) for b in range(a,D+1-a)
        if (a+b)%2==0]
def dual_value(A):
    return sum(v*A.get(k,0) for k,v in dual.items())

def m3(a,b,n,c):
    if (a+b-n-c)%2: return 0
    lo=max(abs(a-b),abs(n-c)); hi=min(a+b,n+c)
    return max(0,(hi-lo)//2+1)

def short_poly(ns,es):
    T=defaultdict(int); L=len(ns)
    if L==0: T[0,0]=1
    elif L==1:
        if ns[0]<=D: T[0,ns[0]]=1
    elif L==2:
        for c in cg(*ns):
            if c<=D: T[0,c]+=2 if c==0 else 1
        if sum(ns)<=D:
            T[min(ns),max(ns)]+=es[0]*(2 if ns[0]==ns[1] else 1)
    else:
        for c in range(sum(ns)%2,D+1,2):
            T[0,c]+=m3(*ns,c)*(2 if c==0 else 1)
        for i,z in enumerate(ns):
            if z>D: continue
            for c in cg(*(ns[:i]+ns[i+1:])):
                if c+z>D: break
                T[min(c,z),max(c,z)]+=es[i]*(2 if c==z else 1)
    return T

dual_checks=0
for L in range(4):
    for ns in combinations_with_replacement(range(1,2*D+2),L):
        if sum(ns)%2: continue
        for es in product((-1,1),repeat=L):
            if prod(es)!=1: continue
            A=short_poly(ns,es)
            assert dual_value(A)>=0,(ns,es)
            if not ns or max(ns)<=4:
                B=character(tuple(n*e for n,e in zip(ns,es)))
                assert all(A.get(k,0)==B.get(k,0) for k in coords)
            dual_checks+=1
assert all(dual.get((0,c),0)>=0 for c in range(0,D+1,2))
assert dual_value(character((-7,-8,-9,-10)))==-62922
print("short-word dual checks:",dual_checks,"target:",-62922)

U=tuple(range(8))
@lru_cache(None)
def invariant(T,r):
    L=len(T); weight=L*r+sum(T)
    if weight%2: return 0
    ans=0
    for k in range(L+1):
        for V in combinations(T,k):
            h=weight//2-k*r-sum(v+1 for v in V)
            if h>=0: ans+=(-1)**k*comb(h+L-2,L-2)
    return ans

def consecutive_g(r):
    out=invariant(U,r)
    for T in combinations(U,3):
        V=tuple(x for x in U if x not in T)
        out-=invariant(T,r)*invariant(V,r)
    for T in combinations(U,4):
        if 0 not in T: continue
        V=tuple(x for x in U if x not in T)
        out+=invariant(T,r)*invariant(V,r)
    return out

# Every inclusion-exclusion branch has a fixed sign for s >= 7.
for parity in (0,1):
    for L in (3,4,5,8):
        for T in combinations(U,L):
            if (L*parity+sum(T))%2: continue
            for k in range(L+1):
                for V in combinations(T,k):
                    a=L-2*k
                    b=(L*parity+sum(T))//2-k*parity-sum(v+1 for v in V)
                    if a>0: assert 7*a+b>=0
                    elif a<0: assert 7*a+b<0

def claimed(s,parity):
    if parity==0:
        return Fraction(32*s**5+360*s**4+1376*s**3+
                        2016*s*s+764*s+291,3)
    return Fraction(32*s**5+440*s**4+2176*s**3+
                    4660*s*s+3984*s+1323,3)

# Degree <= 6, so seven exact values prove each stable polynomial identity.
for parity in (0,1):
    for s in range(7,14):
        assert consecutive_g(2*s+parity)==claimed(s,parity)>0
corners=[consecutive_g(r) for r in range(1,14)]
assert min(corners)>0
for r in (1,3,14,15):
    assert phi(tuple(-n for n in range(r,r+8)))==2*consecutive_g(r)
print("consecutive corners r=1..13:",corners)
print("FM-MECH142 EXACT VERIFIER PASS")
