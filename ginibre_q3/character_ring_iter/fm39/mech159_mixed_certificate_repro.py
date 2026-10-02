import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from fractions import Fraction as Q
from math import comb
argparse.ArgumentParser(description="Exact FM-MECH159 verifier; no files written.").parse_args()

def cg(a,b): return range(abs(a-b),a+b+1,2)

def fusion(ns):
    d={0:1}
    for n in ns:
        e=defaultdict(int)
        for a,v in d.items():
            for b in cg(a,n): e[b]+=v
        d=e
    return d

@lru_cache(None)
def inv(ns):
    if not ns: return 1
    W=sum(ns); k=len(ns)
    if W%2 or 2*max(ns)>W: return 0
    if k==1: return 0
    ans=0
    for mask in range(1<<k):
        t=W//2-sum(ns[i]+1 for i in range(k) if mask>>i&1)
        if t>=0: ans+=(-1)**mask.bit_count()*comb(t+k-2,k-2)
    return ans

def wt(v):
    v=v[:]; h=1
    while h<len(v):
        for i in range(0,len(v),2*h):
            for j in range(i,i+h):
                a,b=v[j],v[j+h]; v[j],v[j+h]=a+b,a-b
        h*=2
    return v

def phi_table(ns):
    full=(1<<len(ns))-1
    m=[inv(tuple(ns[i] for i in range(len(ns)) if s>>i&1))
       for s in range(full+1)]
    return m,wt([m[s]*m[full^s] for s in range(full+1)])

# Derivatives track the multiplier attached to each y choice.
def row(B,x=0,y=0):
    rem=sum(map(abs,B)); d={(0,0):(1,0,0)}
    for z in sorted(B,key=abs,reverse=True):
        n,e=abs(z),1 if z>0 else -1; rem-=n; out={}
        for (a,b),(v,v1,v2) in d.items():
            for c in cg(a,n):
                if abs(c-x)+abs(b-y)<=rem:
                    q=out.setdefault((c,b),[0,0,0])
                    q[0]+=v; q[1]+=v1; q[2]+=v2
            for c in cg(b,n):
                if abs(a-x)+abs(c-y)<=rem:
                    q=out.setdefault((a,c),[0,0,0])
                    q[0]+=e*v; q[1]+=e*(v1+v); q[2]+=e*(v2+2*v1)
        d={k:v for k,v in out.items() if any(v)}
    return d.get((x,y),(0,0,0))

def quarter(L,a,b):
    C=list(L); C.remove(a); C.remove(b)
    return (1 if b>0 else -1)*row(C,abs(a),abs(b))[0]

def top(B):
    i,j=max(((i,j) for i,j in combinations(range(len(B)),2)
             if (B[i]+B[j])%2==0),
            key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                           max(abs(B[ij[0]]),abs(B[ij[1]]))))
    return tuple(z for k,z in enumerate(B) if k not in (i,j))

witnesses=[
    ((-1,)*10+(-2,)+(-3,)*8+(-4,),6,
     (2135745622,3396151040,104711736020)),
    ((-1,)*7+(-2,)*2+(-3,)*6+(-5,6),6,
     (12867724,3037025,-823007452))]
bounds=[]
for B,p,expected in witnesses:
    N=len(B)+1; g,d1,d2=row(B,p); child=row(top(B),p)[0]
    ds=(N-1)*d1-d2
    assert (g,child,ds)==expected
    L=B+((-1)**sum(z<0 for z in B)*p,)
    counts=Counter(L); direct=0
    for a,b in combinations_with_replacement(counts,2):
        count=comb(counts[a],2) if a==b else counts[a]*counts[b]
        if count: direct+=count*quarter(L,a,b)
    assert ds==direct
    bounds.append(Q(-N*(N-1)*(g-child),ds))
assert bounds==[Q(26468513778,5235586801),Q(1504096947,411503726)]
assert bounds[0]>5>4>bounds[1]

# The lemma covers lengths 3 and 4; length 5 is a conjecture.
lowchecks=0
for k in (3,4,5):
    for ns in combinations_with_replacement(range(1,13),k):
        f=fusion(ns); m=f.get(0,0)
        assert inv(ns)==m
        for r in range(min(ns)//2+1):
            assert f.get(2*r,0)>=(r+1)*m
            lowchecks+=1

def choose(n,k): return comb(n,k) if 0<=k<=n else 0
def negcuts(N,q,s):
    return sum(choose(q,j)*choose(N-q,s-j) for j in range(1,s+1,2))
rows={7:(12,Q(160,7),Q(394,8281)),
      8:(22,Q(127,2),Q(7853,82524)),
      9:(32,Q(4571,36),Q(58909,666468))}
certchecks=0
for N,(u,K,rho) in rows.items():
    values=[sum((Q(negcuts(N,q,s),2 if 2*s==N else 1)*
                 (1+Q(s*(N-s),2*N*(N-1)))
                 for s in range(3,N//2+1)),Q(0))
            for q in range(0,N+1,2)]
    assert max(values)==K
    T=(u//2+1)*(u//2+2)//2
    H=2**(N-3)-N+1
    assert 1-K/T-Q(u+1+H,(u+1)**2)==rho>0
    consecutive=list(range(u,u+N))
    consecutive[-1]+=sum(consecutive)%2
    for ns in ((u,u)+(u+2,)*(N-2),tuple(consecutive)):
        m,v=phi_table(ns); full=len(v)-1
        children=[]
        for a,b in combinations(range(N-1),2):
            if (ns[a]+ns[b])%2: continue
            idx=[i for i in range(N) if i not in (a,b)]
            _,cv=phi_table(tuple(ns[i] for i in idx))
            children.append((idx,cv))
        for mask in range(full+1):
            if mask.bit_count()%2: continue
            if any(ns[i]==ns[i+1] and
                   ((mask>>i&1)!=(mask>>(i+1)&1)) for i in range(N-1)):
                continue
            dif=sum(v[mask]-v[mask^(1<<i)^(1<<j)]
                    for i,j in combinations(range(N),2))
            assert dif%4==0
            ds=dif//4
            for idx,cv in children:
                cm=sum(((mask>>i)&1)<<j for j,i in enumerate(idx[:-1]))
                cm|=(cm.bit_count()%2)<<(N-3)
                drop=Q(v[mask]-cv[cm],2)
                assert cv[cm]>=0
                assert drop>=rho*m[full]
                assert drop+Q(ds,2*N*(N-1))>=rho*m[full]
                certchecks+=1
        if ns[0]==ns[1] and len(set(ns))==2:
            beta=Q(0)
            for s in range(1,1<<(N-1)):
                if (s&3).bit_count()%2 and m[s]*m[full^s]:
                    beta+=Q(1,1+max(ns[i] for i in range(N-1) if s>>i&1))
            assert beta=={7:Q(4,3),8:Q(2),9:Q(16,5)}[N]

signchecks=0
for N,u,negative_count in ((8,20,56),(9,28,112)):
    ns=(u,u)+(u+2,)*(N-2)
    m,v=phi_table(ns); T=(u//2+1)*(u//2+2)//2
    for mask in range(len(v)):
        if mask.bit_count()%2: continue
        if any(ns[i]==ns[i+1] and
               ((mask>>i&1)!=(mask>>(i+1)&1)) for i in range(N-1)):
            continue
        assert Q(v[mask],2)>=Q(T-negative_count,T)*m[-1]
        signchecks+=1
print("PASS", "low-channel controls",lowchecks,"sign checks",signchecks,
      "mixed/removal checks",certchecks,"constant bounds",bounds)
