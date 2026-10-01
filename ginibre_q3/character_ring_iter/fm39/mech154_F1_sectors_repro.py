import argparse
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, product
from math import comb
argparse.ArgumentParser(description='Exact FM-MECH154 checks; no files written.').parse_args()

def A(word, x=0, y=0):
    word=sorted(word,key=abs,reverse=True); rem=sum(map(abs,word))
    if rem<x+y or (rem-x-y)%2: return 0
    out={(0,0):1}
    for z in word:
        n,e=abs(z),1 if z>0 else -1; rem-=n; nxt=defaultdict(int)
        for (a,b),v in out.items():
            for c in range(abs(a-n),a+n+1,2):
                if abs(c-x)+abs(b-y)<=rem: nxt[c,b]+=v
            for c in range(abs(b-n),b+n+1,2):
                if abs(a-x)+abs(c-y)<=rem: nxt[a,c]+=e*v
        out={q:v for q,v in nxt.items() if v}
    return out.get((x,y),0)

def D(L,x,y):
    R=list(L); R.remove(x); R.remove(y)
    return (1 if y>0 else -1)*A(R,abs(x),abs(y))

def flips(L):
    c=Counter(L)
    return {(x,y):D(L,x,y) for x,y in combinations_with_replacement(c,2)
            if x!=y or c[x]>=2}

def mu(ns,p=0):
    r=len(ns); W=sum(ns)
    if p<0 or p>W or (W-p)%2: return 0
    if r<2: return int((ns[0] if r else 0)==p)
    d=(W-p)//2; ans=0
    for mask in range(1<<r):
        t=d-sum(ns[i]+1 for i in range(r) if mask>>i&1)
        if t>=0: ans+=(-1)**mask.bit_count()*comb(t+r-2,r-2)
    return ans

B=(1,1,-2)+(3,)*4+(-4,)*3+(5,)*4; L=B+(8,)
v=flips(L); R=list(B); R.remove(5); R.remove(5)
assert len(v)==19 and (min(v.values()),max(v.values()))==(-4567048,-20406)
assert (A(L),A(B,8),A(R,8))==(8150742,4075371,160741)
assert max(flips(tuple(-z if abs(z)%2 else z for z in L)).values())<0

L=(1,)*2+(2,)*5+(-3,)*4+(5,)*8+(-6,)*2+(8,)
assert A(L)==357306925012
assert [D(L,1,z) for z in (1,2,-3,5,-6,8)]==[
    -1117195873496,-1690562261,-39718883522,
    -29334565523,-35819446453,-95417686281]
assert D(L,2,-3)==28541128639
L=(1,)*2+(2,)*5+(-3,)*4+(5,)*4+(8,)
assert [D(L,z,z) for z in (1,2,-3,5)]==[
    -15335632,-72818,-1543816,-1003404]
assert D(L,1,2)==1199234
assert (D((1,1,-5,-5,-5,7,-8),1,1),
        D((1,1,-5,-5,-5,7,-8),1,-8))==(-4,28)

# Small exact controls for the two uniform CG arguments.
for ns in combinations_with_replacement(range(2,8),4):
    for signs in product((-1,1),repeat=len(set(ns))):
        e=dict(zip(sorted(set(ns)),signs))
        L=(1,1)+tuple(e[n]*n for n in ns)
        if sum(z<0 for z in L)%2==0:
            assert max(flips(L).values())>=0
for a,b,c in combinations_with_replacement(range(3,12,2),3):
    h=mu((a,b,c),1)
    if not h: continue
    for P in range(max(5,c),14,2):
        ns=(a,b,c,P); T=0
        for i,j in combinations(range(4),2):
            if ns[i]+ns[j]>=P+1:
                T+=mu((1,1)+tuple(ns[k] for k in range(4) if k not in (i,j)))
        assert mu(ns,P+1)>=ns.count(P)*h+T
for n in range(1,7):
    for p in range(1,16):
        rhs=6*mu((n,)*3,p)-6*(n%2==0)*mu((n,)*2,p)+2*(n+1)*mu((n,),p)
        assert D((-n,)*5+(-p,),-n,-n)==rhs>=0

for r,C in ((2,(-2,-3,5)),(3,(-2,-4,5)),(4,(-3,-3,4)),
            (2,(-5,-5,-5,7,-8)),(5,(2,4,6))):
    L=(1,)*r+C; R=(1,)*(r-2)+C
    rhs=4*(r-2)*A(R)
    for z,count in Counter(C).items():
        rest=list(R); rest.remove(z)
        rhs+=4*count*(abs(z)+1)*(1 if z>0 else -1)*A(rest,1,abs(z)-1)
    assert (sum(map(abs,L))+4)*D(L,1,1)==rhs
print('PASS: counterexamples, partial-sector controls, shifted-label recurrence')
