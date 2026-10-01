from collections import defaultdict
from functools import lru_cache

@lru_cache(None)
def profile(B):
    d={(0,0):1}
    for signed_n in B:
        n=abs(signed_n); eps=1 if signed_n>0 else -1
        nd=defaultdict(int)
        for (a,b),v in d.items():
            for u in range(abs(a-n),a+n+1,2): nd[u,b]+=v
            for u in range(abs(b-n),b+n+1,2): nd[a,u]+=eps*v
        d={key:value for key,value in nd.items() if value}
    return tuple(d.items())

def canon(B): return tuple(sorted(B,key=lambda z:(abs(z),z)))
def gp(B,p): return dict(profile(canon(B))).get((p,0),0)
def fuses(a,n,p): return abs(a-n)<=p<=a+n and (a+n-p)%2==0

def class_pairs(B,same_parity=False):
    labels=sorted(set(B),key=lambda z:(abs(z),z)); out=[]
    for i,x in enumerate(labels):
        n=abs(x); ex=1 if x>0 else -1
        if B.count(x)>=2: out.append((x,x))
        for y in labels[i+1:]:
            m=abs(y); ey=1 if y>0 else -1
            if ex==ey and (not same_parity or (n-m)%2==0):
                out.append((x,y))
    return out

def rest_for(B,x,y):
    C=list(B); C.remove(x); C.remove(y)
    return canon(C)

def cross(C,p,n,m):
    c=dict(profile(canon(C)))
    return (sum(v for (a,b),v in c.items() if b==m and fuses(a,n,p))+
            sum(v for (a,b),v in c.items() if b==n and fuses(a,m,p)))

def channel_values(B,p,x,y):
    C=rest_for(B,x,y); n=abs(x); m=abs(y); out=[]
    for k in range(abs(n-m),n+m+1,2):
        value=2*gp(C,p) if k==0 else gp(C+(k,),p)
        out.append((k,value))
    return out

def positive_moves(B,p):
    out=[]
    for x,y in class_pairs(B):
        for k,value in channel_values(B,p,x,y):
            if k>0: out.append((value,x,y,k))
    return out

def residual(B,p):
    W=sum(map(abs,B)); M=max(map(abs,B)); delta=(W-p)//2
    return (p>=6 and p>=M and (W-p)%2==0 and delta>=8 and M<=delta
            and sum(abs(z)>=3 for z in B)>=2)

def all_backgrounds(maxW):
    def rec(start,rem,B):
        if rem==0: yield B; return
        for n in range(start,rem+1):
            for c in range(1,rem//n+1):
                for eps in (-1,1):
                    yield from rec(n+1,rem-c*n,B+(eps*n,)*c)
    for W in range(1,maxW+1): yield from rec(1,W,())

def screen(maxW):
    cases=moves=0; best=None
    for B in all_backgrounds(maxW):
        if len(B)<2: continue
        W=sum(map(abs,B)); M=max(map(abs,B))
        for p in range(max(6,M),W+1):
            if not residual(B,p): continue
            cases+=1; parent=gp(B,p); vals=positive_moves(B,p); moves+=len(vals)
            if not vals or min(v[0] for v in vals)>parent:
                return cases,moves,(B,p,parent,vals),best
            child,x,y,k=min(vals); drop=parent-child
            if best is None or drop*best[1] < best[0]*parent:
                best=(drop,parent,B,p,x,y,k,child)
    return cases,moves,None,best

# Fixed smallest-pair selector witness.
B85=canon((-11,)*3+(-2,)*2+(6,)*8); p85=11
pairs85=class_pairs(B85,same_parity=True)
distinct=[(abs(x),abs(y),x,y) for x,y in pairs85 if x!=y]
if distinct: _,_,x85,y85=min(distinct)
else: _,_,x85,y85=min((abs(x),abs(y),x,y) for x,y in pairs85 if x==y)
k85=abs(abs(x85)-abs(y85)) if x85!=y85 else 2
fixed_child=gp(rest_for(B85,x85,y85)+(k85,),p85)
assert residual(B85,p85)
assert (gp(B85,p85),x85,y85,k85,fixed_child)==(
    444227708,-2,-2,2,549660463)
best85=min(positive_moves(B85,p85))
assert best85[0]<gp(B85,p85)
print('fixed-selector witness',B85,'p',p85,'parent',gp(B85,p85),
      'chosen',(x85,y85,k85),'child',fixed_child)
print('same-sign fusion menu on that witness: best (child,n,m,k)',best85)

# All-negative correction-sign witness; check the fusion identity exactly.
B22=canon((-1,-2,-3,-3,-4,-4,-5)); p22=6; parent22=gp(B22,p22)
cross_rows=[]
for x,y in class_pairs(B22):
    eps=1 if x>0 else -1
    X=cross(rest_for(B22,x,y),p22,abs(x),abs(y))
    channels=channel_values(B22,p22,x,y)
    assert parent22==sum(value for k,value in channels)+eps*X
    pos=[value for k,value in channels if k>0]
    cross_rows.append((eps*X,x,y,min(pos),max(pos)))
best22=min(positive_moves(B22,p22))
assert residual(B22,p22) and parent22==273
assert max(row[0] for row in cross_rows)==-6
assert best22==(27,-4,-5,1)
print('all-negative-cross witness',B22,'p',p22,'parent',parent22)
print('(eps*X,n,m,min positive-channel child,max positive-channel child)',
      cross_rows)
print('best direct fusion child (value,n,m,k)',best22)

# Exact full-menu screen on the residual through total weight 22.
result=screen(22)
assert result==(4724,68252,None,
    (34002,51874,(1,1,1,1,1,1,1,1,-2,-2,-2,-2,3,3),6,1,3,2,17872))
print('W<=22 residual screen: cases, positive-k class-pair channels,'
      ' failure, least relative drop =',result)
