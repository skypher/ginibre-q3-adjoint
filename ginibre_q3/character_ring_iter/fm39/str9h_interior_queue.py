import sys
if any(x in ('-h','--help') for x in sys.argv[1:]):
    print('FM-STR9h exact interior-cut verifier: labels <= 4, length <= 7, every distinguished pair; also all pairs in three FM-SEC181 witnesses. Python integers.')
    raise SystemExit(0)
from collections import defaultdict, Counter
from datetime import datetime, timezone
from functools import lru_cache
from itertools import combinations_with_replacement
def log(*x): print(datetime.now(timezone.utc).isoformat(timespec='seconds'),*x,flush=True)
def cg(a,b): return range(abs(a-b),a+b+1,2)
@lru_cache(None)
def char(w):
    d={(0,0):1}
    for z in sorted(w,key=lambda q:(-abs(q),q)):
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (r,s),v in d.items():
            for t in cg(r,n): q[t,s]+=v
            for t in cg(s,n): q[r,t]+=e*v
        d={k:v for k,v in q.items() if v}
    return d
def add(P,Q):
    d=defaultdict(int)
    for k,v in P.items():d[k]+=v
    for k,v in Q.items():d[k]+=v
    return {k:v for k,v in d.items() if v}
def mul(d,n,e):
    q=defaultdict(int)
    for (r,s),v in d.items():
        for t in cg(r,n):q[t,s]+=v
        for t in cg(s,n):q[r,t]+=e*v
    return {k:v for k,v in q.items() if v}
def xx(d,a,b,ey):
    q=defaultdict(int)
    for (r,s),v in d.items():
        for c in cg(a,b):
            for t in cg(r,c):q[t,s]+=v
            if ey:
                for t in cg(s,c):q[r,t]+=ey*v
    return {k:v for k,v in q.items() if v}
def xy(d,a,b,ea,eb):
    q=defaultdict(int)
    for (r,s),v in d.items():
        for t in cg(r,a):
            for u in cg(s,b):q[t,u]+=eb*v
        for t in cg(r,b):
            for u in cg(s,a):q[t,u]+=ea*v
    return {k:v for k,v in q.items() if v}
def cut(C):
    A=[];B=[];wa=wb=0
    for z in sorted(C,key=lambda q:(-abs(q),q)):
        if wa<=wb:A.append(z);wa+=abs(z)
        else:B.append(z);wb+=abs(z)
    if wa>wb:A,B=B,A
    return tuple(A),tuple(B)
def fifo(src,dst):
    q=[]
    for h in range(len(src)):
        if src[h]:q.append([h,src[h]])
        need=dst[h]
        while need:
            if not q:return False
            k,v=q[0];assert k<=h
            take=min(need,v);need-=take;v-=take
            if v:q[0][1]=v
            else:q.pop(0)
    return True
def profile(w,i,j):
    u,v=w[i],w[j];a,b=abs(u),abs(v);eu=1 if u>0 else -1;ev=1 if v>0 else -1
    C=tuple(z for k,z in enumerate(w) if k not in (i,j));A,B=cut(C)
    fA,fB=char(A),char(B);tA=char(tuple(abs(z) for z in A));tB=char(tuple(abs(z) for z in B))
    ps=xx(fB,a,b,eu*ev);ms=xy(fB,a,b,eu,ev)
    pt=xx(tB,a,b,1);mt=xy(tB,a,b,1,1)
    bp=mul(mul(fB,a,eu),b,ev);assert bp==add(ps,ms)
    P=defaultdict(int);M=defaultdict(int);PT=defaultdict(int);MT=defaultdict(int)
    for k,x in fA.items():
        h=sum(k);P[h]+=x*ps.get(k,0);M[h]+=x*ms.get(k,0)
    for k,x in tA.items():
        h=sum(k);PT[h]+=x*pt.get(k,0);MT[h]+=x*mt.get(k,0)
    H=max(set(P)|set(M)|set(PT)|set(MT),default=0)
    pos=[];neg=[];supply=[];demand=[];pi=bud=mp=mb=0;badT=None
    for h in range(H+1):
        p,m=P[h],M[h];tp,tm=PT[h],MT[h]
        assert tp>=abs(p) and tm>=abs(m) and (tp+p)%2==0 and (tm+m)%2==0
        z=tp+tm;x=p+m
        assert z>=abs(x) and (z+x)%2==0
        pos.append((z+x)//2);neg.append((z-x)//2)
        supply.append(max(p,0));demand.append(max(-p,0)+max(-m,0))
        pi+=x;bud+=p+min(m,0);mp=min(mp,pi);mb=min(mb,bud)
        if bud<0 and badT is None:badT=h
        assert pi>=0,(w,i,j,h,pi)
    assert fifo(pos,neg)
    bmatch=mb>=0 and fifo(supply,demand)
    assert sum(P.values())+sum(M.values())==char(A+B+(u,v)).get((0,0),0)
    return mp,mb,badT,bmatch
def pairfree(w):
    d=defaultdict(set)
    for z in w:d[abs(z)].add(1 if z>0 else -1)
    return all(len(v)==1 for v in d.values())
def main():
    log('begin interior-cut census; signed multisets; even minus count')
    types=tuple(e*n for n in range(1,5) for e in (-1,1))
    nl=np=bf=pf=0;minp=minb=0
    for L in range(1,8):
        cl=cp=0
        for w in combinations_with_replacement(types,L):
            if sum(z<0 for z in w)%2:continue
            nl+=1;cl+=1
            for i in range(L):
                for j in range(i+1,L):
                    p,b,t,ok=profile(w,i,j);np+=1;cp+=1
                    pf+=p<0;bf+=b<0;minp=min(minp,p);minb=min(minb,b)
                    assert p>=0 and b>=0 and ok
        log('length',L,'lists',cl,'interior pairs',cp)
    assert (nl,np,pf,bf)==(3234,54236,0,0)
    log('small-box PASS lists',nl,'pairs',np,'plain failures',pf,'budget failures',bf,
        'min plain prefix',minp,'min budget',minb)
    ws=((1,-2,-4,-4,-4,1)+(3,)*6+(5,)*4+(6,),
        (1,-2,-4,-4,-4,1)+(3,)*5+(5,)*5+(6,),
        (1,-2,-4,-4,-4,1,1,2)+(3,)*5+(4,4)+(5,5)+(8,))
    expected=(136,136,153)
    for ix,w in enumerate(ws):
        n=len(w)*(len(w)-1)//2;assert n==expected[ix]
        fail=Counter();minp=minb=0;plainbad=0
        for i in range(len(w)):
            for j in range(i+1,len(w)):
                p,b,t,ok=profile(w,i,j);minp=min(minp,p);minb=min(minb,b);plainbad+=p<0
                if b<0:fail[tuple(sorted((w[i],w[j])))]+=1
                else:assert ok
        log('witness',ix+1,'W',sum(map(abs,w)),'factors',len(w),'pairfree',pairfree(w),
            'pairs',n,'plain failures',plainbad,'budget failures',sum(fail.values()),
            'failing pair values',dict(fail),'min plain',minp,'min budget',minb)
        assert plainbad==0
        if ix==0:assert set(fail)=={(-4,-2)}
        else:assert not fail
    A0,B0=cut((-1,-1,-2))
    fa,fb=char(A0),char(B0);assert A0==(-2,) and B0==(-1,-1)
    assert fa[(0,2)]*fb[(0,2)]==-1 and fa[(2,0)]*fb[(2,0)]==1
    cross=profile((-1,-1,-1,-2,3),2,4)
    assert cross[0]>=0 and cross[1]>=0
    log('pair-free cross-channel example','word=(-1,-1,-1,-2,3)','pair=(-1,3)',
        'cut A',A0,'B',B0,'negative channel (0,2)=-1','positive channel (2,0)=+1',
        'height=2')
    A=(1,-2,-4,-4,-4);B=(1,)+(3,)*6+(5,)*4+(6,)
    pi=sum(x*char(B).get(k,0) for k,x in char(A).items() if sum(k)<=11)
    assert pi==-24695910
    log('segregated Pi_11',pi,'interior cuts remain nonnegative')
    log('PASS FM-STR9h')
main()
