import argparse, collections, datetime, functools, heapq, itertools, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

ap = argparse.ArgumentParser(description="Exact F1 complete-block obstruction; memory only.")
ap.add_argument("--dump", action="store_true", help="Print the integer scalar direction.")
args = ap.parse_args()
def log(*a):
    print(datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), *a, flush=True)

@functools.lru_cache(None)
def character(Q):
    if not Q: return {0: 1}
    r = collections.defaultdict(int)
    for a, v in character(Q[:-1]).items():
        for b in range(abs(a-Q[-1]), a+Q[-1]+1, 2): r[b] += v
    return dict(r)
def m(Q): return character(Q).get(0, 0)
def cuts(Q):
    cc = tuple(collections.Counter(Q).items())
    for z in itertools.product(*(range(n+1) for a,n in cc)):
        yield (tuple(a for (a,n),j in zip(cc,z) for _ in range(j)),
               tuple(a for (a,n),j in zip(cc,z) for _ in range(n-j)))
def fusions(Q):
    cc = collections.Counter(Q)
    for a in cc:
        for b in cc:
            if a>b or (a==b and cc[a]<2): continue
            i=Q.index(a); j=Q.index(b,i+1) if a==b else Q.index(b)
            C=Q[:i]+Q[i+1:j]+Q[j+1:]
            yield i,j,tuple(tuple(sorted(C+((c,) if c else ())))
                            for c in range(abs(a-b),a+b+1,2))
@functools.lru_cache(None)
def K(n):
    return np.array([[sum((-1)**h*math.comb(t,h)*math.comb(n-t,j-h)
                         for h in range(max(0,j-n+t),min(j,t)+1))
                      for j in range(n+1)] for t in range(n+1)],dtype=object)
def negative_counts(Q,s):
    return tuple(sum((s>>i)&1 for i,b in enumerate(Q) if b==a)
                 for a in collections.Counter(Q))
def maskof(Q,z):
    seen=collections.Counter(); out=0
    for i,a in enumerate(Q):
        if seen[a]<z.get(a,0): out|=1<<i; seen[a]+=1
    return out
def tensor(Q,v=None):
    cc=collections.Counter(Q); h=[]
    for A,B in cuts(Q):
        h.append(m(A)*m(B) if v is None else v.get(A,0)*m(B)+m(A)*v.get(B,0))
    h=np.array(h,dtype=object).reshape(tuple(n+1 for n in cc.values()))
    for ax,n in enumerate(cc.values()):
        h=np.moveaxis(np.tensordot(K(n),h,axes=(1,ax)),0,ax)
    return h
def child(L,s):
    i,j=max(((i,j) for i in range(len(L)-1) for j in range(i+1,len(L)-1)
             if (L[i]-L[j])%2==0),
            key=lambda ij:(L[ij[0]]+L[ij[1]],max(L[ij[0]],L[ij[1]])))
    keep=[h for h in range(len(L)) if h not in (i,j)]
    C=tuple(L[h] for h in keep)
    t=sum(((s>>h)&1)<<q for q,h in enumerate(keep))
    if ((s>>i)^(s>>j))&1: t^=1<<(len(C)-1)
    return C,t
def actual(L,s):
    h=tensor(L); C,t=child(L,s); hc=tensor(C)
    p=int(h[negative_counts(L,s)]); q=int(hc[negative_counts(C,t)])
    ds=[F(p-int(h[negative_counts(L,s^(1<<i)^(1<<j))]),4)
        for i in range(len(L)) for j in range(i+1,len(L))]
    assert (p-q)%2==0
    return p//2,q//2,(p-q)//2,max(ds)
def sparse(rows,n):
    ii=[]; jj=[]; vv=[]
    for i,r in enumerate(rows):
        for j,v in r.items():
            if v:
                assert -(1<<63)<v<(1<<63)
                ii.append(i); jj.append(j); vv.append(v)
    return coo_matrix((vv,(ii,jj)),shape=(len(rows),n),dtype=np.int64).tocsr()
def mv(A,x):
    return [sum(int(v)*x[j] for j,v in
                zip(A.indices[A.indptr[i]:A.indptr[i+1]],
                    A.data[A.indptr[i]:A.indptr[i+1]])) for i in range(A.shape[0])]
def exact(A,b,approx):
    A=A.tocsr(); eq=[]; bs=[]; inc=[set() for _ in approx]
    for i in range(A.shape[0]):
        r={int(j):F(int(v)) for j,v in
           zip(A.indices[A.indptr[i]:A.indptr[i+1]],A.data[A.indptr[i]:A.indptr[i+1]]) if v}
        assert int(b[i])==b[i]
        z=F(int(b[i]))
        if not r:
            assert not z
            continue
        h=len(eq); eq.append(r); bs.append(z)
        for j in r: inc[j].add(h)
    hp=[(len(r),i) for i,r in enumerate(eq)]; heapq.heapify(hp); saved=[]
    while hp:
        w,i=heapq.heappop(hp); r=eq[i]
        if r is None or w!=len(r): continue
        if not r:
            assert not bs[i]
            eq[i]=None; continue
        p=min(r,key=lambda j:(len(inc[j]),abs(r[j])!=1,j))
        a=r[p]; rr={j:v/a for j,v in r.items() if j!=p}; z=bs[i]/a
        saved.append((p,rr,z))
        for j in r: inc[j].discard(i)
        eq[i]=None
        for h in list(inc[p]):
            t=eq[h]; a=t.pop(p); inc[p].discard(h); bs[h]-=a*z
            for j,v in rr.items():
                q=t.get(j,0)-a*v
                if q:
                    if j not in t: inc[j].add(h)
                    t[j]=q
                elif j in t: del t[j]; inc[j].discard(h)
            heapq.heappush(hp,(len(t),h))
    x=[F(float(v)).limit_denominator(10**6) for v in approx]
    for p,r,z in reversed(saved): x[p]=z-sum(v*x[j] for j,v in r.items())
    return x

L=(1,1,2,3,3,3,3,4,4,4,5,5,5,5,8); s=900
roots={A for A,B in cuts(L)}; words=set(roots); er=[]
for Q in sorted(roots):
    for i,j,rr in fusions(Q):
        r=collections.Counter({Q:1})
        for R in rr: r[R]-=1; words.add(R)
        er.append({A:v for A,v in r.items() if v})
keys=sorted(Q for Q in words if len(Q)>=3); ix={Q:i for i,Q in enumerate(keys)}
E=sparse([{ix[Q]:v for Q,v in r.items() if Q in ix} for r in er],len(keys))
@functools.lru_cache(None)
def derivative(Q,t):
    cc=collections.Counter(Q)
    ng=dict(zip(cc,negative_counts(Q,t)))
    w={a:[int(v) for v in K(n)[ng[a],:]] for a,n in cc.items()}
    r=collections.Counter()
    for A,B in cuts(Q):
        ca=collections.Counter(A); z=math.prod(w[a][ca[a]] for a in cc)
        if A in ix: r[ix[A]]+=z*m(B)
        if B in ix: r[ix[B]]+=z*m(A)
    return {j:v for j,v in r.items() if v}
pa=derivative(L,s); C,t=child(L,s); ob=collections.Counter(pa)
for j,v in derivative(C,t).items(): ob[j]-=v
assert all(v%2==0 for v in ob.values())
ob={j:v//2 for j,v in ob.items() if v}
flips=[]; blocks=[]
for i,j,rr in fusions(L):
    f=derivative(L,s^(1<<i)^(1<<j)); union=pa.keys()|f.keys()
    flips.append({h:pa.get(h,0)-f.get(h,0) for h in union if pa.get(h,0)!=f.get(h,0)})
    blocks.append({h:-(pa.get(h,0)+f.get(h,0)) for h in union if pa.get(h,0)!=-f.get(h,0)})
extra=(16521,31929)
U=sparse(flips+blocks+[{j:-2*v for j,v in derivative(L,t0).items()} for t0 in extra],len(keys))
EE=vstack([E,sparse([ob],len(keys))]).tocsr()
b=np.zeros(EE.shape[0]); b[-1]=-1
log("LP matrices",E.shape,"inequalities",U.shape)
lp=linprog(np.ones(len(keys)),A_eq=EE,b_eq=b,A_ub=U,b_ub=np.zeros(U.shape[0]),
           bounds=(0,None),method="highs",options={"time_limit":60})
assert lp.success,lp.message
ss=np.flatnonzero(lp.x>1e-7); active=np.flatnonzero(abs(lp.ineqlin.residual)<1e-7)
log("rational reconstruction","support",len(ss),"active inequalities",len(active))
xx=exact(vstack([EE[:,ss],U[active,:][:,ss]]),np.r_[b,np.zeros(len(active))],lp.x[ss])
den=math.lcm(*(v.denominator for v in xx)); x=[0]*len(keys)
for j,v in zip(ss,xx): x[j]=int(v*den)
assert min(x)>=0 and mv(EE,x)==[int(z)*den for z in b] and max(mv(U,x))<=0
vt=dict(zip(keys,x))
assert all(not v for Q,v in vt.items() if len(Q)<=3 or sum(Q)%2)
log("exact scalar equations checked",E.shape[0],"direction support",sum(v!=0 for v in x))
cc=collections.Counter(L); labels=tuple(cc); counts=tuple(cc.values())
profiles=list(itertools.product(*(range(n+1) for n in counts)))
index={z:i for i,z in enumerate(profiles)}; edges=set()
for i,z in enumerate(profiles):
    for a in range(len(labels)):
        for b0 in range(a,len(labels)):
            for da in (-1,1):
                for db in (-1,1):
                    zz=list(z); zz[a]+=da; zz[b0]+=db
                    if a==b0:
                        if da!=db:
                            if not 0<z[a]<counts[a]: continue
                        elif da==1 and counts[a]-z[a]<2: continue
                        elif da==-1 and z[a]<2: continue
                    elif not (0<=z[a]+da<=counts[a] and 0<=z[b0]+db<=counts[b0]): continue
                    zz=tuple(zz)
                    if zz in index: edges.add(tuple(sorted((i,index[zz]))))
h=tensor(L,vt); flat=list(h.flat)
h0=tensor(L); flat0=list(h0.flat)
assert min(flat0[i]+flat0[j] for i,j in edges)>=0
assert len(edges)==19592 and min(flat[i]+flat[j] for i,j in edges)>=0
for t0 in extra:
    i=index[negative_counts(L,t0)]; assert (i,i) in edges
dparent=sum(v*x[j] for j,v in pa.items())
dchild=sum(v*x[j] for j,v in derivative(C,t).items())
log("ALL BLOCKS",len(edges),"min derivative of 2Q",min(flat[i]+flat[j] for i,j in edges))
log("LIFTED FEASIBLE POINT Delta",3914630-den)
log("SEPARATOR","dDelta",-den,"dPhi",dparent,"dPhiChild",dchild,
    "max derivative of 4D",max(mv(sparse(flips,len(keys)),x)))
kids=sorted({Q for i,j,rr in fusions(L) for Q in rr})
best=None
for Q in kids:
    hh=tensor(Q,vt); ind=np.unravel_index(np.argmin(hh),hh.shape)
    item=(int(hh[ind]),Q,tuple(int(z) for z in ind))
    if best is None or item<best: best=item
assert best[0]<0
log("PARTIAL CHANNEL WITNESS",best)
selected=[]
base=list(L); base.remove(4); base.remove(8)
for channel in range(4,13,2):
    Q=tuple(sorted(base+[channel]))
    t0=maskof(Q,{2:1,channel:1})
    hd=tensor(Q,vt); hv=tensor(Q)
    selected.append((channel,int(hd[negative_counts(Q,t0)]),int(hv[negative_counts(Q,t0)])))
assert sum(z for c,z,w in selected)>=0
log("SELECTED (4,8) BLOCK: channel, derivative, Phi",selected)
assert len(kids)==58
log("one-fusion children",len(kids),"signed profiles",
    sum(math.prod(n+1 for n in collections.Counter(Q).values()) for Q in kids))
for mask in (17,68,170,255):
    r=actual(tuple(range(1,9)),mask); assert r[2]>=0 and r[3]<0
    log("eight-label signing",mask,"g,child,Delta,maxD",r)
r=actual(tuple(range(1,11))+(13,),340)
assert r==(194087,3086,191001,F(-185))
log("eleven-factor","g,child,Delta,maxD",r)
for q in range(2,7):
    Q=(1,1,2)+(3,)*q+(4,)*(q-1)+(5,)*q+(8,)
    mask=maskof(Q,{2:1,4:q-1,8:q%2})
    r=actual(Q,mask); assert r[2]>=0
    if q==4: assert r==(4075371,160741,3914630,F(-20406))
    log("F1 multiplicity",q,"g,child,Delta,maxD",r)
for k in range(5,13):
    W=k*(k+1)//2; p=k+1+(W-k-1)%2
    Q=tuple(range(1,k+1))+(p,); mask=(1<<k)-1
    if k%2: mask|=1<<k
    r=actual(Q,mask); assert r[2]>=0
    log("all-minus background run",k,"p",p,"g,child,Delta,maxD",r)
if args.dump:
    for Q,v in vt.items():
        if v: print("DIRECTION",Q,v,flush=True)
log("ALL EXACT CHECKS PASS")
