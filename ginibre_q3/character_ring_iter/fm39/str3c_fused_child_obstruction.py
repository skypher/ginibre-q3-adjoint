import argparse,collections,functools,heapq,itertools,math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,vstack

ap=argparse.ArgumentParser(description="Exact depth-one identities and a fused-child necessity test; memory only.")
ap.add_argument("--dump",action="store_true",help="Print every identity coefficient and obstruction entry.")
args=ap.parse_args()

@functools.lru_cache(None)
def character(Q):
    if not Q:return {0:1}
    r=collections.defaultdict(int)
    for a,v in character(Q[:-1]).items():
        for b in range(abs(a-Q[-1]),a+Q[-1]+1,2):r[b]+=v
    return dict(r)
def m(Q):return character(Q).get(0,0)
def cuts(Q):
    for s in range(1<<len(Q)):
        yield tuple(a for i,a in enumerate(Q) if s>>i&1),tuple(a for i,a in enumerate(Q) if not s>>i&1)
def ucuts(Q):
    cc=tuple(collections.Counter(Q).items())
    for z in itertools.product(*(range(c+1) for a,c in cc)):
        yield tuple(a for (a,c),b in zip(cc,z) for _ in range(b)),tuple(a for (a,c),b in zip(cc,z) for _ in range(c-b))
def fused(Q,i,j):
    C=Q[:i]+Q[i+1:j]+Q[j+1:]
    return [tuple(sorted(C+((c,) if c else ()))) for c in range(abs(Q[i]-Q[j]),Q[i]+Q[j]+1,2)]
def ufused(Q):
    cc=collections.Counter(Q)
    for a in cc:
        for b in cc:
            if a>b or (a==b and cc[a]<2):continue
            i=Q.index(a);j=Q.index(b,i+1) if a==b else Q.index(b)
            yield i,j,fused(Q,i,j)
def child(L,s):
    i,j=max(((i,j) for i in range(len(L)-1) for j in range(i+1,len(L)-1) if (L[i]-L[j])%2==0),
            key=lambda ij:(L[ij[0]]+L[ij[1]],max(L[ij[0]],L[ij[1]])))
    keep=[h for h in range(len(L)) if h not in (i,j)]
    Q=tuple(L[h] for h in keep);t=sum(((s>>h)&1)<<q for q,h in enumerate(keep))
    if ((s>>i)^(s>>j))&1:t^=1<<(len(Q)-1)
    return Q,t
@functools.lru_cache(None)
def values(Q):
    z=[m(A)*m(B) for A,B in cuts(Q)]
    for k in range(len(Q)):
        for s in range(1<<len(Q)):
            if not s>>k&1:
                a,b=z[s],z[s|1<<k];z[s]=a+b;z[s|1<<k]=a-b
    return tuple(z)
def actual(L,s):
    Q,t=child(L,s);f=values(L);g=values(Q)
    delta=(f[s]-g[t])//2
    flips=[F(f[s]-f[s^(1<<i)^(1<<j)],4) for i in range(len(L)) for j in range(i+1,len(L))]
    assert max(flips)<0 and delta-m(L)<0
    print("ACTUAL",len(L),s,"g",f[s]//2,"child",g[t]//2,"Delta",delta,"m",m(L),"maxD",max(flips),"depth-zero",delta-m(L),flush=True)
    return delta
def mat(rows,n):
    ii=[];jj=[];vv=[]
    for i,r in enumerate(rows):
        for j,v in r.items():
            if v:
                assert -(1<<63)<v<(1<<63)
                ii.append(i);jj.append(j);vv.append(v)
    return coo_matrix((vv,(ii,jj)),shape=(len(rows),n),dtype=np.int64).tocsr()
def mv(A,x):
    return [sum(int(v)*x[j] for j,v in zip(A.indices[A.indptr[i]:A.indptr[i+1]],A.data[A.indptr[i]:A.indptr[i+1]])) for i in range(A.shape[0])]
def solve_exact(A,b,approx):
    A=A.tocsr();eq=[];bs=[];inc=[set() for _ in approx]
    for i in range(A.shape[0]):
        r={int(j):F(int(v)) for j,v in zip(A.indices[A.indptr[i]:A.indptr[i+1]],A.data[A.indptr[i]:A.indptr[i+1]]) if v}
        assert int(b[i])==b[i]
        z=F(int(b[i]))
        if not r:
            assert not z
            continue
        h=len(eq);eq.append(r);bs.append(z)
        for j in r:inc[j].add(h)
    hp=[(len(r),i) for i,r in enumerate(eq)];heapq.heapify(hp);saved=[]
    while hp:
        w,i=heapq.heappop(hp);r=eq[i]
        if r is None or w!=len(r):continue
        if not r:
            assert not bs[i]
            eq[i]=None;continue
        p=min(r,key=lambda j:(len(inc[j]),abs(r[j])!=1,j))
        a=r[p];rr={j:v/a for j,v in r.items() if j!=p};z=bs[i]/a
        saved.append((p,rr,z))
        for j in r:inc[j].discard(i)
        eq[i]=None
        for h in list(inc[p]):
            t=eq[h];a=t.pop(p);inc[p].discard(h);bs[h]-=a*z
            for j,v in rr.items():
                q=t.get(j,0)-a*v
                if q:
                    if j not in t:inc[j].add(h)
                    t[j]=q
                elif j in t:del t[j];inc[j].discard(h)
            heapq.heappush(hp,(len(t),h))
    x=[F(float(v)).limit_denominator(10**6) for v in approx]
    for p,r,z in reversed(saved):x[p]=z-sum(v*x[j] for j,v in r.items())
    return x

def eight():
    L=tuple(range(1,9));d={A:0 for A,B in cuts(L)}
    for Q in list(d):
        for i in range(len(Q)):
            for j in range(i+1,len(Q)):
                for R in fused(Q,i,j):d.setdefault(R,1)
    coords=[];ix={}
    def var(A,B):
        if len(A)==2:A=() if A[0]==A[1] else (-1,)
        if len(B)==2:B=() if B[0]==B[1] else (-1,)
        if len(A)==1 or len(B)==1:return None
        t=(A,B) if A<=B else (B,A)
        if t not in ix:ix[t]=len(coords);coords.append(t)
        return ix[t]
    terms={}
    for Q in d:
        terms[Q]=[(j,s) for s,(A,B) in enumerate(cuts(Q)) if (j:=var(A,B)) is not None]
    rows=[{var((),()):1}];tags=[("unit",)]
    for j,(A,B) in enumerate(coords):
        for side in (0,1):
            C,D=(A,B) if not side else (B,A)
            if d[C]:continue
            for i in range(len(C)):
                for h in range(i+1,len(C)):
                    r=collections.Counter({j:1})
                    for R in fused(C,i,h):
                        q=var(R,D)
                        if q is not None:r[q]-=1
                    rows.append(r);tags.append((A,B,side,i,h))
    E=mat(rows,len(coords));be=np.zeros(E.shape[0]);be[0]=1
    def phi(Q,s):
        r=collections.Counter()
        for j,t in terms[Q]:r[j]+=(-1 if (s&t).bit_count()%2 else 1)
        return r
    for s,wanted in ((17,143),(68,143),(170,155),(255,155)):
        delta=actual(L,s);Q,t=child(L,s);r=phi(L,s)
        for j,v in phi(Q,t).items():r[j]-=v
        c=np.zeros(len(coords),dtype=np.int64)
        for j,v in r.items():
            assert v%2==0
            c[j]=v//2
        lp=linprog(c,A_eq=E,b_eq=be,bounds=(0,None),method="highs",options={"time_limit":120})
        assert lp.success,lp.message
        support=np.flatnonzero(abs(lp.eqlin.marginals)>1e-8)
        active=np.flatnonzero(abs(lp.lower.marginals)<1e-7)
        V=E[support,:].tocsc()
        lam=solve_exact(V[:,active].T,c[active],lp.eqlin.marginals[support])
        den=math.lcm(*(v.denominator for v in lam));z=[int(v)*den for v in c]
        for j in range(len(coords)):
            for h,v in zip(V.indices[V.indptr[j]:V.indptr[j+1]],V.data[V.indptr[j]:V.indptr[j+1]]):
                z[j]-=int(v)*int(lam[h]*den)
        q=lam[list(support).index(0)] if 0 in support else F(0)
        assert q==wanted and min(z)>=0
        assert q+sum(F(v,den)*m(A)*m(B) for v,(A,B) in zip(z,coords))==delta
        print("EIGHT IDENTITY",s,"constant",q,"fusion",len(support)-(0 in support),"products",sum(v>0 for v in z),"den",den,"shorter=0 flips=0",flush=True)
        if args.dump:
            for j,v in zip(support,lam):
                if j and v:print("FUSION",s,tags[j],v)
            for j,v in enumerate(z):
                if v:print("PRODUCT",s,coords[j],F(v,den))

def maskof(Q,negative):
    seen=collections.Counter();s=0
    for i,a in enumerate(Q):
        if seen[a]<negative.get(a,0):s|=1<<i;seen[a]+=1
    return s
@functools.lru_cache(None)
def K(n):
    return np.array([[sum((-1)**h*math.comb(t,h)*math.comb(n-t,j-h)
                         for h in range(max(0,j-n+t),min(j,t)+1))
                      for j in range(n+1)] for t in range(n+1)],dtype=object)
def transformed(Q,vtab,also_value=False):
    cc=collections.Counter(Q);dims=tuple(c+1 for c in cc.values())
    vals=[];base=[]
    for A,B in ucuts(Q):
        vals.append(vtab.get(A,0)*m(B)+m(A)*vtab.get(B,0))
        if also_value:base.append(m(A)*m(B))
    h=np.array(vals,dtype=object).reshape(dims)
    h0=np.array(base,dtype=object).reshape(dims) if also_value else None
    for ax,c in enumerate(cc.values()):
        h=np.moveaxis(np.tensordot(K(c),h,axes=(1,ax)),0,ax)
        if also_value:h0=np.moveaxis(np.tensordot(K(c),h0,axes=(1,ax)),0,ax)
    return h,h0

def f1_obstruction():
    L=(1,1,2,3,3,3,3,4,4,4,5,5,5,5,8);s=900
    actual(L,s);roots={A for A,B in ucuts(L)};er=[];allwords=set(roots)
    for Q in roots:
        for i,j,rs in ufused(Q):
            r=collections.Counter({Q:1})
            for R in rs:r[R]-=1;allwords.add(R)
            er.append({A:v for A,v in r.items() if v})
    keys=sorted(Q for Q in allwords if len(Q)>=3);ix={Q:i for i,Q in enumerate(keys)}
    E=mat([{ix[Q]:v for Q,v in r.items() if Q in ix} for r in er],len(keys))
    @functools.lru_cache(None)
    def derivative(Q,t):
        cc=collections.Counter(Q);neg=collections.Counter(a for i,a in enumerate(Q) if t>>i&1)
        weights={}
        for a,n in cc.items():weights[a]=[int(v) for v in K(n)[neg[a],:]]
        r=collections.Counter()
        for A,B in ucuts(Q):
            counts=collections.Counter(A);w=math.prod(weights[a][counts[a]] for a in cc)
            if A in ix:r[ix[A]]+=w*m(B)
            if B in ix:r[ix[B]]+=w*m(A)
        return {j:v for j,v in r.items() if v}
    parent=derivative(L,s);Q,t=child(L,s);c=collections.Counter(parent)
    for j,v in derivative(Q,t).items():c[j]-=v
    assert all(v%2==0 for v in c.values())
    c={j:v//2 for j,v in c.items() if v}
    flips=[];seen=set()
    for i in range(len(L)):
        for j in range(i+1,len(L)):
            sig=tuple(sorted(((L[i],(s>>i)&1),(L[j],(s>>j)&1))))
            if sig in seen:continue
            seen.add(sig);r=collections.Counter(parent)
            for h,v in derivative(L,s^(1<<i)^(1<<j)).items():r[h]-=v
            flips.append(dict(r))
    U=list(flips)
    tests="""1 1 -3 -3 -3 -3 -4 -4 -4 -5 -5 -5 -5 -8
1 3 4 4 4 -5 -5 5 13
-1 1 -4 -4 4 4 -5 5 5 5
-1 1 3 -4 -4 4 4 -5 5 5
1 1 2 3 3 -5 -5 -5 -13
-1 1 2 2 4 4 -5 5 5 5
-1 1 1 3 -4 -4 -5 5 5 5
-1 1 2 2 3 3 -5 5 5 5
-1 1 2 2 3 4 4 -5 5 5
-1 1 1 2 3 4 -5 5 5 5
-1 1 1 -3 3 3 4 4 4 8
1 1 4 4 -5 -5 5 13
3 -4 -4 -4 5 5 -13
-1 1 1 2 3 3 4 -5 5 5
-1 1 -2 2 3 3 4 4 5 5
-1 1 1 3 3 3 3 -5 5 5
-1 1 2 2 4 -5 5 5 5
-1 1 1 3 4 -5 5 5 5
-1 1 1 2 3 3 3 4 -5 5
-1 1 1 2 4 4 -5 5 5
-1 1 1 3 3 3 3 4 4 -5
-1 1 2 2 3 4 -5 5 5
-1 1 1 2 3 -5 5 5 5
-1 1 2 2 3 3 3 4 4 -5
2 -3 3 -4 4 4 5 5 5 13
-3 -3 3 3 -4 4 5 5 8 -10
-3 -3 3 3 4 4 5 5 5 13
-1 1 3 4 4 5 5 5 5 -11
2 -3 -3 3 3 4 5 5 5 13
2 -3 -3 3 3 5 5 5 8 9
-1 1 3 3 3 5 5 5 5 -11
-1 1 2 3 4 5 5 5 5 -11
-1 1 2 3 5 5 5 5 7 -8
-1 1 2 4 4 4 5 5 5 -11
1 1 -2 4 4 5 5 5 5 -12
-1 1 2 4 4 5 5 5 7 -8
-1 1 2 4 5 5 5 5 6 -8
-1 1 3 4 5 5 5 5 5 -8
1 2 -3 4 4 5 5 5 5 -12
2 -3 3 4 4 5 5 5 5 -12
2 3 -4 4 4 5 5 5 5 -11
2 -3 -3 3 3 5 5 5 5 12
2 -3 3 -4 4 5 5 5 5 12
-1 1 2 3 3 5 5 5 5 -12
-3 -3 3 3 -4 4 5 5 5 -13
-1 1 2 3 3 5 5 5 5 -10
-1 1 3 3 3 3 5 5 5 -13"""
    for line in tests.splitlines():
        word=tuple(map(int,line.split()));Q=tuple(abs(a) for a in word)
        t=sum((a<0)<<i for i,a in enumerate(word))
        U.append({j:-v for j,v in derivative(Q,t).items()})
    UU=mat(U,len(keys));EE=vstack([E,mat([c],len(keys))]).tocsr()
    be=np.zeros(EE.shape[0]);be[-1]=-1
    lp=linprog(np.ones(len(keys)),A_eq=EE,b_eq=be,A_ub=UU,b_ub=np.zeros(UU.shape[0]),bounds=(0,None),method="highs",options={"time_limit":120})
    assert lp.success,lp.message
    ss=np.flatnonzero(lp.x>1e-7);active=np.flatnonzero(abs(lp.ineqlin.residual)<1e-7)
    xx=solve_exact(vstack([EE[:,ss],UU[active,:][:,ss]]),np.r_[be,np.zeros(len(active))],lp.x[ss])
    den=math.lcm(*(z.denominator for z in xx));x=[0]*len(keys)
    for j,z in zip(ss,xx):x[j]=int(z*den)
    assert min(x)>=0 and mv(EE,x)==[int(z)*den for z in be] and max(mv(UU,x))<=0
    vt=dict(zip(keys,x));count=0
    for Q in allwords:
        if len(Q)==len(L) or (Q not in roots and len(Q)>10):continue
        h,h0=transformed(Q,vt,True)
        assert min(h.flat)>=0 and min(h0.flat)>=0,(Q,min(h.flat))
        count+=h.size
    assert count==700367
    Q=(1,1,2,3,3,3,4,4,4,5,5,5,5,11)
    t=maskof(Q,{3:3,4:3,5:4})
    bad=sum(v*x[j] for j,v in derivative(Q,t).items())
    assert bad<0
    print("F1 EXACT OBSTRUCTION","scalar equations",E.shape[0],"support",sum(bool(z) for z in x),"den",den,"Delta derivative",-den,"max4D",max(mv(mat(flips,len(keys)),x)),"sum",sum(x),flush=True)
    print("DELETION CHILDREN AND FUSED WORDS <=10:",count,"signed profiles checked exactly",flush=True)
    print("FUSED FORM",Q,"negative counts",{3:3,4:3,5:4},"derivative",bad,flush=True)
    if args.dump:
        for Q,z in vt.items():
            if z:print("DIRECTION",Q,z)

eight()
actual(tuple(range(1,11))+(13,),340)
f1_obstruction()
print("ALL EXACT CHECKS PASS",flush=True)
