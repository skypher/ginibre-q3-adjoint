import argparse,collections,functools,heapq,itertools,math,time
from fractions import Fraction as F
from types import SimpleNamespace as NS
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,vstack

ap=argparse.ArgumentParser(description="Exact fusion-depth certificates; no files are written.")
ap.add_argument("--dump",action="store_true",help="Print every rational identity coefficient and ray.")
args=ap.parse_args()

@functools.lru_cache(None)
def mult(Q):
    r={0:1}
    for a in Q:
        s=collections.defaultdict(int)
        for b,x in r.items():
            for c in range(abs(a-b),a+b+1,2):s[c]+=x
        r=s
    return r.get(0,0)

def cuts(Q):
    for s in range(1<<len(Q)):
        yield tuple(a for i,a in enumerate(Q) if s>>i&1),tuple(a for i,a in enumerate(Q) if not s>>i&1)

def fuse(Q,i,j):
    R=Q[:i]+Q[i+1:j]+Q[j+1:]
    return [tuple(sorted(R+((c,) if c else ()))) for c in range(abs(Q[i]-Q[j]),Q[i]+Q[j]+1,2)]

def words_at(L,k):
    d={A:0 for A,B in cuts(L)}
    for h in range(k):
        new={}
        for Q,z in list(d.items()):
            if z!=h:continue
            for i in range(len(Q)):
                for j in range(i+1,len(Q)):
                    for R in fuse(Q,i,j):
                        if R not in d:new[R]=h+1
        d.update(new)
    return d

def character_values(Q):
    f=[mult(A)*mult(B) for A,B in cuts(Q)]
    for k in range(len(Q)):
        for s in range(1<<len(Q)):
            if not s>>k&1:
                a,b=f[s],f[s|1<<k]
                f[s]=a+b;f[s|1<<k]=a-b
    return f

def child(L,mask):
    n=len(L)
    i,j=max(((i,j) for i in range(n-1) for j in range(i+1,n-1)
             if (L[i]-L[j])%2==0),
            key=lambda ij:(L[ij[0]]+L[ij[1]],max(L[ij[0]],L[ij[1]])))
    keep=[h for h in range(n) if h not in (i,j)]
    C=tuple(L[h] for h in keep)
    cm=sum(((mask>>h)&1)<<q for q,h in enumerate(keep))
    if ((mask>>i)^(mask>>j))&1:cm^=1<<(len(C)-1)
    return C,cm,(i,j)

def sparse(rows,n):
    rr=[];cc=[];vv=[]
    for i,row in enumerate(rows):
        for j,v in row.items():
            if v:rr.append(i);cc.append(j);vv.append(v)
    return coo_matrix((vv,(rr,cc)),shape=(len(rows),n),dtype=np.int64).tocsr()

def mv(M,x):
    return [sum(int(v)*x[j] for j,v in
                zip(M.indices[M.indptr[i]:M.indptr[i+1]],
                    M.data[M.indptr[i]:M.indptr[i+1]]))
            for i in range(M.shape[0])]

def model(L,k,mode,pairs=False):
    tick=time.monotonic();d=words_at(L,k);coords=[];index={}
    def var(A,B):
        if pairs and len(A)==2:A=() if A[0]==A[1] else (-1,)
        if pairs and len(B)==2:B=() if B[0]==B[1] else (-1,)
        if len(A)==1 or len(B)==1:return None
        q=(A,B) if A<=B else (B,A)
        if q not in index:index[q]=len(coords);coords.append(q)
        return index[q]
    terms={}
    for Q in d:
        terms[Q]=[(j,s) for s,(A,B) in enumerate(cuts(Q))
                  if (j:=var(A,B)) is not None]
    E=[{var((),()):1}];et=[("unit",)]
    for j,(A,B) in enumerate(coords):
        for side in (0,1):
            C,D=(A,B) if side==0 else (B,A)
            level=d[C] if mode=="component" else d[tuple(sorted(A+B))]
            if level>=k:continue
            assert d[C]<k
            for i in range(len(C)):
                for h in range(i+1,len(C)):
                    row=collections.Counter({j:1})
                    for R in fuse(C,i,h):
                        z=var(R,D)
                        if z is not None:row[z]-=1
                    E.append(row);et.append((A,B,side,i,h))
    def phi(Q,mask):
        r=collections.Counter()
        for j,s in terms[Q]:r[j]+=(-1 if (s&mask).bit_count()%2 else 1)
        return r
    P=[];pt=[]
    for Q in d:
        if len(Q)>=len(L):continue
        for s in range(1<<len(Q)):
            if s.bit_count()%2:continue
            P.append({j:-v for j,v in phi(Q,s).items()})
            pt.append((Q,s))
    N=len(coords)
    M=NS(L=L,k=k,mode=mode,pairs=pairs,d=d,coords=coords,index=index,phi=phi,N=N,
         E=sparse(E,N),P=sparse(P,N),et=et,pt=pt)
    print("MODEL",L,k,mode,"pair anchors",pairs,N,M.E.shape[0],M.P.shape[0],
          "seconds",round(time.monotonic()-tick,2),flush=True)
    return M

def profile(M,mask):
    L=M.L;n=len(L);C,cm,pair=child(L,mask)
    p=M.phi(L,mask);obj=p.copy()
    for j,v in M.phi(C,cm).items():obj[j]-=v
    assert all(v%2==0 for v in obj.values())
    c=np.zeros(M.N,dtype=np.int64)
    for j,v in obj.items():c[j]=v//2
    rows=[];tags=[]
    for i in range(n):
        for j in range(i+1,n):
            r=p.copy()
            for z,v in M.phi(L,mask^(1<<i)^(1<<j)).items():r[z]-=v
            rows.append(r);tags.append((i,j))
    return c,sparse(rows,M.N),tags

def exact_solve(mat,right,approx):
    mat=mat.tocsr();eq=[];bs=[];inc=[set() for _ in approx]
    for i in range(mat.shape[0]):
        row={int(j):F(int(v)) for j,v in
             zip(mat.indices[mat.indptr[i]:mat.indptr[i+1]],
                 mat.data[mat.indptr[i]:mat.indptr[i+1]]) if v}
        assert int(right[i])==right[i]
        b=F(int(right[i]))
        if not row:
            assert not b
            continue
        h=len(eq);eq.append(row);bs.append(b)
        for j in row:inc[j].add(h)
    heap=[(len(r),i) for i,r in enumerate(eq)];heapq.heapify(heap)
    saved=[]
    while heap:
        width,i=heapq.heappop(heap);row=eq[i]
        if row is None or width!=len(row):continue
        if not row:
            assert not bs[i]
            eq[i]=None;continue
        p=min(row,key=lambda j:(len(inc[j]),abs(row[j])!=1,j))
        a=row[p];r={j:v/a for j,v in row.items() if j!=p};b=bs[i]/a
        saved.append((p,r,b))
        for j in row:inc[j].discard(i)
        eq[i]=None
        for h in list(inc[p]):
            rr=eq[h];a=rr.pop(p);inc[p].discard(h);bs[h]-=a*b
            for j,v in r.items():
                z=rr.get(j,0)-a*v
                if z:
                    if j not in rr:inc[j].add(h)
                    rr[j]=z
                elif j in rr:del rr[j];inc[j].discard(h)
            heapq.heappush(heap,(len(rr),h))
    x=[F(float(v)).limit_denominator(10**6) for v in approx]
    for p,r,b in reversed(saved):x[p]=b-sum(v*x[j] for j,v in r.items())
    return x

def check_point(M,U,c,x,wanted):
    b=[0]*M.E.shape[0];b[0]=1
    assert min(x)>=0 and mv(M.E,x)==b and max(mv(U,x),default=0)<=0
    assert sum(int(v)*z for v,z in zip(c,x))==wanted

def certify(M,mask,wanted,use_shorter=True,with_flips=False):
    c,flips,ft=profile(M,mask)
    U=M.P if use_shorter else sparse([],M.N)
    npterms=U.shape[0]
    if with_flips:U=vstack([U,flips]).tocsr()
    b=np.zeros(M.E.shape[0]);b[0]=1
    print("SOLVE_DUAL",M.L,mask,M.k,flush=True)
    r=linprog(c,A_ub=U if U.shape[0] else None,
              b_ub=np.zeros(U.shape[0]) if U.shape[0] else None,
              A_eq=M.E,b_eq=b,bounds=(0,None),method="highs",
              options={"time_limit":180})
    assert r.success,r.message
    iy=np.flatnonzero(abs(r.eqlin.marginals)>1e-8)
    iz=np.flatnonzero(abs(r.ineqlin.marginals)>1e-8)
    V=vstack([M.E[iy,:],U[iz,:]]).tocsc()
    active=np.flatnonzero(abs(r.lower.marginals)<1e-7)
    a=exact_solve(V[:,active].T,c[active],
                  np.r_[r.eqlin.marginals[iy],r.ineqlin.marginals[iz]])
    den=math.lcm(*(v.denominator for v in a))
    z=[int(v)*den for v in c];ai=[int(v*den) for v in a]
    for j in range(M.N):
        for h,v in zip(V.indices[V.indptr[j]:V.indptr[j+1]],
                       V.data[V.indptr[j]:V.indptr[j+1]]):
            z[j]-=int(v)*ai[h]
    q=a[list(iy).index(0)] if 0 in iy else F(0)
    assert q==wanted and min(z)>=0 and all(v<=0 for v in a[len(iy):])
    x=[0]*M.N
    if wanted==0:x[M.index[((),())]]=1
    else:
        x=[mult(A)*mult(B) for A,B in M.coords]
        if sum(int(v)*z for v,z in zip(c,x))!=wanted:
            ss=np.flatnonzero(r.x>1e-8);tt=np.flatnonzero(abs(r.ineqlin.residual)<1e-7)
            xx=exact_solve(vstack([M.E[:,ss],U[tt,:][:,ss]]),
                           np.r_[b,np.zeros(len(tt))],r.x[ss])
            x=[F(0)]*M.N
            for j,v in zip(ss,xx):x[j]=v
    check_point(M,U,c,x,wanted)
    assert max(mv(flips,x),default=0)<=0
    print("EXACT_DUAL",M.L,mask,"bound",q,"fusion terms",
          len(iy)-(0 in iy),"shorter terms",len(iz),
          "product terms",sum(v>0 for v in z),flush=True)
    alpha=[-4*a[len(iy)+q] for q,j in enumerate(iz) if j>=npterms]
    assert not any(alpha)
    print("All flip coefficients are zero.",flush=True)
    if args.dump:
        print("CONSTANT",q)
        for j,v in zip(iy,a[:len(iy)]):
            if j:print("FUSION",M.et[j],v)
        for j,v in zip(iz,a[len(iy):]):
            if j<npterms:print("SHORTER",M.pt[j],-v)
            else:print("FLIP",ft[j-npterms],-4*v)
        for j,v in enumerate(z):
            if v:print("PRODUCT",M.coords[j],F(v,den))
    return q

def verify_ray(M,mask,ray,den):
    c,flips,ft=profile(M,mask);U=vstack([M.P,flips]).tocsr()
    assert min(ray)>=0 and not any(mv(M.E,ray))
    assert max(mv(U,ray),default=0)<=0
    value=F(sum(int(v)*x for v,x in zip(c,ray)),den)
    assert value<0
    print("EXACT_NEGATIVE_RAY",M.L,mask,value,"denominator",den,
          "support",sum(v!=0 for v in ray),flush=True)

def find_ray(M,mask):
    c,flips,ft=profile(M,mask);U=vstack([M.P,flips]).tocsr()
    E=vstack([M.E,coo_matrix(c.reshape(1,-1))]).tocsr()
    b=np.zeros(E.shape[0]);b[-1]=-1
    print("SOLVE_RAY",M.L,mask,flush=True)
    r=linprog(np.ones(M.N),A_ub=U,b_ub=np.zeros(U.shape[0]),
              A_eq=E,b_eq=b,bounds=(0,None),method="highs",
              options={"time_limit":180})
    assert r.success,r.message
    ss=np.flatnonzero(r.x>1e-8);tt=np.flatnonzero(abs(r.ineqlin.residual)<1e-7)
    a=exact_solve(vstack([E[:,ss],U[tt,:][:,ss]]),
                  np.r_[b,np.zeros(len(tt))],r.x[ss])
    den=math.lcm(*(v.denominator for v in a));ray=[0]*M.N
    for j,v in zip(ss,a):ray[j]=int(v*den)
    verify_ray(M,mask,ray,den)
    if args.dump:
        for j,v in enumerate(ray):
            if v:print("RAY",M.coords[j],F(v,den))
    return ray,den

def first_jet():
    L=(1,3,8,9,10,11);mask=63;d=words_at(L,1)
    keys=[Q for Q in d if len(Q)>=2];ix={Q:i for i,Q in enumerate(keys)}
    N=len(keys);E=[]
    for Q,h in d.items():
        if h:continue
        for i in range(len(Q)):
            for j in range(i+1,len(Q)):
                r=collections.Counter()
                if Q in ix:r[ix[Q]]+=1
                for R in fuse(Q,i,j):
                    if R in ix:r[ix[R]]-=1
                E.append(r)
    def derivative(Q,s):
        r=collections.Counter()
        for t,(A,B) in enumerate(cuts(Q)):
            sg=-1 if (s&t).bit_count()%2 else 1
            if A in ix:r[ix[A]]+=sg*mult(B)
            if B in ix:r[ix[B]]+=sg*mult(A)
        return r
    U=[]
    for Q in d:
        if len(Q)>=len(L):continue
        for s in range(1<<len(Q)):
            if not s.bit_count()%2:U.append({j:-v for j,v in derivative(Q,s).items()})
    p=derivative(L,mask);C,cm,pair=child(L,mask);obj=p.copy()
    for j,v in derivative(C,cm).items():obj[j]-=v
    for i in range(len(L)):
        for j in range(i+1,len(L)):
            r=p.copy()
            for z,v in derivative(L,mask^(1<<i)^(1<<j)).items():r[z]-=v
            U.append(r)
    c=np.zeros(N,dtype=np.int64)
    for j,v in obj.items():assert v%2==0;c[j]=v//2
    EE=sparse(E,N);UU=sparse(U,N)
    T=vstack([EE,coo_matrix(c.reshape(1,-1))]).tocsr()
    b=np.zeros(T.shape[0]);b[-1]=-1
    r=linprog(np.ones(N),A_ub=UU,b_ub=np.zeros(UU.shape[0]),
              A_eq=T,b_eq=b,bounds=(0,None),method="highs")
    assert r.success
    x=[F(float(v)).limit_denominator(10**6) for v in r.x]
    assert min(x)>=0 and not any(mv(EE,x)) and max(mv(UU,x))<=0
    assert sum(int(v)*z for v,z in zip(c,x))==-1
    v=dict(zip(keys,x))
    print("EXACT_FIRST_JET","support",sum(bool(z) for z in x),
          "fusion equations",EE.shape[0],"sign tests",UU.shape[0],
          "TopPair derivative -1",flush=True)
    if args.dump:print("FIRST_JET",[(Q,z) for Q,z in v.items() if z])
    return v

def actual(L,mask):
    f=character_values(L);C,cm,pair=child(L,mask);g=character_values(C)
    drop=F(f[mask]-g[cm],2)
    ds=[F(f[mask]-f[mask^(1<<i)^(1<<j)],4)
        for i in range(len(L)) for j in range(i+1,len(L))]
    assert max(ds)<0 and drop>0 and drop-mult(L)<0
    print("ACTUAL",L,mask,"g",f[mask]//2,"child",g[cm]//2,
          "drop",drop,"m",mult(L),"max D",max(ds),flush=True)
    return drop

L6=(1,3,8,9,10,11);L7=(1,2,3,4,7,8,9);L8=tuple(range(1,9))
for L,s in [(L6,63),(L7,15)]+[(L8,s) for s in (17,68,170,255)]:
    actual(L,s)
f=character_values(L8)
mins=[s for s in range(256) if not s.bit_count()%2 and
      all(f[s]<f[s^(1<<i)^(1<<j)] for i in range(8) for j in range(i+1,8))]
assert mins==[17,68,170,255]
v=first_jet()
M=model(L6,1,"component")
ray=[mult(A)*v.get(B,0)+v.get(A,0)*mult(B) for A,B in M.coords]
assert all(z.denominator==1 if isinstance(z,F) else True for z in ray)
verify_ray(M,63,[int(z) for z in ray],1)
certify(model(L6,2,"total"),63,62,use_shorter=False)
certify(model(L7,1,"component"),15,0)
M1=model(L8,1,"component")
ray,den=find_ray(M1,17)
verify_ray(M1,68,ray,den)
certify(M1,170,0);certify(M1,255,0)
M2=model(L8,2,"total")
certify(M2,17,460);certify(M2,68,460)
certify(model(L6,1,"component",pairs=True),63,1,use_shorter=False)
MP=model(L8,1,"component",pairs=True)
certify(MP,17,167,with_flips=True)
certify(MP,68,167,with_flips=True)
actual((1,1,2,3,3,3,3,4,4,4,5,5,5,5,8),900)
actual(tuple(range(1,11))+(13,),340)
print("ALL EXACT CHECKS PASS",flush=True)
