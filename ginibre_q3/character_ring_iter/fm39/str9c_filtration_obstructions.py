import sys
if any(a in ("-h","--help") for a in sys.argv[1:]):
    print("Usage: python3 -u -B verifier.py; exact filtration census and CG obstructions; no file writes.")
    raise SystemExit
from fractions import Fraction as Q
import sympy as sp

import sys
from datetime import datetime,timezone
from functools import lru_cache
from itertools import combinations_with_replacement,product,combinations
from math import comb,prod
def log(*a):print(datetime.now(timezone.utc).isoformat(timespec="seconds"),*a,flush=True)
@lru_cache(None)
def char(w):
    out={(0,0):1}
    for v in w:
        n=abs(v);eps=1 if v>0 else -1;new={}
        for (a,b),z in out.items():
            for c in range(abs(a-n),a+n+1,2):new[c,b]=new.get((c,b),0)+z
            for c in range(abs(b-n),b+n+1,2):new[a,c]=new.get((a,c),0)+eps*z
        out={ab:z for ab,z in new.items() if z}
    return out
@lru_cache(None)
def sing(ns):
    out={0:1}
    for n in ns:
        new={}
        for a,z in out.items():
            for c in range(abs(a-n),a+n+1,2):new[c]=new.get(c,0)+z
        out=new
    return out.get(0,0)
def hdata(w):
    A=char(w[::2]);B=char(w[1::2]);v={ab:z*B.get(ab,0) for ab,z in A.items() if z*B.get(ab,0)}
    return sum(max(z,0) for z in v.values()),sum(max(-z,0) for z in v.values()),v
def cutpoly(w):
    counts=[(v,w.count(v)) for v in dict.fromkeys(w)]
    out=[0]*(len(w)+1)
    for ks in product(*(range(c+1) for v,c in counts)):
        X=tuple(sorted(abs(v) for (v,c),k in zip(counts,ks) for _ in range(c-k)))
        Y=tuple(sorted(abs(v) for (v,c),k in zip(counts,ks) for _ in range(k)))
        z=prod(comb(c,k)*((-1)**k if v<0 else 1) for (v,c),k in zip(counts,ks))
        out[sum(ks)]+=z*sing(X)*sing(Y)
    return out
def pairs(w):
    ans=[]
    for i,j in combinations(range(len(w)),2):
        C=w[:i]+w[i+1:j]+w[j+1:]
        mixed=2*(1 if w[j]>0 else -1)*char(C).get((abs(w[i]),abs(w[j])),0)
        ans.append((i,j,mixed))
    return ans

log("start Euler filtration census")
counts={k:0 for k in ("all","corrected","cut_negative","height_negative","casimir_negative","some_mixed_negative","all_mixed_negative")}
first={}
corrected=[]
for L in range(11):
    for ns in combinations_with_replacement(range(1,5),L):
        labels=tuple(sorted(set(ns)))
        for sig in range(1<<len(labels)):
            signs={a:(-1 if sig>>h&1 else 1) for h,a in enumerate(labels)}
            w=tuple(signs[a]*a for a in ns)
            if sum(v<0 for v in w)%2:continue
            counts["all"]+=1
            he,ho,v=hdata(w)
            if not ho:continue
            counts["corrected"]+=1;corrected.append(w)
            cuts=cutpoly(w)
            assert sum(cuts)==he-ho
            hh={};cc={}
            for (a,b),z in v.items():
                hh[a+b]=hh.get(a+b,0)+z
                key=a*(a+2)+b*(b+2);cc[key]=cc.get(key,0)+z
            pp=pairs(w)
            flags={"cut_negative":min(cuts)<0,
                   "height_negative":min(hh.values())<0,
                   "casimir_negative":min(cc.values())<0,
                   "some_mixed_negative":any(z<0 for i,j,z in pp),
                   "all_mixed_negative":all(z<0 for i,j,z in pp)}
            for key,yes in flags.items():
                if yes:
                    counts[key]+=1
                    first.setdefault(key,(w,(he,ho),cuts,hh,pp))
    log("L",L,counts,"cache",char.cache_info().currsize)
log("COUNTS",counts)

assert counts=={"all":4521,"corrected":531,"cut_negative":463,
               "height_negative":44,"casimir_negative":511,
               "some_mixed_negative":377,"all_mixed_negative":0}
assert first["height_negative"][0]==(1,)*5+(-2,3,-4)
assert first["height_negative"][3][6]==-6
assert first["cut_negative"][0]==(-1,-1,-1,-2,3,4)
assert first["cut_negative"][2][3]==-6
log("PASS exact Euler obstructions on all 531 corrected profiles")

for t in range(6):
    w=(1,)*(2*t)+(-2,4,-6);he,ho,v=hdata(w);hh={}
    for (a,b),z in v.items():hh[a+b]=hh.get(a+b,0)+z
    cp=cutpoly(w)
    log("LAMBDA",t,"E1",(he,ho),"negative cut grades",
        [(j,z) for j,z in enumerate(cp) if z<0],
        "negative heights",[(j,z) for j,z in hh.items() if z<0])
expected=[(40,4),(196,12),(1080,100),(7256,822),(53872,7664)]
for k in range(5,10):
    W=k*(k+1)//2;p=max(k,6)
    while (W-p)%2 or (k%2==0 and p<=k):p+=1
    w=tuple(-n for n in range(1,k+1))+(((-1)**k)*p,)
    he,ho,v=hdata(w);assert (he,ho)==expected[k-5]
    hh={}
    for (a,b),z in v.items():hh[a+b]=hh.get(a+b,0)+z
    pp=pairs(w);cp=cutpoly(w)
    assert min(cp)<0 and min(hh.values())>=0
    log("RUN",k,"E1",(he,ho),"cut minimum",min(cp),
        "mixed Euler range",(min(z for i,j,z in pp),max(z for i,j,z in pp)))
strict=(-1,-2,-3,-4,-5,-6,-7,-8)
assert char(strict).get((0,0),0)==980
pp=pairs(strict)
assert (min(z for i,j,z in pp),max(z for i,j,z in pp))==(-130,-16)
f1=(1,1,-2)+(3,)*4+(-4,)*3+(5,)*4+(8,)
assert char(f1).get((0,0),0)==8150742
pp=pairs(f1)
assert (min(z for i,j,z in pp),max(z for i,j,z in pp))==(-9134096,-40812)
log("PASS every pair has negative mixed Euler on the strict root and F1")

def clean(p): return {v:c for v,c in p.items() if c}
def lower(p,ns):
    q={}
    for v,c in p.items():
        for i,n in enumerate(ns):
            if v[i]<n:
                w=list(v);w[i]+=1;w=tuple(w)
                q[w]=q.get(w,0)+c*(n-v[i])
    return clean(q)
def raise_op(p):
    q={}
    for v,c in p.items():
        for i,h in enumerate(v):
            if h:
                w=list(v);w[i]-=1;w=tuple(w)
                q[w]=q.get(w,0)+c*h
    return clean(q)
def inner(p,q,ns):
    if len(p)>len(q):p,q=q,p
    return sum(c*q.get(v,0)/prod(comb(n,h) for n,h in zip(ns,v)) for v,c in p.items())
@lru_cache(None)
def cg(ns):
    # Entries: final spin label, fusion path, normalized monomial multiplet.
    if not ns:return ((0,(),({():Q(1)},)),)
    old,n=ns[:-1],ns[-1];out=[]
    for a,path,states in cg(old):
        for c in range(abs(a-n),a+n+1,2):
            j=(a+n-c)//2;hv={}
            for h in range(j+1):
                for v,z in states[h].items():
                    w=v+(j-h,)
                    hv[w]=hv.get(w,0)+(-1)**h*comb(j,h)*z
            hv=clean(hv);assert not raise_op(hv)
            st=[hv]
            for h in range(c):
                st.append({v:z/Q(c-h) for v,z in lower(st[-1],ns).items()})
            assert not lower(st[-1],ns)
            out.append((c,path+(c,),tuple(st)))
    assert sum(c+1 for c,_,_ in out)==prod(n+1 for n in ns)
    for i,(a,_,st) in enumerate(out):
        assert inner(st[0],st[0],ns)>0
        for b,_,tt in out[:i]:
            if a==b:assert inner(st[0],tt[0],ns)==0
    return tuple(out)

@lru_cache(None)
def copies(word):
    ns=tuple(map(abs,word));L=len(ns);d={}
    for mask in range(1<<L):
        ix=tuple(i for i in range(L) if mask>>i&1)
        iy=tuple(i for i in range(L) if not mask>>i&1)
        parity=sum(word[i]<0 for i in iy)%2
        for a,pa,_ in cg(tuple(ns[i] for i in ix)):
            for b,pb,_ in cg(tuple(ns[i] for i in iy)):
                d.setdefault((a,b),[[],[]])[parity].append((mask,pa,pb))
    for value in d.values():
        for v in value:v.sort()
    return d

def surplus(word):
    out={}
    for ab,(e,o) in copies(word).items():
        r=min(len(e),len(o))
        if len(e)>r:out[ab]=(0,e[r:])
        if len(o)>r:out[ab]=(1,o[r:])
    return out

def tensors(ns,p,q):
    ans={}
    for a,c in p.items():
        for b,d in q.items():
            v=tuple(x+y for x,y in zip(a,b))
            ans[v]=ans.get(v,0)+c*d
    return clean(ans)
def embed(p,idx,L):
    out={}
    for v,c in p.items():
        w=[0]*L
        for i,h in zip(idx,v):w[i]=h
        out[tuple(w)]=c
    return out
def states_for(ns,path):
    return next(st for a,pa,st in cg(ns) if pa==path)
def invpair(ns,left,right,pl,pr,spin):
    L=len(ns);ls=states_for(tuple(ns[i] for i in left),pl)
    rs=states_for(tuple(ns[i] for i in right),pr);ans={}
    for h in range(spin+1):
        p=tensors(ns,embed(ls[h],left,L),embed(rs[spin-h],right,L))
        c=(-1)**h*comb(spin,h)
        for v,z in p.items():ans[v]=ans.get(v,0)+c*z
    return clean(ans)
def hvector(word,A,B,ab,u,v):
    ns=tuple(map(abs,word));L=len(ns);sx,pa,qa=u;tx,pb,qb=v
    ax=tuple(i for j,i in enumerate(A) if sx>>j&1)
    ay=tuple(i for j,i in enumerate(A) if not sx>>j&1)
    bx=tuple(i for j,i in enumerate(B) if tx>>j&1)
    by=tuple(i for j,i in enumerate(B) if not tx>>j&1)
    x=invpair(ns,ax,bx,pa,pb,ab[0])
    y=invpair(ns,ay,by,qa,qb,ab[1]);ans={}
    mask=sum(1<<i for i in ay+by)
    p=tensors(ns,x,y)
    assert p and not raise_op(p) and not lower(p,ns)
    return mask,p

def split(word):
    return tuple(range(0,len(word),2)),tuple(range(1,len(word),2))

def meta(w,par):
    A,B=split(w);aa=surplus(tuple(w[i] for i in A));bb=surplus(tuple(w[i] for i in B));out=[]
    for ab in sorted(aa.keys()&bb.keys()):
        ep,al=aa[ab];eq,bl=bb[ab]
        if ep^eq!=par:continue
        for u in al:
            for v in bl:
                S,p=hvector(w,A,B,ab,u,v);out.append((S,ab,p))
    return out
def matrix(w):
    log("CG matrix start",w)
    E,O=meta(w,0),meta(w,1);ns=tuple(map(abs,w))
    M=sp.Matrix([[inner(p,q,ns)*prod(h+2 for h in range(len(w)) if (S&T)>>h&1)
                  for S,ab,q in O] for T,cd,p in E])
    log("CG matrix built",M.shape)
    return E,O,M
def signs_of_edges(E,O,M,degree):
    return {(degree(E[i])-degree(O[j])>0)-(degree(E[i])-degree(O[j])<0)
            for i in range(len(E)) for j in range(len(O)) if M[i,j]}
cut=lambda v:v[0].bit_count()
height=lambda v:sum(v[1])
casimir=lambda v:sum(a*(a+2) for a in v[1])

w=(-1,-1,-1,-2,3,4);E,O,M=matrix(w)
assert M.shape==(22,2) and M.rank()==2
assert M[0,0]==630 and M[9,0]==5 and M[7,0]==756
assert signs_of_edges(E,O,M,cut)=={-1,1}
assert signs_of_edges(E,O,M,casimir)=={-1,1}
log("PASS six-factor cut and Casimir incompatibility")

w=(-2,)*3+(-4,)*3;E,O,M=matrix(w)
assert M.shape==(60,2) and M.rank()==2
assert M[0,0]==Q(-75,2) and M[36,0]==105
assert height(O[0])==6 and height(E[0])==2 and height(E[36])==8
assert signs_of_edges(E,O,M,height)=={-1,0,1}
log("PASS height has both directions already at six factors")

w=(1,)*5+(-2,3,-4);E,O,M=matrix(w)
assert M.shape==(176,26)
assert {height(v) for v in O}=={6}
assert signs_of_edges(E,O,M,height)=={-1,0}
E6=[i for i,v in enumerate(E) if height(v)==6]
E4=[i for i,v in enumerate(E) if height(v)==4]
rr,piv=M.extract(E6,range(26)).rref()
assert len(E6)==20 and len(piv)==20
free=[j for j in range(26) if j not in piv]
K=sp.zeros(26,6)
for col,j in enumerate(free):
    K[j,col]=1
    for r,p in enumerate(piv):K[p,col]=-rr[r,j]
assert M.extract(E6,range(26))*K==sp.zeros(20,6)
log("height E1 odd dimension",K.cols)
U=M.extract(E4,range(26))*K
chosen=[28,29,31,32,33,35]
minor=U.extract([E4.index(j) for j in chosen],range(6)).det(method="domain-ge")
assert minor==sp.Rational(-472768994315782089853519649227223484375,
                         69853013611007444)
log("higher differential minor",minor)
assert {h:sum(height(v)==h for v in E) for h in (0,2,4,6)}=={0:4,2:72,4:80,6:20}
bidir=set()
for a,b in combinations(range(8),2):
    flag=lambda v:((v[0]>>a) ^(v[0]>>b))&1
    ss=signs_of_edges(E,O,M,flag)
    if -1 in ss and 1 in ss:bidir.add((a,b))
assert bidir==set(combinations(range(8),2))-{(5,7)}
assert all((((v[0]>>5) ^(v[0]>>7))&1)==1 for v in O)
assert all((((v[0]>>5) ^(v[0]>>7))&1)==0 for v in E)
log("PASS: 27 pair-colour filtrations incompatible; remaining one has 26 odd E1 dimensions")
log("PASS: height E1=(156,6), E2=(150,0)")

def CG(a,b):return range(abs(a-b),a+b+1,2)
def insert(tab,a,eps=1):
    ans={}
    for (r,s),z in tab.items():
        for t in CG(r,a):ans[t,s]=ans.get((t,s),0)+z
        for t in CG(s,a):ans[r,t]=ans.get((r,t),0)+eps*z
    return {k:v for k,v in ans.items() if v}
def defect(tab,n,m):
    return sum(tab.get((c,0),0) for c in CG(n,m))-tab.get((n,m),0)
def tetra(tab,n,m,a):
    T=sum(tab.get((s,0),0) for c in CG(n,m) for s in CG(c,a))
    X=sum(tab.get((c,a),0) for c in CG(n,m))
    Y=sum(tab.get((c,m),0) for c in CG(n,a))
    Z=sum(tab.get((c,n),0) for c in CG(m,a))
    return T,X,Y,Z
R={(0,0):2,(2,0):2,(0,2):2,(1,1):2,(1,3):2,(3,1):2,(3,3):2}
assert all(v>=0 and R.get((s,r),0)==v for (r,s),v in R.items())
assert all(defect(R,n,m)>=0 for n in range(30) for m in range(30))
assert all(defect(R,n,m)>=0 for n,m in R if n and m)
assert tetra(R,1,3,2)==(6,0,4,4)
assert defect(R,1,3)==0 and defect(insert(R,2),1,3)==-2
assert insert(insert(insert(R,2),1,-1),3,-1).get((0,0),0)==-4
H=sp.Matrix([[6,0,4,4],[0,6,4,4],[4,4,6,0],[4,4,0,6]])
v=sp.Matrix([1,1,-1,-1])
assert H*v==-2*v

def addp(*pp):
    out={}
    for p in pp:
        for a,z in p.items():out[a]=out.get(a,0)+z
    return {a:z for a,z in out.items() if z}
def scal(p,c):return {a:c*z for a,z in p.items() if c*z}
def mulp(p,q):
    out={}
    for (i,j),x in p.items():
        for (k,l),y in q.items():out[i+k,j+l]=out.get((i+k,j+l),0)+x*y
    return {a:z for a,z in out.items() if z}
one={(0,0):1};xx={(1,0):1};yy={(0,1):1}
def U(n,x):
    p,q=one,x
    if n==0:return p
    for j in range(1,n):p,q=q,addp(mulp(x,q),scal(p,-1))
    return q
S2=addp(U(2,xx),U(2,yy))
RR=scal(addp(one,S2,mulp(addp(U(1,xx),U(3,xx)),
                         addp(U(1,yy),U(3,yy)))),2)
D1=addp(U(1,xx),scal(U(1,yy),-1))
D3=addp(U(3,xx),scal(U(3,yy),-1))
def moment(n):return comb(n,n//2)//(n//2+1) if n%2==0 else 0
def mean(p):return sum(z*moment(i)*moment(j) for (i,j),z in p.items())
assert mean(mulp(mulp(RR,D1),D3))==0
assert mean(mulp(mulp(mulp(RR,D1),D3),S2))==-4
log("PASS genuine-G countermodel: T,X,Y,Z=(6,0,4,4), parent 0, insertion -4")
checks=0
for L in range(7):
    for w in combinations_with_replacement(range(1,4),L):
        tab=char(w)
        for n,m,a in product(range(1,5),repeat=3):
            T,X,Y,Z=tetra(tab,n,m,a)
            assert defect(insert(tab,a),n,m)==T+X-Y-Z
            checks+=1
assert checks==5376
log("PASS",checks,"exact insertion identities")
log("PASS FM-STR9c: filtration obstructions, pages, census, and genuine-G countermodel")
