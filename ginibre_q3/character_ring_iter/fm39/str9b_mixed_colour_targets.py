import sys,time
if any(x in ("-h","--help") for x in sys.argv[1:]):
    print("Usage: python3 -u -B verifier.py; symbolic mixed-colour family, exact dimensions, and modular rank certificates.")
    raise SystemExit
from datetime import datetime, timezone
from fractions import Fraction as QQ
from math import comb,prod
from functools import lru_cache
from itertools import product
from sympy.polys.fields import field
from sympy.polys.domains import QQ as SQ
import sympy as sp
F,n=field("n",SQ)
def say(*a): print(datetime.now(timezone.utc).isoformat(timespec="seconds"),*a,flush=True)
def ff(x): return F(x.numerator)/x.denominator if isinstance(x,QQ) else F(x)
def add(q,key,x):
    if x:q[key]=q.get(key,F.zero)+x
def cl(p): return {v:c for v,c in p.items() if c}
def low(p,ns):
    q={}
    for v,c in p.items():
        for i,a in enumerate(ns):
            w=list(v);w[i]+=1
            add(q,tuple(w),c*(a-v[i]))
    return cl(q)
@lru_cache(None)
def paths(ns):
    if not ns:return ((0,()),)
    return tuple((c,pa+(c,)) for a,pa in paths(ns[:-1])
                 for c in range(abs(a-ns[-1]),a+ns[-1]+1,2))
@lru_cache(None)
def cgs(ns,pa):
    if not ns:return ({():F.one},)
    n0=ns[-1];c=pa[-1];old=cgs(ns[:-1],pa[:-1]);a=pa[-2] if len(pa)>1 else 0
    j=(a+n0-c)//2
    hv={}
    for h in range(j+1):
        for v,z in old[h].items():add(hv,v+(j-h,),(-1)**h*comb(j,h)*z)
    st=[cl(hv)]
    for h in range(c):st.append({v:z/(c-h) for v,z in low(st[-1],ns).items()})
    return tuple(st)
def surplus(word):
    out={}
    for mask in range(1<<len(word)):
        xx=tuple(abs(word[i]) for i in range(len(word)) if mask>>i&1)
        yy=tuple(abs(word[i]) for i in range(len(word)) if not mask>>i&1)
        par=sum(word[i]<0 for i in range(len(word)) if not mask>>i&1)%2
        for a,pa in paths(xx):
            for b,pb in paths(yy):out.setdefault((a,b),[[],[]])[par].append((mask,pa,pb))
    ans={}
    for ab,eo in out.items():
        for z in eo:z.sort()
        d=len(eo[0])-len(eo[1])
        if d:ans[ab]=(int(d<0),eo[int(d<0)][min(map(len,eo)):])
    return ans
def metadata(n0,par):
    w=(1,)*5+(-2,n0,-n0-1);aa=surplus(w[::2]);bb=surplus(w[1::2]);ans=[]
    for ab in sorted(aa.keys()&bb.keys()):
        ep,al=aa[ab];eq,bl=bb[ab]
        if ep^eq !=par:continue
        for u in al:
            for v in bl:ans.append((ab,u,v))
    return ans
small=(1,1,1,1,1,2)
def emb(p,idx):
    ans={}
    for v,c in p.items():
        w=[0]*6
        for i,h in zip(idx,v):w[i]=h
        add(ans,tuple(w),c)
    return cl(ans)
def mul(p,q):
    ans={}
    for v,c in p.items():
        for w,d in q.items():add(ans,tuple(a+b for a,b in zip(v,w)),c*d)
    return cl(ans)
def cup_small(li,ri,lp,rp,s):
    ls=cgs(tuple(small[i] for i in li),lp)
    rs=cgs(tuple(small[i] for i in ri),rp);ans={}
    for h in range(s+1):
        for v,c in mul(emb(ls[h],li),emb(rs[s-h],ri)).items():
            add(ans,v,(-1)**h*comb(s,h)*c)
    return cl(ans)
def cup_big(li,ri,lp,rp,out0,n0):
    assert li[-1]==6 and ri[-1]==7
    li,ri=li[:-1],ri[:-1];sl=lp[-2] if len(lp)>1 else 0;sr=rp[-2] if len(rp)>1 else 0
    off=out0-n0;out=n+off;jl=(sl-off)//2;jr=(sr+1-off)//2
    ls=cgs(tuple(small[i] for i in li),lp[:-1])
    rs=cgs(tuple(small[i] for i in ri),rp[:-1])
    hv={}
    for h in range(jr+1):
        for v,c in rs[h].items():add(hv,v+(jr-h,),(-1)**h*comb(jr,h)*c)
    st=[cl(hv)]
    for h in range(sl-jl):
        st.append({v:z/(out-h) for v,z in low(st[-1],tuple(small[i] for i in ri)+(n+1,)).items()})
    kr=(sum(small[i] for i in ri)+1-off)//2
    ans={}
    for h in range(sl-jl+1):
        rr={}
        for v,c in st[h].items():
            smallrev=tuple(small[i]-v[j] for j,i in enumerate(ri))
            add(rr,smallrev,c)
        for v,c in mul(emb(ls[jl+h],li),emb(rr,ri)).items():
            add(ans,v,(-1)**(jl+h+kr)*comb(sl-jl,h)*c)
    return cl(ans)
def vector(md,n0=8):
    ab,u,v=md;sx,pa,qa=u;tx,pb,qb=v
    A=(0,2,4,6);B=(1,3,5,7)
    ax=tuple(i for j,i in enumerate(A) if sx>>j&1)
    ay=tuple(i for j,i in enumerate(A) if not sx>>j&1)
    bx=tuple(i for j,i in enumerate(B) if tx>>j&1)
    by=tuple(i for j,i in enumerate(B) if not tx>>j&1)
    S=sum(1<<i for i in ay+by)
    assert (6 in ax)==(7 in bx)
    if 6 in ax:
        p=cup_big(ax,bx,pa,pb,ab[0],n0);q=cup_small(ay,by,qa,qb,ab[1])
    else:
        p=cup_small(ax,bx,pa,pb,ab[0]);q=cup_big(ay,by,qa,qb,ab[1],n0)
    ans=mul(p,q)
    assert all(0<=sum(z)-3<=4 for z in ans)
    return S,ans
def bn(x,h):return prod(x-i for i in range(h))/sp.factorial(h)
metric={}
for z in product(range(2),range(2),range(2),range(2),range(2),range(3)):
    ell=sum(z)-3
    if 0<=ell<=4:metric[z]=(n+1)/(prod(comb(a,b) for a,b in zip(small,z))*bn(n+1,ell))
def ip(p,q):
    return sum((c*q.get(z,0)*metric[z] for z,c in p.items()),F.zero)
def row(ev,od):
    T,p=ev
    return [ip(p,q)*prod(i+2 for i in range(8) if (S&T)>>i&1) for S,q in od]

def canonical(md,n0):
    ab,u,v=md
    def one(w):
        mask,px,py=w;px=list(px);py=list(py)
        if mask&8:px[-1]-=n0
        else:py[-1]-=n0
        return mask,tuple(px),tuple(py)
    aa=list(ab);aa[0 if u[0]&8 else 1]-=n0
    return tuple(aa),one(u),one(v)
reference=[set(canonical(md,8) for md in metadata(8,p)) for p in (0,1)]
for n0 in range(4,13):
    for p in (0,1):
        assert set(canonical(md,n0) for md in metadata(n0,p))==reference[p]
say("matching stability n=4..12 PASS")

say("start symbolic compressed CG")
od=[vector(md) for md in metadata(8,1)]
ev=[vector(md) for md in metadata(8,0)]
assert (len(ev),len(od))==(162,22)
pure=[v for v in ev if v[0] in (0,255)]
say("vectors",len(ev),len(od),"pure",len(pure))
R=[row(z,od) for z in pure]
ks=[
{1:QQ(-17,12),4:QQ(-29,12),7:QQ(187,12),9:QQ(-34,9),10:QQ(221,61),11:QQ(153,244),12:QQ(-17,36),14:-17,17:1},
{1:QQ(29,18),4:QQ(47,18),7:QQ(-235,18),9:QQ(235,27),10:QQ(-1363,183),11:QQ(-47,122),12:QQ(-47,54),14:QQ(47,3),20:1}]
K=[[ff(d.get(i,0)) for i in range(22)] for d in ks]
assert all(not sum((a*b for a,b in zip(r,k)),F.zero) for r in R for k in K)
say("symbolic pure nullspace inclusion PASS")
H=[[sum((a*b for a,b in zip(row(ev[h],od),k)),F.zero) for k in K] for h in (12,13)]
assert H==[[-ff(QQ(506872,549))*(n+2),ff(QQ(1928128,1647))*(n+2)],[-1080*(n+2),940*(n+2)]]
say("mixed 2x2 formula PASS")
det=H[0][0]*H[1][1]-H[0][1]*H[1][0]
assert det==ff(QQ(217666400,549))*(n+2)**2
say("mixed determinant = 217666400*(n+2)^2/549")

norms=[ip(q,q) for S,q in od]
G0=[[QQ(8092446253,2411208),QQ(-11244780959,3616812)],
    [QQ(-11244780959,3616812),QQ(16744235653,5425218)]]
for a in range(2):
    for b in range(2):
        assert sum((K[a][i]*K[b][i]*norms[i] for i in range(22)),F.zero)==(n+2)*ff(G0[a][b])
assert ip(ev[12][1],ev[12][1])==2*(n+1)*(n+2)**2/(3*n*(n-1))
assert ip(ev[13][1],ev[13][1])==4*(n+1)*(n+2)/n
C0=[[QQ(-506872,549),QQ(1928128,1647)],[-1080,940]]
c0=QQ(217666400,549)**2/(5*sum(x*x for rr in C0 for x in rr)*(G0[0][0]+G0[1][1]))
assert c0==QQ(23133561138977057395200,20256219355112208219797) and c0>1
say("uniform mixed energy constant",c0,"> 1 PASS")

# symbolic elimination on pure matrix; certificates tracked as positive factors later
rows=[r[:] for r in R];rank=0;piv=[];pcols=[]
for col in range(22):
    h=next((h for h in range(rank,len(rows)) if rows[h][col]),None)
    if h is None:continue
    rows[rank],rows[h]=rows[h],rows[rank]
    pv=rows[rank][col];piv.append(pv);pcols.append(col)
    rows[rank]=[v/pv for v in rows[rank]]
    for ii in range(rank+1,len(rows)):
        c=rows[ii][col]
        if c:rows[ii]=[x-c*y for x,y in zip(rows[ii],rows[rank])]
    rank+=1

assert rank==20
constants=[-3240,864,216,1080,576,288,-1080,288,72,360,QQ(2,5),
 QQ(305,36),QQ(6,5),QQ(33,2),QQ(-1,9),QQ(-83,360),QQ(27,8),
 QQ(-353,1440),QQ(-17,18),QQ(803,3240)]
for h,(pv,c) in enumerate(zip(piv,constants),1):
    assert pv==ff(c)*(n+2)*((n+3)/(n+1) if h in (1,4,7,16,18,20) else 1)
assert pcols==[h for h in range(22) if h not in (17,20)]
say("all 20 nonzero pivot formulas PASS; pure rank = 20")
say("PASS symbolic family")

def conv(tab,label,sign):
    z={}
    for (a,b),v in tab.items():
        for c in range(abs(a-label),a+label+1,2):z[c,b]=z.get((c,b),0)+v
        for c in range(abs(b-label),b+label+1,2):z[a,c]=z.get((a,c),0)+sign*v
    return {k:v for k,v in z.items() if v}

for n0 in range(4,13):
    w=(1,)*5+(-2,n0,-n0-1)
    for ww in (w,tuple(-v if abs(v)%2 else v for v in w)):
        z={(0,0):1}
        for v in ww:z=conv(z,abs(v),1 if v>0 else -1)
        assert z.get((0,0),0)==140
say("18 direct Phi and parity-reflection checks PASS")

def spin0(ns):
    z={0:1}
    for label in ns:
        zz={}
        for a,v in z.items():
            for c in range(abs(a-label),a+label+1,2):zz[c]=zz.get(c,0)+v
        z=zz
    return z.get(0,0)
def cb(t,j):return comb(t,j) if 0<=j<=t else 0
def exact_G(t,a,b):
    if (t+a+b)%2:return 0
    num=(a+1)*(b+1)*cb(t+2,(t+a+b)//2+2)*cb(t+2,(t+a-b)//2+1)
    den=(t+1)*(t+2);assert num%den==0
    return num//den
G={(0,0):1};checked=0
for t in range(33):
    if t<=12:
        for a in range(t+2):
            for b in range(t+2):
                assert G.get((a,b),0)==exact_G(t,a,b);checked+=1
    if t in (4,8,16,32):
        FA=conv(conv(G,2,-1),6,-1);FB=conv(G,4,1)
        vals=[v*FB.get(k,0) for k,v in FA.items()]
        he=sum(max(z,0) for z in vals);ho=sum(max(-z,0) for z in vals)
        cap=2*spin0((1,)*(2*t)+(2,4,6))
        interval={4:(6,7),8:(13,14),16:(19,20),32:(23,24)}[t]
        assert interval[0]*he<100*ho<interval[1]*he
        say("GROWTH t=",t,"H0=",(he,ho),"pure capacity <=",cap,
            "100*odd/even in",interval)
    G=conv(G,1,1)
xx,yy=sp.symbols("xx yy",real=True)
hh=xx*yy*sp.exp(-(xx**2+yy**2)/2)
LL=lambda h:sp.diff(h,xx,2)-sp.diff(h,yy,2)
assert sp.simplify(LL(LL(hh))/hh-
 ((xx**2-yy**2)**2-4*(xx**2+yy**2)+12))==0
say("binomial identities",checked,"and limiting polynomial PASS")
# Store symbolic vectors reduced at n=4 for an independent full-CG comparison.
prime=65521
def spec(v):
    z=v.numer.evaluate(0,4)/v.denom.evaluate(0,4)
    return int(z.numerator)*pow(int(z.denominator),-1,prime)%prime
symbolic_signatures=[sorted((S,tuple(sorted((z,spec(c)) for z,c in q.items() if spec(c))))
                           for S,q in vv) for vv in (ev,od)]

import sys,time
if any(x in ("-h","--help") for x in sys.argv[1:]):
    print("Usage: python3 -u -B verifier.py; exact finite-field exchange ranks.")
    raise SystemExit
from functools import lru_cache
from itertools import product
from math import comb,prod,isqrt
from datetime import datetime,timezone
import numpy as np
P=65521
assert all(P%d for d in range(2,isqrt(P)+1))
def log(*args):print(datetime.now(timezone.utc).isoformat(),*args,flush=True)
def clean(x):return {v:c%P for v,c in x.items() if c%P}
def lower(p,ns):
    q={}
    for v,c in p.items():
        for i,n in enumerate(ns):
            if v[i]<n:
                w=list(v);w[i]+=1;w=tuple(w);q[w]=(q.get(w,0)+c*(n-v[i]))%P
    return clean(q)
@lru_cache(None)
def cg(ns):
    if not ns:return ((0,(),({():1},)),)
    old,n=ns[:-1],ns[-1];out=[]
    for a,path,states in cg(old):
        for c in range(abs(a-n),a+n+1,2):
            j=(a+n-c)//2;hv={}
            for h in range(j+1):
                for v,z in states[h].items():
                    w=v+(j-h,);hv[w]=(hv.get(w,0)+(-1)**h*comb(j,h)*z)%P
            st=[clean(hv)]
            for h in range(c):
                den=pow(c-h,-1,P);st.append({v:z*den%P for v,z in lower(st[-1],ns).items()})
            assert not lower(st[-1],ns)
            out.append((c,path+(c,),tuple(st)))
    return tuple(out)
@lru_cache(None)
def states(ns,path):return next(st for a,pa,st in cg(ns) if pa==path)
@lru_cache(None)
def copies(w):
    ns=tuple(map(abs,w));out={}
    for mask in range(1<<len(w)):
        ix=tuple(i for i in range(len(w)) if mask>>i&1);iy=tuple(i for i in range(len(w)) if not mask>>i&1);p=sum(w[i]<0 for i in iy)%2
        for a,pa,_ in cg(tuple(ns[i] for i in ix)):
            for b,pb,_ in cg(tuple(ns[i] for i in iy)):
                out.setdefault((a,b),[[],[]])[p].append((mask,pa,pb))
    for ev,od in out.values():ev.sort();od.sort()
    return out
def surplus(w):
    d={}
    for ab,(ev,od) in copies(w).items():
        n=min(len(ev),len(od))
        if len(ev)>n:d[ab]=(0,ev[n:])
        if len(od)>n:d[ab]=(1,od[n:])
    return d
def embed(p,ix,L):
    d={}
    for v,c in p.items():
        w=[0]*L
        for i,h in zip(ix,v):w[i]=h
        d[tuple(w)]=c
    return d
def mul(p,q):
    d={}
    for v,c in p.items():
        for w,z in q.items():
            t=tuple(a+b for a,b in zip(v,w));d[t]=(d.get(t,0)+c*z)%P
    return clean(d)
@lru_cache(None)
def cup(ns,left,right,pl,pr,spin):
    ls=states(tuple(ns[i] for i in left),pl);rs=states(tuple(ns[i] for i in right),pr);d={}
    for h in range(spin+1):
        z=mul(embed(ls[h],left,len(ns)),embed(rs[spin-h],right,len(ns)));sg=(-1)**h*comb(spin,h)
        for v,c in z.items():d[v]=(d.get(v,0)+sg*c)%P
    return clean(d)
def vectors(w,parity,coord):
    ns=tuple(map(abs,w));L=len(w);A=tuple(range(0,L,2));B=tuple(range(1,L,2));aa=surplus(tuple(w[i] for i in A));bb=surplus(tuple(w[i] for i in B));rows=[];masks=[];channels=[]
    for ab in sorted(aa.keys()&bb.keys()):
        ep,al=aa[ab];eq,bl=bb[ab]
        if ep^eq!=parity:continue
        for sx,pa,qa in al:
            ax=tuple(i for j,i in enumerate(A) if sx>>j&1);ay=tuple(i for j,i in enumerate(A) if not sx>>j&1)
            for tx,pb,qb in bl:
                bx=tuple(i for j,i in enumerate(B) if tx>>j&1);by=tuple(i for j,i in enumerate(B) if not tx>>j&1)
                z=mul(cup(ns,ax,bx,pa,pb,ab[0]),cup(ns,ay,by,qa,qb,ab[1]))
                row=np.zeros(len(coord),dtype=np.int64)
                for v,c in z.items():row[coord[v]]=c
                rows.append(row);masks.append(sum(1<<i for i in ay+by));channels.append(ab)
                if len(rows)%200==0:log("vectors",parity,"completed",len(rows))
    return np.array(rows,dtype=np.int64),masks,channels
def rank(a):
    a=a.copy();r=0
    for j in range(a.shape[1]):
        v=np.flatnonzero(a[r:,j])
        if not len(v):continue
        i=r+int(v[0]);a[[i,r]]=a[[r,i]];a[r,:]=a[r,:]*pow(int(a[r,j]),-1,P)%P
        for first in range(r+1,len(a),256):
            end=min(first+256,len(a));a[first:end,:]=(a[first:end,:]-a[first:end,j:j+1]*a[r:r+1,:])%P
        r+=1
        if r%100==0:log("rank pivots",r,"columns",a.shape[1])
        if r==min(a.shape):break
    return r
for w in ((1,)*5+(-2,4,-5),(1,)*7+(-2,3,-4)):
    start=time.monotonic();ns=tuple(map(abs,w));deg=sum(ns);assert deg<P and deg%2==0
    mon=[v for v in product(*(range(n+1) for n in ns)) if sum(v)*2==deg];coord={v:i for i,v in enumerate(mon)}
    assert len(mon)*(P-1)**2<2**63
    norm=np.array([pow(prod(comb(n,h) for n,h in zip(ns,v)),-1,P) for v in mon],dtype=np.int64)
    log("start",w,"monomials",len(mon))
    od,sm,sa=vectors(w,1,coord);ev,tm,ta=vectors(w,0,coord)
    log("homology dimensions",len(ev),len(od))
    if len(w)==8:
        for pp,(aa,mm) in enumerate(((ev,tm),(od,sm))):
            sig=sorted((S,tuple(sorted((v[:6],int(rr[h])) for h,v in enumerate(mon)
                                      if v[6]==0 and rr[h]))) for rr,S in zip(aa,mm))
            assert sig==symbolic_signatures[pp]
        log("all 184 projected vectors match independent full CG modulo",P)
    m=np.zeros((len(ev),len(od)),dtype=np.int64)
    for first in range(0,len(ev),128):
        end=min(first+128,len(ev));z=((ev[first:end,:]*norm)%P)@od.T%P
        for i in range(first,end):
            weights=np.array([prod(j+2 for j in range(len(w)) if (tm[i]&s)>>j&1)%P for s in sm],dtype=np.int64)
            m[i,:]=z[i-first,:]*weights%P
        log("exchange rows",end,"of",len(ev))
    full=(1<<len(w))-1;pure=[i for i,t in enumerate(tm) if t in (0,full)]
    rp=rank(m[pure,:]);rf=rank(m)
    expected=(162,22,47,20) if len(w)==8 else (1678,336,150,150)
    assert (len(ev),len(od),len(pure),rp)==expected and rf==len(od)
    log("RESULT",w,"H0",(len(ev),len(od)),"pure dimension",len(pure),"pure rank mod",rp,"full rank",rf,"seconds",time.monotonic()-start)

log("PASS: symbolic family, growth identities, and full-rank certificates")
