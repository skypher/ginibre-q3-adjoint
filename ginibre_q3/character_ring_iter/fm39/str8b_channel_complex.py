import sys
if any(a in ('-h','--help') for a in sys.argv[1:]):
    print('Usage: python3 -u -B verifier.py; exact channel complex and cut-exchange completion, labels 1..3, at most 8 factors.')
    raise SystemExit
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb,prod
from fractions import Fraction as Q
import sympy as sp

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
    y=invpair(ns,ay,by,qa,qb,ab[1])
    # S denotes the y-colour subset.
    mask=sum(1<<i for i in ay+by)
    p=tensors(ns,x,y)
    assert p and not raise_op(p) and not lower(p,ns)
    return mask,p

def split(word):
    return tuple(range(0,len(word),2)),tuple(range(1,len(word),2))
def hdata(word):
    A,B=split(word)
    a=copies(tuple(word[i] for i in A))
    b=copies(tuple(word[i] for i in B))
    he=ho=ro=ri=0
    for ab in a.keys()&b.keys():
        ea,oa=map(len,a[ab]);eb,ob=map(len,b[ab])
        ra,rb=min(ea,oa),min(eb,ob)
        ha,hb=ea-oa,eb-ob
        he+=max(ha*hb,0);ho+=max(-ha*hb,0)
        ro+=ra*rb+ra*max(hb,0)+max(ha,0)*rb
        ri+=ra*rb+ra*max(-hb,0)+max(-ha,0)*rb
    return he,ho,ro,ri

def hvectors(word,parity,selected=None):
    A,B=split(word)
    aa=surplus(tuple(word[i] for i in A))
    bb=surplus(tuple(word[i] for i in B))
    out=[];index=0
    for ab in sorted(aa.keys()&bb.keys()):
        ep,al=aa[ab];eq,bl=bb[ab]
        if (ep^eq)!=parity:continue
        for u in al:
            for v in bl:
                if selected is None or index in selected:
                    out.append(hvector(word,A,B,ab,u,v))
                index+=1
    return out,index

def phi(word):
    d={(0,0):1}
    for signed in word:
        n=abs(signed);sg=1 if signed>0 else -1;q={}
        for (a,b),v in d.items():
            for c in range(abs(a-n),a+n+1,2):q[c,b]=q.get((c,b),0)+v
            for c in range(abs(b-n),b+n+1,2):q[a,c]=q.get((a,c),0)+sg*v
        d=q
    return d.get((0,0),0)

def koszul_check(m):
    d=sp.zeros(2*m)
    for j in range(m):d[j,m+j]=1
    gamma=sp.diag(*([1]*m+[-1]*m))
    D=sp.kronecker_product(d,sp.eye(2*m))+sp.kronecker_product(gamma,d)
    assert D*D==sp.zeros(4*m*m)
    ev=[];od=[]
    for i in range(2*m):
        for j in range(2*m):
            (od if (i>=m)^(j>=m) else ev).append(i*2*m+j)
    assert D.extract(ev,od).rank()==m*m
    assert D.extract(od,ev).rank()==m*m

# Row selections only shorten the rank verification; they do not define d.
minors={
 (1,1,1,2,2,2,3):[0,1,3,4],
 (1,2,2,2,3,3,3):[2,3],
 (1,1,1,1,1,2,2,3):[0,2,7,11],
 (1,1,1,2,2,3,3,3):[1,2],
 (1,2,2,3,3,3,3,3):[0,1,9,10,93,94,95,102,103,104],
}

koszul_check(1);koszul_check(3)
assert hdata((-2,)*6)==(100,0,10,10)
assert phi((-2,)*6)==100
print('six adjoints: outgoing/incoming ranks 10,10; homology 100,0',flush=True)
total=nonzero=patched=0
for L in range(9):
    nc=nh=0
    for ns in combinations_with_replacement(range(1,4),L):
        for signs in range(8):
            if any(signs>>(n-1)&1 for n in (1,2,3) if n not in ns):continue
            word=tuple(-n if signs>>(n-1)&1 else n for n in ns)
            if sum(n<0 for n in word)%2:continue
            he,ho,ro,ri=hdata(word)
            assert he-ho==phi(word)
            total+=1;nc+=1
            if ro+ri+ho:nonzero+=1
            if not ho:continue
            rows=minors[ns]
            ev,ne=hvectors(word,0,set(rows))
            od,no=hvectors(word,1)
            assert (ne,no)==(he,ho) and len(ev)==len(od)==ho
            M=sp.Matrix([
                [inner(p,q,ns)*prod(i+2 for i in range(L) if (S&T)>>i&1)
                 for S,q in od] for T,p in ev])
            assert M.det()!=0
            print('CORRECTION',word,'H0=',(he,ho),
                  'rank=',ho,'H=',(he-ho,0),flush=True)
            patched+=1;nh+=1
    print('LEVEL',L,'profiles=',nc,'corrections=',nh,flush=True)
assert (total,nonzero,patched)==(481,88,10)
print('PASS: 481 profiles, 88 nonzero odd chain spaces, 10 corrections',flush=True)
