import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import prod
if any(x in ("-h","--help") for x in sys.argv[1:]):
    print("Exact FM-STR4 verifier for two fusion identities, one LP dual, and one Farkas vector.")
    raise SystemExit(0)
def cg(a,b): return range(abs(a-b),a+b+1,2)
def mcount(ns):
    d={0:1}
    for n in ns:
        q=defaultdict(int)
        for a,v in d.items():
            for b in cg(a,n): q[b]+=v
        d=dict(q)
    return d.get(0,0)
def coeff2(L,A,B):
    d={(0,0):1}
    for z in L:
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (x,y),v in d.items():
            for t in cg(x,n): q[t,y]+=v
            for t in cg(y,n): q[x,t]+=e*v
        d=dict(q)
    return d.get((A,B),0)
def fwht(a):
    a=list(a); h=1
    while h<len(a):
        for st in range(0,len(a),2*h):
            for j in range(st,st+h):
                x,y=a[j],a[j+h]; a[j],a[j+h]=x+y,x-y
        h*=2
    return a
@lru_cache(None)
def table(labels):
    n=len(labels); full=(1<<n)-1
    m=tuple(mcount([labels[i] for i in range(n) if (s>>i)&1]) for s in range(1<<n))
    f=tuple(m[s]*m[full^s] for s in range(1<<n))
    return m,f,tuple(fwht(f))
def gp(m,bmask,bminus,pbit=64):
    eps=[-1 if (bminus>>i)&1 else 1 for i in range(6)]
    z=0; sub=bmask
    while True:
        rest=bmask^sub
        z+=prod(eps[i] for i in range(6) if (rest>>i)&1)*m[sub|pbit]*m[rest]
        if sub==0: break
        sub=(sub-1)&bmask
    return z
def phi_child(labels,eps,i,j,c):
    outside=[k for k in range(len(labels)) if k not in (i,j)]
    labs=tuple(labels[k] for k in outside)+(c,)
    ce=[eps[k] for k in outside]+[eps[i]*eps[j]]
    mask=sum(1<<q for q,e in enumerate(ce) if e<0)
    return table(labs)[2][mask]
def get_rows(labels):
    m,f,phi=table(labels); pairs=list(combinations(range(7),2)); rows=[]
    for bm in range(64):
        sig=(-1)**bm.bit_count(); T=bm|(64 if sig<0 else 0); P=phi[T]
        ds=[]
        for i,j in pairs:
            z=P-phi[T^(1<<i)^(1<<j)]; assert z%4==0; ds.append(z//4)
        rows.append((bm,T,P,ds))
    return m,f,phi,pairs,rows
def check_identity(labels,top,cert,noflip):
    m,f,phi,pairs,rows=get_rows(labels); deltas=[]
    for bm,T,P,ds in rows:
        eps=[-1 if (bm>>i)&1 else 1 for i in range(6)]; eps.append(prod(eps))
        signed=[labels[i]*eps[i] for i in range(7)]
        g=gp(m,63,bm); child=gp(m,63^(1<<top[0])^(1<<top[1]),bm)
        assert coeff2(signed,0,0)==P
        assert coeff2(signed[:6],labels[6],0)==g
        assert coeff2([signed[k] for k in range(6) if k not in top],labels[6],0)==child
        delta=g-child; deltas.append(delta); rhs=Fraction(0)
        for i,j,c,beta in cert:
            v=phi_child(labels,eps,i,j,c)
            assert v>=0
            rhs+=beta*v
        assert rhs==delta,(bm,delta,rhs)
    assert [bm for bm,T,P,ds in rows if max(ds)<0]==noflip
    return m,f,phi,pairs,rows,deltas
L1=(1,2,5,9,11,12,14); TP1=(3,4)
C1=((0,1,3,Fraction(212,1331)),(0,4,10,Fraction(271,1331)),
(1,3,11,Fraction(907,2662)),(1,6,12,Fraction(907,2662)),
(3,4,2,Fraction(212,1331)),(3,6,5,Fraction(424,1331)),
(5,6,2,Fraction(153,2662)))
m1,f1,phi1,pairs1,rows1,d1=check_identity(L1,TP1,C1,[6,13,16,27,35,40,53,62])
assert (min(d1),max(d1),m1[-1],f1[0])==(313,329,325,325)
L2=(1,3,4,6,7,9,12); TP2=(4,5)
C2=((0,2,3,Fraction(180,463)),(0,3,7,Fraction(108,463)),
(1,2,1,Fraction(49,1389)),(2,4,3,Fraction(535,1852)),
(2,6,8,Fraction(391,926)),(4,6,5,Fraction(36,463)),
(5,6,3,Fraction(175,926)),(5,6,7,Fraction(175,926)))
m2,f2,phi2,pairs2,rows2,d2=check_identity(L2,TP2,C2,[7,12,18,25,33,42,52,63])
assert (min(d2),max(d2),m2[-1],f2[0])==(208,236,221,221)
assert (phi1[6],gp(m1,63,6),gp(m1,63^(1<<3)^(1<<4),6),d1[6])==(638,319,6,313)
assert [rows1[6][3][k] for k in (15,17,19)]==[-2,-6,-4]

# Verify Q_ij equals the sum of individual CG-channel children.
for labels,phi,rows in ((L1,phi1,rows1),(L2,phi2,rows2)):
    for bm,T,parent,ds in rows:
        eps=[-1 if (bm>>i)&1 else 1 for i in range(6)]; eps.append(prod(eps))
        for i,j in combinations(range(7),2):
            q=parent+phi[T^(1<<i)^(1<<j)]
            assert q%2==0
            assert q//2==sum(phi_child(labels,eps,i,j,c) for c in cg(labels[i],labels[j]))
print("FUSION_CERT_1",min(d1),max(d1),"noflip",[6,13,16,27,35,40,53,62])
print("FUSION_CERT_2",min(d2),max(d2),"noflip",[7,12,18,25,33,42,52,63])
print("BOUNDARY mask 6: Phi=638, g=319, child=6, Delta=313, three D values=(-2,-6,-4)")

# Exact LP dictionary: 21 D columns, every individual pair-channel child,
# 31 nonempty parity-preserving removals, and one constant f-cone column.
specs=[]
for i,j in pairs1:
    for c in cg(L1[i],L1[j]): specs.append(("F",i,j,c))
for rm in range(1,64):
    if sum(L1[i] for i in range(6) if (rm>>i)&1)%2==0: specs.append(("R",rm))
A=[]; b=[]
for bm,T,parent,ds in rows1:
    eps=[-1 if (bm>>i)&1 else 1 for i in range(6)]; eps.append(prod(eps)); vals=[]
    for spec in specs:
        if spec[0]=="F": vals.append(phi_child(L1,eps,*spec[1:]))
        else: vals.append(gp(m1,63^spec[1],bm))
    assert min(vals)>=0
    A.append(ds+[-v for v in vals]+[-1]); b.append(-d1[bm])
sel={("F",0,1,3):Fraction(212,1331),("F",0,4,10):Fraction(271,1331),
("F",1,3,11):Fraction(907,2662),("F",1,6,12):Fraction(907,2662),
("F",3,4,2):Fraction(212,1331),("F",3,6,5):Fraction(424,1331),
("F",5,6,2):Fraction(153,2662)}
x=[Fraction(0)]*(21+len(specs)+1)
for k,spec in enumerate(specs):
    if spec in sel: x[21+k]=sel[spec]
assert all(sum(A[r][j]*x[j] for j in range(len(x)))==b[r] for r in range(64))
obj=sum(x[:-1]); assert obj==Fraction(4205,2662) and sum(v>0 for v in x)==7
y={1:Fraction(49,10648),2:Fraction(1313,10648),5:Fraction(-7,121),
7:Fraction(-813,5324),9:Fraction(1663,10648),10:Fraction(1301,10648),
14:Fraction(-260,1331)}
for j in range(len(x)):
    yc=sum(y.get(r,Fraction(0))*A[r][j] for r in range(64))
    assert yc <= (0 if j==len(x)-1 else 1)
assert sum(y.get(r,Fraction(0))*b[r] for r in range(64))==obj
print("LP",len(specs),"child columns; exact L1 optimum",obj,"support 7")

# Farkas obstruction after replacing each pair's channels by aggregate Q.
rems=[("R2",i,j) for i,j in combinations(range(6),2) if L1[i]%2==L1[j]%2]
rems += [("R1",i,-1) for i in range(6) if L1[i]%2==0]
Agg=[]
for bm,T,parent,ds in rows1:
    vals=[(parent+phi1[T^(1<<i)^(1<<j)])//2 for i,j in pairs1]
    for typ,i,j in rems:
        rm=(1<<i)^(1<<j) if typ=="R2" else (1<<i)
        vals.append(gp(m1,63^rm,bm))
    Agg.append(ds+[-v for v in vals]+[-1])
yf={0:Fraction(321,652),6:Fraction(-1,2),13:Fraction(-5,652)}
cols=[sum(yf.get(r,Fraction(0))*Agg[r][j] for r in range(64))
      for j in range(len(Agg[0]))]
yb=sum(yf.get(r,Fraction(0))*(-d1[r]) for r in range(64))
assert len(cols)==52 and min(cols)==0 and all(v>=0 for v in cols)
assert yb==Fraction(-1003,326)
print("FARKAS",len(cols),"columns; min yA",min(cols),"yb",yb)
print("ALL CHECKS PASS")
