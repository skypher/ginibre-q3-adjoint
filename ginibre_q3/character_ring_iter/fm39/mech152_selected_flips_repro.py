import argparse,re
from collections import Counter,defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,prod
from pathlib import Path

ROOT=Path("/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc")
parser=argparse.ArgumentParser(description="FM-MECH152: exact flip regions and averaging obstructions.")
parser.add_argument("--census",type=Path,nargs="*",default=[ROOT/"fx3_40.log",ROOT/"fx3b_44.log"])
args=parser.parse_args()

@lru_cache(None)
def cg(a,b):return tuple(range(abs(a-b),a+b+1,2))

def entry(L,a=0,b=0):
    L=sorted(L,key=abs,reverse=True)
    rem=sum(map(abs,L));F={(0,0):1}
    for z in L:
        n=abs(z);sg=1 if z>0 else -1;rem-=n;G=defaultdict(int)
        for (i,j),v in F.items():
            for k in cg(i,n):
                if abs(k-a)+abs(j-b)<=rem:G[k,j]+=v
            for k in cg(j,n):
                if abs(i-a)+abs(k-b)<=rem:G[i,k]+=sg*v
        F={q:v for q,v in G.items() if v}
    return F.get((a,b),0)

def gain(L,u,v):
    C=list(L);C.remove(u);C.remove(v)
    return (1 if v>0 else -1)*entry(C,abs(u),abs(v))

def flipped(L,u,v):
    C=list(L);C.remove(u);C.remove(v)
    return C+[-u,-v]

def params(C):
    t=sum(z<0 for z in C);assert t%2==0
    lam=sum((Q((abs(z)-1)*(abs(z)+3),5) if z<0
             else Q(z*(z+2),3)) for z in C)
    return t//2,lam

def criterion(C,a):
    if a<2 or (a+sum(map(abs,C)))%2:return False
    r,lam=params(C);h=a+2*r+4;kap=(2*r+3)*lam
    return (a*a-1)*(h-kap)>=(2*r+1)*(2*r+3)*h

def mul(A,B):
    C=defaultdict(int)
    for (i,j),v in A.items():
        for (k,l),w in B.items():C[i+k,j+l]+=v*w
    return {q:v for q,v in C.items() if v}

@lru_cache(None)
def V(n):
    if n==0:return {(0,0):1}
    if n==1:return {(1,0):1,(0,1):1}
    P=mul({(1,0):1,(0,1):1},V(n-1))
    for q,v in V(n-2).items():P[q]=P.get(q,0)-4*v
    return {q:v for q,v in P.items() if v}

def core_poly(C):
    P={(0,0):1}
    for z in C:
        sg=1 if z>0 else -1
        P=mul(P,{(i,j):v*(1+sg*(-1)**j) for (i,j),v in V(abs(z)).items()
                 if 1+sg*(-1)**j})
    return P

@lru_cache(None)
def J(m,r):
    h=Q(comb(2*m,m)*comb(2*r,r),comb(m+r,m))
    ans=h*h*Q(2*(2*m+1)*(2*r+1),(m+r+1)**2*(m+r+2))
    assert ans.denominator==1
    return ans.numerator

def moment_value(C,a,extra_d=0):
    P=core_poly(tuple(C));ans=0
    for (i,j),v in P.items():
        if (a+i)%2 or (j+extra_d)%2:continue
        ans+=v*J((a+i)//2,(j+extra_d)//2)
    return Q(ans,2**sum(map(abs,C)))

def cat(n):return comb(2*n,n)//(n+1)
def Jdirect(m,r):
    n=2*(m+r);ans=0
    for h in range(0,n+1,2):
        c=sum((-1)**j*comb(2*r,j)*comb(2*m,h-j)
              for j in range(max(0,h-2*m),min(2*r,h)+1))
        ans+=c*cat(h//2)*cat((n-h)//2)
    return ans

for m in range(9):
    for r in range(9):
        assert J(m,r)==Jdirect(m,r)
        if m:
            assert Q(J(m-1,r+1),J(m,r))==Q((2*r+1)*(2*r+3),(2*m-1)*(2*m+1))
print("81 direct moment identities: PASS")

cores=[(-2,-2),(-3,-3),(-2,-3,4),(-3,-5,4),
       (-3,-3,4,6),(-2,-3,-4,-5,7)]
small=0
for C in cores:
    for a in range(2,10):
        if (a+sum(map(abs,C)))%2:continue
        assert moment_value(C,a)==entry([1]*a+list(C))
        assert moment_value(C,a-2,2)==entry([1]*(a-2)+[-1,-1]+list(C))
        small+=1
print("small independent fusion/moment pairs:",small)
for C in cores:
    a=2+sum(map(abs,C))%2
    while not criterion(C,a):a+=2
    r,lam=params(C)
    dim=prod(comb(abs(z)+2,3) if z<0 else 2*(z+1) for z in C)
    for aa in (a,a+2,a+10):
        difference=moment_value(C,aa)-moment_value(C,aa-2,2)
        assert difference>=0
        if aa%2==0:
            bound=dim*J(aa//2,r)*(1-Q((2*r+3)*lam,aa+2*r+4)
                                            -Q((2*r+1)*(2*r+3),aa*aa-1))
            assert difference>=bound>=0
    print("fundamental criterion:",C,"first admissible a =",a)

def psi(i,j):
    if i==j:return {}
    return {(i,j):1,(j,i):-1}
def plus(P,n):
    F=defaultdict(int)
    for (i,j),v in P.items():
        for k in cg(i,n):F[k,j]+=v
        for k in cg(j,n):F[i,k]+=v
    return {q:v for q,v in F.items() if v}
identities=0
for j in range(8):
    for i in range(j+2,j+11,2):
        P=plus(psi(i,j),2)
        P[(i,j)]=P.get((i,j),0)-1
        P[(j,i)]=P.get((j,i),0)+1
        R=defaultdict(int)
        terms=[(i+2,j),(i-2,j),(i,j+2)]
        if j>=2:terms.append((i,j-2))
        if j>=1:terms.append((i,j))
        for u,v in terms:
            for ij,c in psi(u,v).items():R[ij]+=c
        assert {q:v for q,v in P.items() if v}=={q:v for q,v in R.items() if v}
        assert all(v>=0 for (u,w),v in P.items() if u>w)
        identities+=1
print("ordered-cone S2 identities:",identities)

templates=[(8,6,2,(1,1,3,3)),(9,7,3,(1,4)),
           (12,10,4,(3,3,6)),(10,9,5,(3,3,4))]
checks=0
for P,Qn,q,H in templates:
    assert P>=sum(H) and (P-sum(H))%2==0
    assert q<Qn and (q-Qn)%2==0
    for b in range(9):
        L=[-P,-Qn,q]+list(H)+[2]*b
        v=gain(L,q,-Qn)
        assert v>=0
        assert entry(L)-entry(flipped(L,q,-Qn))==4*v
        checks+=1
print("unbounded-S2 theorem sample identities:",checks)

small_witness=[1,1,-2,3,3,-4]
assert entry(small_witness)==20
sv=[gain(small_witness,1,z) for z in small_witness[1:]]
assert sv==[-4,1,1,1,0] and sum(sv)==-1
print("uniform-star obstruction:",sv)

L=[1]*2+[-3]*12+[10]*3+[12]*7
assert entry(L)==10440907700952478
expected={1:-1116871498387752264,-3:1707100719943765,
          10:655179726140625,12:530410045330396}
for z,v in expected.items():
    assert gain(L,1,z)==v
    assert entry(L)-entry(flipped(L,1,z))==4*v
cnt=Counter(L)
label=sum((cnt[z]-(z==1))*abs(z)*v for z,v in expected.items())
casimir=sum((cnt[z]-(z==1))*abs(z)*(abs(z)+2)*v for z,v in expected.items())
assert label==-991206036877804710
assert casimir==-2183709450854208396
print("large-star obstructions:",label,casimir)
print("successful (1,-3) flip difference:",4*expected[-3])

def residual(L,p_signed):
    B=list(L);B.remove(p_signed);p=abs(p_signed);W=sum(map(abs,B))
    signs={}
    for z in L:
        if abs(z) in signs and signs[abs(z)]!=(z>0):return False
        signs[abs(z)]=z>0
    return (sum(z<0 for z in L)%2==0 and (W-p)%2==0
            and p_signed==(-1 if sum(z<0 for z in B)%2 else 1)*p
            and p>=max(6,max(map(abs,B)))
            and (W-p)//2>=max(8,max(map(abs,B)))
            and sum(abs(z)>=3 for z in B)>=2)
assert residual(L,12)

def weighted_star(L,weight):
    R=list(L);R.remove(1);R.sort(key=abs,reverse=True)
    rem=sum(map(abs,R));F={(0,0):(1,0)}
    for z in R:
        n=abs(z);s=1 if z>0 else -1;rem-=n;G=defaultdict(lambda:[0,0])
        for (i,j),(f,h) in F.items():
            for k in cg(i,n):
                if abs(k-1)+j<=rem:
                    G[k,j][0]+=f;G[k,j][1]+=h
            for k in cg(j,n):
                if abs(i-1)+k<=rem:
                    G[i,k][0]+=s*f;G[i,k][1]+=s*(h+weight(n)*f)
        F={q:tuple(v) for q,v in G.items() if any(v)}
    return F.get((1,0),(0,0))[1]

L2=[1]*2+[2]*4+[4]*10+[6]*12+[-7]*12+[-8]*6+[-9]*2+[-11]*12
assert residual(L2,-11)
for w,expected in ((lambda n:n-1,-529311115421531794338728680568405269354888985571690328),
                   (lambda n:n*n-1,-4964123857541026132159734948670969874427790952628165672)):
    value=weighted_star(L2,w);assert value==expected
    print("shifted-star obstruction:",value)
supplier=gain(L2,1,2)
assert supplier==46638554327418058868978815596190232129633997622673
print("successful (1,2) quarter-difference:",supplier)

def covered(L):
    for T in (tuple(L),tuple(-z if abs(z)%2 else z for z in L)):
        if all(z>0 for z in T):return True
        if -1 not in T:
            C=tuple(z for z in T if z!=1);a=T.count(1)
            if a>=2 and criterion(C,a):return True
        neg=[-z for z in T if z<0]
        pos=[z for z in T if z>0]
        if len(neg)!=2:continue
        for P,Qn in (neg,neg[::-1]):
            for q in set(pos):
                if q>=Qn or (Qn-q)%2:continue
                H=pos.copy();H.remove(q);H=[h for h in H if h!=2]
                if P>=sum(H) and (P-sum(H))%2==0:return True
    return False
for path in args.census:
    if not path.exists():
        print("optional census absent:",path);continue
    rows=[]
    for line in path.read_text().splitlines():
        m=re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?) phi=",line)
        if m:
            L=tuple(map(int,m[3].split()))+(int(m[2]),)
            rows.append(L)
    assert rows
    assert not any(covered(L) for L in rows)
    print("census",path.name,"rows",len(rows),"covered no-flip rows",0,
          "max factors",max(map(len,rows)),
          "max fundamentals",max(sum(abs(z)==1 for z in L) for L in rows))
print("PASS")
