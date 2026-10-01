"""FM-MECH97: charging, folding, a uniform obstruction, and one paid chord."""
import argparse
from collections import defaultdict
from functools import lru_cache
from math import comb

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--box",type=int,default=16)
ap.add_argument("--max-b",type=int,default=12)
ap.add_argument("--fold-n",type=int,default=24)
o=ap.parse_args()

def C(n,j):
    return comb(n,j) if 0<=j<=n else 0

@lru_cache(None)
def walks(B):
    z={(0,1):1}
    for _ in range(B):
        w=defaultdict(int)
        for (i,j),v in z.items():
            w[i,j]+=(1 if j==i+1 else 2)*v
            for ij in ((i-1,j-1),(i-1,j+1),(i+1,j+1)):
                w[ij]+=v
            if j>i+1:
                w[i+1,j-1]+=v
        z=dict(w)
    assert all(w>0 and i<j and (j-i)%2 for (i,j),w in z.items())
    return z

def folded(N,eps,k,B):
    z=defaultdict(int)
    for (i,j),w in walks(B).items():
        i+=k
        j+=k
        if i < -1 or j > N+1:
            continue
        u=max(i,N-i)
        v=max(j,N-j)
        if u==v:
            continue
        s=eps**((2*i<N)+(2*j<N))
        if u>v:
            u,v,s=v,u,-s
        z[u,v]+=s*w
    return {ij:w for ij,w in z.items() if w}

def row(a,e):
    N=a+e
    c=[1]
    for j in range(N):
        z,rem=divmod((a-e)*c[j]-(N-j+1)*(c[j-1] if j else 0),j+1)
        assert rem==0
        c.append(z)
    return c

def tools_for(cc):
    N=len(cc)-1
    def c(j):
        return cc[j] if 0<=j<=N else 0
    def D(j):
        return c(j)**2-c(j-1)*c(j+1)
    def W(i,j):
        return c(i)*(c(j-1)+c(j+1))-(c(i-1)+c(i+1))*c(j)
    @lru_cache(None)
    def v(h,j):
        if h<0:
            return 0
        if h==0:
            return c(j)
        return v(h-1,j-1)+v(h-1,j+1)
    @lru_cache(None)
    def T(B,j):
        if B<0:
            return 0
        return sum(C(B,h)*2**(B-h)*(
            v(h,j)**2-v(h,j-1)*v(h,j+1)) for h in range(B+1))
    def F(B,j):
        return T(B,j)-T(B,j+1)
    return c,D,W,v,T,F

expected=[1,0,-4,-34,-216,-1288,-7520,-43538,
          -251608,-1454104,-8413536,-48754088,-282978080]
mins=[]
for B in range(o.max_b+1):
    z=walks(B)
    values=[]
    for l in range(-B,B+1):
        cost=sum(w for (i,j),w in z.items()
                 if j>i+1 and i<=l<j)
        values.append(z.get((l,l+1),0)-cost)
    mins.append(min(values))
    edge_cost=sum(w for (i,j),w in z.items() if i==-B and j>i+1)
    assert z[-B,-B+1]-edge_cost == 2-2**B
    assert z[-B,B+1]==1
    assert all(z[l,l+1]>=1 for l in range(-B,B+1))
    if B<len(expected):
        assert mins[-1]==expected[B]
print("MINIMUM WALK BUDGETS",mins,flush=True)

fold_count=fold_negative=0
for N in range(1,o.fold_n+1):
    for eps in (1,-1):
        for B in range(o.max_b+1):
            for k in range((N+1)//2,N+B+2):
                z=folded(N,eps,k,B)
                fold_count+=1
                fold_negative+=any(w<0 for w in z.values())
assert folded(3,-1,2,1)=={(2,3):1,(2,4):-1,(3,4):1}
print("FOLDED TABLES",fold_count,
      "WITH NEGATIVE WEIGHT",fold_negative,flush=True)

def charging(N,eps,k,B,D,W):
    z=folded(N,eps,k,B)
    value=credit=cost=0
    for (i,j),w in z.items():
        assert 2*i>=N and j<=N+1
        assert D(i)-D(j)>=abs(W(i,j))
        q=w*W(i,j)
        value+=q
        if q>=0:
            credit+=q
        else:
            cost+=abs(w)*(D(i)-D(j))
    return value,credit,cost

profiles=fail=energy_checks=0
first=large=None
for a in range(o.box+1):
    for e in range(a+1):
        N=a+e
        c,D,W,v,T,F=tools_for(row(a,e))
        for B in range(o.max_b+1):
            for k in range(N//2+1,N+B+1):
                p=2*k-N
                d=N+B-k
                value,credit,cost=charging(N,(-1)**e,k,B,D,W)
                assert value==F(B,k)>=0
                assert value==sum(w*W(k+i,k+j)
                                  for (i,j),w in walks(B).items())
                norms=lambda j:sum(C(B,h)*2**(B-h)*v(h,j)**2
                                   for h in range(B+1))
                assert (k+B+2)*value==(
                    (p+1)*T(B,k)+norms(k)-norms(k+1)-2*B*F(B-1,k))
                energy_checks+=1
                profiles+=1
                if credit<cost:
                    fail+=1
                    case=(a,e,B,p,k,d,value,credit,cost)
                    if first is None:
                        first=case
                    if (large is None and min(a,e)>=2 and B>=3
                            and p>=11 and d>=8):
                        large=case
print("ACTUAL PROFILES",profiles,"CHARGING FAILURES",fail,flush=True)
print("FIRST FAILURE",first,flush=True)
print("LARGE-LABEL FAILURE",large,flush=True)
print("ACTUAL ENERGY COMPARISONS",energy_checks,flush=True)

for a,e,B,p in [(8,8,2,4),(40,40,40,82)]:
    N=a+e
    c,D,W,v,T,F=tools_for(row(a,e))
    k=(N+p)//2
    value,credit,cost=charging(N,(-1)**e,k,B,D,W)
    assert value==F(B,k)>0 and credit-cost<0
    if B==2:
        assert (value,credit,cost)==(5012,12348,21896)
    print("ACTUAL CONTROL",(a,e,B,p),"F",value,
          "CHARGING LOWER BOUND",credit-cost,flush=True)

def cat(j):
    return C(2*j,j)//(j+1)

def invariant(m):
    return (
        sum(C(m,2*j)*2**(m-2*j)*cat(j+1)*cat(j)
            for j in range(m//2+1))
        -sum(C(m,2*j+1)*2**(m-2*j-1)*cat(j+1)**2
             for j in range((m-1)//2+1)))

obstructions=0
for N in range(4,34,2):
    cc=[0]*(N+1)
    for j in range(N//2+1):
        cc[2*j]=(-1)**j
    c,D,W,v,T,F=tools_for(cc)
    assert all(D(j)==1 for j in range(N+1))
    assert D(N+1)==0
    for i in range(N//2,N+1):
        for j in range(i+1,N+2):
            assert D(i)-D(j)>=abs(W(i,j))
    for B in range(2,min(o.max_b,N-2)+1):
        I=invariant(B-2)
        assert I>=1
        assert F(B,N-1)==-B*(B-1)*I//2<0
        assert F(B,N-1)==sum(w*W(N-1+i,N-1+j)
                            for (i,j),w in walks(B).items())
        assert 2*c(2)+N*c(0)==N-2>0
        obstructions+=1
print("UNIFORM ABSTRACT OBSTRUCTION CHECKS",obstructions,flush=True)
print("SP4 INVARIANT COUNTS",[invariant(m) for m in range(11)],flush=True)

N=16
cc=[(-1)**(j//2) if j%2==0 else 0 for j in range(N+1)]
c,D,W,v,T,F=tools_for(cc)
assert F(7,15)==-567
assert F(2,15)==-1
print("ABSTRACT RESIDUAL-SIZE CONTROL",
      "N=16, B=7, p=14, d=8, F=-567",flush=True)

def old_band(a,e,B,p):
    S=a+e+2*B+2
    R=p*(p+2)*(S-p)*(S+p+2)
    A=abs(a-e)*(p+1)
    H=B*(S+3*p+4)
    Z=2*R-2*A*A-H*H
    return 0<p<S and Z>=0 and Z*Z>=8*A*A*H*H

def criterion89(a,e,B,p):
    n=(a+e-p)//2
    if n<2*B:
        return False
    q=p*(p+2)
    A=abs(a-e)*(p+1)
    C2=32*B*B*(q+(a-e)**2)
    Z=4*q*(n+1-2*B)**2-A*A-C2
    return Z>=0 and Z*Z>=4*A*A*C2

a,e,B,p=180,20,8,118
N=a+e
k=(N+p)//2
L,U=k-B,k+B+1
c,D,W,v,T,F=tools_for(row(a,e))
z=walks(B)
assert 2*L>=N and U<=N+1
assert W(L,U)<0
assert all(W(k+i,k+j)>=0 for (i,j) in z
           if (i,j)!=(-B,B+1))
paid=D(L)-D(U)+W(L,U)
rest=sum((z[l-k,l+1-k]-1)*W(l,l+1) for l in range(L,U))
rest+=sum(w*W(k+i,k+j) for (i,j),w in z.items()
          if j>i+1 and (i,j)!=(-B,B+1))
assert paid>=0 and rest>=0 and paid+rest==F(B,k)>0
assert not old_band(a,e,B,p)
assert not criterion89(a,e,B,p)
assert (p-2*B)**2<4*(e-1)*(N-e+2)
print("PROVED CHARGING REGION CONTROL",(a,e,B,p),
      "WINDOW",(L,U),"OUTSIDE 86/89/92",flush=True)

def direct(a,e,B):
    N=a+e
    degree=N+2*B
    cc=row(a,e)
    mon=[0]*(degree+1)
    for j in range(0,N+1,2):
        for l in range(B+1):
            for h in range(l+1):
                mon[N-j+2*(l-h)]+=(
                    cc[j]*C(B,l)*(-2)**(B-l)*C(l,h)*cat(j//2+h))
    return [sum(mon[m]*(C(m,(m-p)//2)-C(m,(m-p)//2-1))
                for m in range(p,degree+1,2))
            for p in range(degree+1)]

c,D,W,v,T,F=tools_for(row(180,20))
assert direct(180,20,8)[118]==F(8,159)>0
print("NEW REGION DIRECT CATALAN CHECK PASS",flush=True)

bridges=0
for a in range(7):
    for e in range(7):
        N=a+e
        c,D,W,v,T,F=tools_for(row(a,e))
        for B in range(5):
            for p,q in enumerate(direct(a,e,B)):
                if (N+p)%2:
                    assert q==0
                else:
                    assert q==F(B,(N+p)//2)
                bridges+=1
print("CATALAN BRIDGES",bridges,flush=True)
print("PASS",flush=True)