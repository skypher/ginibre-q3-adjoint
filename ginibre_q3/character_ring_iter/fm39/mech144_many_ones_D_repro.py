import argparse, time
argparse.ArgumentParser(description='Exact verifier for the FM-MECH144 Jacobi descent criterion.').parse_args()
started=time.monotonic()
from fractions import Fraction as F
from math import comb, factorial
from functools import lru_cache
from collections import defaultdict
def add(A,B):
 C=dict(A)
 for k,v in B.items():C[k]=C.get(k,0)+v
 return {k:v for k,v in C.items() if v}
def mul(A,B):
 C=defaultdict(int)
 for (i,j),a in A.items():
  for (k,l),b in B.items():C[i+k,j+l]+=a*b
 return {k:v for k,v in C.items() if v}
@lru_cache(None)
def V(n):
 if n==0:return {(0,0):1}
 if n==1:return {(1,0):1,(0,1):1}
 return add(mul({(1,0):1,(0,1):1},V(n-1)),{k:-4*v for k,v in V(n-2).items()})
def feature(n,e):
 return {k:v*(1+e*(-1)**k[1]) for k,v in V(n).items() if 1+e*(-1)**k[1]}
@lru_cache(None)
def kernel(C,p):
 P={(0,0):1};t=sum(z<0 for z in C)
 for z in C:P=mul(P,feature(abs(z),1 if z>0 else -1))
 P=mul(P,feature(p,(-1)**t))
 assert all(j%2==0 for i,j in P)
 return P
def norm(C,p,a):
 S=sum(map(abs,C));eta=(S+p)%2
 assert (a+S-p)%2==0
 m=(a+eta)//2
 @lru_cache(None)
 def ratio(i,h):
  if i:
   q=m+i-1
   return ratio(i-1,h)*F(4*(2*q+1)*(2*q+3),(q+h+2)*(q+h+3))
  if h:
   q=h-1
   return ratio(0,q)*F(4*(2*q+1)*(2*q+3),(m+q+2)*(m+q+3))
  return F(1)
 value=sum(F(v)*ratio((i-eta)//2,j//2) for (i,j),v in kernel(tuple(C),p).items())/2**(S+p+1)
 return value,m
def par(C,p):
 t=sum(z<0 for z in C);h=(t+1)//2
 alp=lambda n:F(n*(n+1)*(n+2),6)
 lam=sum((F(abs(z)*(abs(z)+2),12) if z>0 else F((abs(z)-1)*(abs(z)+3),20)) for z in C)
 c=F(1)
 for z in C:c*=2*(z+1) if z>0 else alp(-z)
 if t%2:c*=alp(p)/2;lam+=F((p-1)*(p+3),20)
 else:c*=p+1;lam+=F(p*(p+2),12)
 return h,lam,c
def eps(C,p,a):
 h,L,c=par(C,p)
 return L*F(4)*(1-F((a+1)*(a+3),(a+2*h+4)*(a+2*h+6)))
def score_minus(C,p,a):
 h,L,c=par(C,p);Cr=tuple(C[2:]);hr,Lr,cr=par(Cr,p)
 assert h==hr+1
 rJ=F(16*(2*h-1)*(2*h+1),(a+2*h+2)*(a+2*h+4))
 return c/cr*rJ*(1-eps(C,p,a)),c/cr*rJ,1-eps(Cr,p,a)

from random import Random
def catalan(m):return comb(2*m,m)//(m+1)
def J(m,h):
 u=F(comb(2*m,m)*comb(2*h,h),comb(m+h,m))
 q=F(2*(2*m+1)*(2*h+1),(m+h+1)**2*(m+h+2))*u*u
 assert q.denominator==1
 return q.numerator
def Jdirect(m,h):
 D=2*(m+h);ans=0
 for j in range(0,D+1,2):
  coeff=sum((-1)**k*comb(2*h,k)*comb(2*m,j-k)
            for k in range(max(0,j-2*m),min(2*h,j)+1))
  ans+=coeff*catalan(j//2)*catalan((D-j)//2)
 return ans
def fusion(B,p):
 st={(0,0):1};rem=sum(map(abs,B))
 for z in B:
  n=abs(z);e=1 if z>0 else -1;rem-=n;out=defaultdict(int)
  for (i,j),v in st.items():
   if j<=rem:
    lo=max(abs(i-n),p-rem);lo+=(i+n-lo)%2
    for h in range(lo,min(i+n,p+rem)+1,2):out[h,j]+=v
   if abs(i-p)<=rem:
    for h in range(abs(j-n),min(j+n,rem)+1,2):out[i,h]+=e*v
  st={ij:v for ij,v in out.items() if v}
 return st.get((p,0),0)
def Q(C,p,a):
 h,L,c=par(C,p);K=60*L+1;assert K.denominator==1;b=2*h+3
 return 15*a*a-K*(2*b*a+b*b)-15
for m in range(9):
 for h in range(9):
  assert J(m,h)==Jdirect(m,h)
  assert F(J(m+1,h),J(m,h))==F(4*(2*m+1)*(2*m+3),(m+h+2)*(m+h+3))
print("81 moment identities pass",flush=True)
rng=Random(144);total=0
for _ in range(20):
 C=tuple(n*rng.choice((-1,1)) for n in rng.sample(range(2,8),rng.randint(2,4)))
 p=max(map(abs,C))+rng.randint(0,3)
 h,L,c=par(C,p);P=kernel(C,p)
 den=2**(sum(map(abs,C))+p+1)*c
 assert sum(v*4**i for (i,j),v in P.items() if j==2*h)==den
 assert sum(i*v*4**(i-1) for (i,j),v in P.items() if j==2*h and i)/den==L
 for a in range(2,10):
  if (a+sum(map(abs,C))-p)%2:continue
  z,m=norm(C,p,a)
  assert z*J(m,0)==fusion([1]*a+list(C),p)
  total+=1
print(total,"independent fusion/Jacobi comparisons pass",flush=True)
for C,p in [((-3,4,7,-10),12),((-3,5,-7),9),((3,-4,7),8),((-3,-6,7,8),10),((-3,-4,-7,6),9)]:
 a=2+(2+sum(map(abs,C))-p)%2
 while Q(C,p,a)<0:a+=2
 for aa in (a,a+2,a+20):
  x,m=norm(C,p,aa);y,_=norm(C,p,aa-2)
  base=F(4*(2*m-1)*(2*m+1),(m+1)*(m+2))
  assert x>0 and x>=y/base
  h,L,c=par(C,p);b=2*h+3
  assert eps(C,p,aa)<1
  if aa%2==0:
   actual=x*J(m,0);child=y*J(m-1,0)
   bound=c*J(aa//2-1,h)*Q(C,p,aa)/((aa+b)**2-1)
   assert actual-child>=bound>=0
   assert c*J(aa//2,h)*(1-eps(C,p,aa))<=actual<=c*J(aa//2,h)
 print("certified C,p,a0:",C,p,a,"Q",Q(C,p,a),flush=True)
for a in (5536,20000):
 C=(-21,-21);p=22
 score,upper,childlow=score_minus(C,p,a)
 x,m=norm(C,p,a);y,_=norm((),p,a);z,_=norm(C,p,a-2)
 base=F(4*(2*m-1)*(2*m+1),(m+1)*(m+2))
 assert x>z/base
 assert (x>y)==(a==5536)
 assert score>=1 if a==5536 else upper<childlow
 print("switch a",a,"safe-score",score,"upper",upper,"child-lower",childlow,flush=True)
C=tuple(-n for n in range(3,16,2) for _ in range(2))
x,m=norm(C,16,2);y,_=norm(C,16,0)
assert x*J(m,0)==214784436431 and y==296201364755
assert Q(C,16,2)<0
print("old small-removal obstruction recovered",flush=True)
C=tuple([-11]*5+[-54,-2]+[n for n in range(3,54,2) if n!=11 for _ in range(54//n)])
h,L,c=par(C,55)
assert (h,L)==(4,F(47201,15))
assert eps(C,55,54)==F(1935241,510)
assert Q(C,55,54)==-247102020
assert Q(C,55,276918)<0<Q(C,55,276920)
print("WM core threshold: a >= 276920, a even; epsilon at 54 =",eps(C,55,54),flush=True)
print("ALL CHECKS PASS; seconds =",round(time.monotonic()-started,3))