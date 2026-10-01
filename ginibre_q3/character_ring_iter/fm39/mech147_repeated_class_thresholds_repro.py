import argparse, time
argparse.ArgumentParser(description="FM-MECH147 exact moment-tail verifier; no files written.").parse_args()
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb, prod
from random import Random
start=time.monotonic()

def add(A,B):
 C=dict(A)
 for k,v in B.items(): C[k]=C.get(k,0)+v
 return {k:v for k,v in C.items() if v}
def mul(A,B):
 C=defaultdict(int)
 for (i,j),a in A.items():
  for (k,l),b in B.items(): C[i+k,j+l]+=a*b
 return {k:v for k,v in C.items() if v}
@lru_cache(None)
def V(n):
 if n==0: return {(0,0):1}
 if n==1: return {(1,0):1,(0,1):1}
 return add(mul(V(1),V(n-1)),{k:-4*v for k,v in V(n-2).items()})
def feat(n,e):
 return {k:v*(1+e*(-1)**k[1]) for k,v in V(n).items() if 1+e*(-1)**k[1]}
@lru_cache(None)
def kernel(C,p):
 P={(0,0):1}
 for z in C: P=mul(P,feat(abs(z),1 if z>0 else -1))
 P=mul(P,feat(p,(-1)**sum(z<0 for z in C)))
 assert all(i%2==j%2==0 for i,j in P)
 return {(i//2,j//2):v for (i,j),v in P.items()}

def pa(A,B):
 C=A[:]+[0]*max(0,len(B)-len(A))
 for i,v in enumerate(B): C[i]+=v
 while len(C)>1 and C[-1]==0: C.pop()
 return C
def pm(A,B):
 C=[0]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  if a:
   for j,b in enumerate(B):
    if b: C[i+j]+=a*b
 while len(C)>1 and C[-1]==0: C.pop()
 return C
def ps(A,v): return [v*x for x in A]
def shift(A,v):
 C=[0]
 for x in reversed(A): C=pa(pm(C,[v,1]),[x])
 return C
def pe(A,v):
 z=0
 for x in reversed(A): z=z*v+x
 return z
@lru_cache(None)
def numerator(C,p,r,delta):
 K=dict(kernel(C,p))
 if delta:
  L={ij:-v for ij,v in K.items()}
  di,dj=(1,0) if r==1 else (1,1)
  for (i,j),v in K.items(): L[i+di,j+dj]=L.get((i+di,j+dj),0)+v
  K={ij:v for ij,v in L.items() if v}
 D=max(i+j for i,j in K)
 A=[[1]]
 for i in range(D):
  A.append(ps(pm(pm(A[-1],[2*i+1,2]),[2*i+3,2]),4))
 fill=[None]*(D+1); fill[D]=[1]
 for d in range(D-1,-1,-1):
  fill[d]=pm(pm(fill[d+1],[d+2,r]),[d+3,r])
 T=[[0] for _ in range(D+1)]
 for (i,j),v in K.items():
  R=ps(A[i],A[j][0]) if r==1 else pm(A[i],A[j])
  T[i+j]=pa(T[i+j],ps(R,v))
 out=[0]
 for d in range(D+1): out=pa(out,pm(T[d],fill[d]))
 assert out[-1]>0
 return out,D
def cert(C,p,r):
 G,_=numerator(C,p,r,False); P,_=numerator(C,p,r,True)
 def good(q): return min(shift(G,q))>=0 and min(shift(P,q-1))>=0
 lo,hi=0,1
 while not good(hi): lo,hi=hi,2*hi
 while hi-lo>1:
  mid=(lo+hi)//2
  if good(mid): hi=mid
  else: lo=mid
 return 2*hi
def cat(n): return comb(2*n,n)//(n+1)
@lru_cache(None)
def J(m,h):
 u=F(comb(2*m,m)*comb(2*h,h),comb(m+h,m))
 z=F(2*(2*m+1)*(2*h+1),(m+h+1)**2*(m+h+2))*u*u
 assert z.denominator==1
 return z.numerator
def value(C,p,r,q,delta=False):
 P,D=numerator(C,p,r,delta)
 den=prod((r*q+i+2)*(r*q+i+3) for i in range(D))
 z=F(pe(P,q)*J(q,0 if r==1 else q),den*2**(sum(map(abs,C))+p+1))
 assert z.denominator==1
 return z.numerator
def fusion(B,p):
 st={(0,0):1}; rem=sum(map(abs,B))
 for z in B:
  n=abs(z);e=1 if z>0 else -1;rem-=n;out=defaultdict(int)
  for (i,j),v in st.items():
   if j<=rem:
    lo=max(abs(i-n),p-rem);lo+=(i+n-lo)%2
    for h in range(lo,min(i+n,p+rem)+1,2): out[h,j]+=v
   if abs(i-p)<=rem:
    for h in range(abs(j-n),min(j+n,rem)+1,2): out[i,h]+=e*v
  st={ij:v for ij,v in out.items() if v}
 return st.get((p,0),0)
def par(C,p):
 t=sum(z<0 for z in C); h=(t+1)//2
 L=sum(F(abs(z)*(abs(z)+2),12) if z>0 else F((abs(z)-1)*(abs(z)+3),20) for z in C)
 L+=F(p*(p+2),12) if t%2==0 else F((p-1)*(p+3),20)
 return h,L
def old_threshold(C,p):
 h,L=par(C,p); K=60*L+1; nu=2*h+3
 a=2+(sum(map(abs,C))-p)%2
 while 15*a*a-K*nu*(2*a+nu)-15<0: a+=2
 return a

rng=Random(147); tests=0
for _ in range(12):
 C=tuple(n*rng.choice((-1,1)) for n in rng.sample(range(3,7),2))
 p=max(map(abs,C)); p+=(sum(map(abs,C))-p)%2
 for r in (1,2):
  for q in range(3):
   label=1 if r==1 else -2
   parent=fusion(C+(label,)*(2*q),p)
   childstep=fusion(C+(label,)*(2*q+2),p)-parent
   assert value(C,p,r,q)==parent
   assert value(C,p,r,q,True)==childstep
   tests+=2
print(tests,"independent fusion/numerator identities PASS",flush=True)

rows=[((-3,4,7,-10),12,1112,2),
      ((-3,5,-7),9,594,2),
      ((3,-4,7),8,460,2),
      ((-3,-6,7,8),10,994,2),
      ((-3,-4,-7,6),9,759,3),
      ((-3,)*64,6,22770,8)]
for C,p,old,new in rows:
 eta=(sum(map(abs,C))-p)%2
 assert old_threshold(C,p)==old
 assert cert(C+(1,)*eta,p,1)+eta==new
 print("fundamentals:",len(C),"core factors; p",p,"old/new",old,new,flush=True)

for C,p in [(rows[i][0],rows[i][1]) for i in range(3)]:
 assert cert(C,p,2)==2
 assert cert(C+(-2,),p,2)+1==3
print("three -2 families: every multiplicity >= 2 PASS",flush=True)
run=tuple(-n for n in range(3,16,2) for _ in range(2))
assert value(run,16,1,1)==214784436431
assert value(run,16,1,0)==296201364755
assert cert(run,16,1)==4
C=(-1,)+run
assert value(C,15,2,1)==703572078442
assert value(C,15,2,0)==1087917863463
assert cert(C,15,2)==4
C=(-3,)*64
assert 0<400*value(C,6,1,2)<value(C,6,1,1)
assert 0<value((-3,)*62,6,1,2)<value(C,6,1,2)
print("old counterexamples and the constant-4 obstruction PASS",flush=True)

def moment3(q,h):
 n=2*q; cur=J(q,h); bc=(-8)**n; total=0
 for j in range(n+1):
  z=cur;v=bc
  for i in range(n-j+1):
   total+=z*v
   if i<n-j:
    z=z*4*(2*(q+i)+1)*(2*(q+i)+3)//((q+i+h+j+2)*(q+i+h+j+3))
    v=v*(n-j-i)//(-8*(i+1))
  if j<n:
   cur=cur*4*(2*(h+j)+1)*(2*(h+j)+3)//((q+h+j+2)*(q+h+j+3))
   bc=bc*(n-j)*3//(-8*(j+1))
 d=16**q
 assert total%d==0
 return total//d
for q in range(4):
 assert moment3(q,1)==fusion((3,)*(2*q)+(-1,-1),0)
Q=[8,2,-24,24,-8,1]
R=[v*2**i for i,v in enumerate(Q)]
bern=[sum(F(R[j]*comb(i,j),comb(5,j)) for j in range(i+1)) for i in range(6)]
assert bern==list(map(F,[8,F(44,5),0,F(4,5),F(24,5),12]))
assert shift(Q,2)==[12,18,8,0,2,1]
assert par((-4,-4),6)==(1,F(61,10))
assert 3923*moment3(247,1)>246032*moment3(246,1)
print("+/-3 multiplicity >= 494, even, C=(-4,-4), p=6 PASS",flush=True)

def z_moments(N):
 C=[1];B=[F(1)]
 for n in range(N):
  c=lambda j:C[j] if j>=0 else 0
  b=lambda j:B[j] if j>=0 else F(0)
  num=2*n*(3*n+5)*c(n)+4*n*(n+5)*c(n-1)-24*n*(n-1)*c(n-2)
  den=(n+2)*(n+3)
  assert num%den==0
  C.append(num//den)
  B.append(((12*n*n+40*n+24)*b(n)+(8*n*n+28*n)*b(n-1)-48*n*(n-1)*b(n-2))/((2*n+5)*(n+4)))
 return C,B
Z,Beta=z_moments(267)
u=[F(1)]
for i in range(1,10): u.append(u[-1]*F(8*i,2*i+3))
for n in range(9):
 direct=lambda v:sum(comb(n,k)*(-2)**(n-k)*sum(comb(k,i)*v[i]*v[k-i] for i in range(k+1)) for k in range(n+1))
 assert Z[n]==direct([cat(i) for i in range(10)])
 assert Beta[n]==direct(u)
def atan_bounds(k,N=30):
 v=sum(F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(N))
 w=v+F((-1)**N,(2*N+1)*k**(2*N+1))
 return min(v,w),max(v,w)
l5,h5=atan_bounds(5);l239,h239=atan_bounds(239)
pi_lo,pi_hi=16*l5-4*h239,16*h5-4*l239
assert 3<pi_lo<pi_hi<4
def plus2_moment(q):
 n=2*q;A=F(Z[n+1]+2*Z[n]);D=F(128,9)*Beta[n]
 assert D>0
 return A-D/pi_lo**2,A-D/pi_hi**2
assert par((-3,5),6)==(1,F(173,30))
assert par((3,-5),6)==(1,F(51,10))
lo,hi=plus2_moment(133);pl,ph=plus2_moment(132)
assert 1259*lo>43644*ph
print("+2 multiplicity >= 266, even, C=(-3,+5), p=6 PASS",flush=True)
print("ALL CHECKS PASS; seconds =",round(time.monotonic()-start,3))
