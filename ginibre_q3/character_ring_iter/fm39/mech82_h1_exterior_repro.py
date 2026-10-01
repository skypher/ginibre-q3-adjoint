"""FM-MECH82: one-label combined rows and a uniform exterior theorem.
Run from /home/yang/q3adjoint. No files are written.
"""
import argparse
from math import comb, prod
from functools import lru_cache
from fractions import Fraction as Qf
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-r",type=int,default=44)
ap.add_argument("--max-m",type=int,default=128)
args=ap.parse_args()

def choose(n,k):
 return comb(n,k) if 0<=k<=n else 0

@lru_cache(None)
def row(a,e):
 N=a+e
 c=[1]
 for k in range(N):
  z,rem=divmod(
   (a-e)*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
  assert rem==0
  c.append(z)
 return tuple(c)

def data(a,e):
 N=a+e
 rr=row(a,e)
 def c(k):
  return rr[k] if 0<=k<=N else 0
 @lru_cache(None)
 def v(h,k):
  if h<0:return 0
  return sum(comb(h,j)*c(k+h-2*j) for j in range(h+1))
 @lru_cache(None)
 def dh(h,k):
  return v(h,k)**2-v(h,k-1)*v(h,k+1)
 @lru_cache(None)
 def T(b,k):
  return sum(comb(b,h)*2**(b-h)*dh(h,k)
             for h in range(b+1))
 def W(b,j,i):
  return sum(comb(b,h)*2**(b-h)*
   (v(h,j)*v(h+1,i)-v(h+1,j)*v(h,i))
   for h in range(b+1))
 def F(b,k):
  return T(b,k)-T(b,k+1)
 return N,c,v,dh,T,W,F

# Universal adjacent insertion identity.
zm2,zm1,z0,z1,z2,z3=S.symbols("zm2 zm1 z0 z1 z2 z3")
d0=z0*z0-zm1*z1
d1=z1*z1-z0*z2
db0=(zm1+z1)**2-(zm2+z0)*(z0+z2)
db1=(z0+z2)**2-(zm1+z1)*(z1+z3)
rhs=zm1**2-zm2*z0-z2**2+z1*z3
rhs+=zm1*(z1+z3)-(zm2+z0)*z2
assert S.expand(2*(d0-d1)+db0-db1-rhs)==0

# Uniform root bound and escaping-family polynomial.
R,h,l,t=S.symbols("R h l t")
assert S.expand(
 h*(R+2)-l*(R+h+2-l)-(h-l)*(R+2-l))==0
assert S.expand((t+4)**4-12*((t+4)**2+(t+4)-1))==(
 t**4+16*t**3+84*t**2+148*t+28)
print("SYMBOLIC IDENTITIES: PASS",flush=True)

def Kraw(n,d,X):
 if d==0:return 1
 prev,cur=1,X
 for j in range(1,d):
  prev,cur=cur,X*cur-j*(n-j+1)*prev
 return cur

@lru_cache(None)
def Q(R,h,l):
 return sum((-1)**j*comb(h,j)*choose(R,l-j)
            for j in range(h+1))

duals=0
for R0 in range(19):
 for h0 in range(R0+1):
  n=R0+h0
  den=prod(range(n-h0+1,n+1))
  for l0 in range(n+1):
   assert den*Q(R0,h0,l0)==(
    (-1)**h0*comb(n,l0)*Kraw(n,h0,2*l0-n))
   duals+=1
print("KRAWTCHOUK DUALITY",duals,flush=True)

recurrences=insertions=0
for a in range(9):
 for e in range(9):
  N,c,v,dh,T,W,F=data(a,e)
  for h0 in range(5):
   for k in range(-h0,N+h0+2):
    assert (k+h0+1)*v(h0,k+1)==(
     (a-e)*v(h0,k)-(N+h0-k+1)*v(h0,k-1)
     +4*h0*v(h0-1,k))
    recurrences+=1
  for b in range(1,5):
   for k in range((N+1)//2,N+b+1):
    assert F(b,k)==(
     T(b-1,k-1)-T(b-1,k+2)+W(b-1,k-1,k+2))
    insertions+=1
print("COUPLED RECURRENCES",recurrences,
      "INSERTIONS",insertions,flush=True)

# Parity identities, including points outside the proved region.
parity_checks=0
for R0 in range(13):
 N,c,v,dh,T,W,F=data(R0,R0)
 N1,c1,v1,dh1,T1,W1,F1=data(R0+1,R0)
 for h0 in range(7):
  for k in range(2*R0+h0+4):
   if (k+h0)%2==0:
    l0=(k+h0)//2
    bal=Q(R0,h0,l0)*(Q(R0,h0,l0)-Q(R0,h0,l0+1))
    adj=Q(R0,h0,l0)*(Q(R0,h0,l0-1)-Q(R0,h0,l0+1))
   else:
    l0=(k+h0+1)//2
    bal=Q(R0,h0,l0)*(Q(R0,h0,l0-1)-Q(R0,h0,l0))
    l1=(k+h0-1)//2
    adj=Q(R0,h0,l1)**2-Q(R0,h0,l1+1)**2
   assert dh(h0,k)-dh(h0,k+1)==bal
   assert dh1(h0,k)-dh1(h0,k+1)==adj
   parity_checks+=2
print("PARITY FACTORIZATIONS",parity_checks,flush=True)

# Bounded vetting of the uniform exterior theorem.
balanced=adjacent=0
for R0 in range(1,args.max_r+1):
 N,c,v,dh,T,W,F=data(R0,R0)
 for b in range(R0+1):
  for s in range(R0+b+1):
   if s*s<4*b*(R0+2):continue
   k=R0+s
   for h0 in range(b+1):
    assert dh(h0,k)-dh(h0,k+1)>=0
    balanced+=1
   assert F(b,k)>=0
for R0 in range(1,max(1,args.max_r-3)):
 N,c,v,dh,T,W,F=data(R0+1,R0)
 for b in range(R0+1):
  for s in range(R0+b+1):
   if s*s<4*b*(R0+2):continue
   k=R0+s+1
   for h0 in range(b+1):
    assert dh(h0,k)-dh(h0,k+1)>=0
    adjacent+=1
   assert F(b,k)>=0
print("EXTERIOR CHANNELS",balanced,adjacent,flush=True)

# Independent Catalan moments and ballot conversion.
@lru_cache(None)
def cat(j):
 return comb(2*j,j)//(j+1)

def direct(a,e,b):
 N=a+e
 D=N+2*b
 mon=[0]*(D+1)
 cc=row(a,e)
 for j in range(0,N+1,2):
  for l0 in range(b+1):
   for h0 in range(l0+1):
    mon[N-j+2*(l0-h0)]+=(
     cc[j]*comb(b,l0)*(-2)**(b-l0)*
     comb(l0,h0)*cat(j//2+h0))
 return [
  sum(mon[d]*(choose(d,(d-p)//2)-choose(d,(d-p)//2-1))
      for d in range(p,D+1,2))
  for p in range(D+1)]

bridges=0
for a in range(9):
 for e in range(9):
  for b in range(5):
   if a+e+2*b>24:continue
   N,c,v,dh,T,W,F=data(a,e)
   actual=direct(a,e,b)
   for p,val in enumerate(actual):
    if (p+N)%2:
     assert val==0
    else:
     assert val==F(b,(N+p)//2)
    bridges+=1
print("CATALAN BRIDGES",bridges,flush=True)

# The escaping family: additional checks using four channels.
def bv(R0,h0,k):
 if (k+h0)%2:return 0
 l0=(k+h0)//2
 return (-1)**l0*Q(R0,h0,l0)

def bd(R0,h0,k):
 return bv(R0,h0,k)**2-bv(R0,h0,k-1)*bv(R0,h0,k+1)

family=0
example=None
for m in range(4,args.max_m+1):
 R0=m*m+m-3
 s=m*m
 k=R0+s
 assert R0%2==1 and s*s>=12*(R0+2)
 parts=[bd(R0,h0,k)-bd(R0,h0,k+1) for h0 in range(4)]
 assert min(parts)>=0
 chosen=0 if m%2 else 1
 assert parts[chosen]>0
 value=sum(comb(3,h0)*2**(3-h0)*parts[h0]
           for h0 in range(4))
 assert value>0
 family+=1
 if m==23:example=2*value
print("ESCAPING FAMILY",family,"PASS",flush=True)
print("RAW EXPECTATION AT m=23",example,flush=True)

# Global channelwise positivity fails.
N,c,v,dh,T,W,F=data(10,10)
parts=[dh(h0,16)-dh(h0,17) for h0 in range(4)]
assert parts==[1575,1400,560,-546]
assert F(3,16)==32214
print("ROOT-CROSSING CHANNELS",parts,
      "TOTAL",F(3,16),flush=True)

# Fixed-output Newton screen: an unproved strengthening.
newton_checks=0
for a in range(13):
 for e in range(a+1):
  if a+e==0:continue
  N,c,v,dh,T,W,F=data(a,e)
  for k in range((N+1)//2,N+19):
   layer=[F(b,k) for b in range(19)]
   for order in range(1,19):
    layer=[y-x for x,y in zip(layer,layer[1:])]
    assert min(layer)>=0,(a,e,2*k-N,order)
    newton_checks+=len(layer)
print("FIXED-OUTPUT NEWTON ENTRIES",
      newton_checks,"NO FAILURE",flush=True)

x,theta=S.symbols("x theta")
assert S.expand(
 x**3-(1+theta)*x-(x**3-2*x)-(1-theta)*x)==0
print("SHIFT-1 NECESSARY BOUND: PASS",flush=True)

# Shifted transfer and the whole-U-positive-cone obstruction.
def clean(P):
 return {k:v for k,v in P.items() if v}

def add1(*terms):
 out={}
 for c,P in terms:
  for k,v in P.items():
   out[k]=out.get(k,0)+c*v
 return clean(out)

def xp1(P):
 return add1((1,P),(1,{k+2:v for k,v in P.items()}))

def resolvent(P,D):
 return clean({k:Qf(v,D+4-k) for k,v in P.items()})

def shifted_mon(a,e,b):
 N=a+e
 out={}
 for j in range(0,N+1,2):
  for l0 in range(b+1):
   for h0 in range(l0+1):
    k=N-j+2*(l0-h0)
    z=(row(a,e)[j]*comb(b,l0)*(-3)**(b-l0)*comb(l0,h0)
       *choose(2*(j//2+h0),j//2+h0)//(j//2+h0+1))
    out[k]=out.get(k,0)+z
 return clean(out)

shift_checks=0
for a in range(6):
 for e in range(6):
  if a+e==0:continue
  for b in range(5):
   D=a+e+2*b
   current=shifted_mon(a,e,b)
   previous=shifted_mon(a,e,b-1) if b else {}
   brace=add1((b,xp1(previous)),(-(b+2),current))
   image=add1((1,xp1(current)),(6,resolvent(brace,D)))
   assert image==shifted_mon(a,e,b+1)
   shift_checks+=1

seed={3:1,1:-2}
assert add1((1,xp1(seed)),(-12,resolvent(seed,3)))=={
 5:1,3:-4,1:2}
print("SHIFTED TRANSFER",shift_checks,"PASS",flush=True)
print("UNIVERSAL-SEED CONTROL: U_3 -> U_5 - U_1",flush=True)
print("PASS",flush=True)