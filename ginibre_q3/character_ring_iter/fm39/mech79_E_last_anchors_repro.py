"""FM-MECH79: uniform completion of the larger-anchor range of (E).
Run from the repository root with python3 -u; no files are written.
The finite certificates prove parameter bounds. The row census is diagnostic.
"""
import argparse
from fractions import Fraction as F
from math import comb, factorial, isqrt
from functools import lru_cache
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n",type=int,default=140,
                help="upper N in the bounded row diagnostic (default: 140)")
args=ap.parse_args()

# The coefficient of delta in the frozen-coordinate phase estimate.
rho,alpha,sigma,xx=S.symbols("rho alpha sigma xx",positive=True)
kap=rho*rho*sigma/(2*alpha)
VV=rho*rho*sigma*sigma*(1-xx*xx)+2*sigma*(1-rho*xx)
assert S.cancel(alpha*alpha*VV/(rho*rho*sigma*sigma*(1-rho*rho))
 -(alpha*alpha*(1-xx*xx)+alpha*(1-rho*xx)/kap)/(1-rho*rho))==0
print("PHASE NORMALIZATION: PASS",flush=True)

# If h<=1/50 and L=kappa*h<=11/(10*sqrt(h)), the phase budget fails.
A=F(65,64); b=F(25,49); c=F(1,64); q=F(1,7); lam=F(11,10)
upper=A*A/(1-b*b)*(
 (1+c)*lam*lam+(3*(1+c)*b+1)*lam*q
 +(3*(1+c)*b*b+2*b)*F(1,50)
 +((1+c)*b**3+b*b)/15)
assert upper<F(12,5)<F(157,100)**2
print("PHASE-BUDGET CONSTANT: PASS",upper)

# t=h^(1/4) <= 19/50; checkpoint delta0=9/(4t).
t=F(19,50)
assert t**4>F(1,50)
remaining=F(11,10)/t**2-F(9,4)/t-b
assert remaining>1
print("CHECKPOINT REMAINING LENGTH > 1:",remaining)

# Angular coefficient A_x <= 7t^2/4.
ang=(2*A*A+10*A*t*t/11)/(1-b*b)
assert ang<F(7,4)**2
turn=F(7,4)*t*t*(F(9,4)/t+b)+F(1,98)
assert turn<F(157,50)
assert F(1,98)-F(1,98)**3/6>F(1,99)
print("CHECKPOINT PRECEDES EVERY INNER LONG ENDPOINT: PASS")

# w(delta0) < 1/10.
def ex(q,n=36):
 return sum(q**k/F(factorial(k)) for k in range(n+1))
pref=F(100,99)/(t**8*(1-F(45,22)*t))
assert ex(F(441,100)/t)>10*pref
assert F(441,100)-8*t>0
assert F(1,15)/(F(99,50))-F(49,25)+1<0
margin=F(99,100)**2-4*F(33,32)**2/F(10)
assert margin>0
print("MOVING CHECKPOINT MARGIN",margin)

# A stronger frozen crossing bound valid for x>=7/10.
r=S.symbols("r")
pol=F(107,100)**2*(1-r*r)-A*A*F(51,100)-A*F(1,15)*(1-F(7,10)*r)
assert S.Poly(pol,r).coeff_monomial(r*r)<0
assert pol.subs(r,0)>0 and pol.subs(r,F(5,7))>0
assert F(119,1000)-F(119,1000)**3/6>F(2,17)
assert F(14,5)*F(107,100)+F(119,1000)<F(157,50)
print("LARGER-ANCHOR CROSSING delta > 14/5: PASS",pol.subs(r,F(5,7)))

# (lower endpoint numerator / 100, kappa cutoff, gamma numerator / 10000).
TABLE=[
(70,18,9525),(71,18,9525),(72,18,9525),(73,19,9539),
(74,19,9539),(75,20,9552),(76,20,9552),(77,22,9586),
(78,24,9624),(79,26,9655),(80,28,9682),(81,30,9705),
(82,32,9726),(83,35,9751),(84,39,9778),(85,43,9800),
(86,47,9818),(87,52,9837),(88,58,9854),(89,67,9875),
(90,76,9890),(91,89,9907),(92,106,9922),(93,128,9936),
(94,162,9950),(95,210,9961),(96,292,9972),(97,400,9980)]

def ceil(q):
 return -((-q.numerator)//q.denominator)

def amin(e,lo):
 b=e+1; p,q=lo.numerator,lo.denominator
 A=max(b+2,ceil(F(1024,b)),ceil(F(15*(b+1),2*b-15)))
 good=lambda a:(a+b)**2*q*q>16*p*p*a*b
 if good(A):
  return A
 left=A; right=2*A
 while not good(right):
  right*=2
 while right-left>1:
  mid=(left+right)//2
  if good(mid):
   right=mid
  else:
   left=mid
 return right

def minors(e,A,lo):
 N=A+e-1; C=4*A*(e+1); p,q=lo.numerator,lo.denominator
 old,new=1,p
 ans=[old,new]
 for l in range(1,e):
  nxt=p*(C if l%2 else 1)*new-l*(N-l+1)*q*q*old
  ans.append(nxt)
  old,new=new,nxt
 return ans

assert F(699,1000)**2<F(24,49)
assert F(1699,1000)*236*F(472,473)>400
roots=minor_count=0
for k,K,G in TABLE:
 lo=F(k,100)
 for e in range(7,235):
  A=amin(e,lo)
  if F(2*A*(e+1),A+e+2)>=K:
   continue
  mm=minors(e,A,lo)
  assert min(mm)>0,(k,K,e,A)
  roots+=1
  minor_count+=len(mm)
print("ROOT CERTIFICATES",roots,"POSITIVE MINORS",minor_count,flush=True)

def ex(q,n=40):
 return sum(q**h/F(factorial(h)) for h in range(n+1))

for k,K,G in TABLE:
 lo,hi=F(k,100),F(k+1,100)
 gamma=F(G,10000)
 step=min(F(1,32),1/(2*lo*K))
 assert hi+F(14,5*K)<gamma<gamma+step<1
 for z in (gamma,gamma+step):
  power=K*(z*z-hi*hi)
  pref=(1+z)/((1+hi)*(1-hi)*(1-z))
  assert ex(power)>5*pref,(k,K,G,z)
 assert 4*gamma/(1-gamma*gamma)**2>2*K
print("CHECKPOINT ENDPOINTS",2*len(TABLE),"PASS",flush=True)

x,t,U,T=S.symbols("x t U T")
delta=S.Rational(14,5)
y=x+delta*t
power=2*delta*x+delta*delta*t
P=S.expand(
 (1-y)*(1-x*x)*sum(power**h/S.factorial(h) for h in range(13))
 -(1+y)*(1-x*x+t))
counts=0
least=None
for k,K,G in TABLE:
 lo=S.Rational(k,100); hi=S.Rational(k+1,100)
 pp=S.Poly(
  S.expand(P.subs({x:lo+(hi-lo)*U,t:T/S.Integer(K)})),U,T)
 du,dt=pp.degree(U),pp.degree(T)
 terms=dict(pp.terms())
 cc=[]
 for i in range(du+1):
  for j in range(dt+1):
   cc.append(sum(
    z*S.Rational(comb(i,h),comb(du,h))*
      S.Rational(comb(j,l),comb(dt,l))
    for (h,l),z in terms.items() if h<=i and l<=j))
 assert (du,dt,len(cc))==(15,13,224)
 assert min(cc)>0,(k,K,min(cc))
 counts+=len(cc)
 least=min(cc) if least is None else min(least,min(cc))
assert least>S.Rational(1,25)
assert (roots,minor_count,counts)==(930,54863,6272)
print("BE COEFFICIENTS",counts,"MINIMUM > 1/25",flush=True)
assert F(1,5)*F(3,10)*(1-F(9525,10000))<F(1,100)
assert F(1,5)*(1+F(1,15))<1
assert F(99,100)**2-4*F(33,32)**2/F(5)==F(20691,160000)

def bounds(C,X):
 den=1<<80
 z=isqrt(C*den*den)
 ol=F(z,den); oh=F(z+1,den)
 return ol,oh,1-X/ol,1-X/oh

def fourth_bounds(q):
 den=1<<80
 z=isqrt(isqrt(q.numerator*den**4//q.denominator))
 return F(z,den),F(z+1,den)

def data(a,e):
 N=a+e; si=N+2; de=a-e; C=4*(a+1)*(e+1); cc=[1]
 for k in range(N):
  val,rem=divmod(
   de*cc[k]-(N-k+1)*(cc[k-1] if k else 0),k+1)
  assert not rem
  cc.append(val)
 def c(k):
  return cc[k] if 0<=k<=N else 0
 @lru_cache(None)
 def B(k):
  return c(k-1)+c(k+1)
 @lru_cache(None)
 def D(k):
  return c(k)**2-c(k-1)*c(k+1)
 @lru_cache(None)
 def H(k):
  return si*(c(k)**2+c(k-1)**2)-2*de*c(k)*c(k-1)
 return N,si,de,C,c,B,D,H

samples={
 (N-e,e)
 for N in range(80,args.max_n+1)
 for e in range(7,(N-2)//2+1)}
large={
 (2000,250),(4000,250),(3200,400),
 (6400,400),(5600,700),(11200,700)}
samples|=large
counts=[0]*8

for a,e in sorted(samples):
 N=a+e; si=N+2; de=a-e; C=4*(a+1)*(e+1)
 kap=F(C,2*(si+1))
 if C<4096 or kap<15:
  continue
 K=(N+isqrt(C))//2
 while (2*K-N)**2<C:
  K+=1
 if (a,e) in large:
  js=set()
  for p in (705,745,795,845,895,945,975,980,981,985,990,995):
   x0=isqrt(C*p*p//10**6)
   j0=(N+x0)//2
   js.update((j0,j0+1))
 else:
  js=range(N//2+1,K)
 js=[
  j for j in js
  if 100*(2*j-N)**2>49*C
  and 2*(2*j-N)<=si and (2*j-N)**2<C]
 if not js:
  continue
 N,si,de,C,c,B,D,H=data(a,e)
 for j in sorted(js):
  counts[0]+=1
  X=2*j-N
  s=K-j
  V=C-X*X+2*(si-X)
  SV=s*V-2*(X+1)*s*(s-1)-2*s*(s-1)*(2*s-1)//3
  m=None
  early=False
  short=False
  edge=2500*X*X>=2401*C
  if edge:
   ol,oh,hl,hh=bounds(C,X)
   if kap*kap*hh**3<=F(121,100):
    counts[1]+=1
    short=True
    assert 10000*s*SV<24649*de*de
   else:
    assert kap*kap*hl**3>F(121,100)
    counts[2]+=1
    tl=fourth_bounds(hl)[0]
    th=fourth_bounds(hh)[1]
    targetl=X+ol*F(9,4)/th/kap
    targeth=X+oh*F(9,4)/tl/kap
    ml=ceil((N+targetl)/2)
    mh=ceil((N+targeth)/2)
    assert ml==mh
    m=ml
    assert j<m<K
  else:
   k=isqrt(10000*X*X//C)
   if 10000*X*X==k*k*C:
    k-=1
   lo,cut,G=TABLE[k-70]
   assert lo==k
   if kap<cut:
    counts[3]+=1
    short=True
    assert all((-1)**e*c(l)>0 for l in range(j,N+1))
   else:
    counts[4]+=1
    early=True
    gamma=F(G,10000)
    m=j+1
    while gamma.denominator**2*(2*m-N)**2<gamma.numerator**2*C:
     m+=1
    assert m<K

  seen=False
  for i in list(range(j+1,K+1))+[N+1]:
   Y=2*i-N
   W=c(j)*B(i)-B(j)*c(i)
   seen|=W<0 or (W==0 and c(j)*c(i)+B(j)*B(i)<0)
   drop=D(j)-D(i)
   assert drop>=abs(W),(a,e,j,i,"E")
   counts[5]+=1
   if short:
    assert not seen,(a,e,j,i,"SHORT")
   if edge and m is not None and i<=m:
    assert not seen,(a,e,j,i,"EDGE CROSSING")
   if m is not None and i>=m:
    assert C*(si+Y)**2*drop**2>=4*Y*Y*H(j)*H(i),(a,e,j,i,"J")
    counts[6]+=1
   if early and seen and i<m:
    P=C-X*X
    VV=P+2*(si-X)
    bj,bi=comb(si,j+1),comb(si,i+1)
    aa=P*bj-VV*bi
    bb=Y*(P*bj+VV*bi)
    assert aa>=0 and C*aa*aa>=bb*bb
  if edge and seen:
   counts[7]+=1

print("VETTING anchors,edge-short,edge-checkpoint,root-short,"
      "compact-checkpoint,pairs,J,edge-long")
print(counts)

# Definition-level Catalan-moment bridges for the E normalization.
@lru_cache(None)
def mom(power,label):
 if (power+label)%2:
  return 0
 return sum(
  (-1)**h*comb(label-h,h)*
  comb(power+label-2*h,(power+label-2*h)//2)//
  ((power+label-2*h)//2+1)
  for h in range(label//2+1))

@lru_cache(None)
def kernel(e,a,p,q):
 return sum(
  (-1)**h*comb(e,h)*comb(a,k)*
  mom(e+a-h-k,p)*mom(h+k,q)
  for h in range(e+1) for k in range(a+1))

bridges=0
for e in range(3,9):
 for a in range(e+2,e+6):
  N,si,de,C,c,B,D,H=data(a,e)
  j=N//2+2
  i=j+4
  p,q=i+j-N-1,i-j-1
  assert kernel(e,a,p,q)==c(j)*B(i)-B(j)*c(i)
  assert sum(
   kernel(e,a,l,0)
   for l in range(p-q,p+q+1,2))==D(j)-D(i)
  bridges+=1
print("CATALAN BRIDGES",bridges)
print("PASS")