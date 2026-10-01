import argparse
from fractions import Fraction as F
from math import comb
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap=argparse.ArgumentParser(description="FM-MECH84 exact certificates and bounded searches.")
ap.add_argument("--limit",type=int,default=100)
ap.add_argument("--b2-limit",type=int,default=140)
ap.add_argument("--large",action="store_true")
args=ap.parse_args()

# Universal local identities.
t,u,d,p,q=S.symbols("t u d p q")
sig=t+u; X=t-u; C=sig**2-d**2; den=(t+1)*(u+1)
cp=(d*p-u*q)/t
cm2=(d*q-(t-1)*p)/(u+1)
cp2=(d*cp-(u-1)*p)/(t+1)
D=p*p-q*cp
T=p*p+q*q+cp*cp-p*(cm2+cp2)-cm2*cp2
U=(t*t+u*u+d*d-2)/den
H0=sig*(p*p+q*q)-2*d*p*q
H1=sig*(p*p+cp*cp)-2*d*p*cp
K=t*t+t+u*u+u
assert S.cancel(T-U*D-(2*t*t*C*p*p+K*(sig*q-d*p)**2)/(sig*t*t*den))==0
assert S.cancel(T-U*D-((X+1)*H0+(X-1)*H1)/(X*den))==0
assert S.cancel(H0-H1-X*(sig*q-d*p)**2/t**2)==0
cc=dict(zip(range(-2,5),S.symbols("c0:7")))
def AB(c,k): return (-c[k-2]-c[k+2],c[k-1]+c[k+1])
def tv(c,k):
 return c[k]**2+c[k-1]**2+c[k+1]**2-c[k]*(c[k-2]+c[k+2])-c[k-2]*c[k+2]
def wd(c,j,i):
 a,b=AB(c,j);c0,d0=AB(c,i)
 return a*d0-b*c0
assert S.expand(tv(cc,0)-tv(cc,1)-wd(cc,0,1))==0
for eps in (-1,1):
 gg={k:cc[k]+eps*cc[k-1] for k in range(-1,5)}
 assert S.expand(tv(gg,1)-tv(gg,2)-tv(cc,0)+tv(cc,2)-eps*wd(cc,0,2))==0
print("local identities PASS",flush=True)

# Formal normalized Sonin certificate.
R,n,x,d,o=ring("n,x,d,o",QQ);j=n+x;sig=2*n+x+2
de={-2:(n+1)*(n+2),-1:n+1,0:R.one,1:R.one,2:j+2}
nu={-2:(d*d-j*(n+1),-d*(j+1)),-1:(d,-j-1),
    0:(R.one,R.zero),1:(R.zero,R.one),2:(-n,d)}
for h in (2,3):
 de[h+1]=de[h]*(j+h+1)
 nu[h+1]=tuple(d*nu[h][l]-(n-h+1)*(j+h)*nu[h-1][l] for l in (0,1))
H=(n+1)**2*(n+2)*(j+2)**2*(j+3)
def dot(h,l):
 z,r=divmod(H,de[h]*de[l]);assert not r
 a,b=nu[h];c,f=nu[l]
 return [z*a*c,z*(a*f+b*c),z*b*f]
def tv0(k):
 out=[R.zero]*3
 for a,b,sg in ((0,0,1),(-1,-1,1),(1,1,1),(-2,0,-1),(0,2,-1),(-2,2,-1)):
  out=[z+sg*w for z,w in zip(out,dot(k+a,k+b))]
 return out
def red(poly):
 out=R.zero
 for (a,b,h,z),v in poly.items():
  out+=v*n**a*x**b*d**(h%2)*o**z*(sig*sig-o*o)**(h//2)
 return out
K0=(j+1)**2+(n+1)**2+d*d-2; L0=(j+2)*(n+2)
K1=(j+2)**2+n*n+d*d-2; L1=(j+3)*(n+1)
v=[red(a*L0*K1*(n+1)*(o+x+2)-b*L1*K0*(j+2)*(o+x))
   for a,b in zip(tv0(0),tv0(1))]
gc=v[0].gcd(v[1]).gcd(v[2])
assert gc==(n+1)*(n+2)*(j+2)*(j+3)
A,B,C=[z.exquo(gc) for z in v]
FD,rem=divmod(red(4*A*C-B*B),(n+1)**2*o**2*(j+2)**2)
assert not rem
P,u,vv=ring("u,v",QQ)
totals=[]
for name,poly in (("A",A),("C",C),("det",FD)):
 deg=poly.degree(o);parts=[];ss=2*u+vv+3
 for h in range(deg+1):
  parts.append(sum((v*u**i*(vv+1)**k for (i,k,dd,l),v in poly.items()
                    if l==h and dd==0),P.zero)*ss**h)
 assert all(dd==0 for i,k,dd,l in poly)
 total=0
 for k in range(deg+1):
  z=sum((QQ(comb(k,h),comb(deg,h))*parts[h] for h in range(k+1)),P.zero)
  assert z and all(v>0 for v in z.values())
  total+=len(z)
 totals.append(total)
assert totals==[213,214,298]
print("Sonin positive coefficients",totals,flush=True)

# Homogeneous resultant excludes consecutive flat steps.
H=(n+1)**2*(n+2)*(j+2)**2*(j+3)**2*(j+4)**2
def wedge(k):
 out=[R.zero]*3
 for h in (-2,2):
  for l in (0,2):out=[a-b for a,b in zip(out,dot(k+h,k+l))]
 for h in (-1,1):
  for l in (-1,3):out=[a+b for a,b in zip(out,dot(k+h,k+l))]
 g=out[0].gcd(out[1]).gcd(out[2])
 assert g==((j+3)*(j+4)**2 if k==0 else (n+1)*(n+2)*(j+4))
 return [z.exquo(g) for z in out]
A,B,C=wedge(0);D,E,F0=wedge(1)
res=(A*F0-C*D)**2-(A*E-B*D)*(B*F0-C*E)
N=2*n+x
G=2*n*n+2*n*x+4*n+x*x+4*x+2
F1=-2*n*n-2*n*x-8*n+x*x+2*x
P4=6*n*n+6*n*x+20*n+5*x*x+24*x+16
P2=(12*n**4+24*n**3*x+80*n**3+20*n*n*x*x+128*n*n*x+164*n*n
    +8*n*x**3+88*n*x*x+244*n*x+160*n+5*x**4+52*x**3+166*x*x+180*x+64)
target=((n+1)**2*((N+2)**2-d*d)*((N+4)**2-d*d)*(j+3)**2
        *(x+2)**2*(j+2)**4*(d*d+G)*(d**6+P4*d**4+P2*d*d+F1*F1*G))
assert res==target
assert all(v>0 for poly in (G,P4,P2) for v in poly.values())
print("flat resultant PASS",flush=True)

# Root-bound algebra and balanced Jacobi characteristic polynomial.
a,e,l=S.symbols("a e l",positive=True)
entry=l*(a+e-l+1)/(4*(a+1)*(e+1))
assert S.cancel(entry-l*(1+(e-l)/(a+1))/(4*(e+1)))==0
xx=S.Symbol("x")
for e0 in range(1,21):
 p0,p1=S.Integer(1),xx
 for h in range(2,e0+1):
  p0,p1=p1,S.expand(xx*p1-(h-1)*(2*e0-h+2)*p0)
 assert S.expand(p1-S.prod(xx-(2*e0-2-4*r) for r in range(e0)))==0
NN,xx0=S.symbols("N X")
tt=(NN+xx0+2)/2;uu=(NN-xx0+2)/2
assert S.cancel((tt-1)*(xx0-2)/uu+(uu-1)*(xx0+2)/tt
                -xx0*(NN*NN-2*NN-4+xx0*xx0)/(2*tt*uu))==0
print("root-bound algebra PASS",flush=True)

# Exact actual-row searches.
def row(a,e):
 N=a+e;out=[1]
 for k in range(N):
  z,r=divmod((a-e)*out[-1]-(N-k+1)*(out[-2] if k else 0),k+1)
  assert r==0
  out.append(z)
 return out
def mat(N,d,k):
 t=k+1;u=N-k+1;X=2*k-N
 den=(t+1)*(u+1);L=t*t+u*u+d*d-2;v=d*(N+3)
 return (t*L-v*d,-v*X,d*den,X*den),t*den
def energy(N,d,j,i,R0,Hj,Hi,ms):
 (a,b,c,f),zj=ms[j];(g,h,l,m),zi=ms[i]
 a,b,c,f=a*l-c*g,a*m-c*h,b*l-f*g,b*m-f*h
 sig=N+2;C=sig*sig-d*d
 tau=sig*sig*(a*a+b*b+c*c+f*f)+2*sig*d*(a*c+b*f+a*b+c*f)+2*d*d*(a*f+b*c)
 qdet=a*f-b*c;Hp=Hj*Hi;E=C*C*zj*zj*zi*zi*R0*R0
 return 2*E>=tau*Hp and E*E-tau*Hp*E+C*C*qdet*qdet*Hp*Hp>=0
short=long=0
for N in range(14,args.limit+1):
 for e in range(1,(N+1)//2):
  d=N-2*e;cc=row(N-e,e)
  c=lambda k:cc[k] if 0<=k<=N else 0
  aa=[-c(k-2)-c(k+2) for k in range(N+2)]
  bb=[c(k-1)+c(k+1) for k in range(N+2)]
  tt=[c(k)**2+c(k-1)**2+c(k+1)**2-c(k)*(c(k-2)+c(k+2))-c(k-2)*c(k+2) for k in range(N+2)]
  hh=[(N+2)*(c(k)**2+c(k-1)**2)-2*d*c(k)*c(k-1) for k in range(N+2)]
  ms=[mat(N,d,k) for k in range(N+2)]
  for j in range((N+1)//2,N-3):
   if 2*j<=N:continue
   seen=False
   for i in range(j+1,N+2):
    W=aa[j]*bb[i]-bb[j]*aa[i]
    seen|=W<0 or (W==0 and aa[j]*aa[i]+bb[j]*bb[i]<0)
    if i-j<5:continue
    R0=tt[j]-tt[i]
    assert R0>=abs(W)
    if seen:
     long+=1
     assert energy(N,d,j,i,R0,hh[j],hh[i],ms)
    else:short+=1
print("b1 short,long",short,long,"no failures",flush=True)
if args.limit==100:assert (short,long)==(194508,1111266)

if args.large:
 from math import isqrt
 count=long=fail=0
 for N in (200,500,1000,2000):
  for e in sorted({1,2,3,6,10,30,N//8,N//4,N//2-1}):
   if e<=0 or 2*e>=N:continue
   d=N-2*e;cc=row(N-e,e);C=(N+2)**2-d*d
   c=lambda k:cc[k] if 0<=k<=N else 0
   aa=[-c(k-2)-c(k+2) for k in range(N+2)]
   bb=[c(k-1)+c(k+1) for k in range(N+2)]
   tt=[c(k)**2+c(k-1)**2+c(k+1)**2-c(k)*(c(k-2)+c(k+2))-c(k-2)*c(k+2) for k in range(N+2)]
   hh=[(N+2)*(c(k)**2+c(k-1)**2)-2*d*c(k)*c(k-1) for k in range(N+2)]
   ms=[mat(N,d,k) for k in range(N+2)]
   om=isqrt(C)
   js=sorted({N//2+1,*((N+(v*om)//100)//2 for v in (5,20,40,60,75,85,95,99,105))})
   for j in js:
    if 2*j<=N or j+5>N-1:continue
    seen=False
    for i in range(j+1,N):
     W=aa[j]*bb[i]-bb[j]*aa[i]
     seen|=W<0 or (W==0 and aa[j]*aa[i]+bb[j]*bb[i]<0)
     if i-j<5 or (i-j>40 and (i-j)%max(1,N//100)):continue
     R0=tt[j]-tt[i];assert R0>=abs(W)
     count+=1
     if seen:
      long+=1
      fail+=not energy(N,d,j,i,R0,hh[j],hh[i],ms)
  print("large N",N,"pairs",count,"long",long,"failures",fail,flush=True)
 assert fail==0

count=0
for N in range(10,args.b2_limit+1):
 for e in range(1,N//2+1):
  if 2*e==N:continue
  cc=row(N-e,e);c=lambda k:cc[k] if 0<=k<=N else 0
  v=[{k:sum(comb(h,r)*c(k+h-2*r) for r in range(h+1))
      for k in range(-1,N+4)} for h in range(4)]
  def TB(k):
   return sum(comb(2,h)*2**(2-h)*(v[h][k]**2-v[h][k-1]*v[h][k+1]) for h in range(3))
  for j in range(N//2+1,N-3):
   i=j+4;R0=TB(j)-TB(i)
   W=sum(comb(2,h)*2**(2-h)*(v[h][j]*v[h+1][i]-v[h+1][j]*v[h][i]) for h in range(3))
   assert R0>abs(W)
   count+=1
print("b2 gap4 positive",count,flush=True)
if args.b2_limit==140:assert count==209304

# Failed free-state extension at b=4, distance 7.
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sc(a,x):return tuple(x*z for z in a)
def mul(a,b):return a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[1]*b[1]
def local(N,d,k,b):
 t=k+1;u=N-k+1;c={-1:(F(0),F(1)),0:(F(1),F(0))}
 for h in range(1,b+1):
  c[-h-1]=sc(add(sc(c[-h],d),sc(c[-h+1],-(t-h))),F(1,u+h))
 for h in range(b+2):
  c[h+1]=sc(add(sc(c[h],d),sc(c[h-1],-(u-h))),F(1,t+h))
 def v(h,k0):
  z=(F(0),F(0))
  for r in range(h+1):z=add(z,sc(c[k0+h-2*r],comb(h,r)))
  return z
 def tb(k0):
  z=(F(0),)*3
  for h in range(b+1):
   z=add(z,sc(add(mul(v(h,k0),v(h,k0)),sc(mul(v(h,k0-1),v(h,k0+1)),-1)),comb(b,h)*2**(b-h)))
  return z
 A,B,C=tb(0)
 return tb(0),tb(1),-F(t,d*(2*k-N))*((N+2)*B+2*d*C)
T0,T1,U0=local(48,40,45,4);_,_,U1=local(48,40,46,4)
A,B,C=add(sc(T0,F(37,423)/U0),sc(T1,-1/U1))
r=F(13,9);w=A+B*r+C*r*r
assert w==-F(727416812969514908039917081,134084600436422613358531417928310)
cc=row(44,4)
assert A*cc[45]**2+B*cc[45]*cc[44]+C*cc[44]**2>0
print("b4 free-state obstruction",w,flush=True)

# Independent semicircle-moment consumer checks.
def pmul(p,q):
 r={}
 for (i,j),v in p.items():
  for (k,l),w in q.items():r[i+k,j+l]=r.get((i+k,j+l),0)+v*w
 return {m:v for m,v in r.items() if v}
def moment(a,e,b,n,m,eps):
 p={(0,0):1}
 for count,fac in ((a,{(1,0):1,(0,1):1}),(e,{(1,0):1,(0,1):-1}),
                   (b,{(2,0):1,(0,2):1,(0,0):-2})):
  for _ in range(count):p=pmul(p,fac)
 for k,sg in ((n,eps),(m,eps*(-1)**e)):
  q={}
  for h in range(k//2+1):
   v=(-1)**h*comb(k-h,h);power=k-2*h
   q[power,0]=q.get((power,0),0)+v
   q[0,power]=q.get((0,power),0)+sg*v
  p=pmul(p,q)
 return sum(v*(comb(i,i//2)//(i//2+1))*(comb(j,j//2)//(j//2+1))
            for (i,j),v in p.items() if i%2==j%2==0)
def consumer(a,e,b,n,m):
 N=a+e;cs=row(a,e)
 c=lambda k:cs[k] if 0<=k<=N else 0
 v=lambda h,k:sum(comb(h,r)*c(k+h-2*r) for r in range(h+1))
 t=lambda k:sum(comb(b,h)*2**(b-h)*(v(h,k)**2-v(h,k-1)*v(h,k+1)) for h in range(b+1))
 j=(N+n-m)//2;i=j+m+1
 R0=t(j)-t(i)
 W=sum(comb(b,h)*2**(b-h)*(v(h,j)*v(h+1,i)-v(h+1,j)*v(h,i)) for h in range(b+1))
 return sorted([2*(R0-W),2*(R0+W)])
for tup in ((10,5,1,5,4),(8,6,1,6,4),(3,10,1,5,4),(9,6,2,4,3)):
 a,e,b,n,m=tup
 vals=sorted(moment(a,e,b,n,m,eps) for eps in (-1,1))
 assert vals==consumer(*tup)
 print("consumer",tup,vals)
print("ALL CLAIMED CERTIFICATES AND BOUNDED TESTS PASS",flush=True)