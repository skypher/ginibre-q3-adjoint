"""FM-MECH91 exact certificates and diagnostics.
Run from /home/yang/q3adjoint; writes no files.
"""
import argparse
from math import comb
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--limit",type=int,default=60,
                help="N cutoff for the supplemental exact row census")
args=ap.parse_args()

# Universal outer matrix certificate.
R,t,u,d=ring("t,u,d",QQ)
def covnum(t,u,d):
 s=t+u;X=t-u;den=(t+1)*(u+1)
 L=t*t+u*u+d*d-2;v=d*(s+1)
 a,b,c,f=t*L-v*d,-v*X,d*den,X*den
 return (s*(a*a+b*b)+2*d*a*b,
         s*(a*c+b*f)+d*(a*f+b*c),
         s*(c*c+f*f)+2*d*c*f)

C0=covnum(t,u,d);C1=covnum(t+1,u-1,d)
V=[t*u*u*(t+2)**2*p-u*t*t*(u+1)**2*q
   for p,q in zip(C0,C1)]
assert V[0].gcd(V[1]).gcd(V[2])==t*u
A,B,D=[p.exquo(t*u) for p in V]
det=(A*D-B*B).exquo(t*u*(u+1)**2)

P,l,r,h=ring("l,r,h",QQ)
rr=r+1;ll=l+2
tt=(rr+ll)**2+h;uu=rr**2;dd=2*rr*(rr+ll)+h
def subst(poly):
 powers=[[P.one]+[v**i for i in range(1,poly.degree(z)+1)]
         for v,z in zip((tt,uu,dd),(t,u,d))]
 return sum((c*powers[0][i]*powers[1][j]*powers[2][k]
             for (i,j,k),c in poly.items()),P.zero)

counts=[]
for p in (A,D,det):
 q=subst(p)
 assert q and all(c>0 for c in q.values())
 counts.append(len(q))
assert counts==[907,755,3280]
print("OUTER MATRIX: positive coefficients",counts,flush=True)

# Identities consumed by the proof.
t,u,d,p,q=S.symbols("t u d p q")
s=t+u;X=t-u;C=s*s-d*d;den=(t+1)*(u+1)
cp=(d*p-u*q)/t
cm2=(d*q-(t-1)*p)/(u+1)
cp2=(d*cp-(u-1)*p)/(t+1)
H=s*(p*p+q*q)-2*d*p*q
Hn=s*(p*p+cp*cp)-2*d*p*cp
T=p*p+q*q+cp*cp-p*(cm2+cp2)-cm2*cp2
DT=p*p-q*cp;U=(t*t+u*u+d*d-2)/den
assert S.cancel(H-Hn-X*(s*q-d*p)**2/t**2)==0
assert S.cancel(
 T-U*DT-((X+1)*H+(X-1)*Hn)/(X*den))==0
assert S.expand(4*(t*t+u*u+d*d-2-den)
                -((s-6)*(s+2)+3*X*X+4*d*d))==0
alpha=(t-1)/(u+1)+(u-1)/(t+1)-d*d/(t*(t+1))
beta=-d/(u+1)+d*u/(t*(t+1))
assert S.cancel(-cm2-cp2-alpha*p-beta*q)==0
assert S.cancel(q+cp-d*p/t-X*q/t)==0
r0=S.symbols("r0")
hh=s*(r0*r0+1)-2*d*r0
assert S.cancel(S.diff((s-d*r0)**2/hh,r0)
                +2*C*r0*(s-d*r0)/hh**2)==0
eta=d*(s+1)/(s*s+2*s+d*d)
assert S.cancel(
 -cm2-cp2+2*eta*(q+cp)-U*(p-eta*(q+cp)))==0
print("LOCAL IDENTITIES: PASS",flush=True)

def row(N,e):
 d=N-2*e;c=[1]
 for k in range(N):
  v,rem=divmod(
   d*c[-1]-(N-k+1)*(c[-2] if k else 0),k+1)
  assert rem==0
  c.append(v)
 return c

def data(N,e):
 cc=row(N,e);s=N+2;d=N-2*e
 c=lambda k:cc[k] if 0<=k<=N else 0
 A=[-c(k-2)-c(k+2) for k in range(N+2)]
 B=[c(k-1)+c(k+1) for k in range(N+2)]
 T=[c(k)**2+c(k-1)**2+c(k+1)**2
    -c(k)*(c(k-2)+c(k+2))-c(k-2)*c(k+2)
    for k in range(N+2)]
 H=[s*(c(k)**2+c(k-1)**2)-2*d*c(k)*c(k-1)
    for k in range(N+2)]
 return A,B,T,H

def matrix(N,d,k):
 t=k+1;u=N-k+1;X=t-u;den=(t+1)*(u+1)
 L=t*t+u*u+d*d-2;v=d*(N+3)
 return (t*L-v*d,-v*X,d*den,X*den),t*den

def criterion(N,d,j,i,R,Hj,Hi):
 (a,b,c,f),zj=matrix(N,d,j)
 (g,h,l,m),zi=matrix(N,d,i)
 a,b,c,f=a*l-c*g,a*m-c*h,b*l-f*g,b*m-f*h
 s=N+2;C=s*s-d*d
 tau=(s*s*(a*a+b*b+c*c+f*f)
      +2*s*d*(a*c+b*f+a*b+c*f)+2*d*d*(a*f+b*c))
 P=Hj*Hi;E=C*C*zj*zj*zi*zi*R*R
 return (2*E>=tau*P
         and E*E-tau*P*E+C*C*(a*f-b*c)**2*P*P>=0)

# Original curve short, combined curve long.
aa,bb,tt,hh=data(14,2);cc=row(14,2)
assert min(cc[9:15])>0
assert aa[9]*bb[14]-bb[9]*aa[14]==-252
assert all(aa[k]*bb[k+1]-bb[k]*aa[k+1]>0
           for k in range(9,14))
print("PHASE TRANSFER OBSTRUCTION: PASS",flush=True)

# Supplemental census, including i=N and i=N+1.
pairs=longs=outer=0
for N in range(14,args.limit+1):
 for e in range(2,(N+1)//2):
  d=N-2*e;C=(N+2)**2-d*d
  A,B,T,H=data(N,e)
  for k in range(N//2+1,N+1):
   if (2*k-N)**2>=C:
    assert (k+1)*H[k+1]<=(N-k+1)*H[k]
    outer+=1
  for j in range(N//2+1,N-3):
   X=2*j-N
   if X>2 and (X-2)**2*(e+1)**2>=(e-1)**2*C:
    continue
   seen=False
   for i in range(j+1,N+2):
    W=A[j]*B[i]-B[j]*A[i]
    seen|=W<0 or (W==0 and A[j]*A[i]+B[j]*B[i]<0)
    if i-j<5:
     continue
    R=T[j]-T[i]
    assert R>=abs(W)
    pairs+=1
    if seen:
     assert criterion(N,d,j,i,R,H[j],H[i])
     longs+=1
print("ROW CENSUS: pairs,long,outer steps",
      pairs,longs,outer,flush=True)
if args.limit==60:
 assert (pairs,longs,outer)==(130729,125016,2438)

# Exact example of the uniform tail sector.
N,e,j,l=1000,100,501,700;s=N+2;d=N-2*e
assert 0<2*j-N<2*l-N and (2*l-N)**2<s*s-d*d
assert 16*s**9*comb(s,l+1)<=comb(s,j+1)
A,B,T,H=data(N,e)
assert 4*s**8*T[l]<=T[j]
for i in range(l,N+2):
 assert criterion(N,d,j,i,T[j]-T[i],H[j],H[i])
 assert 2*(T[j]-T[i]-abs(A[j]*B[i]-B[j]*A[i]))>=T[j]
print("TAIL SECTOR: 302 endpoints PASS",flush=True)

# Corrected b=2 energy, with scalar channel relations imposed.
N,e,b,j,i=150,2,2,96,100
d=N-2*e;cc=row(N,e);w=(4,4,1)

def vc(h,k):
 return sum(
  comb(h,r)*(cc[k+h-2*r] if 0<=k+h-2*r<=N else 0)
  for r in range(h+1))

def TB(k):
 return sum(w[h]*(vc(h,k)**2-vc(h,k-1)*vc(h,k+1))
            for h in range(3))

def EB(k):
 P=[vc(h,k)+vc(h,k-1) for h in range(3)]
 Q=[vc(h,k)-vc(h,k-1) for h in range(3)]
 return (
  sum(w[h]*((e+1)*P[h]**2+(N-e+1)*Q[h]**2)
      for h in range(3))
  +sum(w[h]*(2-h)*(
    (2*P[h]-P[h+1])**2+(2*Q[h]+Q[h+1])**2)//2
       for h in range(2)))

def pullback(k):
 cs={k:S.Matrix([1,0]),k-1:S.Matrix([0,1])}
 for l in range(k,k+3):
  cs[l+1]=(d*cs[l]-(N-l+1)*cs[l-1])/(l+1)
 for l in range(k-1,k-3,-1):
  cs[l-1]=(d*cs[l]-(l+1)*cs[l+1])/(N-l+1)
 def v(h,l):
  return sum((comb(h,r)*cs[l+h-2*r]
              for r in range(h+1)),S.zeros(2,1))
 V=[v(h,k) for h in range(4)]
 P=[v(h,k)+v(h,k-1) for h in range(3)]
 Q=[v(h,k)-v(h,k-1) for h in range(3)]
 G=sum((w[h]*(
  (e+1)*P[h]*P[h].T+(N-e+1)*Q[h]*Q[h].T)
        for h in range(3)),S.zeros(2))
 for h in range(2):
  pp=2*P[h]-P[h+1];qq=2*Q[h]+Q[h+1]
  G+=S.Rational(w[h]*(2-h),2)*(pp*pp.T+qq*qq.T)
 assert G[0,0]>0 and G.det()>0
 return G,V

Gj,Vj=pullback(j);Gi,Vi=pullback(i)
Q=sum((w[h]*(Vj[h]*Vi[h+1].T-Vj[h+1]*Vi[h].T)
       for h in range(3)),S.zeros(2))
zj=S.Matrix([cc[j],cc[j-1]])
zi=S.Matrix([cc[i],cc[i-1]])
assert (zj.T*Gj*zj)[0]==EB(j)
assert (zi.T*Gi*zi)[0]==EB(i)
W=sum(w[h]*(vc(h,j)*vc(h+1,i)-vc(h+1,j)*vc(h,i))
      for h in range(3))
assert (zj.T*Q*zi)[0]==W
R=TB(j)-TB(i)
tau=(Gj.inv()*Q*Gi.inv()*Q.T).trace()
ratio=S.cancel(2*R*R/(tau*EB(j)*EB(i)))
assert 0<ratio<S.Rational(1,10)
assert ratio==S.Rational(
 1617024255773736151655927508203784264127587702624130035060389668901,
 19462293570447534407236891093722356673792884741856061270982167191936)
assert R>abs(W)
print("b2 ENDPOINT-NORM OBSTRUCTION:",ratio,flush=True)
print("b2 ACTUAL SLACK:",R-abs(W),flush=True)
print("PASS",flush=True)