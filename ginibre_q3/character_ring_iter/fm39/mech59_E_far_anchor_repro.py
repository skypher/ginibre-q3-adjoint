"""FM-MECH59: far-anchor (E), exact certificates and Catalan bridges."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n",type=int,default=90)
args=ap.parse_args()
if args.max_n<53:
    ap.error("--max-n must be at least 53 for the finite proof remainder.")

# Exact Sonin identity; omega^2=sigma^2-Delta^2 makes the remainder a square.
sig,X,de,omega,p,q=S.symbols("sig X de omega p q")
cn=(2*de*p-(sig-X)*q)/(sig+X)
cnn=(2*de*cn-(sig-X-2)*p)/(sig+X+2)
D=p*p-q*cn
Dn=cn*cn-p*cnn
rate=(sig-X)*(omega+X+2)/((sig+X+2)*(omega+X))
rem=2*((sig-omega)*cn**2-2*de*cn*p+(sig+omega)*p*p)/(
    (sig+X+2)*(omega+X))
assert S.cancel(rate*D-Dn-rem)==0
assert S.expand((sig-omega)*(sig+omega)-de**2
                -(sig**2-omega**2-de**2))==0

z,v,w=S.symbols("z v w",nonnegative=True)
P=3*X*(sig+X+2)-8*(sig-X)*(X+2)
assert S.expand(4*P.subs(X,(sig+1)/2+z)
    -((sig-4)**2+39+24*sig*z+44*z*z+132*z))==0
P=3*X*(sig+X+2)-10*(sig-X)*(X+2)
pc=S.Poly(S.expand(25*P.subs(X,3*sig/5+w).subs(sig,10+v)),v,w)
assert all(c>0 for c in pc.coeffs())
assert pc.as_expr()==12*v*v+215*v*w+130*v+325*w*w+2800*w+100
rho=F(3,8)
assert 1-rho**5-5*(1+rho)*rho**2==F(845,32768)
assert F(6,5)**2*rho<1
rho=F(3,10)
assert (1-rho**4)**2-16*(1+rho)**2*rho**3==F(25378561,100000000)
print("SONIN IDENTITY, RATE BOUNDS, LD SCALARS PASS")

# Root-free coverage above the recurrence discriminant threshold.
O=sig**2-X**2-6*X-2*sig-8
G=S.expand(25*O-(3*X+2*sig-8)**2)
assert S.diff(G,X)==-68*X-12*sig-102
assert S.expand(25*G.subs(X,3*sig/5).subs(sig,56+v))==(
    39*v*v+2388*v+4824)
print("ROOT-FREE COVERAGE POLYNOMIAL PASS")

# Gap four, n=N-j>=4, X>=2n+3, Delta^2<=4(n-1)(n+X+2).
R,sig,X,de,n,u,v=ring("sig,X,de,n,u,v",QQ)
def sq(t):return [t[0]*t[0],2*t[0]*t[1],t[1]*t[1]]
def cross(t,w):return [
    t[0]*w[0],t[0]*w[1]+t[1]*w[0],t[1]*w[1]]
def add(*xs):return [sum(z,R.zero) for z in zip(*xs)]
def scale(a,xs):return [a*z for z in xs]
Y=X+8
ns=[(R.one,R.zero),(2*de,-(sig-X))]
den=sig+X
for h in range(1,5):
    fac=(sig-X-2*h)*(sig+X+2*h-2)
    ns.append(tuple(2*de*a-fac*b for a,b in zip(ns[-1],ns[-2])))
    if h<4:den*=sig+X+2*h
fac=sig+X+6
Dh=add(scale(sig+Y,sq(ns[4])),scale(-fac,cross(ns[3],ns[5])))
T=add(scale((sig+Y)*den*den,[sig+X,-2*de,sig-X]),
      scale(-(sig+X),Dh))
Bi=tuple(fac*(sig+Y)*a+b for a,b in zip(ns[3],ns[5]))
W=add(scale((sig+X)*den,cross((R.one,R.zero),Bi)),
      scale(-(sig+Y)*den,cross((2*de,2*X),ns[4])))
A,B,C=add(T,W)
minus=add(T,scale(-1,W))
assert minus[0]==A.compose(de,-de)
assert minus[1]==-B.compose(de,-de)
assert minus[2]==C.compose(de,-de)
det=4*A*C-B*B
def evenodd(P):
    pp=P.compose(de,-de)
    return (P+pp)/2,(P-pp)/2
Ae,Ao=evenodd(A)
De,Do=evenodd(det)
tests=[("Aeven",Ae,(4,330,QQ(1024))),
       ("Aproduct",Ae*Ae-Ao*Ao,(8,2079,QQ(1048576))),
       ("deteven",De,(4,1130,QQ(1048576))),
       ("detproduct",De*De-Do*Do,(8,7605,QQ(4398046511104,7)))]
bound=4*(n-1)*(n+X+2)
for label,pol,expected in tests:
    assert all(mon[2]%2==0 for mon in pol)
    degree=max(mon[2]//2 for mon in pol)
    powers=[R.zero for _ in range(degree+1)]
    for mon,co in pol.items():
        h=mon[2]//2
        mm=list(mon)
        mm[2]=0
        powers[h]+=R.from_dict({tuple(mm):co})
    bb=bound.compose(X,2*n+3+u).compose(n,4+v)
    for h in range(degree+1):
        powers[h]=(powers[h].compose(sig,2*n+X+2)
                   .compose(X,2*n+3+u).compose(n,4+v)*bb**h)
    count=0
    minimum=None
    for j in range(degree+1):
        P=sum((powers[h]*QQ(comb(j,h),comb(degree,h))
               for h in range(j+1)),R.zero)
        cs=list(P.values())
        assert all(c>0 for c in cs),(label,j)
        count+=len(cs)
        minimum=min(cs) if minimum is None else min(minimum,min(cs))
    assert (degree,count,minimum)==expected
    print("MATRIX CERTIFICATE",label,degree,count,minimum,flush=True)

# Independent recurrence evaluation of the cleared quadratic form.
def evaluate(P,si,xx,dd):
    total=F(0)
    for mon,co in P.items():
        assert mon[3:]==(0,0,0)
        total+=F(int(co.numerator),int(co.denominator))*(
            si**mon[0]*xx**mon[1]*dd**mon[2])
    return total
identities=0
for nn in range(4,7):
    for extra in (0,2):
        xx=2*nn+3+extra
        si=2*nn+xx+2
        jj=nn+xx
        for dd in (0,1,2*nn):
            aa,bb,cc=[evaluate(P,si,xx,dd) for P in (A,B,C)]
            denominator=(si+xx)*(si+xx+8)
            for h in range(4):
                denominator*=(si+xx+2*h)**2
            for pp,qq in ((1,0),(0,1),(1,1)):
                c={-1:F(qq),0:F(pp)}
                for h in range(5):
                    c[h+1]=(dd*c[h]-(nn-h+1)*c[h-1])/(jj+h+1)
                dv=lambda h:c[h]**2-c[h-1]*c[h+1]
                bv=lambda h:c[h-1]+c[h+1]
                value=dv(0)-dv(4)+c[0]*bv(4)-bv(0)*c[4]
                assert aa*pp*pp+bb*pp*qq+cc*qq*qq==denominator*value
                identities+=1
print("QUADRATIC IDENTITY CHECKS",identities)

def data(a,e):
    N=a+e
    sig,de=N+2,a-e
    C=sig*sig-de*de
    c=[1]
    for k in range(N):
        z,rem=divmod(de*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
        assert rem==0
        c.append(z)
    val=lambda k:c[k] if 0<=k<=N else 0
    B=[val(k-1)+val(k+1) for k in range(N+2)]
    D=[val(k)**2-val(k-1)*val(k+1) for k in range(N+2)]
    H=[sig*(val(k)**2+val(k-1)**2)-2*de*val(k)*val(k-1)
       for k in range(N+2)]
    return sig,de,C,val,B,D,H

counts=dict(largegap=0,rate4=0,matrix4=0,RF4=0,finite4=0,terminal4=0)
smallmin=None
checks=0
for N in range(8,args.max_n+1):
    for e in range(3,(N-2)//2+1):
        a=N-e
        sig,de,C,val,B,D,H=data(a,e)
        for j in range(N//2+1,N-2):
            X=2*j-N
            if 2*X<=sig:continue
            assert 8*D[j+1]<=3*D[j]
            if 5*X>=3*sig:assert 10*D[j+1]<=3*D[j]
            for i in range(j+4,N+2):
                s=i-j
                W=val(j)*B[i]-B[j]*val(i)
                slack=D[j]-D[i]-abs(W)
                assert slack>=0,(a,e,j,i)
                checks+=1
                if s>=5:
                    kind="largegap"
                elif i==N+1:
                    kind="terminal4"
                    z=abs(val(j))
                    assert 4*(N-2)*D[j]>=(N+1)*z*z
                    assert D[j]>=z
                elif 5*X>=3*sig:
                    kind="rate4"
                else:
                    nn=N-j
                    threshold=4*(nn-1)*(j+2)
                    if de*de<=threshold:
                        kind="matrix4"
                    elif sig>=56:
                        kind="RF4"
                        assert X*X>=4*(e-1)*(a+4)
                    else:
                        kind="finite4"
                        assert N<=53
                        item=slack,(a,e,j,i)
                        if smallmin is None or item<smallmin:smallmin=item
                counts[kind]+=1
assert counts["finite4"]==100
assert smallmin==(719376,(19,3,18,22))
print("FAR COVERAGE",checks,counts)
print("FINITE BOX MINIMUM",smallmin)

# Independent definition-level Catalan evaluation.
@lru_cache(None)
def moment(power,label):
    if (power+label)%2:return 0
    return sum((-1)**h*comb(label-h,h)
               *comb(power+label-2*h,(power+label-2*h)//2)
               //((power+label-2*h)//2+1)
               for h in range(label//2+1))
@lru_cache(None)
def kernel(e,a,p,q):
    return sum((-1)**h*comb(e,h)*comb(a,k)
               *moment(e+a-h-k,p)*moment(h+k,q)
               for h in range(e+1) for k in range(a+1))
bridges=0
for N in range(8,23):
    for e in range(3,(N-2)//2+1):
        a=N-e
        sig,de,C,val,B,D,H=data(a,e)
        for j in range(N//2+1,N-2):
            if 2*(2*j-N)<=sig:continue
            for i in range(j+4,N+2):
                p,q=i+j-N-1,i-j-1
                W=val(j)*B[i]-B[j]*val(i)
                assert kernel(e,a,p,q)==W
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2))==D[j]-D[i]
                bridges+=1
print("CATALAN REDUCTION PAIRS",bridges)

# A remaining unbalanced central long arc: the separate-chord bound fails.
a,e,j,i=4000,4,2016,2051
sig,de,C,val,B,D,H=data(a,e)
X,Y=2*j-a-e,2*i-a-e
assert 4*de>sig and 2*X<=sig and Y*Y<C
assert any(val(j)*B[k]-B[j]*val(k)<0 for k in range(j+1,i+1))
s=i-j
drop=D[j]-D[i]
assert C*(sig+Y)**2*drop**2>=4*Y*Y*H[j]*H[i]
assert 10*drop*drop<s*s*D[j]*D[i-1]
print("UNBALANCED WITNESS",(a,e,j,i),"J>=0; LD fails by a factor >sqrt(10)")
print("PASS")

