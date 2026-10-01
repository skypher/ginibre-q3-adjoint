"""FM-MECH73: central anchors, exact certificates and row vetting.
Run from the repository root with python3 -u. No file writes.
"""
import argparse
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-a",type=int,default=220,
                help="bounded diagnostic grid; the proof uses the symbolic certificates")
args=ap.parse_args()

# Exact transfer and the fixed map to psi.
si,de,om,X,h,l=S.symbols("si de om X h l")
p=h/om
q=(l+de*p)/si
pn=(2*de*p-(si-X)*q)/(si+X)
rr=(si-X)/(si+X)
assert S.cancel(om*pn-(de*h-rr*om*l)/si)==0
err=S.cancel(si*p-de*pn-(om*h+rr*de*l)/si)
assert S.cancel(err.subs(om**2,si**2-de**2))==0
assert S.cancel(q+pn-2*de*p/si-2*X*l/(si*(si+X)))==0
print("TRANSFER IDENTITIES PASS")

# Rational lower sine bounds and pi > 157/50.
sin3=lambda z:z-z**3/6
assert 16*(F(1,5)-F(1,375))-F(4,239)>F(157,50)
assert sin3(F(11,15))>F(2,3)
assert sin3(F(17,16)*F(11,20))>F(11,20)
assert sin3(F(53,50)*F(8,15))>F(8,15)
assert sin3(F(101,100)*F(68,375))>F(68,375)
assert sin3(F(5,8))>F(7,12)
assert sin3(F(21,80))>F(7,27)
assert F(512,195)>F(5,2)
assert F(97,96)*F(43,80)<F(11,20)
assert F(65,64)*F(21,40)<F(8,15)
R=F(101,100)
anchor=F(97,96)*(F(3,4)*F(11,10)+R*F(3,4)**2/30)
small=F(97,96)*(F(13,5)*F(11,10)+R*F(13,5)**2/30)
middle=F(1291,1280)*(F(53,20)*F(17,16)+R*F(53,20)**2/30)
assert anchor<F(157,100)
assert small<F(157,50) and middle<F(157,50)
angles=[]
for z,dist in ((F(3,4),F(5,8)),(F(3),F(21,80))):
    bound=F(121,120)*F(21,10)*(
        F(53,50)+R*(2*z+F(21,10))/30)+dist
    assert bound<F(157,50)
    angles.append(bound)
print("ANGLE BOUNDS",anchor,small,middle,angles)

def bernstein(poly):
    ds=tuple(poly.degree(z) for z in poly.gens)
    a={idx:F(co) for idx,co in poly.terms()}
    for ax,d in enumerate(ds):
        b={}
        for idx,co in a.items():
            h=idx[ax]
            for i in range(h,d+1):
                ix=list(idx);ix[ax]=i;ix=tuple(ix)
                b[ix]=b.get(ix,F(0))+co*F(comb(i,h),comb(d,h))
        a=b
    return ds,[a.get(ix,F(0))
               for ix in product(*(range(d+1) for d in ds))]

# Four entire rectangles, with no subdivision.
u,v,t,z=S.symbols("u v t z")
tests=(
    ("gap5/2",F(0),F(3),"gap",F(5,2),132),
    ("central13/5",F(0),F(1,2),"end",F(13,5),209),
    ("central53/20",F(1,2),F(3,4),"end",F(53,20),209),
    ("gap21/10",F(3,4),F(3),"gap",F(21,10),132))
for name,lo,hi,kind,c,count in tests:
    c=S.Rational(c.numerator,c.denominator)
    w=z+c if kind=="gap" else c
    x=z*t;y=w*t;lam=t*(w*w-z*z)
    poly=S.expand(
        (1-y)*(1-x*x)*sum(lam**h/factorial(h) for h in range(9))
        -(1+y)*(1-x*x+t))
    assert poly.subs(t,0)==0
    poly=S.cancel(poly/t)
    poly=S.Poly(S.expand(poly.subs({
        z:S.Rational(lo.numerator,lo.denominator)+
          S.Rational((hi-lo).numerator,(hi-lo).denominator)*u,
        t:v/15})),u,v)
    ds,coeffs=bernstein(poly)
    assert len(coeffs)==count and min(coeffs)>=0
    print("BE CERTIFICATE",name,"DEGREES",ds,
          "COEFFICIENTS",len(coeffs),"ALL NONNEGATIVE")

def row(a,e):
    N=a+e;sig=N+2;de=a-e;C=4*(a+1)*(e+1)
    c=[1]
    for k in range(N):
        h,rem=divmod(de*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
        assert rem==0;c.append(h)
    v=lambda k:c[k] if 0<=k<=N else 0
    B=[v(k-1)+v(k+1) for k in range(N+2)]
    D=[v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    beta=[comb(sig,k+1) for k in range(N+2)]
    return N,sig,de,C,v,B,D,beta

def be(sig,C,j,i,beta):
    X=2*j-sig+2;Y=2*i-sig+2
    P=C-X*X;V=P+2*(sig-X)
    h=P*beta[j]-V*beta[i]
    z=Y*(P*beta[j]+V*beta[i])
    return h>=0 and C*h*h>=z*z

def crossing_case(sig,C,X,Y,s):
    if 9*C>=4*sig*sig:
        assert 4*C*s*s>25*(sig+1)**2
        return 0
    if C*X*X<=(sig+1)**2:
        assert 25*C*Y*Y>676*(sig+1)**2
        return 1
    if 4*C*X*X<=9*(sig+1)**2:
        assert 100*C*Y*Y>2809*(sig+1)**2
        return 2
    assert 100*C*s*s>441*(sig+1)**2
    return 3

# The stronger constant-gap assertion is false.
N,sig,de,C,v,B,D,beta=row(103,10)
j,i=58,62
W=v(j)*B[i]-B[j]*v(i)
assert W<0
gap2=F(C*(i-j)**2,(sig+1)**2)
assert gap2==F(4576,841) and gap2<F(25,4)
assert be(sig,C,j,i,beta) and D[j]-D[i]>=abs(W)
print("5/2 CONTROL",(103,10,58,62),"SQUARED GAP",gap2,
      "BE AND E PASS")

counts=[0,0,0,0];rows=anchors=pairs=0
for a in range(9,args.max_a+1):
    for e in range(7,a-1):
        N=a+e;sig=N+2;C=4*(a+1)*(e+1)
        if 4*(a-e)<=sig or C<4096 or C<30*(sig+1):continue
        rows+=1
        N,sig,de,C,v,B,D,beta=row(a,e)
        base=(N+1)//2
        if e%2 and N%2==0:
            scale=1 if v(base-1)>0 else -1
        elif e%2:
            scale=-1 if v(base)>0 else 1
        else:
            scale=1 if v(base)>0 else -1
        for j in range(N//2+1,N-3):
            X=2*j-N
            if C*X*X>=36*(sig+1)**2:break
            anchors+=1;seen=False
            if 9*C<4*sig*sig and 4*C*X*X<=9*(sig+1)**2:
                hj=scale*v(j)
                lj=scale*(sig*v(j-1)-de*v(j))
                assert lj>=0 and ((hj>=0) if e%2==0 else (hj<=0))
            for i in range(j+1,N+2):
                Y=2*i-N
                if 4*Y*Y>=C:break
                W=v(j)*B[i]-B[j]*v(i)
                seen |= W<0 or W==0 and v(j)*v(i)+B[j]*B[i]<0
                if i-j<4 or not seen:continue
                pairs+=1
                assert D[j]-D[i]>=abs(W),(a,e,j,i)
                assert be(sig,C,j,i,beta),(a,e,j,i,"BE")
                counts[crossing_case(sig,C,X,Y,i-j)]+=1
if args.max_a==220:
    assert (rows,anchors,pairs)==(12322,44062,1542704)
    assert counts==[1368370,22185,18873,133276]
print("CENTRAL ROWS",rows,"ANCHORS",anchors,
      "LONG PAIRS",pairs,"CASES",counts)

extra=0
for e in (7,8,15,30,60):
    for a in (400,1000,4000):
        N,sig,de,C,v,B,D,beta=row(a,e)
        if 4*de<=sig or C<4096 or C<30*(sig+1):continue
        for j in range(N//2+1,N-3):
            X=2*j-N
            if C*X*X>=36*(sig+1)**2:break
            seen=False
            for i in range(j+1,N+2):
                Y=2*i-N
                if 4*Y*Y>=C:break
                W=v(j)*B[i]-B[j]*v(i)
                seen |= W<0 or W==0 and v(j)*v(i)+B[j]*B[i]<0
                if i-j<4 or not seen:continue
                assert be(sig,C,j,i,beta)
                assert D[j]-D[i]>=abs(W)
                crossing_case(sig,C,X,Y,i-j)
                extra+=1
assert extra==13513
print("LARGE CENTRAL PAIRS",extra)

# J at the turning point cannot be required on a short-arc anchor.
N,sig,de,C,v,B,D,beta=row(160,9)
j,K=120,125
X,Y=2*j-N,2*K-N
assert (X,Y,C)==(71,81,6440)
assert (Y-2)**2<C<=Y*Y
assert 25*X*X>9*C and 2*X<=sig and X*X<4*(9-1)*(160+4)
energy=lambda k:sig*(v(k)**2+v(k-1)**2)-2*de*v(k)*v(k-1)
J=C*(sig+Y)**2*(D[j]-D[K])**2-4*Y*Y*energy(j)*energy(K)
assert J<0
assert all(D[k]>D[k+1] for k in range(j,N+1))
for i in range(j+1,N+2):
    W=v(j)*B[i]-B[j]*v(i)
    assert W>=0
    assert W>0 or v(j)*v(i)+B[j]*B[i]>=0
    assert D[j]-D[i]>=abs(W)
print("TURNING CONTROL",(160,9,j,K),"J NEGATIVE",
      "SHORT ENDPOINTS",N+1-j,"E PASS")
print("PASS")

