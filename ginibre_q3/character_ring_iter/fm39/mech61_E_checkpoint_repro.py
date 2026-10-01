"""FM-MECH61: exact obstruction, checkpoint theorems, and fixed remainder."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, prod
import sympy as S

argparse.ArgumentParser(description=__doc__).parse_args()

# Endpoint ratio and infinite R3 obstruction.
ss, om, y = S.symbols("ss om y")
den = (ss+y+2)*(om-y-2)*(om+y)
num = (ss-y)*(om+y+2)*(om-y)
assert S.expand(den-num-2*((y+1)*(om**2+1-(y+1)**2)
                         -2*om*(ss+1))) == 0
h,z = S.symbols("h z")
aa = 80*h*h+60*h+11
sig = aa+6
xx,yy = 18*h+7,40*h+15
cc = 20*(aa+1)
K2 = lambda x: x*x-sig
K4 = lambda x: x**4-(6*sig-8)*x*x+3*sig*(sig-2)
G = (24*sig*(sig-1)*(sig-2)*(sig-3)
     +12*(sig-2)*(sig-3)*K2(xx)*K2(yy)+K4(xx)*K4(yy))
assert all(v > 0 for v in S.Poly(S.expand(-G/960),h).all_coeffs())
assert S.expand(cc-yy*yy) == 15
assert S.expand(25*xx*xx-4*cc) == 1700*h*h+1500*h+265
assert S.expand(cc-4*xx*xx) == 304*h*h+192*h+44
assert S.expand((18*(sig-yy)-(yy*yy-xx*xx)).subs(h,z+4)) == (
    164*z*z+724*z+132)
up = sum(F(9)**k/factorial(k) for k in range(51))
up += F(9)**51/factorial(51)/(1-F(9,52))
assert up < 8104
assert F(15,4*175**2) < F(1,8104)
print("ENDPOINT IDENTITY AND INFINITE BE FAILURE: PASS")

# Exact checkpoint constants.
def exp_lower(t,degree):
    return sum(t**k/factorial(k) for k in range(degree+1))

assert exp_lower(F(297,64),16) > 100
for t in (F(7,8),F(29,32)):
    assert 18*t**3-18*t+2 < 0
margin4 = F(79,80)**2-4*F(33,32)**2/F(5)
assert margin4 == F(199,1600)
gamma1 = F(31,4)*(F(9,10)**2-F(4,9)**2)
gamma2 = F(31,4)*(F(149,160)**2-F(4,9)**2)
assert exp_lower(gamma1,16) > 112
assert exp_lower(gamma2,18) > 170
v3,w3 = F(171,14560),F(1539,7280)
assert F(25029,715)/170 < w3
assert 4*F(9,10)/(1-F(9,10)**2)**2-F(31,2) > 0
margin3 = (1-v3)**2-4*F(33,32)**2*w3
assert margin3 == F(65606479,847974400)

# Upper boundaries of R1/R2 and the new checkpoints.
assert exp_lower(F(45,32),3) > F(18,5) > F(257,75)
assert exp_lower(F(1107,640),4) > F(845,168) > F(949,189)
assert F(2)/(1-F(17,32)**2) < F(15,2)
assert F(2)/(1-F(21,32)**2) < F(75,8)
assert F(31,27)*15/100 < 1
assert F(2339,2015)*19/112 < 1
assert F(2339,2015)*F(309,11)/170 < 1
print("CHECKPOINT MARGINS",margin4,margin3)

def data(a,e):
    N=a+e
    si,dd=N+2,a-e
    C=4*(a+1)*(e+1)
    c=[1]
    for k in range(N):
        q,r=divmod(dd*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
        assert r == 0
        c.append(q)
    v=lambda k: c[k] if 0 <= k <= N else 0
    B=[v(k-1)+v(k+1) for k in range(N+2)]
    D=[v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H=[si*(v(k)**2+v(k-1)**2)-2*dd*v(k)*v(k-1)
       for k in range(N+2)]
    return N,si,C,v,B,D,H

def be(si,C,j,i):
    N=si-2
    x,y=2*j-N,2*i-N
    P=C-x*x
    V=P+2*(si-x)
    bj,bi=comb(si,j+1),comb(si,i+1)
    A=y*(P*bj+V*bi)
    B=P*bj-V*bi
    return B >= 0 and B*B*C >= A*A

def joint_parts(si,C,Y,drop,Hj,Hi):
    return C*(si+Y)**2*drop*drop,4*Y*Y*Hj*Hi

def million_interval(t):
    k=t.numerator*10**6//t.denominator
    return k,k+1

# Actual-row witnesses, including failure of the local metric.
for hh in (4,8):
    a=80*hh*hh+60*hh+11
    N,si,C,v,B,D,H=data(a,4)
    j=40*hh*hh+39*hh+11
    i=40*hh*hh+50*hh+15
    X,Y=2*j-N,2*i-N
    W=v(j)*B[i]-B[j]*v(i)
    drop=D[j]-D[i]
    assert W < 0 and not be(si,C,j,i)
    assert drop >= abs(W)
    local=F((C-Y*Y)*drop*drop,4*Y*Y*D[j]*D[i])
    jnum,jden=joint_parts(si,C,Y,drop,H[j],H[i])
    assert jnum >= jden
    if hh == 8:
        assert local < 1
    print("R3 WITNESS",(a,4,j,i),"LOCAL",million_interval(local),
          "JOINT",million_interval(F(jnum,jden)),"units=1e-6")

# Entire fixed remainder omega < 64, plus the J checkpoint claims.
rows=pairs=inner_long=0
r12=[0,0]
tail={3:0,4:0}
least=None
for e in range(3,31):
    for a in range(e+2,256):
        si=a+e+2
        if 4*(a-e) <= si or 4*(a+1)*(e+1) >= 4096:
            continue
        rows += 1
        N,si,C,v,B,D,H=data(a,e)
        for j in range(N//2+1,N-3):
            X=2*j-N
            if 2*X > si or X*X >= C:
                continue
            crossed=False
            for i in range(j+1,N+2):
                Y=2*i-N
                W=v(j)*B[i]-B[j]*v(i)
                if W < 0 or (W == 0 and v(j)*v(i)+B[j]*B[i] < 0):
                    crossed=True
                if i-j < 4:
                    continue
                drop=D[j]-D[i]
                slack=drop-abs(W)
                assert slack >= 0,(a,e,j,i)
                pairs += 1
                item=slack,(a,e,j,i)
                if least is None or item < least:
                    least=item
                typ=0
                if e == 3 and X*X < 3*si-2 and 100*Y*Y >= 81*C:
                    typ=3
                if e >= 4 and 4*X*X <= C and 64*Y*Y >= 49*C:
                    typ=4
                if typ:
                    jnum,jden=joint_parts(si,C,Y,drop,H[j],H[i])
                    assert jnum >= jden,(a,e,j,i)
                    tail[typ] += 1
                if not crossed or Y*Y >= C:
                    continue
                inner_long += 1
                region=-1
                if 16*X*X <= C and 4*Y*Y < C:
                    region=0
                elif 16*X*X > C and 25*X*X <= 4*C and 64*Y*Y < 25*C:
                    region=1
                if region >= 0:
                    assert be(si,C,j,i),(a,e,j,i)
                    r12[region] += 1
assert (rows,pairs,least) == (1494,1105010,(67,(6,3,5,9)))
assert tail == {3:133792,4:317222}
assert inner_long == 235571 and r12 == [31903,11519]
print("FIXED REMAINDER",rows,pairs,"MINIMUM",least)
print("FINITE J REGIONS",tail,"R1/R2 BE",r12)

# Independent Catalan bridges and RF normalizations.
@lru_cache(None)
def moment(power,label):
    if (power+label)%2:
        return 0
    return sum((-1)**h*comb(label-h,h)*
               (comb(power+label-2*h,(power+label-2*h)//2)//
                ((power+label-2*h)//2+1))
               for h in range(label//2+1))

@lru_cache(None)
def kernel(e,a,p,q):
    return sum((-1)**h*comb(e,h)*comb(a,k)*
               moment(e+a-h-k,p)*moment(h+k,q)
               for h in range(e+1) for k in range(a+1))

bridges=0
for e in (3,4):
    for a in range(e+2,19):
        N,si,C,v,B,D,H=data(a,e)
        for j in range(N//2+1,N):
            X=2*j-N
            for i in range(j+1,N+2):
                Y=2*i-N
                def kp(x):
                    ans=[1,x]
                    for l in range(1,e):
                        ans.append(x*ans[-1]-l*(si-l+1)*ans[-2])
                    return ans
                px,py=kp(X),kp(Y)
                summ=sum(F(px[l]*py[l],
                           factorial(l)*prod(range(si-l+1,si+1)))
                         for l in range(e%2,e+1,2))
                eta=F(factorial(e),4*prod(range(si-e-1,si+1)))
                rhs=(eta*(Y*Y-X*X)*comb(si,j+1)*comb(si,i+1)*summ)
                W=v(j)*B[i]-B[j]*v(i)
                assert rhs == W
                p,q=i+j-N-1,i-j-1
                assert kernel(e,a,p,q) == W
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2)) == D[j]-D[i]
                bridges += 1
assert bridges == 966
print("RF AND CATALAN BRIDGES",bridges)
print("PASS")
