"""FM-MECH66: exact proofs and bounded vetting. No file writes."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, prod
import sympy as S

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n",type=int,default=120,
                help="upper N for the diagnostic band census")
args=ap.parse_args()

# Rotation-contraction and the fixed change from psi.
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

# Rational certificates for the angle constants.
sin3=lambda x:x-x**3/6
t=F(1815,1024)
assert F(11,8)*(1-t/6+t*t/120-t**3/5040)>1
for ang,upper in ((F(9,25),F(103,300)),
                  (F(1,5),F(7,37)),
                  (F(16,25),F(89,150)),
                  (F(4,25),F(14,89)),
                  (F(77,100),F(52,75)),
                  (F(43,500),F(7,82))):
    assert sin3(ang)>upper
# Machin's identity supplies the rational lower bound for pi.
assert 16*(F(1,5)-F(1,375))-F(4,239)>F(157,50)
angles=(F(91,64)*F(347,200)+F(1,5),
        F(91,64)*F(403,200)+F(4,25),
        F(91,64)*F(429,200)+F(43,500))
assert all(z<F(157,50) for z in angles)
print("ANGLE BOUNDS",angles)

exp_lower=lambda t,n:sum(t**k/factorial(k) for k in range(n+1))
assert 23*F(103,300)**2<3
assert F(7,10)*3-F(307,150)==F(4,75)
assert exp_lower(F(623,750),7)>F(21359,9456)
assert exp_lower(F(679,375),12)>F(6731,1104)
for y in (F(1,4),F(29,32)):
    assert 30*y**3-30*y+2<0
assert exp_lower(F(1947,320),22)>432
assert F(75,64)/432<F(1,360)
assert F(375,16)/432<F(1,18)
margin=F(359,360)**2-4*F(33,32)**2/F(18)
assert margin==F(786023,1036800)
aa,ee=S.symbols("aa ee")
assert S.expand(4*(aa+1)*(ee+1)-36*(aa+ee-7)
                -4*(aa-8)*(ee-8))==0
assert 36*(64-9)>30*(64+1)
print("BINOMIAL AND TAIL CONSTANTS PASS; MARGIN",margin)

# The e=3 crossing bound.
s,z=S.symbols("s z")
A=6*(s-1)*(s-2)
poly=S.Poly(S.expand((25*A-144*(s+1)**2).subs(s,260+z)),z)
assert poly.all_coeffs()==[6,2382,213876]
assert exp_lower(F(12,5),6)>F(7*2339,2015)
assert exp_lower(F(14663,5184),6)>F(7*2339,2015)
assert 4*F(3,4)/(1-F(3,4)**2)**2>F(31,2)
assert exp_lower(F(43049111,8294400),18)>170
assert F(2339,2015)*F(309,11)/170<1
print("E3 CROSSING CONSTANTS PASS")

def data(a,e):
    N=a+e;si=N+2;dd=a-e;C=4*(a+1)*(e+1)
    c=[1]
    for k in range(N):
        val,rem=divmod(dd*c[k]-(N-k+1)*(c[k-1] if k else 0),k+1)
        assert not rem
        c.append(val)
    v=lambda k:c[k] if 0<=k<=N else 0
    B=[v(k-1)+v(k+1) for k in range(N+2)]
    D=[v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H=[si*(v(k)**2+v(k-1)**2)-2*dd*v(k)*v(k-1)
       for k in range(N+2)]
    beta=[si]
    for k in range(N+1):
        val,rem=divmod(beta[-1]*(si-k-1),k+2)
        assert not rem
        beta.append(val)
    return N,si,C,v,B,D,H,beta

def be(si,C,j,i,beta):
    N=si-2;X=2*j-N;Y=2*i-N
    P=C-X*X;V=P+2*(si-X)
    b=P*beta[j]-V*beta[i]
    A=Y*(P*beta[j]+V*beta[i])
    return b>=0 and C*b*b>=A*A

# Bounded e=3 census; the unbounded proof is analytic.
total=cross=0
for a in range(5,101):
    N,si,C,v,B,D,H,beta=data(a,3)
    M=3*si-2;A=6*(si-1)*(si-2)
    norm=4*prod(range(si-4,si+1))*si*(si-1)*(si-2)
    for j in range(N//2+1,N):
        X=2*j-N
        for i in range(j+1,N+2):
            Y=2*i-N
            G=A+(X*X-M)*(Y*Y-M)
            W=v(j)*B[i]-B[j]*v(i)
            rhs=(Y*Y-X*X)*beta[j]*beta[i]*X*Y*G
            assert W*norm==rhs
            assert D[j]-D[i]>=abs(W)
            total+=1
            if W<=0 and Y*Y<C and 100*Y*Y<81*C:
                assert X*X<M and (M-X*X)*(Y*Y-M)>=A
                assert be(si,C,j,i,beta)
                cross+=1
assert (total,cross)==(48104,2084)
print("E3 PAIRS",total,"CROSSING BE",cross)

# Samples beyond the earlier omega<64 finite remainder.
large=0
for a in (255,1000,4000):
    N,si,C,v,B,D,H,beta=data(a,3)
    for j in range(N//2+1,N):
        X=2*j-N
        if X*X>=3*si-2:break
        for i in range(j+4,N+2):
            Y=2*i-N
            if Y*Y>=C:break
            W=v(j)*B[i]-B[j]*v(i)
            assert D[j]-D[i]>=abs(W)
            if W<=0 and 100*Y*Y<81*C:
                assert be(si,C,j,i,beta)
            if 100*Y*Y>=81*C:
                drop=D[j]-D[i]
                assert C*(si+Y)**2*drop**2>=4*Y*Y*H[j]*H[i]
            large+=1
print("E3 LARGE-ROW SAMPLES",large)

# Fixed exception: e=7, omega>=64, kappa<15.
fixed=0
fleast=None
for a in range(127,134):
    N,si,C,v,B,D,H,beta=data(a,7)
    assert C>=4096 and C<30*(si+1)
    for j in range(N//2+1,N-3):
        for i in range(j+4,N+2):
            slack=D[j]-D[i]-abs(v(j)*B[i]-B[j]*v(i))
            assert slack>=0
            fixed+=1
            item=slack,(a,7,j,i)
            if fleast is None or item<fleast:fleast=item
assert fixed==15341
assert fleast==(13847976002664,(127,7,130,134))
print("FIXED E7 BOX",fixed,"MINIMUM",fleast)

# Uniform-region census: vetting, not a finite-box proof.
rows=pairs=long_inner=tail=0
least=None
for N in range(8,args.max_n+1):
    for e in range(4,(N-2)//2+1):
        a=N-e;si=N+2;C=4*(a+1)*(e+1)
        if 4*(a-e)<=si or C<4096 or C<30*(si+1):
            continue
        rows+=1
        N,si,C,v,B,D,H,beta=data(a,e)
        for j in range(N//2+1,N-3):
            X=2*j-N
            if C*X*X<36*(si+1)**2 or 25*X*X>9*C or 2*X>si:
                continue
            seen=False
            for i in range(j+1,N+2):
                Y=2*i-N
                W=v(j)*B[i]-B[j]*v(i)
                if W<0 or W==0 and v(j)*v(i)+B[j]*B[i]<0:
                    seen=True
                if i-j<4:continue
                drop=D[j]-D[i]
                slack=drop-abs(W)
                assert slack>=0,(a,e,j,i)
                pairs+=1
                item=slack,(a,e,j,i)
                if least is None or item<least:least=item
                if 64*Y*Y>=49*C:
                    assert C*(si+Y)**2*drop**2>=4*Y*Y*H[j]*H[i]
                    tail+=1
                elif seen:
                    assert 25*C*(Y-X)**2>196*(si+1)**2
                    assert be(si,C,j,i,beta),(a,e,j,i)
                    long_inner+=1
print("BAND",rows,pairs,"INNER LONG",long_inner,"TAIL",tail)
print("BAND MINIMUM",least)

# Independent definition-level Catalan checks.
@lru_cache(None)
def moment(power,label):
    if (power+label)%2:return 0
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
for e in range(3,9):
    for a in range(e+2,16):
        N,si,C,v,B,D,H,beta=data(a,e)
        for j in range(N//2+1,N-3):
            for i in range(j+4,N+2):
                p,q=i+j-N-1,i-j-1
                W=v(j)*B[i]-B[j]*v(i)
                assert kernel(e,a,p,q)==W
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2))==D[j]-D[i]
                bridges+=1
print("CATALAN BRIDGES",bridges)
print("PASS")
