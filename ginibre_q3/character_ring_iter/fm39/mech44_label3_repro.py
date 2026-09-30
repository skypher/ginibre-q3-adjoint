"""FM-MECH44 (astra_max_ceres): the {1,2} sector plus one label 3 (phi_r(h_2 h_1^a hat S_2^b), phi_r(hat S_3 h_1^a hat S_2^b) >= 0) via radial IBP;
general-n elimination; raywise-positivity counterexamples."""
import argparse
from functools import lru_cache
from fractions import Fraction as Q
from math import comb
import sympy as sp

argparse.ArgumentParser(
    description="FM-MECH44: radial certificate and uniform IBP elimination"
).parse_args()

t,z=sp.symbols("t z", positive=True)
NN,kk=sp.symbols("N kappa")
kap=(1-t)**2/(1+t)**2
beta=8*(1+t*t)/(1+t)**2
alpha=4*(1+3*t*t)/(1+t)**2
xx=2*sp.sqrt(z)
yy=2*sp.sqrt(z)*(1-t)/(1+t)
assert sp.factor(xx**2+yy**2-2-(beta*z-2))==0
assert sp.factor(xx**2+yy**2-xx*yy-2-(alpha*z-2))==0
assert sp.factor(beta-4-4*(t-1)**2/(1+t)**2)==0
assert sp.factor(alpha-3-(3*t-1)**2/(1+t)**2)==0
jac=sp.diff(xx,z)*sp.diff(yy,t)-sp.diff(xx,t)*sp.diff(yy,z)
assert sp.factor(jac+4/(1+t)**2)==0
assert sp.factor((4-xx**2)*(4-yy**2)-16*(1-z)*(1-kap*z))==0
logder=NN/z-1/(2*(1-z))-kk/(2*(1-kk*z))
ratio=1-2*z+z*(1-z)*logder
claimed=NN+1-(NN+3)*z+(1-kk)*z/(2*(1-kk*z))
assert sp.factor(ratio-claimed)==0
assert sp.factor(sp.diff((1-kk*z)/(1-z),z)-(1-kk)/(1-z)**2)==0
print("Radial coordinates, normalization, and IBP: PASS")

def rising(a,j):
    out=Q(1)
    for k in range(j):
        out*=a+k
    return out

radial_checks=0
for N in range(1,13):
    alpha=Q(2*(N+3),N+1)
    for tail in (Q(3,2),Q(2)):
        for beta_value in (4,5,8):
            for b in range(17):
                value=sum(
                    Q(comb(b,j))*beta_value**j*(-2)**(b-j)*
                    (alpha*rising(Q(N+1),j+1)/
                     rising(Q(N+1)+tail,j+1)
                     -2*rising(Q(N+1),j)/
                     rising(Q(N+1)+tail,j))
                    for j in range(b+1))
                assert value>=0
                radial_checks+=1
print("Positive radial IBP: symbolic identity and",
      radial_checks,"exact endpoint cases")

def add(p,q):
    out=p.copy()
    for mon,c in q.items():
        out[mon]=out.get(mon,0)+c
    return {mon:c for mon,c in out.items() if c}

def scale(p,c):
    return {mon:c*v for mon,v in p.items() if c*v}

def mul(p,q):
    out={}
    for (i,j),v in p.items():
        for (k,l),w in q.items():
            mon=(i+k,j+l)
            out[mon]=out.get(mon,0)+v*w
    return {mon:c for mon,c in out.items() if c}

ONE={(0,0):1}
S={(1,0):1,(0,1):1}
D={(1,0):1,(0,1):-1}
Z={(2,0):1,(0,2):1,(0,0):-2}
P={(1,1):1}

@lru_cache(None)
def power(which,k):
    if k==0:
        return ONE
    return mul(power(which,k-1),{"s":S,"d":D,"z":Z,"p":P}[which])

@lru_cache(None)
def base(A,E,b):
    return mul(mul(power("s",A),power("d",E)),power("z",b))

@lru_cache(None)
def cm(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

def moment(poly,k=0):
    return sum(v*cm(i+k)*cm(j+k) for (i,j),v in poly.items())

@lru_cache(None)
def R(A,E,k,b):
    return moment(base(A,E,b),k)

@lru_cache(None)
def coeff(A,E,k,b):
    if k<0 or b<0:
        return ()
    if k==0:
        return (Q(0),)*b+(Q(1),)
    terms=[(4*(A-E),coeff(A,E,k-1,b)),
           (4*(k-1),coeff(A,E,k-2,b+1)),
           (8*(k-1),coeff(A,E,k-2,b))]
    if b:
        terms.append((12*b,coeff(A,E,k,b-1)))
    length=max(len(v) for _,v in terms)
    out=[Q(0)]*length
    for c,v in terms:
        for j,x in enumerate(v):
            out[j]+=c*x
    return tuple(x/Q(A+E+2*b+2*k+4) for x in out)

checks=0
for A in range(0,17,2):
    for E in range(0,17,2):
        for b in range(21):
            assert R(A,E,0,b+1)>=abs(R(A,E,1,b))
            checks+=1
print("Label-3 moment inequalities:",checks,"failures: 0")

checks=0
for A in range(0,9,2):
    for E in range(0,9,2):
        for k in range(1,6):
            for b in range(6):
                got=sum(c*R(A,E,0,j)
                        for j,c in enumerate(coeff(A,E,k,b)))
                assert got==R(A,E,k,b)
                checks+=1
print("IBP elimination identities:",checks)

s,p,zz=sp.symbols("s p Z")
Apol=[sp.Integer(2),s]
Bpol=[sp.Integer(1),s]
for k in range(2,13):
    Apol.append(sp.expand(s*Apol[-1]-p*Apol[-2]))
    Bpol.append(sp.expand(s*Bpol[-1]-p*Bpol[-2]))

def external(n,kind):
    eps=n%2 if kind=="plus" else (n-1)%2
    raw=0
    for j in range(n//2+1):
        m=n-2*j
        c=(-1)**j*comb(n-j,j)
        if kind=="plus":
            raw+=c*Apol[m]
        elif m:
            raw+=c*Bpol[m-1]
    out=0
    for (i,k),c in sp.Poly(raw,s,p).terms():
        assert i%2==eps
        out+=c*(zz+2+2*p)**((i-eps)//2)*p**k
    return eps,{
        (k,l):int(c)
        for (k,l),c in sp.Poly(sp.expand(out),p,zz).terms()
    }

@lru_cache(None)
def U(n,axis):
    return {
        ((n-2*j,0) if axis==0 else (0,n-2*j)):
        (-1)**j*comb(n-j,j)
        for j in range(n//2+1)
    }

def direct_external(n,kind):
    if kind=="plus":
        return add(U(n,0),U(n,1))
    out={}
    for j in range(n):
        out=add(out,mul(U(j,0),U(n-1-j,1)))
    return out

h2=direct_external(3,"minus")
assert h2==add(Z,P)
assert moment(mul(base(0,2,0),h2))==0
assert moment(mul(base(0,2,0),add(h2,ONE)))==2
print("Corrected h_2 = Z + P; phi_1(h_2) = 0; adding 1 gives 1")

bridges=0
for n in range(3,13):
    for kind in ("plus","minus"):
        eps,F=external(n,kind)
        expanded={}
        for (k,l),c in F.items():
            expanded=add(
                expanded,
                scale(mul(power("p",k),power("z",l)),c))
        expanded=mul(power("s",eps),expanded)
        direct=direct_external(n,kind)
        assert expanded==direct
        for r in range(1,4):
            for a in range(8):
                if (a+eps)%2:
                    continue
                for b in range(5):
                    got=sum(
                        Q(c)*v*R(a+eps,2*r,0,j)
                        for (k,l),c in F.items()
                        for j,v in enumerate(
                            coeff(a+eps,2*r,k,b+l)))
                    want=moment(mul(base(a,2*r,b),direct))
                    assert got==want
                    assert want>=0
                    bridges+=1
print("General-label direct bridges:",bridges,"failures: 0")

ray=4*Q(8,11)-3
full=Q(moment(
    mul(base(3,2,0),direct_external(4,"minus"))),2)
assert ray==Q(-1,11) and full==1
print("Label-4 radial control:",ray,
      "; full phi_1(h_3 h_1^3):",full)

ray6=(1024*rising(Q(3),4)/rising(Q(5),4)
      -1536*rising(Q(3),3)/rising(Q(5),3)
      +704*rising(Q(3),2)/rising(Q(5),2)
      -112*rising(Q(3),1)/rising(Q(5),1)+4)
full6=Q(moment(
    mul(base(0,4,1),direct_external(6,"plus"))),2)
assert ray6==Q(-36,35) and full6==1
print("Label-6 b=1 radial control:",ray6,
      "; full phi_2(S_6 S_2):",full6)
print("PASS")