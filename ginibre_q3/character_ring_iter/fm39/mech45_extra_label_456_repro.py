import argparse
from functools import lru_cache
from math import comb
from fractions import Fraction as Q
import sympy as sp

argparse.ArgumentParser(
    description="FM-MECH45: exact cutoff and finite-box certificates"
).parse_args()

q,v,f,s,p,Z=sp.symbols("q v f s p Z")
c=q/(1+q*q)
u=c*c
B=2+4*q*q
beta=B+2
certificates={
 "h3":(6,f-1,sp.Integer(1)),
 "S4":(12,(1-2*u)*f*f+(1-8*u)*f-8*u,
       2*(1-2*u)*B+(1-4*u)**2),
 "h4/S5":(24,(1+c-u)*f*f-4*u*f-(2*c+1)**2,
          2*(1+c-u)*B),
 "h5":(25,(1-u)*f*f-(1+4*u)*f-4*u,
        2*(1-u)*B),
 "S6":(51,(1-3*u)*f**3+(1-8*u)*f*f+(-2+4*u)*f-2+16*u,
        3*(1-3*u)*B*B+2*(1-4*u)**2*B)
}

for name,(K,F,D) in certificates.items():
    num,den=sp.fraction(
        sp.factor((K+3)*F.subs(f,B)-2*beta*D))
    a=sp.Poly(sp.expand(num.subs(q,2*v-1)),v)
    bern=[
        sum(a.nth(j)*sp.Rational(comb(i,j),comb(64,j))
            for j in range(min(i,a.degree())+1))
        for i in range(65)
    ]
    assert min(bern)>=0
    assert sp.expand(
        den-(1+q*q)**(sp.degree(den,q)//2))==0
    print("cutoff",name,K,"Bernstein minimum",min(bern))

C=[sp.Integer(2),s]
T=[sp.Integer(1),s]
for j in range(2,8):
    C.append(sp.expand(s*C[-1]-p*C[-2]))
    T.append(sp.expand(s*T[-1]-p*T[-2]))

@lru_cache(None)
def external(n,kind):
    eps=n%2 if kind=="plus" else (n-1)%2
    raw=sum(
        (-1)**j*comb(n-j,j)*
        (C[n-2*j] if kind=="plus" else
         (T[n-2*j-1] if n-2*j else 0))
        for j in range(n//2+1))
    out=0
    for (i,k),a in sp.Poly(raw,s,p).terms():
        assert i%2==eps
        out+=a*(Z+2+2*p)**((i-eps)//2)*p**k
    out=sp.Poly(sp.expand(out),p,Z)
    return eps,{(k,l):int(a) for (k,l),a in out.terms()}

families=[
    (4,"minus",6,"h3"),(4,"plus",12,"S4"),
    (5,"minus",24,"h4/S5"),(5,"plus",24,"h4/S5"),
    (6,"minus",25,"h5"),(6,"plus",51,"S6")
]

for n,kind,K,name in families:
    eps,F=external(n,kind)
    ray=sp.expand(sum(
        a*(c*(f+2))**k*f**l for (k,l),a in F.items()))
    target=certificates[name][1]
    if n==5 and kind=="plus":
        target=target.subs(q,-q)
    assert sp.factor(ray-target)==0
print("All six ray-polynomial identities: PASS")

@lru_cache(None)
def cat(j):
    return comb(2*j,j)//(j+1)

@lru_cache(None)
def M0(m,r):
    A,E=2*m,2*r
    N=A+E
    row=[1]
    if N:
        row.append(A-E)
    for k in range(1,N):
        num=(A-E)*row[k]-(N-k+1)*row[k-1]
        assert num%(k+1)==0
        row.append(num//(k+1))
    return sum(
        row[k]*cat(k//2)*cat((N-k)//2)
        for k in range(0,N+1,2))

@lru_cache(None)
def moment(m,r,b):
    if m>r:
        return moment(r,m,b)
    if b==0:
        return M0(m,r)
    top=moment(m+1,r,b-1)+moment(m,r+1,b-1)
    assert top%2==0
    return top//2-2*moment(m,r,b-1)

def mixed(m,r,k,b):
    num=sum(
        (-1)**(k-j)*comb(k,j)*moment(m+j,r+k-j,b)
        for j in range(k+1))
    assert num%(4**k)==0
    return num//(4**k)

total=0
for n,kind,K,name in families:
    eps,F=external(n,kind)
    count=zeros=0
    minimum=None
    for N in range(2,K):
        for r in range(1,N+1):
            m=N-r
            if m<eps:
                continue
            for b in range(K-N):
                val=sum(
                    a*mixed(m,r,k,b+l)
                    for (k,l),a in F.items())
                assert val%2==0
                val//=2
                assert val>=0,(n,kind,r,2*m-eps,b,val)
                count+=1
                zeros+=val==0
                if val and (minimum is None or val<minimum):
                    minimum=val
    total+=count
    print("finite box",n,kind,"count",count,
          "zeros",zeros,"least positive",minimum)
print("TOTAL finite consumer cases:",total)

# A separate evaluator expands the original x,y integrand.
x,y=sp.symbols("x y")

def U(n,z):
    return sum(
        (-1)**j*comb(n-j,j)*z**(n-2*j)
        for j in range(n//2+1))

def direct(n,kind,r,a,b):
    W=U(n,x)+U(n,y) if kind=="plus" else sum(
        U(j,x)*U(n-1-j,y) for j in range(n))
    poly=sp.Poly(sp.expand(
        (x-y)**(2*r)*(x+y)**a*(x*x+y*y-2)**b*W),x,y)
    num=sum(
        int(c)*cat(i//2)*cat(j//2)
        for (i,j),c in poly.terms() if i%2==j%2==0)
    return Q(num,2)

bridges=0
for n,kind,K,name in families:
    eps,F=external(n,kind)
    for r in range(1,3):
        for a in range(4):
            if (a+eps)%2:
                continue
            for b in range(3):
                actual=direct(n,kind,r,a,b)
                m=(a+eps)//2
                reduced=Q(sum(
                    c*mixed(m,r,k,b+l)
                    for (k,l),c in F.items()),2)
                assert actual==reduced
                bridges+=1
print("Independent direct Catalan bridges:",bridges)

def rising(a,k):
    out=Q(1)
    for j in range(k):
        out*=a+j
    return out

z=sp.symbols("z")

def beta_average(poly,N,tail):
    return sum(
        Q(int(c))*rising(Q(N+1),k)/rising(Q(N+1)+tail,k)
        for (k,),c in sp.Poly(sp.expand(poly),z).terms())

bounds=0
for N in range(2,9):
    for beta_value,tail in ((4,Q(3,2)),(8,Q(2))):
        J=[
            beta_average((beta_value*z-2)**b,N,tail)
            for b in range(16)
        ]
        BB=beta_value-2
        for b in range(13):
            assert J[b]>=0
            for j in range(1,4):
                err=BB**j*J[b]-J[b+j]
                assert 0<=err<=Q(
                    2*beta_value*j*BB**(j-1),N+b+3)*J[b]
                bounds+=1
print("Exact radial moment bounds:",bounds)

for kind,qq,aa in (("minus",-1,2),("plus",1,1)):
    eps,F=external(7,kind)
    print("Label-7 polynomial",kind,
          sp.expand(sum(a*p**k*Z**l for (k,l),a in F.items())))
    cc=Q(qq,1+qq*qq)
    ray=sp.expand(sum(
        a*(cc*8*z)**k*(8*z-2)**l for (k,l),a in F.items()))
    radial=beta_average(ray*(8*z-2)**3,2,Q(2))
    true=direct(7,kind,1,aa,3)
    assert radial==Q(-32,105) and true==4
    print("Label-7 control",kind,"ray",radial,"consumer",true)
print("PASS")
