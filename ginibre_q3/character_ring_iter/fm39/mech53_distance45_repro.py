import argparse
from fractions import Fraction
from functools import lru_cache
from math import comb,factorial
import random
import sympy as s

argparse.ArgumentParser(
    description="FM-MECH53: distance 4 and 5 certificates"
).parse_args()
T,u,b,x,v,w=s.symbols("T u b x v w")

def C(n,k):
    return comb(n,k) if 0<=k<=n else 0
def choose(z,k):
    if k<0:return s.Integer(0)
    return s.Rational(1,factorial(k))*s.prod(z-h for h in range(k))
@lru_cache(None)
def cat(j):
    return C(2*j,j)//(j+1)
def mu(j,l):
    return sum(C(l,h)*(-2)**(l-h)*cat(j+h) for h in range(l+1))

c=[s.Integer(1),x]
for j in range(1,12):
    c.append(s.expand((x*c[j]-(T-j+1)*c[j-1])/(j+1)))
ce=[]
for j in range(7):
    ce.append(s.expand(sum(
        co*u**(mon[0]//2)
        for mon,co in s.Poly(c[2*j],x).terms()
    )))

G={}
for d in range(7):
    D=T+2*b
    raw=sum(
        ce[j]*choose(b,l)*mu(j,l)*
        (choose(D-2*j-2*l,d-j-l)
         -choose(D-2*j-2*l,d-j-l-1))
        for j in range(d+1) for l in range(d-j+1)
    )
    G[d]=s.Poly(
        s.expand(factorial(d)*factorial(d+1)*raw),u,T,b
    )
    assert not G[d].as_expr().atoms(s.Float)
    assert all(co.q==1 for co in G[d].coeffs())
    assert s.Poly(G[d].as_expr(),u).LC()==1
    if d:
        A=d*(-(d-2)*T+2*(2*d-1)*b
             +s.Rational((2*d-1)*(d-7),3))
        assert s.expand(
            s.Poly(G[d].as_expr(),u).nth(d-1)-A
        )==0
print("Distance polynomials and general penultimate coefficient: PASS")

def positive(expr,variables):
    pp=s.Poly(s.expand(expr),*variables)
    assert all(co>=0 for co in pp.coeffs())
    return len(pp.terms()),min(pp.coeffs())

def bernstein(expr,N):
    pp=s.Poly(s.expand(expr),v)
    assert pp.degree()<=N
    return [
        s.expand(sum(
            pp.nth(j)*s.Rational(C(i,j),C(N,j))
            for j in range(min(i,pp.degree())+1)
        ))
        for i in range(N+1)
    ]

R={}
for d,B0,N in ((4,4,4),(5,6,8)):
    pp=s.Poly(G[d].as_expr(),u)
    A=pp.nth(d-1)
    q=u*u+A*u/2+(10*T*T-80*T-960*b if d==5 else 0)
    R[d]=s.Poly(
        s.expand(pp.as_expr()-(u if d==5 else 1)*q*q),u
    )
    outer=7*b-3 if d==4 else 6*b-1
    for k in range(d-1):
        cert=R[d].nth(k).subs(T,outer+v).subs(b,3+w)
        count,minimum=positive(cert,(v,w))
        print("outer",d,k,"terms",count,"minimum",minimum)
    assert R[d].nth(0).subs(
        {T:int(outer.subs(b,3)),b:3}
    )>0
    hi=7*b-4 if d==4 else 6*b-2
    all_coefficients=[]
    polynomial_count=0
    for k in range(d):
        aa=pp.nth(k).subs(T,4+(hi-4)*v)
        for bb in bernstein(aa,N):
            cc=s.Poly(s.expand(bb.subs(b,B0+w)),w)
            assert all(co>=0 for co in cc.coeffs())
            if k==0:assert cc.nth(0)>0
            all_coefficients+=cc.coeffs()
            polynomial_count+=1
    print("inner",d,"Bernstein degree",N,"b >=",B0,
          "polynomials",polynomial_count,
          "coefficient minimum",min(all_coefficients))
assert s.expand(
    R[5].nth(3)-s.Rational(55,4)*(T-6*b+2)**2-128
)==0

def row(e,a):
    c=[1]
    for sign,n in ((1,a),(-1,e)):
        for _ in range(n):
            c=[
                (c[j] if j<len(c) else 0)
                +sign*(c[j-1] if j else 0)
                for j in range(len(c)+1)
            ]
    return c
def ballot(m,p):
    if m<p or (m-p)%2:return 0
    h=(m-p)//2
    return C(m,h)-C(m,h-1)
def direct(e,a,b,p):
    c=row(e,a)
    t=e+a
    return sum(
        c[j]*C(b,l)*(-2)**(b-l)*C(l,h)*
        cat(j//2+h)*ballot(t-j+2*(l-h),p)
        for j in range(0,t+1,2)
        for l in range(b+1) for h in range(l+1)
    )
def exact_G(d,t,uu,bb):
    return sum(
        int(co)*uu**i*t**j*bb**k
        for (i,j,k),co in G[d].terms()
    )

for d,B0 in ((4,4),(5,6)):
    count=0
    least=None
    for bb in range(3,B0):
        lo=max(4,2*d+7-2*bb)
        hi=7*bb-4 if d==4 else 6*bb-2
        for t in range(lo,hi+1):
            p=t+2*bb-2*d
            for e in range(2,t-1):
                a=t-e
                num=exact_G(d,t,(a-e)**2,bb)
                den=factorial(d)*factorial(d+1)
                assert num%den==0
                val=num//den
                assert val==direct(e,a,bb,p)>0
                count+=1
                if least is None or val<least[0]:
                    least=(val,e,a,bb,p)
    print("finite remainder",d,"count",count,"least",least)

rng=random.Random(53)
bridges=0
for _ in range(128):
    d=rng.choice((4,5))
    e=rng.randrange(2,25)
    a=rng.randrange(2,25)
    bb=rng.randrange(3,14)
    p=e+a+2*bb-2*d
    if p<7:continue
    val=direct(e,a,bb,p)
    assert exact_G(d,e+a,(a-e)**2,bb)==(
        factorial(d)*factorial(d+1)*val
    )
    assert val>0
    bridges+=1
print("Additional direct Catalan bridges:",bridges)

# Gaussian leading part, checked against both closed descriptions.
H=[s.Integer(1),x]
for j in range(1,6):
    H.append(s.expand(x*H[j]-j*T*H[j-1]))
for d in range(7):
    top=sum(
        co*u**i*T**j*b**k
        for (i,j,k),co in G[d].terms() if i+j+k==d
    )
    raw=sum(
        T**k*(b-T/2)**l*u**j*cat(j+l)/
        (s.Integer(factorial(k))*factorial(l)*factorial(2*j))
        for j in range(d+1) for l in range(d-j+1)
        for k in (d-j-l,)
    )
    assert s.expand(
        top-factorial(d)*factorial(d+1)*raw
    )==0
    sos=factorial(d)*sum(
        T**(d-j)*H[j]**2/s.Integer(factorial(j))
        for j in range(d+1)
    )
    sos=sum(
        co*u**(mon[0]//2)
        for mon,co in s.Poly(s.expand(sos),x).terms()
    )
    assert s.expand(top.subs(b,0)-sos)==0
print("Gaussian generating function and Hermite-square identities: PASS")

# The distance-six continuation with lower-degree corrections fails.
pp=s.Poly(G[6].as_expr(),u)
A=pp.nth(5)
q0=u**3+A*u*u/2+35*T*T*u
bad=s.Poly(s.expand(pp.as_expr()-q0*q0),u).nth(4)
assert s.expand(
    bad.subs(T,s.Rational(11,2)*b)
    +(5005*b*b+26100*b-1628)/4
)==0
print("Degree-six template obstruction:",
      s.factor(bad.subs(T,s.Rational(11,2)*b)))

# The natural mixed-term repair on its uniform outer region.
q,rem=s.div(pp.nth(3),A,T)
rr=s.Poly(
    s.expand(pp.as_expr()-(u**3+A*u*u/2+q*u)**2),u
)
for k in range(5):
    cert=rr.nth(k).subs(
        T,s.Rational(11,2)*b-s.Rational(1,2)+v
    ).subs(b,6+w)
    positive(cert,(v,w))
print("Degree-six mixed-term repair: outer b>=6 PASS")

# Exact checks on the family excluding a fixed Gaussian margin.
for p in (8,10,12):
    bb=p*p
    t=p+2
    d=bb+1
    assert d<=2*p*p-p//2+1
    actual=direct(2,p,bb,p)
    lower=Fraction((bb-t//2)**d*cat(d),factorial(d))
    ratio=Fraction(abs(actual),1)/lower
    exponent=0
    while ratio < Fraction(1,10**(exponent+1)):
        exponent+=1
    assert exponent>=1
    print("Gaussian comparison family p",p,
          "actual positive",actual>0,
          "|F|/leading <= lower-test ratio < 10^-"+str(exponent))
print("PASS")
