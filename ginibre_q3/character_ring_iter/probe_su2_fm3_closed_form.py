#!/usr/bin/env python3
"""Check the Andreief closed form of I_r(F)=<F, D_1^{2r}> in tau-Schur coordinates.
usage: probe_su2_fm3_closed_form.py RMAX"""
import sys, sympy as sp
from fractions import Fraction
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
RMAX=int(sys.argv[1])
A,B,t,u=sp.symbols('A B t u')  # A=chi1(x), B=chi1(y); t=tau, u=tau'
def U(n,X):
    a,b=sp.Integer(1),X
    if n==0: return a
    for _ in range(n-1): a,b=b,sp.expand(X*b-a)
    return b
def S(p): return U(p,A)+U(p,B)
def H(q): return sp.cancel((U(q,A)-U(q,B))/(A-B))
def doubled(F):
    # coefficients F_ab of chi_a(x)chi_b(y) from polynomial in A,B
    P=sp.Poly(sp.expand(F),A,B); deg=max(P.degree(A),P.degree(B))
    out={}; rem=sp.expand(F)
    for a in range(deg,-1,-1):
        for b in range(deg,-1,-1):
            c=sp.Poly(rem,A,B).coeff_monomial(A**a*B**b)
            if c!=0:
                out[(a,b)]=int(c); rem=sp.expand(rem-c*U(a,A)*U(b,B))
    assert rem==0
    return out
def cg(a,p): return range(abs(a-p),a+p+1,2)
def ell(F,r):
    Fd=doubled(F); D={(1,0):1,(0,1):-1}; R=Fd
    for _ in range(2*r):
        out={}
        for (a,b),v in R.items():
            for (c,d),w in D.items():
                for e in cg(a,c):
                    for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+v*w
        R={k:v for k,v in out.items() if v}
    return R.get((0,0),0)
def tau_schur(F):
    # substitute A+B=4cc', AB=4(t+u-1); keep even part in cc' -> (A+B)^2=16 t u
    P,W=sp.symbols('P W')
    G=sp.expand(F.subs({A:P/2+W,B:P/2-W},simultaneous=True))
    PG=sp.Poly(G,P,W)
    ee=0
    for (i,j),c in PG.terms():
        if j%2: continue
        if i%2: continue
        # P=4cc' -> P^i = 4^i (tu)^(i/2); W^2 = (A-B)^2/4 = 4(1-t)(1-u)
        ee+=c*4**i*(t*u)**(i//2)*(4*(1-t)*(1-u))**(j//2)
    ee=sp.expand(ee)
    # Schur expansion in (t,u): s_{l1,l2} = (t^(l1+1)u^l2 - u^(l1+1)t^l2)/(t-u)
    num=sp.Poly(sp.expand(ee*(t-u)),t,u)
    cl={}
    for (i,j),c in num.terms():
        if i>j: cl[(i-1,j)]=sp.Rational(c)
    return cl
def pochh(x,n):
    r=sp.Integer(1)
    for k in range(n): r*=(x+k)
    return r
def R_closed(cl,r):
    h=sp.Rational(1,2); tot=0
    for (l1,l2),c in cl.items():
        tot+=c*(l1-l2+1)*pochh(h,l1+1)*pochh(h,l2)/(pochh(r+1,l1+2)*pochh(r+1,l2+1))
    return sp.nsimplify(tot)
tests={'H3':H(3),'S2':S(2),'S1^2*H3':sp.expand(S(1)**2*H(3)),'H3*H4*S1':sp.expand(H(3)*H(4)*S(1)),'S3*S1':sp.expand(S(3)*S(1)),'H5':H(5)}
one_cl={(0,0):1}
rr=sp.symbols('r')
for name,F in tests.items():
    cl=tau_schur(F)
    print(name,"tau-Schur coeffs:",dict(cl), " R_F(r)/R_1(r) =",sp.factor(R_closed(cl,rr)/R_closed(one_cl,rr)))
    for r in range(0,RMAX+1):
        direct=Fraction(ell(F,r),ell(1,r))
        cf=R_closed(cl,r)/R_closed(one_cl,r)
        assert sp.Rational(direct.numerator,direct.denominator)==cf, (name,r,direct,cf)
print("closed form agrees with direct doubled-ring computation for r<=%d"%RMAX)
