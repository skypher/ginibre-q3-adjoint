# FM49 verifier: Q3 for [q^-, 1^-, p^+, 1^+ x a], p > q >= 1 (q >= p is FM39, p = 1 is FM11)
import sympy as sp, time
t0=time.time()
p,q,h,R=sp.symbols('p q h R', positive=True)
S,D=p+q,p-q; N=S+2*h
# (A) symbolic closed form: beta(t)=(2t/N)binom(N,N/2+t); ratios to binom(N,h); R = rho = binom(N,p+h)/binom(N,h)
Bn={-1:(S+h)/(h+1), 0:sp.Integer(1), 1:h/(S+h+1), 2:h*(h-1)/((S+h+1)*(S+h+2)), 3:h*(h-1)*(h-2)/((S+h+1)*(S+h+2)*(S+h+3))}
Pn={-1:(p+h)/(q+h+1), 0:sp.Integer(1), 1:(q+h)/(p+h+1), 2:(q+h)*(q+h-1)/((p+h+1)*(p+h+2))}
def bS(k): return (S+2*k)/N*Bn[k]          # beta(S/2+k)/binom(N,h)
def bD(k): return (D+2*k)/N*R*Pn[k]        # beta(D/2+k)/binom(N,h)
def g(b,k): return b(k)**2-b(k+1)*b(k-1)   # g(v) = beta(v)^2 - beta(v+1)beta(v-1)
def al(b,k): return b(k-1)+b(k)+b(k+1)
Delta_over = g(bD,0) - g(bS,1) - (al(bS,1)*bD(0) - bS(1)*al(bD,0))
Psi=((D*D+S+2*h)*(S+h+1)**2*(S+h+2)*R**2 - D*(p+1)*(q+1)*(S+2)*(S+h+1)*(S+2*h+1)*R
     - h*(p+h+1)*(q+h+1)*(S*S+5*S+2*h+4))
claimed = Psi/(N*(S+h+1)**2*(S+h+2)*(p+h+1)*(q+h+1))
print("(A) closed form (3) symbolic identity:", sp.simplify(sp.together(Delta_over-claimed))==0, f"{time.time()-t0:.0f}s", flush=True)
# (A') telescoping of m_a(A,0) and regrouped determinant form vs Brauer formula (1), exact, small grid
from math import comb
from fractions import Fraction as Fr
def Cb(n,k):
    if k.denominator!=1: return 0
    k=int(k); return comb(n,k) if 0<=k<=n else 0
def c(a,P,Q): return Cb(a,Fr(a+P+Q,2))*Cb(a,Fr(a-P+Q,2))
def m1(a,A,B): return (c(a,A,B)-c(a,A+4,B)-c(a,A,B+2)+c(a,A+4,B+2)-c(a,A+1,B-1)+c(a,A+3,B-1)+c(a,A+1,B+3)-c(a,A+3,B+3))
def Bb(a,t): return Cb(a,Fr(a,2)+t)
def be(a,t): return Bb(a,t-Fr(1,2))-Bb(a,t+Fr(1,2))
def alp(a,t): return Bb(a,t-Fr(3,2))-Bb(a,t+Fr(3,2))
bad=0
for a in range(0,25):
  for A in range(0,a+2):
    for B in range(0,A+1):
      u,v=Fr(A+B+3,2),Fr(A-B+1,2)
      if m1(a,A,B)!=alp(a,u)*be(a,v)-be(a,u)*alp(a,v): bad+=1
      for t in (u,v):
        if alp(a,t)!=be(a,t-1)+be(a,t)+be(a,t+1): bad+=1
        if be(a,t)*(a+1)!=2*t*Cb(a+1,Fr(a+1,2)+t): bad+=1
print("(A') regrouped form, alpha=beta+beta+beta, Pascal (2): mismatches", bad, flush=True)
# (B) q=1,2 factorisations and discriminant
P_,H_=sp.symbols('P_ H_')
Q1=12*H_**2+(-4*P_**2+4*P_+36)*H_+P_**4+2*P_**3-5*P_**2+2*P_+24
Q2=12*H_**2+(-4*P_**2+4*P_+60)*H_+P_**4+2*P_**3-9*P_**2+6*P_+72
PsiPH=Psi.subs({p:P_,h:H_},simultaneous=True)
e1=sp.simplify(PsiPH.subs({q:1,R:(P_+H_+1)/(H_+1)})-(P_+1)*(P_+H_+1)*(P_+2*H_+2)/(H_+1)**2*Q1)
e2=sp.simplify(PsiPH.subs({q:2,R:(P_+H_+1)*(P_+H_+2)/((H_+1)*(H_+2))})-(P_+1)*(P_+H_+1)*(P_+2*H_+2)*(P_+2*H_+3)*(P_+2*H_+4)/((H_+1)**2*(H_+2)**2)*Q2)
d1,d2=sp.factor(sp.discriminant(Q1,H_)),sp.factor(sp.discriminant(Q2,H_))
print("(B) q=1,2 factorisations:",e1==0 and e2==0,"; disc:",d1,d2, flush=True)
# 2p^4+8p^3+2p^2-12p-9 > 0 for p>=2: value 71 at 2, derivative 8p^3+24p^2+4p-12>0
# (C) q>=3: Jensen rho >= (1+y)^q >= B3, y = 2p/(2h+q+1); Psi(B3) >= 0 on the orthant
P,Q,H=sp.symbols('P Q H', nonnegative=True)
sub={q:3+Q, p:4+Q+P, h:H}
y=2*p/(2*h+q+1)
B3=1+q*y+q*(q-1)/2*y**2+q*(q-1)*(q-2)/6*y**3
num=sp.expand(sp.cancel(sp.together((Psi.subs(R,B3)*(2*h+q+1)**6*36).subs(sub,simultaneous=True))))
poly=sp.Poly(num,P,Q,H); d=dict(poly.terms())
negs={m:cf for m,cf in d.items() if cf<0}
print("(C) terms",len(d),"negative",len(negs), f"{time.time()-t0:.0f}s", flush=True)
used=set(); ok=True
for (i,j,k),cf in sorted(negs.items()):
    pair={(2,7):((0,8),(4,6)),(3,6):((1,7),(5,5)),(3,7):((1,8),(5,6))}[(i,k)]
    m1_,m2_=(pair[0][0],j,pair[0][1]),(pair[1][0],j,pair[1][1])
    c1,c2=d.get(m1_,0),d.get(m2_,0)
    good = c1>0 and c2>0 and cf*cf <= 4*c1*c2 and m1_ not in used and m2_ not in used
    used |= {m1_,m2_}; ok &= good
    print("   ",(i,j,k),cf,"<-",m1_,c1,m2_,c2,"AM-GM ok" if good else "FAIL")
print("(C) AM-GM certificate (disjoint pairs, midpoint exponents):",ok, flush=True)
# (D) Jensen: f(i)=log(1+p/(h+i)) convex in i: f''= 1/(h+i)^2 - 1/(h+i+p)^2 > 0 ; binomial expansion of (1+y)^q has nonneg terms.
# (E) Psi quadratic in R, leading coeff > 0, Psi(0) = -h(...) <= 0  =>  {R>=0: Psi>=0} = [x_+, inf); B3>0, Psi(B3)>=0 => rho>=B3>=x_+.
print("done", f"{time.time()-t0:.0f}s")
