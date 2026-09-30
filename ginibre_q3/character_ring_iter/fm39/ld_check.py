# Independent check: for random real-rooted P and windows, Q(t) = sum_j C(n,j) c_(x-1+j) t^j (n = C+2) is real-rooted
# (Sturm count over exact rationals via sympy), and (i) holds exactly.
import random, sympy as sp
from fractions import Fraction as Fr
from math import comb
random.seed(21); t=sp.symbols('t'); n_rr=n_tot=0; viol=0
for _ in range(150):
    roots=[Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,5)) for _ in range(random.randint(3,10))]
    c=[Fr(1)]
    for rho in roots: c=[(c[k] if k<len(c) else 0)+rho*(c[k-1] if k>=1 else 0) for k in range(len(c)+1)]
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    for x in range(0,N+1):
        for Cw in range(1,N-x+1):
            n=Cw+2; coeffs=[comb(n,j)*C_(x-1+j) for j in range(n+1)]
            Q=sp.Poly(sum(sp.Rational(v.numerator,v.denominator)*t**j for j,v in enumerate(coeffs)),t)
            if Q.is_zero: continue
            n_tot+=1
            deg=Q.degree(); real=sp.polys.polytools.count_roots(Q) if deg>0 else 0
            # count real roots with multiplicity via real_roots
            rr=len(sp.real_roots(Q)) if deg>0 else 0
            n_rr+= rr==deg
            D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
            viol+= T*T>(Cw+1)**2*D(x)*D(x+Cw)
print('windows',n_tot,' Q real-rooted (with multiplicity):',n_rr,' (i) violations:',viol)
