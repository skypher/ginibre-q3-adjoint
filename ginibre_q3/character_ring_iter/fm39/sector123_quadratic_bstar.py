# Main agent: largest B with E[s^(2m) d^(2k) prod (x^2 +- 2B xy + y^2 - 2)] >= 0 at all sign vertices (argument: max number of factors).
import sympy as sp, sys
from fractions import Fraction as Q
from math import factorial
Bs,u,v=sp.symbols('B u v')
def mu(m,k): return sp.Rational(2*factorial(2*m)*factorial(2*m+1)*factorial(2*k)*factorial(2*k+1), factorial(m)**2*factorial(k)**2*factorial(m+k+1)*factorial(m+k+2))
Lp=(1+Bs)/2*u+(1-Bs)/2*v-2; Lm=(1-Bs)/2*u+(1+Bs)/2*v-2
N=int(sys.argv[1]); out=[]
for n in range(2,N+1):
  for m in range(0,5):
    for k in range(0,5):
      if k<m: continue
      worst=None
      for npos in range(n+1):
        P=sp.Poly(sp.expand(u**m*v**k*Lp**npos*Lm**(n-npos)),u,v)
        f=sp.expand(sum(c*mu(i,j) for (i,j),c in P.terms()))
        f=sp.Poly(f,Bs)
        cands=[r for r in sp.real_roots(f) if r>0]
        crit=None
        for r in sorted(cands,key=float):
            if f.eval(r+sp.Rational(1,10**6))<0: crit=float(r); break
        if crit is not None and (worst is None or crit<worst[0]): worst=(crit,npos)
      if worst: out.append((worst[0],n,2*m,2*k,worst[1])); 
  mn=min(o for o in out if o[1]==n); print("n",n,"min B*",round(mn[0],4),"at A,E,npos",mn[2:],flush=True)
