# Positivity certificate for Mars's polynomial Q(p,s) (e1_q0_closedform.py) on p >= 0, s >= 5 (integers; we prove it for reals):
# substitute s = 5 + S and check the coefficients of Q(p, 5+S) (all >= 0 would prove Q > 0 for p, S >= 0).
# If some are negative, multiply by (1 + p + S)^k (Polya) and retry.
import sympy as sp
src=open('e1_q0_closedform.py').read()
src=src.split("# Boundary values at s=3,4,5")[0]+src[src.index("# Symbolic normalized closed form"):].split("assert poly.degree(p) == 12")[0]
exec(src)
S=sp.symbols('S',nonnegative=True)
Qs=sp.Poly(sp.expand(poly.as_expr().subs(s,5+S)),p,S)
neg=[(m,c) for m,c in zip(Qs.monoms(),Qs.coeffs()) if c<0]
print('Q(p,5+S): terms',len(Qs.coeffs()),'negative coefficients',len(neg), neg[:8])
for k in range(1,8):
    if not neg: break
    Pk=sp.Poly(sp.expand(Qs.as_expr()*(1+p+S)**k),p,S)
    neg=[(m,c) for m,c in zip(Pk.monoms(),Pk.coeffs()) if c<0]
    print('Polya k=%d: negative coefficients %d'%(k,len(neg)), neg[:5])
# Numerical scan of Q on reals and integers
import numpy as np
f=sp.lambdify((p,S),Qs.as_expr(),'numpy')
Pg,Sg=np.meshgrid(np.linspace(0,60,601),np.linspace(0,60,601))
vals=f(Pg,Sg)
# normalize by a positive majorant to see relative sign
k=np.unravel_index(np.argmin(vals),vals.shape); print('min over real grid [0,60]^2 of Q(p,5+S):',vals[k],'at p=%.2f S=%.2f'%(Pg[k],Sg[k]))
# directions at infinity: Q(t cos a, t sin a) leading homogeneous part
lead=sp.Poly(sp.expand(Qs.as_expr()),p,S); d=lead.total_degree()
top=sum(c*p**m[0]*S**m[1] for m,c in zip(lead.monoms(),lead.coeffs()) if m[0]+m[1]==d)
print('total degree',d,'top homogeneous part:',sp.factor(top))
ints=[(pp,ss) for pp in range(0,401,1) for ss in range(0,396,5)]
mn=min((int(Qs.eval({p:pp,S:ss})),pp,ss+5) for pp in range(0,201,7) for ss in range(0,196,7))
print('min over integer sample p<=200, s<=200:',mn)
