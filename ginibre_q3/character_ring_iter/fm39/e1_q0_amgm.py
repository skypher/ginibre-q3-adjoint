# AM-GM (SONC-type) certificate for Q(p, 5+S) >= 0 on p, S >= 0: every negative monomial x^m is dominated by
# lam x^u + (1-lam) x^v >= x^m (weighted AM-GM, m = lam u + (1-lam) v, u, v positive monomials of Q), with the positive
# coefficients shared out by an LP; the solution is then verified in exact rational arithmetic.
import sympy as sp, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
src=open('e1_q0_Qpositivity.py').read().split("S=sp.symbols('S',nonnegative=True)")[0]
exec(src)
S=sp.symbols('S',nonnegative=True)
Qs=sp.Poly(sp.expand(poly.as_expr().subs(s,5+S)),p,S)
coef={m:int(c) for m,c in zip(Qs.monoms(),Qs.coeffs())}
negs=[m for m,c in coef.items() if c<0]; poss=[m for m,c in coef.items() if c>0]
# candidate pairs: u, v positive, m on the open segment between them with rational lam
cands=[]
for m in negs:
    for u,v in itertools.combinations(poss,2):
        du=(u[0]-m[0],u[1]-m[1]); dv=(v[0]-m[0],v[1]-m[1])
        if du[0]*dv[1]-du[1]*dv[0]!=0: continue          # not collinear
        # m = lam u + (1-lam) v with 0<lam<1  <=>  du and dv point in opposite directions
        if du[0]*dv[0]+du[1]*dv[1]>=0: continue
        # lam = |dv| / (|du|+|dv|) along the line
        nu=abs(du[0])+abs(du[1]); nv=abs(dv[0])+abs(dv[1])
        lam=Fr(nv,nu+nv)
        cands.append((m,u,v,lam))
print('negative monomials',len(negs),'candidate AM-GM pairs',len(cands))
# LP: variables w_c >= 0; for each negative m: sum_{c for m} w_c >= |coef_m|; for each positive u: sum lam_c w_c (u side) + (1-lam_c) w_c (v side) <= coef_u
nv=len(cands); A=[]; b=[]
for m in negs:
    row=[-1.0 if c[0]==m else 0.0 for c in cands]; A.append(row); b.append(-float(-coef[m]))
for u in poss:
    row=[float(c[3]) if c[1]==u else (float(1-c[3]) if c[2]==u else 0.0) for c in cands]; A.append(row); b.append(float(coef[u]))
scale=max(abs(x) for x in b)
res=linprog(np.zeros(nv),A_ub=np.array(A),b_ub=np.array(b)/scale,bounds=(0,None),method='highs')
print('LP status',res.status,res.message)
if res.status==0:
    w=[Fr(x*scale).limit_denominator(10**6) for x in res.x]
    # exact verification: rescale w up slightly to satisfy >= exactly, then check budgets
    from collections import defaultdict
    need=defaultdict(Fr); use=defaultdict(Fr)
    for wc,c in zip(w,cands):
        need[c[0]]+=wc; use[c[1]]+=wc*c[3]; use[c[2]]+=wc*(1-c[3])
    # fix up: scale each negative's allocations to cover exactly |coef|
    ok=True
    fac={m:(Fr(-coef[m])/need[m] if need[m]>0 else None) for m in negs}
    if any(v is None for v in fac.values()): ok=False
    use=defaultdict(Fr)
    for wc,c in zip(w,cands):
        f=fac[c[0]]; use[c[1]]+=wc*f*c[3]; use[c[2]]+=wc*f*(1-c[3])
    over=[(u,float(use[u]),coef[u]) for u in poss if use[u]>coef[u]]
    print('exact check: every negative covered:',ok,' positive budgets exceeded:',len(over),over[:5])
    if ok and not over: print('CERTIFICATE: Q(p,5+S) >= 0 for all real p, S >= 0 (weighted AM-GM), hence Q(p,s) > 0 for p >= 0, s >= 5.')
