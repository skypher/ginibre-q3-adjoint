# AM-GM with pairs and triangles: x^m <= sum_t lam_t x^(u_t) when m = sum lam_t u_t (lam >= 0, sum 1), u_t positive monomials.
# LP allocation of positive coefficients; exact rational verification.
import sympy as sp, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from collections import defaultdict
src=open('e1_q0_Qpositivity.py').read().split("S=sp.symbols('S',nonnegative=True)")[0]
exec(src)
S=sp.symbols('S',nonnegative=True)
Qs=sp.Poly(sp.expand(poly.as_expr().subs(s,5+S)),p,S)
coef={m:int(c) for m,c in zip(Qs.monoms(),Qs.coeffs())}
negs=[m for m,c in coef.items() if c<0]; poss=[m for m,c in coef.items() if c>0]
def bary(m,u,v,w):
    # solve m = a u + b v + c w, a+b+c = 1
    M=sp.Matrix([[u[0],v[0],w[0]],[u[1],v[1],w[1]],[1,1,1]])
    if M.det()==0: return None
    sol=M.LUsolve(sp.Matrix([m[0],m[1],1]))
    lam=[Fr(int(sp.fraction(x)[0]),int(sp.fraction(x)[1])) for x in sol]
    return lam if all(l>=0 for l in lam) else None
cands=[]
for m in negs:
    near=[u for u in poss if abs(u[0]-m[0])<=6 and abs(u[1]-m[1])<=8]
    for u,v in itertools.combinations(near,2):
        du=(u[0]-m[0],u[1]-m[1]); dv=(v[0]-m[0],v[1]-m[1])
        if du[0]*dv[1]-du[1]*dv[0]==0 and du[0]*dv[0]+du[1]*dv[1]<0:
            nu=abs(du[0])+abs(du[1]); nv=abs(dv[0])+abs(dv[1]); lam=Fr(nv,nu+nv)
            cands.append((m,((u,lam),(v,1-lam))))
    for u,v,w in itertools.combinations(near,3):
        lam=bary(m,u,v,w)
        if lam and all(l>0 for l in lam): cands.append((m,((u,lam[0]),(v,lam[1]),(w,lam[2]))))
print('negative monomials',len(negs),'candidates',len(cands),flush=True)
A=[];b=[]
for m in negs:
    A.append([-1.0 if c[0]==m else 0.0 for c in cands]); b.append(float(coef[m]))
for u in poss:
    A.append([float(sum(l for (uu,l) in c[1] if uu==u)) for c in cands]); b.append(float(coef[u]))
scale=max(abs(x) for x in b)
res=linprog(np.zeros(len(cands)),A_ub=np.array(A),b_ub=np.array(b)/scale,bounds=(0,None),method='highs')
print('LP status',res.status,res.message,flush=True)
if res.status==0:
    w=[Fr(x*scale).limit_denominator(10**9) for x in res.x]
    need=defaultdict(Fr)
    for wc,c in zip(w,cands): need[c[0]]+=wc
    fac={m:Fr(-coef[m])/need[m] for m in negs}
    use=defaultdict(Fr)
    for wc,c in zip(w,cands):
        for (u,l) in c[1]: use[u]+=wc*fac[c[0]]*l
    over=[(u,float(use[u]),coef[u]) for u in poss if use[u]>coef[u]]
    print('exact check: positive budgets exceeded:',len(over),over[:5])
    if not over: print('CERTIFICATE: Q(p,5+S) >= 0 for all real p, S >= 0; hence the family phi_2(h_(p+2) h_2^2 h_1^(p+2s+2)) > 0 for s >= 5.')
