#!/usr/bin/env python3
"""Test sign pattern: v_lam(M) = <Res(M (x) V_lam), D_1^{2r}>/2 has sign (-1)^{lam_2} for Sym-products M.
Sp(4) irreps restricted to K via Weyl (sympy), products in doubled ring.  usage: probe_su2_fm3_signcone.py L N R LAMMAX"""
import sys, itertools, sympy as sp
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
L,NF,R,LM=map(int,sys.argv[1:5])
z,w=sp.symbols('z w')
def a(e1,e2): return (z**e1-z**-e1)*(w**e2-w**-e2)
def resV(l1,l2):
    num=sp.expand(a(l1+2,l2+1)-a(l2+1,l1+2))
    q=sp.expand(sp.cancel(num/(z+1/z-w-1/w)))
    P=sp.Poly(sp.expand(q*z**80*w**80),z,w)
    return {(i-81,j-81):int(c) for (i,j),c in P.terms() if i-80>=1 and j-80>=1}
def cg(x,p): return range(abs(x-p),x+p+1,2)
def cmul(F,G):
    out={}
    for (p,q),u in F.items():
        for (c,d),v in G.items():
            for e in cg(p,c):
                for f in cg(q,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
H=lambda k:{(i,k-i):1 for i in range(k+1)}
D1={(1,0):1,(0,1):-1}
lams=[(l1,l2) for l1 in range(LM+1) for l2 in range(l1+1)]
RV={lam:resV(*lam) for lam in lams}
viol=0; tot=0; zeros_wrong=0
for n in range(0,NF+1):
    for ks in itertools.combinations_with_replacement(range(1,L+1),n):
        F={(0,0):1}
        for k in ks: F=cmul(F,H(k))
        for r in range(1,R+1):
            Fd=F
            for _ in range(2*r): Fd=cmul(Fd,D1)
            for lam in lams:
                val=sum(Fd.get(key,0)*c for key,c in RV[lam].items())  # <F*D^{2r}, Res V_lam>
                tot+=1
                s=(-1)**lam[1]
                if s*val<0:
                    viol+=1
                    if viol<=10: print("VIOL ks=%s r=%d lam=%s val=%d"%(ks,r,lam,val))
print("tested=%d sign violations=%d"%(tot,viol))
