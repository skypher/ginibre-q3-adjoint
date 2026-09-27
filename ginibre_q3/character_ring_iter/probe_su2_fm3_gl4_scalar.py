#!/usr/bin/env python3
"""Scalar test: ell_r(Res S_lambda(C^4)) = <S_lambda(C^2_x+C^2_y), D_1^{2r}>_{SU2xSU2} >= 0 ?
usage: probe_su2_fm3_gl4_scalar.py MAXDEG RMAX"""
import sys, sympy as sp
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
MAXDEG,RMAX=map(int,sys.argv[1:3])
z,w=sp.symbols('z w'); X=[z,1/z,w,1/w]
def schur4(lam):
    lam=list(lam)+[0]*(4-len(lam))
    num=sp.Matrix(4,4,lambda i,j: X[j]**(lam[i]+3-i)).det(); den=sp.Matrix(4,4,lambda i,j: X[j]**(3-i)).det()
    return sp.expand(sp.cancel(num/den))
def to_chars(expr):
    q=sp.expand(expr*(z-1/z)*(w-1/w)); P=sp.Poly(sp.expand(q*z**80*w**80),z,w)
    return {(i-81,j-81):int(c) for (i,j),c in P.terms() if i-80>=1 and j-80>=1}
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
def parts(n,k=4,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    if k==0: return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,k-1,a): yield (a,)+r
D2=cmul({(1,0):1,(0,1):-1},{(1,0):1,(0,1):-1})
bad=0
for d in range(0,MAXDEG+1,2):
    for lam in (parts(d) if d else [()]):
        F=to_chars(schur4(lam)) if lam else {(0,0):1}
        vals=[]; R=F
        for r in range(0,RMAX+1):
            if r: R=cmul(R,D2)
            vals.append(R.get((0,0),0))
        flag='' if all(v>=0 for v in vals) else '  <-- NEGATIVE'
        if flag: bad+=1
        print("lam=%-14s ell_r, r=0..%d: %s%s"%(lam,RMAX,vals,flag),flush=True)
print("negative GL(4) irreps:",bad)
