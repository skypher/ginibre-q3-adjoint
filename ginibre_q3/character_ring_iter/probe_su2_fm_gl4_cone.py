#!/usr/bin/env python3
"""Admissibility (column 0 of F*D_1^m >= 0, m<=M) of GL(4) Schur functors S_lam(C^2_x + C^2_y)
restricted to SU(2)xSU(2), alone and multiplied by S_p. usage: probe_su2_fm_gl4_cone.py MAXDEG M"""
import sys, itertools, sympy as sp
if any(a in ('-h','--help') for a in sys.argv[1:]): print(__doc__); sys.exit(0)
MAXDEG, M = map(int, sys.argv[1:3])
z,w=sp.symbols('z w')
X=[z,1/z,w,1/w]
def schur4(lam):
    lam=list(lam)+[0]*(4-len(lam))
    num=sp.Matrix(4,4,lambda i,j: X[j]**(lam[i]+3-i)).det()
    den=sp.Matrix(4,4,lambda i,j: X[j]**(3-i)).det()
    return sp.expand(sp.cancel(num/den))
def to_chars(expr):
    q=sp.expand(expr*(z-1/z)*(w-1/w))
    P=sp.Poly(sp.expand(q*z**80*w**80),z,w)
    return {(i-81,j-81):int(c) for (i,j),c in P.terms() if i-80>=1 and j-80>=1}
def cg(a,p): return range(abs(a-p), a+p+1, 2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
D1={(1,0):1,(0,1):-1}
def first_fail(F):
    R=F
    for m in range(M+1):
        if m: R=cmul(R,D1)
        neg={a:v for (a,b),v in R.items() if b==0 and v<0}
        if neg: return (m,neg)
    return None
def parts(n,k=4,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    if k==0: return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,k-1,a): yield (a,)+r
for d in range(1,MAXDEG+1):
    for lam in parts(d):
        F=to_chars(schur4(lam))
        r=first_fail(F)
        extra=[]
        for p in (1,2,3):
            rr=first_fail(cmul(F,{(p,0):1,(0,p):1}))
            extra.append('S%d:%s'%(p,'ok' if rr is None else 'm=%d'%rr[0]))
        print("lam=%s  alone:%s  %s"%(lam,'ok' if r is None else 'FAIL m=%d %s'%r, ' '.join(extra)),flush=True)
