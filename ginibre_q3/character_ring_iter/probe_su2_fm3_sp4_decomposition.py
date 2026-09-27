#!/usr/bin/env python3
"""Decompose swap-symmetric SU(2)xSU(2) characters into Sp(4) characters (virtual).
usage: probe_su2_fm3_sp4_decomposition.py PMAX"""
import sys, sympy as sp
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
PMAX=int(sys.argv[1])
z,w=sp.symbols('z w')
def a(e1,e2): return (z**e1-z**-e1)*(w**e2-w**-e2)
def sp4num(l1,l2):
    m1,m2=l1+2,l2+1
    return sp.expand(a(m1,m2)-a(m2,m1))
def gg_num(F):
    # F: dict (a,b)->coef of chi_a(z)chi_b(w); return F*(z-1/z)(w-1/w)*(z+1/z-w-1/w) i.e. times full Sp4 denominator
    tot=0
    for (A,B),c in F.items():
        tot+=c*(z**(A+1)-z**-(A+1))*(w**(B+1)-w**-(B+1))
    return sp.expand(tot*(z+1/z-w-1/w))
def decompose(F):
    num=gg_num(F)
    P=sp.Poly(sp.expand(num*z**200*w**200),z,w)
    terms={(i-200,j-200):int(c) for (i,j),c in P.terms()}
    out={}
    # dominant chamber for numerator: exponent (m1,m2) with m1>m2>0 -> lambda=(m1-2,m2-1)
    for (i,j),c in terms.items():
        if i>j>0: out[(i-2,j-1)]=c
    # verify
    chk=sp.expand(sum(c*sp4num(l1,l2) for (l1,l2),c in out.items())-num)
    assert chk==0, "decomposition check failed"
    return out
for p in range(1,PMAX+1):
    print("S_%d ="%p, decompose({(p,0):1,(0,p):1}))
print("V1xV1 =",decompose({(1,1):1}))
print("D_1^2 =",decompose({(2,0):1,(0,2):1,(0,0):2,(1,1):-2}))
