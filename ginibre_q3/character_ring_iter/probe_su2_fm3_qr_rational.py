#!/usr/bin/env python3
"""Rational form of Q_r(t_1..t_J) = int_K prod_i det(1-t_i g)^{-1} D_1(g)^{2r} dg via constant terms.
usage: probe_su2_fm3_qr_rational.py J R"""
import sys, sympy as sp
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
J,R=map(int,sys.argv[1:3])
z,w=sp.symbols('z w'); t=sp.symbols('t1:%d'%(J+1))
# SU(2) Haar: int f = CT_z[(1-z^2) f(z)] for class functions; do residue computation symbolically
def su2_int(expr,var):
    # expr rational in var; integral = CT[(1-var^2)*expr] over |var|=1 with |t|<1: sum of residues inside unit circle of (1-var^2)*expr/var
    f=sp.together((1-var**2)*expr/var)
    num,den=sp.fraction(f)
    # poles inside: var=0 and var=t_i (from 1/(1-t_i/var) = var/(var-t_i))
    poles=[0]+list(t)
    tot=0
    for p0 in poles:
        tot+=sp.residue(f,var,p0)
    return sp.simplify(tot)
Pz=1
for ti in t: Pz*=1/((1-ti*z)*(1-ti/z))
Pw=Pz.subs(z,w)
D=(z+1/z-w-1/w)**(2*R)
# expand D and integrate termwise: D = sum_k binom (-1)^k A^(2R-k) B^k
A=z+1/z; B=w+1/w
tot=0
Mz=[su2_int(sp.together(Pz*A**k),z) for k in range(2*R+1)]
for k in range(2*R+1):
    tot+=sp.binomial(2*R,k)*(-1)**k*Mz[2*R-k]*Mz[k]   # Mw = Mz by symmetry
Q=sp.factor(sp.simplify(tot))
print("Q_%d(t) for J=%d:"%(R,J)); print(Q)
num,den=sp.fraction(sp.together(Q))
print("numerator:",sp.expand(num)); print("denominator:",sp.factor(den))
