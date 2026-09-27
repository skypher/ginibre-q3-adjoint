#!/usr/bin/env python3
"""Test candidate inductive invariants on Sp(4)-multiplicities y_lambda of products of S_p-hat.
Sp(4) multiplicities computed from the doubled character via Res^{-1} stencils.
usage: probe_su2_fm3_row_domination.py L N"""
import sys, itertools
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
L,N=map(int,sys.argv[1:3])
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
def sp4mult(F):
    """Res^{-1}: symmetric doubled char -> dict lambda->mult, via D_1*F antisymmetric = sum x_lam (V_{l1+1}xV_{l2}-V_{l2}xV_{l1+1})"""
    D=cmul(F,{(1,0):1,(0,1):-1})
    y={}
    for (c,d),v in D.items():
        if c>d: y[(c-1,d)]=y.get((c-1,d),0)+v
    return {k:v for k,v in y.items() if v}
fails={'I_p (b=1)':0,'GND all odd b':0,'GND even b':0}
ex={}
total=0
for n in range(1,N+1):
    for combo in itertools.combinations_with_replacement(range(1,L+1),n):
        F={(0,0):1}
        for p in combo: F=cmul(F,{(p,0):1,(0,p):1})
        y=sp4mult(F); total+=1
        g=lambda a,b: y.get((a,b),0) if a>=b>=0 else 0
        amax=max([a for a,b in y]+[0])+2
        for a in range(0,amax+1):
            for b in range(0,a+1):
                if b==0: continue
                lhs=g(a,b); rhs=g(a+1,b-1)+g(a-1,b-1)
                key='I_p (b=1)' if b==1 else ('GND all odd b' if b%2==1 else 'GND even b')
                if lhs>rhs:
                    fails[key]+=1
                    ex.setdefault(key,(combo,(a,b),lhs,rhs))
print("products=%d"%total)
for k,v in fails.items(): print(k,"violations:",v,"example:",ex.get(k))
