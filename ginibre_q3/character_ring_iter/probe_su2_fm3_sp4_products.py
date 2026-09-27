#!/usr/bin/env python3
"""Sp(4)-decomposition sign pattern of products of S_p-hat (virtual Sp(4) characters).
usage: probe_su2_fm3_sp4_products.py L N"""
import sys, itertools, sympy as sp
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
L,N=map(int,sys.argv[1:3])
z,w=sp.symbols('z w')
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
def decompose(F):
    tot=0
    for (A,B),c in F.items(): tot+=c*(z**(A+1)-z**-(A+1))*(w**(B+1)-w**-(B+1))
    num=sp.expand(tot*(z+1/z-w-1/w))
    P=sp.Poly(sp.expand(num*z**200*w**200),z,w)
    out={}
    for (i,j),c in P.terms():
        i-=200;j-=200
        if i>j>0: out[(i-2,j-1)]=int(c)
    return out
bad_even=0; bad_odd_pos=0; total=0
for n in range(1,N+1):
    for combo in itertools.combinations_with_replacement(range(1,L+1),n):
        F={(0,0):1}
        for p in combo: F=cmul(F,{(p,0):1,(0,p):1})
        D=decompose(F); total+=1
        neg_even={k:v for k,v in D.items() if v<0 and k[1]%2==0}
        pos_odd={k:v for k,v in D.items() if v>0 and k[1]%2==1}
        if neg_even: bad_even+=1
        if pos_odd: bad_odd_pos+=1
        if n<=2 or neg_even:
            print("S%s -> triv=%d  negatives=%s"%(list(combo),D.get((0,0),0),{k:v for k,v in sorted(D.items()) if v<0}))
print("products=%d  with negative even-lambda2 coeff: %d ; with positive odd-lambda2 coeff: %d"%(total,bad_even,bad_odd_pos))
