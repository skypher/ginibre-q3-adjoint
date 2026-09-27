#!/usr/bin/env python3
"""For symmetric 0/1 kernels F supported on an antidiagonal a+b=n, test admissibility:
column 0 of F*D_1^m >= 0 for m<=M. usage: probe_su2_fm_antidiagonal.py NMAX M"""
import sys, itertools
if any(a in ('-h','--help') for a in sys.argv[1:]): print(__doc__); sys.exit(0)
NMAX, M = map(int, sys.argv[1:3])
def cg(a,p): return range(abs(a-p), a+p+1, 2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
D1={(1,0):1,(0,1):-1}
def admissible(F):
    R=F
    for m in range(M+1):
        if m: R=cmul(R,D1)
        if any(v<0 for (a,b),v in R.items() if b==0): return m
    return None
for n in range(1,NMAX+1):
    orbits=[(a,n-a) for a in range(n,(n-1)//2,-1)]  # a>=b
    for mask in range(1,1<<len(orbits)):
        F={}
        sel=[]
        for i,(a,b) in enumerate(orbits):
            if mask>>i&1:
                sel.append((a,b)); F[(a,b)]=1; F[(b,a)]=1
        r=admissible(F)
        print("n=%d orbits=%s : %s"%(n,sel,"OK" if r is None else "fails at m=%d"%r))
