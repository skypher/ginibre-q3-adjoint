#!/usr/bin/env python3
"""H-only FM3: ell_r(prod_j Sym^{k_j} C^4) for all multisets k with sum<=S, parts<=KMAX.
Reports zeros (with nonzero parity) and minimal normalized values.  usage: probe_su2_fm3_honly_support.py S KMAX RMAX"""
import sys, itertools
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
S,KMAX,RMAX=map(int,sys.argv[1:4])
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
H=lambda k:{(a,k-a):1 for a in range(k+1)}   # Sym^k C^4 = H_{k+1}
D2=cmul({(1,0):1,(0,1):-1},{(1,0):1,(0,1):-1})
def parts(n,mx):
    if n==0: yield (); return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,a): yield (a,)+r
zeros={r:[] for r in range(RMAX+1)}; neg=0; cnt=0
for tot in range(0,S+1,2):
    for k in parts(tot,KMAX):
        F={(0,0):1}
        for kj in k: F=cmul(F,H(kj))
        R=F
        vals=[]
        for r in range(RMAX+1):
            if r: R=cmul(R,D2)
            vals.append(R.get((0,0),0))
        cnt+=1
        for r,v in enumerate(vals):
            if v<0: neg+=1; print("NEG",k,r,v)
            if v==0 and r>=1: zeros[r].append(k)
print("multisets=%d negatives=%d"%(cnt,neg))
for r in range(1,RMAX+1):
    print("r=%d zeros (%d): %s"%(r,len(zeros[r]),zeros[r][:40]))
