#!/usr/bin/env python3
"""n-wise refinement: <M,E_n> = (1/4) [u^{-n}] CT_{z,w} chi_M(t)|Delta_SU4(t)|^2, t=(z,u/z,w,1/(uw)),
for M = tensor Sym^{k_j}(C^4). Tests nonnegativity of every n-coefficient.  usage: probe_su2_fm3_nwise.py MAXSUM MAXPART MAXF"""
import sys, itertools
from collections import defaultdict
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
MS,MP,MF=map(int,sys.argv[1:4])
# Laurent polys in (u,z,w) as dict (a,b,c)->coef
def mul(P,Q):
    R=defaultdict(int)
    for k1,v1 in P.items():
        for k2,v2 in Q.items(): R[(k1[0]+k2[0],k1[1]+k2[1],k1[2]+k2[2])]+=v1*v2
    return {k:v for k,v in R.items() if v}
t=[(0,1,0),(1,-1,0),(0,0,1),(-1,0,-1)]   # exponents of u,z,w for t1..t4
def mono(e): return {e:1}
def h(k):
    # complete homogeneous of degree k in t1..t4
    P=defaultdict(int)
    for c in itertools.product(range(k+1),repeat=4):
        if sum(c)!=k: continue
        e=tuple(sum(c[i]*t[i][d] for i in range(4)) for d in range(3))
        P[e]+=1
    return dict(P)
Delta2={(0,0,0):1}
for i in range(4):
    for j in range(4):
        if i!=j:
            e=tuple(t[i][d]-t[j][d] for d in range(3))
            Delta2=mul(Delta2,{(0,0,0):1, e:-1})
def parts(n,mx,ml):
    if n==0: yield (); return
    if ml==0: return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,a,ml-1): yield (a,)+r
H={}
bad=0; tot=0
for s in range(0,MS+1,2):
    for ks in parts(s,MP,MF):
        P=Delta2
        for k in ks:
            if k not in H: H[k]=h(k)
            P=mul(P,H[k])
        coeffs=defaultdict(int)
        for (a,b,c),v in P.items():
            if b==0 and c==0: coeffs[-a]+=v
        coeffs={n:v for n,v in coeffs.items() if v}
        assert all(v%4==0 for v in coeffs.values()), (ks,coeffs)
        cs={n:v//4 for n,v in coeffs.items()}
        tot+=1
        if any(v<0 for v in cs.values()):
            bad+=1
            if bad<=10: print("NEG",ks,dict(sorted(cs.items())))
        if len(ks)<=2 and s<=4: print(ks,dict(sorted(cs.items())),"sum",sum(cs.values()))
print("tested",tot,"with a negative n-coefficient:",bad)
