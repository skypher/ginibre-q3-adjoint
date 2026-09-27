#!/usr/bin/env python3
"""Sp(4) multiplicities (m00,m11,m20) of Sym-products; test candidate linear inequalities.
usage: probe_su2_fm3_ratios.py MAXSUM MAXPART MAXFACTORS"""
import sys, itertools
from fractions import Fraction
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
MS,MP,MF=map(int,sys.argv[1:4])
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
H=lambda k:{(a,k-a):1 for a in range(k+1)}
D1={(1,0):1,(0,1):-1}
rows=[]
def parts(n,mx,maxlen):
    if n==0: yield (); return
    if maxlen==0: return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,a,maxlen-1): yield (a,)+r
for tot in range(0,MS+1,2):
    for ks in parts(tot,MP,MF):
        F={(0,0):1}
        for k in ks: F=cmul(F,H(k))
        Psi=cmul(F,D1)   # m_mu = Psi[(mu1+1,mu2)]
        m00=Psi.get((1,0),0); m11=Psi.get((2,1),0); m20=Psi.get((3,0),0)
        rows.append((ks,m00,m11,m20))
print("rows",len(rows))
viol={}
tests={'3m11<=m20+5m00':lambda a,b,c:3*b<=c+5*a,
       'm11<=m20':lambda a,b,c:b<=c,
       'm11<=m00+m20':lambda a,b,c:b<=a+c,
       '2m11<=m20+2m00':lambda a,b,c:2*b<=c+2*a,
       '3m11<=m20+3m00':lambda a,b,c:3*b<=c+3*a,
       '3m11<=m20+4m00':lambda a,b,c:3*b<=c+4*a,
       '3m11<=2m20+2m00':lambda a,b,c:3*b<=2*c+2*a,
       'm11<=2m00':lambda a,b,c:b<=2*a}
for name,f in tests.items():
    v=[(ks,a,b,c) for ks,a,b,c in rows if not f(a,b,c)]
    print("%-18s violations %4d  e.g. %s"%(name,len(v),v[:2]))
# extremal ratios
best=None
for ks,a,b,c in rows:
    if b>0:
        r=Fraction(c+5*a,3*b)
        if best is None or r<best[0]: best=(r,ks,a,b,c)
print("min (m20+5m00)/(3m11) =",best)
