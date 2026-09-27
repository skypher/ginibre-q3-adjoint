#!/usr/bin/env python3
"""Test whether restrictions of Sp(4) irreps to SU(2)xSU(2) are admissible plus
factors: column 0 of F * D_1^m >= 0 for F = products of generators from
{S_p} u {Res V_lam^{Sp4}}.  usage: probe_su2_fm_sp4_cone.py L N M"""
import sys, itertools
if any(a in ('-h','--help') for a in sys.argv[1:]):
    print(__doc__); sys.exit(0)
L, N, M = map(int, sys.argv[1:4])
def lmul(A,B):
    C={}
    for (a1,b1),c1 in A.items():
        for (a2,b2),c2 in B.items():
            k=(a1+a2,b1+b2); C[k]=C.get(k,0)+c1*c2
    return {k:v for k,v in C.items() if v}
def sp4_char_times_denominator(l1,l2):
    # numerator of Weyl formula for C2: det[t_j^{m_i} - t_j^{-m_i}], m=(l1+2,l2+1); t1=z, t2=w
    m1,m2=l1+2,l2+1
    def a(e1,e2):  # (z^e1 - z^-e1)(w^e2 - w^-e2)
        return {(e1,e2):1,(-e1,e2):-1,(e1,-e2):-1,(-e1,-e2):1}
    num = a(m1,m2); 
    for k,v in a(m2,m1).items(): num[k]=num.get(k,0)-v
    return {k:v for k,v in num.items() if v}
import sympy as _sp
_z,_w=_sp.symbols('z w')
def restrict(l1,l2):
    num = sp4_char_times_denominator(l1,l2)
    N=sum(v*_z**i*_w**j for (i,j),v in num.items())
    # chi*(z-1/z)(w-1/w) = num / (z+1/z-w-1/w)
    q=_sp.cancel(N/(_z+1/_z-_w-1/_w))
    q=_sp.expand(q)
    P=_sp.Poly(_sp.expand(q*_z**60*_w**60),_z,_w)
    out={}
    for (i,j),c in P.terms():
        i-=60; j-=60
        if i>=1 and j>=1: out[(i-1,j-1)]=int(c)
    return out
def cg(a,p): return range(abs(a-p), a+p+1, 2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
gens=[('S%d'%p,{(p,0):1,(0,p):1}) for p in range(1,L+1)]
for l1 in range(0,L+1):
    for l2 in range(0,l1+1):
        if l1==0: continue
        R=restrict(l1,l2)
        assert all(v>0 for v in R.values()), (l1,l2,R)
        gens.append(('Sp(%d,%d)'%(l1,l2),R))
print("sanity Sym^2 = ", restrict(2,0), " Lambda2_0 =", restrict(1,1), flush=True)
D1={(1,0):1,(0,1):-1}
fails=0; tested=0
for n in range(1,N+1):
    for combo in itertools.combinations_with_replacement(range(len(gens)),n):
        F={(0,0):1}
        for i in combo: F=cmul(F,gens[i][1])
        R=F
        for m in range(0,M+1):
            if m: R=cmul(R,D1)
            tested+=1
            neg={a:v for (a,b),v in R.items() if b==0 and v<0}
            if neg:
                fails+=1
                if fails<=12: print("FAIL m=%d F=%s neg=%s"%(m,[gens[i][0] for i in combo],neg),flush=True)
print("done gens=%d tested=%d fails=%d"%(len(gens),tested,fails),flush=True)
