#!/usr/bin/env python3
"""J-fraction (Jacobi) coefficients of the moment sequence M_k(s)=<exp(s chi_p) chi_1^k>/<exp(s chi_p)>
as power series in s; checks coefficientwise nonnegativity of beta_n, gamma_n.
usage: probe_su2_fm3_jfraction.py P SORDER NLEVELS"""
import sys
from fractions import Fraction
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
P,SO,NL=map(int,sys.argv[1:4])
KMAX=2*NL+2
def cg(a,p): return range(abs(a-p),a+p+1,2)
# mult[n][a] = multiplicity of V_a in V_p^{n}
mult=[{0:1}]
for n in range(SO):
    d={}
    for a,v in mult[-1].items():
        for c in cg(a,P): d[c]=d.get(c,0)+v
    mult.append(d)
fact=[1]
for n in range(1,SO+1): fact.append(fact[-1]*n)
def ballot(a,k):
    # multiplicity of V_0 in V_a (x) V_1^k
    d={a:1}
    for _ in range(k):
        nd={}
        for c,v in d.items():
            for e in cg(c,1): nd[e]=nd.get(e,0)+v
        d=nd
    return d.get(0,0)
# series ops
def sm(a,b):
    c=[Fraction(0)]*SO
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:SO-i]): c[i+j]+=x*y
    return c
def sinv(a):
    c=[Fraction(0)]*SO; c[0]=1/a[0]
    for n in range(1,SO): c[n]=-sum(a[k]*c[n-k] for k in range(1,n+1))/a[0]
    return c
M=[]
for k in range(KMAX+1):
    M.append([Fraction(sum(v*ballot(a,k) for a,v in mult[n].items()),fact[n]) for n in range(SO)])
inv0=sinv(M[0]); m=[sm(Mk,inv0) for Mk in M]
# Chebyshev algorithm (modified via monic orthogonal polynomials with series coefficients)
# represent polynomials as lists of series coefficients in x
def pmul_x(p): return [[Fraction(0)]*SO]+p
def padd(p,q,sg=1):
    n=max(len(p),len(q)); z=[Fraction(0)]*SO
    return [[ (p[i][j] if i<len(p) else 0)+sg*(q[i][j] if i<len(q) else 0) for j in range(SO)] for i in range(n)]
def pscal(p,c): return [sm(pi,c) for pi in p]
def funct(p,shift=0):  # L[x^shift p(x)]
    tot=[Fraction(0)]*SO
    for i,pi in enumerate(p): tot=[t+u for t,u in zip(tot,sm(pi,m[i+shift]))]
    return tot
one=[Fraction(1)]+[Fraction(0)]*(SO-1)
pprev=[[Fraction(0)]*SO]; pcur=[one]
norm_prev=None
for n in range(NL):
    h=funct(pscal(pcur,[Fraction(1)]+[0]*(SO-1)) if False else pcur,0)
    # h_n = L[p_n^2] ; compute via L[p_n * p_n] using x-multiplication
    # compute L[p_n(x)^2]: expand
    sq=[[Fraction(0)]*SO for _ in range(2*len(pcur)-1)]
    for i,a in enumerate(pcur):
        for j,b in enumerate(pcur):
            sq[i+j]=[u+v for u,v in zip(sq[i+j],sm(a,b))]
    hn=funct(sq); 
    xsq=[[Fraction(0)]*SO]+sq
    gam=sm(funct(xsq),sinv(hn))
    beta=sm(hn,sinv(norm_prev)) if norm_prev is not None else None
    ok=lambda s: all(x>=0 for x in s)
    print("n=%d gamma_n: %s  %s"%(n,[str(x) for x in gam[:6]],'' if ok(gam) else 'NEG'))
    if beta is not None: print("     beta_n : %s  %s"%([str(x) for x in beta[:6]],'' if ok(beta) else 'NEG'))
    # recurrence p_{n+1} = (x-gam) p_n - beta p_{n-1}
    nxt=padd(pmul_x(pcur),pscal(pcur,gam),-1)
    if beta is not None: nxt=padd(nxt,pscal(pprev,beta),-1)
    pprev,pcur=pcur,nxt; norm_prev=hn
