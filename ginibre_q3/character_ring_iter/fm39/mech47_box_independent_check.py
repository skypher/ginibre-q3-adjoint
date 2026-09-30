# Main agent: independent (u,v)-expansion evaluator vs FM-MECH47 IBP table on random box entries.  Usage: python3 mech47_box_independent_check.py SEED COUNT
import random, sys, importlib.util
from fractions import Fraction as Q
from math import factorial, comb
from functools import lru_cache
@lru_cache(None)
def mu(m,k): return Q(2*factorial(2*m)*factorial(2*m+1)*factorial(2*k)*factorial(2*k+1), factorial(m)**2*factorial(k)**2*factorial(m+k+1)*factorial(m+k+2))
def pmul(a,b):
    r={}
    for (i,j),c in a.items():
        for (k,l),d in b.items(): r[(i+k,j+l)]=r.get((i+k,j+l),0)+c*d
    return r
@lru_cache(None)
def ppow(name,n):
    base={'Z':{(1,0):Q(1,2),(0,1):Q(1,2),(0,0):Q(-2)},'H':{(1,0):Q(3,4),(0,1):Q(1,4),(0,0):Q(-2)},'G':{(1,0):Q(1,4),(0,1):Q(3,4),(0,0):Q(-2)}}[name]
    if n==0: return {(0,0):Q(1)}
    h=ppow(name,n//2); r=pmul(h,h)
    return pmul(r,base) if n%2 else r
def M(m,k,b,al,ga):
    P=pmul(pmul(ppow('Z',b),ppow('H',al)),ppow('G',ga))
    return sum(c*mu(m+i,k+j) for (i,j),c in P.items())
# load ceres table evaluator
src=open('/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39/mech47_sector123_repro.py').read()
pre=src.split("coeffs={")[0]
ns={}; _sv=sys.argv; sys.argv=[sys.argv[0]]; exec(pre.split("large_L=")[0]+"\n"+pre[pre.index("@lru_cache(None)\ndef cat"):],ns)
sys.argv=_sv; random.seed(int(sys.argv[1])); cnt=0; worst=None
for _ in range(int(sys.argv[2])):
    N=random.randint(2,19); L=random.randint(2,min(21,2*N)); al=random.randint(0,L); ga=L-al; b=random.randint(0,41); m=random.randint(0,N); r=N-m
    T=ns['mixed_table'](m,r,min(21,2*N),62)
    cc=[(j,v) for j,v in enumerate(ns['row'](al,ga)) if v]
    val=sum(v*T[j][b+L-j] for j,v in cc)
    ind=M(m,r,b,al,ga)
    assert ind==val,(2*m,2*r,b,al,ga,ind,val)
    ratio=ind/mu(m,r)
    if worst is None or ratio<worst[0]: worst=(ratio,2*m,2*r,b,al,ga)
    cnt+=1
print("independent agreements",cnt,"min value/mu",float(worst[0]),worst[1:])
