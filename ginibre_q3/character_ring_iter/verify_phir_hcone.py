#!/usr/bin/env python3
# phi_r(kappa) = (1/2) E_{sc x sc}[(x-y)^{2r} prod_i h_{kappa_i}(x,y)] >= 0 on the whole H-cone (FM3 note, Euler-characteristic section).
# usage: verify_phir_hcone.py Nmax r1,r2,...   (exact integers; prints per-level counts of negatives)
# phi_r(kappa) = 1/2 E_{sc x sc}[(x-y)^{2r} prod h_{kappa_i}(x,y)], h_k = sum_{a+b=k} U_a(x)U_b(y)
import sys
from math import comb
from functools import lru_cache
from collections import defaultdict
@lru_cache(None)
def ballot(a,n):  # E[U_a x^n] for semicircle = # walks 0->a length n on N
    if n<a or (n-a)%2: return 0
    return comb(n,(n-a)//2)-comb(n,(n-a)//2-1) if (n-a)//2>=1 else 1
@lru_cache(None)
def pair_ab(a,b,r):
    return sum(comb(2*r,j)*(-1)**j*ballot(a,2*r-j)*ballot(b,j) for j in range(2*r+1))
def mul_h(F,k):
    G=defaultdict(int)
    for (a,b),c in F.items():
        for i in range(k+1):
            i2=k-i
            for s in range(min(a,i)+1):
                for t in range(min(b,i2)+1):
                    G[(a+i-2*s,b+i2-2*t)]+=c
    return G
def partitions(n,maxpart=None):
    if maxpart is None: maxpart=n
    if n==0: yield (); return
    for k in range(min(n,maxpart),0,-1):
        for rest in partitions(n-k,k): yield (k,)+rest
R=[int(v) for v in sys.argv[2].split(',')]; Nmax=int(sys.argv[1])
cache={():{(0,0):1}}
def G(kap):
    if kap in cache: return cache[kap]
    g=mul_h(G(kap[1:]),kap[0]); cache[kap]=g; return g
for r in R:
    neg=[];cnt=0;zero=0
    for N in range(0,Nmax+1):
        for kap in partitions(N):
            F=G(kap)
            v=sum(c*pair_ab(a,b,r) for (a,b),c in F.items() if a+b<=2*r)
            cnt+=1
            if v<0: neg.append((kap,v//2))
            if v==0: zero+=1
    print('r',r,'kappas',cnt,'zeros',zero,'neg',len(neg), neg[:8],flush=True)
