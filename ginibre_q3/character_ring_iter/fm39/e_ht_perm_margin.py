# True (HT) margin from the closed form (HT*): rho* = (2V(D_j - D_i) / (|m_j n_i| + |n_j m_i|))^2 over case-(c) pairs
# (psi-sweep >= pi, q >= 2), rows a >= e+2, e >= 3.  Exact integers (2m, 2n), ratio in floats.  Replaces the finite-t
# search of e_ht_margin.py, which missed the limit t -> 0 at X_j = 0 (FM-SEC54).
import sys, math
from math import comb
from collections import defaultdict
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
AM=int(sys.argv[1])
jobs=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [(a,e) for a in range(5,AM+1) for e in range(3,a-1)]
best=defaultdict(lambda:(math.inf,None)); n=0
for (a,e) in jobs:
    c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
    D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
    Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
    M2=[2*sg*C_(k)-d*Bv[k] for k in range(N+3)]; L2=[(2*k-N)*Av[k] for k in range(N+3)]
    for j in range((N+1)//2,N+1):
        cj,Bj=C_(j),Bv[j]; imax=j
        for k in range(j+1,N+2):
            if cj*Bv[k]-Bj*C_(k)>0: imax=k
            else: break
        for i in range(imax+1,N+2):
            if i-j-1<=1: continue
            n+=1; lhs=8*V*(D[j]-D[i]); rhs=abs(M2[j]*L2[i])+abs(L2[j]*M2[i])
            r=min(lhs/rhs,1e100)**2 if rhs else math.inf
            cl=('e odd' if e%2 else 'e even')+(', X_j=0' if 2*j==N else '')
            if r<best[cl][0]: best[cl]=(r,(a,e,j,i,i-j-1))
print('case-(c) pairs',n)
for k,v in sorted(best.items()): print('%-16s min (2VE/perm)^2 = %.4f at (a,e,j,i,q) = %s'%(k,v[0],v[1]))
