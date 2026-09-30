# EM2 at level r: phi_r(kappa) >= phi_r(kappa - e_{s-1} - e_s)  (two smallest parts), H-only words.
# phi_r(kappa) = (1/2) [U0 x U0] ( h_kappa * (x-y)^(2r) ) in the SU(2)xSU(2) character ring.
import sys, time
from functools import lru_cache
from collections import defaultdict
def mul1(A,B):
    # A,B: dict (a,b)->c ; Clebsch-Gordan in each factor
    C=defaultdict(int)
    for (a,b),x in A.items():
        for (c,d),y in B.items():
            xy=x*y
            for e in range(abs(a-c),a+c+1,2):
                for f in range(abs(b-d),b+d+1,2): C[(e,f)]+=xy
    return {k:v for k,v in C.items() if v}
def hk(k): return {(p,k-p):1 for p in range(k+1)}
R=int(sys.argv[1]); NMAX=int(sys.argv[2]); LMAX=int(sys.argv[3]) if len(sys.argv)>3 else 99
X={(1,0):1,(0,1):-1}
Q={(0,0):1}
for _ in range(2*R): Q=mul1(Q,X)
cacheP={}
def prod(kap):
    kap=tuple(sorted([k for k in kap if k],reverse=True))
    if kap in cacheP: return cacheP[kap]
    if not kap: P={(0,0):1}
    else: P=mul1(prod(kap[:-1]),hk(kap[-1]))
    cacheP[kap]=P; return P
def phi(kap):
    P=prod(kap)
    # [U0xU0](P*Q) = sum_{(a,b)} P[a,b]*Q[a,b]  (U's self-dual, orthonormal)
    return sum(v*Q.get(k,0) for k,v in P.items())  # = 2*phi_r
def parts(n,m,l):
    if n==0: yield (); return
    if l==0: return
    for k in range(min(n,m),0,-1):
        for r in parts(n-k,k,l-1): yield (k,)+r
fails=[];tot=0;minm=None;neg=[]
for n in range(2,NMAX+1,2):
    for lam in parts(n,n,LMAX):
        v=phi(lam)
        if v<0: neg.append(lam)
        if len(lam)<2: continue
        k2=list(lam); k2[-1]-=1; k2[-2]-=1
        d=v-phi(k2); tot+=1
        if d<0: fails.append((lam,d//2))
        elif minm is None or d<minm[0]: minm=(d//2,lam)
    print(f"[{time.strftime('%H:%M:%S')}] r={R} size {n}: EM2 tests {tot}, failures {len(fails)} {fails[:6]}; min margin {minm}; negative phi {len(neg)}",flush=True)
