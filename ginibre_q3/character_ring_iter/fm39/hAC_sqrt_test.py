# Test: is p = IWHT(sqrt(WHT(f))) >= 0 for f(S) = m(S) m(S^c)?  (canonical H_AC square root)
import itertools, numpy as np, sys, random
from functools import lru_cache
def mult_trivial(labels):
    # number of SU(2) invariants in tensor of V_n, n in labels (dims n+1)
    d={0:1}
    for n in labels:
        e={}
        for j,c in d.items():
            for k in range(abs(j-n), j+n+1, 2):
                e[k]=e.get(k,0)+c
        d=e
    return d.get(0,0)
def wht(v):
    v=v.astype(float).copy(); h=1; n=len(v)
    while h<n:
        for i in range(0,n,2*h):
            a=v[i:i+h].copy(); b=v[i+h:i+2*h].copy()
            v[i:i+h]=a+b; v[i+h:i+2*h]=a-b
        h*=2
    return v
def test(lam):
    L=len(lam); N=1<<L
    m=np.zeros(N)
    for S in range(N):
        m[S]=mult_trivial([lam[i] for i in range(L) if S>>i&1])
    f=np.array([m[S]*m[(N-1)^S] for S in range(N)])
    F=wht(f)
    if F.min()< -1e-9: return ('FM3FAIL',F.min())
    p=wht(np.sqrt(np.maximum(F,0)))/N
    return ('ok',p.min(), p.max())
random.seed(1)
worst=[]
cases=0
for L in range(2,11):
    for lam in itertools.combinations_with_replacement(range(1,7),L):
        if sum(lam)%2: continue
        if L>=8 and random.random()>0.05: continue
        r=test(lam); cases+=1
        if r[0]!='ok' or r[1]< -1e-9:
            worst.append((r[1],lam))
worst.sort()
print("cases",cases,"negative sqrt",len(worst))
print(worst[:15])
