# Necessary condition for Hypothesis B (FM-MECH30): the sign-pattern function F(T) = sum_S (-1)^|S cap T| m(S) m(S^c)
# (FM3 value up to the factor 2 for the word with minus set T; F = 0 at odd T) must be a nonnegative combination of
# subgroup indicators on the even patterns, hence F(T1+T2) + F(0) >= F(T1) + F(T2) for all T1, T2 (and F >= 0).
# Walsh-Hadamard transform of g(S) = m(S) m(S^c); exact integers.  Random label lists beyond the FM-MECH30 screens.
import sys, random, numpy as np
from functools import lru_cache
@lru_cache(None)
def fusion(labels):
    b={0:1}
    for n in labels:
        out={}
        for p,v in b.items():
            for q in range(abs(p-n),p+n+1,2): out[q]=out.get(q,0)+v
        b=out
    return b
def inv(labels): return fusion(tuple(sorted(labels))).get(0,0)
def F_of(labels):
    L=len(labels); g=np.zeros(1<<L,dtype=object)
    for S in range(1<<L):
        A=[labels[i] for i in range(L) if S>>i&1]; B=[labels[i] for i in range(L) if not S>>i&1]
        g[S]=inv(A)*inv(B)
    F=g.copy(); h=1
    while h<(1<<L):
        for i in range(0,1<<L,2*h):
            a=F[i:i+h].copy(); b=F[i+h:i+2*h].copy(); F[i:i+h]=a+b; F[i+h:i+2*h]=a-b
        h*=2
    return F
seed,NL,LMIN,LMAX,NMAX=[int(x) for x in sys.argv[1:6]]
rng=random.Random(seed); bad=[]; tested=0; negF=0
for t in range(NL):
    L=rng.randint(LMIN,LMAX); labels=tuple(sorted(rng.randint(1,NMAX) for _ in range(L)))
    if sum(labels)%2: labels=labels[:-1]+(labels[-1]+1,)
    F=F_of(labels)
    even=[T for T in range(1<<L) if bin(T).count('1')%2==0]
    Fe=np.array([int(F[T]) for T in even],dtype=object)
    if min(Fe)<0: negF+=1
    idx={T:k for k,T in enumerate(even)}
    F0=int(F[0]); worst=None
    for T1 in even:
        for T2 in even:
            if T2<T1: continue
            v=int(F[T1^T2])+F0-int(F[T1])-int(F[T2])
            if worst is None or v<worst[0]: worst=(v,T1,T2)
    tested+=1
    if worst[0]<0 and len(bad)<8: bad.append((labels,worst,int(F[worst[1]]),int(F[worst[2]]),int(F[worst[1]^worst[2]]),F0))
    print('list',labels,'min slack F(T1+T2)+F(0)-F(T1)-F(T2) =',worst[0],flush=True)
print('lists tested',tested,'with negative F:',negF,'violations:',len(bad)); print(bad)
