# "B-lite": is f_lambda(S) = m(S) m(S^c) a nonnegative combination of indicators of the special subspaces
#   W_T = { S : S or its complement is a subset X of T with |X cap Odd| even }  (Odd = positions of odd labels)   (T subset of [L]; W_T is the span of {e_i + e_j : i, j in T}
#   together with <1>, taken mod <1>)?  If yes for all lists, B holds with a transparent family.  LP (scipy) + exact recheck.
import sys, random, numpy as np
from functools import lru_cache
from fractions import Fraction as Fr
from scipy.optimize import linprog
@lru_cache(None)
def fusion(labels):
    b={0:1}
    for n in labels:
        out={}
        for p,v in b.items():
            for q in range(abs(p-n),p+n+1,2): out[q]=out.get(q,0)+v
        b=out
    return b
def inv(l): return fusion(tuple(sorted(l))).get(0,0)
def test(labels):
    L=len(labels); full=(1<<L)-1; ODD=sum(1<<i for i in range(L) if labels[i]%2)
    reps=[S for S in range(1<<L) if S< (full^S)]            # one rep per complement class
    idx={S:k for k,S in enumerate(reps)}
    f=[inv([labels[i] for i in range(L) if S>>i&1])*inv([labels[i] for i in range(L) if not S>>i&1]) for S in reps]
    cols=[]
    for T in range(1<<L):
        col=np.zeros(len(reps))
        for S in reps:
            for X in (S, full^S):
                if (X & ~T)==0 and bin(X & ODD).count('1')%2==0: col[idx[S]]=1; break
        cols.append(col)
    A=np.array(cols).T
    res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=np.array(f,dtype=float),bounds=(0,None),method='highs')
    return res.status==0, f
seed,N,LMIN,LMAX,NMAX=[int(x) for x in sys.argv[1:6]]
rng=random.Random(seed); ok=bad=0; badl=[]
for t in range(N):
    L=rng.randint(LMIN,LMAX); lab=tuple(sorted(rng.randint(1,NMAX) for _ in range(L)))
    if sum(lab)%2: lab=lab[:-1]+(lab[-1]+1,)
    feas,f=test(lab)
    if max(f)==0: continue
    if feas: ok+=1
    else:
        bad+=1
        if len(badl)<8: badl.append(lab)
print('B-lite feasible',ok,'infeasible',bad,'examples of infeasible:',badl)
