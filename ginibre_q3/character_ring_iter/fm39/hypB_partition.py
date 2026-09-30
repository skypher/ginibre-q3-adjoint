# "B-partition": is f_lambda a nonnegative combination of partition-code indicators
#   W_pi = { S : sum of labels in S cap B even for every block B },  pi a set partition of [L] into blocks of even total label
# (W_pi contains the all-ones vector, so it is a subspace of G_L = F_2^L/<1>)?  LP test on random lists.
import sys, random, numpy as np
from functools import lru_cache
from scipy.optimize import linprog
exec(open('hypB_lite.py').read().split("def test(labels):")[0])
def even_partitions(elems, lab=None):
    # all set partitions whose blocks have even total label (labels via global LAB)
    if not elems: yield []; return
    first=elems[0]; rest=elems[1:]
    from itertools import combinations
    for k in range(0,len(rest)+1):
        for comb in combinations(rest,k):
            if (LAB[first]+sum(LAB[c] for c in comb))%2: continue
            block=(first,)+comb
            remaining=[x for x in rest if x not in comb]
            for p in even_partitions(remaining): yield [block]+p
def test(labels):
    global LAB; LAB=labels
    L=len(labels); full=(1<<L)-1
    reps=[S for S in range(1<<L) if S<(full^S)]; idx={S:k for k,S in enumerate(reps)}
    f=[inv([labels[i] for i in range(L) if S>>i&1])*inv([labels[i] for i in range(L) if not S>>i&1]) for S in reps]
    cols=[]
    for pi in even_partitions(list(range(L))):
        masks=[sum(1<<i for i in B) for B in pi]
        col=np.zeros(len(reps))
        for S in reps:
            if all(sum(labels[i] for i in range(L) if (S&mk)>>i&1)%2==0 for mk in masks): col[idx[S]]=1
        cols.append(col)
    A=np.array(cols).T
    res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=np.array(f,dtype=float),bounds=(0,None),method='highs')
    return res.status==0, f, len(cols)
seed,N,LMIN,LMAX,NMAX=[int(x) for x in sys.argv[1:6]]
rng=random.Random(seed); ok=bad=0; badl=[]
for t in range(N):
    L=rng.randint(LMIN,LMAX)
    lab=tuple(sorted(rng.randint(1,NMAX) for _ in range(L)))
    if sum(lab)%2: lab=lab[:-1]+(lab[-1]+1,)
    feas,f,nc=test(lab)
    if max(f)==0: continue
    if feas: ok+=1
    else:
        bad+=1
        if len(badl)<8: badl.append(lab)
print('B-partition feasible',ok,'infeasible',bad,'examples:',badl)
