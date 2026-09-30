from math import comb
from fractions import Fraction as Q
import numpy as np, itertools, sys
from scipy.optimize import linprog
def riordan(K):
    st={0:1}; out=[1]
    for _ in range(K):
        nx={}
        for h,v in st.items():
            for j in (h+2,h,h-2):
                if j<0 or (j==h and h==0): continue
                nx[j]=nx.get(j,0)+v
        st=nx; out.append(st.get(0,0))
    return out
R=riordan(40)
def partitions_into_blocks(n,minb=2):
    # integer partitions (block sizes >= minb) of m<=n
    res=[]
    def rec(rem,maxp,cur):
        res.append(tuple(cur))
        for p in range(min(rem,maxp),minb-1,-1):
            rec(rem-p,p,cur+[p])
    rec(n,n,[]); return res
def code_profile(M,blocks,evenA):
    # coordinates: first evenA coords carry even-weight code; then disjoint blocks with block indicators
    gens=[]
    pos=0
    for i in range(evenA-1): gens.append((1<<pos)|(1<<(pos+1))); pos+=1
    if evenA>0: pos+=1
    for b in blocks:
        if pos+b>M: return None
        gens.append(((1<<b)-1)<<pos); pos+=b
    H={0}
    for g in gens: H|={x^g for x in H}
    full=(1<<M)-1
    Hq={min(x,x^full) for x in H}
    cnt=[0]*(M//2+1)
    for x in Hq:
        w=bin(x).count('1'); cnt[min(w,M-w)]+=1
    sizes=[(comb(M,w) if 2*w!=M else comb(M,w)//2) for w in range(M//2+1)]
    return tuple(Q(c,s) for c,s in zip(cnt,sizes))
for M in range(3,15):
    t=[Q(R[w]*R[M-w]) for w in range(M//2+1)]
    cols=set()
    for evenA in [0]+list(range(2,M+1)):
        for blocks in partitions_into_blocks(M-evenA):
            pr=code_profile(M,list(blocks),evenA)
            if pr is not None: cols.add(pr)
    cols=sorted(cols)
    A=np.array([[float(c[i]) for c in cols] for i in range(len(t))]); b=np.array([float(x) for x in t])
    r=linprog(np.zeros(len(cols)),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    print(M,'cols',len(cols),'feasible' if r.status==0 else 'INFEASIBLE',flush=True)
