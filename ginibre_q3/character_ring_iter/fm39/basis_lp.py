# Is the three-factor quadratic form (closed form in c) a nonnegative combination of basic forms
# {OL: D_k - D_(k+1), 2k>N} U {W: P_C(x)} U {(E): D_j - D_i +/- W_ij}, as forms on anti-reciprocal sequences (c_(N-k) = -c_k)?
import sys, itertools
import numpy as np
from scipy.optimize import linprog
exec(open('recip_test.py').read().split("random.seed(11)")[0].split("exec(open('split3.py')")[0])
def fold(N,sign):
    # free variables: c_k for k < N/2 (and middle forced 0 if sign=-1 and N even)
    free=[k for k in range(N+1) if k < N-k]
    mid=[k for k in range(N+1) if k==N-k]
    idx={}
    for t,k in enumerate(free): idx[k]=(t,1); idx[N-k]=(t,sign)
    nf=len(free)
    if mid and sign==1: idx[mid[0]]=(nf,1); nf+=1
    elif mid: idx[mid[0]]=None
    return idx,nf
def form(pairs,N,idx,nf):   # pairs: list of (coef, i, j) meaning coef*c_i*c_j
    Mx=np.zeros((nf,nf))
    for co,i,j in pairs:
        if not(0<=i<=N and 0<=j<=N): continue
        a=idx.get(i); b=idx.get(j)
        if a is None or b is None: continue
        v=co*a[1]*b[1]; Mx[a[0],b[0]]+=v/2; Mx[b[0],a[0]]+=v/2
    return Mx
def Dp(k): return [(1,k,k),(-1,k-1,k+1)]
def Tp(x,C): return [(1,x,x+C),(-1,x-1,x+C+1)]
def Wp(i,j): return [(1,i-1,j),(1,i+1,j),(-1,i,j-1),(-1,i,j+1)]
def neg(p): return [(-c,i,j) for c,i,j in p]
def P_C(x,C): return sum((Dp(k) for k in range(x,x+C+1)),[])+neg(Tp(x,C))
def phi3_pairs(N,A,B,C):
    al=(N+A-B-C)//2; tau=al+B+1
    Ac=lambda x:[(1,x),(-1,x+C)]; Bc=lambda x:[(1,x-1),(-1,x+C+1)]
    def prod(u,v,s): return [(s*a*b,i,j) for a,i in u for b,j in v]
    return P_C(al,C)+neg(P_C(tau,C))+neg(prod(Ac(al),Bc(tau),1))+prod(Bc(al),Ac(tau),1)
from collections import Counter
st=Counter(); ex=[]
for r in (2,3):
    e=2*r-3
    for a in range(2,10):
        N=a+e; idx,nf=fold(N,-1)
        basics=[]
        for k in range(N+1):
            if 2*k>N: basics.append(form(Dp(k)+neg(Dp(k+1)),N,idx,nf))
        for x in range(0,N+1):
            for C in range(1,N-x+1): basics.append(form(P_C(x,C),N,idx,nf))
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                for s in (1,-1): basics.append(form(Dp(j)+neg(Dp(i))+[(s*co,p,q) for co,p,q in Wp(i,j)],N,idx,nf))
        iu=np.triu_indices(nf); Amat=np.array([Bm[iu] for Bm in basics]).T
        for C in range(3,N+2):
            for B in range(C,N+3):
                for A in range(B,B+C+1):
                    if (N+A-B-C)%2: continue
                    al=(N+A-B-C)//2
                    if al+B>N: continue   # gamma<=N branch only
                    tg=form(phi3_pairs(N,A,B,C),N,idx,nf)[iu]
                    res=linprog(np.zeros(Amat.shape[1]),A_eq=Amat,b_eq=tg,bounds=[(0,None)]*Amat.shape[1],method='highs')
                    st['feasible' if res.status==0 else 'infeasible']+=1
                    if res.status!=0 and len(ex)<5: ex.append((r,a,A,B,C))
print(dict(st),'first infeasible (r,a,A,B,C):',ex)
