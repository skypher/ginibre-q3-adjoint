# LP: is (E)_(q,+/-) at (j, i=j+q+1) a nonnegative combination of OL differences
#   delta^R_(j+t) = D_(j+t)(R c) - D_(j+t+1)(R c),  R = (1+z)^al (1-z)^be,  2t > deg R  (OL range for every j >= N/2)?
# Identity required as quadratic forms in a free sequence c (indices relative to j = 0).
import sys, itertools
import numpy as np
from scipy.optimize import linprog
def pmul(a,b):
    o=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): o[i+j]+=x*y
    return o
def Rpoly(al,be):
    R=[1]
    for _ in range(al): R=pmul(R,[1,1])
    for _ in range(be): R=pmul(R,[1,-1])
    return R
def lin(R,m): return {m-s:R[s] for s in range(len(R)) if R[s]}
def prod(A,B,sg):
    o={}
    for i,x in A.items():
        for j,y in B.items():
            k=(min(i,j),max(i,j)); o[k]=o.get(k,0)+sg*x*y
    return o
def addq(*qs):
    o={}
    for q in qs:
        for k,v in q.items(): o[k]=o.get(k,0)+v
    return {k:v for k,v in o.items() if v}
def Dform(R,k): return addq(prod(lin(R,k),lin(R,k),1),prod(lin(R,k-1),lin(R,k+1),-1))
def deltaform(R,k): return addq(Dform(R,k),{kk:-v for kk,v in Dform(R,k+1).items()})
def Etarget(q,sg):
    j=0; i=q+1; c=lambda m:{m:1}
    Bi=addq(c(i-1),c(i+1)); Bj=addq(c(j-1),c(j+1))
    W=addq(prod(Bi,c(j),1),prod(c(i),Bj,-1))
    return addq(Dform([1],j),{k:-v for k,v in Dform([1],i).items()},{k:-sg*v for k,v in W.items()})
q=int(sys.argv[1]); dmax=int(sys.argv[2]); tmax=int(sys.argv[3]) if len(sys.argv)>3 else q+dmax+2
for sg in (1,-1):
    cols=[]; names=[]
    for d in range(0,dmax+1):
        for al in range(d+1):
            be=d-al; R=Rpoly(al,be)
            for t in range(d//2+1, tmax+1):
                if 2*t<=d: continue
                cols.append(deltaform(R,t)); names.append((al,be,t))
    tg=Etarget(q,sg)
    keys=sorted(set(k for cq in cols for k in cq)|set(tg))
    A=np.array([[float(cq.get(k,0)) for cq in cols] for k in keys]); b=np.array([float(tg.get(k,0)) for k in keys])
    res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
    print('q=%d sign=%+d dmax=%d: %s'%(q,sg,dmax,'FEASIBLE' if res.status==0 else 'infeasible'), end='')
    if res.status==0:
        print('  ', [(n,round(w,4)) for w,n in zip(res.x,names) if w>1e-9])
    else: print()
