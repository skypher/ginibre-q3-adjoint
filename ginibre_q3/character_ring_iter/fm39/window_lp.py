# LP search: P_C(0) = sum_i w_i D_(k_i)(R_i c), w_i >= 0, R_i real-rooted multipliers, as an identity of quadratic forms in c.
# c is a free sequence on indices -L..C+1+L (window at x=0 uses c_(-1)..c_(C+1)); R*c shifts indices up by deg R (we allow index offset).
import itertools, sys
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from math import comb
def polymul(a,b):
    o=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): o[i+j]+=x*y
    return o
def mults(maxdeg,ts):
    out=[]
    for j in range(maxdeg+1):
        for l in range(maxdeg+1-j):
            R=[1]
            for _ in range(j): R=polymul(R,[1,-1])
            for _ in range(l): R=polymul(R,[1,1])
            out.append((('(1-z)^%d(1+z)^%d'%(j,l)),R))
    base=list(out)
    for t in ts:
        for nm,R in base:
            if len(R)-1<maxdeg: out.append(('(1+%sz)'%t+nm,polymul(R,[1,t])))
    return out
def quad_D(R,k,lo,hi):
    # D_k(R c) as dict over pairs (i<=j) of c indices; (Rc)_m = sum_s R_s c_(m-s)
    def lin(m): return {m-s:R[s] for s in range(len(R)) if R[s]!=0}
    def prod(A,B,sg):
        o={}
        for i,x in A.items():
            for j,y in B.items():
                key=(min(i,j),max(i,j)); o[key]=o.get(key,0)+sg*x*y
        return o
    q={}
    for d in (prod(lin(k),lin(k),1),prod(lin(k-1),lin(k+1),-1)):
        for kk,v in d.items(): q[kk]=q.get(kk,0)+v
    return {kk:v for kk,v in q.items() if v!=0}
def target(C):
    # P_C(0) = sum_{k=0}^{C} (c_k^2 - c_(k-1)c_(k+1)) - c_0 c_C + c_(-1) c_(C+1)
    q={}
    def add(i,j,v):
        key=(min(i,j),max(i,j)); q[key]=q.get(key,0)+v
    for k in range(C+1): add(k,k,1); add(k-1,k+1,-1)
    add(0,C,-1); add(-1,C+1,1)
    return {kk:v for kk,v in q.items() if v!=0}
C=int(sys.argv[1]); maxdeg=int(sys.argv[2]); ts=[Fr(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 and sys.argv[3] else []
lo,hi=-1,C+1
EXT=int(sys.argv[4]) if len(sys.argv)>4 else 0
cols=[]; names=[]
for name,R in mults(maxdeg,ts):
    for k in range(lo-EXT, hi+EXT+len(R)):
        q=quad_D(R,k,lo,hi)
        if all(lo-EXT<=i and j<=hi+EXT for (i,j) in q) and q: cols.append(q); names.append((name,k))
keys=sorted(set(k for q in cols for k in q)|set(target(C)))
A=np.array([[float(q.get(key,0)) for q in cols] for key in keys]); b=np.array([float(target(C).get(key,0)) for key in keys])
res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
print('C=%d maxdeg=%d ts=%s: columns %d, status %s'%(C,maxdeg,ts,len(cols),res.status), 'FEASIBLE' if res.status==0 else 'infeasible')
if res.status==0:
    for w,nm in zip(res.x,names):
        if w>1e-9: print('   %.6g * D_%d[%s c]'%(w,nm[1],nm[0]))
