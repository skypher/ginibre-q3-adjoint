# LP search for Newton certificates of (E): D_j - D_i +/- W(p,q) >= 0 on reciprocal rows (c_{N-k} = c_k, e even),
# as an identity of quadratic forms in the free coefficients c_0..c_floor(N/2).
import sys
import numpy as np
from scipy.optimize import linprog
exec(open('window_lp.py').read().split("C=int(sys.argv[1])")[0])
def sym_index(k,N):
    if k<0 or k>N: return None
    return min(k,N-k)
def quad_sym(q,N):
    o={}
    for (i,j),v in q.items():
        a,b=sym_index(i,N),sym_index(j,N)
        if a is None or b is None: continue
        key=(min(a,b),max(a,b)); o[key]=o.get(key,0)+v
    return {k:v for k,v in o.items() if v!=0}
def E_target(N,p,q,sg):
    j=(N+p-q)//2; i=(N+p+q)//2+1; o={}
    def add(x,y,v):
        key=(min(x,y),max(x,y)); o[key]=o.get(key,0)+v
    for (k,s) in ((j,1),(i,-1)): add(k,k,s); add(k-1,k+1,-s)
    # W = (c_{i-1}+c_{i+1})c_j - c_i(c_{j-1}+c_{j+1})
    for (x,y,v) in ((i-1,j,1),(i+1,j,1),(i,j-1,-1),(i,j+1,-1)): add(x,y,sg*v)
    return quad_sym(o,N)
maxdeg=int(sys.argv[1]); ts=[Fr(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 and sys.argv[2] else []
from collections import Counter
st=Counter(); ex=[]
for N in range(4,13,2):
    cols=[]; names=[]
    for name,R in mults(maxdeg,ts):
        for k in range(0,N+len(R)):
            q=quad_sym(quad_D(R,k,0,N),N)
            if q: cols.append(q); names.append((name,k))
    for p in range(0,N+1):
        for q in range(0,p+1):
            if (N+p+q)%2 or p+q>N: continue
            for sg in (1,-1):
                tg=E_target(N,p,q,sg)
                if not tg: continue
                keys=sorted(set(k for cq in cols for k in cq)|set(tg))
                A=np.array([[float(cq.get(key,0)) for cq in cols] for key in keys]); b=np.array([float(tg.get(key,0)) for key in keys])
                res=linprog(np.zeros(len(cols)),A_eq=A,b_eq=b,bounds=[(0,None)]*len(cols),method='highs')
                st['feasible' if res.status==0 else 'infeasible']+=1
                if res.status!=0 and len(ex)<6: ex.append((N,p,q,sg))
print('maxdeg',maxdeg,'ts',ts,dict(st),'first infeasible (N,p,q,sign):',ex)
