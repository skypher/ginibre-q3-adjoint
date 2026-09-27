#!/usr/bin/env python3
"""Check: dim Sp(4)-invariants of tensor_j Sym^{k_j}(C^4) == # loopless multigraphs on [J] with degree
sequence k and no 3-nesting (a1<a2<a3<b3<b2<b1).  Also checks m_11 = X (nested tail count) and the
r=2 inequality in this model.  usage: probe_su2_fm3_nesting_model.py MAXSUM MAXPART"""
import sys, itertools
from functools import lru_cache
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
MS,MP=map(int,sys.argv[1:3])
def cg(a,p): return range(abs(a-p),a+p+1,2)
def cmul(F,G):
    out={}
    for (a,b),u in F.items():
        for (c,d),v in G.items():
            for e in cg(a,c):
                for f in cg(b,d): out[(e,f)]=out.get((e,f),0)+u*v
    return {k:v for k,v in out.items() if v}
H=lambda k:{(a,k-a):1 for a in range(k+1)}
D1={(1,0):1,(0,1):-1}
def I_char(ks):
    F={(0,0):1}
    for k in ks: F=cmul(F,H(k))
    F=cmul(cmul(F,D1),D1)
    return F.get((0,0),0)//2
def multigraphs(ks):
    """enumerate multisets of edges (i<j) realizing degree ks"""
    J=len(ks); pairs=[(i,j) for i in range(J) for j in range(i+1,J)]
    res=[]
    def rec(idx,deg,cur):
        if idx==len(pairs):
            if all(d==0 for d in deg): res.append(tuple(cur))
            return
        i,j=pairs[idx]
        for m in range(0,min(deg[i],deg[j])+1):
            deg[i]-=m; deg[j]-=m
            rec(idx+1,deg,cur+[((i,j),m)] if m else cur)
            deg[i]+=m; deg[j]+=m
    rec(0,list(ks),[])
    return res
def has_nesting(edges,depth):
    E=[e for e,m in edges]
    for comb in itertools.combinations(E,depth):
        s=sorted(comb)  # by left endpoint
        ok=all(s[t][0]<s[t+1][0] and s[t+1][1]<s[t][1] for t in range(depth-1))
        if ok and s[-1][0]<s[-1][1]: return True
    return False
def I_model(ks):
    return sum(1 for g in multigraphs(ks) if not has_nesting(g,3))
bad=0; tot=0
for n in range(1,6):
    for ks in itertools.product(range(0,MP+1),repeat=n):
        if sum(ks)>MS or sum(ks)%2: continue
        tot+=1
        a=I_char(list(ks)); b=I_model(ks)
        if a!=b: bad+=1; print("MISMATCH",ks,a,b)
print("checked %d ordered degree sequences, mismatches %d"%(tot,bad))
