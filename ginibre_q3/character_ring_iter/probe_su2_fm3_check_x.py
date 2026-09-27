#!/usr/bin/env python3
"""Verify in the 3-nesting-free model: m_20 = sum_{p<=q} #{G on k-e_p-e_q: alpha(G)<=p},
m_11 = X = sum_{p<q} #{G on k-e_p-e_q: alpha(G)<=p, beta(G)<=q}; compare with character values.
usage: probe_su2_fm3_check_x.py MAXSUM MAXPART"""
import sys, itertools
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
MS,MP=map(int,sys.argv[1:3])
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__),'probe_su2_fm3_nesting_model.py')).read().split("bad=0; tot=0")[0].replace('if len(sys.argv)<2','if False'))
def nested2_right_of(edges,p):
    E=[e for e,m in edges if e[0]>p]
    for e1 in E:
        for e2 in E:
            if e1[0]<e2[0] and e2[1]<e1[1]: return True
    return False
def alpha(edges,J):
    for p in range(-1,J):   # vertices 0..J-1 ; p=-1 means everything to the right
        if not nested2_right_of(edges,p): return p
    return J
def beta(edges):
    return max([e[0] for e,m in edges],default=-1)
def mults(ks):
    F={(0,0):1}
    for k in ks: F=cmul(F,H(k))
    Psi=cmul(F,D1)
    return Psi.get((1,0),0),Psi.get((2,1),0),Psi.get((3,0),0)
bad=0; tot=0
for n in range(1,5):
    for ks in itertools.product(range(0,MP+1),repeat=n):
        if sum(ks)>MS or sum(ks)%2: continue
        J=len(ks); m00,m11,m20=mults(list(ks))
        X=0; Y=0
        for p in range(J):
            for q in range(p,J):
                d=list(ks); d[p]-=1; d[q]-=1
                if min(d)<0: continue
                for g in multigraphs(d):
                    if has_nesting(g,3): continue
                    if alpha(g,J)<=p:
                        Y+=1
                        if p<q and beta(g)<=q: X+=1
        tot+=1
        if (X,Y)!=(m11,m20): bad+=1; print("MISMATCH",ks,(m11,m20),(X,Y))
print("checked",tot,"sequences; mismatches",bad)
