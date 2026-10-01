import argparse, ast
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import linprog
import sympy as S

argparse.ArgumentParser(description="ADV-1b exact Gram and reindex tests").parse_args()
path=Path("ginibre_q3/character_ring_iter/fm39/mech54_genfun_obstructions_repro.py")
nodes=[n for n in ast.parse(path.read_text()).body
       if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
lib={}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),"exec"),lib)
fp,kr,profile=lib["fpoly"],lib["craw"],lib["profile"]

def square(c,d):
    ans=[Q(0)]*(d+1)
    for i,a in enumerate(c):
        for j,b in enumerate(c):
            if (i+j)%2==0: ans[(i+j)//2]+=a*b
    return ans

def certify(T,b,d):
    target=fp(T,b,d)
    basis=kr(T-2*b,d)
    columns=[square(c,d) for c in basis]
    for j in range(d-1):
        for exponent in range(-10,11):
            for sign in (-1,1):
                q=sign*Q(2)**exponent
                c=[Q(0)]*(j+3)
                for h,v in enumerate(basis[j]): c[h]+=v
                for h,v in enumerate(basis[j+2]): c[h]+=q*v
                columns.append(square(c,d))
    A=np.array(columns,dtype=float).T
    v=np.array(target,dtype=float)
    rs=np.where(abs(v)>0,abs(v),np.max(abs(A),axis=1))
    cs=np.max(abs(A/rs[:,None]),axis=0)
    cs[cs==0]=1
    result=linprog(np.ones(len(columns)),A_eq=A/rs[:,None]/cs,
                   b_eq=v/rs,bounds=(0,None),method="highs",
                   options={"primal_feasibility_tolerance":1e-9,
                            "dual_feasibility_tolerance":1e-9})
    assert result.success
    ids=[j for j,x in enumerate(result.x) if x>1e-10]
    exact=S.Matrix([[S.Rational(columns[j][i]) for j in ids]
                    for i in range(d+1)])
    rhs=S.Matrix([S.Rational(x) for x in target])
    weights,parameters=exact.gauss_jordan_solve(rhs)
    assert parameters.rows==0
    assert all(x>=0 for x in weights)
    assert exact*weights==rhs
    return len(ids)

cases=[(T,b,d) for T in (13,14) for b in range(3,10)
       for d in (6,7) if T+2*b-2*d>=7]
for case in cases: certify(*case)
print("Requested rows: 26 exact band-two certificates")
for case in [(40,3,6),(50,3,7),(100,3,7),
             (48,18,10),(36,14,14),(77,9,15)]:
    print("Additional certificate:",case,certify(*case))

for e,a,b,d in ((3,11,3,6),(6,8,3,6),(3,11,4,7),(6,8,4,7)):
    p=e+a+2*b-2*d
    G=profile(e,a,b)[p]
    H=profile(e-1,a+1,b)[p]
    parent=profile(e-1,a,b)
    A=parent[p-1]+(parent[p+1] if p+1<len(parent) else 0)
    assert G+H==2*A
    print("Profile, F, paired F:",(e,a,b,d),G,H)

xis=[1,3,5,7,9]
dual=[Q(1),Q(0),Q(-4,21),Q(-1,30),Q(-1,2310)]
for TT in range(9,40,2):
    for p in range(TT%2,TT+1,2):
        column=[profile((TT-x)//2,(TT+x)//2,0)[p] for x in xis]
        assert sum(y*z for y,z in zip(dual,column))>=0
target=[profile((13-x)//2,(13+x)//2,3)[7] for x in xis]
assert sum(y*z for y,z in zip(dual,target))==Q(-346112,1155)
print("Output-reindexed dictionary separator:",Q(-346112,1155))
