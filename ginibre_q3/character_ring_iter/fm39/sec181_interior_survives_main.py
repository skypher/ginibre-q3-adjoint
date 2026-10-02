import sys
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from str7e_budget_runs_main import table, interior_cut, mul_xx, mul_xy, toppair
from collections import defaultdict
def check(L):
    N=len(L); res=[]; cache={}
    def tab(X):
        X=tuple(sorted(X))
        if X not in cache: cache[X]=table(list(X))
        return cache[X]
    pb=bb=0; tot=0
    for i in range(N):
        for j in range(i+1,N):
            zu,zv=L[i],L[j]; a,b=abs(zu),abs(zv); eu=1 if zu>0 else -1; ev=1 if zv>0 else -1
            C=[z for q,z in enumerate(L) if q not in (i,j)]; A,Bb=interior_cut(C)
            rA=tab(A); FB=tab(Bb); Gp=mul_xx(FB,a,b,eu*ev); Gm=mul_xy(FB,a,b,eu,ev)
            P=defaultdict(int); M=defaultdict(int)
            for (r,s),v in rA.items(): P[r+s]+=v*Gp.get((r,s),0); M[r+s]+=v*Gm.get((r,s),0)
            bud=0; plain=0; okb=okp=True
            for t in sorted(set(P)|set(M)):
                bud+=P[t]+min(M[t],0); plain+=P[t]+M[t]; okb&=bud>=0; okp&=plain>=0
            tot+=1; pb+=not okp; bb+=not okb
    return tot,pb,bb
# FM-SEC181 main-agent check: on the three HPP counterexample lists, the interior cut of every pair has nonnegative plain prefixes.
for L in ([1,1,-2,3,3,3,3,3,3,-4,-4,-4,5,5,5,5,6],[1,1,-2,3,3,3,3,3,-4,-4,-4,5,5,5,5,5,6],[1,1,1,-2,2,3,3,3,3,3,-4,-4,-4,4,4,5,5,8]):
    tot,pb,bb=check(L)
    B=L[:-1]; i,j=toppair(B)
    print(L,"pairs",tot,"plain-prefix failures",pb,"budget failures",bb,"TopPair",(B[i],B[j]))
    assert pb==0
print('FM-SEC181 INTERIOR CUT SURVIVES PASS')
