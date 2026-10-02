# FM-STR12b (main agent): which half-plane cutoffs keep sum_{lambda.(r,s) <= T} f_A f_B >= 0 for every split?
# Pair-free lists with labels <= 4, length 3..7, even minus count, every split, every T.
# Expected: no failure when lambda_1 + lambda_2 >= 0; failures for the listed directions with lambda_1 + lambda_2 < 0,
# and for non-half-plane cutoffs (box, max, band |r-s|, Casimir disk).
import itertools, sys, time
from collections import defaultdict
def cg(a,b): return range(abs(a-b),a+b+1,2)
def table(word):
    d={(0,0):1}
    for z in sorted(word,key=abs,reverse=True):
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n): q[c,b]+=v
            for c in cg(b,n): q[a,c]+=e*v
        d={k:v for k,v in q.items() if v}
    return d
GOOD=[(1,1),(1,2),(1,3),(2,3),(1,10),(10,1),(1,0),(0,1),(1,-1),(-1,1),(2,-1),(3,-2),(10,-9),(-9,10),(11,-10)]
BAD=[(-1,-1),(-1,0),(1,-2),(-1,-2)]
def first_bad(prod,lam):
    lay=defaultdict(int)
    for (r,s),v in prod.items(): lay[lam[0]*r+lam[1]*s]+=v
    run=0
    for t in sorted(lay):
        run+=lay[t]
        if run<0: return t
    return None
labs=[s*n for n in range(1,5) for s in (1,-1)]
fails=defaultdict(int); splits=0
for N in range(3,8):
    for w in itertools.combinations_with_replacement(sorted(labs),N):
        if sum(1 for z in w if z<0)%2 or sum(map(abs,w))%2 or any(-z in w for z in w): continue
        for m in range(1,1<<(N-1)):
            A=[w[i] for i in range(N) if m>>i&1]; B=[w[i] for i in range(N) if not m>>i&1]
            fa=table(A); fb=table(B); prod={k:v*fb.get(k,0) for k,v in fa.items() if fb.get(k,0)}
            if not prod: continue
            splits+=1
            for lam in GOOD+BAD:
                if first_bad(prod,lam) is not None: fails[lam]+=1
    print(time.strftime('%H:%M:%S'),'length',N,'splits',splits,flush=True)
print('failures by direction:',dict(fails))
assert all(fails[l]==0 for l in GOOD), 'a direction with lambda_1+lambda_2 >= 0 failed'
assert all(fails[l]>0 for l in BAD)
print('FM-STR12b HALF-PLANE DIRECTIONS PASS')
