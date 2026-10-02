# FM-STR12d (main agent): every distinct split of the pair-free list F1 = (+1)^2(-2)(+3)^4(-4)^3(+5)^4(+8).
# Height prefixes (r+s <= T) pass on all 599 splits; other directions fail (pair-free counterexamples to FM-STR12b).
import itertools
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
w=[1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8]
cls=sorted(set(w)); mult=[w.count(c) for c in cls]
dirs=[(1,1),(1,0),(1,-1),(1,2),(1,3),(2,3),(1,10)]
bad=defaultdict(int); first={}; n=0; cache={}
def tab(X):
    X=tuple(sorted(X))
    if X not in cache: cache[X]=table(list(X))
    return cache[X]
for take in itertools.product(*[range(m+1) for m in mult]):
    A=[c for c,t in zip(cls,take) for _ in range(t)]; B=[c for c,m,t in zip(cls,mult,take) for _ in range(m-t)]
    if not A or not B or tuple(sorted(A))>tuple(sorted(B)): continue
    fa=tab(A); fb=tab(B); P={k:v*fb[k] for k,v in fa.items() if k in fb}; n+=1
    assert sum(P.values())==8150742
    for d in dirs:
        lay=defaultdict(int)
        for (r,s),v in P.items(): lay[d[0]*r+d[1]*s]+=v
        run=0
        for t in sorted(lay):
            run+=lay[t]
            if run<0:
                bad[d]+=1; first.setdefault(d,(tuple(A),tuple(B),t,run)); break
print("F1 distinct splits",n,"failures by direction",dict(bad))
for d,v in first.items(): print("first failure",d,v)
assert n==599 and bad[(1,1)]==0
assert bad[(1,0)]==37 and bad[(1,-1)]==13 and bad[(1,2)]==5
print("FM-STR12d PASS: height survives every F1 split; (1,0), (1,-1), (1,2), (1,3), (2,3), (1,10) fail")
