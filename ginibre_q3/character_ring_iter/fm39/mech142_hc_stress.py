# Independent stress test of h_c >= 0 for Minerva's quartet identity at the boundary q = D + M + 1.
import sys
from collections import defaultdict
from random import Random
def cg(a,b): return range(abs(a-b),a+b+1,2)
def fusion(ns):
    A={0:1}
    for n in ns:
        T=defaultdict(int)
        for a,v in A.items():
            for c in cg(a,n): T[c]+=v
        A=T
    return A
PARTS=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
def hc(ns,D):
    row=fusion(ns); N=[0]*(D+1); Z=[0]*(D+1)
    for I,J in PARTS:
        FA=fusion(tuple(ns[i] for i in I)); FB=fusion(tuple(ns[j] for j in J))
        for a in FA:
            for b in FB:
                if a+b>D: continue
                if a==0 or b==0: Z[a+b]+=1
                else:
                    for c in cg(a,b):
                        if c<=D: N[c]+=1
    return [row.get(c,0)+Z[c]-N[c] for c in range(D%2,D+1,2)]
rng=Random(int(sys.argv[1])); worst=None; n=0
for _ in range(int(sys.argv[2])):
    D=rng.randrange(0,60); M=rng.randrange(0,16)
    t=[0]+sorted(rng.randrange(M+1) for _ in range(2))+[M]
    if (D+sum(t))%2: continue
    for q in (D+M+1, D+M+2):
        h=hc(tuple(q+x for x in t),D); n+=1
        m=min(h)
        if worst is None or m<worst[0]: worst=(m,D,M,t,q)
        assert m>=0,(D,M,t,q,h)
print("cases",n,"min h_c",worst[0],"at",worst[1:])
