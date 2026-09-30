# F_R(T) = sum_S eps^S #{(M1 NC on S-legs, M2 NC on S^c-legs): every mutual crossing's 4-set of legs lies in R}
import itertools, sys, random
from collections import defaultdict
def matchings(legs, owner):
    if not legs: yield []; return
    a=legs[0]
    for j in range(1,len(legs)):
        b=legs[j]
        if owner[a]==owner[b]: continue
        for m in matchings(legs[1:j]+legs[j+1:],owner): yield [(a,b)]+m
def cr(e,f):
    a,b=e; c,d=f
    return (a<c<b<d) or (c<a<d<b)
def nc(M): return all(not cr(e,f) for e,f in itertools.combinations(M,2))
def configs(labels,T):
    owner=[]
    for i,n in enumerate(labels): owner+=[i]*n
    L=len(labels); out=[]  # list of (sign, frozenset of crossing 4-sets)
    for S in range(1<<L):
        sign=(-1)**bin(S&T).count('1')
        l1=[k for k in range(len(owner)) if S>>owner[k]&1]; l2=[k for k in range(len(owner)) if not S>>owner[k]&1]
        if len(l1)%2: continue
        M1s=[M for M in matchings(l1,owner) if nc(M)]; M2s=[M for M in matchings(l2,owner) if nc(M)]
        for M1 in M1s:
            for M2 in M2s:
                X=frozenset(tuple(sorted(e+f)) for e in M1 for f in M2 if cr(e,f))
                out.append((sign,X))
    return out, len(owner)
rng=random.Random(1)
tot=0; neg=0; worst=None
for L in range(2,7):
    for lab in itertools.combinations_with_replacement(range(1,4),L):
        D=sum(lab)
        if D%2 or D>8: continue
        for T in range(1<<L):
            if bin(T).count('1')%2: continue
            cf,D=configs(lab,T)
            quads=sorted(set().union(*[X for s,X in cf])) if cf else []
            nq=len(quads)
            Rs = range(1<<nq) if nq<=14 else [rng.getrandbits(nq) for _ in range(20000)]
            for Rm in Rs:
                R={quads[k] for k in range(nq) if Rm>>k&1}
                v=sum(s for s,X in cf if X<=R)
                tot+=1
                if v<0:
                    neg+=1
                    if worst is None or v<worst[0]: worst=(v,lab,T,sorted(R))
print("R-profiles",tot,"negative",neg,"worst",worst)
