# F_{q,s}(T) = sum_S eps^S sum_{M1 on S-legs, M2 on S^c-legs} q^{cr(M1)+cr(M2)} s^{cr(M1,M2)}
# legs in linear order, block i legs consecutive; no chord inside a block.
import itertools, sys
from collections import defaultdict
from functools import lru_cache
def matchings(legs, owner):
    if not legs: yield []; return
    a=legs[0]
    for j in range(1,len(legs)):
        b=legs[j]
        if owner[a]==owner[b]: continue
        rest=legs[1:j]+legs[j+1:]
        for m in matchings(rest,owner): yield [(a,b)]+m
def cr(e,f):
    a,b=e; c,d=f
    return (a<c<b<d) or (c<a<d<b)
def Fqs(labels,T):
    owner=[]; 
    for i,n in enumerate(labels): owner+= [i]*n
    L=len(labels); poly=defaultdict(int)
    for S in range(1<<L):
        sign=(-1)**bin(S&T).count('1')
        l1=[k for k in range(len(owner)) if S>>owner[k]&1]
        l2=[k for k in range(len(owner)) if not S>>owner[k]&1]
        if len(l1)%2: continue
        M1s=list(matchings(l1,owner)); M2s=list(matchings(l2,owner))
        if not M1s or not M2s: continue
        for M1 in M1s:
            c1=sum(cr(e,f) for e,f in itertools.combinations(M1,2))
            for M2 in M2s:
                c2=sum(cr(e,f) for e,f in itertools.combinations(M2,2))
                c12=sum(cr(e,f) for e in M1 for f in M2)
                poly[(c1+c2,c12)]+=sign
    return {k:v for k,v in poly.items() if v}
def ev(p,q,s): return sum(v*q**i*s**j for (i,j),v in p.items())
import random
random.seed(3)
lists=[]
for L in range(2,7):
    for lab in itertools.combinations_with_replacement(range(1,4),L):
        if sum(lab)%2==0 and sum(lab)<=12: lists.append(lab)
grid=[i/20 for i in range(-20,21)]
viol_diag=0; minsq=(1e9,None); negcoef_s1=0; cnt=0; neg_region=defaultdict(int)
for lab in lists:
    L=len(lab)
    for T in range(1<<L):
        if bin(T).count('1')%2: continue
        p=Fqs(lab,T); cnt+=1
        # region |s|<=q
        for q in grid:
            for s in grid:
                v=ev(p,q,s)
                if v< -1e-9:
                    neg_region[(q,s)]+=1
                    if abs(s)<=q+1e-12: viol_diag+=1
                if 0<=q and 0<=s and v<minsq[0]: minsq=(v,(lab,T,q,s))
        # s=1 line: q-coefficients
        qc=defaultdict(int)
        for (i,j),v in p.items(): qc[i]+=v
        if any(v<0 for v in qc.values()): negcoef_s1+=1
print("profiles",cnt,"violations in |s|<=q:",viol_diag,"neg q-coeff on s=1:",negcoef_s1)
print("min over [0,1]^2 grid:",minsq)
pts=sorted(neg_region.items())
print("grid points with some negative profile:",len(pts))
print(pts[:40])
print("negative points with q>=0:",[k for k,v in pts if k[0]>=0])
print("max q with a negative:",max(k[0] for k,v in pts))
