# Local route to MP_2 (H-only level 2, any number of factors):  phi_2(h_kappa) = sum_lambda K(lambda,kappa) w(lambda).
# Negative shapes (w = -2): lambda = (k+b+2, k+b+1, k+1, k).  Their positive neighbours:
#   lambda + e1 - e3 (w = 1),  lambda + e2 - e4 (w = 1, needs k >= 1),  lambda + e2 - e3 (w = 2);
# each positive shape of weight 1 or 2 is such a neighbour of exactly one negative shape (checked below).
# LOCAL inequality:  2 K(lambda,kappa) <= K(lambda+e1-e3,kappa) + K(lambda+e2-e4,kappa) + 2 K(lambda+e2-e3,kappa).
# If it holds for all negative lambda and all kappa, MP_2 follows (w = 5 shapes are extra slack).
import sys
from functools import lru_cache
from collections import Counter
exec(open('mp2_kostka_check.py').read().split("bad=0; n=0; ratios=set()")[0].split("src=open('mech28_eval.py').read()")[0])
exec(open('mp2_kostka_check.py').read().split("def w(lam):")[0].split("from itertools import combinations_with_replacement as cwr")[1])
def w(lam):
    l=list(lam)+[0]*(4-len(lam)); a,b,c=l[0]-l[1],l[1]-l[2],l[2]-l[3]
    if a==b==c==0: return 5
    if a==0 and c==0 and b>=1: return 2
    if (a,c) in ((2,0),(0,2)): return 1
    if a==1 and c==1: return -2
    return 0
def norm(l): return tuple(x for x in l if x>0)
def valid(l): return all(l[i]>=l[i+1] for i in range(3)) and l[3]>=0
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 14
st=Counter(); ex=[]; owner=Counter()
for size in range(4,NMAX+1):
    negs=[]
    for k in range(0,size//4+1):
        for b in range(0,size):
            if 4*k+2*b+4==size: negs.append((k+b+2,k+b+1,k+1,k))
    # partner check: each positive weight-1/2 shape used by one negative
    for lam in negs:
        for (i,j) in ((0,2),(1,3),(1,2)):
            mu=list(lam); mu[i]+=1; mu[j]-=1
            if valid(mu): owner[tuple(mu)]+=1
    for kappa in partitions(size, maxlen=size):
        for lam in negs:
            K0=kostka(norm(lam),kappa)
            if K0==0: continue
            nb=0
            for (i,j),wt in (((0,2),1),((1,3),1),((1,2),2)):
                mu=list(lam); mu[i]+=1; mu[j]-=1
                if valid(mu): nb+=wt*kostka(norm(mu),kappa)
            st['(lambda,kappa) tested']+=1
            if 2*K0<=nb: st['local holds']+=1
            else:
                st['local FAILS']+=1
                if len(ex)<10: ex.append((lam,kappa,K0,nb))
print(dict(st)); print('max times a positive shape is used as a neighbour:',max(owner.values()) if owner else None)
print('failures (lambda,kappa,K,neighbour sum):',ex)
