# Check the note's MP_2 weights: phi_2(prod h_(kappa_i)) (no bound on the number of factors) versus
# sum_lambda K(lambda, kappa) w(lambda), l(lambda) <= 4, w by part differences (a,b,c) = (l1-l2, l2-l3, l3-l4):
#   w = 5 (a=b=c=0), 2 (a=c=0, b>=1), 1 ((a,c) in {(2,0),(0,2)}), -2 (a=c=1), 0 otherwise.
src=open('mech28_eval.py').read()
exec(src)
from functools import lru_cache
from itertools import combinations_with_replacement as cwr
def partitions(n, maxlen=4, maxpart=None):
    if maxpart is None: maxpart=n
    if n==0: yield (); return
    if maxlen==0: return
    for p in range(min(n,maxpart),0,-1):
        for rest in partitions(n-p,maxlen-1,p): yield (p,)+rest
@lru_cache(None)
def kostka(lam, kappa):
    # number of SSYT of shape lam and content kappa: strip off the largest entry as a horizontal strip
    lam=tuple(x for x in lam if x>0)
    if not kappa: return 1 if not lam else 0
    k=kappa[-1]; rest=kappa[:-1]
    if sum(lam)!=sum(kappa): return 0
    tot=0
    def strips(i, mu, left):
        nonlocal tot
        if i==len(lam):
            if left==0: tot_add(tuple(mu))
            return
        lo=lam[i+1] if i+1<len(lam) else 0
        for m in range(lam[i], lo-1, -1):
            if lam[i]-m>left: continue
            strips(i+1, mu+[m], left-(lam[i]-m))
    acc=[]
    def tot_add(mu): acc.append(mu)
    strips(0,[],k)
    return sum(kostka(mu,rest) for mu in acc)
def w(lam):
    l=list(lam)+[0]*(4-len(lam)); a,b,c=l[0]-l[1],l[1]-l[2],l[2]-l[3]
    if a==b==c==0: return 5
    if a==0 and c==0 and b>=1: return 2
    if (a,c) in ((2,0),(0,2)): return 1
    if a==1 and c==1: return -2
    return 0
bad=0; n=0; ratios=set()
for size in range(2,11):
    for kappa in partitions(size, maxlen=size):
        kap=tuple(kappa)
        lhs=phi(2, word(kap,()))
        rhs=sum(kostka(lam,kap)*w(lam) for lam in partitions(size))
        n+=1
        if rhs==0 and lhs==0: continue
        ratios.add(lhs/rhs if rhs else None)
        if (lhs<0)!=(rhs<0): bad+=1
print('kappa tested',n,'sign mismatches',bad,'ratios phi_2 / sum K w:',sorted(r for r in ratios if r is not None)[:5], len(ratios))
