# q-strengthening test for the H-only level-2 sector: replace Kostka numbers by Kostka-Foulkes polynomials
#   K_(lambda,mu)(q) = sum over SSYT T of shape lambda, content mu (mu a partition) of q^charge(T),
# and test whether P_mu(q) = sum_lambda K_(lambda,mu)(q) w(lambda) (MP_2 weights) has nonnegative coefficients.
# charge via Lascoux-Schutzenberger on the reading word (standard subword decomposition), for partition content.
import sys
from itertools import combinations
from collections import defaultdict
def partitions(n, maxlen=None, maxpart=None):
    if maxpart is None: maxpart=n
    if maxlen is None: maxlen=n
    if n==0: yield (); return
    if maxlen==0: return
    for p in range(min(n,maxpart),0,-1):
        for rest in partitions(n-p,maxlen-1,p): yield (p,)+rest
def ssyt(shape, content):
    # generate SSYT row by row via horizontal strips: fill numbers 1..len(content)
    res=[]
    def rec(i, cur_shape, rows):
        if i==len(content):
            if tuple(cur_shape)==tuple(shape): res.append([r[:] for r in rows])
            return
        k=content[i]
        # add a horizontal strip of size k to cur_shape staying inside shape
        L=len(shape)
        cs=list(cur_shape)+[0]*(L-len(cur_shape))
        def strips(r, left, new):
            if r==L:
                if left==0:
                    newrows=[row[:] for row in rows]+[[] for _ in range(L-len(rows))]
                    for rr in range(L):
                        newrows[rr]=newrows[rr]+[i+1]*(new[rr]-cs[rr])
                    rec(i+1, new, newrows)
                return
            upper=shape[r] if r==0 else min(shape[r], cs[r-1])   # horizontal strip: new[r] <= old[r-1]
            for add in range(0, min(left, upper-cs[r])+1):
                strips(r+1, left-add, new+[cs[r]+add])
        strips(0, k, [])
    rec(0, [0]*len(shape), [[] for _ in shape])
    return res
def reading_word(T):
    w=[]
    for row in reversed(T): w.extend(row)   # rows bottom to top, left to right
    return w
def charge_word(w):
    # Lascoux-Schutzenberger charge for words with partition content: extract standard subwords
    w=list(w); total=0
    idx=list(range(len(w)))
    while w:
        n=max(w)
        # extract standard subword: scan cyclically from the right for 1, then 2 to its left cyclically, ...
        pos=[]; used=set(); cur=len(w)
        c=0; charges=[]
        for letter in range(1,n+1):
            # find letter scanning leftward cyclically from cur-1
            L=len(w); found=None
            for step in range(1,L+1):
                j=(cur-step)%L
                if j in used: continue
                if w[j]==letter: found=j; break
            if found is None: break
            if letter>1 and found>cur: c+=1          # wrapped around: index increments
            charges.append(c); used.add(found); cur=found
        total+=sum(charges)
        w=[w[j] for j in range(len(w)) if j not in used]
    return total
def w_weight(lam):
    l=list(lam)+[0]*(4-len(lam)); a,b,c=l[0]-l[1],l[1]-l[2],l[2]-l[3]
    if a==b==c==0: return 5
    if a==0 and c==0 and b>=1: return 2
    if (a,c) in ((2,0),(0,2)): return 1
    if a==1 and c==1: return -2
    return 0
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 10
bad=[]; n=0
for size in range(2,NMAX+1):
    for mu in partitions(size):
        P=defaultdict(int)
        for lam in partitions(size, maxlen=4):
            wt=w_weight(lam)
            if wt==0: continue
            for T in ssyt(lam, mu):
                P[charge_word(reading_word(T))]+=wt
        n+=1
        neg={k:v for k,v in P.items() if v<0}
        if neg and len(bad)<10: bad.append((mu,dict(sorted(P.items()))))
        if neg: pass
print('partitions mu tested:',n,'with a negative q-coefficient:',len(bad)); 
for b in bad[:6]: print(b)
