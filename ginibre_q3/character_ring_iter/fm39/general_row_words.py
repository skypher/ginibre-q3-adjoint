# Generalized FM3 test: replace the row of (1+z)^a(1-z)^e by any reciprocal real-rooted row
# c = (1+z)^a (1-z)^e prod(1+t z+z^2), |t|>=2, and evaluate the split identity with the quadratic W_c.
import random, time, sys
from fractions import Fraction as Fr
from itertools import product
exec(open('wformula_check.py').read().split("n=0")[0])
exec(open('recip_test.py').read().split("def tools")[0].split("import random")[1].replace("from fractions import Fraction as Fr\nfrom math import comb\n",""))
from collections import Counter
def cgmul(A,l):  # A: Counter of labels -> product with U_l
    out=Counter()
    for p,mlt in A.items():
        for s in range(abs(p-l),p+l+1,2): out[s]+=mlt
    return out
def phi_row(c,Ls,Ps,eps):
    tot=Fr(0); m=len(Ls)
    for mask in range(1<<m):
        for tmask in range(1<<len(Ps)):
            X=Counter({0:1}); Y=Counter({0:1})
            for i,l in enumerate(Ls): X,Y=(cgmul(X,l),Y) if mask>>i&1 else (X,cgmul(Y,l))
            for j,p in enumerate(Ps): X,Y=(cgmul(X,p),Y) if tmask>>j&1 else (X,cgmul(Y,p))
            sgn=(-1)**(m-bin(mask).count('1'))
            for p,mp in X.items():
                for q,mq in Y.items(): tot+=sgn*mp*mq*Wq(c,p,q,eps)
    return tot/2
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 21)
def randt():
    s=random.choice([-1,1]); return s*(2+Fr(random.randint(0,12),random.randint(1,6)))
stats=Counter(); fails=[]
t0=time.time(); last=t0
# control: base rows must reproduce the direct evaluator for a few words
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
for (hs,ss,r,a) in [((2,3,4),(),3,2),((2,2,3,5),(),3,1),((3,4),(2,),2,3),((2,2,2,2,2),(),4,1)]:
    m=len(hs); e=2*r-m; c=mkrow(a,e,[])
    assert phi_row(c,[k+1 for k in hs],list(ss),e%2)==phi_kernel(core(hs,ss),r,a),(hs,ss,r,a)
print('control: split evaluation with quadratic W reproduces fm3kern on base rows')
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 400):
    m=random.randint(2,6); Ls=[random.randint(3,7) for _ in range(m)]
    Ps=[random.randint(1,5) for _ in range(random.choice([0,0,1,2]))]
    e=random.choice([x for x in range(0,6) if x%2==m%2]); a=random.randint(0,4)
    ts=[randt() for _ in range(random.randint(1,3))]
    c=mkrow(a,e,ts); N=len(c)-1
    if (N+sum(Ls)+sum(Ps))%2: a+=1; c=mkrow(a,e,ts); N+=1
    v=phi_row(c,Ls,Ps,e%2); stats['words']+=1
    if v<0: stats['neg']+=1; fails.append((m,Ls,Ps,a,e,[str(t) for t in ts],float(v)))
    if v==0: stats['zero']+=1
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat',dict(stats),flush=True)
print(dict(stats),'elapsed',round(time.time()-t0,1)); print('first failures:',fails[:5])
