# GFM3 stress screen on the hardest words: clustered small labels, many factors, minimal e, several hat S.
import random, time, sys
from fractions import Fraction as Fr
exec(open('general_row_words.py').read().split("random.seed(")[0])
random.seed(int(sys.argv[1])); count=int(sys.argv[2])
def randt():
    s=random.choice([-1,1]); return s*(2+Fr(random.randint(0,8),random.randint(1,8)))
from collections import Counter
stats=Counter(); fails=[]; t0=time.time(); last=t0
for trial in range(count):
    m=random.randint(3,8); base=random.randint(3,5)
    Ls=[base+random.randint(0,1) for _ in range(m)]            # clustered labels
    Ps=[random.randint(1,4) for _ in range(random.choice([0,1,2,3]))]
    e=m%2                                                      # minimal e (FM2 boundary m = 2r or 2r-1)
    a=random.randint(0,3); ts=[randt() for _ in range(random.randint(1,4))]
    c=mkrow(a,e,ts); N=len(c)-1
    if (N+sum(Ls)+sum(Ps))%2: a+=1; c=mkrow(a,e,ts)
    v=phi_row(c,Ls,Ps,e%2); stats['words']+=1
    if v<0: stats['neg']+=1; fails.append((Ls,Ps,a,e,[str(t) for t in ts],float(v)))
    if v==0: stats['zero']+=1
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat',trial+1,'/',count,dict(stats),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(stats),'elapsed',round(time.time()-t0,1)); print('failures:',fails[:8])
