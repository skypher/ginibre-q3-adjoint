# Is Conjecture W implied by Newton's inequalities alone?  Positive ultra-log-concave rows w_k = C(n,k) u_k, u log-concave,
# and plain log-concave rows (control).
import random
from fractions import Fraction as Fr
from math import comb
exec(open('window_f.py').read().split("random.seed(7)")[0])
def Pmin_row(c):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    T=lambda j,k: C_(j)*C_(k)-C_(j-1)*C_(k+1)
    worst=None
    for x in range(N+1):
        for Cw in range(1,N-x+1):
            P=sum(T(k,k) for k in range(x,x+Cw+1))-T(x,x+Cw)
            if worst is None or P<worst[0]: worst=(P,x,Cw)
    return worst
random.seed(5)
for name,ulc in [('ultra-log-concave (Newton only), not nec. real-rooted',True),('log-concave only (control)',False)]:
    neg=0; tot=0; ex=None; nonrr=0
    for _ in range(600):
        n=random.randint(3,14)
        # log-concave u: u_k = exp(-convex) via decreasing ratios
        ratios=sorted([Fr(random.randint(1,40),10) for _ in range(n)],reverse=True)
        u=[Fr(1)]
        for q in ratios: u.append(u[-1]*q)
        c=[comb(n,k)*u[k] for k in range(n+1)] if ulc else u
        w=Pmin_row(c); tot+=1
        if w[0]<0: neg+=1; ex=ex or (n,[str(v) for v in c],w)
    print(name,': rows',tot,'rows with negative window',neg,'example',ex)
print('--- mixed signs ---')
random.seed(6)
def newton_ok(c):
    n=len(c)-1
    return all(c[k]**2*k*(n-k) >= c[k-1]*c[k+1]*(k+1)*(n-k+1) for k in range(1,n))
for name in ['ULC magnitudes with random signs (Newton holds with signs)','ULC magnitudes, one sign change']:
    neg=tot=0; ex=None
    while tot<600:
        n=random.randint(3,14)
        ratios=sorted([Fr(random.randint(1,40),10) for _ in range(n)],reverse=True)
        u=[Fr(1)]
        for q in ratios: u.append(u[-1]*q)
        c=[comb(n,k)*u[k] for k in range(n+1)]
        if name.startswith('ULC magnitudes with random'): c=[x*random.choice([-1,1]) for x in c]
        else:
            j=random.randint(1,n); c=[x if k<j else -x for k,x in enumerate(c)]
        if not newton_ok(c): continue
        w=Pmin_row(c); tot+=1
        if w[0]<0: neg+=1; ex=ex or (n,w)
    print(name,': rows',tot,'rows with negative window',neg,ex)
