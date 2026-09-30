# S(w_0..w_M) = sum_{j=1}^{M-1} w_j^2 + sum_{j=1}^{M-2} w_j w_{j+1} - w_0 sum_{s=2}^{M} w_s - w_M sum_{s=1}^{M-2} w_s
# Claim tested: S>=0 on every window of consecutive coefficients (zero-extended) of every real-rooted polynomial.
import random
from fractions import Fraction as Fr
exec(open('window_realrooted.py').read().split("random.seed(1)")[0])
def S(w):
    M=len(w)-1
    return sum(w[j]**2 for j in range(1,M))+sum(w[j]*w[j+1] for j in range(1,M-1))-w[0]*sum(w[2:M+1])-w[M]*sum(w[1:M-1])
def windows(c,pad=2):
    z=[Fr(0)]*pad+list(c)+[Fr(0)]*pad
    for lo in range(len(z)):
        for M in range(2,len(z)-lo):
            yield z[lo:lo+M+1]
random.seed(7)
# sanity: S on window f_x..f_{x+C+1} of f=(1-z)P equals P_C(x) of P
def Pwin(c,x,Cw):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    T=lambda j,k: C_(j)*C_(k)-C_(j-1)*C_(k+1)
    return sum(T(k,k) for k in range(x,x+Cw+1))-T(x,x+Cw)
for _ in range(50):
    c=row([Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(2,9))])
    f=[ (c[k] if k<len(c) else 0)-(c[k-1] if k>=1 else 0) for k in range(len(c)+1)]
    F_=lambda k: f[k] if 0<=k<len(f) else 0
    for x in range(len(c)):
        for Cw in range(1,len(c)-x):
            assert S([F_(k) for k in range(x,x+Cw+2)])==Pwin(c,x,Cw)
print('f-form identity S(f window)=P_C(x) confirmed')
for name,gen in [('arbitrary real roots',lambda: [Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(1,12))]),
                 ('positive coefficients (neg roots)',lambda: [Fr(random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(1,12))]),
                 ('(1+z)^a(1-z)^e',lambda: [1]*random.randint(0,12)+[-1]*random.randint(0,12))]:
    neg=tot=0; ex=None
    for _ in range(400):
        c=row(gen())
        for w in windows(c):
            tot+=1; v=S(w)
            if v<0: neg+=1; ex=ex or ([str(t) for t in w],str(v))
    print(name,': windows',tot,'negative',neg,ex)
