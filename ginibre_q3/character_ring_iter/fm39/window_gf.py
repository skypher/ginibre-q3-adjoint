# Is G_C(s) = sum_x P_C(x) s^x real-rooted for real-rooted P?  (numpy roots, tolerance on imaginary parts)
import numpy as np, random
from fractions import Fraction as Fr
exec(open('window_realrooted.py').read().split("random.seed(1)")[0])
def Pw(c,x,Cw):
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    T=lambda j,k: C_(j)*C_(k)-C_(j-1)*C_(k+1)
    return sum(T(k,k) for k in range(x,x+Cw+1))-T(x,x+Cw)
random.seed(9); res={}
for fam in ['(1+z)^a(1-z)^e','random real roots']:
    for _ in range(300):
        if fam.startswith('(1'):
            a=random.randint(0,12); e=random.randint(0,12); c=row([1]*a+[-1]*e)
        else:
            c=row([Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(3,12))])
        N=len(c)-1
        for Cw in range(1,min(N,6)+1):
            co=[float(Pw(c,x,Cw)) for x in range(-1,N+2)]
            while co and abs(co[-1])<1e-300: co.pop()
            while co and abs(co[0])<1e-300: co.pop(0)
            if len(co)<2: continue
            rts=np.roots(co[::-1]); scale=max(1,max(abs(rts)))
            rr=all(abs(z.imag)<=1e-7*max(1,abs(z)) for z in rts)
            key=(fam,'C=%d'%Cw); s=res.setdefault(key,[0,0]); s[0]+=1; s[1]+=rr
for k in sorted(res): print(k,'real-rooted G_C in',res[k][1],'of',res[k][0])
