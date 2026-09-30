# (i) |T(x,x+C)| <= (C+1) sqrt(D_x D_(x+C)) and (ii) sum D >= (C+1) sqrt(D_x D_(x+C)) on general real-rooted rows.
import math, random
from fractions import Fraction as Fr
from collections import Counter
exec(open('window_realrooted.py').read().split("random.seed(1)")[0])
random.seed(8); st=Counter(); ex={}
for fam in ['mixed real roots','positive coefficients','reciprocal real-rooted (e odd)']:
    for _ in range(500):
        if fam=='mixed real roots': rts=[Fr(random.choice([-1,1])*random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(3,14))]
        elif fam=='positive coefficients': rts=[Fr(random.randint(1,9),random.randint(1,9)) for _ in range(random.randint(3,14))]
        else:
            rts=[-1]*random.choice([1,3,5])+[1]*random.randint(0,5)
            for _ in range(random.randint(1,3)):
                t=Fr(random.randint(20,40),10)*random.choice([-1,1]); 
                # roots of 1+tz+z^2 are real: include as quadratic factor via two roots rho,1/rho -> use exact quadratic multiply below
                rts.append(('q',t))
        c=[Fr(1)]
        for rho in rts:
            if isinstance(rho,tuple):
                t=rho[1]; c=[(c[k] if k<len(c) else 0)+t*(c[k-1] if 1<=k<=len(c) else 0)+(c[k-2] if 2<=k<=len(c)+1 else 0) for k in range(len(c)+2)]
            else: c=[(c[k] if k<len(c) else 0)+rho*(c[k-1] if k>=1 else 0) for k in range(len(c)+1)]
        N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+1)]
        for x in range(N+1):
            for Cw in range(1,N-x+1):
                n=Cw+1; T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1); g2=D[x]*D[x+Cw]
                st[(fam,'windows')]+=1
                if T*T>n*n*g2: st[(fam,'(i) fails')]+=1; ex.setdefault((fam,'i'),(x,Cw,float(T*T/(n*n*g2))))
                S=sum(D[x:x+Cw+1])
                if S*S<n*n*g2: st[(fam,'(ii) fails')]+=1; ex.setdefault((fam,'ii'),(x,Cw,float(S*S/(n*n*g2))))
for k,v in sorted(st.items()): print(k,v)
print(ex)
