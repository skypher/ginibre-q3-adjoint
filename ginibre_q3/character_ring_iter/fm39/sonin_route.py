# Route test on base rows: (i) |T(x,x+C)| <= n sqrt(D_x D_(x+C)), n = C+1;  (ii) sum D >= n sqrt(D_x D_(x+C));
# also log-concavity of D_k.
exec(open('split3.py').read())
import math
from collections import Counter
st=Counter(); ex={}
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(0,N+1)]
        lc=all(D[k]**2>=D[k-1]*D[k+1] for k in range(1,N)); st['rows']+=1; st['D log-concave rows']+=lc
        if not lc: ex.setdefault('nonLC',(r,a))
        for x in range(0,N+1):
            for Cw in range(1,N-x+1):
                n=Cw+1; T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1); g=math.sqrt(D[x]*D[x+Cw])
                st['windows']+=1
                if abs(T)>n*g*(1+1e-12): st['(i) fails']+=1; ex.setdefault('i',(r,a,x,Cw,T/(n*g)))
                if sum(D[x:x+Cw+1])<n*g*(1-1e-12): st['(ii) fails']+=1; ex.setdefault('ii',(r,a,x,Cw))
print(dict(st)); print(ex)
