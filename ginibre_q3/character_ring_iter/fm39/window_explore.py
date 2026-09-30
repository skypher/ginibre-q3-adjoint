# Explore the one-row window inequality P_C(x) = sum_{k=x}^{x+C} T(k,k) - T(x,x+C) >= 0,
# T(j,k)=c_j c_k - c_{j-1} c_{k+1}, c = row of (1+z)^a (1-z)^e; also shifted family and the C-increment.
from math import comb
def cv(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return c
res={}
for a in range(0,31):
    for e in range(0,31):
        c=cv(a,e); N=a+e; C_=lambda k: c[k] if 0<=k<=N else 0
        T=lambda j,k: C_(j)*C_(k)-C_(j-1)*C_(k+1)
        for x in range(0,N+1):
            for Cw in range(1,N-x+1):
                P=sum(T(k,k) for k in range(x,x+Cw+1))-T(x,x+Cw)
                key=('P','odd e' if e%2 else 'even e'); res.setdefault(key,[0,0]); res[key][0]+=1; res[key][1]+= P<0
                inc=T(x+Cw,x+Cw) - (C_(x)*(C_(x+Cw)-C_(x+Cw-1)) - C_(x-1)*(C_(x+Cw+1)-C_(x+Cw)))
                key=('C-increment',); res.setdefault(key,[0,0]); res[key][0]+=1; res[key][1]+= inc<0
                for s in range(1,4):
                    Ps=sum(T(k,k+s) for k in range(x,x+Cw+1))-T(x,x+Cw+s)
                    key=('shift',s); res.setdefault(key,[0,0]); res[key][0]+=1; res[key][1]+= Ps<0
for k,v in sorted(res.items(),key=str): print(k,'tested',v[0],'negative',v[1])
