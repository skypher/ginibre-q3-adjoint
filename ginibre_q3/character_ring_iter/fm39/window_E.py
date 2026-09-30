# gamma <= N branch: test (a) P_C(l) >= P_C(l-A-1) and (b) P_C(l) - P_C(l-A-1) >= |(g_l + g_(tau-A-1)) ^ (g_tau + g_(l-A-1))|.
exec(open('split3.py').read())
from collections import Counter
st=Counter(); ex={}
for r in range(2,7):
    e=2*r-3
    for a in range(0,30):
        N=a+e; c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
        g=lambda k:(C_(k),C_(k-1)); wedge=lambda p,q:p[0]*q[1]-p[1]*q[0]; add=lambda p,q:(p[0]+q[0],p[1]+q[1])
        def P(x,C): return sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+C+1))-C_(x)*C_(x+C)+C_(x-1)*C_(x+C+1)
        for C in range(3,N+2):
            for B in range(C,N+3):
                for A in range(B,B+C+1):
                    if (N+A-B-C)%2: continue
                    l=(N+A-B-C)//2
                    if l+B>N or l<0: continue
                    tau=l+B+1; tp=l-A-1
                    d=P(l,C)-P(tp,C); X=wedge(add(g(l),g(tau-A-1)),add(g(tau),g(l-A-1)))
                    st['words']+=1
                    if d<0: st['(a) fails']+=1; ex.setdefault('a',(r,a,A,B,C,d))
                    if d<abs(X): st['(b) fails']+=1; ex.setdefault('b',(r,a,A,B,C,d,X))
                    if d-X<0: st['phi<0 ?!']+=1
print(dict(st)); print(ex)
