# Mirror form of the gamma <= N branch (anti-reciprocal rows): P_C(tau) = P_C(l-A-1) and
# phi = P_C(l) - P_C(l-A-1) - (g_l + g_(tau-A-1)) ^ (g_tau + g_(l-A-1)),  l = alpha, tau = gamma+1.
exec(open('split3.py').read())
import random
random.seed(2); n=0; ok1=ok2=0
for _ in range(3000):
    r=random.randint(2,6); a=random.randint(0,20); e=2*r-3; N=a+e
    w=random.randint(2,8); v=random.randint(w,12); u=random.randint(v,v+w)
    if (a+u+v+w)%2: continue
    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2; ga=al+B
    if not(al>=0 and ga<=N): continue
    c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
    g=lambda k:(C_(k),C_(k-1)); wedge=lambda p,q:p[0]*q[1]-p[1]*q[0]
    def P(x): return sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+C+1))-C_(x)*C_(x+C)+C_(x-1)*C_(x+C+1)
    l=al; tau=ga+1; taup=N-tau-C
    val=phi3(r,a,u,v,w); n+=1
    ok1+= (P(tau)==P(taup)) and taup==l-A-1
    add=lambda p,q:(p[0]+q[0],p[1]+q[1])
    ok2+= P(l)-P(taup)-wedge(add(g(l),g(tau-A-1)),add(g(tau),g(l-A-1)))==val
print('gamma<=N cases',n,' mirror window identity:',ok1,' mirror formula == phi:',ok2)
