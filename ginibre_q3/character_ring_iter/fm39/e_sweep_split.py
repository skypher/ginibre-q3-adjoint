# Two-case split of (E):  sweep(psi, j..i) < pi  -> convex-polygon argument from EQ1;
#                         sweep >= pi            -> does a metric bound (M_i), (M_j) or the general t-bound hold?
# The general bound 4t(4V-t) W^2 <= F_j(t) F_i(t), F_k(t) = 4t D_k + (X_k^2 - t) A_k^2, is tested at
# t in {X_i^2, X_j^2, 2V} and, if those fail, on a grid of 399 values of t in (0, 4V).  Exact integers.
# usage: python3 e_sweep_split.py AMAX [a:e,...]
import sys, time, math
from math import comb
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
AM=int(sys.argv[1])
jobs=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [(a,e) for a in range(3,AM+1) for e in range(3,a-1)]
st=Counter(); ex=[]; t0=time.time(); last=t0
for n,(a,e) in enumerate(jobs):
    c=row(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
    D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
    Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
    for j in range((N+1)//2,N+1):
        cj,Bj=C_(j),Bv[j]; imax=j
        for k in range(j+1,N+2):
            if cj*Bv[k]-Bj*C_(k)>0: imax=k
            else: break
        for i in range(imax+1,N+2):
            q=i-j-1
            if q<=1: continue
            st['sweep>=pi pairs']+=1
            E=D[j]-D[i]; x=2*j-N; y=2*i-N
            def Mt(t): return 0<t<4*V and 4*t*(4*V-t)*E*E>=(4*t*D[j]+(x*x-t)*Av[j]**2)*(4*t*D[i]+(y*y-t)*Av[i]**2)
            if Mt(y*y) or Mt(x*x) or Mt(2*V): st['metric at X_i^2, X_j^2 or 2V']+=1; continue
            ok=False
            for s in range(1,400):
                if Mt((4*V*s)//400): ok=True; break
            if ok: st['metric at some other t']+=1; continue
            st['NO METRIC']+=1
            if len(ex)<20: ex.append((a,e,j,i,q,x,y))
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat %d/%d rows'%(n+1,len(jobs)),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed %.1f'%(time.time()-t0)); print('no-metric examples (a,e,j,i,q,X_j,X_i):',ex)
