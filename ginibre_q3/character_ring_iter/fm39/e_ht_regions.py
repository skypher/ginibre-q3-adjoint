# Two-region plan for (HT*) on case-(c) pairs (rows a >= e+2, e >= 3):
#   inner (X_i^2 <= 2V): the rigorous chain (b): sum_k sin(dth_k^(M_k))/sqrt(M_k) >= 1/sqrt(mu_i), M_k = mu_k (M_j = mu_(j+1)
#                         when X_j = 0, where n_j = 0);  then (HT*) holds (fm39/e_ht_sinsum.py for the derivation);
#   (b') same with the exact ratio perm/(rho_j rho_i) in place of 1/sqrt(mu_i) (rho^2 = m^2 + mu n^2 = 4V D / sc^2);
#   outer (X_i^2 > 2V):  (HT) at t = 2V:  16V^2 (D_j - D_i)^2 >= F_j(2V) F_i(2V)  (Q_(2V) unimodular).
# Reports failures of each half.  usage: python3 e_ht_regions.py AMAX  |  python3 e_ht_regions.py 0 a:e,...
import sys, math, time
from math import comb
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
def ang(m,n,M): return math.atan2(math.sqrt(M)*n, m)
AM=int(sys.argv[1])
jobs=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [(a,e) for a in range(5,AM+1) for e in range(3,a-1)]
import os; HB=float(os.environ.get('HB','30'))
st=Counter(); exI=[]; exO=[]; t0=time.time(); last=t0
for n_,(a,e) in enumerate(jobs):
    trow=time.time()
    c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
    D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
    Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
    sc=max(abs(x) for x in c) or 1
    mm=[(sg*C_(k)-d*Bv[k]/2)/sc for k in range(N+2)]; nn=[((2*k-N)*Av[k]/2)/sc for k in range(N+2)]
    mu=[(4*V/(2*k-N)**2-1) if 2*k!=N else math.inf for k in range(N+2)]
    for j in range((N+1)//2, N+1):
        cj,Bj=C_(j),Bv[j]; imax=j
        for k in range(j+1,N+2):
            if cj*Bv[k]-Bj*C_(k)>0: imax=k
            else: break
        if cj==0 and Bj==0: continue          # v_j = 0: W = 0, (E) is OL
        Sb=0.0; kk=j
        for i in range(imax+1,N+2):
            if i-j-1<=1: continue
            Xi2=(2*i-N)**2
            if Xi2<=2*V:
                st['inner']+=1
                # accumulate chain sum over steps j..i-1
                Sb=0.0
                for k in range(j,i):
                    M=mu[k] if 2*k!=N else mu[k+1]
                    dt=(ang(mm[k+1],nn[k+1],M)-ang(mm[k],nn[k],M))%(2*math.pi)
                    Sb+=math.sin(dt)/math.sqrt(M)
                rj=math.sqrt(max(mm[j]**2+(0 if 2*j==N else mu[j]*nn[j]**2),0)); ri=math.sqrt(max(mm[i]**2+mu[i]*nn[i]**2,0))
                permr=(abs(mm[j]*nn[i])+abs(nn[j]*mm[i]))/(rj*ri) if rj*ri>0 else 0.0
                if Sb>=1/math.sqrt(mu[i]): st['inner: (b) holds']+=1
                elif Sb>=permr: st["inner: (b') exact ratio holds"]+=1
                else:
                    st['inner: (b) FAILS']+=1
                    if len(exI)<8: exI.append((a,e,j,i,i-j-1,round(Sb*math.sqrt(mu[i]),3)))
            else:
                st['outer']+=1
                E=D[j]-D[i]; x=2*j-N; y=2*i-N; t=2*V
                if 4*t*(4*V-t)*E*E>=(4*t*D[j]+(x*x-t)*Av[j]**2)*(4*t*D[i]+(y*y-t)*Av[i]**2): st['outer: t=2V holds']+=1
                else:
                    st['outer: t=2V FAILS']+=1
                    if len(exO)<8: exO.append((a,e,j,i,i-j-1,round(x*x/(4*V),3),round(Xi2/(4*V),3)))
    if time.time()-last>HB: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat rows %d/%d done, last row (a,e)=(%d,%d) took %.2fs, elapsed %.0fs'%(n_+1,len(jobs),a,e,time.time()-trow,time.time()-t0),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed %.1f'%(time.time()-t0))
print('inner (b) failures (a,e,j,i,q,ratio):',exI); print('outer 2V failures (a,e,j,i,q,X_j^2/4V,X_i^2/4V):',exO)
