# Probe of the log-spiral heuristic for (HT): fix t, R_k = Q_t(v_k) = F_k(t)/(4t), theta_k = Q_t-polar angle of v_k.
# For case-(c) pairs report: is R monotone on [j, i]?  min over steps of log(R_k/R_(k+1)) / dtheta_k versus 2/sqrt(det Q_t),
# and the summed version log(R_j/R_i) / (theta_i - theta_j).
import sys, math
from math import comb
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
import numpy as np
for (a,e) in [(40,5),(60,7),(200,3),(200,21),(300,150)]:
    c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
    B=lambda k: C_(k-1)+C_(k+1); A=lambda k: C_(k-1)-C_(k+1)
    stats=[]; nonmono=0; tot=0; worst_step=9e9; worst_sum=9e9
    for j in range((N+1)//2,N+1):
        cj,Bj=C_(j),B(j); imax=j
        for k in range(j+1,N+2):
            if cj*B(k)-Bj*C_(k)>0: imax=k
            else: break
        for i in range(imax+1,min(N+2,imax+6)):
            if i-j-1<=1: continue
            for t in [ (2*i-N)**2, 2*V ]:
                if not (0<t<4*V): continue
                M=np.array([[1+d*d/t, -d*sg/(2*t)],[-d*sg/(2*t), -0.25+sg*sg/(4*t)]])   # Q_t as a form in (c, B)
                w,Uv=np.linalg.eigh(M); Mh=Uv@np.diag(np.sqrt(w))@Uv.T
                sc=max(abs(C_(k)) for k in range(j,i+1)) or 1
                us=[Mh@np.array([C_(k)/sc,B(k)/sc]) for k in range(j,i+1)]
                R=[float(u@u) for u in us]; th=np.unwrap([math.atan2(u[1],u[0]) for u in us])
                lam=2/math.sqrt(np.linalg.det(M))
                tot+=1
                if any(R[k+1]>R[k] for k in range(len(R)-1)): nonmono+=1
                for k in range(len(R)-1):
                    dth=th[k+1]-th[k]
                    if dth>1e-12 and R[k+1]>0: worst_step=min(worst_step, math.log(R[k]/R[k+1])/dth/lam)
                if th[-1]-th[0]>0 and R[-1]>0: worst_sum=min(worst_sum, math.log(R[0]/R[-1])/(th[-1]-th[0])/lam)
    print('(a,e)=(%d,%d): %d (pair,t) cases; R non-monotone in %d; min step rate/lambda = %.3f; min summed rate/lambda = %.3f'%(a,e,tot,nonmono,worst_step,worst_sum),flush=True)
