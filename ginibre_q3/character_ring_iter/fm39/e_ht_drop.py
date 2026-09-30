# Sufficient conditions for (HT*) on case-(c) pairs, via p_k = (m_k, n_k), rho_k^2 = 4V D_k = m_k^2 + mu_k n_k^2.
#  perm = |m_j n_i| + |n_j m_i| <= rho_j rho_i / sqrt(mu_i)  (mu_i <= mu_j, bulk), so (HT*) follows from the DROP LEMMA
#     D_i / D_j <= r*(mu_i)^2,   r*(mu) = sqrt(1 + 1/mu) - 1/sqrt(mu).
#  Also: per step, sinh(dx_k) >= sin(dtheta_k)/sqrt(mu_j) (x = log rho, theta = elliptic angle in diag(1, mu_j)); with all
#  steps <= pi/2 this gives sinh(x_j - x_i) >= 2/sqrt(mu_j), enough when mu_j <= 4 mu_i.
# Report: drop-lemma failures (bulk case-(c) pairs), max elliptic step angle, and the fraction with mu_j <= 4 mu_i.
import sys, math
from math import comb
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
AM=int(sys.argv[1]); st=Counter(); ex=[]; maxstep=0; exstep=None
for a in range(5,AM+1):
    for e in range(3,a-1):
        c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
        Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
        m=lambda k: sg*C_(k)-d*Bv[k]/2; n=lambda k: (2*k-N)*Av[k]/2
        for j in range((N+1)//2,N+1):
            if 2*j==N: continue                      # X_j = 0: mu_j infinite, separate line
            cj,Bj=C_(j),Bv[j]; imax=j
            for k in range(j+1,N+2):
                if cj*Bv[k]-Bj*C_(k)>0: imax=k
                else: break
            muj=4*V/(2*j-N)**2-1
            for i in range(imax+1,N+2):
                if i-j-1<=1: continue
                Xi2=(2*i-N)**2
                if Xi2>=4*V: st['tail (mu_i <= 0)']+=1; continue
                st['bulk case-(c)']+=1
                mui=4*V/Xi2-1
                rs=math.sqrt(1+1/mui)-1/math.sqrt(mui)
                if D[i]/D[j] <= rs*rs: st['drop lemma holds']+=1
                else:
                    st['drop lemma FAILS']+=1
                    if len(ex)<8: ex.append((a,e,j,i,round(D[i]/D[j],4),round(rs*rs,4),round(muj,2),round(mui,2)))
                if muj<=4*mui: st['mu_j <= 4 mu_i']+=1
                # elliptic step angles in metric diag(1, mu_j)
                sc=max(abs(m(k))+abs(n(k)) for k in range(j,i+1)) or 1
                th=[math.atan2(math.sqrt(muj)*n(k)/sc, m(k)/sc) for k in range(j,i+1)]
                for k in range(len(th)-1):
                    dt=(th[k+1]-th[k])%(2*math.pi)
                    if dt>maxstep: maxstep=dt; exstep=(a,e,j,i,j+k)
for k,v in st.items(): print(k,v)
print('max elliptic step angle (metric j): %.4f rad at (a,e,j,i,k) = %s'%(maxstep,exstep))
print('drop-lemma failures (a,e,j,i,D_i/D_j,r*^2,mu_j,mu_i):',ex)
