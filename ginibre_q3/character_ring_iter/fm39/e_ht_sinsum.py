# Rigorous sufficient conditions for (HT*) on bulk case-(c) pairs (0 < X_j < X_i, X_i^2 < 4V), from
#   4V D_k = m_k^2 + mu_k n_k^2,  2V delta_k = det(p_k, p_(k+1)),  mu_k = 4V/X_k^2 - 1 decreasing, OL (delta >= 0):
#   per step, in any fixed metric diag(1, M) with M >= mu_k, mu_(k+1):  sinh(dx_k) >= sin(dth_k^(M)) / sqrt(M),  x = log sqrt(D);
#   sinh is superadditive, and perm <= rho_j rho_i / sqrt(mu_i).  Hence (HT*) follows from
#     (a) sum_k sin(dth_k^(mu_j)) >= sqrt(mu_j / mu_i)              [one metric, mu_j]
#     (b) sum_k sin(dth_k^(mu_k)) / sqrt(mu_k) >= 1 / sqrt(mu_i)      [local metric at each step]
# Coverage of (a), (b), and of the drop lemma itself, on the grid.
import sys, math
from math import comb
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
def ang(m,n,M):   # elliptic angle of (m,n) in metric diag(1,M)
    return math.atan2(math.sqrt(M)*n, m)
AM=int(sys.argv[1]); st=Counter(); ex=[]; fails_by=Counter()
jobs=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [(a,e) for a in range(5,AM+1) for e in range(3,a-1)]
for (a,e) in jobs:
    if True:
        c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
        Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
        sc=max(abs(x) for x in c) or 1
        mm=[(sg*C_(k)-d*Bv[k]/2)/sc for k in range(N+2)]; nn=[((2*k-N)*Av[k]/2)/sc for k in range(N+2)]
        mu=[(4*V/(2*k-N)**2-1) if 2*k!=N else math.inf for k in range(N+2)]
        for j in range((N+1)//2+1 if N%2==0 else (N+1)//2, N+1):
            if 2*j==N: continue
            cj,Bj=C_(j),Bv[j]; imax=j
            for k in range(j+1,N+2):
                if cj*Bv[k]-Bj*C_(k)>0: imax=k
                else: break
            for i in range(imax+1,N+2):
                if i-j-1<=1 or (2*i-N)**2>=4*V: continue
                st['bulk case-(c)']+=1
                Sa=0.0; Sb=0.0
                for k in range(j,i):
                    da=(ang(mm[k+1],nn[k+1],mu[j])-ang(mm[k],nn[k],mu[j]))%(2*math.pi)
                    db=(ang(mm[k+1],nn[k+1],mu[k])-ang(mm[k],nn[k],mu[k]))%(2*math.pi)
                    Sa+=math.sin(da); Sb+=math.sin(db)/math.sqrt(mu[k])
                okA = Sa>=math.sqrt(mu[j]/mu[i]); okB = Sb>=1/math.sqrt(mu[i])
                st['(a) one metric']+=okA; st['(b) local metrics']+=okB; st['(a) or (b)']+=okA or okB
                if not(okA or okB):
                    fails_by[(e, 'a<=10' if a<=10 else 'a<=20' if a<=20 else 'a<=40' if a<=40 else 'a>40')]+=1
                    if len(ex)<8: ex.append((a,e,j,i,i-j-1,round(Sb*math.sqrt(mu[i]),3)))
for k,v in st.items(): print(k,v)
print('neither (a,e,j,i,q,Sb*sqrt(mu_i)):',ex); print('failures by (e, a-range):',sorted(fails_by.items()))
