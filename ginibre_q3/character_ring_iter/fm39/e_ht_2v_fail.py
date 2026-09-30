# Where does t = 2V fail on case-(c) pairs?  Q_(2V) is unimodular; HT(2V) <=> D_j - D_i >= 2 sqrt(Dt_j Dt_i),
# Dt_k = Q_(2V)(v_k) = D_k + (X_k^2 - 2V) A_k^2/(8V).  Report the location (X_j^2/4V, X_i^2/4V) of failures, and whether
# the failure has Dt_k > D_k (tail side X_k^2 > 2V) at i or j.
import sys
from math import comb
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
AM=int(sys.argv[1]); st=Counter(); ex=[]
for a in range(5,AM+1):
    for e in range(3,a-1):
        c=row(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
        Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
        for j in range((N+1)//2,N+1):
            cj,Bj=C_(j),Bv[j]; imax=j
            for k in range(j+1,N+2):
                if cj*Bv[k]-Bj*C_(k)>0: imax=k
                else: break
            for i in range(imax+1,N+2):
                if i-j-1<=1: continue
                E=D[j]-D[i]; x=2*j-N; y=2*i-N; t=2*V
                Fj=4*t*D[j]+(x*x-t)*Av[j]**2; Fi=4*t*D[i]+(y*y-t)*Av[i]**2
                st['pairs']+=1
                if 4*t*(4*V-t)*E*E>=Fj*Fi: continue
                st['2V fails']+=1
                st[('X_j^2 vs 2V', 'X_j^2>2V' if x*x>2*V else 'X_j^2<=2V', 'X_i^2 vs 4V','X_i^2>=4V' if y*y>=4*V else 'X_i^2>2V' if y*y>2*V else 'X_i^2<=2V')]+=1
                st[('q',min(i-j-1,6))]+=1
                if len(ex)<8: ex.append((a,e,j,i,i-j-1,round(x*x/(4*V),3),round(y*y/(4*V),3)))
for k,v in sorted(st.items(),key=lambda z:str(z[0])): print(k,v)
print('examples (a,e,j,i,q,X_j^2/4V,X_i^2/4V):',ex)
