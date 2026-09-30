# Is the HT margin rho(j,i) = max_t 4t(4V-t)E^2/(F_j F_i) nondecreasing in i (fixed j) over case-(c) pairs?
# Also for fixed t in {2V}: is 4t(4V-t)E_(j,i)^2 - F_j(t)F_i(t) nondecreasing in i?  If so, HT reduces to the first
# half-turn index i0(j) = imax(j) + 1 (or i0 = j+3 if the half-turn completes earlier).
import sys, math
from math import comb
from collections import Counter
exec(open('e_ht_margin.py').read().split("AM=int(sys.argv[1])")[0])
def rho_of(D,Av,V,N,j,i):
    E=D[j]-D[i]; x=2*j-N; y=2*i-N
    uj,ui=4*D[j]-Av[j]**2,4*D[i]-Av[i]**2; vj,vi=x*x*Av[j]**2,y*y*Av[i]**2
    sc=max(abs(uj),abs(ui),vj,vi,abs(E),1); Uj,Ui,Vj,Vi=uj/sc,ui/sc,vj/sc,vi/sc; E2=(E/sc)**2
    def ratio(t):
        den=(Uj*t+Vj)*(Ui*t+Vi); return 4*t*(4*V-t)*E2/den if den>0 else -1.0
    ts=[y*y,x*x,2*V]; p2=Uj*Ui; p1=Uj*Vi+Ui*Vj; p0=Vj*Vi; c2=-16*V*p2-4*p1; c1=-8*p0; c0=16*V*p0
    if abs(c2)>1e-300:
        disc=c1*c1-4*c2*c0
        if disc>=0: ts+=[(-c1+s*math.sqrt(disc))/(2*c2) for s in (1,-1)]
    return max([ratio(t) for t in ts if 0<t<4*V]+[-1.0])
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
            i0=max(imax+1,j+3)
            prev=None
            for i in range(i0,N+2):
                r=rho_of(D,Av,V,N,j,i); st['pairs']+=1
                if prev is not None:
                    if r<prev*(1-1e-9):
                        st['rho decreases']+=1
                        if len(ex)<8: ex.append((a,e,j,i,round(prev,3),round(r,3)))
                prev=r
print(dict(st)); print('decrease examples (a,e,j,i,rho(i-1),rho(i)):',ex)
