# Margin map of Conjecture HT on case-(c) pairs (sweep >= pi, q >= 2), rows a >= e+2, e >= 3:
#   rho = max over 0 < t < 4V of 4t(4V-t)(D_j-D_i)^2 / (F_j(t) F_i(t)),  F_k(t) = u_k t + v_k, u_k = 4D_k - A_k^2, v_k = X_k^2 A_k^2.
# (HT) holds at the pair iff rho >= 1.  The max is taken over t in {X_i^2, X_j^2, 2V} and the stationary points of the
# ratio (roots of a quadratic, found numerically in floats after exact big-int scaling).  Reports min rho per row-class.
# usage: python3 e_ht_margin.py AMAX   or   python3 e_ht_margin.py 0 a:e,a:e,...
import sys, time, math
from math import comb
from collections import defaultdict
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
AM=int(sys.argv[1])
jobs=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [(a,e) for a in range(5,AM+1) for e in range(3,a-1)]
worst=[]; t0=time.time(); last=t0; npairs=0
byclass=defaultdict(lambda:(9e9,None))
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
            if i-j-1<=1: continue
            npairs+=1
            E=D[j]-D[i]; x=2*j-N; y=2*i-N
            uj,ui=4*D[j]-Av[j]**2,4*D[i]-Av[i]**2; vj,vi=x*x*Av[j]**2,y*y*Av[i]**2
            sc=max(abs(uj),abs(ui),vj,vi,abs(E),1)
            # scale everything to floats (ratio is homogeneous of degree 0 in (E^2, F_j F_i))
            f=lambda z: z/sc
            Uj,Ui,Vj,Vi=f(uj),f(ui),f(vj),f(vi); E2=f(E)**2
            def ratio(t):
                den=(Uj*t+Vj)*(Ui*t+Vi)
                return 4*t*(4*V-t)*E2/den if den>0 else (math.inf if den==0 and E2>0 else -1.0)
            ts=[y*y,x*x,2*V]
            # stationary points of 4t(4V-t)/((Uj t+Vj)(Ui t+Vi)): numerator of derivative is a quadratic in t
            # d/dt [ (16Vt - 4t^2) / (Uj Ui t^2 + (Uj Vi + Ui Vj) t + Vj Vi) ] = 0
            p2=Uj*Ui; p1=Uj*Vi+Ui*Vj; p0=Vj*Vi
            # derivative numerator: (16V-8t)(p2 t^2+p1 t+p0) - (16Vt-4t^2)(2 p2 t + p1); expansion:
            # (16V-8t)(p2 t^2 + p1 t + p0) = 16V p2 t^2 + 16V p1 t + 16V p0 - 8p2 t^3 - 8p1 t^2 - 8p0 t
            # (16Vt-4t^2)(2p2 t + p1) = 32V p2 t^2 + 16V p1 t - 8 p2 t^3 - 4 p1 t^2
            # difference: t^2 (16V p2 - 8 p1 - 32V p2 + 4 p1) + t (16V p1 - 8 p0 - 16 V p1) + 16 V p0
            c2=-16*V*p2-4*p1; c1=-8*p0; c0=16*V*p0
            if abs(c2)>1e-300:
                disc=c1*c1-4*c2*c0
                if disc>=0:
                    for sgn in (1,-1):
                        tr=(-c1+sgn*math.sqrt(disc))/(2*c2)
                        if 0<tr<4*V: ts.append(tr)
            elif abs(c1)>1e-300:
                tr=-c0/c1
                if 0<tr<4*V: ts.append(tr)
            rho=max(ratio(t) for t in ts if 0<t<4*V) if any(0<t<4*V for t in ts) else -1.0
            key=(a,e)
            cl=('e<=5' if e<=5 else 'e<=20' if e<=20 else 'e>20')+(' a-e<=4' if a-e<=4 else '')+(' e odd' if e%2 else ' e even')
            if rho<byclass[cl][0]: byclass[cl]=(rho,(a,e,j,i,i-j-1,x,y))
            if rho<1.02 and len(worst)<30: worst.append((round(rho,5),a,e,j,i))
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat %d/%d rows, pairs %d'%(n+1,len(jobs),npairs),{k:round(v[0],4) for k,v in byclass.items()},flush=True)
print(time.strftime('%H:%M:%S'),'done: case-(c) pairs',npairs,'elapsed %.1f'%(time.time()-t0))
for k,v in sorted(byclass.items()): print('class %-14s min rho = %.5f at (a,e,j,i,q,X_j,X_i) = %s'%(k,v[0],v[1]))
print('pairs with rho < 1.02:',worst)
