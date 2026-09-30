# (E) on extreme rows (small e, large a; a >= e suffices since (E) is invariant under a <-> e), exact integers.
# Criteria, all proved for every a, e:
#   q <= 1 (EQ1); outer endpoint (N-1, N+1); short psi-arc (sweep < pi: W_(j,k) > 0 for j < k <= i, from EQ1 by the
#   convex-polygon argument); q = 2 with W >= 0; LD energy drop; FM-MECH26 (M_i), (M_j); Theorem RF.
# The float metric chords and the binary forms are not used, so the residual is conservative.
# Also FM-MECH26's general metric bound at t = 2V (closes the j = N/2, even-e tail; see e_center_mt.py).
# usage: python3 e_extreme_residual.py a1,a2,.. e1,e2,..   or   python3 e_extreme_residual.py pairs a:e,a:e,..   (heartbeat every 30 s)
import sys, time, math
from math import comb
exec(open('e_residual_sample2.py').read().split("from collections import Counter")[0].split("exec(open('e_residual3.py')")[0])
import numpy as np
def kraw_roots(l,n):
    if l==0: return []
    off=[math.sqrt(j*(n-j+1)) for j in range(1,l)]
    J=np.diag(off,1)+np.diag(off,-1)
    return list(np.linalg.eigvalsh(J)) if l>1 else [0.0]
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
def ld(E,s,p,q):
    z=E*E-s*s*(p+q); return z>=0 and z*z>=4*s**4*p*q
from collections import Counter
if sys.argv[1]=='pairs': jobs0=[tuple(int(y) for y in x.split(':')) for x in sys.argv[2].split(',')]
else: jobs0=[(int(a),int(e)) for a in sys.argv[1].split(',') for e in sys.argv[2].split(',')]
t0=time.time(); last=t0; tot=Counter(); ex=[]
jobs=[(a,e) for (a,e) in jobs0 if a>=e+2 and e>=3]
for n_done,(a,e) in enumerate(jobs):
    tr=time.time(); st=Counter(); byq=Counter()
    c=row(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
    D=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]+[0]
    Bv=[C_(k-1)+C_(k+1) for k in range(N+3)]; Av=[C_(k-1)-C_(k+1) for k in range(N+3)]
    roots=sorted(rt for l in range(e%2,e+1,2) for rt in kraw_roots(l,N+2))
    for j in range((N+1)//2,N+1):
        cj,Bj=C_(j),Bv[j]
        imax=j
        for k in range(j+1,N+2):
            if cj*Bv[k]-Bj*C_(k)>0: imax=k
            else: break
        for i in range(j+1,N+2):
            q=i-j-1; s=i-j; st['pairs']+=1
            if q<=1 or (j,i)==(N-1,N+1): st['q<=1/endpoint']+=1; continue
            if i<=imax: st['short psi-arc']+=1; continue
            W=cj*Bv[i]-Bj*C_(i)
            if q==2 and W>=0: st['q=2, W>=0']+=1; continue
            E=D[j]-D[i]
            if ld(E,s,D[j]*D[i-1],D[j+1]*D[i]): st['LD']+=1; continue
            x=2*j-N; y=2*i-N
            if y*y<4*V and (4*V-y*y)*E*E>=D[i]*(4*y*y*D[j]-(y*y-x*x)*Av[j]**2): st['M_i']+=1; continue
            if 0<x*x<4*V and (4*V-x*x)*E*E>=D[j]*(4*x*x*D[i]+(y*y-x*x)*Av[i]**2): st['M_j']+=1; continue
            t=2*V
            if 4*t*(4*V-t)*E*E>=(4*t*D[j]+(x*x-t)*Av[j]**2)*(4*t*D[i]+(y*y-t)*Av[i]**2): st['M_t(2V)']+=1; continue
            if not any(min(x,y)+1e-9<rt<max(x,y)-1e-9 for rt in roots): st['RF']+=1; continue
            st['RESIDUAL']+=1; byq[q]+=1
            if E<abs(W): st['(E) FAILS']+=1
            if len(ex)<40: ex.append((a,e,j,i,q,x,y))
        if time.time()-last>30:
            last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat row %d/%d (a,e)=(%d,%d) j=%d/%d pairs %d residual %d'%(n_done+1,len(jobs),a,e,j,N,st['pairs'],st['RESIDUAL']),flush=True)
    tot.update(st)
    print(time.strftime('%H:%M:%S'),'row (a,e)=(%d,%d) done in %.1fs:'%(a,e,time.time()-tr),dict(st),'residual by q',sorted(byq.items()),flush=True)
print(time.strftime('%H:%M:%S'),'ALL DONE',dict(tot),'elapsed %.1f'%(time.time()-t0)); print('residual examples (a,e,j,i,q,X_j,X_i):',ex)
