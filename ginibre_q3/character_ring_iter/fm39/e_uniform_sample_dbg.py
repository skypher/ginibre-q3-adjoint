# (E) residual after the UNIFORM regions only (no binary forms): strips, min<=2, q<=1 (EQ1), outer endpoint, LD energy drop,
# metric chords (e_metric), FM-MECH26 regions (M_i), (M_j), and Theorem RF.  Distribution by q.
import sys
exec(open('e_residual_sample2.py').read().split("seed,rows,AMAX=")[0])
import numpy as np
def kraw_roots(l,n):
    # monic Krawtchouk K_0=1, K_1=X, K_(j+1) = X K_j - j(n-j+1) K_(j-1): roots = eigenvalues of the Jacobi matrix
    if l==0: return []
    off=[math.sqrt(j*(n-j+1)) for j in range(1,l)]
    J=np.diag(off,1)+np.diag(off,-1)
    return list(np.linalg.eigvalsh(J)) if l>1 else [0.0]
def rf_ok(a,e,j,i):
    A,E=max(a,e),min(a,e); N=A+E; n=N+2
    Xj,Xi=2*j-N,2*i-N
    lo,hi=min(Xj,Xi),max(Xj,Xi)
    for l in range(E%2,E+1,2):
        for rt in kraw_roots(l,n):
            if lo+1e-9<rt<hi-1e-9: return False
    return True
from collections import Counter
st=Counter(); byq=Counter(); maxq=None; RES=[]
import random, time
random.seed(int(sys.argv[1])); ROWS=int(sys.argv[2]); AM=int(sys.argv[3]); t0=time.time(); last=t0
for rr in range(ROWS):
    a=random.randint(3,AM); e=random.randint(3,AM)
    if abs(a-e)<=1: continue
    if True:
        c=cv(a,e); N=a+e; d=a-e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
        Dint=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]; Dm=max(Dint); Dsc=[x/Dm for x in Dint]
        Dd=lambda k: Dsc[k] if 0<=k<=N+1 else 0.0
        Ak=lambda k: C_(k-1)-C_(k+1)
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                q=i-j-1; s=i-j
                if q<=1 or (j,i)==(N-1,N+1): continue
                st['pairs']+=1
                ld=s*(math.sqrt(Dd(j)*Dd(i-1))+math.sqrt(Dd(j+1)*Dd(i)))
                b=min(ld, chord(j,i-1,Dd,N,V)+chord(j+1,i,Dd,N,V))
                if Dd(j)-Dd(i)>=b*(1+1e-9): continue
                x=2*j-N; y=2*i-N; Dl=Dint[j]-Dint[i]
                Mi = y*y<4*V and (4*V-y*y)*Dl*Dl >= Dint[i]*(4*y*y*Dint[j]-(y*y-x*x)*Ak(j)**2)
                Mj = 0<x*x<4*V and (4*V-x*x)*Dl*Dl >= Dint[j]*(4*x*x*Dint[i]+(y*y-x*x)*Ak(i)**2)
                if Mi or Mj: st['metric M']+=1; continue
                if rf_ok(a,e,j,i): st['RF']+=1; continue
                st["RESIDUAL"]+=1; byq[q]+=1; RES.append((a,e,j,i,round((Dint[j]-Dint[i]+ (C_(i-1)+C_(i+1))*C_(j)-C_(i)*(C_(j-1)+C_(j+1)))/max(1,Dint[j]-Dint[i]),4)))
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat rows %d/%d'%(rr+1,ROWS),dict(st),flush=True)
RES.sort(); print('residual (a,e,j,i,(plus-sign value)/(D_j-D_i)):',RES[:40]); print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed',round(time.time()-t0,1)); print('residual by q:',sorted(byq.items()))
