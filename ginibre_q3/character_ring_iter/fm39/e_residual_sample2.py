# Sampled (E) coverage beyond the grid: random rows (a, e) up to AMAX, all pairs, all criteria (binary forms exact).
import sys, time, random, math
exec(open('e_residual3.py').read().split("st=Counter(); byq=Counter(); ex=[]")[0])
seed,rows,AMAX=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); random.seed(seed)
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
st=Counter(); res=[]; t0=time.time(); last=t0
for rr in range(rows):
    a=random.randint(3,AMAX); e=random.randint(3,AMAX)
    c=cv(a,e); N=a+e; d=a-e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
    Dint=[C_(k)**2-C_(k-1)*C_(k+1) for k in range(N+2)]; Dmax=max(Dint)
    Dsc=[x/Dmax for x in Dint]          # homogeneous rescaling (exact big-int true division -> float)
    Dd=lambda k: Dsc[k] if 0<=k<=N+1 else 0.0
    for j in range((N+1)//2,N+1):
        for i in range(j+1,N+2):
            st['pairs']+=1; q=i-j-1; s=i-j
            if abs(a-e)<=1 or q<=1 or (j,i)==(N-1,N+1): continue
            ld=s*(math.sqrt(Dd(j)*Dd(i-1))+math.sqrt(Dd(j+1)*Dd(i)))
            b=min(ld, chord(j,i-1,Dd,N,V)+chord(j+1,i,Dd,N,V))
            if Dd(j)-Dd(i)>=b*(1+1e-9): continue
            ok=True
            for sg in (1,-1):
                if q==2 and sg==1: continue
                A,B,Cq=form(d,N,j,q,sg)
                if not(B*B-4*A*Cq<0 and A>0): ok=False; break
            if ok: st['definite']+=1; continue
            if rf_ok(a,e,j,i): st['RF']+=1; continue
            st['RESIDUAL']+=1
            if len(res)<12: res.append((a,e,j,i))
    st['rows']+=1
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat rows %d/%d pairs %d residual %d'%(rr+1,rows,st['pairs'],st['RESIDUAL']),flush=True)
print(time.strftime('%H:%M:%S'),'done',dict(st),'elapsed',round(time.time()-t0,1)); print('residual examples (a,e,j,i):',res)
