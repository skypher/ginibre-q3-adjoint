# Fixed-e limit of the (E) coverage menu (see e_continuum_limit.py for the limit of D, W, A):
#   short psi-arc: W(th,tau) > 0 for all tau in (th, et];   M_i, M_j: limits of FM-MECH26 with 4V/N -> 4(e+1);
#   RF: no zero of He_l (l = e mod 2, l <= e) in (th, et);   LD fails in this scaling (i - j ~ sqrt N).
# Reports, per e, the fraction of the (th, et) grid left uncovered and sample uncovered points.
import math, numpy as np, sys
from scipy.special import eval_hermitenorm as He
def FF(e,t):
    g=np.exp(-t*t/2); return He(e,t)*g, -He(e+1,t)*g
es=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else list(range(0,21))
for e in es:
    T=math.sqrt(4*e+4)+3; M=1200; ts=np.linspace(0,T,M+1)
    F,Fp=FF(e,ts); Gv=Fp*Fp+ts*F*Fp+(e+1)*F*F
    zs=sorted(set(float(r) for l in range(e%2,e+1,2) for r in (np.polynomial.hermite_e.hermeroots([0]*l+[1]) if l>0 else []) if r>1e-12))
    unc=0; tot=0; ex=[]
    for p in range(M):
        th=ts[p]; et=ts[p+1:]
        W=th*Fp[p]*F[p+1:]-et*F[p]*Fp[p+1:]
        pos=W>0; arc=np.cumprod(pos).astype(bool)          # short arc: W>0 on all of (th, et]
        E=Gv[p]-Gv[p+1:]; Gi=Gv[p+1:]
        Mi=(et*et<4*(e+1)) & ((4*(e+1)-et*et)*E*E >= Gi*(4*et*et*Gv[p]-(et*et-th*th)*4*Fp[p]**2))
        Mj=(0<th*th<4*(e+1)) & ((4*(e+1)-th*th)*E*E >= Gv[p]*(4*th*th*Gi+(et*et-th*th)*4*Fp[p+1:]**2))
        rf=np.ones_like(et,dtype=bool)
        for z in zs: rf &= ~((th<z)&(z<et))
        cov=arc|Mi|Mj|rf
        # ignore the far tail where both G values are negligible
        rel=(Gv[p]>1e-9*Gv.max())
        if not rel: continue
        tot+=len(et); u=np.where(~cov)[0]; unc+=len(u)
        if len(u) and len(ex)<3: ex.append((round(th,3),round(float(et[u[0]]),3),round(float(et[u[-1]]),3)))
    print('e=%2d grid pairs %d uncovered %d (%.3f%%)  samples (th, et_first, et_last): %s'%(e,tot,unc,100*unc/max(tot,1),ex),flush=True)
