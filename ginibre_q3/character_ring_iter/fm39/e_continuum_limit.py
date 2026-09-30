# Fixed-e, a -> infinity limit of (E) (scaling X_k = 2k - N = t sqrt(N)): the row tends to F = D^e exp(-t^2/2)
# (F'' + t F' + (e+1) F = 0), D_k/S -> G = F'^2 - F F'' = F'^2 + t F F' + (e+1) F^2,
# W_ij/S -> F(th)F''(et) - F''(th)F(et) = th F'(th) F(et) - et F(th) F'(et).
# (E_inf,e):  G(th) - G(et) >= |th F'(th) F(et) - et F(th) F'(et)|,  0 <= th < et.
# Part 1: check the limit against exact rows (e = 3, 5; a large).  Part 2: scan (E_inf,e) for e <= 40.
import math, numpy as np
from math import comb
from scipy.special import eval_hermitenorm as He
def FF(e,t):  # F, F' with F = He_e(t) exp(-t^2/2), F' = -He_(e+1)(t) exp(-t^2/2)
    g=np.exp(-t*t/2); return He(e,t)*g, -He(e+1,t)*g
def G(e,t):
    F,Fp=FF(e,t); return Fp*Fp+t*F*Fp+(e+1)*F*F
def Wc(e,th,et):
    F1,F1p=FF(e,th); F2,F2p=FF(e,et); return th*F1p*F2-et*F1*F2p
# Part 1
for (a,e) in [(4000,3),(4000,5),(16000,5)]:
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    C_=lambda k: c[k] if 0<=k<=N else 0
    Dk=lambda k: C_(k)**2-C_(k-1)*C_(k+1); B=lambda k: C_(k-1)+C_(k+1)
    j=(N+1)//2+int(0.6*math.sqrt(N)/2)*1; i=j+int(1.1*math.sqrt(N)/2)
    th=(2*j-N)/math.sqrt(N); et=(2*i-N)/math.sqrt(N)
    lhs=(Dk(j)-Dk(i)); W=C_(j)*B(i)-B(j)*C_(i)
    r_disc=W/lhs; r_cont=Wc(e,th,et)/(G(e,th)-G(e,et))
    print('limit check (a,e)=(%d,%d) th=%.3f et=%.3f: W/(D_j-D_i) exact %.5f, continuum %.5f'%(a,e,th,et,r_disc,r_cont))
# Part 2
worst=[]
for e in []:
    T=math.sqrt(4*e+4)+6; ts=np.linspace(0,T,1601)
    Gv=G(e,ts); F,Fp=FF(e,ts)
    mn=(np.inf,None)
    for p in range(len(ts)-1):
        th=ts[p]; et=ts[p+1:]
        W=th*Fp[p]*F[p+1:]-et*F[p]*Fp[p+1:]
        sl=(Gv[p]-Gv[p+1:]-np.abs(W))/Gv.max()
        k=np.argmin(sl)
        if sl[k]<mn[0]: mn=(sl[k],(th,et[k]))
    worst.append((e,mn[0],mn[1]))
    print('e=%2d  min normalized slack of (E_inf) = %+.3e at (th,et)=(%.3f,%.3f)'%(e,mn[0],mn[1][0],mn[1][1]),flush=True)
# Part 3: worst ratio |W| / (G(th) - G(et)) (must be <= 1), excluding pairs where G(th) - G(et) < 1e-12 G(0)-scale
print('Part 3: max ratio |W|/(G(th)-G(et)) over 0 <= th < et <= sqrt(4e+4)+4')
for e in list(range(0,41))+[60,80,120]:
    T=math.sqrt(4*e+4)+4; ts=np.linspace(0,T,2401)
    Gv=G(e,ts); F,Fp=FF(e,ts); sc=Gv.max()
    best=(0,None)
    for p in range(len(ts)-1):
        th=ts[p]; et=ts[p+1:]
        W=th*Fp[p]*F[p+1:]-et*F[p]*Fp[p+1:]; den=Gv[p]-Gv[p+1:]
        ok=den>1e-10*sc
        if not ok.any(): continue
        r=np.where(ok,np.abs(W)/np.where(ok,den,1),0); k=np.argmax(r)
        if r[k]>best[0]: best=(r[k],(th,et[k]))
    print('e=%3d  max |W|/(G(th)-G(et)) = %.6f at (th,et)=(%.3f,%.3f)'%(e,best[0],best[1][0],best[1][1]),flush=True)
