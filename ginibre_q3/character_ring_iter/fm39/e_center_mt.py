# The only fixed-e continuum gap of the menu: th = 0 (j = N/2, N even), e even, et^2 >= 4(e+1), where (M_i) needs
# X_i^2 < 4V and (M_j) needs X_j != 0.  FM-MECH26's general bound 4t(4V-t) W^2 <= F_j(t) F_i(t), 0 < t < 4V,
# F_k(t) = 4t D_k + (X_k^2 - t) A_k^2, with t = 2V closes it.  (e odd: c_(N/2) = B_(N/2) = 0, so W = 0 there.)
# Part 1: continuum check.  Part 2: exact check on rows with e even, a large, all i at j = N/2.
import math, numpy as np
from math import comb
from scipy.special import eval_hermitenorm as He
for e in range(0,41,2):
    et=np.linspace(2*math.sqrt(e+1),2*math.sqrt(e+1)+12,4001); g0=1.0
    F0=He(e,0.0); Fp0=0.0; G0=Fp0**2+(e+1)*F0**2
    ge=np.exp(-et*et/2); F=He(e,et)*ge; Fp=-He(e+1,et)*ge; Ge=Fp*Fp+et*F*Fp+(e+1)*F*F
    tau=2*(e+1)
    ok=(4*(e+1)-tau)*(G0-Ge)**2 >= 4*(e+1)*F0**2*(tau*Ge+(et*et-tau)*Fp*Fp)
    print('continuum e=%2d: M_t(tau=2(e+1)) holds on all of th=0, et in [2sqrt(e+1), +12]: %s'%(e,ok.all()))
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
for (a,e) in [(500,4),(1000,4),(2000,6),(2000,10),(4000,4),(4000,20)]:
    c=row(a,e); N=a+e; V=(a+1)*(e+1); C_=lambda k: c[k] if 0<=k<=N else 0
    D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); A=lambda k: C_(k-1)-C_(k+1); B=lambda k: C_(k-1)+C_(k+1)
    j=N//2; bad=0; n=0
    for i in range(j+2,N+2):
        X=2*i-N
        if X*X<4*V: continue
        n+=1; t=2*V; E=D(j)-D(i)
        Fj=4*t*D(j)+(0-t)*A(j)**2; Fi=4*t*D(i)+(X*X-t)*A(i)**2
        W=C_(j)*B(i)-B(j)*C_(i)
        if not (4*t*(4*V-t)*E*E >= Fj*Fi and E>=abs(W)): bad+=1
    print('exact (a,e)=(%d,%d), j=N/2, all i with X_i^2 >= 4V: %d pairs, M_t(t=2V) fails on %d'%(a,e,n,bad))
