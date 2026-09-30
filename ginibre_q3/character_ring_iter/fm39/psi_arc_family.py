# Jupiter's growing-q family (e = 3, a -> infinity): does the short-psi-arc criterion cover it?
# psi_k = (c_k, B_k); left turns = (E) q=1 minus (EQ1, proved); ccw = OL.  Arc sweep < pi => (E) both signs (convex polygon);
# pi <= sweep <= 2pi => W <= 0 side trivial for one sign.
import math
from math import comb
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        for k in range(u,u+a+1): c[k]+=(-1)**u*comb(e,u)*comb(a,k-u)
    return c
for (a,j,i) in [(2000,1040,1050),(4000,2056,2070),(8000,4078,4099),(16000,8111,8139)]:
    e=3; c=row(a,e); N=a+e; C_=lambda k: c[k] if 0<=k<=N else 0
    B=lambda k: C_(k-1)+C_(k+1)
    ps=[(C_(k),B(k)) for k in range(j,i+1)]
    sc=max(max(abs(p[0]),abs(p[1])) for p in ps)
    psf=[(p[0]/sc,p[1]/sc) for p in ps]
    ang=lambda p,q: math.atan2(p[0]*q[1]-p[1]*q[0], p[0]*q[0]+p[1]*q[1])
    sw=sum(ang(psf[k],psf[k+1]) for k in range(len(psf)-1))
    W=B(i)*C_(j)-C_(i)*B(j); Dl=(C_(j)**2-C_(j-1)*C_(j+1))-(C_(i)**2-C_(i-1)*C_(i+1))
    print('a=%d j=%d i=%d q=%d: psi-sweep = %.4f rad (pi=%.4f); sign W = %s; (E) slack/(D_j-D_i) = %.4f'%(a,j,i,i-j-1,sw,math.pi,'+' if W>0 else '-', (Dl-abs(W))/Dl))
