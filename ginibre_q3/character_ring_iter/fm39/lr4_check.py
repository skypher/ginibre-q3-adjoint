# Independent checks of LR4: (a) boundary lower bounds c on a fine y-grid, labels <= 30 (floating, dense);
# (b) the 1/25 conclusion at the cutoff, and LS with the corrected normalization at odd a, via direct phi.
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
import numpy as np, itertools, math
from fractions import Fraction as Fr
def Uv(k,y):
    a,b=np.ones_like(y),y.copy()
    if k==0: return a
    for _ in range(k-1): a,b=b,y*b-a
    return b
ys=np.linspace(-2,2,200001); worst={}
for k in range(1,31):
    j=k//2
    if k%2==0:
        hb=sum((Uv(l,ys)+(Uv(l-1,ys) if l>=1 else 0))**2 for l in range(j+1)); c=(j+1)/4
        sb=(k+1)+Uv(k,ys); cs=(k+1)/2
        worst[('h',k)]=hb.min()/c; worst[('s',k)]=sb.min()/cs
    else:
        hb=sum(Uv(l,ys)**2 for l in range(j+1)); c=(j+1)/4; worst[('h',k)]=hb.min()/c
        yy=-ys; sb=sum((Uv(l,yy)+(Uv(l-1,yy) if l>=1 else 0))**2 for l in range(j+1)); worst[('s',k)]=sb.min()/c
print('(a) min over y of boundary/claimed lower bound (must be >=1):',round(min(worst.values()),4),'at',min(worst,key=worst.get))
def ellH(k): return Fr(2*k*k*(k+1)*(k+3),3) if k%2==0 else Fr((k-1)**2*(k+2)*(k+3),6)
def ellS(p): return Fr(2*p*p) if p%2==0 else ellH(p-1)
bad=0; cnt=0
for hs in [(1,),(2,),(3,2),(1,1,2),(2,2)]:
    for ss in [(),(2,),(3,)]:
        w=ONE
        for k in hs: w=mul(w,h(k))
        for p in ss: w=mul(w,shat(p))
        t=sum(k%2 for k in hs)+sum(p%2 for p in ss)
        L=sum(ellH(k) for k in hs)+sum(ellS(p) for p in ss)
        for a in range(t%2,t%2+6,2):
            rmin=max(1,math.ceil((7*L-a-t-6)/2))
            for r in (rmin,rmin+1):
                if r>12: continue
                val=phi(r,w,a); base=phi(r,ONE,a+t); cnt+=1
                if val < base/25: bad+=1; print('LR4 FAILS',hs,ss,r,a,val,base)
print('(b1) LR4 1/25 bound at cutoff (r<=12):',cnt,'cases, failures',bad)
# LS corrected normalization for odd a: Z=(1/2)E[(x-y)^{2r}|x+y|^a] numerically
g,wts=np.polynomial.legendre.leggauss(400); th=(g+1)*np.pi/2; dens=(2/np.pi)*np.sin(th)**2*wts*np.pi/2; X=2*np.cos(th)
XX,YY=np.meshgrid(X,X); D=np.outer(dens,dens)
bad=0; cnt=0
for hs,ss in [((1,),()),((2,1),()),((3,),()),((1,),(2,)),((2,),(3,))]:
    w=ONE
    for k in hs: w=mul(w,h(k))
    for p in ss: w=mul(w,shat(p))
    degw=sum(hs)+sum(ss); dim=math.prod(math.comb(k+3,3) for k in hs)*math.prod(2*(p+1) for p in ss)
    Lam=sum(k*(k+4) for k in hs)/5+sum(p*(p+2) for p in ss)/3
    for r in (1,2,3):
        for a in range(1,30,2):
            if (degw+a)%2: continue
            Z=0.5*(D*(XX-YY)**(2*r)*np.abs(XX+YY)**a).sum()
            val=float(phi(r,w,a)); b=dim*Z*(1-(2*r+3)*Lam/(a+2*r+4)); cnt+=1
            if val < b*(1-1e-9)-1e-9: bad+=1; print('LS odd-a fails',hs,ss,r,a,val,b) if bad<4 else None
print('(b2) LS with corrected Z at odd a:',cnt,'cases, failures',bad)
