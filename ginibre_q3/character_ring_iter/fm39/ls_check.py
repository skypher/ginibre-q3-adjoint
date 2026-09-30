# Independent checks of FM-MECH12: (a) E_{r,a}[4-|s|] <= 4(2r+3)/(a+2r+4) numerically (weight (x-y)^{2r}|s|^a on semicircles);
# (b) the quantitative bound (3)/(4) phi_r(w h_1^a) >= d_w Z_{r,a} [1 - (2r+3) Lambda/(a+2r+4)] exactly, for words with h and hat S;
# (c) identities (7),(8).
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
import itertools, math
from fractions import Fraction as Fr
# (b) exact: Z_{r,a} = (1/2)E[(x-y)^{2r}|s|^a] = phi_r(h_1^a) for a even; for a odd use |s|: skip odd a (words with deg w + a even and a even)
def dimw(hs,ss): return math.prod(math.comb(k+3,3) for k in hs)*math.prod(2*(p+1) for p in ss)
bad=0; cnt=0; tight=None
for r in range(1,5):
    for hs in itertools.combinations_with_replacement(range(1,5),2):
        for ss in [(),(2,),(3,),(2,2)]:
            w=ONE
            for k in hs: w=mul(w,h(k))
            for p in ss: w=mul(w,shat(p))
            degw=sum(hs)+sum(ss)
            Lam=Fr(sum(k*(k+4) for k in hs),5)+Fr(sum(p*(p+2) for p in ss),3)
            for a in range(0,40,2):
                if (degw+a)%2: continue
                val=phi(r,w,a); Z=phi(r,ONE,a)
                bound=dimw(hs,ss)*Z*(1-Fr(2*r+3,a+2*r+4)*Lam)
                cnt+=1
                if val<bound: bad+=1; print('BOUND FAILS',r,hs,ss,a,val,bound) if bad<4 else None
                if bound>=0 and val<0: print('THEOREM FAILS',r,hs,ss,a)
print('(b) quantitative bound (3)/(4):',cnt,'cases, failures',bad)
# (c) identities (7),(8) as polynomials: h_{2j+1}(x,y) = (x+y)[h_j(x^2-2,y^2-2)+h_{j-1}(x^2-2,y^2-2)]; h_{2j}(2,y)=sum_{l<=j}(U_l(y)+U_{l-1}(y))^2; h_{2j+1}(2,y)=(2+y) sum U_l(y)^2
import sympy as sp
x,y,t=sp.symbols('x y t')
def hp(k,X,Y): return sp.expand(sp.series(1/((1-X*t+t**2)*(1-Y*t+t**2)),t,0,k+1).removeO().coeff(t,k))
def Up(l,Y): return sp.chebyshevu(l,Y/2) if l>=0 else 0
ok=True
for j in range(0,5):
    lhs=hp(2*j+1,x,y); rhs=sp.expand((x+y)*(hp(j,x**2-2,y**2-2)+(hp(j-1,x**2-2,y**2-2) if j>=1 else 0)))
    ok&= sp.expand(lhs-rhs)==0
    ok&= sp.expand(hp(2*j,2,y)-sum((Up(l,y)+Up(l-1,y))**2 for l in range(j+1)))==0
    ok&= sp.expand(hp(2*j+1,2,y)-(2+y)*sum(Up(l,y)**2 for l in range(j+1)))==0
print('(c) identities (7),(8) for j<=4:',ok)
# (a) numeric check of E[4-|s|] bound via exact moments: E[(4-s) s^a (x-y)^{2r}] / E[s^a (x-y)^{2r}] on s>0 half = (4*Z_a - Z_{a+1}')...
# use symmetry: for a even, E[|s|^a d^{2r}] = E[s^a d^{2r}], and E[|s|^{a+1} d^{2r}] computed via numerical quadrature
import numpy as np
def Ez(r,a,n=1200):
    g,wts=np.polynomial.legendre.leggauss(n); th=(g+1)*np.pi/2; wt=wts*np.pi/2
    X=2*np.cos(th); dens=(2/np.pi)*np.sin(th)**2*wt
    XX,YY=np.meshgrid(X,X); D=np.outer(dens,dens)
    S=XX+YY; W=D*(XX-YY)**(2*r)*np.abs(S)**a
    return (W*(4-np.abs(S))).sum()/W.sum()
worst=0
for r in range(1,5):
    for a in range(0,30,3):
        e=Ez(r,a); b=4*(2*r+3)/(a+2*r+4); worst=max(worst,e/b)
print('(a) max ratio E[4-|s|]/bound over r<=4, a<=27:',round(worst,4))
