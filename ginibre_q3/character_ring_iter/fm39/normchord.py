# Test the norm-chord route: T_ij >= eta (z_i - z_j) sqrt(Gamma(i,i) Gamma(j,j)) (squared form to stay exact).
from math import comb, factorial
from fractions import Fraction as Fr
def cvec(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return c
def K(l,x,n):
    if l==0: return 1
    o,c=1,x
    for q in range(1,l): o,c=c,x*c-q*(n-q+1)*o
    return c
fails=0; cnt=0; worst=None
for e in range(3,13):
    for a in range(e+2,e+40):
        N=a+e; n=N+2; c=cvec(a,e); C=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C(k)**2-C(k-1)*C(k+1)
        eta=Fr(factorial(e),4*(factorial(n)//factorial(n-e-2)))
        nu={l:factorial(l)*(factorial(n)//factorial(n-l)) for l in range(e%2,e+1,2)}
        G=lambda u: Fr(comb(n,u+1)**2)*sum(Fr(K(l,2*u-N,n)**2,nu[l]) for l in nu)
        for j in range(N//2+1,N+1):
            for i in range(j+2,N+2):
                T=D(j)-D(i); z=lambda k:(2*k-N)**2
                rhs2=(eta*(z(i)-z(j)))**2*G(i)*G(j); cnt+=1
                if T<0 or T*T<rhs2:
                    fails+=1
                    r=float(T*T/rhs2) if rhs2 else 0
                    if worst is None or r<worst[0]: worst=(r,a,e,j,i)
print('norm-chord pairs',cnt,'failures',fails,'worst ratio (T^2/rhs^2)',worst)
