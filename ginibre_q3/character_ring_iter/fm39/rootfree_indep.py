# Independent checks for FM-MECH8: (a) identities (1),(2) with own code; (b) Corollary 5 against direct phi (Catalan moments).
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
from math import comb, factorial
from fractions import Fraction as Fr
def cvec(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        for j in range(a+1): c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    return c
def K(l,x,n):
    o,c=1,x
    if l==0: return 1
    for q in range(1,l): o,c=c,x*c-q*(n-q+1)*o
    return c
bad=0; cnt=0
for e in range(3,11):
    for a in range(e,e+15):
        N=a+e; n=N+2; c=cvec(a,e); C=lambda k: c[k] if 0<=k<=N else 0
        B=lambda k: C(k-1)+C(k+1); D=lambda k: C(k)**2-C(k-1)*C(k+1)
        eta=Fr(factorial(e),4*(factorial(n)//factorial(n-e-2)))
        def Gam(u,v): return comb(n,u+1)*comb(n,v+1)*sum(Fr(K(l,2*u-N,n)*K(l,2*v-N,n),factorial(l)*(factorial(n)//factorial(n-l))) for l in range(e%2,e+1,2))
        z=lambda k:(2*k-N)**2
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                W=B(i)*C(j)-C(i)*B(j); T=D(j)-D(i); cnt+=1
                if W!=eta*(z(i)-z(j))*Gam(i,j): bad+=1
                if T!=eta*sum((z(k+1)-z(k))*Gam(k,k+1) for k in range(j,i)): bad+=1
print('identities (1),(2) own code:',cnt,'pairs, mismatches',bad,flush=True)
# Corollary 5: a >= (r-1)(s+t+4)^2 - 4  =>  phi_r(h_s h_t h_1^a) >= 0 ; check at the threshold with direct phi
bad=0; cnt=0
for r in (3,4):
    for s_,t_ in ((1,1),(2,1),(2,2),(3,1)):
        a0=(r-1)*(s_+t_+4)**2-4
        for a in (a0,a0+1):
            if (a+s_+t_)%2: continue
            v=phi(r,mul(h(s_),h(t_)),a); cnt+=1
            if v<0: bad+=1
            print(f'  r={r} s={s_} t={t_} a={a}: phi={v}',flush=True)
print('Corollary 5 threshold checks:',cnt,'negatives',bad)
