import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial as fac
import sympy as S

argparse.ArgumentParser(
    description="FM-MECH62 grouped-positivity obstruction"
).parse_args()

def C(n,k):
    return comb(n,k) if 0<=k<=n else 0
def cat(n):
    return C(2*n,n)//(n+1)
def ballot(N,p):
    if N<p or (N-p)%2: return 0
    k=(N-p)//2
    return C(N,k)-C(N,k-1)

@lru_cache(None)
def beta(m,g):
    return Q(2*fac(2*m)*fac(2*g+2),
             4**(m+g+1)*fac(m)*fac(g+1)*fac(m+g+1))

@lru_cache(None)
def R(e,a,b,p,l,j):
    assert (a+p-l)%2==0 and (e+l)%2==0
    g=(e+l)//2+j
    return sum(
        Q((-1)**k*fac(p-k)*2**(p-l-2*k),
          fac(l)*fac(k)*fac(p-l-2*k))
        *C(b-j,h)*4**h*(-1)**(b-j-h)
        *beta((a+p-l-2*k+2*h)//2,g)
        for k in range((p-l)//2+1)
        for h in range(b-j+1))

@lru_cache(None)
def I(e,l,j):
    return sum(
        Q((-1)**k*fac(2*l-2*k),
          2**l*fac(k)*fac(l-k)*fac(l-2*k))
        *Q(C(j,h)*3**h*(-1)**(j-h),
           2**j*(e+l-2*k+2*h+1))
        for k in range(l//2+1) for h in range(j+1))

def A(p,l):
    return Q(4**l*fac(l)**2*(2*l+1)*fac(p-l),
             fac(p+l+1))
def term(e,a,b,p,l,j):
    return A(p,l)*I(e,l,j)*R(e,a,b,p,l,j)*R(e,a+2,b,p,l,j)
def group(e,a,b,p,j):
    return sum(term(e,a,b,p,l,j)
               for l in range(e%2,min(p,e+2*j)+1,2))

# Independent expansion of (x-y)^e (x+y)^a (x^2+y^2-2)^b.
def direct(e,a,b,p):
    T=e+a
    row=[sum((-1)**h*C(e,h)*C(a,k-h)
             for h in range(e+1)) for k in range(T+1)]
    return sum(
        row[k]*C(b,u)*C(b-u,v)*(-2)**(b-u-v)
        *cat(k//2+v)*ballot(T-k+2*u,p)
        for k in range(0,T+1,2)
        for u in range(b+1) for v in range(b-u+1))

e=a=2; b=15; p=22
G=[group(e,a,b,p,j) for j in range(b+1)]
assert G[9]==Q(-1244088089047,58367203565363422822400)
assert [j for j,x in enumerate(G) if x<0]==[9,10,11,12]
H={l:sum(C(b,j)*8**j*term(e,a,b,p,l,j)
         for j in range(b+1) if l<=e+2*j)
   for l in range(e%2,p+1,2)}
assert H[6]==Q(-49400854497767270673,12447435037081600)
assert [l for l,x in H.items() if x<0]==[6,8,10,12,14,16]

factor=4**(e+a+1)*Q(2,3)**b
F=factor*sum(C(b,j)*8**j*G[j] for j in range(b+1))
assert F==direct(e,a,b,p)==504656
assert factor*C(b,9)*8**9*G[9]==Q(
    -517540645043552,15456714364935)
print("Balanced profile: G9 =",G[9],"H6 =",H[6],"F =",F)

# Independent polynomial integration of every channel in G9.
z=S.symbols("z")
def moment(poly,semicircle):
    return sum(
        co*(Q(cat(n//2),4**(n//2)) if semicircle else Q(1,n+1))
        for (n,),co in S.Poly(poly,z).terms() if n%2==0)
for l in range(0,21,2):
    for aa in (2,4):
        integrand=(
            z**aa*(1-z*z)**(1+l//2+9)*(4*z*z-1)**6
            *S.gegenbauer(22-l,l+1,z))
        assert moment(integrand,True)==R(2,aa,15,22,l,9)
    assert moment(
        z**2*S.legendre(l,z)*S.legendre(2,z)**9,False
    )==I(2,l,9)

# A b=3 obstruction, including the positive overall prefactor.
e,a,b,p=23,2,3,19
factor=4**(e+a+1)*Q(2,3)**b
pieces=[factor*C(b,j)*8**j*group(e,a,b,p,j)
        for j in range(b+1)]
assert pieces==[
    Q(129908066,2925),Q(-55104044146,2925),
    Q(-235766538154,2925),Q(1239529442878,975)]
assert sum(pieces)==direct(e,a,b,p)==1171913728
print("b=3 pieces:",pieces,"F =",sum(pieces))

# Valid square identity after summing output labels.
e=a=2; b=15; j=9; D=e+a+2*b
lhs=sum((p+1)*group(e,a,b,p,j)
        for p in range(D%2,D+1,2))
rhs=sum(
    C(2*(b-j),h)*4**h*(-1)**h*beta(a+1+h,e+2*j)
    for h in range(2*(b-j)+1))
assert lhs==rhs==Q(40666219,70368744177664)
assert sum((2*l+1)*I(e,l,j)
           for l in range(e%2,e+2*j+1,2))==1
for p in range(D%2,D+1,2):
    for l in range(e%2,min(p,e+2*j)+1,2):
        n=p-l
        norm=Q(fac(n+2*l+1),
               4**l*fac(n)*(n+l+1)*fac(l)**2)
        assert A(p,l)*norm==Q(2*l+1,p+1)
print("Output-summed square:",rhs)
