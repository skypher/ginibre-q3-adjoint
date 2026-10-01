import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
import sympy as sp

argparse.ArgumentParser(
    description="FM-MECH52: quadratic cutoff and exact obstructions"
).parse_args()

h,u,z,nu,b = sp.symbols("h u z nu b")
s0 = 1/(1+h)
v = s0-u
P = lambda t:t*(1-t)*(1-h*t)
assert sp.factor(P(v)-P(u)-h*(v-u)*(s0*s0-u*v)) == 0
beta = 4*(1+h)
dlog = nu/z-1/(2*(1-z))-h/(2*(1-h*z))+b*beta/(beta*z-2)
rhs = (nu*(1-z)/z+b*beta*(1-z)/(beta*z-2)
       -sp.Rational(3,2)-h*(1-z)/(2*(1-h*z)))
assert sp.factor((1-z)*dlog-1-rhs) == 0
print("Reflection and IBP identities: PASS")

pair_checks=0
for iq in range(9):
    h0=Q(iq,8)**2
    s=1/(1+h0)
    for j in range(1,8):
        x=s*Q(j,16)
        y=s-x
        for N in range(2,6):
            ratio=(x/y)**(2*N)*(1-x)*(1-h0*x)/((1-y)*(1-h0*y))
            assert ratio <= (x/y)**(2*N-1)
            assert ratio <= ((1-y)/y)**2
            for ell in range(3):
                if 2*N < ell+3:
                    continue
                err=ratio*((1-x)/(1-y))**2*(y/x)**ell
                assert err <= 1
                pair_checks+=1
print("Exact weighted reflection checks:",pair_checks)

def U(n,x):
    a,b=Q(1),x
    if n==0:return a
    for j in range(2,n+1):
        a,b=b,x*b-a
    return b

norm_checks=0
for p in range(3,25):
    for iq in range(-4,5):
        q=Q(iq,4)
        for sign in (-1,1):
            den=p+1+sign*U(p,2*q)
            if den==0:
                continue
            ell=1 if p%2 else (0 if sign==1 else 2)
            L=Q(p*p-1) if ell==1 else Q(p*p,2)-(ell==2)
            for it in range(1,9):
                t=Q(it,8)
                R=(U(p,2*t)+sign*U(p,2*q*t))/den
                W=R/t**ell
                assert abs(R)<=1
                assert abs(W-1)<=L*(1-t*t)/t**ell
                norm_checks+=1
print("Exact normalized-factor checks:",norm_checks)

def C(n,k):
    return comb(n,k) if 0<=k<=n else 0

@lru_cache(None)
def cat(j):
    return C(2*j,j)//(j+1)

def ballot(m,n):
    if m<n or (m-n)%2:return 0
    q=(m-n)//2
    return C(m,q)-C(m,q-1)

def profiles(e,a,B):
    T=e+a
    c=[1]
    if T:c.append(a-e)
    for j in range(1,T):
        z=(a-e)*c[j]-(T-j+1)*c[j-1]
        assert z%(j+1)==0
        c.append(z//(j+1))
    now=[0]*(T+1)
    prev=[]
    for j in range(0,T+1,2):
        now[T-j]=c[j]*cat(j//2)
    for b in range(B+1):
        D=T+2*b
        yield [sum(now[m]*ballot(m,n) for m in range(n,D+1,2))
               for n in range(D+1)]
        get=lambda w,j:w[j] if 0<=j<len(w) else 0
        nxt=[]
        for j in range(D+3):
            num=(2*(T-2-j)*get(now,j)
                 +4*b*(get(prev,j-2)+2*get(prev,j)))
            assert num%(D+4-j)==0
            nxt.append(get(now,j-2)+num//(D+4-j))
        prev,now=now,nxt

pp=list(profiles(1,1,16))
qq=list(profiles(2,1,16))
state={(1,0):1,(0,1):-1}
walk_checks=0

def neighbours(i):
    return (1,) if i==0 else (i-1,i,i+1)

for b in range(17):
    assert all(v>=0 for (i,j),v in state.items() if i>j)
    assert all(state.get((j,i),0)==-v for (i,j),v in state.items())
    for i in range(b+2):
        assert state.get((i,0),0)==pp[b][2*i]
        walk_checks+=1
    for n,val in enumerate(qq[b]):
        expected=(pp[b][n-1] if 0<=n-1<len(pp[b]) else 0)
        expected+=(pp[b][n+1] if n+1<len(pp[b]) else 0)
        assert val==expected
    nxt=defaultdict(int)
    for (i,j),v in state.items():
        for k in neighbours(i):nxt[k,j]+=v
        for k in neighbours(j):nxt[i,k]+=v
    state={ij:v for ij,v in nxt.items() if v}
print("Chamber-walk coefficient bridges:",walk_checks)

state={(1,0):1,(0,1):-1}
sector_bridges=0
for b in range(13):
    current={(2*i,2*j):v for (i,j),v in state.items()}
    for a in range(1,9):
        assert all(v>=0 for (i,j),v in current.items() if i>j)
        assert all(current.get((j,i),0)==-v
                   for (i,j),v in current.items())
        co=list(profiles(1,a,b))[b]
        folded=list(profiles(a,1,b))[b]
        assert co==folded
        for n,val in enumerate(co):
            assert current.get((n,0),0)==val>=0
            sector_bridges+=1
        nxt=defaultdict(int)
        for (i,j),v in current.items():
            for k in ((1,) if i==0 else (i-1,i+1)):
                nxt[k,j]+=v
            for k in ((1,) if j==0 else (j-1,j+1)):
                nxt[i,k]+=v
        current={ij:v for ij,v in nxt.items() if v}
    nxt=defaultdict(int)
    for (i,j),v in state.items():
        for k in neighbours(i):nxt[k,j]+=v
        for k in neighbours(j):nxt[i,k]+=v
    state={ij:v for ij,v in nxt.items() if v}
print("All-suffix min(a,e)<=1 bridges:",sector_bridges)

@lru_cache(None)
def base_coefficient(e,a,j):
    return sum((-1)**h*C(e,h)*C(a,j-h) for h in range(e+1))

@lru_cache(None)
def mu(j,l):
    return sum(C(l,h)*(-2)**(l-h)*cat(j+h) for h in range(l+1))

def distance_value(e,a,b,p):
    D=e+a+2*b
    if p>D or (D-p)%2:return 0
    d=(D-p)//2
    return sum(
        base_coefficient(e,a,2*j)*C(b,l)*mu(j,l)*
        (C(D-2*j-2*l,d-j-l)-C(D-2*j-2*l,d-j-l-1))
        for j in range(min(d,(e+a)//2)+1)
        for l in range(min(b,d-j)+1))

def direct(e,a,b,p):
    T=e+a
    return sum(
        base_coefficient(e,a,j)*C(b,l)*(-2)**(b-l)*C(l,h)*
        cat(j//2+h)*ballot(T-j+2*(l-h),p)
        for j in range(0,T+1,2)
        for l in range(b+1)
        for h in range(l+1))

bridges=0
for e in range(1,7):
    for a in range(7):
        for b,co in enumerate(profiles(e,a,5)):
            for p in range((e+a)%2,len(co),2):
                assert co[p]==distance_value(e,a,b,p)==direct(e,a,b,p)
                assert co[p]>=0
                bridges+=1
print("Distance formula / recurrence / Catalan bridges:",bridges)

def ray(p,N,b):
    assert p%2==0
    return sum(
        Q((-1)**j*C(p-j,j)*2**(p-2*j)*C(b,h)*8**h*(-2)**(b-h),
          (p+1)*(N+p//2-j+h+1)*(N+p//2-j+h+2))
        for j in range(p//2+1) for h in range(b+1))

bad_ray=ray(10,6,3)
assert bad_ray==Q(-17093,495495)
assert direct(2,10,3,10)==11128
print("Negative ray:",bad_ray,"; consumer:",11128)

K=20
exp_upper=sum(Q(5**j,factorial(j)) for j in range(K+1))
exp_upper+=Q(5**(K+1),factorial(K+1))/(1-Q(5,K+2))
int_lower=sum(Q(5**j,factorial(j)*(2*j+1)) for j in range(K+1))
gaussian_upper=(exp_upper-9*int_lower)/(2*exp_upper)
assert gaussian_upper<0
print("Quadratic-scale negative-ray limit: upper bound",gaussian_upper)

rr=list(profiles(3,3,5))
assert (rr[4][6],rr[4][8],rr[4][10],rr[5][8])==(54,34,18,100)
assert rr[5][8]-rr[4][6]-rr[4][8]-rr[4][10]==-6
ss=list(profiles(15,16,3))
assert ss[3][7]-2*ss[2][7]==-424528
print("Actual U2-transfer failure:",-6)
print("Actual doubling failure:",ss[2][7],ss[3][7],-424528)

def cutoff(p,kind):
    if p%2:return 4*p*p-2
    return 2*p*p+(2 if kind=="plus" else -2)

for p in (7,8,20,100):
    print("p =",p,"K_plus =",cutoff(p,"plus"),
          "K_minus =",cutoff(p,"minus"))
print("PASS")
