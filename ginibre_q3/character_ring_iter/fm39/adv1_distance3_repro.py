import argparse
from math import comb
from fractions import Fraction as Q
import mpmath as mp

argparse.ArgumentParser(description="ADV-1 exact and Bessel tests; no files").parse_args()
mp.mp.dps = 60

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def c(a, e, j):
    return sum((-1)**h*C(e,h)*C(a,j-h) for h in range(e+1))

def ballot(k, n):
    return C(k,(k-n)//2)-C(k,(k-n)//2-1) if k >= n and (k-n)%2 == 0 else 0

def cat(k):
    return C(k,k//2)//(k//2+1) if k >= 0 and k%2 == 0 else 0

def terms(a, e, b, n):
    out = []
    for i in range(b+1):
        z = lambda j: sum(C(i,h)*c(a,e,j-2*h) for h in range(i+1))
        D = lambda j: z(j)**2-z(j-1)*z(j+1)
        j = (a+e+2*i+n)//2
        out.append(D(j)-D(j+1))
    return out

def F(a, e, b, n):
    return sum(C(b,i)*2**(b-i)*v for i,v in enumerate(terms(a,e,b,n)))

def pair(a, e, p, q):
    N = a+e
    D = lambda j: c(a,e,j)**2-c(a,e,j-1)*c(a,e,j+1)
    B = lambda j: c(a,e,j-1)+c(a,e,j+1)
    j, i = (N+p-q)//2, (N+p+q)//2+1
    return D(j)-D(i), B(i)*c(a,e,j)-c(a,e,i)*B(j)

print("b=2 witness:", terms(5,4,2,7), F(5,4,2,7))
for m in (8,16,32,64):
    a, e, n = 3*m, m, 2*m
    z2 = [(4,0,1),(2,2,2),(0,4,1),(2,0,-4),(0,2,-4),(0,0,4)]
    direct = sum(c(a,e,h)*w*ballot(a+e-h+dx,n)*cat(h+dy)
                 for h in range(a+e+1) for dx,dy,w in z2)
    assert F(a,e,2,n) == direct > 0
    p, q, N = m+2, m, a+e
    T, W = pair(a,e,p,q)
    cross = sum(c(a,e,h)*ballot(N-h,p)*ballot(h,q) for h in range(N+1))
    same = sum(c(a,e,h)*cat(h)*sum(ballot(N-h,j)
               for j in range(abs(p-q),p+q+1,2)) for h in range(N+1))
    assert (T,W) == (same,cross) and T >= abs(W)
    print("m, H2 relative margin:", m, float(Q(T-abs(W),T)))

count = 0
for t in range(25):
    for b in range(2,11):
        n = t+2*b-6
        if n < 7: continue
        A = 30*b-3*t-20
        B = 180*b*b-36*b*t-300*b+9*t*t-6*t+64
        K = 3*(40*b**3-12*b*b*t-120*b*b+6*b*t*t
               +48*b*t+32*b+3*t**3-18*t*t)
        for e in range(t+1):
            u = (t-2*e)**2
            assert 144*F(t-e,e,b,n) == u**3+A*u*u+B*u+K > 0
            count += 1
print("all-b distance-three checks:", count)

for t in (6,8,10,14,20):
    for e in range(t+1):
        u = (t-2*e)**2
        P = u*((Q(2)*u-3*t-8)**2+27*t*t-72*t)/4+9*t*t*(t-2)
        for m in (7,19,101):
            T,W = pair(t-e,e,m+t-6,m)
            assert W == 0 and 144*T == P
print("H2 distance-three stabilization: PASS")

def real(q):
    return mp.mpf(q.numerator)/q.denominator

def J(x):
    u = x*x
    return 16*(u*u-2*u+3)*mp.exp(-u)

def GP(x,y):
    lo, hi = (x-y)**2, (x+y)**2
    T = 4/(x*y)*((lo*lo+3)*mp.exp(-lo)-(hi*hi+3)*mp.exp(-hi))
    W = 16*mp.exp(-x*x-y*y)*(3-2*(x*x+y*y)+(x*x-y*y)**2)
    return (T+W)/96, (T-W)/96

for n in (10,20,40):
    a = (n+1)**2-1
    x = mp.mpf(n+1)/mp.sqrt(a)
    h1 = real(Q(F(a,2,2,n),(n+1)*F(a,2,2,0)))
    T,W = pair(a,2,2*n,n)
    den = 2*(2*n+1)*(n+1)*F(a,2,0,0)
    gp,gm = GP(mp.mpf(2*n+1)/mp.sqrt(a),x)
    errors = (a*(h1-J(x)/48),
              a*(real(Q(T+W,den))-gp),
              a*(real(Q(T-W,den))-gm))
    print("n,a; scaled errors:", n,a,*[mp.nstr(z,10) for z in errors])
