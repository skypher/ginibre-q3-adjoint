import argparse
from fractions import Fraction as Q
from math import comb
import mpmath as mp
argparse.ArgumentParser(description="ADV-2 exact and high-precision checks").parse_args()
mp.mp.dps = 70

def poly(n):
    p = [Q(0)]*(n//2+1)
    for j in range(n//2+1):
        p[n//2-j] = Q((-1)**j*comb(n-j,j)*2**(n-2*j))
    return p

def val(p,z):
    out = Q(0)
    for a in reversed(p):
        out = out*z+a
    return out

count = 0
for n in range(3,65):
    M = n+1
    ep = n%2
    p = poly(n)
    dp = [j*p[j] for j in range(1,len(p))]
    zs = [Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)]
    zs += [1-Q(c,M*M) for c in
           (Q(1,100),Q(1,2),Q(1),Q(4),Q(16))]
    for q in [Q(j,8) for j in range(9)]:
        for sg in (-1,1):
            den = M+sg*q**ep*val(p,q*q)
            for z in zs:
                if not 0 <= z <= 1:
                    continue
                num = val(p,z)+sg*q**ep*val(p,q*q*z)
                dd = den
                if not dd:
                    num = ep*val(p,z)+2*z*val(dp,z)
                    dd = Q(n*(n+1)*(n+2),3)
                R2 = z**ep*(num/dd)**2
                X2 = M*M*(1-z)
                assert R2*(1+X2/4096) <= 1
                if X2 >= 16:
                    assert R2*X2 <= 16
                count += 1

T0 = 2**21
assert Q(249,250)**2*(1+Q(25,4096)) < 1
assert Q(16,25)+Q(1,256) < 1
assert Q(512*(48*4096)**2,2**64) == Q(9,2**23)
assert (T0-12)//1024 == 2047
assert 2**83*T0*T0*Q(1,2)**2047 == Q(1,2**1922)
print("pair checks",count,"T0",T0)

def ray(h,n):
    p = [0]*(3*h+2)
    for j in range(2*h+1):
        p[h+1+j] = (-1)**j*comb(2*h,j)*2**j
    u = poly(n)
    moment = lambda j: Q(1,(j+1)*(j+2))
    den = sum((a*moment(i) for i,a in enumerate(p)),Q(0))
    num = sum((a*b*moment(i+j)
               for i,a in enumerate(p)
               for j,b in enumerate(u)),Q(0))
    return num/den

def fm(h,n):
    rows = [{0:1}]
    for _ in range(2*h):
        c = {}
        for j,a in rows[-1].items():
            for k in range(abs(j-3),j+4,2):
                c[k] = c.get(k,0)+a
        rows.append(c)
    total = 0
    for j in range(2*h+1):
        a,b = rows[j],rows[2*h-j]
        total += 2*comb(2*h,j)*(
            (a.get(n-2,0)+2*a.get(n,0)+a.get(n+2,0))*b.get(0,0)
            -2*(a.get(n-1,0)+a.get(n+1,0))*b.get(1,0)
            +a.get(n,0)*(b.get(0,0)+b.get(2,0)))
    return total

assert ray(4,20) == -Q(3277672,5980345)
assert fm(4,20) == 406
print("exact witness",ray(4,20),fm(4,20))

term = Q(1)
upper = term
for j in range(40):
    term *= -Q(20*(j+2),(2*j+3)*(2*j+2))
    upper += term
assert upper < -Q(1,50)
print("negative limit",
      mp.nstr(mp.hyp1f1(2,mp.mpf(3)/2,-5),30))

for m in (2,4,10,16):
    h,n = m*m,10*m
    M = n+1
    a = mp.mpf(5)*h/2
    t = mp.mpf(M*M)/(10*h)
    f = fm(h,n)
    scaled = mp.pi*a**5*mp.mpf(f)/(mp.mpf(8)**(2*h)*M)
    predicted = mp.exp(-t)*(t*t-2*t+3)
    assert ray(h,n) < 0 and f > 0
    print("m",m,"full/Gaussian",mp.nstr(scaled/predicted,16))

for h in (2,3,4,8,16,32):
    assert fm(h,6*h-4) == (2*h-1)*(4*h*h+26*h+6)//3
    assert fm(h,6*h+2) == 2

for K in (8,32,128):
    for M in (4,K//2,K):
        for u in map(mp.mpf,("0.1","0.5","1","3")):
            v = u/K
            finite = mp.sin(M*v)/(M*mp.sin(v))
            sphere = mp.sin(M*v)/(M*v)
            assert abs(finite-sphere) <= u*u/(5*K*K)
print("distance-three identities and finite-character errors: PASS")
