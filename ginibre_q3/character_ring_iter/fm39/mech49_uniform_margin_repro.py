import argparse, random
from fractions import Fraction as Q
from math import ceil, log, comb

p = argparse.ArgumentParser(
    description="FM-MECH49: exact constants, edge bounds, "
                "thresholds and consumer bridges")
p.parse_args()

def sin_poly(x, last):
    term = x
    out = term
    for j in range(1, last+1):
        term *= -x*x/Q((2*j)*(2*j+1))
        out += term
    return out

upper = sin_poly(Q(9,4), 4)
lower = sin_poly(Q(9,16), 1)
assert upper < Q(3,2)*lower
assert Q(832,1195) < Q(7,10)
assert Q(4,990) > Q(1,250)
assert Q(1,220) > Q(1,250)
assert Q(250,251)*Q(502,501)**2 < 1
assert 2+Q(249,500)+Q(249*248,6*250*250) > Q(5,2)
assert Q(251,250)**2 < 2
assert 2**111 < 5**48 and 2**30 < 5**16
print("Sine certificate slack:", Q(3,2)*lower-upper)
print("Scalar constants: PASS")

for M in range(4,501):
    L = M//4
    coeff = [0]*M
    for j in range(L):
        coeff[j+L] += 5
        coeff[j+2*L] += 5
    for j in range(M-L,M):
        coeff[j-L] += 5
        coeff[j-2*L] += 5
    assert max(coeff) <= 10
    assert all(coeff[j] == 0
               for j in list(range(L))+list(range(M-L,M)))
    for j in range(L,M-L):
        assert Q(M*M-1-(M-1-2*j)**2,2) >= Q(M*M,4)
print("Middle-weight combinatorics:", 497, "dimensions")

def coefficient(k,t):
    K = k+1
    delta = Q(1,2*K*K)
    rho = 1-delta
    J = (rho*rho-Q(4,25))/2
    raw = Q(k*k*K,2)
    b0 = 3*K*K-2
    A = Q(33*(t+1),32)+1+delta+b0
    return raw*raw/delta**2*(Q(1,8)+J)*A*(A+1)

def good(k,t):
    c = coefficient(k,t)
    return c.numerator*250**(t-2) < c.denominator*251**(t-2)

def threshold(k):
    lo, hi = 500, 1000
    def approx(t):
        c = coefficient(k,t)
        return (log(c.numerator)-log(c.denominator)
                -(t-2)*log(Q(251,250)))
    while approx(hi) >= 0:
        hi *= 2
    while hi > lo+1:
        mid = (lo+hi)//2
        if approx(mid) < 0:
            hi = mid
        else:
            lo = mid
    t = hi
    while not good(k,t):
        t += 1
    while t > 500 and good(k,t-1):
        t -= 1
    assert good(k,t)
    assert t == 500 or not good(k,t-1)
    H = 1+ceil(Q((k+1)**2*(t-1),16))
    return t,H

for k in (5,6,7,10,20,100,1000,10**6):
    t,H = threshold(k)
    assert t <= 4000*k.bit_length()
    print("k =",k,"T0 =",t,"H0 =",H)

count = 0
for n in range(3,65):
    M = n+1
    h = [Q(0)]*(n//2+1)
    for j in range(n//2+1):
        h[n//2-j] = Q((-1)**j*comb(n-j,j)*2**(n-2*j))
    eps = n%2

    def val(z):
        ans = Q(0)
        for a in reversed(h):
            ans = ans*z+a
        return ans

    def deriv(z):
        ans = Q(0)
        for j in range(len(h)-1,0,-1):
            ans = ans*z+j*h[j]
        return ans

    zs = [Q(0),Q(1,4),Q(1,2),Q(3,4)]
    zs += [1-Q(c,M*M) for c in
           (Q(1,100),Q(1,2),Q(1),Q(2),Q(4),Q(16))]
    for q in [Q(i,8) for i in range(9)]:
        end = q**eps*val(q*q)
        for z in zs:
            if not 0 <= z <= 1:
                continue
            bound = 1-min(Q(M*M)*(1-z),Q(1,2))/125
            for sign in (-1,1):
                den = M+sign*end
                if den:
                    num = val(z)+sign*q**eps*val(q*q*z)
                else:
                    num = eps*val(z)+2*z*deriv(z)
                    den = Q(n*(n+1)*(n+2),3)
                assert z**eps*num*num <= bound*bound*den*den, (
                    n,q,z,sign)
                count += 1
print("Exact pair-profile checks:", count)

def U(n,x):
    a,b = Q(1),x
    if n == 0:
        return a
    for j in range(2,n+1):
        a,b = b,x*b-a
    return b

def hcoeff(n):
    h = [Q(0)]*(n//2+1)
    for j in range(n//2+1):
        h[n//2-j] = Q((-1)**j*comb(n-j,j)*2**(n-2*j))
    return h

def residual(n,sign,z,q):
    h = hcoeff(n)
    if n%2:
        return sum(c*(1+sign*q**(2*j+1))*z**j
                   for j,c in enumerate(h))/(2*(1+sign*q))
    if sign == 1:
        return sum(c*(1+q**(2*j))*z**j
                   for j,c in enumerate(h))
    return sum(h[j]*(1-q**(2*j))*z**(j-1)
               for j in range(1,len(h)))/(4*(1-q*q))

rng = random.Random(49)
for rep in range(400):
    factors = [(rng.randrange(3,11),rng.choice((-1,1)))
               for _ in range(rng.randrange(1,7))]
    k = max(n for n,sign in factors)
    K = k+1
    m = sum(sign == -1 for n,sign in factors)
    r = max(1,(m+1)//2)+rng.randrange(3)
    b = rng.randrange(4)
    a = rng.randrange(5)
    sextra = sum(
        (sign == -1 and n%2 == 0)
        or (sign == 1 and n%2 == 1)
        for n,sign in factors)
    if (a+sextra)%2:
        a += 1
    A = a+sextra
    E = 2*r
    N = (A+E)//2
    nus = [
        Q(1) if sign == -1 and n%2 == 0
        else Q(1,2) if n%2 else Q(0)
        for n,sign in factors]
    R = max(1,ceil(sum(nus)))
    T = sum(Q((n+1)**2,K*K) for n,sign in factors)
    delta = Q(1,2*K*K)
    rho = 1-delta
    assert N >= R and T >= 1
    assert delta*R <= T/32+delta

    q = Q(rng.randrange(-6,7),7)
    t = Q(2,3)
    z = t*t
    x,y = 2*t,2*q*t
    s,d = x+y,x-y
    Z = x*x+y*y-2
    lhs = d**E*s**a*Z**b
    rhs = d**E*s**A*Z**b
    prod = Q(1)
    for n,sign in factors:
        pair = U(n,x)+sign*U(n,y)
        lhs *= pair/d if sign == -1 else pair
        f = residual(n,sign,z,q)
        rhs *= f
        boundary = residual(n,sign,Q(1),q)
        assert boundary > 0
        prod *= f/boundary
        ze = 1-delta/2
        assert (residual(n,sign,ze,q)/boundary
                >= 1-2*(n+1)**2*(1-ze))
    assert lhs == rhs
    C = Q(k*k*K,2)
    S = C*C*Q(250,251)**(T.numerator//T.denominator-2)
    assert z**N*abs(prod) <= S*rho**(N-R)*z
print("Consumer ray/resource bridges:",400)

def moment(j):
    return 0 if j < 0 or j%2 else comb(j,j//2)//(j//2+1)

def row(a,e):
    c = [1]
    for sign,times in ((1,a),(-1,e)):
        for _ in range(times):
            c = [
                (c[j] if j < len(c) else 0)
                +sign*(c[j-1] if j else 0)
                for j in range(len(c)+1)]
    return c

def direct(e,a,n):
    c = row(a,e)
    N = a+e
    ans = 0
    for j,cj in enumerate(c):
        for h in range(n//2+1):
            x = N-j+n-2*h
            uj = (-1)**h*comb(n-h,h)
            ans += cj*uj*(
                moment(x+2)*moment(j)
                +moment(x)*moment(j+2)
                -2*moment(x)*moment(j))
    return ans

count = 0
for e in range(1,9):
    for a in range(11):
        c = row(a,e)
        N = a+e
        def cv(j):
            return c[j] if 0 <= j <= N else 0
        def B(j):
            return cv(j-1)+cv(j+1)
        def D(j):
            return cv(j)**2-cv(j-1)*cv(j+1)
        for n in range(3,11):
            if (N+n)%2:
                continue
            j = (N+n-2)//2
            i = j+3
            ans = D(j)-D(i)+B(i)*cv(j)-cv(i)*B(j)
            assert direct(e,a,n) == ans, (e,a,n)
            assert ans >= 0
            count += 1
print("Direct q=2 consumer bridges:",count)
print("PASS")

