import argparse
from fractions import Fraction as Q
from math import comb, factorial, isqrt

argparse.ArgumentParser(
    description="FM-MECH75 exact concentration-model verifier"
).parse_args()

# Constants in both analytic estimates.
c = Q(99, 100)
assert c**7 > Q(1, 4)
assert Q(256, 15) + Q(64, 81) < 18
assert (Q(64, 45) + Q(2, 3)
        + Q(115712, 3375) + Q(32, 15)) < 43
assert Q(5829, 128) < 46
assert 17*2**17 + 1300 < 2**22
assert Q(46,256)+Q(6,32)+Q(factorial(10),2**30) < Q(1,2)
assert Q(323,2048)+Q(6,32)+Q(factorial(10),2**52) < Q(1,2)
for ell in range(2, 100):
    assert 16*ell <= 3*4**ell
print("Relative-error constants: PASS")

# Critical quartic identity, with v=y^2.
v = [Q(1), -Q(1,2), Q(1,16)]
prod = [v[j] if j < 3 else Q(0) for j in range(4)]
for j in range(1, 4):
    prod[j] += v[j-1]/2
assert prod == [1, 0, -Q(3,16), Q(1,32)]
print("Critical quartic identity: PASS")

def C(n, j):
    return comb(n, j) if 0 <= j <= n else 0

def exact_word(e, a, b, p):
    assert (e+a+p) % 2 == 0
    T = e+a
    row = [
        sum((-1)**i*C(e,i)*C(a,j-i)
            for i in range(max(0,j-a), min(e,j)+1))
        for j in range(T+1)
    ]
    cats = [C(2*j,j)//(j+1) for j in range(T//2+b+1)]
    ballot = [
        C(t,(t-p)//2)-C(t,(t-p)//2-1)
        if t >= p and (t-p) % 2 == 0 else 0
        for t in range(T+2*b+1)
    ]
    return 2*sum(
        row[j]*C(b,l)*(-2)**(b-l)
        *sum(C(l,h)*cats[j//2+h]*ballot[T-j+2*(l-h)]
             for h in range(l+1))
        for j in range(0,T+1,2)
        for l in range(b+1)
    )

# Rational enclosures only.
def sqrt_iv(x):
    unit = 2**100
    j = isqrt((x.numerator*unit*unit)//x.denominator)
    return Q(j,unit), Q(j+1,unit)

def expminus_iv(x):
    assert 0 <= x <= 2
    term = Q(1)
    s = term
    for j in range(1,16):
        term *= -x/Q(j)
        s += term
    return s, s+term*(-x)/16

def atan_iv(d):
    s = sum((Q((-1)**j,(2*j+1)*d**(2*j+1))
             for j in range(24)), Q(0))
    return s, s+Q(1,49*d**49)

p5, p239 = atan_iv(5), atan_iv(239)
pi_iv = (16*p5[0]-4*p239[1], 16*p5[1]-4*p239[0])
pi32 = (sqrt_iv(pi_iv[0]**3)[0],
        sqrt_iv(pi_iv[1]**3)[1])

def critical_model_iv(b, p, N=512):
    # a=e=2b; return G/32**b.
    M = p+1
    lo = hi = Q(0)
    step = Q(2,N)
    for j in range(N):
        y0, y1 = j*step, (j+1)*step
        v0, v1 = y0*y0, y1*y1
        H0 = (1-Q(3,16)*v0*v0+v0**3/32)**b
        H1 = (1-Q(3,16)*v1*v1+v1**3/32)**b
        Amin = Q(8*b)/(4-v0)+Q(4*b)/(2+v1)
        upper_root = sqrt_iv((4-v0)/Amin**3)[1]
        if j+1 == N:
            high = H0*upper_root
            low = Q(0)
        else:
            Amax = Q(8*b)/(4-v1)+Q(4*b)/(2+v0)
            lower_root = sqrt_iv((4-v1)/Amax**3)[0]
            elo = expminus_iv(Q(M*M)/(4*Amin))[0]
            ehi = expminus_iv(Q(M*M)/(4*Amax))[1]
            low = H1*lower_root*elo
            high = H0*upper_root*ehi
        lo += step*low
        hi += step*high
    return 2*M*lo/pi32[1], 2*M*hi/pi32[0]

for b, p in ((8,12),(16,16),(32,24),(64,32)):
    F = exact_word(2*b,2*b,b,p)
    Glo, Ghi = critical_model_iv(b,p)
    normalized = Q(F,32**b)
    assert Q(9,10)*Ghi < normalized < Q(11,10)*Glo
    print("Exact model ratio in (9/10,11/10):", (2*b,2*b,b,p))

expected = {
    6: 6344,
    7: 33408,
    8: 182816,
    12: 180399520,
    16: 190050665248,
    24: 234674980036583712,
}
for m, answer in expected.items():
    assert exact_word(m,m,m,2*m) == answer
    M = 2*m+1
    ell = (M.bit_length()-1)//2
    assert 3*M*M > 16*ell*m
print("Residual-family exact values:", expected)

p = 1024
M = p+1
b = 2**16
a = e = 2*b
ell = (M.bit_length()-1)//2
B = Q(a+e,2)+b
distance = B-Q(p,2)
assert M >= 2**10 and 3*M*M <= 16*ell*b
assert b <= min(a,e) <= max(a,e) <= 4*b
assert B < 2*p*p+2 and B < 384*M*M
assert distance == 196096
print("New proved profile:", (e,a,b,p), "distance", distance)
print("PASS")
