
import argparse
from fractions import Fraction as F
from math import comb
import sympy as S

parser = argparse.ArgumentParser(
    description="Check the uniform one-label inequalities."
)
parser.parse_args()

N,t,s,n,p,d,r,x,y = S.symbols("N t s n p d r x y")
kN,nN = (N+p)/2,(N-p)/2
prev = (d*x-(kN+1)*y)/(nN+1)
nxt2 = (d*y-nN*x)/(kN+2)
qxy = p*(kN+2)*y*y-(p+1)*d*x*y+(p+2)*(nN+1)*x*x
assert S.factor(
    (x*x-prev*y-y*y+x*nxt2)*(nN+1)*(kN+2)-qxy
) == 0

lam = (n-1)*(n+p+2)
Q = p*(n+p+2)*r*r-(p+1)*d*r+(p+2)*(n+1)
assert S.expand(
    (n-1)*Q-((p-n+1)*d*r+2*n*n-p-2)
    -p*(lam*r*r-n*d*r+n*n)
) == 0

Disc = s**4-((4*t+6)*N-4*t*t+10)*s*s+(N+3)**2
original = s*s*(N-2*t)**2-(s*s-1)*((N+3)**2-s*s)
assert S.expand(Disc-original) == 0
assert S.factor(Disc.subs(N,2*t)) == (
    (s-1)*(s+1)*(s-2*t-3)*(s+2*t+3)
)

nS,kS = (N-s+1)/2,(N+s-1)/2
lamS = (nS-1)*(kS+2)
assert S.expand(
    Disc-4*(kS+2)*(nS+1-2*s*s)
    -s*s*((N-2*t)**2-4*lamS)
) == 0

q = ((s-1)*(N+s+3)*r*r/2-s*(N-2*t)*r
     +(s+1)*(N-s+3)/2)
U = (s+1)/(s-1)
b = (4*t+2)/3
C = s*s+(2*t+4)*s+3
assert S.factor(q.subs(r,1)-2*(t+1)*s) == 0
assert S.factor(q.subs(r,U)-2*(t+2)*s*U) == 0
assert S.factor(S.diff(q,r).subs(r,U)-(N+C)) == 0
margin = S.factor(
    2*(t+2)-b*(S.Rational(4,3)+S.Rational(5,4)/t)
)
assert margin == (4*t*t+26*t-15)/(18*t)
print("symbolic identities: 9 passed; odd margin =",margin)

def coefficient(a,t,j):
    if j < 0 or j > a+t:
        return 0
    lo,hi = max(0,j-a),min(t,j)
    if lo > hi:
        return 0
    choose = comb(a,j-lo)
    total = 0
    for i in range(lo,hi+1):
        total += (-1)**i*comb(t,i)*choose
        if i < hi:
            numerator = choose*(j-i)
            denominator = a-j+i+1
            assert numerator % denominator == 0
            choose = numerator//denominator
    return total

def low_coefficients(t,N):
    older,current = [1,0,0,0],[0,1,0,0]
    for j in range(1,t):
        aj = j*(N-j+1)
        newer = [
            (current[k-1] if k else 0)-aj*older[k]
            for k in range(4)
        ]
        older,current = current,newer
    return current

root_checks = 0
for tt in (4,5,6,7,16,17,32,33,64,65):
    for NN in (2*tt,12*tt,147*tt):
        low = low_coefficients(tt,NN)
        if tt % 2 == 0:
            actual = -F(low[2],low[0])
            bound = F(tt,2*(NN-tt+2))
        else:
            actual = -F(low[3],low[1])
            bound = F(tt-1,4*(NN-tt+3))
        assert 0 < actual <= bound
        root_checks += 1
print("reciprocal-root-sum checks:",root_checks)

points = set()
for tt in (4,5,16,17,32,33,64,65):
    for pp in (1,2,6,7):
        ss = pp+1
        for base in (2*tt,12*tt,3*tt*ss*ss,6*tt*ss*ss):
            NN = base+(pp-base)%2
            if 3*pp < NN-2:
                points.add((tt,NN,pp))
    NN = 128*tt
    pp = NN//4
    pp += (NN-pp)%2
    points.add((tt,NN,pp))

counts = dict(band=0,outer=0,central_even=0,central_odd=0)
central_checks = 0
for tt,NN,pp in sorted(points):
    ss = pp+1
    nn,kk,dd = (NN-pp)//2,(NN+pp)//2,NN-2*tt
    cm,ck,ck1,ck2 = [
        coefficient(NN-tt,tt,kk+i) for i in (-1,0,1,2)
    ]
    delta = ck*ck-cm*ck1-ck1*ck1+ck*ck2
    assert delta >= 0
    numerator = (
        pp*(kk+2)*ck1*ck1
        -ss*dd*ck*ck1+(pp+2)*(nn+1)*ck*ck
    )
    assert delta*(nn+1)*(kk+2) == numerator

    ff = ss*ss*dd*dd-4*pp*(pp+2)*(nn+1)*(kk+2)
    ll = (nn-1)*(kk+2)
    if ff <= 0:
        counts["band"] += 1
    elif dd*dd >= 4*ll:
        counts["outer"] += 1
        ratio = F(ck1,ck)
        assert ratio > 0
        assert ll*ratio*ratio-nn*dd*ratio+nn*nn >= 0
        assert ratio <= F(nn*dd,2*ll)
    else:
        assert NN > 3*tt*ss*ss
        key = "central_odd" if tt % 2 else "central_even"
        counts[key] += 1

    if NN >= 3*tt*ss*ss:
        assert ck != 0
        ratio = F(ck1,ck)
        qq = (
            pp*(kk+2)*ratio*ratio-ss*dd*ratio
            +(pp+2)*(nn+1)
        )
        if tt % 2 == 0:
            assert 0 < ratio <= 1
            assert qq >= 2*(tt+1)*ss
        else:
            U = F(ss+1,ss-1)
            b = F(4*tt+2,3)
            R = U*(1-b*F(ss,NN))
            margin = F(4*tt*tt+26*tt-15,18*tt)
            assert ratio >= R
            assert qq >= U*ss*margin
        central_checks += 1

print("direct coefficient points:",len(points))
print("coverage counts:",counts)
print("central ratio and margin checks:",central_checks)

N,t,s,z,w = S.symbols("N t s z w")
disc_z = z*z-((4*t+6)*N-4*t*t+10)*z+(N+3)**2
upper = (
    (3*t*t+18*t-1)*w*w
    +(20*t*t+126*t+2)*w+32*t*t+216*t+15
)
assert S.expand(
    -disc_z.subs(N,3*t*z).subs(z,4+w)-upper
) == 0
assert S.expand(
    disc_z.subs(N,4*z)+(8*t+7)*z*z-14*z-9
    -4*t*z*(t-2*z)
) == 0

H = 4*(N-t+3)-(t-1)*(s-1)**2
assert S.expand(
    H-3*N-(N-3*t*s*s)
    -((2*t+1)*(s*s-4)+2*(t-1)*(s-2)+7*t+13)
) == 0
assert S.expand(
    2*(N-t+2)-t*(s+1)**2
    -2*(N-3*t*s*s)-t*(s-1)*(5*s+3)-4
) == 0
assert S.expand(
    N-s*s-(2*t+2)*s+3
    -(N-3*t*s*s)-(2*t-2)*s*s-(t+1)*s*(s-2)-3
) == 0
b = (4*t+2)/3
assert S.expand(
    N-b*s*(s+1)-(N-3*t*s*s)
    -s*((5*t-2)*(s-2)+6*(t-1))/3
) == 0
print("endpoint and central-bound identities: 6 passed")

def monic_kraw(degree,value,size):
    if degree == 0:
        return 1
    older,current = 1,value
    for j in range(1,degree):
        older,current = (
            current,value*current-j*(size-j+1)*older
        )
    return current

dual_checks = reflection_checks = 0
for tt,NN,pp in sorted(points):
    kk,nn = (NN+pp)//2,(NN-pp)//2
    ck = coefficient(NN-tt,tt,kk)
    ck1 = coefficient(NN-tt,tt,kk+1)
    falling = 1
    for j in range(tt):
        falling *= NN-j
    assert ck*falling == (
        (-1)**tt*comb(NN,kk)*monic_kraw(tt,pp,NN)
    )
    assert ck == (-1)**tt*coefficient(NN-tt,tt,nn)
    assert ck1 == (-1)**tt*coefficient(NN-tt,tt,nn-1)
    dual_checks += 1
    reflection_checks += 2
print("duality checks:",dual_checks)
print("reflection checks:",reflection_checks)
