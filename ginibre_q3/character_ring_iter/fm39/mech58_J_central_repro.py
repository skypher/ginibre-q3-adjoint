"""FM-MECH58: uniform central-anchor region for the joint-energy criterion."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n", type=int, default=80)
args = ap.parse_args()
if args.max_n < 53:
    ap.error("--max-n must be at least 53: the finite transition proof needs it.")

# Energy, endpoint, and two-step contraction identities.
sig, de, X, p, q = S.symbols("sig de X p q")
C = sig**2-de**2
nxt = (2*de*p-(sig-X)*q)/(sig+X)
H = sig*(p*p+q*q)-2*de*p*q
Hn = sig*(nxt*nxt+p*p)-2*de*nxt*p
ell = sig*q-de*p
D = p*p-q*nxt
assert S.cancel(H-Hn-4*X*ell**2/(sig+X)**2) == 0
assert S.cancel((sig+X)*D-H-X*(p*p-q*q)) == 0
assert S.expand(H*H-C*(p*p-q*q)**2
                -(de*(p*p+q*q)-2*sig*p*q)**2) == 0
r, c = S.symbols("r c")
aa = 1-r*r
lam = r*r+c*aa
tr = 2*r*r+c*c*aa*aa
assert S.expand(lam*lam-tr*lam+r**4-c*c*(1-c)*aa**3) == 0
x = S.symbols("x")
assert S.expand((1-x/4)**2-(1-x+x*x)-x*(8-15*x)/16) == 0
assert S.expand((1-x/4)*(1+5*x/4)-(1+x)+5*x*x/16) == 0
m = S.symbols("m")
assert S.cancel(5*m*(m-1)/(8*m+4)-S.Rational(15,14)
                -5*(m-3)*(7*m+2)/(28*(2*m+1))) == 0
assert S.cancel(m*(m-1)/(2*m+1)-S.Rational(4,3)
                -(m-4)*(3*m+1)/(3*(2*m+1))) == 0
print("ENERGY AND CONTRACTION IDENTITIES PASS")

def bernstein(poly, variables):
    P = S.Poly(S.expand(poly), *variables)
    d, e = P.degree_list()
    out = []
    for i in range(d+1):
        for j in range(e+1):
            out.append(sum(
                co*S.Rational(comb(i,h),comb(d,h))
                  *S.Rational(comb(j,k),comb(e,k))
                for (h,k),co in P.terms() if h <= i and k <= j))
    return (d,e), out

# Two fixed certificates, uniform in the gap.
x,y,U,V = S.symbols("x y U V")
beta = S.Rational(24,25)
specs = [
    ("SIX", S.Rational(1,2),
     1+S.Rational(75,28)*x+S.Rational(15,14)*y
       +S.Rational(75,16)*x*x, 313600),
    ("EIGHT", S.Integer(1),
     1+S.Rational(8,3)*x+S.Rational(4,3)*y+6*x*x, 225)]
for name, ymax, K, denominator in specs:
    G = (beta-x)*(1+y)*K*K-(1+x)*(2*y*K+beta+y)
    num, den = S.fraction(S.factor(G))
    assert den == denominator
    P = num.subs(y,x+(ymax-x)*V).subs(x,U/2)
    degrees, co = bernstein(P, (U,V))
    assert degrees == (6,3) and len(co) == 28
    assert min(co) == 0 and sum(z == 0 for z in co) == 1
    assert all(z >= 0 for z in co)
    print(name, "BERNSTEIN", degrees, len(co), "MIN", min(co))

# The accepted MECH55 band certificate also proves J itself.
f = lambda z: (1+S.Rational(3,4)*z-z*z)/(1+z)**2
eta = f(x)*f((x+y)/2)
G = ((19-20*x)*(1+y)/(1+x)
     -(19+20*y)*eta**2-40*y*eta)
P = S.cancel(64*(1+x)**4*(x+y+2)**4*G)
degrees, co = bernstein(
    P.subs(y,x+(1-x)*V).subs(x,(1+U)/4), (U,V))
assert degrees == (9,5) and len(co) == 60
assert min(co) == S.Rational(462163,32)
print("BAND BERNSTEIN", degrees, len(co), "MIN", min(co))

# Gaps four and five: a quadratic form retaining endpoint directions.
R, sig, X, de, v, w = ring("sig,X,de,v,w", QQ)
def sq(n):
    return [n[0]*n[0], 2*n[0]*n[1], n[1]*n[1]]
def cross(n,m):
    return [n[0]*m[0], n[0]*m[1]+n[1]*m[0], n[1]*m[1]]
def scale(t,a):
    return [t*x for x in a]
def add(*arrays):
    return [sum(z,R.zero) for z in zip(*arrays)]
expected = {
    4: {"A":(4,825,142), "C":(3,660,142),
        "det":(4,2205,75615)},
    5: {"A":(5,1404,190), "C":(4,1170,190),
        "det":(5,3750,135375)}}
for gap in (4,5):
    Y = X+2*gap
    ns = [(R.one,R.zero), (2*de,-(sig-X))]
    den = sig+X
    for h in range(1,gap+1):
        factor = (sig-X-2*h)*(sig+X+2*h-2)
        ns.append(tuple(2*de*u-factor*z
                        for u,z in zip(ns[-1],ns[-2])))
        if h < gap:
            den *= sig+X+2*h
    nh, prev, nxt = ns[gap], ns[gap-1], ns[gap+1]
    fac = sig+X+2*gap-2
    D0 = [sig+X,-2*de,sig-X]
    H0 = [sig,-2*de,sig]
    Dh = add(scale(sig+Y,sq(nh)), scale(-fac,cross(prev,nxt)))
    Hh = add(scale(sig,sq(nh)), scale(sig*fac*fac,sq(prev)),
             scale(-2*de*fac,cross(nh,prev)))
    # This is 25*(sig+X)*den^2 times
    # (24/25)*sig*(sig+Y)*(D_0-D_gap)-Y*(H_0+H_gap).
    A,B,Q = add(
        scale(24*sig*(sig+Y)*den*den,D0),
        scale(-24*sig*(sig+X),Dh),
        scale(-25*Y*(sig+X),add(scale(den*den,H0),Hh)))
    for label, pol in (("A",A),("C",Q),("det",4*A*Q-B*B)):
        assert all(mon[2] % 2 == 0 for mon in pol)
        degree = max(mon[2]//2 for mon in pol)
        blocks = [R.zero for _ in range(degree+1)]
        for mon, co in pol.items():
            h = mon[2]//2
            mm = list(mon)
            mm[2] = 0
            mm[0] += 2*h
            base = R.from_dict({tuple(mm):co/16**h})
            for i in range(h,degree+1):
                blocks[i] += base*QQ(comb(i,h),comb(degree,h))
        coefficients = []
        for block in blocks:
            # 1 <= X <= 2*gap, integer X; sig >= 2Y.
            for xx in range(1,2*gap+1):
                pp = block.compose(sig,2*Y+v).compose(X,R(xx))
                coefficients.extend(pp.values())
            # X >= 2*gap; sig >= 4X.
            pp = block.compose(sig,4*X+v).compose(X,2*gap+w)
            coefficients.extend(pp.values())
        got = degree,len(coefficients),min(coefficients)
        assert got == expected[gap][label]
        assert all(z > 0 for z in coefficients)
        print("SHORT",gap,label,"DEG/COUNT/MIN",got)

def data(a,e):
    N = a+e
    sig, de = N+2, a-e
    C = sig*sig-de*de
    cc = [1]
    for k in range(N):
        z, rem = divmod(de*cc[k]-(N-k+1)*(cc[k-1] if k else 0), k+1)
        assert rem == 0
        cc.append(z)
    value = lambda k: cc[k] if 0 <= k <= N else 0
    B = [value(k-1)+value(k+1) for k in range(N+2)]
    D = [value(k)**2-value(k-1)*value(k+1) for k in range(N+2)]
    H = [sig*(value(k)**2+value(k-1)**2)
         -2*de*value(k)*value(k-1) for k in range(N+2)]
    return sig,de,C,value,B,D,H

rows = pairs = longs = newlongs = transition = 0
least = first = transition_least = None
for N in range(8,args.max_n+1):
    for e in range(3,(N-2)//2+1):
        a = N-e
        sig,de,C,val,B,D,H = data(a,e)
        rows += 1
        for j in range(N//2+1,N-3):
            X = 2*j-N
            crossed = False
            for i in range(j+1,N+1):
                W = val(j)*B[i]-B[j]*val(i)
                if W < 0 or (W == 0 and val(j)*val(i)+B[j]*B[i] < 0):
                    crossed = True
                gap = i-j
                if gap < 4:
                    continue
                Y = 2*i-N
                drop = D[j]-D[i]
                num = C*(sig+Y)**2*drop*drop
                den = 4*Y*Y*H[j]*H[i]
                region = 4*de <= sig and 2*X <= sig
                if region:
                    assert drop >= 0 and num >= den and drop >= abs(W)
                    pairs += 1
                    longs += crossed
                    ratio = F(num,den)
                    item = ratio,(a,e,j,i)
                    if least is None or ratio < least[0]:
                        least = item
                    m = gap//2
                    old = (sig <= 4*X <= 2*sig
                           or 3*m*(X+2*m-2) >= 16*sig)
                    newlongs += bool(crossed and not old)
                    if gap <= 7 and 4*X < sig and 2*Y > sig:
                        assert N <= 53
                        transition += 1
                        if transition_least is None or ratio < transition_least[0]:
                            transition_least = item
                elif crossed and Y*Y < C and first is None:
                    first = ((a,e,j,i), X,Y,C,F(num,den),drop-abs(W))
assert transition == 475
assert transition_least == (F(182329,26550),(8,5,7,11))
print("ROWS",rows,"REGION PAIRS",pairs,"LONG",longs,
      "BEYOND MECH55",newlongs)
print("REGION MINIMUM",least)
print("FINITE TRANSITIONS",transition,"MINIMUM",transition_least)
print("FIRST OUTSIDE REGION",first)

# Definition-level checks of the two-label reduction.
@lru_cache(None)
def moment(power,label):
    if (power+label)%2:
        return 0
    return sum((-1)**h*comb(label-h,h)
               *comb(power+label-2*h,(power+label-2*h)//2)
               // ((power+label-2*h)//2+1)
               for h in range(label//2+1))
@lru_cache(None)
def kernel(e,a,p,q):
    return sum((-1)**u*comb(e,u)*comb(a,v)
               *moment(e+a-u-v,p)*moment(u+v,q)
               for u in range(e+1) for v in range(a+1))
direct = 0
for N in range(8,17):
    for e in range(3,(N-2)//2+1):
        a = N-e
        sig,de,C,val,B,D,H = data(a,e)
        if 4*de > sig:
            continue
        for j in range(N//2+1,N-3):
            if 2*(2*j-N) > sig:
                continue
            for i in range(j+4,N+1):
                p,q = i+j-N-1,i-j-1
                W = val(j)*B[i]-B[j]*val(i)
                assert kernel(e,a,p,q) == W
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2)) == D[j]-D[i]
                direct += 1
print("CATALAN REDUCTION PAIRS",direct)

# Outer propagation still needs its initial value.
a,e,j = 21,3,16
sig,de,C,val,B,D,H = data(a,e)
crossed = False
for i in range(j+1,a+e+1):
    W = val(j)*B[i]-B[j]*val(i)
    if W < 0 or (W == 0 and val(j)*val(i)+B[j]*B[i] < 0):
        crossed = True
    if crossed and i-j >= 4:
        break
assert i == 22 and i < a+e-1
assert (2*i-a-e)**2 == 400 and C == 352
print("INITIAL-OUTER EXAMPLE",(a,e,j,i),"X_i^2",400,"C",352)

# Recorded endpoint-envelope obstruction, outside the new domain.
a,e,j,i = 48,8,29,33
sig,de,C,val,B,D,H = data(a,e)
assert C == 42**2 and 4*de > sig
X,Y = 2*j-a-e,2*i-a-e
h = F(H[i],H[j])
left = F(42-X,sig+X)-h*F(42+Y,sig+Y)
deficit = left*left-F(2*Y,sig+Y)**2*h
assert left > 0 and deficit == -F(23639809903556039,5171578111388160000)
assert C*(sig+Y)**2*(D[j]-D[i])**2 >= 4*Y*Y*H[j]*H[i]
print("ENDPOINT-ENVELOPE DEFICIT",deficit,"ACTUAL J NONNEGATIVE")
print("PASS")

