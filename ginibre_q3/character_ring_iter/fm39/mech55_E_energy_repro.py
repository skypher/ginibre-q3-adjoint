"""FM-MECH55: exact verifier."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import sympy as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n", type=int, default=60)
ap.add_argument("--direct-n", type=int, default=14)
args = ap.parse_args()

# Fixed energy and transfer identities.
sig, de, X, p, q = S.symbols("sig de X p q")
t, u = (sig+X)/2, (sig-X)/2
C = sig**2-de**2
pn = (de*p-u*q)/t
H = sig*(p*p+q*q)-2*de*p*q
Hn = sig*(pn*pn+p*p)-2*de*pn*p
ell = sig*q-de*p
assert S.cancel(H-Hn-X*ell**2/t**2) == 0
assert S.cancel(2*t*(p*p-q*pn)-H-X*(p*p-q*q)) == 0
assert S.expand(H**2-C*(p*p-q*q)**2
                -(de*(p*p+q*q)-2*sig*p*q)**2) == 0
assert S.cancel(pn-de*p/sig+u*ell/(sig*t)) == 0
assert S.cancel(sig*p-de*pn-C*p/sig-de*u*ell/(sig*t)) == 0

# Operator-norm estimates.
r, c = S.symbols("r c")
A = 1-r*r
lam = r*r+c*A
trace = 2*r*r+c*c*A*A
assert S.expand(lam*lam-trace*lam+r**4-c*c*(1-c)*A**3) == 0
assert S.expand((2*r+A/8)**2-(4*r*r+A*A/16)
                -A*(32*r-3*A)/64) == 0
assert sum(F(4**k, factorial(k)) for k in range(5)) > 25
print("SYMBOLIC IDENTITIES: PASS")

# Fixed Bernstein certificate.
x, y, U, V = S.symbols("x y U V")
f = lambda z: (1+S.Rational(3,4)*z-z*z)/(1+z)**2
eta = f(x)*f((x+y)/2)
G = (19-20*x)*(1+y)/(1+x)-(19+20*y)*eta**2-40*y*eta
num = S.cancel(64*(1+x)**4*(x+y+2)**4*G)
P = S.Poly(S.expand(num.subs(y, x+(1-x)*V)
                        .subs(x, (1+U)/4)), U, V)
du, dv = P.degree(U), P.degree(V)
bern = []
for i in range(du+1):
    for j in range(dv+1):
        bern.append(sum(
            co*S.Rational(comb(i,h),comb(du,h))
              *S.Rational(comb(j,k),comb(dv,k))
            for (h,k),co in P.terms() if h <= i and k <= j))
assert (du,dv,len(bern),min(bern)) == (
    9,5,60,S.Rational(462163,32))
assert all(z > 0 for z in bern)
print("BERNSTEIN: degrees", (du,dv), "coefficients", len(bern),
      "minimum", min(bern))

@lru_cache(None)
def row(a,e):
    return tuple(
        sum((-1)**h*comb(e,h)*comb(a,k-h)
            for h in range(max(0,k-a), min(e,k)+1))
        for k in range(a+e+1))

def data(a,e):
    N=a+e
    cc=row(a,e)
    v=lambda k: cc[k] if 0 <= k <= N else 0
    B=[v(k-1)+v(k+1) for k in range(N+2)]
    D=[v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H=[(N+2)*(v(k)**2+v(k-1)**2)-2*(a-e)*v(k)*v(k-1)
       for k in range(N+2)]
    return v,B,D,H

# Bounded checks supplement the uniform proofs.
rows=energy=band=long=union=0
least=None
for N in range(args.max_n+1):
    for a in range(N+1):
        e=N-a
        v,B,D,H=data(a,e)
        sigma=N+2
        delta=a-e
        rows += 1
        for k in range(N+2):
            assert ((k+1)*v(k+1)
                    == delta*v(k)-(N-k+1)*v(k-1))
        for k in range((N+1)//2,N+1):
            xx=2*k-N
            assert ((H[k]-H[k+1])*(k+1)**2
                    == xx*(delta*v(k)-sigma*v(k-1))**2)
            assert (2*(k+1)*D[k]
                    == H[k]+xx*(v(k)**2-v(k-1)**2))
            energy += 1
        if 4*abs(delta) > sigma:
            continue
        for j in range((N+1)//2,N-2):
            xx=2*j-N
            for i in range(j+4,N+2):
                m=(i-j)//2
                in_band=(sigma <= 4*xx and 2*xx <= sigma)
                in_long=(2*xx <= sigma
                         and 3*m*(xx+2*m-2) >= 16*sigma)
                if not (in_band or in_long):
                    continue
                slack=D[j]-D[i]-abs(B[i]*v(j)-v(i)*B[j])
                assert slack >= 0, (a,e,j,i,slack)
                band += in_band
                long += in_long
                union += 1
                least=slack if least is None else min(least,slack)
print("ROWS", rows, "ENERGY CHECKS", energy)
print("REGIONS: band", band, "long", long, "union", union,
      "least slack", least)

# Independent definition-level Catalan moments.
@lru_cache(None)
def cat(h):
    return comb(2*h,h)//(h+1)

@lru_cache(None)
def moment(power,label):
    if (power+label)%2:
        return 0
    return sum(
        (-1)**h*comb(label-h,h)*cat((power+label-2*h)//2)
        for h in range(label//2+1))

@lru_cache(None)
def kernel(e,a,p,q):
    return sum(
        (-1)**u*comb(e,u)*comb(a,v)
        *moment(e+a-u-v,p)*moment(u+v,q)
        for u in range(e+1) for v in range(a+1))

direct=0
for N in range(args.direct_n+1):
    for a in range(N+1):
        e=N-a
        v,B,D,H=data(a,e)
        for j in range((N+1)//2,N-2):
            for i in range(j+4,N+2):
                p=i+j-N-1
                q=i-j-1
                W=B[i]*v(j)-v(i)*B[j]
                T=D[j]-D[i]
                assert kernel(e,a,p,q) == W
                assert sum(kernel(e,a,l,0)
                           for l in range(p-q,p+q+1,2)) == T
                direct += 1
print("CATALAN REDUCTION PAIRS", direct)

# Endpoint-envelope obstruction, using the actual energy ratio.
a,e,j,i=48,8,29,33
v,B,D,H=data(a,e)
sigma=a+e+2
xx,yy=2*j-a-e,2*i-a-e
root=42
assert root*root == sigma*sigma-(a-e)**2
h=F(H[i],H[j])
left=F(root-xx,sigma+xx)-h*F(root+yy,sigma+yy)
right2=F(2*yy,sigma+yy)**2*h
deficit=left*left-right2
W=B[i]*v(j)-v(i)*B[j]
assert left > 0 and deficit < 0
assert all(D[k] >= D[k+1] for k in range(j,i))
assert W < 0
assert D[j]-D[i] >= abs(W)
assert h == F(26305021,44590400)
assert deficit == -F(23639809903556039,5171578111388160000)
print("ENVELOPE WITNESS", (a,e,j,i), "H ratio", h)
print("ENVELOPE SQUARED DEFICIT", deficit)
print("ACTUAL DROP", D[j]-D[i], "CHORD", abs(W),
      "SLACK", D[j]-D[i]-abs(W))
print("PASS")
