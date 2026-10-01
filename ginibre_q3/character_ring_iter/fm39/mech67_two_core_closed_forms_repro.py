"""FM-MECH67: exact symbolic and character verifications; no file writes."""
import argparse
from math import comb, factorial
from fractions import Fraction as F
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--degree", type=int, default=12)
ap.add_argument("--dump", action="store_true")
args = ap.parse_args()
J = args.degree
R, N, b, t = ring("N,b,t", QQ)
M = N - 2*b

def bn(x, k):
    if k < 0:
        return R.zero
    ans = R.one
    for h in range(k):
        ans = ans*(x-h)/(h+1)
    return ans

def cat(k):
    return comb(2*k, k)//(k+1)

c = {-1:R.zero, 0:R.one}
q = c.copy()
for k in range(J+1):
    c[k+1] = (t*c[k]-(N-k+1)*c[k-1])/(k+1)
    q[k+1] = (t*q[k]-(M-k+1)*q[k-1])/(k+1)

def row(k, h):
    if k < 0:
        return R.zero
    return sum((bn(b-h,u)*c[k-2*u]
                for u in range(k//2+1)), R.zero)

def det(k, h):
    return row(k,h)**2-row(k-1,h)*row(k+1,h)

def td(d):
    if d < 0:
        return R.zero
    return sum((2**h*bn(b,h)*det(d-h,h)
                for h in range(d+1)), R.zero)

def wd(d, m):
    out = R.zero
    for h in range(d+1):
        j = d-h
        i = j-m-1
        out += 2**h*bn(b,h)*(
            (row(i-1,h)+row(i+1,h))*row(j,h)
            -row(i,h)*(row(j-1,h)+row(j+1,h)))
    return out

def vv(d, l):
    out = R.zero
    for s in range((d-l)//2+1):
        js = sum((-1)**j*comb(s,j)*2**(s-j)*cat(l+s+j)
                 for j in range(s+1))
        inner = sum((
            bn(b-s,u)*bn(N-2*l-2*s,d-l-2*s-2*u)
            for u in range((d-l-2*s)//2+1)), R.zero)
        out += bn(b,s)*inner*QQ(js,comb(2*l,l))
    return out

AA = {}
checks = 0
for d in range(J+1):
    v = [vv(d,l) for l in range(d+1)]
    A = [
        sum((v[l]*(bn(2*l-M,l-k)
                   -2*bn(2*l-M-1,l-k-1))
             for l in range(k,d+1)), R.zero)
        for k in range(d+1)]
    assert sum((A[k]*q[k]**2 for k in range(d+1)),
               R.zero) == td(d)
    if d >= 2:
        k = d-2
        num = (N*N+(2*(k+1)*b-3*k)*N
               +2*(k+1)*(k+2)*b*b
               -2*k*(6*k+5)*b+2*k*k)
        assert A[k] == num/((k+1)*(k+2)*(k+3))
        assert A[k].compose(N,R(2*d-2)).compose(
            b,R.one) == QQ(-6*(d-4),d*(d+1))
    AA[d] = A
    checks += d+1
    if args.dump:
        print(d, [str(z) for z in A])
print("diagonal coefficient identities", checks)

def basecoef(d, k, M0):
    if k < 0 or d < k:
        return R.zero
    if d == k:
        return R.one/(k+1)
    v = R(factorial(k))*(M0-2*k)/factorial(d+1)
    for h in range(d-k-1):
        v *= M0-k-h
    return v

for d in range(J+1):
    for k in range(d+1):
        assert AA[d][k].compose(b,R.zero) == basecoef(d,k,N)
        rhs = (basecoef(d,k,N-2)+2*basecoef(d-1,k,N-2)
               +2*basecoef(d-3,k,N-2)+basecoef(d-4,k,N-2)
               +24*basecoef(d-1,k+1,N)
               -30*basecoef(d,k+2,N+2))
        assert AA[d][k].compose(b,R.one) == rhs
print("b=0 and b=1 closed identities", 2*checks)

def alpha(h, r):
    return sum(((-1)**(r-u)*bn(b-h,u)*bn(h,r-2*u)
                for u in range(r//2+1)), R.zero)

def cross(d, m):
    out = {}
    def add(i, j, v):
        if i < 0 or j < 0 or not v:
            return
        if i > j:
            i,j = j,i
        if j == d+1:
            add(d,i+1,v*(i+1)/(d+1))
            add(d,i-1,v*(M-i+1)/(d+1))
            add(d-1,i,-v*(M-d+1)/(d+1))
        else:
            out[i,j] = out.get((i,j),R.zero)+v
    def product(i, j, h, v):
        if i < 0 or j < 0:
            return
        for u in range(i//2+1):
            for z in range(j//2+1):
                add(i-2*u,j-2*z,v*alpha(h,u)*alpha(h,z))
    for h in range(d+1):
        j = d-h
        i = j-m-1
        v = 2**h*bn(b,h)
        product(i-1,j,h,v)
        product(i+1,j,h,v)
        product(i,j-1,h,-v)
        product(i,j+1,h,-v)
    return out

checks = 0
for d in range(3,J+1):
    for m in range(3,d+1):
        C = cross(d,m)
        assert sum((z*q[i]*q[j] for (i,j),z in C.items()),
                   R.zero) == wd(d,m)
        assert all((i+j-m)%2 == 0
                   for (i,j),z in C.items() if z)
        checks += 1
print("bilinear cross identities", checks)

import sympy as S
p = S.symbols("p0:18")

def at(c, j):
    return c[j] if 0 <= j < len(c) else 0

def D(c, j):
    return at(c,j)**2-at(c,j-1)*at(c,j+1)

def K(c, j):
    return at(c,j)**2-at(c,j-2)*at(c,j+2)

r = [at(p,j)+at(p,j-2) for j in range(len(p)+2)]
for j in range(15):
    assert S.expand(D(r,j)+2*D(p,j-1)
                    -D(p,j)-D(p,j-2)-K(p,j-1)) == 0
print("stride-two insertion identities", 15)

if J >= 5:
    A = [S.Rational(str(z.compose(N,R(7)).compose(b,R(3))))
         for z in AA[5]]
    diag = A[:]
    diag[0] -= S.Rational(13,2)
    diag[1] -= S.Rational(1,2)
    assert all(z > 0 for z in diag)
    assert diag[2]*diag[5]-S.Rational(1,3)**2 == -S.Rational(31,360)
    G = S.diag(*diag)
    def add(i,j,z):
        z = S.sympify(z)
        G[i,j] += z/2
        G[j,i] += z/2
    add(2,5,-S.Rational(2,3))
    add(1,4,-4*(A[4]-S.Rational(1,10)))
    add(0,3,-4*(A[3]-3))
    add(0,1,18)
    H = G.copy()
    H[2,2] += 1
    H[1,3] -= S.Rational(3,4)
    H[3,1] -= S.Rational(3,4)
    H[1,1] += S.Rational(1,4)
    H[0,0] -= S.Rational(1,4)
    minors = [H[:j,:j].det() for j in range(1,7)]
    assert all(z > 0 for z in minors)
    qq = S.Matrix([
        S.sympify(str(q[k].compose(N,R(7)).compose(b,R(3))))
        for k in range(6)])
    assert S.expand((qq.T*(H-G)*qq)[0]) == 0
    print("failed minor", -S.Rational(31,360),
          "repaired leading minors", minors)

def poly(a,e,l):
    out = [
        sum((-1)**h*comb(e,h)*comb(a,j-h)
            for h in range(max(0,j-a),min(e,j)+1))
        for j in range(a+e+1)]
    for _ in range(l):
        out = [at(out,j)+at(out,j-2)
               for j in range(len(out)+2)]
    return out

def pairrow(c,j,m):
    i = j-m-1
    return (D(c,j)-D(c,i),
            (at(c,i-1)+at(c,i+1))*at(c,j)
            -at(c,i)*(at(c,j-1)+at(c,j+1)))

def pair(a,e,B,d,m):
    out = [0,0]
    for h in range(min(B,d)+1):
        v = pairrow(poly(a,e,B-h),d-h,m)
        for z in (0,1):
            out[z] += 2**h*comb(B,h)*v[z]
    return tuple(out)

p = poly(3,1,0)
full = pair(3,1,1,3,3)
base1 = pairrow(p,3,3)
base2 = pairrow(p,1,3)
remainder = tuple(full[z]-base1[z]-base2[z] for z in (0,1))
assert (full,base1,base2,remainder) == (
    (9,0),(4,-2),(4,0),(1,2))

for d in range(6,66,2):
    C = comb(d,d//2)
    R0,W0 = pair(d,d,1,d,d-2)
    T = F(2*d,d+2)*C*C+F(4,(d+2)**2)*C*C
    assert R0 == T-d-1
    assert W0 == (-1)**(d//2+1)*F(2*d,d+2)*C
    assert R0 > abs(W0)

def cc(d,k,M0):
    if k < 0 or d < k:
        return F(0)
    if d == k:
        return F(1,k+1)
    out = F(factorial(k)*(M0-2*k),factorial(d+1))
    for h in range(d-k-1):
        out *= M0-k-h
    return out

def ab1(d,k,N0):
    M0 = N0-2
    return (cc(d,k,M0)+2*cc(d-1,k,M0)
            +2*cc(d-3,k,M0)+cc(d-4,k,M0)
            +24*cc(d-1,k+1,M0+2)
            -30*cc(d,k+2,M0+4))

d = 24
qv = [(-1)**(k//2)*comb(d-1,k//2) if k%2 == 0 else 0
      for k in range(d+1)]
group = sum(ab1(d,k,48)*qv[k]**2 for k in range(d-6,d+1))
assert group == -F(358324309078,5)
assert pair(24,24,1,24,22) == (13543194541047,-4992288)
print("negative seven-term group", group,
      "full pair", pair(24,24,1,24,22))

checks = 0
for a in range(9):
    for e in range(a+1):
        for B in range(5):
            for d in range(a+e+2*B+1):
                T = pair(a,e,B,d,d+1)[0]
                assert T >= abs(at(poly(a,e,B),d))
                checks += 1
print("boundary checks", checks)

E,x = ring("x",QQ)

def bnf(n,j):
    if j < 0:
        return F(0)
    z = F(1)
    for h in range(j):
        z = z*(n-h)/(h+1)
    return z

def product(seq):
    z = E.one
    for v in seq:
        z *= v
    return z

def edge(B,r,g):
    out = E.zero
    for j in range(r+1):
        bb = (bnf(2*j-2*r+4*B-g,j)
              -2*bnf(2*j-2*r+4*B-g-1,j-1))
        for s in range((r-j)//2+1):
            v = bb*bnf(B,s)*sum(
                bnf(B-s,u)*bnf(2*r-2*B+g-2*j-2*s,
                              r-j-2*s-2*u)
                for u in range((r-j-2*s)//2+1))
            for t0 in range(s+1):
                h = s+t0
                pp = product(
                    x+k for k in range(1,r+2)
                    if not j+1 <= k <= j+h+1)
                pp *= product(x+j+F(1,2)+k for k in range(h))
                out += v*(-1)**t0*comb(s,t0)*2**(s-t0)*4**h*pp
    return out

quartic = 1+2*x-6*x*x+2*x**3+x**4
checks = 0
for B in range(1,9):
    for r in range(13):
        for gap in (0,2,9):
            p = edge(B,r,gap)
            assert p.degree() <= r
            assert p.get((r,),0) == (quartic**B).get((r,),0)
            checks += 1
    r = 2*B if B%2 else 2*B-1
    p = edge(B,r,2)
    lc = p.get((r,),0)
    assert lc < 0
    bound = 1+sum(abs(v/lc) for mon,v in p.items() if mon[0] < r)
    cutoff = bound.numerator//bound.denominator+1
    assert p.evaluate(x,cutoff) < 0
print("edge polynomial checks", checks)

from collections import defaultdict

def cg(i,j):
    return range(abs(i-j),i+j+1,2)

def multiply(F0,G0):
    out = defaultdict(int)
    for (i,j),v in F0.items():
        for (k,l),w in G0.items():
            for u in cg(i,k):
                for z in cg(j,l):
                    out[u,z] += v*w
    return {key:v for key,v in out.items() if v}

checks = 0
for a,e,B in ((2,2,1),(3,1,1),(7,5,2),
              (6,6,1),(4,3,3),(0,4,4)):
    ff = {(0,0):1}
    for fac,power in (
        ({(1,0):1,(0,1):1},a),
        ({(1,0):1,(0,1):-1},e),
        ({(2,0):1,(0,2):1},B)):
        for _ in range(power):
            ff = multiply(ff,fac)
    degree = a+e+2*B
    for n in range(3,degree+4):
        for m in range(3,n+1):
            if (degree+m-n)%2:
                continue
            d = (degree+m-n)//2
            if d < 0:
                continue
            value = (sum(ff.get((j,0),0) for j in cg(n,m)),
                     ff.get((n,m),0))
            assert value == pair(a,e,B,d,m)
            if d <= m:
                assert value[0] >= abs(value[1])
            checks += 1
print("independent character bridges", checks)

def rows(N0,t0,B,d):
    a = (N0+t0)//2
    e = N0-a
    p = [
        sum((-1)**h*comb(e,h)*comb(a,j-h)
            for h in range(max(0,j-a),min(e,j)+1))
        for j in range(d+2)]
    return [
        [sum(comb(B-h,u)*p[j-2*u]
             for u in range(min(B-h,j//2)+1))
         for j in range(d+2)]
        for h in range(min(B,d)+1)]

def quickpair(rr,B,d,m):
    ans = [0,0]
    for h in range(min(B,d)+1):
        value = pairrow(rr[h],d-h,m)
        for z in (0,1):
            ans[z] += 2**h*comb(B,h)*value[z]
    return ans

count = 0
least = None
zeros = 0

def check(N0,t0,B,d,ms):
    global count,least,zeros
    rr = rows(N0,t0,B,d)
    for m in ms:
        r0,w0 = quickpair(rr,B,d,m)
        v = r0-abs(w0)
        assert v >= 0, (N0,t0,B,d,m,r0,w0)
        count += 1
        zeros += int(v == 0)
        if v > 0:
            least = v if least is None else min(least,v)

for d in range(3,13):
    for B in range(1,7):
        low = max(0,2*d-2*B)
        for N0 in range(low,low+9):
            for t0 in range(N0%2,N0+1,2):
                check(N0,t0,B,d,range(3,d+2))

for d in (16,24,32,48,64):
    for B in (1,2,3,4,6,8):
        for gap in (0,1,2,3,6,15):
            N0 = 2*d-2*B+gap
            ts = {N0%2,N0%2+2,N0%2+4,N0-2,N0-4,N0}
            ts.update(t0 for t0 in (N0//2,N0//3,2*N0//3)
                      if 0 <= t0 <= N0 and (N0-t0)%2 == 0)
            for t0 in sorted(ts):
                if 0 <= t0 <= N0:
                    check(N0,t0,B,d,
                          sorted({3,4,5,d//2,d-2,d-1,d,d+1}))

assert (count,zeros,least) == (39502,0,5)
print("targeted pairs", count, "zeros", zeros, "least slack", least)
print("PASS")
