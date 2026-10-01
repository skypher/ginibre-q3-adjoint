import argparse, random, time
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial
import sympy as S

argparse.ArgumentParser(description="FM-MECH68 exact certificates").parse_args()
started = time.monotonic()
t,b,k,c,f,x,y,z = S.symbols("t b k c f x y z")
Q = S.Rational

def choose(v,j):
    return S.sympify(S.prod(v-i for i in range(j)))/factorial(j)

A = {0:S.Integer(1), 1:x}
P = {0:S.Integer(1), 1:x}
for j in range(1,8):
    A[j+1] = S.expand((x*A[j]-(t-j+1)*A[j-1])/(j+1))
for j in range(1,4):
    P[j+1] = S.expand((x*P[j]-(t-2*b-j+1)*P[j-1])/(j+1))

def base_symbolic(d):
    out = 0
    for i in range(d+1):
        for v in range((d-i)//2+1):
            for u in range(d-i-2*v+1):
                r = d-i-u-2*v
                coeff = sum(choose(k+h-2,h)*choose(t-2*i,r-h)
                            for h in range(r+1))
                out += A[2*i]*choose(b,u)*choose(b-u,v)*coeff \
                       *S.binomial(2*(i+u),i+u)/(i+u+1)
    return S.expand(out)

B = {d:base_symbolic(d) for d in (3,4)}

def diagonal(poly,d):
    out = {}
    rem = poly
    for j in range(d,-1,-1):
        out[j] = S.factor(rem.coeff(x,2*j)*factorial(j)**2)
        rem = S.expand(rem-out[j]*P[j]**2)
    assert rem == 0
    return out

a3 = diagonal(B[3],3)
a4 = diagonal(B[4],4)
F3 = S.expand(B[3]+(y*y-c)/2+y*(A[3]+b*x))
N = t+b+k

# Independent enumeration of the distance-three contraction types.
table3 = choose(N+1,3)-t*(N-1)-b \
    +A[2]*(choose(N-2,2)-t+2)+2*(N-5)*A[4]+5*A[6] \
    +(N-3)*choose(b,2)+choose(b,3)+(y*y-c)/2 \
    +(N-4)*A[2]*b+3*A[4]*b+2*A[2]*choose(b,2) \
    +(x*b+A[3])*y
assert S.expand(F3-table3) == 0

sq3 = (P[3]+2*y)**2/4+a3[2]*P[2]**2+a3[1]*x*x \
      +a3[0]-(y*y+c)/2
assert S.expand(F3-sq3) == 0

alpha = x*x+2*b+k-3
Y = A[3]*(N+3*b-5)+x*b*(N+b-5)+4*A[5]
Z = A[4]+b*A[2]+choose(b,2)
F4 = S.expand(B[4]-c+alpha*(y*y-c)/2
              +(z*z-f)/2+y*Y+z*(Z+x*y))
L = 8*b+5*k+t-13
D1 = a4[1]-c/2
D0 = a4[0]-c*(2*b+k-1)/2-f/2
sq4 = (P[4]+2*x*y+Q(5,2)*z)**2/5 \
    +L*(P[3]+2*y)**2/20+P[3]**2/10 \
    +a4[2]*P[2]**2+D1*x*x+D0 \
    -Q(3,10)*x*x*y*y-(6*b+5*k+2*t-11)*y*y/10 \
    -Q(3,4)*z*z-x*y*z-2*b*x*y-b*z
assert S.expand(F4-sq4) == 0

R3 = S.expand(a3[0]-k*(k+1)/2)
R41 = S.expand(a4[1]-b-k-Q(3,10)*k*k)
R40 = S.expand(a4[0]-2*b*k-Q(5,4)*k*k-k/2
               -(16*b+10*k+2*t)*k*k/10)

for d,polys,tails in (
    (3,(a3[2],a3[1],R3),((t,11),(b,6),(k,4))),
    (4,(L,a4[2],R41,R40),((t,27),(b,11),(k,24)))):
    for var,H in tails:
        assert all(min(S.Poly(p.subs(var,var+H),t,b,k).coeffs()) >= 0
                   for p in polys)
    print("symbolic identities and tails:",d,"PASS",flush=True)

def C(n,j):
    return comb(n,j) if 0 <= j <= n else 0

def cat(n):
    return C(2*n,n)//(n+1)

def row(T,X,n=8):
    out = [1,X]
    for j in range(1,n):
        v = X*out[j]-(T-j+1)*out[j-1]
        assert v % (j+1) == 0
        out.append(v//(j+1))
    return out

def K(r,n,T,B):
    return (1,n-1,C(n,2)-T,
            C(n+1,3)-T*(n-1)-B,
            C(n+2,4)-T*C(n,2)-B*(n-1)+C(T,2))[r]

@lru_cache(None)
def mult(i,j):
    return sum((-1)**(j-v)*C(j,v)*cat(i+v)
               for v in range(j+1))

def baseline(d,T,B,K0,aa):
    n = T+B+K0
    return sum(aa[2*i]*C(B,j)*mult(i,j)
               *K(d-i-j,n-2*i-j,T-2*i,B-j)
               for i in range(min(d,T//2)+1)
               for j in range(min(B,d-i)+1))

def data4(T,B,K0,X):
    aa = row(T,X)
    n = T+B+K0
    return (2*baseline(4,T,B,K0,aa),
            X*X+2*B+K0-3,
            aa[3]*(n+3*B-5)+X*B*(n+B-5)+4*aa[5],
            aa[4]+B*aa[2]+C(B,2))

def doubled(data,X,C0,F0,Y0,Z0):
    BB,al,YY,ZZ = data
    return BB+al*Y0*Y0+2*YY*Y0+Z0*Z0 \
           +2*(ZZ+X*Y0)*Z0-(al+2)*C0-F0

n3 = 0
minimum3 = None
for T in range(11):
    for B0 in range(6):
        for K0 in range(4):
            for X in range(-T,T+1,2):
                aa = row(T,X)
                BB = 2*baseline(3,T,B0,K0,aa)
                for C0 in range(K0+1):
                    if T+2*B0+3*C0+4*(K0-C0) < 6:
                        continue
                    for Y0 in range(-C0,C0+1,2):
                        value = BB+Y0*Y0-C0+2*Y0*(aa[3]+B0*X)
                        assert value >= 0
                        n3 += 1
                        minimum3 = value if minimum3 is None \
                                   else min(minimum3,value)
assert n3 == 7869 and minimum3 == 0
print("distance 3 finite profiles:",n3,
      "minimum EVEN:",minimum3,flush=True)

corners = quadratics = 0
minimum4 = None
for T in range(27):
    for B0 in range(11):
        for K0 in range(24):
            # F4(-X,-Y0,Z0) = F4(X,Y0,Z0).
            for X in range(T%2,T+1,2):
                data = data4(T,B0,K0,X)
                BB,al,YY,ZZ = data
                beta = al+2
                if T+2*B0+3*K0 < 8 or beta < 1:
                    for C0 in range(K0+1):
                        for F0 in range(K0-C0+1):
                            if T+2*B0+3*C0+4*F0+5*(K0-C0-F0) < 8:
                                continue
                            for Y0 in range(-C0,C0+1,2):
                                for Z0 in range(-F0,F0+1,2):
                                    value = doubled(data,X,C0,F0,Y0,Z0)
                                    assert value >= 0
                                    corners += 1
                                    minimum4 = value if minimum4 is None \
                                               else min(minimum4,value)
                    continue

                for Y0 in range(-K0,K0+1):
                    M = K0-abs(Y0)
                    DD = ZZ+X*Y0
                    for q in (0,1):
                        if q > M:
                            continue
                        r = (M-q)%2
                        J = (M-q)//2
                        for sg in (-1,1):
                            lin = 2*sg*DD+beta-1
                            constant = BB+al*Y0*Y0+2*YY*Y0 \
                                       -beta*K0+beta*r
                            # Exact minimizer of w^2+lin*w,
                            # w=q+2*j, 0 <= j <= J.
                            j = max(0,min(J,(2-lin-2*q)//4))
                            w = q+2*j
                            value = constant+w*w+lin*w
                            assert value == doubled(
                                data,X,K0-w-r,w,Y0,sg*w)
                            assert value >= 0
                            quadratics += 1
                            minimum4 = value if minimum4 is None \
                                       else min(minimum4,value)
    if T%6 == 0:
        print("distance 4 through t =",T,flush=True)

assert (corners,quadratics,minimum4) == (115,4764410,0)
print("distance 4 small profiles / quadratic minima:",
      corners,quadratics,"minimum EVEN:",minimum4,flush=True)

def character_value(word):
    out = {(0,0):1}
    for tok in word:
        n = abs(tok)
        sg = 1 if tok > 0 else -1
        nxt = defaultdict(int)
        for (i,j),v in out.items():
            for q in range(abs(i-n),i+n+1,2):
                nxt[q,j] += v
            for q in range(abs(j-n),j+n+1,2):
                nxt[i,q] += sg*v
        out = {ij:v for ij,v in nxt.items() if v}
    return out.get((0,0),0)

def formula_value(word):
    rest = list(word)
    ii = max(range(len(rest)),key=lambda i:abs(rest[i]))
    n = abs(rest.pop(ii))
    d = (sum(map(abs,rest))-n)//2
    T = sum(abs(v)==1 for v in rest)+2*rest.count(-2)
    X = rest.count(1)-rest.count(-1)
    B0 = rest.count(2)
    K0 = sum(abs(v)>=3 for v in rest)
    C0 = rest.count(3)+rest.count(-3)
    Y0 = rest.count(3)-rest.count(-3)
    F0 = rest.count(4)+rest.count(-4)
    Z0 = rest.count(4)-rest.count(-4)
    if d == 3:
        aa = row(T,X)
        return 2*baseline(3,T,B0,K0,aa)+Y0*Y0-C0 \
               +2*Y0*(aa[3]+B0*X)
    assert d == 4
    return doubled(data4(T,B0,K0,X),X,C0,F0,Y0,Z0)

bridges = 0
for length in range(1,7):
    for word in combinations_with_replacement(
            (-1,1,-2,2,-3,3,-4,4),length):
        W = sum(map(abs,word))
        if W%2 or sum(v<0 for v in word)%2:
            continue
        d = W//2-max(map(abs,word))
        if d not in (3,4):
            continue
        assert formula_value(word) == character_value(word)
        bridges += 1
assert bridges == 324
print("all signed multisets: labels <= 4, length <= 6, d=3,4:",
      bridges,flush=True)

rng = random.Random(6804)
large = 0
while large < 80:
    d = rng.choice((3,4))
    rest = [rng.choice((-1,1))*rng.randrange(1,18)
            for _ in range(rng.randrange(3,8))]
    n = sum(map(abs,rest))-2*d
    if n < max(map(abs,rest)):
        continue
    sg = (-1)**sum(v<0 for v in rest)
    word = rest+[sg*n]
    assert formula_value(word) == character_value(word)
    large += 1

for word,expected in (
    ((-1,-1,-1,-1,-1,1,1,-3,-3,-5),118),
    ((-1,-1,-1,-1,-2,-3,-4,-5),142)):
    assert formula_value(word) == character_value(word) == expected
    print("atlas bridge:",word,expected,flush=True)

print("large-label bridges:",large,"PASS")
print("ALL CHECKS PASS; seconds =",round(time.monotonic()-started,2))
