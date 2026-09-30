
import argparse, random, sympy as s
from math import comb
from functools import lru_cache

ap = argparse.ArgumentParser()
ap.add_argument("--rows", type=int, default=300)
args = ap.parse_args()

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def mul(c, p):
    q = [0] * (len(c) + len(p) - 1)
    for i, a in enumerate(c):
        for j, b in enumerate(p):
            q[i+j] += a*b
    return q

def D(c, k):
    return at(c,k)**2 - at(c,k-1)*at(c,k+1)

def T(c, j, k):
    return at(c,j)*at(c,k) - at(c,j-1)*at(c,k+1)

def W(c, x, C):
    return sum(D(c,k) for k in range(x,x+C+1)) - T(c,x,x+C)

def L(c, x, C):
    y = x+C
    return (sum(T(c,k-1,k) for k in range(x,y+1))
            - at(c,x)*at(c,y-1) + at(c,x-2)*at(c,y+1))

# Symbolic identity (1).
for C in range(1,9):
    c = s.symbols("c:"+str(C+3))
    rhs = 0
    for i in range((C+1)//2):
        R = [0]*(C+1)
        R[i], R[C-i] = 1, -1
        rhs += D(mul(c,R), C+1)
    assert s.expand(W(c,1,C)-rhs) == 0
print("multiplier identities: 8")

# Unique covariance in (2).
c = s.symbols("c:6")
pairs = [(i,j) for i in range(4) for j in range(i,4)]
unknowns = s.symbols("G:"+str(len(pairs)))
G = s.zeros(4)
for (i,j), h in zip(pairs,unknowns):
    G[i,j] = G[j,i] = h
expr = sum(G[i,j]*(at(c,4-i)*at(c,4-j)
                   - at(c,3-i)*at(c,5-j))
           for i in range(4) for j in range(4))
sol = list(s.linsolve(s.Poly(expr-W(c,1,3),*c).coeffs(),
                     unknowns))
assert len(sol) == 1
G = G.subs(dict(zip(unknowns,sol[0])))
J = s.zeros(4)
for i in range(4):
    J[i,3-i] = 1
assert G == s.eye(4)-J
print("C=3 unique covariance:", G.tolist())
print("forced mean quadratic discriminant:",
      G[1,1]+2*G[0,1]-3*G[0,0])

# Bounded vetting of the unproved induction inequality.
rng = random.Random(2312)
checks = 0
for trial in range(args.rows):
    c = [1]
    for _ in range(rng.randrange(0,13)):
        c = mul(c, [rng.choice([-5,-3,-2,-1,1,2,3,5]),
                    rng.choice([1,2,3])])
    for C in range(1,9):
        for x in range(-C-1,len(c)+1):
            A, B, E = W(c,x,C), L(c,x,C), W(c,x-1,C)
            assert A >= 0 and E >= 0 and B*B <= 4*A*E
            for rho in (-2,0,3):
                assert W(mul(c,[1,rho]),x,C) == A+rho*B+rho*rho*E
            checks += 1
print("linear-factor identities and PSD tests:", checks, "passed")

# Definition-level Catalan evaluation.
def pmul(p, q):
    out = {}
    for (i,j), b in p.items():
        for (k,l), d in q.items():
            out[i+k,j+l] = out.get((i+k,j+l),0)+b*d
    return {k:v for k,v in out.items() if v}

def hp(k):
    return {(k-2*d-j,j): (-1)**d*comb(k+1-d,d)
            for d in range(k//2+1) for j in range(k+1-2*d)}

@lru_cache(None)
def cat(k):
    return 0 if k % 2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r, i, j):
    return sum((-1)**d*comb(2*r,d)*cat(i+2*r-d)*cat(j+d)
               for d in range(2*r+1))

def phi(r, p):
    twice = sum(v*moment(r,i,j) for (i,j),v in p.items())
    assert twice % 2 == 0
    return twice//2

p = {(0,0):1}
for k in (21,4,2):
    p = pmul(p,hp(k))
p = pmul(p,{(15-j,j):comb(15,j) for j in range(16)})
A = phi(3,p)
B = phi(3,pmul(p,{(1,1):1}))
C = phi(3,pmul(p,{(2,0):1,(0,2):1,(0,0):-4}))
assert (A,B,C) == (1843,9097,16813)
print("direct Catalan (A,B,C):", (A,B,C))
print("4A, |B|:", 4*A, abs(B))
print("endpoint numerators:", 4*A-2*B+C, 4*A+2*B+C)
print("4AC-B^2:", 4*A*C-B*B)

# Local-inequality obstruction.
z = s.symbols("z")
c = [4,-3,1,1,-2,-4]
p = s.Poly(sum(v*z**k for k,v in enumerate(c)),z)
print("real-root count:", p.count_roots(-s.oo,s.oo),
      "degree:", p.degree())
print("cubic slacks:",
      [4*D(c,k)*D(c,k+1)-T(c,k,k+1)**2 for k in range(5)])
print("W_3(1):", W(c,1,3))
