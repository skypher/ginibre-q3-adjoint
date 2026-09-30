
import argparse
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, prod

p = argparse.ArgumentParser(
    description="FM-MECH10: exact local splitting and deletion checks")
p.add_argument("--max-n", type=int, default=22,
               help="largest N=a+e in the bounded screen")
args = p.parse_args()

def det(A):
    A = [row[:] for row in A]
    m = len(A)
    if not m:
        return 1
    sign, old = 1, 1
    for k in range(m-1):
        if not A[k][k]:
            h = next((h for h in range(k+1, m) if A[h][k]), None)
            if h is None:
                return 0
            A[k], A[h] = A[h], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1, m):
            for j in range(k+1, m):
                num = pivot*A[i][j] - A[i][k]*A[k][j]
                assert num % old == 0
                A[i][j] = num // old
            A[i][k] = 0
        old = pivot
    return sign*A[-1][-1]

def H(z, mu, q):
    if q == 0:
        return 1
    moments = [
        sum(w*x**h for x, w in zip(z, mu))
        for h in range(2*q-1)
    ]
    return det([[moments[i+j] for j in range(q)] for i in range(q)])

def table(z, mu, g, q):
    return {
        (u,v): (z[v]-z[u])*g[u]*g[v] *
        H(z, [w*abs((x-z[u])*(x-z[v]))
              for x,w in zip(z,mu)], q)
        for u,v in combinations(range(len(z)), 2)
    }

def slack(A, u, v):
    return sum(A[k,k+1] for k in range(u,v)) - A[u,v]

def physical(a, e):
    N, n, eps = a+e, a+e+2, e % 2
    ks = [k for k in range((N+1)//2, N+2)
          if eps == 0 or 2*k != N]
    X = [2*k-N for k in ks]
    b = [comb(n,k+1) for k in ks]
    z = [x*x for x in X]
    mu = [(1 if x == 0 else 2)*v*x**(2*eps)
          for x,v in zip(X,b)]
    g = [v*x**eps for x,v in zip(X,b)]
    q = e//2
    scale = F(factorial(a)*factorial(e)*2**n,
              4*factorial(n)*H(z,mu,q+1))
    return ks, X, z, mu, g, q, scale

def coeffs(a,e):
    return [
        sum((-1)**h*comb(e,h)*comb(a,k-h)
            for h in range(e+1) if 0 <= k-h <= a)
        for k in range(a+e+1)
    ]

def delta(z, R):
    return prod(z[v]-z[u] for u,v in combinations(R,2))

adj_checks = local_checks = grouped_checks = 0
for N in range(8, args.max_n+1):
    for e in range(3, min(10,(N-2)//2)+1):
        a = N-e
        ks,X,z,mu,g,q,scale = physical(a,e)
        U = table(z,mu,g,q)
        c = coeffs(a,e)
        C = lambda k: c[k] if 0 <= k < len(c) else 0
        D = lambda k: C(k)**2-C(k-1)*C(k+1)

        for h in range(len(z)-1):
            assert scale*U[h,h+1] == D(ks[h])-D(ks[h+1])
            adj_checks += 1

        for h in range(1,len(z)-1):
            u,v,w = h-1,h,h+1
            if X[u] == 0:
                continue
            rem = [s for s in range(len(z))
                   if s not in (u,v,w)]
            d1,d2 = z[v]-z[u],z[w]-z[v]
            total = 0
            for size in (q,q-1):
                for R in combinations(rem,size):
                    fu,fv,fw = [
                        prod(abs(z[x]-z[s]) for s in R)
                        for x in (u,v,w)
                    ]
                    base = prod(mu[s] for s in R)*delta(z,R)**2
                    if size == q:
                        term = base*(
                            d1*g[u]*g[v]*fu*fv
                            + d2*g[v]*g[w]*fv*fw
                            - (d1+d2)*g[u]*g[w]*fu*fw)
                    else:
                        term = (
                            base*d1*d2*(d1+d2)*fu*fv*fw *
                            (g[u]*g[v]*mu[w]*fw
                             + g[v]*g[w]*mu[u]*fu
                             - g[u]*g[w]*mu[v]*fv))
                    assert term >= 0
                    total += term
                    grouped_checks += 1
            assert total == U[u,v]+U[v,w]-U[u,w]
            local_checks += 1

print("adjacent coefficient links:", adj_checks)
print("local identities / nonnegative grouped terms:",
      local_checks, grouped_checks)

# Negative occupied branch on the physical lattice.
a,e,j,h,i = 17,6,12,14,16
ks,X,z,mu,g,q,scale = physical(a,e)
U = table(z,mu,g,q)
J,Hh,I = [ks.index(k) for k in (j,h,i)]
keep = [s for s in range(len(z)) if s != Hh]
zz = [z[s] for s in keep]
E = table(zz, [mu[s] for s in keep],
          [g[s] for s in keep], q)
O = table(
    zz, [mu[s]*(z[s]-z[Hh])**2 for s in keep],
    [g[s]*abs(z[s]-z[Hh]) for s in keep], q-1)
empty = scale*slack(E,J,I-1)
occupied = scale*mu[Hh]*slack(O,J,I-1)
local = scale*(U[Hh-1,Hh]+U[Hh,Hh+1]-U[Hh-1,Hh+1])
whole = scale*slack(U,J,I)
assert whole == empty+occupied+local
assert occupied < 0 < whole
print("physical deletion (empty, occupied, local, whole):")
print(tuple(str(x) for x in (empty,occupied,local,whole)))

# Local concavity alone does not imply the full inequality.
z = list(range(5))
mu, g = [1,4,6,4,1], [1]*5
U = table(z,mu,g,1)
dual = [
    F(g[k], mu[k]*prod(abs(z[k]-x) for x in z if x != z[k]))
    for k in range(5)
]
assert dual == [F(1,24)]*5
print("flat primal/dual weights: local slacks, path, chord:",
      [slack(U,k,k+2) for k in range(3)],
      sum(U[k,k+1] for k in range(4)), U[0,4])

# Independent Catalan-moment check.
def Ucoeff(p):
    return {p-2*h: (-1)**h*comb(p-h,h)
            for h in range(p//2+1)}

def moment(d):
    return 0 if d % 2 else comb(d,d//2)//(d//2+1)

def phi_two(r,s,t,a):
    e = 2*r-2
    total = 0
    for h in range(e+1):
        for k in range(a+1):
            wt = (-1)**h*comb(e,h)*comb(a,k)
            bx,by = e-h+a-k,h+k
            for p,up in Ucoeff(s+1).items():
                for q,uq in Ucoeff(t+1).items():
                    total += wt*up*uq*(
                        moment(bx+p+q)*moment(by)
                        + moment(bx)*moment(by+p+q)
                        - moment(bx+p)*moment(by+q)
                        - moment(bx+q)*moment(by+p))
    assert total % 2 == 0
    return total//2

c = coeffs(a,e)
C = lambda k: c[k] if 0 <= k < len(c) else 0
D = lambda k: C(k)**2-C(k-1)*C(k+1)
B = lambda k: C(k-1)+C(k+1)
T = D(j)-D(i)
W = B(i)*C(j)-C(i)*B(j)
value = phi_two(4,3,2,17)
assert value == T-W
print("Catalan phi_4(h_3 h_2 h_1^17), T, W:", value,T,W)
