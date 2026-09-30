# FM-SEC71 (luna_max_mars) printed screen: Phi(C,p,q,s) - Phi(3,p,0,s) >= 0 on 22,680 e = 1 cases.
from math import comb
from functools import lru_cache

@lru_cache(None)
def c_row(a, k):
    if k < 0 or k > a+1:
        return 0
    return ((comb(a,k) if 0 <= k <= a else 0)
            -(comb(a,k-1) if 0 <= k-1 <= a else 0))

def g(a, k):
    return c_row(a,k), c_row(a,k-1)

def wedge(u, v):
    return u[0]*v[1]-u[1]*v[0]

def P(a, x, C):
    return (sum(c_row(a,k)**2-c_row(a,k-1)*c_row(a,k+1)
                for k in range(x,x+C+1))
            -c_row(a,x)*c_row(a,x+C)
            +c_row(a,x-1)*c_row(a,x+C+1))

def phi(C, p, q, s):
    N = C+p+q+2*s
    a = N-1
    alpha, tau = p+s, N-s+1
    ga = g(a,alpha)
    gb = g(a,q+s)
    gc = g(a,tau)
    gd = g(a,s-C-1)
    U = ga[0]+gb[0], ga[1]+gb[1]
    V = gc[0]+gd[0], gc[1]+gd[1]
    return P(a,alpha,C)-P(a,s-C-1,C)-wedge(U,V)

count = 0
minimum = None
argmin = None
for C in range(3,13):
    for p in range(16):
        for q in range(min(p,8)+1):
            for s in range(21):
                difference = phi(C,p,q,s)-phi(3,p,0,s)
                count += 1
                if minimum is None or difference < minimum:
                    minimum, argmin = difference, (C,p,q,s)

assert count == 22680
assert minimum == 0 and argmin == (3,0,0,0)
print("count, minimum, argmin =", count, minimum, argmin)
