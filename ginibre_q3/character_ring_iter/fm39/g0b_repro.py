
import argparse
from fractions import Fraction

ap = argparse.ArgumentParser()
ap.add_argument("--rmax", type=int, default=12)
ap.add_argument("--amax", type=int, default=50)
args = ap.parse_args()

def row(a, e):
    N, d = a+e, a-e
    c, previous = [1], 0
    for k in range(N):
        num = d*c[-1] - (N-k+1)*previous
        assert num % (k+1) == 0
        previous, cnext = c[-1], num//(k+1)
        c.append(cnext)
    assert c[-1] == (-1)**e
    return c

def half(v):
    return 0 if v[1] > 0 or (v[1] == 0 and v[0] > 0) else 1

def lt(v, w):
    if half(v) != half(w):
        return half(v) < half(w)
    return v[0]*w[1] - v[1]*w[0] > 0

def data(c):
    N = len(c)-1
    at = lambda k: c[k] if 0 <= k <= N else 0
    D = [at(k)**2-at(k-1)*at(k+1) for k in range(N+1)]
    g = [(at(k),at(k-1)) for k in range(N+2)]
    E, turns = [0], [0]
    for k in range(N+1):
        assert D[k] > 0
        E.append(E[-1]+D[k])
        turns.append(turns[-1]+lt(g[k+1],g[k]))
    return at, D, g, E, turns

checks = long_checks = strip_long = 0
first = {}
for r in range(2,args.rmax+1):
    e = 2*r-3
    for a in range(args.amax+1):
        N, d = a+e, a-e
        at,D,g,E,turns = data(row(a,e))
        assert D == D[::-1]
        assert all(D[k] >= D[k+1]
                   for k in range((N+1)//2,N))
        for x in range(N+1):
            for y in range(x+3,N+1):
                C, S = y-x, E[y+1]-E[x]
                T = at(x)*at(y)-at(x-1)*at(y+1)
                wraps = turns[y+1]-turns[x]
                long = (wraps >= 2 or
                        (wraps == 1 and lt(g[x],g[y+1])))
                if long and T > 0:
                    if C in (3,4) and C not in first:
                        first[C] = (r,a,e,x,C,S,T)
                    if abs(a-e) == 2:
                        strip_long += 1
                if d*d > 2*(N+1):
                    continue
                product = D[x]*D[y]
                assert S*S >= 4*C*product
                # T² <= (6+4√2) product, exactly.
                b = T*T-6*product
                assert b <= 0 or b*b <= 32*product*product
                assert S > abs(T)
                checks += 1
                long_checks += long and T > 0

print("first long-sweep C=3 and C=4:", first)
print("proved-band windows checked:", checks)
print("of these, long-sweep T>0:", long_checks)
print("a=e+-2 long-sweep T>0 windows:", strip_long)

r,a,x,C = 150,307,300,4
at,D,g,E,turns = data(row(a,2*r-3))
y = x+C
S = E[y+1]-E[x]
T = at(x)*at(y)-at(x-1)*at(y+1)
wraps = turns[y+1]-turns[x]
assert T > 0 and (wraps >= 2 or
                  (wraps == 1 and lt(g[x],g[y+1])))
assert (a-2*r+2)**2 <= 8*r-9
b = at(301)**2
print("r=150 example (a,parts):", a, (304,304,3))
print("(S,T,S-T)/c_301^2:",
      Fraction(S,b), Fraction(T,b), Fraction(S-T,b))
