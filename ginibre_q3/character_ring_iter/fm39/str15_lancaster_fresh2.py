from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb

def expand(fs):
    F = {(0, 0): 1}
    for eps, n in fs:
        G = {}
        for (a, b), c in F.items():
            for r in range(abs(a-n), a+n+1, 2):
                G[r, b] = G.get((r, b), 0) + c
            for s in range(abs(b-n), b+n+1, 2):
                G[a, s] = G.get((a, s), 0) + eps*c
        F = {k: v for k, v in G.items() if v}
    return F

def upolys(N):
    U = [[1], [0, 1]]  # U_0=1, U_1=x; U_{n+1}=xU_n-U_{n-1}
    for n in range(1, N):
        z = [0]*(n+2)
        for j, c in enumerate(U[n]):
            z[j+1] += c
        for j, c in enumerate(U[n-1]):
            z[j] -= c
        U.append(z)
    return U

def qpoly(F):
    N = max((r for r, s in F if s == 0), default=0)
    U, q = upolys(N), [0]*(N+1)
    for r in range(N+1):
        for j, c in enumerate(U[r]):
            q[j] += F.get((r, 0), 0)*c
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q

def bernstein(p, left, right):
    d = len(p)-1
    a = [Fraction(0)]*(d+1)
    for j, c in enumerate(p):
        for k in range(j+1):
            a[k] += c*comb(j, k)*Fraction(left)**(j-k)*Fraction(right-left)**k
    return [sum(a[k]*Fraction(comb(i, k), comb(d, k))
                for k in range(i+1)) for i in range(d+1)]

def bisect(B):
    w, L, R = list(B), [B[0]], [B[-1]]
    while len(w) > 1:
        w = [(w[i]+w[i+1])/2 for i in range(len(w)-1)]
        L.append(w[0])
        R.append(w[-1])
    return L, list(reversed(R))

def certify(B, depth=0, cap=18):
    if all(v >= 0 for v in B):
        return True, depth
    if B[0] < 0 or B[-1] < 0:
        return False, depth
    if depth == cap:
        return None, depth
    L, R = bisect(B)
    a, b = certify(L, depth+1, cap), certify(R, depth+1, cap)
    if False in (a[0], b[0]):
        return False, max(a[1], b[1])
    if None in (a[0], b[0]):
        return None, max(a[1], b[1])
    return True, max(a[1], b[1])

def evalp(p, x):
    v = Fraction(0)
    for c in reversed(p):
        v = v*x+c
    return v

cases = {
    "1": [(1,1)]*2 + [(1,2)]*3 + [(-1,3)]*14,
    "2": [(-1,n) for n in range(1,9)],
    "3": [(1,1)]*2 + [(-1,2)] + [(1,3)]*4
         + [(-1,4)]*3 + [(1,5)]*4 + [(1,8)],
    "4": [(-1,1)]*6 + [(-1,2)]*3 + [(-1,3)]*3
         + [(-1,4), (-1,7)],
}

for name, fs in cases.items():
    F = expand(fs)
    d = max((r for r, s in F if r == s), default=0)
    D = [F.get((r, r), 0) for r in range(d+1)]
    q = qpoly(F)
    shifts = []
    for i, (eps, n) in enumerate(fs):
        g = list(fs)
        g[i] = (eps, n+2)
        shifts.append(expand(g).get((0,0), 0) - F.get((0,0), 0))
    print(name, "Phi", F.get((0,0),0),
          "D degree", d, "min Bernstein", min(bernstein(D,-1,1)),
          "D certificate", certify(bernstein(D,-1,1)),
          "q certificate", certify(bernstein(q,-2,2)),
          "minimum label shift", min(shifts))
    if name == "4":
        print("q(-29/16) =", evalp(q, Fraction(-29,16)))

tokens = [(eps,n) for eps in (1,-1) for n in (1,2,3)]
small = direct = 0
for k in (2,4,6):
    for fs in combinations_with_replacement(tokens, k):
        if sum(eps < 0 for eps, n in fs) % 2:
            continue
        F = expand(fs)
        d = max((r for r,s in F if r == s), default=0)
        D = [F.get((r,r),0) for r in range(d+1)]
        ok, depth = certify(bernstein(D,-1,1))
        assert ok is True
        small += 1
        direct += (depth == 0)
print("small D checks", small, "direct", direct, "after subdivision", small-direct)

q = qpoly(expand([(-1,1),(-1,2)]))
print("small q witness", evalp(q, Fraction(-15,8)))
