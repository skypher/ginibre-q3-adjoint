
import argparse
from math import comb
from fractions import Fraction
from functools import lru_cache
from collections import defaultdict
from itertools import combinations_with_replacement as triples

argparse.ArgumentParser(
    description="Exact adjacent-strip bounds and Catalan bridges."
).parse_args()

def b(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def q4(x, y, z, t):
    return (x+t)*(y+z)-x*t-y*z

@lru_cache(None)
def endpoint(n, mode, parts):
    N = 2*n + int(mode != 0)
    def c(k):
        if not 0 <= k <= N or (mode == 0 and k % 2):
            return 0
        return (-1)**(k//2)*b(n, k//2)*(mode if k % 2 else 1)
    def D(k):
        return c(k)**2-c(k-1)*c(k+1)
    A, B, C = sorted((v+1 for v in parts), reverse=True)
    assert (N+A-B-C) % 2 == 0
    al = (N+A-B-C)//2
    be, ga, de = al+C, al+B, al+B+C
    M = sum(D(k) for k in range(al, be+1))
    M -= sum(D(k) for k in range(ga+1, de+2))
    R = q4(c(al), c(be), c(ga), c(de+2))
    R -= q4(c(al-1), c(be+1), c(ga+1), c(de+1))
    return M, R, al

def diagonal_upper(n, parts):
    P, Q, S = sorted(((v+1)//2 for v in parts), reverse=True)
    al = n+P-Q-S
    h = al//2
    x, t = b(n,h), b(n,h+Q+S+1)
    y = b(n,h+S+al%2)
    z = b(n,h+Q+al%2)
    I, J, K = x*y-z*t, x*z-y*t, y*z-x*t
    return I+J+(K if al%2 else -K)

@lru_cache(None)
def even_data(n, P, Q, S):
    M, R, al = endpoint(n, 1 if n%2 else -1, (2*S,2*Q,2*P))
    if al%2:
        assert R == 0
        return M, R, None, None
    h = al//2
    x, X = b(n,h), b(n,h-1)
    y, Y = b(n,h+S), b(n,h+S+1)
    z, Z = b(n,h+Q), b(n,h+Q+1)
    T, t = b(n,h+Q+S+1), b(n,h+Q+S+2)
    I = x*y+z*t-X*Y-Z*T
    J = x*z+y*t-X*Z-Y*T
    K = x*t+y*z-X*T-Y*Z
    assert I >= K >= 0 and J >= K
    U = I+J-K
    assert R <= U
    return M, R, U, (h,x,X,y,Y,z,Z,T,t)

diag_count = even_count = inc_count = 0
for n in range(25):
    for P in range(1,15):
        for Q in range(1,P+1):
            for S in range(1,Q+1):
                parts = (2*S-1,2*Q-1,2*P-1)
                M0, R0, _ = endpoint(n,0,parts)
                assert R0 <= diagonal_upper(n,parts) <= M0
                diag_count += 1
                M, R, U, data = even_data(n,P,Q,S)
                assert M >= R
                even_count += 1
                if U is None:
                    continue
                assert M >= U
                if P > Q:
                    m1, _, u1, _ = even_data(n,P,Q+1,S+1)
                    h,x,X,y,Y,z,Z,T,t = data
                    L, f = b(n,h-2), b(n,h+Q+S+3)
                    delta = (X-t)*(2*(X+t)+L+x+T+f)
                    delta += (L-X+t-T)*(Y+Z)
                    delta -= (X-x+f-t)*(y+z)
                    assert m1-u1-(M-U) == delta >= 0
                    inc_count += 1
print("Bounds:", diag_count, even_count, "increments:", inc_count)

transfer = 0
for n in range(0,25,2):
    for parts in triples(range(1,15),3):
        if sum(v%2 for v in parts) != 2:
            continue
        i = next(i for i,v in enumerate(parts) if v%2 == 0)
        lo, hi = list(parts), list(parts)
        lo[i] -= 1
        hi[i] += 1
        lo, hi = tuple(sorted(lo)), tuple(sorted(hi))
        M,R,_ = endpoint(n,-1,parts)
        m0,_,_ = endpoint(n,0,lo)
        m1,_,_ = endpoint(n,0,hi)
        U = diagonal_upper(n,lo)+diagonal_upper(n,hi)
        assert M == m0+m1 and R <= U <= M
        transfer += 1
print("Lower transfers:", transfer)

def mul(f,g):
    out = defaultdict(int)
    for (i,j),x in f.items():
        for (k,l),y in g.items():
            out[i+k,j+l] += x*y
    return {p:v for p,v in out.items() if v}

@lru_cache(None)
def hp(k):
    return {(k-2*d-j,j): (-1)**d*comb(k+1-d,d)
            for d in range(k//2+1) for j in range(k+1-2*d)}

@lru_cache(None)
def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r,i,j):
    return sum((-1)**d*comb(2*r,d)*cat(i+2*r-d)*cat(j+d)
               for d in range(2*r+1))

@lru_cache(None)
def direct(r,a,parts):
    f = {(a-j,j): comb(a,j) for j in range(a+1)}
    for k in parts:
        f = mul(f,hp(k))
    return Fraction(sum(v*moment(r,i,j)
                        for (i,j),v in f.items()),2)

bridge = upper_transfer = 0
for r in range(2,7):
    e = 2*r-3
    for a in (e-1,e+1):
        n, mode = min(a,e), (1 if a>e else -1)
        for parts in triples(range(1,9),3):
            ans = direct(r,a,parts)
            if sum(parts)%2:
                assert ans == 0
            else:
                M,R,_ = endpoint(n,mode,parts)
                assert ans == M-R >= 0
                if a == e+1 and sum(v%2 for v in parts) == 2:
                    i = next(i for i,v in enumerate(parts) if v%2 == 0)
                    lo,hi = list(parts),list(parts)
                    lo[i] -= 1
                    hi[i] += 1
                    assert ans == (
                        direct(r,e,tuple(sorted(lo))) +
                        direct(r,e,tuple(sorted(hi))))
                    upper_transfer += 1
            bridge += 1
print("Catalan bridges:", bridge, "upper transfers:", upper_transfer)
