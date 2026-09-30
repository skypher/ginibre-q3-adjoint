# Exact Catalan-moment evaluator of phi_r (FM-MECH28, astra_max_ceres, printed verifier; helper part only).
# phi(r, f) = (1/2) E[(x-y)^(2r) f], f a polynomial dict {(i,j): coeff}; h(k) = H_k, hat(k) = hat S_k, word(hs, ss).
import argparse, random
from math import comb
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement as cwr


def add(f, g):
    h = dict(f)
    for ij, v in g.items():
        h[ij] = h.get(ij, 0) + v
    return {ij:v for ij,v in h.items() if v}

def scale(f, n):
    return {ij:n*v for ij,v in f.items() if n*v}

def mul(f, g):
    h = {}
    for (i,j),v in f.items():
        for (k,l),w in g.items():
            h[i+k,j+l] = h.get((i+k,j+l),0) + v*w
    return {ij:v for ij,v in h.items() if v}

one = {(0,0):1}
s = {(1,0):1, (0,1):1}
xy = {(1,1):1}

@lru_cache(None)
def h(k):
    return {(k-2*b-j,j):(-1)**b*comb(k+1-b,b)
            for b in range(k//2+1)
            for j in range(k+1-2*b)}

@lru_cache(None)
def hat(k):
    out = {}
    for b in range(k//2+1):
        j = k-2*b
        v = (-1)**b*comb(k-b,b)
        out[j,0] = out.get((j,0),0) + v
        out[0,j] = out.get((0,j),0) + v
    return out

@lru_cache(None)
def cat(k):
    return 0 if k < 0 or k % 2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def moment(r, i, j):
    return sum((-1)**b*comb(2*r,b)*cat(i+2*r-b)*cat(j+b)
               for b in range(2*r+1))

def phi(r, f):
    return Q(sum(v*moment(r,i,j) for (i,j),v in f.items()),2)

def power(f, n):
    z = one
    for _ in range(n):
        z = mul(z,f)
    return z

def jets(f):
    out = [f]
    for n in range(1,max((i+j for i,j in f),default=0)+1):
        z = {}
        for (i,j),v in out[-1].items():
            if i:
                z[i-1,j] = z.get((i-1,j),0) + i*v
            if j:
                z[i,j-1] = z.get((i,j-1),0) + j*v
        assert all(v % n == 0 for v in z.values())
        out.append({ij:v//n for ij,v in z.items() if v})
    return out

@lru_cache(None)
def tilted(r, a, i, j):
    return sum(comb(a,b)*moment(r,i+a-b,j+b)
               for b in range(a+1))

def profile(jet, r, a):
    # Coefficients of 2P; positive scaling preserves all tested signs.
    return [sum(v*tilted(r,a,i,j) for (i,j),v in f.items())
            for f in jet]

def hurwitz(cs):
    while cs and cs[0] == 0:
        cs = cs[1:]
    while cs and cs[-1] == 0:
        cs = cs[:-1]
    if cs and cs[-1] < 0:
        return ("leading",cs[-1])
    if len(cs) <= 1:
        return None
    aa = cs[::-1]
    n = len(aa)-1
    H = [[aa[2*j-i+1] if 0 <= 2*j-i+1 <= n else 0
          for j in range(n)] for i in range(n)]
    prev = 1
    for k in range(n):
        pivot = H[k][k]
        if pivot <= 0:
            return (k+1,pivot)
        if k == n-1:
            break
        for i in range(k+1,n):
            for j in range(k+1,n):
                v = pivot*H[i][j] - H[i][k]*H[k][j]
                assert v % prev == 0
                H[i][j] = v//prev
        prev = pivot
    return None

def word(hs, ss):
    w = one
    for k in hs:
        w = mul(w,h(k))
    for k in ss:
        w = mul(w,hat(k))
    return w


