import argparse
from functools import lru_cache
from fractions import Fraction as Q
from math import comb

argparse.ArgumentParser(
    description="FM-MECH46: exact insertion obstruction and two-label checks"
).parse_args()

def clean(P):
    return {j:v for j,v in P.items() if v}

def add(*terms):
    R = {}
    for c,P in terms:
        for j,v in P.items():
            R[j] = R.get(j,0) + c*v
    return clean(R)

def mul2(P,Qq):
    R = {}
    for (i,j),v in P.items():
        for (k,l),w in Qq.items():
            z = (i+k,j+l)
            R[z] = R.get(z,0) + v*w
    return clean(R)

GEN = {
    "d": {(1,0):1,(0,1):-1},
    "s": {(1,0):1,(0,1):1},
    "z": {(2,0):1,(0,2):1,(0,0):-2},
}

@lru_cache(None)
def power2(which,n):
    return {(0,0):1} if n == 0 else mul2(power2(which,n-1),GEN[which])

@lru_cache(None)
def background(k,a,b):
    return mul2(mul2(power2("d",k),power2("s",a)),power2("z",b))

@lru_cache(None)
def cat(j):
    return comb(2*j,j)//(j+1)

@lru_cache(None)
def g(k,a,b):
    P = {}
    for (i,j),v in background(k,a,b).items():
        if j % 2 == 0:
            P[i] = P.get(i,0) + v*cat(j//2)
    return clean(P)

@lru_cache(None)
def ballot(i,n):
    if i < n or (i-n) % 2:
        return 0
    j = (i-n)//2
    return comb(i,j) - (comb(i,j-1) if j else 0)

def uc(P):
    return clean({
        n:sum(v*ballot(i,n) for i,v in P.items())
        for n in range(max(P,default=0)+1)
    })

def U(n):
    return {n-2*j:(-1)**j*comb(n-j,j) for j in range(n//2+1)}

def x2plus2(P):
    return add((1,{j+2:v for j,v in P.items()}),(2,P))

def resolvent(P,D):
    assert all(j < D+4 for j in P)
    return {j:Q(v,D+4-j) for j,v in P.items()}

def insertion(P,L):
    return add((1,x2plus2(P)),(-12,resolvent(P,L)))

recurrences = 0
for k in range(1,9):
    for a in range(9):
        for b in range(5):
            if k+a+2*(b+1) > 18:
                continue
            D = k+a+2*b
            prev = g(k,a,b-1) if b else {}
            brace = add((b,x2plus2(prev)),(-(b+3),g(k,a,b)))
            rhs = add((1,x2plus2(g(k,a,b))),(4,resolvent(brace,D)))
            assert rhs == g(k,a,b+1),(k,a,b)
            assert all(v >= 0 for v in uc(rhs).values())
            recurrences += 1
print("IBP recurrences:",recurrences)

for L in (1,2):
    for n in range(L%2,L+1,2):
        assert all(v >= 0 for v in uc(insertion(U(n),L)).values())
assert uc(insertion(U(3),3)) == {1:-1,3:1,5:1}
assert add((Q(5,4),g(2,1,0)),(Q(-1,4),g(3,0,0))) == U(3)
assert uc(g(2,1,1)) == {1:1,3:2,5:1}
assert uc(g(3,0,1)) == {1:9,3:6,5:1}
for L in range(4,21):
    assert uc(insertion(U(L),L))[L-4] == -Q((L-3)*(L+2),4)
print("Forced-map witness:",uc(insertion(U(3),3)))
print("Actual seed images:",uc(g(2,1,1)),uc(g(3,0,1)))
print("Uniform coefficient formula: L=4..20")
print("Distance-3 witness:",uc(insertion(U(12),12))[8])

def mul1(P,Qq):
    R = {}
    for i,v in P.items():
        for j,w in Qq.items():
            R[i+j] = R.get(i+j,0) + v*w
    return clean(R)

@lru_cache(None)
def laurent(k,a,i):
    P = {0:1}
    for count,factor in ((k,{1:1,-1:-1}),
                         (a,{1:1,-1:1}),
                         (i,{2:1,-2:1})):
        for _ in range(count):
            P = mul1(P,factor)
    return P

det_checks = 0
for k in range(1,7):
    for a in range(7):
        for b in range(4):
            if k+a+2*b > 16:
                continue
            C = uc(g(k,a,b))
            for n in range(k+a+2*b+1):
                value = 0
                for i in range(b+1):
                    u = laurent(k,a,i)
                    at = lambda j: u.get(j,0)
                    term = at(n)*(at(n)+at(n+4))
                    term -= at(n+2)*(at(n-2)+at(n+2))
                    value += comb(b,i)*2**(b-i)*term
                assert value == C.get(n,0),(k,a,b,n)
                det_checks += 1
print("Determinant/direct bridges:",det_checks)

@lru_cache(None)
def table(k,a,b):
    C = {}
    for (i,j),v in background(k,a,b).items():
        for n in range(i%2,i+1,2):
            for m in range(j%2,j+1,2):
                C[n,m] = C.get((n,m),0) + v*ballot(i,n)*ballot(j,m)
    return clean(C)

pairs = 0
least = None
for k in range(6):
    for a in range(6):
        for b in range(5):
            C = table(k,a,b)
            for n in range(1,9):
                for m in range(n,9):
                    R = sum(C.get((j,0),0)
                            for j in range(abs(n-m),n+m+1,2))
                    slack = R-abs(C.get((n,m),0))
                    assert slack >= 0,(k,a,b,n,m)
                    if slack and (least is None or slack < least):
                        least = slack
                    pairs += 1
print("Two-label inequalities:",pairs,"least positive slack:",least)

def upoly(n,axis):
    return {((j,0) if axis == 0 else (0,j)):v for j,v in U(n).items()}

def moment(P):
    return sum(v*cat(i//2)*cat(j//2)
               for (i,j),v in P.items() if i%2 == j%2 == 0)

bridges = 0
for k in range(4):
    for a in range(3):
        for b in range(3):
            C = table(k,a,b)
            for n,m in ((1,1),(2,3),(3,4),(4,4)):
                R = sum(C.get((j,0),0)
                        for j in range(abs(n-m),n+m+1,2))
                for eta in (-1,1):
                    eps = eta*(-1)**k
                    Fn = add((1,upoly(n,0)),(eps,upoly(n,1)))
                    Fm = add((1,upoly(m,0)),(eta,upoly(m,1)))
                    actual = Q(moment(mul2(background(k,a,b),mul2(Fn,Fm))),2)
                    assert actual == R+eta*C.get((n,m),0)
                    bridges += 1
print("Two-label/direct bridges:",bridges)
print("PASS")
