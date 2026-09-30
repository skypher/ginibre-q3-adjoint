
import argparse
from math import comb
from functools import lru_cache
import sympy as sp
argparse.ArgumentParser(description='Exact checks and q=1 factor verification for (E).').parse_args()

@lru_cache(None)
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for r in range(e+1):
        for s in range(a+1): c[r+s]+=(-1)**r*comb(e,r)*comb(a,s)
    return c

def at(c,k): return c[k] if 0<=k<len(c) else 0

def e_data(a,e,j,i):
    c=row(a,e)
    D=lambda k: at(c,k)**2-at(c,k-1)*at(c,k+1)
    B=lambda k: at(c,k-1)+at(c,k+1)
    W=B(i)*at(c,j)-at(c,i)*B(j)
    return D(j)-D(i),W,D(j)-D(i)-abs(W)

# q=1 finite screen on the open coefficient-row region.
least=None; count=fail=0
for a in range(3,41):
    for e in range(3,41):
        if abs(a-e)<=1: continue
        N=a+e
        for j in range((N+1)//2,N):
            i=j+2
            L,W,slack=e_data(a,e,j,i); count+=1
            if least is None or slack<least[0]: least=(slack,(a,e,j,i,N,L,W))
            if slack<0: fail+=1
assert count==28314 and fail==0 and least[0]==4
print('open q=1 screen',count,'failures',fail,'least slack',least)

print('least-slack exact',e_data(3,5,7,9))
print('cross-root q=1 exact',(5,3,5,7),e_data(5,3,5,7))
N=8
X=lambda k: 2*k-N
K3=lambda x: x*(x*x-28)
assert (X(5),X(7))==(2,6) and K3(2)==-48 and K3(6)==48

# Center-to-edge boundary j=N/2, i=N+1, for even N.
center_count=center_fail=0; center_min=None
for a in range(3,41):
    for e in range(3,41):
        if abs(a-e)<=1: continue
        N=a+e
        if N%2: continue
        L,W,slack=e_data(a,e,N//2,N+1); center_count+=1
        center_fail+=slack<0
        if center_min is None or slack<center_min[0]:
            center_min=(slack,(a,e,N//2,N+1,L,W))
assert center_count==684 and center_fail==0
print('center-to-edge boundary',center_count,'failures',center_fail,'least slack',center_min)

# q=1 along |a-e|=2 and 3, e=3..200.
for d in (2,3):
    total=bad=0; least_d=None
    for e in range(3,201):
        a=e+d; N=a+e
        for j in range((N+1)//2,N):
            L,W,slack=e_data(a,e,j,j+2); total+=1
            if slack<0: bad+=1
            if least_d is None or slack<least_d[0]:
                least_d=(slack,(a,e,j,j+2,L,W))
    assert total==20295 and bad==0
    print('folded q=1 ray',d,'pairs',total,'failures',bad,'least slack',least_d)

# Outer endpoint j=N-1, i=N+1, within its consumer range N>=2.
endpoint_count=endpoint_bad=0
for a in range(41):
    for e in range(41):
        N=a+e
        if N<2: continue
        endpoint_count+=1
        endpoint_bad+=e_data(a,e,N-1,N+1)[2]<0
assert endpoint_count==1678 and endpoint_bad==0
print('outer endpoint j=N-1,i=N+1',endpoint_count,'failures',endpoint_bad)

# Symbolically verify the q=1 factorizations.
e,t,b=sp.symbols('e t b', integer=True, positive=True)
bm1=t*b/(e-t+1)
bm2=t*(t-1)*b/((e-t+1)*(e-t+2))
bp1=(e-t)*b/(t+1)
bp2=(e-t)*(e-t-1)*b/((t+1)*(t+2))
def signed_slacks(seq):
    cm1,c0,cp1,cp2,cp3=seq
    D0=c0**2-cm1*cp1; D2=cp2**2-cp1*cp3
    W=(cp1+cp3)*c0-cp2*(cm1+cp1)
    return sp.factor((D2-D0-W)/b**2),sp.factor((D2-D0+W)/b**2)

d2_even=(-2*bm1,b-bm1,2*b,b-bp1,-2*bp1)
d2_odd=(b-bm1,2*b,b-bp1,-2*bp1,bp2-bp1)
d2_even_expected=(e-2*t)*(e+1)**2*(e+2)/((t+1)**2*(e-t+1)**2)
d2_odd_minus=(e+1)*(e+2)*(e+5)*(e-2*t-1)/((t+1)**2*(t+2)*(e-t+1))
d2_odd_plus=(e+1)**2*(e+2)*(e-2*t-1)/((t+1)**2*(t+2)*(e-t+1))
assert sp.simplify(signed_slacks(d2_even)[0]-d2_even_expected)==0
assert sp.simplify(signed_slacks(d2_even)[1]-d2_even_expected)==0
assert sp.simplify(signed_slacks(d2_odd)[0]-d2_odd_minus)==0
assert sp.simplify(signed_slacks(d2_odd)[1]-d2_odd_plus)==0

d3_even=(-3*bm1+bm2,b-3*bm1,3*b-bm1,3*b-bp1,b-3*bp1)
d3_odd=(b-3*bm1,3*b-bm1,3*b-bp1,b-3*bp1,bp2-3*bp1)
den_even=(t+1)**2*(e-t+2)*(e-t+1)**2
den_odd=(t+1)**2*(t+2)*(e-t+1)**2
d3_even_minus=(e+1)*(e+2)*(e-2*t+1)*(e**2+3*e+8*t+10)/den_even
d3_even_plus=(e-2*t)*(e+1)**2*(e+2)*(e+3)/den_even
d3_odd_minus=(e+1)*(e+2)*(e-2*t-1)*(e**2+11*e-8*t+10)/den_odd
d3_odd_plus=(e-2*t)*(e+1)**2*(e+2)*(e+3)/den_odd
assert sp.simplify(signed_slacks(d3_even)[0]-d3_even_minus)==0
assert sp.simplify(signed_slacks(d3_even)[1]-d3_even_plus)==0
assert sp.simplify(signed_slacks(d3_odd)[0]-d3_odd_minus)==0
assert sp.simplify(signed_slacks(d3_odd)[1]-d3_odd_plus)==0
print('symbolic q=1 factorizations verified for d=2,3')
