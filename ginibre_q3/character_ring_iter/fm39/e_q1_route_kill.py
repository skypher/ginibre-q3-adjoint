
import argparse
from fractions import Fraction as Q
from math import comb
argparse.ArgumentParser(description='Exact recurrence-state quadratic-form check.').parse_args()

def f(N,d,j,x,y):
    z=Q(d*y-(N-j+1)*x,j+1)
    w=Q(d*z-(N-j)*y,j+2)
    v=Q(d*w-(N-j-1)*z,j+3)
    Dj=y*y-x*z; Di=w*w-z*v
    W=y*(z+v)-w*(x+z)
    return Dj-Di-W

N,d,j=18,-12,17
A=f(N,d,j,Q(1),Q(0)); C=f(N,d,j,Q(0),Q(1))
B=f(N,d,j,Q(1),Q(1))-A-C
Delta=A*C-B*B/4
print('L-W matrix coefficients A,B,C',A,B,C,'det',Delta)

c=[0]*(N+1)
for h in range(16):
    for k in range(4):
        c[h+k]+=(-1)**h*comb(15,h)*comb(3,k)
x,y=Q(c[j-1]),Q(c[j])
assert (x,y)==(-63,12) and f(N,d,j,x,y)==93
print('actual (c[j-1],c[j])',x,y,'L-W',f(N,d,j,x,y))
