import argparse, sympy as S
from itertools import product
argparse.ArgumentParser(description="Exact spin-block obstruction").parse_args()

I=S.eye(3)
J=[S.Matrix(3,3,lambda r,c:-S.I*S.LeviCivita(k,r,c))
   for k in range(3)]
K=S.zeros(27)
for j in J:
    T=sum((S.kronecker_product(
        *[j if i==h else I for i in range(3)])
        for h in range(3)),S.zeros(27))
    K+=T*T

Id=S.eye(27)
eigen=[0,2,6,12]
Ps=[]
for v in eigen:
    P=Id
    for w in eigen:
        if v!=w: P=P*(K-w*Id)/(v-w)
    Ps.append(P)

digits=list(product(range(3),repeat=3))
def transpose_leg(A,h):
    B=S.zeros(27)
    for r,rr in enumerate(digits):
        for c,cc in enumerate(digits):
            aa=list(rr); bb=list(cc)
            aa[h],bb[h]=bb[h],aa[h]
            B[r,c]=A[9*aa[0]+3*aa[1]+aa[2],
                     9*bb[0]+3*bb[1]+bb[2]]
    return B
def project(A):
    for h in range(3): A=(A+transpose_leg(A,h))/2
    return A

M0=project(Ps[0])
M2=project(Ps[1])/3
assert ((M0*M0).trace(),(M0*M2).trace(),(M2*M2).trace())==(
    S.Rational(1,2),-S.Rational(1,4),S.Rational(13,20))
blocks=[(P*(M0*M0+M0*M2)).trace() for P in Ps]
assert blocks==[S.Rational(1,8),-S.Rational(1,16),
                S.Rational(3,16),0]
print(blocks)
