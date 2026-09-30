import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

argparse.ArgumentParser(
    description="FM-MECH47: exact proof certificate for the {1,2,3} consumer sector"
).parse_args()

def coarse(L,N):
    K=Q(N,3)+L+1
    return Q(481,50)*Q(2,3)**(N+L//2-2)*K*(K+1)

def bbound(L,N,b):
    K=Q(N,3)+L+1
    low=Q(factorial(N),2**(N+1))
    for j in range(1,N+2):
        low/=b+j
    high=(Q(2,3)**(N+1)-Q(2,5)**(N+1))/Q(N+1)*Q(5,9)**b
    return 9*2**L*(K+b)*(K+b+1)*(low+high)

large_L=Q(481,50)*Q(32,81)**10*Q(80,3)*Q(83,3)
assert large_L<1
assert Q(2,3)*Q(23,22)**2<1
large_N=max((coarse(L,20),L) for L in range(2,22))
assert large_N[0]<1
b_values=[(bbound(L,N,42),L,N)
          for L in range(2,22)
          for N in range(max(2,(L+1)//2),20)]
large_b=max(b_values)
assert large_b[0]<1
print("Large-L bound:",large_L)
print("Large-N worst:",large_N)
print("Large-b checks:",len(b_values),"worst parameters:",large_b[1:])
print("Large-b worst bound:",large_b[0])

@lru_cache(None)
def cat(j):
    return comb(2*j,j)//(j+1)

@lru_cache(None)
def row(A,E):
    n=A+E
    out=[1]
    if n:
        out.append(A-E)
    for j in range(1,n):
        num=(A-E)*out[j]-(n-j+1)*out[j-1]
        assert num%(j+1)==0
        out.append(num//(j+1))
    return tuple(out)

@lru_cache(None)
def M0(m,r):
    c=row(2*m,2*r)
    return sum(c[j]*cat(j//2)*cat(m+r-j//2)
               for j in range(0,2*(m+r)+1,2))

@lru_cache(None)
def M(m,r,b):
    if m>r:
        return M(r,m,b)
    if b==0:
        return M0(m,r)
    top=M(m+1,r,b-1)+M(m,r+1,b-1)
    assert top%2==0
    return top//2-2*M(m,r,b-1)

def mixed_table(m,r,J,B):
    out=[[M(m,r,b) for b in range(B+1)]]
    for j in range(1,J+1):
        rr=[]
        for b in range(B-j+1):
            num=8*(m-r)*out[j-1][b]
            if j>=2:
                num+=4*(j-1)*out[j-2][b+1]+8*(j-1)*out[j-2][b]
            if b:
                num+=12*b*rr[-1]
            den=2*(m+r+b+j+2)
            assert num%den==0,(m,r,j,b)
            rr.append(num//den)
        out.append(rr)
    return out

def shifted_R(m,r,j,b):
    num=sum((-1)**(j-h)*comb(j,h)*M(m+h,r+j-h,b)
            for h in range(j+1))
    assert num%4**j==0
    return num//4**j

coeffs={(L,al):tuple((j,v) for j,v in enumerate(row(al,L-al)) if v)
        for L in range(2,22) for al in range(L+1)}

count=zeros=checks=0
least=None
for N in range(2,20):
    for m in range(N+1):
        r=N-m
        J=min(21,2*N)
        T=mixed_table(m,r,J,62)
        for j in (0,1,J):
            for bb in (0,62-j):
                assert T[j][bb]==shifted_R(m,r,j,bb)
                checks+=1
        for L in range(2,J+1):
            for al in range(L+1):
                cc=coeffs[L,al]
                for b in range(42):
                    value=sum(v*T[j][b+L-j] for j,v in cc)
                    assert value>=0,(2*m,2*r,b,al,L-al,value)
                    count+=1
                    zeros+=value==0
                    if value and (least is None or value<least[0]):
                        least=(value,2*m,2*r,b,al,L-al)
    print("Finite remainder N =",N,"cumulative entries =",count)
assert count==1848168
print("Finite remainder:",count,"zeros:",zeros,"least:",least)
print("IBP / shifted-moment bridges:",checks)

def add(*terms):
    out={}
    for c,P in terms:
        for key,v in P.items():
            out[key]=out.get(key,0)+c*v
    return {key:v for key,v in out.items() if v}

def mul(P,R):
    out={}
    for (i,j),v in P.items():
        for (k,l),w in R.items():
            key=(i+k,j+l)
            out[key]=out.get(key,0)+v*w
    return {key:v for key,v in out.items() if v}

def U(n,axis):
    return {((n-2*j,0) if axis==0 else (0,n-2*j)):
            (-1)**j*comb(n-j,j) for j in range(n//2+1)}

ONE={(0,0):1}
S={(1,0):1,(0,1):1}
D={(1,0):1,(0,1):-1}
Z={(2,0):1,(0,2):1,(0,0):-2}
P={(1,1):1}
H=add((1,U(2,0)),(1,U(2,1)),(1,mul(U(1,0),U(1,1))))
S3=add((1,U(3,0)),(1,U(3,1)))
assert H==add((1,Z),(1,P))
assert S3==mul(S,add((1,Z),(-1,P)))
cross=add((1,mul(Z,Z)),(-1,mul(P,P)))
assert cross==add((1,ONE),(1,U(4,0)),(1,U(4,1)),
                  (1,mul(U(2,0),U(2,1))))
raw=mul(mul(D,D),mul(H,H))
even={key:v for key,v in raw.items() if key[1]%2==0}
boundary={}
for j in range(4):
    boundary=add((1,boundary),(1,U(2*j,0)),(1,U(2*j,1)))
assert even==boundary
print("Character boundary identities: PASS")

GEN={"s":S,"d":D,"z":Z,"h":H,"S3":S3}
@lru_cache(None)
def power(which,n):
    return ONE if n==0 else mul(power(which,n-1),GEN[which])

def moment(poly):
    return sum(v*cat(i//2)*cat(j//2)
               for (i,j),v in poly.items() if i%2==j%2==0)

bridges=0
for r in range(1,4):
    for a in range(4):
        for b in range(3):
            for al in range(min(3,2*r)+1):
                for ga in range(4):
                    poly=mul(mul(power("d",2*r),power("s",a)),
                             mul(power("z",b),
                                 mul(power("h",al),power("S3",ga))))
                    actual=moment(poly)
                    A=a+ga
                    if A%2:
                        expected=0
                    else:
                        L=al+ga
                        expected=sum(v*shifted_R(A//2,r,j,b+L-j)
                                     for j,v in enumerate(row(al,ga)) if v)
                    assert actual==expected,(r,a,b,al,ga)
                    assert actual>=0 and actual%2==0
                    bridges+=1
print("Direct consumer / Catalan bridges:",bridges)
print("PASS")
