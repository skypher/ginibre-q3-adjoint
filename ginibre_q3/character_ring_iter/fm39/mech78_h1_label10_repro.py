from fractions import Fraction as Q
from functools import lru_cache
from math import comb
import sympy as S

x,y,u,v,z,t=S.symbols("x y u v z t")

def U(n,x):
    return sum((-1)**j*comb(n-j,j)*x**(n-2*j)
               for j in range(n//2+1))

def bern(P,N,J):
    terms=[(ij,Q(int(c.p),int(c.q))) for ij,c in P.terms()]
    return [[sum(c*Q(comb(i,k),comb(N,k))*Q(comb(j,l),comb(J,l))
                 for (k,l),c in terms if k<=i and l<=j)
             for j in range(J+1)] for i in range(N+1)]

for n,N,L in ((7,32,41),(9,64,53),(10,64,42)):
    div=x+y if n%2 else 1
    P=S.Poly(S.cancel((U(n,x)+U(n,y))/div),x,y)
    F=S.expand(sum(c*2**(i+j)*(2*t-1)**j*z**((i+j)//2)
                   for (i,j),c in P.terms()))
    B=bern(S.Poly(F.subs(z,1),t,z),N,16)
    D=bern(S.Poly(S.diff(F,z),t,z),N,16)
    assert min(row[0] for row in B)>0
    assert min(L*B[i][j]-abs(D[i][j])
               for i in range(N+1) for j in range(17))>=0
    print("Cutoff certificate:",n,4*L+2)

def core(n,sign,div):
    P=S.Poly(S.cancel(
        (U(n,(u+v)/2)+sign*U(n,(u-v)/2))/div),u,v)
    den=int(S.ilcm(*[c.q for ij,c in P.terms()]))
    out=[]
    for (i,j),c in P.terms():
        assert i%2==j%2==0
        out.append((i//2,j//2,int(c*den)))
    return den,out

jobs=[
 ("7 both",166,2,1,4,False,core(7,1,u),695485,4096),
 ("8 plus",130,1,1,4,True,core(8,1,1),168638,5770),
 ("8 minus",126,2,2,5,True,core(8,-1,u*v),145768,4850),
 ("9 both",214,2,1,5,False,core(9,1,u),1521464,8350),
 ("10 plus",170,1,1,5,True,core(10,1,1),384572,10120),
 ("10 minus",198,2,2,6,True,core(10,-1,u*v),594382,10348),
]
cats=[comb(2*j,j)//(j+1) for j in range(225)]

@lru_cache(None)
def base(m,r):
    A,E=2*m,2*r
    N=A+E
    row=[1]
    if N:row.append(A-E)
    for k in range(1,N):
        num=(A-E)*row[k]-(N-k+1)*row[k-1]
        assert num%(k+1)==0
        row.append(num//(k+1))
    return sum(row[2*j]*cats[j]*cats[N//2-j]
               for j in range(N//2+1))

def moment(m,r,b):
    if m>r:m,r=r,m
    return cached(m,r,b)

@lru_cache(None)
def cached(m,r,b):
    if b==0:return base(m,r)
    q=moment(m+1,r,b-1)+moment(m,r+1,b-1)
    assert q%2==0
    return q//2-2*moment(m,r,b-1)

total=0
for name,K,ml,rl,shift,sym,(den,terms),expected,minimum in jobs:
    count=0
    least=None
    for L in range(7+shift,K):
        for b in range(3,L-ml-rl+1):
            N=L-b
            for m in range(ml,N-rl+1):
                r=N-m
                if sym and m<r:continue
                val=sum(c*moment(m+i,r+j,b) for i,j,c in terms)
                assert val%den==0
                val//=den
                assert val>0,(name,m,r,b,val)
                count+=1
                least=val if least is None else min(least,val)
    assert (count,least)==(expected,minimum)
    total+=count
    print(name,count,least)
assert total==3510309
print("Full H=1 closure through label 10: PASS")