
import argparse
from fractions import Fraction as Q
from math import comb
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement as cwr
import sympy as sy
argparse.ArgumentParser(
    description="Exact verifier for the additive quartic cutoff."
).parse_args()

def add(*terms):
    out=defaultdict(int)
    for scale,f in terms:
        for ij,c in f.items(): out[ij]+=scale*c
    return {ij:c for ij,c in out.items() if c}

def mul(f,g):
    out=defaultdict(int)
    for (i,j),c in f.items():
        for (k,l),d in g.items(): out[i+k,j+l]+=c*d
    return {ij:c for ij,c in out.items() if c}

one={(0,0):1}
X={(1,0):1}
A={(0,0):2,(1,0):1,(0,1):-1}
# X=(x+y)^2/4, Y=(x-y)^2/4.
# E_j=h_(2j), O_j=h_(2j+1)/(x+y).
E=[one]; O=[one]
for j in range(1,9):
    ep=E[j-2] if j>=2 else {}
    op=O[j-2] if j>=2 else {}
    E.append(add((4,mul(X,O[-1])),(-1,mul(A,E[-1])),
                 (4,mul(X,op)),(-1,ep)))
    O.append(add((1,E[-1]),(-1,mul(A,O[-1])),
                 (1,E[-2]),(-1,op)))

def U(k,y):
    if k<0: return Q(0)
    a,b=Q(1),y
    for _ in range(k): a,b=b,y*b-a
    return a

def ell(tag,k):
    if tag=="S":
        return Q(2*k*k) if k%2==0 else ell("H",k-1)
    if k%2: return Q((k-1)**2*(k+2)*(k+3),6)
    return Q(2*k*k*(k+1)*(k+3),3)

def star(tag,k):
    j=k//2
    if tag=="H": return O[j] if k%2 else E[j]
    if k%2: return {(v,u):c for (u,v),c in E[j].items()}
    return add((2,E[j]),(2,E[j-1]),(-4,mul(X,O[j-1])))

def at(f,x,y):
    return sum(c*x**i*y**j for (i,j),c in f.items())

@lru_cache(None)
def raw(tag,k):
    f=defaultdict(int)
    if tag=="H":
        for b in range((k+1)//2+1):
            d=k+1-2*b
            for j in range(d):
                f[d-1-j,j]+=(-1)**b*comb(k+1-b,b)
    else:
        for b in range(k//2+1):
            d=k-2*b
            c=(-1)**b*comb(k-b,b)
            f[d,0]+=c
            f[0,d]+=c
    return dict(f)

@lru_cache(None)
def cat(n):
    return 0 if n%2 else comb(n,n//2)//(n//2+1)

@lru_cache(None)
def kernel(r,i,j):
    return sum((-1)**b*comb(2*r,b)*cat(i+2*r-b)*cat(j+b)
               for b in range(2*r+1))

@lru_cache(None)
def scaled(r,m,i,j):
    # M(r+j,m+i)/(4**(i+j)*M(r,m)).
    value=Q(1)
    for _ in range(j):
        value*=Q((2*r+1)*(2*r+3),(r+m+2)*(r+m+3))
        r+=1
    for _ in range(i):
        value*=Q((2*m+1)*(2*m+3),(r+m+2)*(r+m+3))
        m+=1
    return value

def moment(r,m):
    return Q(1,2)*4**(r+m)*scaled(0,0,m,r)

def ratio(word,r,a):
    t=sum(k%2 for tag,k in word)
    assert (a+t)%2==0
    f=one
    for tag,k in word: f=mul(f,star(tag,k))
    return sum(c*scaled(r,(a+t)//2,i,j)
               for (i,j),c in f.items())

def direct(word,r,a):
    f={(a-j,j):comb(a,j) for j in range(a+1)}
    for tag,k in word: f=mul(f,raw(tag,k))
    return Q(sum(c*kernel(r,i,j) for (i,j),c in f.items()),2)

nb=nr=0
for tag,labels in (("H",range(1,17)),("S",range(2,17))):
    for k in labels:
        j=k//2
        for iq in range(-8,9):
            q=Q(iq,8)
            b=at(star(tag,k),(1+q)**2,(1-q)**2)
            if tag=="H" and k%2:
                expected=sum(U(l,2*q)**2 for l in range(j+1))
            elif tag=="S" and k%2==0:
                expected=k+1+U(k,2*q)
            else:
                y=-2*q if tag=="S" else 2*q
                expected=sum((U(l,y)+U(l-1,y))**2
                             for l in range(j+1))
            lower=Q(k+1,2) if tag=="S" and k%2==0 else Q(j+1,4)
            assert b==expected and b>=max(Q(1),lower)
            nb+=1
            for iz in range(9):
                z=Q(iz,8)
                f=at(star(tag,k),z*(1+q)**2,z*(1-q)**2)
                assert abs(f/b-1)<=ell(tag,k)*(1-z)
                nr+=1
print("Boundary identities/bounds:",nb,"; relative-error checks:",nr)

nt=0
specs=[
    ([("H",k) for k in range(1,7)],5,6,2),
    ([(tag,k) for tag in ("H","S") for k in range(2,6)],4,4,1)
]
for alphabet,rstop,astop,number in specs:
    bridges=0
    for count in range(1,4):
        for word in cwr(alphabet,count):
            t=sum(k%2 for tag,k in word)
            for r in range(1,rstop):
                for a in range(t%2,astop,2):
                    assert direct(word,r,a)==(
                        moment(r,(a+t)//2)*ratio(word,r,a))
                    bridges+=1
            L=sum((ell(tag,k) for tag,k in word),Q(0))
            for offset in range(number):
                a=t%2+2*offset
                bound=(7*L-a-t-6)/2
                r=max(1,(count+1)//2,
                      -((-bound.numerator)//bound.denominator))
                assert 2*r+a+t+6>=7*L
                assert ratio(word,r,a)>=Q(1,25)
                nt+=1
    print("Direct Catalan bridges:",bridges)
print("Quantitative cutoff checks:",nt)

word=(("H",3),("H",3),("H",2))
print("New cutoff for (3,3,2):",7*sum(ell(*z) for z in word))
for n in range(4,8):
    word=(("H",2*n+1),("H",2*n+1),("H",2*n))
    value=ratio(word,n,2*n)
    assert value>0
    print("Uncovered-family normalized value:",n,value)

n,v=sy.symbols("n v")
L=8*n**2*(n+1)*(2*n+3)
deficits=[
    (2*n+3)*(2*(2*n+1)*(2*n+5)+2*n*(2*n+4))-5*(4*n+4),
    (2*n+5)*(6*(2*n+1)*(2*n+5)+5*(2*n+1)*(2*n+3))-15*(4*n+3),
    7*L-(4*n+8)
]
for f in deficits:
    coefficients=sy.Poly(f.subs(n,v+1),v).all_coeffs()
    assert all(c>0 for c in coefficients)
    print("Coverage deficit:",sy.expand(f),
          "; at n=1+v:",coefficients)
print("All exact checks passed.")
