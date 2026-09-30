
from math import comb
from fractions import Fraction
import sympy as sp

# C=1 cubic-discriminant identity.
a,b,c,d,t=sp.symbols('a b c d t', real=True)
Q=a+3*b*t+3*c*t**2+d*t**3
Dx=b*b-a*c
Dy=c*c-b*d
T=b*c-a*d
assert sp.factor(sp.discriminant(Q,t)-27*(4*Dx*Dy-T*T)) == 0
print('C=1 symbolic discriminant identity: PASS')

# Root-sum/Cauchy identities underlying the general coefficient lemma.
def poly_mul(p,q):
    out=[Fraction(0)]*(len(p)+len(q)-1)
    for i,u in enumerate(p):
        for j,v in enumerate(q):
            out[i+j]+=u*v
    return out

for n in range(2,9):
    roots=[Fraction((-1)**i*(i+2),i+1) for i in range(n)]
    q=[Fraction(1)]
    for r0 in roots:
        q=poly_mul(q,[-r0,Fraction(1)])
    aa=[q[j]/comb(n,j) for j in range(n+1)]
    mu=sum(roots,Fraction(0))/n
    nu=sum((1/r0 for r0 in roots),Fraction(0))/n
    s1=sum(((1/roots[i]-1/roots[j])**2
            for i in range(n) for j in range(i+1,n)),Fraction(0))
    s2=sum(((roots[i]-roots[j])**2
            for i in range(n) for j in range(i+1,n)),Fraction(0))
    cross=sum(((roots[i]-roots[j])**2/(roots[i]*roots[j])
               for i in range(n) for j in range(i+1,n)),Fraction(0))
    assert aa[1]**2-aa[0]*aa[2] == aa[0]**2*s1/(n*n*(n-1))
    assert aa[n-1]**2-aa[n-2]*aa[n] == aa[n]**2*s2/(n*n*(n-1))
    assert aa[1]*aa[n-1]-aa[0]*aa[n] == aa[0]*aa[n]*(mu*nu-1)
    assert mu*nu-1 == cross/(n*n)
    u=[(roots[i]-roots[j])/(roots[i]*roots[j])
       for i in range(n) for j in range(i+1,n)]
    v=[roots[i]-roots[j] for i in range(n) for j in range(i+1,n)]
    assert sum((u[i]*v[j]-u[j]*v[i])**2
               for i in range(len(u)) for j in range(i+1,len(u))) == (
        sum(x*x for x in u)*sum(x*x for x in v)
        - sum(u[i]*v[i] for i in range(len(u)))**2)
print('root-sum formulas and exact Cauchy-Gram identity, degrees 2..8: PASS')

# Exact binomial row formulas and their stronger ratio bound.
def at(v,k):
    return v[k] if 0 <= k < len(v) else 0

for N in range(31):
    cb=[comb(N,k) for k in range(N+1)]
    for C in range(N+1):
        for x in range(N-C+1):
            y=x+C
            D0=at(cb,x)**2-at(cb,x-1)*at(cb,x+1)
            D1=at(cb,y)**2-at(cb,y-1)*at(cb,y+1)
            T0=at(cb,x)*at(cb,y)-at(cb,x-1)*at(cb,y+1)
            assert D0 == Fraction(comb(N,x)**2*(N+1),
                                  (x+1)*(N-x+1))
            assert D1 == Fraction(comb(N,y)**2*(N+1),
                                  (y+1)*(N-y+1))
            assert T0 == Fraction(comb(N,x)*comb(N,y)*(N+1)*(C+1),
                                  (N-x+1)*(y+1))
            assert (N-x+1)*(y+1)-(x+1)*(N-y+1) == C*(N+2)
            assert (x+1)*(N-y+1) <= (N-x+1)*(y+1)
print('binomial formulas and stronger ratio bound, 0 <= N <= 30: PASS')

# Exact base-row grid: all a=0..39, odd e=1..39, and every support window.
def row(a0,e0):
    N0=a0+e0
    return [
        sum(((-1)**i)*comb(e0,i)*comb(a0,k-i)
            for i in range(e0+1) if 0 <= k-i <= a0)
        for k in range(N0+1)
    ]

def Drow(v,k):
    return at(v,k)**2-at(v,k-1)*at(v,k+1)

def Trow(v,x,C):
    return at(v,x)*at(v,x+C)-at(v,x-1)*at(v,x+C+1)

keys=('all','C1','C2','edge','center','a0')
counts={k:0 for k in keys}
zeros={k:0 for k in keys}
maxratio={}
maxwitness={}

for a0 in range(40):
    for e0 in range(1,40,2):
        cv=row(a0,e0)
        N0=a0+e0
        assert cv[0]==1 and cv[-1]==(-1)**e0
        for C0 in range(N0+1):
            for x0 in range(N0-C0+1):
                y0=x0+C0
                dx=Drow(cv,x0)
                dy=Drow(cv,y0)
                tv=Trow(cv,x0,C0)
                slack=(C0+1)**2*dx*dy-tv*tv
                assert dx>=0 and dy>=0 and slack>=0, (
                    a0,e0,x0,C0,dx,dy,tv,slack)
                cats=['all']
                if C0==1: cats.append('C1')
                if C0==2: cats.append('C2')
                if x0==0 or y0==N0: cats.append('edge')
                if 2*x0+C0==N0: cats.append('center')
                if a0==0: cats.append('a0')
                for key in cats:
                    counts[key]+=1
                    if slack==0: zeros[key]+=1
                    if dx*dy:
                        ratio=Fraction(tv*tv,(C0+1)**2*dx*dy)
                        if key not in maxratio or ratio>maxratio[key]:
                            maxratio[key]=ratio
                            maxwitness[key]=(a0,e0,x0,C0)
                if x0==0:
                    assert dx==1 and tv==cv[C0]
                if y0==N0:
                    assert dy==cv[N0]**2 and tv==cv[x0]*cv[N0]
                if 2*x0+C0==N0:
                    assert dx==dy
                    assert tv==at(cv,x0-1)**2-at(cv,x0)**2
                if a0==0:
                    cb=[comb(e0,k) for k in range(N0+1)]
                    assert Drow(cv,x0)==Drow(cb,x0)
                    assert Drow(cv,y0)==Drow(cb,y0)
                    assert abs(Trow(cv,x0,C0))==abs(Trow(cb,x0,C0))

print('base-row grid a=0..39, e=1,3,..,39: PASS')
print('window counts:',counts)
print('zero-slack counts:',zeros)
for key in keys:
    if key in maxratio:
        print(key,'max squared ratio=',maxratio[key],
              'witness (a,e,x,C)=',maxwitness[key],
              'decimal=',float(maxratio[key]))
print('support-edge, center, and a=0 twist identities: PASS')
