
from math import comb
from itertools import product

# Exact bivariate polynomials: {(x-degree,y-degree): integer coefficient}.
def add(P,Q):
    R=P.copy()
    for m,c in Q.items():
        R[m]=R.get(m,0)+c
        if R[m]==0: del R[m]
    return R

def scale(P,c): return {m:c*v for m,v in P.items() if c*v}
def mul(P,Q):
    R={}
    for (i,j),a in P.items():
        for (k,l),b in Q.items():
            m=(i+k,j+l); R[m]=R.get(m,0)+a*b
    return {m:c for m,c in R.items() if c}
def powp(P,n):
    R={(0,0):1}
    for _ in range(n): R=mul(R,P)
    return R

X={(1,0):1}; Y={(0,1):1}
U=[{(0,0):1},{(1,0):1}]
for n in range(2,40): U.append(add(mul(X,U[-1]),scale(U[-2],-1)))
def up(n): return {(i,0):c for (i,_),c in U[n].items()}
def uy(n): return {(0,i):c for (i,_),c in U[n].items()}
def cat(n): return comb(2*n,n)//(n+1)
def moment(n): return 0 if n%2 else cat(n//2)
def expectation(P): return sum(c*moment(i)*moment(j) for (i,j),c in P.items())
def coeffs(a,e):
    c=[0]*(a+e+1)
    for i in range(a+1):
        for j in range(e+1): c[i+j]+=comb(a,i)*comb(e,j)*((-1)**j)
    return c
def cget(c,k): return c[k] if 0<=k<len(c) else 0
def fusion(p,q): return range(abs(p-q),p+q+1,2)
def Kpoly(a,e): return mul(powp(add(X,scale(Y,-1)),e),powp(add(X,Y),a))
def F(p,q,a,e):
    N=a+e
    if (N+p+q)%2: return 0
    c=coeffs(a,e); i=(N+p+q)//2+1; j=(N+p-q)//2
    return cget(c,i-1)*cget(c,j)-cget(c,i)*cget(c,j+1)
def Wformula(p,q,a,e):
    N=a+e
    if (N+p+q)%2: return 0
    c=coeffs(a,e); i=(N+p+q)//2+1; j=(N+p-q)//2
    return (cget(c,i-1)+cget(c,i+1))*cget(c,j)-cget(c,i)*(cget(c,j-1)+cget(c,j+1))
def D(c,k): return cget(c,k)**2-cget(c,k-1)*cget(c,k+1)
def Delta(p,q,a,e):
    N=a+e
    if (N+p+q)%2: return 0
    c=coeffs(a,e); j=(N+abs(p-q))//2; i=(N+p+q)//2+1
    return D(c,j)-D(c,i)
def C2(x,y,z,t,a,e):
    return sum(F(p,abs(z-t),a,e)-F(p,z+t+2,a,e) for p in fusion(x,y))
def C3(x,y,z,q,a,e):
    return sum(F(q,abs(x-v),a,e)-F(q,x+v+2,a,e) for v in fusion(y,z))
def hpoly(k):
    R={}
    for i in range(k+1): R=add(R,mul(up(i),uy(k-i)))
    return R
def phi_direct(parts,a,r):
    P={(0,0):1}
    for k in parts: P=mul(P,hpoly(k))
    num=expectation(mul(mul(powp(add(X,scale(Y,-1)),2*r),powp(add(X,Y),a)),P))
    assert num%2==0
    return num//2
def split_direct(parts,a,r):
    As=[k+1 for k in parts]; e=2*r-4; K=Kpoly(a,e)
    P={(0,0):1}
    for A in As: P=mul(P,up(A))
    main=expectation(mul(K,P)); cross31=0
    for i in range(4):
        P={(0,0):1}
        for j,A in enumerate(As):
            if j!=i: P=mul(P,up(A))
        cross31+=expectation(mul(K,mul(P,uy(As[i]))))
    cross22=0
    for (i,j),(k,l) in (((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))):
        P=mul(up(As[i]),up(As[j])); Q=mul(uy(As[k]),uy(As[l]))
        cross22+=expectation(mul(K,mul(P,Q)))
    return main-cross31+cross22
def closed4_terms(parts,a,r):
    A,B,C,Dd=sorted((k+1 for k in parts),reverse=True); e=2*r-4
    main=sum(Delta(p,q,a,e) for p in fusion(A,B) for q in fusion(C,Dd))
    corr31=(C3(B,C,Dd,A,a,e)+C3(A,C,Dd,B,a,e)+
            C3(A,B,Dd,C,a,e)+C3(A,B,C,Dd,a,e))
    corr22=(C2(A,B,C,Dd,a,e)+C2(A,C,B,Dd,a,e)+C2(A,Dd,B,C,a,e))
    return main,corr31,corr22

wchecks=support=boundary=0
for a,e in product(range(7),range(0,7,2)):
    N=a+e; c=coeffs(a,e); K=Kpoly(a,e)
    for p,q in product(range(8),repeat=2):
        wd=expectation(mul(K,mul(up(p),uy(q))))
        assert wd==Wformula(p,q,a,e)
        assert wd==F(p,q,a,e)-F(p,q+2,a,e)
        assert wd==expectation(mul(K,mul(up(q),uy(p))))
        wchecks+=1
        if p+q>N: assert wd==0; support+=1
        if p+q==N: assert wd==cget(c,q)==cget(c,p); boundary+=1

splitcases=closedcases=0
for r,a in product(range(2,5),range(4)):
    for parts in product(range(5),repeat=4):
        value=phi_direct(parts,a,r)
        assert split_direct(parts,a,r)==value; splitcases+=1
        terms=closed4_terms(parts,a,r)
        assert terms[0]-terms[1]+terms[2]==value; closedcases+=1

thresholdcases=0
for N in range(9):
    for A in range(1,13):
        for B in range(1,A+1):
            for C in range(1,B+1):
                for Dd in range(1,C+1):
                    mu1=Dd+max(0,A-B-C); mu2=A-B+C-Dd
                    if mu1>N and mu2>N:
                        for sing,(x,y,z) in ((A,(B,C,Dd)),(B,(A,C,Dd)),
                                             (C,(A,B,Dd)),(Dd,(A,B,C))):
                            for p in fusion(x,y):
                                for q in fusion(p,z): assert q+sing>N
                        for x,y,z,t in ((A,B,C,Dd),(A,C,B,Dd),(A,Dd,B,C)):
                            for p in fusion(x,y):
                                for q in fusion(z,t): assert p+q>N
                        thresholdcases+=1

print('W coefficient/symmetry checks:',wchecks)
print('W=F(p,q)-F(p,q+2) checks:',wchecks)
print('support-zero checks p+q>N:',support)
print('boundary checks W(p,q)=c_q=c_p:',boundary)
print('four-factor split vs direct Catalan checks:',splitcases)
print('endpoint closed form vs direct Catalan checks:',closedcases)
print('strict-threshold fusion tuples checked:',thresholdcases)
for parts,a,r in (((0,0,0,0),0,2),((5,1,1,1),0,2),((2,0,0,0),2,2)):
    m,x,y=closed4_terms(parts,a,r)
    print('example',parts,'a=',a,'r=',r,'direct=',phi_direct(parts,a,r),
          'main=',m,'R31=',x,'R22=',y)
print('all exact assertions passed')
