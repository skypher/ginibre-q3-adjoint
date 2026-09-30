
import argparse
from math import comb, factorial, prod
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

parser=argparse.ArgumentParser(description="FM-MECH9 exact bounded reproducer")
parser.parse_args()

def det(A):
    A=[row[:] for row in A]; n=len(A)
    if n==0: return 1
    sign=1; previous=1
    for k in range(n-1):
        if A[k][k]==0:
            h=next((h for h in range(k+1,n) if A[h][k]),None)
            if h is None: return 0
            A[k],A[h]=A[h],A[k]; sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=pivot*A[i][j]-A[i][k]*A[k][j]
                assert numerator%previous==0
                A[i][j]=numerator//previous
            A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]

def Z(z,w,q):
    if q==0: return 1
    moments=[sum(a*x**h for x,a in zip(z,w)) for h in range(2*q-1)]
    return det([[moments[i+j] for j in range(q)] for i in range(q)])

def falling(n,l): return factorial(n)//factorial(n-l)

@lru_cache(None)
def coeff(a,e):
    N=a+e; c=[1]
    for k in range(N):
        numerator=(a-e)*c[k]-(N-k+1)*(c[k-1] if k else 0)
        assert numerator%(k+1)==0
        c.append(numerator//(k+1))
    return tuple(c)

def TW(a,e,j,i):
    c=coeff(a,e);N=a+e
    def C(k): return c[k] if 0<=k<=N else 0
    def D(k): return C(k)**2-C(k-1)*C(k+1)
    def B(k): return C(k-1)+C(k+1)
    return D(j)-D(i),B(i)*C(j)-C(i)*B(j)

def ensemble(a,e):
    N=a+e;n=N+2;eps=e%2;q=e//2
    ks=list(range((N+1)//2,N+2))
    x=[2*k-N for k in ks];z=[h*h for h in x]
    mu=[(1 if h==0 else 2)*comb(n,k+1)*h**(2*eps)
        for k,h in zip(ks,x)]
    g={k:comb(n,k+1)*(2*k-N)**eps for k in ks}
    norm=prod(2**n*factorial(2*h+eps)*falling(n,2*h+eps)
              for h in range(q+1))
    assert Z(z,mu,q+1)==norm
    return z,mu,g,norm,F(factorial(e)*2**n,4*falling(n,e+2))

links=0; failures=[]; parameters=0
for N in (20,30,40,60,100,200):
    for e in (3,4,5,6,8,12,16,24):
        a=N-e
        if a<e+2: continue
        z,mu,g,norm,scale=ensemble(a,e);parameters+=1
        q=e//2
        for j in sorted({N//2+h for h in (1,2,4,8)}|{3*N//4,N-3}):
            if not N/2<j<N: continue
            for i in sorted({j+h for h in (2,4,8,16)}|{N,N+1}):
                if i>N+1: continue
                x=(2*j-N)**2;y=(2*i-N)**2
                weights=[w*(s-x)*(s-y) for s,w in zip(z,mu)]
                signed=Z(z,weights,q)
                absolute=Z(z,[abs(w) for w in weights],q)
                factor=scale*(y-x)*g[j]*g[i]/norm
                T,W=TW(a,e,j,i)
                assert W==factor*signed
                assert T>=abs(W)
                if T<factor*absolute: failures.append((a,e,j,i))
                links+=1
print("Norm products:",parameters)
print("Exact Hankel links:",links)
print("Positive-ensemble failures:",failures)

def pmul(P,Q):
    R={}
    for (a,b),v in P.items():
        for (c,d),w in Q.items():
            key=(a+c,b+d);R[key]=R.get(key,0)+v*w
    return {key:value for key,value in R.items() if value}

@lru_cache(None)
def U(k):
    if k==0: return {0:1}
    if k==1: return {1:1}
    out={h+1:v for h,v in U(k-1).items()}
    for h,v in U(k-2).items():out[h]=out.get(h,0)-v
    return {h:v for h,v in out.items() if v}

@lru_cache(None)
def H(k):
    out={}
    for p in range(k+1):
        for u,v in U(p).items():
            for s,t in U(k-p).items():
                out[u,s]=out.get((u,s),0)+v*t
    return {h:v for h,v in out.items() if v}

def direct(r,s,t,a):
    P=pmul(H(s),H(t))
    P=pmul(P,{(h,a-h):comb(a,h) for h in range(a+1)})
    P=pmul(P,{(h,2*r-h):(-1)**h*comb(2*r,h) for h in range(2*r+1)})
    def cat(h):return comb(2*h,h)//(h+1)
    return F(sum(v*cat(u//2)*cat(w//2)
                 for (u,w),v in P.items() if u%2==w%2==0),2)

checks=0
for r in range(2,6):
    for s in range(6):
        for t in range(6):
            for a in range(6):
                value=direct(r,s,t,a)
                if (a+s+t)%2:
                    assert value==0
                else:
                    e=2*r-2;N=a+e;p,q=sorted((s+1,t+1),reverse=True)
                    j=(N+p-q)//2;i=(N+p+q)//2+1
                    T,W=TW(a,e,j,i)
                    assert value==T-W
                checks+=1
print("Independent Catalan links:",checks)

a,e,j,i=8,6,8,13
z,mu,g,norm,scale=ensemble(a,e)
x=(2*j-a-e)**2;y=(2*i-a-e)**2
absolute=Z(z,[w*abs((s-x)*(s-y)) for s,w in zip(z,mu)],e//2)
T,W=TW(a,e,j,i)
bound=scale*(y-x)*g[j]*g[i]*absolute/norm
print("At (8,6,8,13): T, W, unsigned chord =",T,W,bound)
assert T>=bound>=abs(W)

N=14;j=8;i=12;S=(9,11,14)
P=lambda k:prod((2*k-N)**2-(2*h-N)**2 for h in S)
adjacent=[P(k)*P(k+1) for k in range(j,i)]
print("Single-subset obstruction:",adjacent,P(j)*P(i))
assert adjacent==[0]*4 and P(j)*P(i)>0

duality=0
for n in range(8,17):
    for eps in (0,1):
        xs=[x for x in range(n%2,n+1,2) if x or eps==0]
        z=[x*x for x in xs]
        other=eps if n%2==0 else 1-eps
        mu=[(1 if x==0 else 2)*comb(n,(n+x)//2)*x**(2*eps) for x in xs]
        target=[(1 if x==0 else 2)*comb(n,(n+x)//2)*x**(2*other) for x in xs]
        ratios=[]
        for h,x in enumerate(z):
            derivative=prod(x-y for k,y in enumerate(z) if k!=h)
            ratios.append(F(1,mu[h]*derivative**2*target[h]))
        assert len(set(ratios))==1
        duality+=1
print("Particle-hole weight identities:",duality)

def unsigned(a,e,j,i):
    z,mu,g,norm,scale=ensemble(a,e)
    N=a+e;x=(2*j-N)**2;y=(2*i-N)**2
    A=Z(z,[w*abs((s-x)*(s-y)) for s,w in zip(z,mu)],e//2)
    return scale*(y-x)*g[j]*g[i]*A/norm

fold_tests=((8,6,8,13),(8,3,7,11),(9,4,8,12),(6,8,8,12))
for a,e,j,i in fold_tests:
    assert unsigned(a,e,j,i)==unsigned(e,a,j,i)
print("Unsigned fold identities:",len(fold_tests))

z=(0,1,10,19,20)
def minor(u,v):
    indices=[h for h in range(5) if h not in (u,v)]
    return prod(z[t]-z[s] for s,t in combinations(indices,2))
left=minor(1,2)+minor(2,3);right=minor(1,3)
assert (left,right)==(760,2000)
print("Arbitrary-grid relaxation, adjacent versus endpoint:",left,right)
print("All assertions passed.")
