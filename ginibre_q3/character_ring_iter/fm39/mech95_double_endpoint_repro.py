from collections import defaultdict
from fractions import Fraction as Q
from math import comb
from random import Random

def ddiag(n,j):
    a=[0]*(n+1)
    for l in range(min(j,n-j)+1):
        for v in range(l+1):
            a[n-2*l+2*v]+=(-1)**(l+v)*comb(j,l)*comb(n-j,l)*comb(l,v)
    return a

def multiply(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out

cnt=0
for n in range(1,25):
    for j in range(n+1):
        a=ddiag(n,j)
        assert sum(Q(2*x,k+2) for k,x in enumerate(multiply(a,a)))==Q(1,n+1)
        cnt+=1
print("Exact rotation L2 norms:",cnt,"PASS",flush=True)

def tensor(word):
    T={(0,0):1}
    for n,e in word:
        R=defaultdict(int)
        for (j,k),v in T.items():
            for a in range(abs(j-n),j+n+1,2): R[a,k]+=v
            for a in range(abs(k-n),k+n+1,2): R[j,a]+=e*v
        T={jk:v for jk,v in R.items() if v}
    return T

def pcoeff(word):
    D=sum(n for n,e in word)
    T={(0,0):1}
    for n,e in word:
        B=defaultdict(int)
        for (a,j),v in T.items():
            for k in range(min(n,(2*D-a)//2)+1): B[a+2*k,j]+=v
            if a+n<=2*D:
                for k in range(abs(j-n),j+n+1,2): B[a+n,k]+=e*v
        T={ij:v for ij,v in B.items() if v}
    return [T.get((2*j,0),0) for j in range(D+1)]

def cg(ns,j):
    if not ns: return int(j==0)
    if len(ns)==1: return int(j==ns[0])
    w=sum(ns)
    if j<0 or (w-j)%2 or j>w: return 0
    a=(w-j)//2;L=len(ns)
    out=0
    for mask in range(1<<L):
        z=a-sum(n+1 for i,n in enumerate(ns) if mask>>i&1)
        if z>=0: out+=(-1)**mask.bit_count()*comb(z+L-2,L-2)
    return out

def fused(ns):
    t={0:1}
    for n in ns:
        r=defaultdict(int)
        for j,a in t.items():
            for k in range(abs(j-n),j+n+1,2): r[k]+=a
        t=r
    return t

rng=Random(9501)
for _ in range(80):
    ns=[rng.randrange(1,10) for _ in range(rng.randrange(1,5))]
    row=fused(ns)
    assert all(cg(ns,j)==row.get(j,0) for j in range(sum(ns)+2))
print("80 inclusion-exclusion / fusion comparisons: PASS",flush=True)

legs=[]
for j in range(13):
    d=ddiag(2*j,j)
    assert all(d[k]==0 for k in range(1,len(d),2))
    legs.append(d[::2])
for i,a in enumerate(legs):
    for j,b in enumerate(legs):
        prod=multiply(a,b)
        val=sum(Q(v,k+1) for k,v in enumerate(prod))
        assert val==(Q(1,2*i+1) if i==j else 0)
print("169 exact shifted-Legendre inner products: PASS",flush=True)

def trace_polynomial(B):
    K={(0,0,0):1}
    for n,e in B:
        f=defaultdict(int)
        for j in range(n+1):
            f[j,j,0]+=1
            for k,v in enumerate(ddiag(n,j)): f[j,n-j,k]+=e*v
        R=defaultdict(int)
        for (i,j,k),a in K.items():
            for (x,y,z),b in f.items():
                if b: R[i+x,j+y,k+z]+=a*b
        K={t:a for t,a in R.items() if a}
    D=sum(n for n,e in B)
    ans=[0]*(D//2+1)
    for (i,j,k),v in K.items():
        if i==j:
            assert k%2==0 and k<=D
            ans[k//2]+=v
    return ans

for _ in range(24):
    B=[(rng.randrange(1,5),rng.choice((-1,1)))
       for j in range(rng.randrange(1,6))]
    D=sum(n for n,e in B);M=D//2
    f=trace_polynomial(B)
    P1=sum(pcoeff(B))
    assert sum(Q(v,j+1) for j,v in enumerate(f))==P1
    expansion=[Q(0)]*(M+1)
    for j in range(M+1):
        g=multiply(f,legs[j])
        coef=(2*j+1)*sum(Q(v,k+1) for k,v in enumerate(g))
        for k,v in enumerate(legs[j]): expansion[k]+=coef*v
    assert expansion==f
    for t in (Q(0),Q(1,9),Q(1,2),Q(9,25),Q(1)):
        v=sum(a*t**j for j,a in enumerate(f))
        assert 0<=v<=(M+1)**2*P1
print("24 exact trace-degree / norm profiles: PASS",flush=True)

def direct_with_cores(B,cores):
    t=tensor(B);D=sum(n for n,e in B)
    vals={}
    for mask in range(1<<len(cores)):
        ns=[n for i,(n,e) in enumerate(cores) if mask>>i&1]
        vals[mask]=[cg(ns,j) for j in range(D+1)]
    out=0;full=(1<<len(cores))-1
    for mask in range(1<<len(cores)):
        eps=1
        for i,(n,e) in enumerate(cores):
            if mask>>i&1: eps*=e
        out+=eps*sum(v*vals[full^mask][j]*vals[mask][k]
                     for (j,k),v in t.items())
    return out,t

def endpoint_test(B,d,s1,s2):
    D=sum(n for n,e in B)
    assert D>=d and max([n for n,e in B]+[0])<=D
    c=pcoeff(B)
    A=[1]+[0]*d
    for n,e in B:
        for j in range(d,n-1,-1): A[j]+=e*A[j-n]
    ad=A[d]
    assert c[d]*(d+1)>=ad*ad and all(x>=1 for x in c)
    F=sum(c[:d+1])+(s1+s2)*ad+s1*s2
    ep=s1*s2*(-1)**sum(e<0 for n,e in B)
    V,_=direct_with_cores(B,[(d,s1),(d,s2),(D,ep)])
    assert V==2*F and F>=0
    return F,ad

for _ in range(50):
    B=[(rng.randrange(1,5),rng.choice((-1,1)))
       for j in range(rng.randrange(2,8))]
    D=sum(n for n,e in B)
    d=rng.randrange(max(n for n,e in B),D+1)
    for s1,s2 in ((1,1),(1,-1),(-1,1),(-1,-1)):
        endpoint_test(B,d,s1,s2)
print("200 double-endpoint sign kernels: PASS",flush=True)

for r in (3,13,20):
    B=[(1,1)]*r+[(1,-1)]*r+[(2,1)]*r
    d=4*r-4;e=(-1)**r
    F,A=endpoint_test(B,d,e,e)
    assert A==(-1)**(r-1)*r and F>=2*r-2
    print("Endpoint family",r,"delta",d,"F",F,"A",A,flush=True)

def quartet_test(B,ns,p,signs):
    D=sum(n for n,e in B)
    a,b,c=sorted(ns)
    ns=(a,b,c)
    assert a>D and p>=c
    assert (D+a+b+c-p)%2==0
    d=(D+a+b+c-p)//2
    assert d>=c
    ep=(-1)**sum(e<0 for n,e in B)
    for e in signs: ep*=e
    cores=list(zip((a,b,c),signs))+[(p,ep)]
    V,T=direct_with_cores(B,cores)
    cc=pcoeff(B);P1=sum(cc)
    A=min(d,p)-D+1
    K=(D//2+1)**2
    geom=[cg((a,b,c),p-D+2*i) for i in range(D+1)]
    assert min(geom)==A
    F0=sum(x*y for x,y in zip(geom,cc))
    direct0=sum(v*cg((a,b,c,p),j) for (j,k),v in T.items() if k==0)
    assert F0==direct0
    pair_sum=0
    for i,j,l in ((0,1,2),(0,2,1),(1,2,0)):
        jj=sum(v*cg((ns[l],p),x)*cg((ns[i],ns[j]),y)
               for (x,y),v in T.items())
        assert abs(jj)<=K*P1
        pair_sum+=signs[i]*signs[j]*jj
    assert V==2*(F0+pair_sum)
    assert abs(V//2-F0)<=3*K*P1
    assert V//2>=(A-3*K)*P1
    if A>=3*K: assert V>=0
    return V//2,F0,P1,A,K

tests=0
for _ in range(60):
    B=[(rng.randrange(1,5),rng.choice((-1,1)))
       for j in range(rng.randrange(1,7))]
    D=sum(n for n,e in B)
    a=D+1+rng.randrange(30)
    b=a+rng.randrange(D//3+1)
    c=b+rng.randrange(D//3+1)
    ps=[p for p in range(c,D+a+b-c+1) if (D+a+b+c-p)%2==0]
    if not ps: continue
    p=rng.choice(ps)
    quartet_test(B,(a,b,c),p,tuple(rng.choice((-1,1)) for _ in range(3)))
    tests+=1
print("Quartet error bounds:",tests,"PASS",flush=True)

for r in (1,3,13,20):
    D=4*r;n=12*r*r+16*r+3
    B=[(1,1)]*r+[(1,-1)]*r+[(2,1)]*r
    F,F0,P1,A,K=quartet_test(B,(n,n,n),n+2,(1,1,(-1)**r))
    assert A-3*K>=1 and F>=P1
    C0=96**r
    assert n*P1<Q(3*(D+1)+24,8)*C0
    print("Quartet family",r,"D",D,"n",n,"F",F,
          "P1",P1,"margin",A-3*K,flush=True)

for s in range(1,13):
    B=[(1,1)]*(2*s)
    T=tensor(B)
    J=sum(v for (j,k),v in T.items() if j%2==k%2==0)
    P1=sum(pcoeff(B))
    assert J==comb(2*s,s)**2
    assert P1==comb(4*s+2,2*s+1)//(2*s+2)
for s in list(range(1,101))+[1000]:
    J=comb(2*s,s)**2
    P1=comb(4*s+2,2*s+1)//(2*s+2)
    assert 64*s*s*J*J >= (s+1)**2*(6*s+4)*P1*P1
    if s in (10,100,1000):
        print("Pair-error family s =",s,"floor(J/P1) =",J//P1,flush=True)
print("Constant-error obstruction: PASS",flush=True)
print("ALL CHECKS PASS")