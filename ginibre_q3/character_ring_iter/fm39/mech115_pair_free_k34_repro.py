import argparse
from math import comb, factorial, lcm
from functools import lru_cache
from itertools import combinations_with_replacement, product
from random import Random

ap=argparse.ArgumentParser(description="FM-MECH115 exact certificates; memory only")
ap.add_argument("--cores", type=int, choices=(3,4), default=3)
ap.add_argument("--distance", type=int, choices=(8,9,10), default=8)
ap.add_argument("--sector", choices=("plus","minus","both"), default="both")
args=ap.parse_args()
k,d=args.cores,args.distance
SCALE=factorial(2*d)

def mul(a,b,cap=None):
    out={}
    for (i,j),v in a.items():
        for (h,l),w in b.items():
            if cap is None or i+h+2*(j+l)<=cap:
                m=(i+h,j+l);out[m]=out.get(m,0)+v*w
    return {m:v for m,v in out.items() if v}

def add(a,b,c=1,dx=0,dt=0):
    for (i,j),v in b.items():
        m=(i+dx,j+dt);a[m]=a.get(m,0)+c*v
    return {m:v for m,v in a.items() if v}

def block(n,e):
    a={(0,j):1 for j in range(min(n,d)+1)}
    for j in range(n//2+1):
        m=(n-2*j,j)
        if n<=2*d:a[m]=a.get(m,0)+e*(-1)**j*comb(n-j,j)
    return {m:v for m,v in a.items() if v}

def newton(ns,es,tau=1):
    P={(0,0):1,(0,1):-1}
    for n,e in zip(ns,es):P=mul(P,block(n,e),2*d)
    u={(0,1):1,(1,0):1}
    v={(0,1):1-tau,(0,2):1,(2,0):tau}
    C=[[0]*(d+1) for _ in range(2*d+1)]
    for j in range(d+1):
        R=P
        for i in range(2*(d-j)+1):
            C[i][j]=sum(w*(comb(h,h//2)//(h//2+1))
                for (h,l),w in R.items() if h%2==0 and h//2+l==d)
            R=mul(R,u,2*d)
        P=mul(P,v,2*d)
    return C

def shift_a(C):
    for i in range(len(C)-1):
        for j in range(d+1):C[i][j]+=C[i+1][j]

def shift_b(row):
    for j in range(d):row[j]+=row[j+1]

def check_plus(C,minimum_weight):
    H=2*d+4;checks=0
    for a in range(H):
        row=C[0][:]
        for b in range(H):
            if a+2*b>=minimum_weight:
                assert row[0]>=0
                checks+=1
            shift_b(row)
        assert min(row)>=0
        shift_a(C)
    assert min(map(min,C))>=0
    return checks

# Monic Krawtchouk polynomials in (X,T).
K=[{(0,0):1},{(1,0):1}]
for j in range(1,2*d):
    p=add({},K[j],dx=1)
    p=add(p,K[j-1],-j,dt=1)
    K.append(add(p,K[j-1],j*(j-1)))
BB=[mul(K[j//2],K[(j+1)//2]) for j in range(2*d+1)]

@lru_cache(None)
def falling(j,r):
    a=[1]
    for i in range(r):
        b=[0]*(len(a)+1)
        for h,v in enumerate(a):b[h]-=(j+i)*v;b[h+1]+=v
        a=b
    return a

@lru_cache(None)
def moment(j,word):
    a={0:1}
    for n in (1,)*j+word:
        b={}
        for h,v in a.items():
            for ell in range(abs(h-n),h+n+1,2):b[ell]=b.get(ell,0)+v
        a=b
    return a.get(0,0)

def gb(n,r):
    if n<0:return (-1)**r*comb(r-n-1,r)
    return comb(n,r) if n>=r else 0

@lru_cache(None)
def balanced(D,word):
    if D<0 or sum(word)>2*D:return {}
    r=len(word);w=sum(word);p={}
    for j in range(w%2,2*D-w+1,2):
        M=moment(j,word)
        if not M:continue
        h=D-(w+j)//2
        for v in range(h+1):
            den=factorial(j)*factorial(h-v)
            assert SCALE%den==0
            c=M*gb(k-r+v-2,v)*(SCALE//den)
            if c:
                for z,b in enumerate(falling(j,h-v)):
                    if b:p=add(p,K[j],c*b,dt=z)
    out={}
    for j in range(2*D,-1,-1):
        row=[p.get((j,t),0) for t in range(d+1)]
        if any(row):
            out[j]=row
            for t,v in enumerate(row):
                if v:p=add(p,BB[j],-v,dt=t)
    assert not p
    return out

def coefficients(ns,es):
    out=[[0]*(d+1) for _ in range(2*d+1)]
    for modes in product(range(3),repeat=k):
        word=tuple(n for n,m in zip(ns,modes) if m==1)
        spent=sum(n+1 for n,m in zip(ns,modes) if m==2)
        if sum(word)>2*(d-spent):continue
        sg=1
        for e,m in zip(es,modes):
            if m==1:sg*=e
            if m==2:sg=-sg
        for j,row in balanced(d-spent,word).items():
            for i,v in enumerate(row):out[j][i]+=sg*v
    return out

def shifted(row,H):
    return [sum(row[j]*comb(j,i)*H**(j-i)
                for j in range(i,d+1)) for i in range(d+1)]

def tail_ok(B,H):
    rows=[shifted(row,H) for row in B]
    scale=2*lcm(*range(1,d+1))**2
    for j in range(d+1):
        for i in range(d+1):
            v=scale*rows[2*j][i]
            if j<d:v-=scale//2*(j+1)**2*abs(rows[2*j+1][i])
            if j:v-=scale//(2*j*j)*abs(rows[2*j-1][i])
            if v<0:return False
    return True

def at_t(B,t):
    out=[]
    for row in B:
        v=0
        for c in reversed(row):v=v*t+c
        out.append(v)
    return out

def at_x(row,t,x):
    a=[1,x]
    for j in range(1,d):a.append(x*a[-1]-j*(t-j+1)*a[-2])
    v=sum(c*a[j//2]*a[(j+1)//2] for j,c in enumerate(row))
    assert v%SCALE==0
    return v//SCALE

def two_spin(word,p):
    a={(0,0):1}
    for n,e in word:
        b={}
        for (j,h),v in a.items():
            for r in range(abs(j-n),j+n+1,2):b[r,h]=b.get((r,h),0)+v
            for r in range(abs(h-n),h+n+1,2):b[j,r]=b.get((j,r),0)+e*v
        a={m:v for m,v in b.items() if v}
    return a.get((p,0),0)

# Independent consumer controls, including sign normalization.
rng=Random(11500+10*k+d)
for _ in range(12):
    ns=tuple(sorted(rng.randrange(3,d+1) for _ in range(k)))
    em={n:rng.choice((-1,1)) for n in ns}
    es=tuple(em[n] for n in ns);sigma=rng.choice((-1,1))
    ee=tuple(e*sigma**n for n,e in zip(ns,es))
    b=rng.randrange(5);a=max(rng.randrange(5),2*d-sum(ns[:-1])-2*b)
    p=a+2*b+sum(ns)-2*d
    for tau in (-1,1):
        word=[(1,sigma)]*a+[(2,tau)]*b+list(zip(ns,es))
        if tau<0:value=at_x(at_t(coefficients(ns,ee),a+2*b),a+2*b,a)
        else:
            C=newton(ns,ee)
            value=sum(C[i][j]*comb(a,i)*comb(b,j)
                for i in range(min(a,2*d)+1) for j in range(min(b,d)+1))
        assert value==two_spin(word,p)
        ep=1
        for n,e in word:ep*=e
        assert two_spin(word+[(p,ep)],0)==2*value
print("24 independent consumer controls PASS",flush=True)

profiles=plus_checks=minus_checks=0;tails={}
for ns in combinations_with_replacement(range(3,d+1),k):
    labels=sorted(set(ns));minimum_weight=2*d-sum(ns[:-1])
    for signs in product((-1,1),repeat=len(labels)):
        em=dict(zip(labels,signs));es=tuple(em[n] for n in ns)
        profiles+=1
        if args.sector!="minus":
            plus_checks+=check_plus(newton(ns,es),minimum_weight)
        if args.sector!="plus":
            B=coefficients(ns,es)
            H=next(h for h in (0,8,16,24,32,48,64,96,128,192)
                   if tail_ok(B,h))
            tails[H]=tails.get(H,0)+1
            for t in range(max(0,minimum_weight),H):
                row=at_t(B,t)
                for x in range(t%2,t+1,2):
                    assert at_x(row,t,x)>=0
                    minus_checks+=1
print("PASS", "k",k,"d",d,"profiles",profiles,
      "plus finite",plus_checks,"minus finite",minus_checks,
      "tail thresholds",tails,flush=True)