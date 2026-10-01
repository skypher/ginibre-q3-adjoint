from fractions import Fraction as Q
from math import comb
from collections import defaultdict
from functools import lru_cache
from random import Random

def rotation_diagonal(n,j):
    a=[0]*(n+1)
    for l in range(min(j,n-j)+1):
        for k in range(l+1):
            a[n-2*l+2*k]+=(-1)**(l+k)*comb(j,l)*comb(n-j,l)*comb(l,k)
    return a

def value(a,c):
    return sum(v*c**i for i,v in enumerate(a))

identities=0
for n in range(1,19):
    for j in range(n+1):
        a=[0]*(n+1)
        for l in range(min(j,n-j)+1):
            a[n-2*l]=(-1)**l*comb(n-l,l)*comb(n-2*l,j-l)
        assert a==rotation_diagonal(n,j)
        assert a==rotation_diagonal(n,n-j)
        for c in (Q(0),Q(1,3),Q(3,5),Q(1)):
            x=value(a,c)
            assert abs(x)<=1
            if 0<c<1: assert abs(x)<1
        identities+=1
print("Rotation identities and rational matrix eigenvalues:",
      identities,"PASS",flush=True)

def pcoeff(word,cap=None):
    D=sum(n for n,e in word)
    if cap is None: cap=D
    T={(0,0):1}
    for n,e in word:
        B=defaultdict(int)
        for (a,j),v in T.items():
            for k in range(min(n,(2*cap-a)//2)+1):
                B[a+2*k,j]+=v
            if a+n<=2*cap:
                for k in range(abs(j-n),j+n+1,2):
                    B[a+n,k]+=e*v
        T={ij:v for ij,v in B.items() if v}
    return [T.get((2*j,0),0) for j in range(cap+1)]

def integrated_gram(word):
    # Keep c as a formal variable; integrate its powers exactly.
    K={(0,0,0):1}
    D=sum(n for n,e in word)
    for n,e in word:
        f=defaultdict(int)
        for j in range(n+1):
            f[j,j,0]+=1
            for k,v in enumerate(rotation_diagonal(n,j)):
                f[j,n-j,k]+=e*v
        H=defaultdict(int)
        for (i,j,k),a in K.items():
            for (x,y,z),b in f.items():
                if b: H[i+x,j+y,k+z]+=a*b
        K={t:a for t,a in H.items() if a}
    out=[Q(0) for j in range(D+1)]
    for (i,j,k),a in K.items():
        if i==j: out[i]+=Q(2*a,k+2)
    return out

rng=Random(9301)
for case in range(48):
    word=[(rng.randrange(1,6),rng.choice((-1,1)))
          for _ in range(rng.randrange(1,6))]
    pp=pcoeff(word)
    assert pp==integrated_gram(word)
    assert all(v>=1 for v in pp)
print("48 complete integrated-Gram / fusion profiles: PASS",flush=True)

@lru_cache(None)
def multiplicity(ns):
    if not ns: return 1
    if sum(ns)%2 or 2*max(ns)>sum(ns): return 0
    row={0:1}
    for n in ns:
        nxt=defaultdict(int)
        for j,a in row.items():
            for k in range(abs(j-n),j+n+1,2): nxt[k]+=a
        row=nxt
    return row.get(0,0)

def even(word):
    ans=0
    for mask in range(1<<len(word)):
        a=[];b=[];sg=1
        for i,(n,e) in enumerate(word):
            if mask>>i&1: a.append(n);sg*=e
            else: b.append(n)
        ans+=sg*multiplicity(tuple(sorted(a)))*multiplicity(tuple(sorted(b)))
    return ans

bridges=0
for case in range(60):
    bg=[(rng.randrange(1,8),rng.choice((-1,1)))
        for _ in range(rng.randrange(0,5))]
    D=sum(n for n,e in bg)
    h=1+case%3
    delta=rng.randrange(D//2+1) if h==1 else rng.randrange(13)
    B=2*delta+max([n for n,e in bg]+[0])+4
    rest=bg+[(B,rng.choice((-1,1))) for _ in range(h)]
    p=sum(n for n,e in rest)-2*delta
    ep=(-1)**sum(e<0 for n,e in rest)
    assert p>=max(n for n,e in rest)
    c=pcoeff(bg,delta)
    F=c[delta] if h==1 else sum(
        c[j]*comb(delta-j+h-2,h-2) for j in range(delta+1))
    lower=1 if h==1 else sum(
        comb(delta-j+h-2,h-2) for j in range(min(delta,D)+1))
    assert even(rest+[(p,ep)])==2*F
    assert F>=lower>=1
    bridges+=1
print("Saturation bridges and integer lower bounds:",
      bridges,"PASS",flush=True)

# Independent binomial/fusion evaluation of the r-family.
R=80
powers=[[1]]
for s in range(2*R):
    out=[0]*min(R+1,len(powers[-1])+3)
    for i,a in enumerate(powers[-1]):
        for k in range(min(4,len(out)-i)): out[i+k]+=a
    powers.append(out)

moments=[]
row={0:1}
for j in range(R//3+1):
    moments.append((row.get(0,0),row.get(4,0)))
    for _ in range(2):
        nxt=defaultdict(int)
        for spin,a in row.items():
            for k in range(abs(spin-3),spin+4,2): nxt[k]+=a
        row=nxt

results={}
for r in range(R+1):
    out=[0]*(r+1)
    for j in range(r//3+1):
        a,b=moments[j]
        sg=(-1)**j*comb(r,j)
        for i,v in enumerate(powers[2*(r-j)][:r-3*j+1]):
            z=i+3*j
            for k in range(min(5,r-z+1)): out[z+k]+=sg*a*v
            if z+2<=r: out[z+2]+=sg*b*v
    assert all(v>=1 for v in out)
    if r<=12:
        assert out==pcoeff([(3,1)]*r+[(3,-1)]*r+[(4,1)],r)
    results[r]=sum(out)
    assert results[r]>=r+1
for r,F in ((9,5329441),(11,203647209),(15,309937000389),
            (21,19541886292558037),(31,2131529581761524209753993)):
    assert results[r]==F
    print("r =",r,"F_r =",F,flush=True)
print("r = 0..80: all requested coefficients positive; prefix bound PASS",
      flush=True)

# The previous negative frequency channel is still negative.
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out

a=[0]*4+[1]
for f in ([1,0,0,0,-1],[1,-1],[1,0,0,1]):
    a=mul(a,f)
f=[0]*(6*13+1)
for j in range(14): f[6*j]=(-1)**j*comb(13,j)
a=mul(a,f)
energy=sum(a[j]**2-(a[j-1] if j else 0)*a[j+1] for j in range(16))
assert energy==-756 and results[15]>0
print("Old single-channel value -756; complete F_15 positive: PASS")
print("ALL CHECKS PASS")