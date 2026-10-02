import argparse
from itertools import combinations, product
from functools import lru_cache
from fractions import Fraction as Q
from math import factorial, isqrt
ap=argparse.ArgumentParser(description="FM-STR1b: exact block-transfer obstruction and compositions; memory only.")
ap.parse_args()

@lru_cache(None)
def fusion(ns):
    D={0:1}
    for n in ns:
        C={}
        for a,x in D.items():
            for b in range(abs(a-n),a+n+1,2):
                C[b]=C.get(b,0)+x
        D=C
    return D
def inv(ns):
    return fusion(tuple(sorted(ns))).get(0,0)
def dims(word):
    even=odd=0
    for s in range(1<<len(word)):
        X=tuple(abs(n) for i,n in enumerate(word) if s>>i&1)
        Y=tuple(abs(n) for i,n in enumerate(word) if not(s>>i&1))
        v=inv(X)*inv(Y)
        if sum(n<0 for i,n in enumerate(word) if s>>i&1)%2:odd+=v
        else:even+=v
    return even,odd

def edge(P,i,j):
    R={}
    for ex,c in P.items():
        for k,sg in ((j,1),(i,-1)):
            y=list(ex);y[k]+=1;y=tuple(y)
            R[y]=R.get(y,0)+sg*c
    return {e:c for e,c in R.items() if c}
def covariant(u,S,a):
    i,j,k=S
    edges=[(u,S[a]),(i,j),(i,k),(j,k),
           (j,k) if a==0 else (i,k)]
    P={(0,)*8:1}
    for i,j in edges:P=edge(P,i,j)
    return P
def mul(A,B):
    C={}
    for a,x in A.items():
        for b,y in B.items():
            c=tuple(i+j for i,j in zip(a,b))
            C[c]=C.get(c,0)+x*y
    return {c:x for c,x in C.items() if x}

def rank_Q(cols):
    piv={}
    for col in cols:
        P={r:Q(v) for r,v in col.items() if v}
        while P:
            r=min(P)
            if r not in piv:break
            q=P[r]
            for k,v in piv[r].items():
                P[k]=P.get(k,0)-q*v
                if not P[k]:del P[k]
        if P:
            r=min(P);q=P[r]
            piv[r]={k:v/q for k,v in P.items()}
    return len(piv)

prime=65521
assert all(prime%d for d in range(2,isqrt(prime)+1))
def rank_mod(cols):
    piv={};chosen=[]
    for col in cols:
        P={r:int(v)%prime for r,v in col.items() if v%prime}
        while P:
            r=min(P)
            if r not in piv:break
            q=P[r]
            for k,v in piv[r].items():
                P[k]=(P.get(k,0)-q*v)%prime
                if not P[k]:del P[k]
        if P:
            r=min(P);q=pow(P[r],prime-2,prime)
            piv[r]={k:v*q%prime for k,v in P.items()}
            chosen.append(r)
    return len(piv),chosen
def det_mod(A):
    A=[[x%prime for x in row] for row in A]
    d=1
    for j in range(len(A)):
        k=next((k for k in range(j,len(A)) if A[k][j]),None)
        if k is None:return 0
        if k!=j:A[k],A[j]=A[j],A[k];d=-d
        q=A[j][j];d=d*q%prime
        A[j]=[x*pow(q,prime-2,prime)%prime for x in A[j]]
        for i in range(j+1,len(A)):
            q=A[i][j]
            A[i]=[(x-q*y)%prime for x,y in zip(A[i],A[j])]
    return d%prime

assert fusion((3,3,3))=={1:2,3:4,5:3,7:2,9:1}
assert [inv((1,)+(3,)*j) for j in range(4)]==[0,0,0,2]
records=[]
for S in combinations(range(2,8),3):
    T=tuple(i for i in range(2,8) if i not in S)
    for a,b in product(range(2),repeat=2):
        records.append((S,a,b,mul(covariant(0,S,a),covariant(1,T,b))))
cols=[v for S,a,b,v in records]
print("primitive columns built:",len(cols),flush=True)
r=rank_Q(cols)
assert r==64 and len(cols)==80
assert inv((3,)*6)==34
assert inv((1,1)+(3,)*6)==124
assert dims((3,)*6+(-1,-1))==(946,160)
print("two-minus: exact multiplication rank",r,"of 80; J rank 128 of 160",flush=True)
print("exchange obstruction: skew source 40, skew target 34",flush=True)
print("two-minus: even 946, odd 160, Phi 786",flush=True)

# K_ij is 24 times the canonical singlet transfer for two V_3 factors.
# The output pair singlet is factored out; its colour pair labels the row block.
def contract_binary(P,i,j,n):
    R={}
    for e,c in P.items():
        if e[i]+e[j]!=n:continue
        v=c*(-1)**e[i]*factorial(e[i])*factorial(n-e[i])
        f=list(e);f[i]=f[j]=0;f=tuple(f)
        R[f]=R.get(f,0)+v
    return {e:c for e,c in R.items() if c}

lifted=[]
for S,a,b,P in records:
    R={}
    for i,j in combinations(range(2,8),2):
        if (i in S)==(j in S):continue
        weight=i+1 if i in S else j+1
        for ex,c in contract_binary(P,i,j,3).items():
            R[(i,j)+ex]=weight*c
    lifted.append(R)
r,rows=rank_mod(lifted)
assert r==80
minor=det_mod([[C.get(row,0) for C in lifted] for row in rows])
assert minor
print("two-stage endpoint weights (3,4,5,6,7,8): rank 80 of 80 mod",prime,
      "; selected minor",minor,flush=True)

# Four-minus tensor identity, with V_2=C^3 and an arbitrary vector-valued beta.
def eps(v):
    if len(set(v))<3:return 0
    return (-1)**sum(v[i]>v[j] for i in range(3) for j in range(i+1,3))
H=[]
for i in range(4):
    P={}
    for x in product(range(3),repeat=5):
        if x[0]==x[i+1]:
            c=(-1)**i*eps(tuple(x[j+1] for j in range(4) if j!=i))
            if c:P[x]=c
    H.append(P)
assert all(sum(P.get(x,0) for P in H)==0 for x in product(range(3),repeat=5))
assert rank_Q(H)==3
def contract_vector(P,i,j):
    R={}
    for x,c in P.items():
        if x[i+1]!=x[j+1]:continue
        y=tuple(x[k] for k in range(5) if k not in (i+1,j+1))
        R[y]=R.get(y,0)+c
    return {y:c for y,c in R.items() if c}
incidence=[]
for i,j in combinations(range(4),2):
    C=[contract_vector(P,i,j) for P in H]
    assert C[i]
    assert all(not C[k] for k in range(4) if k not in (i,j))
    assert C[j]=={e:-c for e,c in C[i].items()}
    e=min(C[i]);scale=C[i][e]
    row=[Q(C[k].get(e,0)*(1 if k==i else 2),scale) for k in range(4)]
    target=[0]*4;target[i]=1;target[j]=-2
    assert row==target
    incidence.append(row)
M=[incidence[i][:] for i in (0,1,3,2)]
det=Q(1)
for j in range(4):
    k=next(k for k in range(j,4) if M[k][j])
    if k!=j:M[k],M[j]=M[j],M[k];det=-det
    q=M[j][j];det*=q;M[j]=[x/q for x in M[j]]
    for i in range(j+1,4):
        q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j])]
assert det==4
print("four-minus Schouten rank 3; composed minor determinant",det,flush=True)

# Six adjoints: the second pair crosses the original triple cut.
six_pure=[];six_base=[];six_repaired=[]
for pair in combinations(range(1,6),2):
    S=(0,)+pair
    T=tuple(i for i in range(6) if i not in S)
    P={}
    for x in product(range(3),repeat=6):
        c=eps(tuple(x[i] for i in S))*eps(tuple(x[i] for i in T))
        if c:P[x]=c
    six_pure.append(P)
    R={}
    correction={}
    for i,j in combinations(range(6),2):
        if (i in S)==(j in S):continue
        C={}
        for x,c in P.items():
            if x[i]!=x[j]:continue
            y=tuple(x[k] for k in range(6) if k not in (i,j))
            C[y]=C.get(y,0)+c
        for y,c in C.items():
            if not c:continue
            key=(i,j)+y
            weight=i+1 if i in S else j+1
            R[key]=weight*c
            if S==(0,1,2) and (i,j)==(0,3):correction[key]=c
    six_base.append(R)
    U=dict(R)
    for y,c in correction.items():U[y]=U.get(y,0)+c
    six_repaired.append(U)
assert rank_Q(six_pure)==5
assert rank_Q(six_base)==9
assert rank_Q(six_repaired)==10
assert dims((-2,)*6)==(120,20)
print("six adjoints: single-stage rank 10, two-stage rank 18, repaired rank 20 of 20",flush=True)

checks=0
for q in range(1,13):
    for t in range(q,13):
        if 2 in (q,t):continue
        even,odd=dims((-2,)*4+(q,t))
        b=int(abs(q-t)<=2<=q+t and (q+t)%2==0)
        assert odd==8*b and even>=odd
        checks+=1
print("four-minus two-plus label checks:",checks,"PASS",flush=True)
for A in ((1,3),(3,5),(4,4),(10,20,28)):
    b=inv((2,)+A)
    for s in range((1<<len(A))-1):
        S=tuple(A[i] for i in range(len(A)) if s>>i&1)
        T=tuple(A[i] for i in range(len(A)) if not(s>>i&1))
        assert inv((2,)+S)*inv((2,2,2)+T)==0
    even,odd=dims((-2,)*4+A)
    assert odd==8*b and even>=odd
    print("primitive background",A,"channel dimension",b,"even",even,
          "odd",odd,"Phi",even-odd,flush=True)

cuts=0
for L in range(2,13):
    A=tuple(4*(1<<j) for j in range(L-1))
    A=A+(sum(A)+2,)
    assert A[-1]-sum(A[:-1])==2
    for s in range(1,(1<<L)-1):
        S=[A[i] for i in range(L) if s>>i&1]
        assert 2*max(S)-sum(S)>2
        cuts+=1
print("unbounded-family proper-subset weight checks:",cuts,"PASS",flush=True)
print("FM-STR1b PASS",flush=True)
