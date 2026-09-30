"""FM-SEC113 (luna_max_venus): H_AC_q screen L = 8 (all labels <= 3), L = 9, 10 structured lists; dictionary separator at (1^8,2^2), q^4."""
from itertools import combinations, combinations_with_replacement
from functools import lru_cache
from fractions import Fraction as Q
from scipy.optimize import linprog
from sympy import Matrix, Rational

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:
        p.pop()
    return tuple(p)

def padd(a,b):
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                 for i in range(max(len(a),len(b)))])

def psub(a,b):
    return trim([(a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0)
                 for i in range(max(len(a),len(b)))])

def pmul(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            z[i+j]+=x*y
    return trim(z)

def pshift(a,k):
    return (0,)*k+tuple(a)

@lru_cache(None)
def qbinom(n,k):
    if not 0<=k<=n:
        return (0,)
    if k==0 or k==n:
        return (1,)
    return padd(qbinom(n-1,k),pshift(qbinom(n-1,k-1),n-k))

@lru_cache(None)
def qfactorial(n):
    z=(1,)
    for j in range(1,n+1):
        z=pmul(z,(1,)*j)
    return z

@lru_cache(None)
def linearization(a,b,k):
    return pmul(pmul(qbinom(a,k),qbinom(b,k)),qfactorial(k))

@lru_cache(None)
def moment(labels):
    states={0:(1,)}
    for n in labels:
        nxt={}
        for j,p in states.items():
            for k in range(min(j,n)+1):
                degree=j+n-2*k
                term=pmul(p,linearization(j,n,k))
                nxt[degree]=padd(nxt.get(degree,(0,)),term)
        states=nxt
    return states.get(0,(0,))

def qtable(labels):
    out=[]
    for S in range(1<<len(labels)):
        left=tuple(sorted(labels[i] for i in range(len(labels)) if S>>i&1))
        right=tuple(sorted(labels[i] for i in range(len(labels)) if not(S>>i&1)))
        out.append(pmul(moment(left),moment(right)))
    return out

def profile(S,labels):
    labs=sorted(set(labels))
    w=tuple(sum(labels[i]==v and bool(S>>i&1) for i in range(len(labels)))
            for v in labs)
    total=tuple(labels.count(v) for v in labs)
    comp=tuple(total[i]-w[i] for i in range(len(labs)))
    return min(w,comp)

def coefficient(p,k):
    return p[k] if k<len(p) else 0

def wht(polys):
    v=list(polys)
    h=1
    while h<len(v):
        for i in range(0,len(v),2*h):
            for j in range(i,i+h):
                v[j],v[j+h]=padd(v[j],v[j+h]),psub(v[j],v[j+h])
        h*=2
    return v

# Exact Walsh screen.
for L in (8,9,10):
    words=[w for w in combinations_with_replacement((1,2,3),L)
           if sum(w)%2==0]
    entries=negative=zeros=0
    maxdeg=0
    minpos=None
    for word in words:
        F=wht(qtable(word))
        maxdeg=max(maxdeg,max(map(len,F))-1)
        for T in range(1<<L):
            if T.bit_count()%2:
                continue
            for k,x in enumerate(F[T]):
                entries+=1
                if x<0:
                    negative+=1
                elif x==0:
                    zeros+=1
                elif minpos is None or x<minpos[0]:
                    minpos=(x,word,T,k)
    print("FOURIER",L,len(words),entries,maxdeg,negative,zeros,minpos)

def subspaces(d):
    # RREF basis enumeration; members are elements of each subspace.
    for k in range(d+1):
        for piv in combinations(range(d),k):
            free=[(i,j) for i,p in enumerate(piv)
                  for j in range(p+1,d) if j not in piv]
            for mask in range(1<<len(free)):
                rows=[1<<p for p in piv]
                for h,(i,j) in enumerate(free):
                    if mask>>h&1:
                        rows[i]|=1<<j
                members=[0]
                for row in rows:
                    members += [x^row for x in members]
                yield members

def subspace_columns(labels):
    d=len(labels)-1
    G=1<<d
    profiles=sorted({profile(q,labels) for q in range(G)})
    ix={p:i for i,p in enumerate(profiles)}
    qclass=[ix[profile(q,labels)] for q in range(G)]
    sizes=[qclass.count(i) for i in range(len(profiles))]
    cols={}
    nsub=0
    for members in subspaces(d):
        nsub+=1
        hist=[0]*len(profiles)
        for q in members:
            hist[qclass[q]]+=1
        h=len(members)
        col=tuple(Q(h*hist[i],sizes[i]) for i in range(len(profiles)))
        if any(col):
            cols.setdefault(col,1)
    return list(cols),nsub

def binary_profile_columns(labels,max_support):
    d=len(labels)-1
    G=1<<d
    profiles=sorted({profile(q,labels) for q in range(G)})
    ix={p:i for i,p in enumerate(profiles)}
    P=len(profiles)
    reps=[next(q for q in range(G) if ix[profile(q,labels)]==i)
          for i in range(P)]
    K=[[[0]*P for _ in range(P)] for __ in range(P)]
    for o,q in enumerate(reps):
        for x in range(G):
            K[o][ix[profile(x,labels)]][ix[profile(x^q,labels)]]+=1
    cols={}
    for n in range(1,min(P,max_support)+1):
        for supp in combinations(range(P),n):
            chosen=set(supp)
            col=tuple(sum(K[o][a][b] for a in chosen for b in chosen)
                      for o in range(P))
            if any(col):
                cols.setdefault(col,1)
    return list(cols)

def exact_lp(cols,y):
    res=linprog([1.0]*len(cols),
        A_eq=[[float(cols[j][i]) for j in range(len(cols))]
              for i in range(len(y))],
        b_eq=list(map(float,y)),bounds=(0,None),method="highs")
    if not res.success:
        raise AssertionError(("LP failed",res.message))
    support=[j for j,x in enumerate(res.x) if x>1e-8]
    def rat(x):
        if isinstance(x,Q):
            return Rational(x.numerator,x.denominator)
        return Rational(int(x))
    A=Matrix([[rat(cols[j][i]) for j in support] for i in range(len(y))])
    sol,params=A.gauss_jordan_solve(Matrix(y))
    if params.rows:
        raise AssertionError("LP support has free parameters")
    coeff=[Q(int(v.p),int(v.q)) for v in sol]
    if min(coeff)<0:
        raise AssertionError("negative exact coefficient")
    got=[sum(coeff[h]*cols[j][i] for h,j in enumerate(support))
         for i in range(len(y))]
    if got!=list(map(Q,y)):
        raise AssertionError("exact re-expansion failed")
    return [(j,coeff[h]) for h,j in enumerate(support)],res.fun

def certify(labels,kind,max_support=4):
    tab=qtable(labels)
    if kind=="subspace":
        cols,nsub=subspace_columns(labels)
    else:
        cols=binary_profile_columns(labels,max_support)
        nsub=None
    d=len(labels)-1
    profiles=sorted({profile(q,labels) for q in range(1<<d)})
    ix={p:i for i,p in enumerate(profiles)}
    reps=[next(q for q in range(1<<d) if ix[profile(q,labels)]==i)
          for i in range(len(profiles))]
    degree=max(map(len,tab))-1
    tested=term_sum=max_support_seen=0
    min_obj=None
    for k in range(degree+1):
        y=[coefficient(tab[q],k) for q in reps]
        if not any(y):
            continue
        terms,obj=exact_lp(cols,y)
        got=[sum(c*cols[j][i] for j,c in terms)
             for i in range(len(profiles))]
        assert got==list(map(Q,y))
        for q in range(1<<d):
            assert got[ix[profile(q,labels)]]==coefficient(tab[q],k)
        tested+=1
        term_sum+=len(terms)
        max_support_seen=max(max_support_seen,len(terms))
        if min_obj is None or obj<min_obj[0]:
            min_obj=(obj,k)
    print("CERT",labels,kind,"subspaces",nsub,
          "degree",degree,"vectors",tested,
          "columns",len(cols),"max terms",max_support_seen,
          "sum terms",term_sum,"min LP objective",min_obj)

# All length-eight even-total lists.
L8=[w for w in combinations_with_replacement((1,2,3),8) if sum(w)%2==0]
for w in L8:
    certify(w,"subspace")
print("L8 SUMMARY",len(L8),
      "coefficient vectors",sum(max(map(len,qtable(w))) for w in L8),
      "certificates use a fixed dictionary per list")

# Longer samples and the named structured lists.
certify((1,)*8+(2,),"subspace")
certify((1,)*6+(2,)*3,"subspace")
certify((1,)*10,"binary",6)
certify((1,)*9+(3,),"binary",4)
certify((1,)*10+(4,),"binary",4)

from random import Random

def fuse_row(row,n):
    out={}
    for j,v in row.items():
        for k in range(abs(j-n),j+n+1,2):
            out[k]=out.get(k,0)+v
    return out

def subset_rows(labels):
    rows=[{0:1}]
    for n in labels:
        old=rows[:]
        for row in old:
            rows.append(fuse_row(row,n))
    return rows

L=(1,)*8+(2,2)
d=len(L)-1
G=1<<d
profiles=sorted({profile(q,L) for q in range(G)})
ix={p:i for i,p in enumerate(profiles)}
P=len(profiles)
reps=[next(q for q in range(G) if ix[profile(q,L)]==i)
      for i in range(P)]
K=[[[0]*P for _ in range(P)] for __ in range(P)]
for o,q in enumerate(reps):
    for x in range(G):
        K[o][ix[profile(x,L)]][ix[profile(x^q,L)]]+=1

def autocorr_profile(p):
    return tuple(sum(K[o][a][b]*p[a]*p[b]
                    for a in range(P) for b in range(P))
                 for o in range(P))

base=[tuple(int(i==j) for j in range(P)) for i in range(P)]
rows=subset_rows(L)
full=(1<<len(L))-1
for a in range(sum(L)+1):
    for b in range(a,sum(L)+1):
        vals=[]
        for q in range(G):
            x=rows[q].get(a,0)*rows[full^q].get(b,0)
            y=rows[q].get(b,0)*rows[full^q].get(a,0)
            vals.append(x if a==b else x+y)
        p=tuple(vals[reps[o]] for o in range(P))
        if any(p):
            base.append(p)

ratios=(Q(1,4),Q(1,2),Q(1),Q(2),Q(4))
candidates=list(base)
for i,j in combinations(range(len(base)),2):
    for t in ratios:
        candidates.append(tuple(Q(base[i][s])+t*Q(base[j][s])
                                for s in range(P)))
rng=Random(113102)
for _ in range(15000):
    prob=rng.choice((.15,.25,.4,.6,.8,1.0))
    p=tuple(rng.randrange(1,9) if rng.random()<prob else 0
            for _ in range(P))
    if any(p):
        candidates.append(p)

atoms={}
for p in candidates:
    atoms.setdefault(autocorr_profile(p),p)
cols=list(atoms)
assert len(cols)==17520
tab=qtable(L)
for k in range(4):
    y=[coefficient(tab[q],k) for q in reps]
    terms,obj=exact_lp(cols,y)
    assert len(terms)==7
    print("q",k,"exact support",len(terms))

target=[coefficient(tab[q],4) for q in reps]
assert target==[1330,0,30,0,0,0,145,29,19,0,0,0,38,19]
eta=[0,8263955,712320,197584,338480,466400,0,
     -3205440,-712320,89040,-801360,-890400,1187200,2088730]
assert sum(eta[i]*target[i] for i in range(P))==-322770
assert all(sum(eta[i]*col[i] for i in range(P))>=0 for col in cols)
f=[coefficient(tab[q],4) for q in range(G)]
ft=[sum((-1)**((s&t).bit_count()%2)*v for t,v in enumerate(f))
    for s in range(G)]
assert min(ft)==814 and ft.index(min(ft))==15
print("q4 finite-dictionary separator verified; Fourier minimum",min(ft))
