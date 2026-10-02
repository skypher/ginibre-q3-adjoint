import sys
if any(x in ('-h','--help') for x in sys.argv[1:]):
    print('Exact independent FM-STR9 verifier; Python 3 standard library only.'); raise SystemExit
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
from collections import Counter
from random import Random
from math import comb, prod
from datetime import datetime, timezone

def stamp():
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), end=" ", flush=True)

@lru_cache(None)
def fuse(ns):
    d={0:1}
    for n in ns:
        z={}
        for a,c in d.items():
            for b in range(abs(a-n),a+n+1,2):
                z[b]=z.get(b,0)+c
        d=z
    return d

def clean(d): return {k:v for k,v in d.items() if v}

@lru_cache(None)
def signed_table(word):
    d={(0,0):1}
    for sn in word:
        n=abs(sn); eps=1 if sn>0 else -1
        z={}
        for (a,b),c in d.items():
            for t in range(abs(a-n),a+n+1,2):
                z[t,b]=z.get((t,b),0)+c
            for t in range(abs(b-n),b+n+1,2):
                z[a,t]=z.get((a,t),0)+eps*c
        d=clean(z)
    return d

def G_table(ns):
    out={}
    L=len(ns)
    for mask in range(1<<L):
        x=tuple(ns[i] for i in range(L) if mask>>i&1)
        y=tuple(ns[i] for i in range(L) if not mask>>i&1)
        for a,ca in fuse(x).items():
            for b,cb in fuse(y).items():
                out[a,b]=out.get((a,b),0)+ca*cb
    return clean(out)

def LR(table,n,left=True):
    out={}
    for (a,b),v in table.items():
        if left:
            for t in range(abs(n-a),n+a+1,2):
                out[t,b]=out.get((t,b),0)+v
        else:
            for t in range(abs(n-b),n+b+1,2):
                out[a,t]=out.get((a,t),0)+v
    return clean(out)

def add_scaled(dst,src,k):
    for ab,v in src.items():
        dst[ab]=dst.get(ab,0)+k*v
    return dst

def verify_prop1(word):
    A=tuple(word[::2]); B=tuple(word[1::2])
    na=[-v for v in A if v<0]; nb=[-v for v in B if v<0]
    fa=signed_table(A); fb=signed_table(B)
    if len(na)==1 and len(nb)==1:
        ca=tuple(v for v in A if v>0); cb=tuple(v for v in B if v>0)
        ga,gb=G_table(ca),G_table(cb)
        xa=add_scaled(LR(ga,na[0],True),LR(ga,na[0],False),-1)
        xb=add_scaled(LR(gb,nb[0],True),LR(gb,nb[0],False),-1)
        assert clean(xa)==fa and clean(xb)==fb
        for t,v in fa.items(): assert fa.get((t[1],t[0]),0)==-v
        for t,v in fb.items(): assert fb.get((t[1],t[0]),0)==-v
        odd={(a,b) for a,b in set(fa)|set(fb)
             if fa.get((a,b),0)*fb.get((a,b),0)<0}
        assert all((b,a) in odd for a,b in odd)
        return 1
    assert (len(na),len(nb)) in ((2,0),(0,2))
    side=A if len(na)==2 else B
    other=B if len(na)==2 else A
    neg=na if len(na)==2 else nb
    positives=tuple(v for v in side if v>0)
    g=G_table(positives)
    n,m=neg
    x={}
    for c in range(abs(n-m),n+m+1,2):
        add_scaled(x,LR(g,c,True),1)
        add_scaled(x,LR(g,c,False),1)
    add_scaled(x,LR(LR(g,m,False),n,True),-1)
    add_scaled(x,LR(LR(g,n,False),m,True),-1)
    actual=signed_table(side)
    assert clean(x)==actual
    assert all(actual.get((b,a),0)==v for (a,b),v in actual.items())
    gd=G_table(tuple(other))
    assert gd==signed_table(tuple(other))
    odd={(a,b) for a,b in set(actual)|set(gd)
         if actual.get((a,b),0)<0 and gd.get((a,b),0)>0}
    assert all(actual.get((b,a),0)==actual.get((a,b),0) for a,b in odd)
    return 2

rng=Random(900112)
prop_counts=[0,0]
for trial in range(800):
    L=rng.randint(4,12)
    mags=sorted(rng.randint(1,12) for _ in range(L))
    cnt=Counter(mags)
    choices=[(n,) for n,c in cnt.items() if c==2]
    choices+=list(combinations([n for n,c in cnt.items() if c==1],2))
    if not choices: continue
    neg=choices[rng.randrange(len(choices))]
    word=tuple(-n if n in neg else n for n in mags)
    prop_counts[verify_prop1(word)-1]+=1
assert min(prop_counts)>0
stamp(); print("Prop1 random exact tables",prop_counts,"PASS",flush=True)

@lru_cache(None)
def no_inv(ns):
    for mask in range(1,1<<len(ns)):
        part=tuple(ns[i] for i in range(len(ns)) if mask>>i&1)
        if fuse(part).get(0,0): return False
    return True

def tc_support(plus,n,m):
    support=[]; full=(1<<len(plus))-1
    for mask in range(full+1):
        x=tuple(plus[i] for i in range(len(plus)) if mask>>i&1)
        y=tuple(plus[i] for i in range(len(plus)) if not mask>>i&1)
        if fuse(x).get(n,0)*fuse(y).get(m,0):
            support.append(mask)
            if len(support)>2: return []
    if len(support)==2 and support[0]^support[1]!=full: return []
    return support

def det_small(A):
    n=len(A)
    if n==0:return F(1)
    total=F(0)
    for p in permutations(range(n)):
        inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term=F(-1 if inv&1 else 1)
        for i,j in enumerate(p): term*=A[i][j]
        total+=term
    return total

def matmul(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0))
             for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return [list(x) for x in zip(*A)]

def pure_gram(active,rho):
    p=len(active); N=2*p
    Z=[[F(0) for _ in range(N)] for _ in range(N)]
    for j,(a,b) in enumerate(active):
        Z[j][2*j]=Z[j][2*j+1]=1
        Z[p+j][2*j]=a; Z[p+j][2*j+1]=b
    H=[[F(0) for _ in range(N)] for _ in range(N)]
    for block in range(2):
        for i in range(p):
            for j in range(p):
                H[block*p+i][block*p+j]=F(1) if i==j else rho
    return matmul(matmul(tr(Z),H),Z)

def psd_by_minors(A):
    n=len(A)
    for k in range(1,n+1):
        for ids in combinations(range(n),k):
            M=[[A[i][j] for j in ids] for i in ids]
            if det_small(M)<0:return False
    return True

def energy_check(word,plus_positions,support):
    negpos=[i for i,x in enumerate(word) if x<0]
    assert len(negpos)==2
    n,m=(-word[i] for i in negpos)
    active=[]
    for mask in support:
        S={negpos[0]}
        S.update(plus_positions[i] for i in range(len(plus_positions))
                 if mask>>i&1)
        allpos=set(range(len(word)))
        a=prod(i+2 for i in S)
        b=prod(i+2 for i in allpos-S)
        assert a!=b
        active.append((a,b))
    if not active:return False
    rho=F(1,n+1) if n==m else F(0)
    mu=[F((a-b)**2,2+a*a+b*b) for a,b in active]
    c=(F(n,n+1) if n==m else F(1))*min(mu)
    G=pure_gram(active,rho)
    Gc=[[G[i][j]-(c if i==j else 0) for j in range(len(G))]
        for i in range(len(G))]
    assert c>0 and psd_by_minors(Gc)
    return True

def d0_dims(word):
    A=tuple(word[::2]); B=tuple(word[1::2])
    fa=signed_table(A); fb=signed_table(B)
    he=ho=0
    for ab in set(fa)|set(fb):
        x=fa.get(ab,0)*fb.get(ab,0)
        he+=max(x,0); ho+=max(-x,0)
    return he,ho

MAXLEN=8
candidate_profiles=0; tc_ns=0; nonzero=0; energy_pass=0
mode_counts=Counter()
for L in range(2,MAXLEN+1):
    profs=0; local_tc=0; local_odd=0
    for mags in combinations_with_replacement(range(1,13),L):
        cnt=Counter(mags)
        if max(cnt.values())>2: continue
        if not no_inv(mags[::2]) or not no_inv(mags[1::2]): continue
        profs+=1
        if profs%5000==0:
            stamp(); print("screen progress length",L,"NS profiles",profs,flush=True)
        choices=[(n,) for n,c in cnt.items() if c==2]
        choices+=list(combinations([n for n,c in cnt.items() if c==1],2))
        for neg in choices:
            word=tuple(-n if n in neg else n for n in mags)
            plus=tuple(x for x in mags if x not in neg)
            plus_positions=tuple(i for i,x in enumerate(word) if x>0)
            n,m=neg if len(neg)==2 else (neg[0],neg[0])
            support=tc_support(plus,n,m)
            if not support: continue
            local_tc+=1; tc_ns+=1
            he,ho=d0_dims(word)
            assert he-ho==signed_table(word).get((0,0),0)
            if not ho: continue
            local_odd+=1; nonzero+=1
            assert energy_check(word,plus_positions,support)
            energy_pass+=1
            mode_counts[tuple(sorted((n,m)))]+=1
    candidate_profiles+=profs
    stamp(); print("energy-screen length",L,"NS magnitude profiles",profs,
                   "TC+NS words",local_tc,"nonzero H_odd",local_odd,
                   "cumulative energy checks",energy_pass,flush=True)
assert energy_pass==nonzero and tc_ns>=nonzero
stamp(); print("energy screen summary",{
    "max_length":MAXLEN,"magnitude_profiles":candidate_profiles,
    "TC_NS_two_minus_words":tc_ns,"nonzero_Hodd":nonzero,
    "exact_Gram_passes":energy_pass,"minus_pair_counts":dict(mode_counts)},
    flush=True)

def family_word(low,M,ell):
    high=tuple(sorted(tuple(M*(2**i) for i in range(ell-1))
                      +(M*(2**(ell-1)-1)-1,)))
    return (-1,-1)+tuple(low)+high,high

def family_checks(low,M):
    r0=fuse(tuple(low)).get(1,0)
    assert r0 in (1,2) and M>sum(low)+2
    checks=0
    for ell in range(2,7):
        word,high=family_word(low,M,ell)
        assert fuse(high).get(1,0)==ell-1
        plus=tuple(abs(x) for x in word if x>0)
        support=tc_support(plus,1,1)
        lowmask=sum(1<<i for i,x in enumerate(plus) if x in low)
        himask=((1<<len(plus))-1)^lowmask
        assert sorted(support)==sorted((lowmask,himask))
        assert no_inv(tuple(abs(x) for x in word[::2]))
        assert no_inv(tuple(abs(x) for x in word[1::2]))
        he,ho=d0_dims(word)
        assert he-ho==signed_table(word).get((0,0),0)
        assert ho>=2*r0*(ell-1)
        K=prod(range(7,ell+7)); b=r0*(ell-1)
        rho=F(1,2)
        active=[(240,3*K),(2*K,360)]
        G=pure_gram(active,rho)
        formula=F(3,4)**(2*b)*F(((3*K-240)*(2*K-360))**(2*b))
        assert det_small(G)**b==formula
        c=F(1,2)*min(F((a-z)**2,2+a*a+z*z) for a,z in active)
        assert psd_by_minors([[G[i][j]-(c if i==j else 0)
                              for j in range(4)] for i in range(4)])
        checks+=1
        stamp(); print("family",low,"ell",ell,"r0",r0,"Hodd",ho,
                       "lower",2*r0*(ell-1),"determinant",formula,flush=True)
    return checks

for q in range(1,13):
    cup_a={}; cup_b={}
    for i in range(q+1):
        for j in range(q+1):
            z=F((-1)**(i+j),q+1)
            cup_a[(i,j,q-i,q-j)]=z
            cup_b[(i,j,q-j,q-i)]=z
    assert sum(z*z for z in cup_a.values())==1
    assert sum(z*z for z in cup_b.values())==1
    overlap=sum(z*cup_b.get(k,F(0)) for k,z in cup_a.items())
    assert overlap==F(1,q+1)
stamp(); print("normalized cup overlaps q=1..12 exact PASS",flush=True)
fam_checks=family_checks((3,5,7),18)+family_checks((3,5,9),20)
assert fam_checks==10

w=(2,2,2,3,3,-4,5,5,5,5,6,6,7,-9)
he,ho=d0_dims(w)
cap=2*fuse(tuple(abs(x) for x in w)).get(0,0)
phi=signed_table(w).get((0,0),0)
assert (he,ho,cap,phi)==(56676672,17240634,15925636,39436038)
assert he-ho==phi and ho-cap==1314998
stamp(); print("Prop4 exact",{"H_even":he,"H_odd":ho,
                              "pure_target_capacity":cap,
                              "kernel_lower_bound":ho-cap,"Phi":phi},flush=True)
print("FM-STR9 independent exact checks PASS",flush=True)
