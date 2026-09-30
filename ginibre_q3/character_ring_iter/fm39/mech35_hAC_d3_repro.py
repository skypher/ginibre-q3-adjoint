"""FM-MECH35 (astra_max_ceres), extracted from its run log: mech35_hAC_d3_repro.py"""
import argparse, random
from fractions import Fraction as Q
from itertools import combinations
from math import comb
argparse.ArgumentParser(description='Exact d=3 certificate for t<=3 or v<=1').parse_args()
def C(n,k): return comb(n,k) if 0<=k<=n else 0
def add(*ps):
    out={}
    for p in ps:
        for x,a in p.items(): out[x]=out.get(x,Q(0))+a
    return out
def scale(p,a): return {x:a*b for x,b in p.items()}
def square(p):
    out={}
    for x,a in p.items():
        for y,b in p.items(): out[x^y]=out.get(x^y,Q(0))+a*b
    return out
def sphere(points,k): return {sum(J):Q(1) for J in combinations(points,k)}
def fusion(mu,n):
    rows=[{0:1}]
    for x in mu:
        for row in rows[:]:
            q={}
            for j,a in row.items():
                for k in range(abs(j-x),j+x+1,2): q[k]=q.get(k,0)+a
            rows.append(q)
    full=len(rows)-1
    return [rows[S].get(0,0)*rows[S^full].get(n,0) for S in range(full+1)]
small={
 (2,2):[(Q(1),(0,15)),(Q(1),(0,3,12,15))],
 (3,2):[(Q(3),(0,22)),(Q(1,2),(0,14,21,27)),(Q(1,2),(0,6,19,21)),
        (Q(1,4),(0,6,11,13,19,21,24,30)),(Q(7,4),(0,3,5,6,24,27,29,30))],
 (3,3):[(Q(1),(0,30,45,51)),(Q(15,2),(0,22,37,51)),(Q(9,2),(0,3,5,6)),
        (Q(1),(0,14,21,27,35,45,54,56)),(Q(3,2),(0,24,40,48)),
        (Q(3,2),(0,6,24,30,40,46,48,54))]}
def code_average(t,v,terms):
    out={}; mask=(1<<t)-1
    for a,H in terms:
        assert a>=0 and all(x^y in H for x in H for y in H)
        counts={}
        for x in H:
            key=((x&mask).bit_count(),(x>>t).bit_count())
            counts[key]=counts.get(key,0)+1
        for x in range(1<<(t+v)):
            i,j=(x&mask).bit_count(),(x>>t).bit_count()
            out[x]=out.get(x,Q(0))+a*Q(counts.get((i,j),0),C(t,i)*C(v,j))
    return out
def certificate(mu):
    L=len(mu); t=mu.count(1); v=mu.count(2); u=mu.count(3); b=L-t
    assert t<=3 or v<=1
    if L==t+v and (t,v) in small: return code_average(t,v,small[t,v])
    I=[1<<i for i,x in enumerate(mu) if x==1]
    II=[1<<i for i,x in enumerate(mu) if x==2]
    III=[1<<i for i,x in enumerate(mu) if x==3]
    E=sphere(I,1); A=sphere(I,2); H=sphere(I,3)
    B=sphere(II,1); K=sphere(III,1); J={i^j:Q(1) for i in I for j in II}
    atoms=[]
    def atom(a,p):
        assert a>=0,(mu,a)
        if a and p: atoms.append((a,p))
    if t<=3:
        for i,j,k in combinations(II,3): atom(Q(1,2),{i:Q(1),j^k:Q(1)})
        if t==0:
            atom(Q(max(L-3,0),2),B); atom(Q(1,2),K)
        elif t==1:
            atom(Q(1,2),add(B,{I[0]^k:Q(1) for k in III}))
            atom(Q(max(L-4,0),2),B)
        elif t==2 and L<=3: pass
        elif t==2 and L==4:
            atom(Q(1,2),add(J,K)); atom(Q(1-v,2),E)
        elif t==3 and L==4: atom(Q(1,2),add(H,K))
        elif t==3 and L==5:
            atom(Q(1,2),add(H,J,K)); atom(Q(2-v,2),E)
        else:
            atom(Q(1,2),add(J,K,H if t==3 else {}))
            atom(Q(L-3-t,2),add(A,B))
            for j in II: atom(Q(1,2),add(A,{j:Q(1)}))
            beta=C(L-2,2)-v if t==2 else C(L-2,2)-L+5-2*v
            atom(Q(beta,2),E)
    elif v==0 and b==0 and t in (6,7):
        even={x:Q(1) for x in range(1<<t) if x.bit_count()%2==0}
        if t==6:
            atom(Q(2,len(even)),even); atom(Q(3,2),{0:Q(1),(1<<t)-1:Q(1)})
        else:
            atom(Q(4,len(even)),even)
            atom(Q(1,2),add(E,{(1<<t)-1:Q(1)}))
    elif v==1 and b==1 and t in (4,5):
        # Change the distinguished coordinate from a label 0 or 1 to the label 2.
        if t==4:
            p={x:Q(1) for x in range(1<<L) if ((x&15)^(15 if x>>4 else 0)).bit_count()==1}
            z={x:Q(1) for x in range(1<<L) if ((x&15)^(15 if x>>4 else 0))==0}
            atom(Q(1,4),p); atom(Q(1,2),z)
        else:
            def transform(x):
                z=x>>5
                return ((x&31)^(31 if z else 0))|(z<<5)
            for k,a in [(2,Q(1,3)),(1,Q(1,6))]:
                atom(a,{x:Q(1) for x in range(1<<L) if transform(x).bit_count()==k})
    else:
        # A nonnegative 3-by-3 Gram factorization; sample a single label-3 point.
        samples=[{k:Q(u)} for k in III] if u else [{}]
        for K0 in samples:
            norm=Q(1,len(samples))
            atom(norm/144,add(scale(H,5),scale(J,6),scale(K0,12)))
            atom(norm/18,add(H,scale(J,3)))
            atom(norm/48,add(H,scale(K0,4)))
        atom(Q(1,2),K)
        alpha=Q(t+4*b-8,12)
        h=Q(t+4*b-10,4)
        if v: atom(alpha,add(A,scale(B,h/(2*alpha))))
        else: atom(alpha,A)
        gamma=Q(6*b*b+4*b*t-14*b+t*t-7*t-6*v*v-12*v+10,12)
        atom(gamma/2,E)
    D=C(L+1,3)-t*(L-1)-v
    rho=Q(D)-sum(a*sum(z*z for z in p.values()) for a,p in atoms)
    atom(rho,{0:Q(1)})
    out={}
    for a,p in atoms:
        for x,z in square(p).items(): out[x]=out.get(x,Q(0))+a*z
    return out
cases=[]
for L in range(1,11):
    for t in range(L+1):
        for v in range(L-t+1):
            if t>3 and v>1: continue
            for u in range(L-t-v+1):
                mu=(1,)*t+(2,)*v+(3,)*u+(7,)*(L-t-v-u)
                if sum(mu)>=6: cases.append(mu)
rng=random.Random(3535)
for _ in range(100):
    L=rng.randrange(4,13)
    if rng.randrange(2):
        t=rng.randrange(4); rest=tuple(rng.randrange(2,10) for _ in range(L-t))
    else:
        t=rng.randrange(L+1); v=rng.randrange(min(1,L-t)+1)
        rest=(2,)*v+tuple(rng.randrange(3,10) for _ in range(L-t-v))
    mu=tuple(sorted((1,)*t+rest))
    if sum(mu)>=6: cases.append(mu)
entries=0
for mu in cases:
    f=fusion(mu,sum(mu)-6); out=certificate(mu)
    assert all(out.get(S,0)==a for S,a in enumerate(f)),mu
    entries+=len(f)
print('d=3, t<=3 or v<=1:',len(cases),'lists;',entries,'exact entries; PASS')
for mu in [(1,)*8,(1,)*6+(2,),(2,)*4,(1,)*3+(3,)]:
    f=fusion(mu,sum(mu)-6); out=certificate(mu)
    assert all(out.get(S,0)==a for S,a in enumerate(f))
    print('boundary',mu,'n',sum(mu)-6,'origin',f[0],'PASS')