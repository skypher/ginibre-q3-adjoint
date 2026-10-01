import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

ap=argparse.ArgumentParser(description="FM-MECH109 exact certificates; memory only")
ap.parse_args()
r=(19,66,125,167,190,168,137,78,42,8)
R=sum(r);D=sum(2*j*w for j,w in enumerate(r,1))
S=sum((2*j)**2*w for j,w in enumerate(r,1))
A2=S+2*D+R
assert (R,D,S,A2)==(1000,10386,123140,144912)
K=20;scale=1<<K

P=[[1],[-1,1]]
for j in range(2,11):
    v=[0]*(j+1)
    for i,a in enumerate(P[-1]):v[i]-=2*a;v[i+1]+=a
    for i,a in enumerate(P[-2]):v[i]-=a
    P.append(v)
for j in range(11):
    assert P[j]==[(-1)**(j-k)*comb(j+k,j-k) for k in range(j+1)]

root=tuple(tuple(sum(P[j][h]*4**h*comb(i,h)*factorial(h)*factorial(j-h)
                     for h in range(i+1)) for i in range(j+1))
           for j in range(1,11))

@lru_cache(None)
def bernstein(d,i):
    if d==0:return root
    out=[]
    for j,b in enumerate(bernstein(d-1,i//2),1):
        if i%2==0:
            out.append(tuple(sum(comb(k,h)*b[h] for h in range(k+1))<<(j-k)
                             for k in range(j+1)))
        else:
            out.append(tuple(sum(comb(j-k,h)*b[k+h] for h in range(j-k+1))<<k
                             for k in range(j+1)))
    return tuple(out)

@lru_cache(None)
def ranges(d,i):
    return tuple((min(b),max(b)) for b in bernstein(d,i))

# Independent endpoint checks against the character polynomials.
for d,i in ((0,0),(2,1),(5,21),(9,502),(12,4095)):
    for j,b in enumerate(bernstein(d,i),1):
        den=factorial(j)<<(j*d)
        for index,x in ((0,Q(4*i,1<<d)),(-1,Q(4*(i+1),1<<d))):
            assert Q(b[index],den)==sum(a*x**h for h,a in enumerate(P[j]))

def certificate(negative_only):
    stack=[((0,0),(0,0))]
    seen=leaves=sign_leaves=depth_max=0;area=Q(0)
    while stack:
        A,B=stack.pop();seen+=1
        da,ia=A;db,ib=B;depth=max(da,db)
        depth_max=max(depth_max,depth)
        aa=ranges(da,ia);bb=ranges(db,ib)
        intervals=[]
        for j,((l1,u1),(l2,u2)) in enumerate(zip(aa,bb),1):
            s1=j*(depth-da);s2=j*(depth-db)
            intervals.append(((l1<<s1)-(u2<<s2),
                              (u1<<s1)-(l2<<s2)))
        sign_ok=False
        if negative_only:
            l6,u6=intervals[2];l8,u8=intervals[3]
            sign_ok=(l6>=0 and l8>=0) or (u6<=0 and u8<=0)
        ok=sign_ok
        if not ok:
            num=1
            for j,((lo,hi),w) in enumerate(zip(intervals,r),1):
                den=factorial(j)*(2*j+1)<<(j*depth)
                q=(scale*max(-lo,hi)+den-1)//den
                num*=q**w
            power=64 if negative_only else 1
            ok=(num<<power)<=1<<(K*R)
        if ok:
            leaves+=1;sign_leaves+=sign_ok
            area+=Q(8,1<<(2*da)) if A==B else Q(16,1<<(da+db))
        elif A==B:
            d,i=A;left=(d+1,2*i);right=(d+1,2*i+1)
            stack.extend(((left,left),(right,left),(right,right)))
        elif da<=db:
            stack.extend((((da+1,2*ia),B),((da+1,2*ia+1),B)))
        else:
            stack.extend(((A,(db+1,2*ib)),(A,(db+1,2*ib+1))))
    assert area==8
    return seen,leaves,sign_leaves,depth_max

global_cert=certificate(False)
negative_cert=certificate(True)
assert global_cert==(2161,1083,0,12)
assert negative_cert==(39,22,8,4)
print("Global certificate:",global_cert,"PASS")
print("Negative-region certificate:",negative_cert,"PASS")

H0=Q(1)
for j,w in enumerate(r,1):
    H0*=Q(2*j+1-(-1)**j,2*j+1)**w
assert H0>Q(1,32)
assert 5*S<=4*(D+8)**2
assert 4*(D+8)**2>=1024
assert (1<<116)>1008*32768*(D+8)**6
for ell in (1,2,4,8,16,32,64,128):
    L=ell*D+8
    assert (1<<(116*ell))>1008*32768*L**6
print("Uniform positive-corner comparison: PASS")

# Derivative identity used in the local estimates.
U=[[1],[0,1]]
for n in range(2,21):
    v=[0]+U[-1][:]
    for j,a in enumerate(U[-2]):v[j]-=a
    U.append(v)
for n in range(1,21):
    derivative=[(j+1)*U[n][j+1] for j in range(n)]
    rhs=[0]*n
    for j in range(n-1,-1,-2):
        for k,a in enumerate(U[j]):rhs[k]+=(j+1)*a
    assert derivative==rhs

def fusion(ns):
    row={0:1}
    for n in ns:
        nxt=defaultdict(int)
        for j,a in row.items():
            for h in range(abs(j-n),j+n+1,2):nxt[h]+=a
        row=dict(nxt)
    return row
assert fusion([2,2,6,8])[0]==2

def tensor(word):
    row={(0,0):1}
    for n,e in word:
        nxt=defaultdict(int)
        for (j,k),a in row.items():
            for h in range(abs(j-n),j+n+1,2):nxt[h,k]+=a
            for h in range(abs(k-n),k+n+1,2):nxt[j,h]+=e*a
        row={key:a for key,a in nxt.items() if a}
    return row

for p in (20,22,25):
    B=[(n,-1) for n in (2,2,4,4,6,6,8,8,10,10,6,8)]
    bg=tensor(B+[(p,-1)])
    full=tensor(B+[(p,-1)]*2)
    assert full.get((0,0),0)==2*bg.get((p,0),0)>0
print("Derivative, fusion, and consumer normalization: PASS")

# The integrand of the new family takes both signs.
P3=sum(a*Q(1,2)**j for j,a in enumerate(P[3]))
P4=sum(a*Q(1,2)**j for j,a in enumerate(P[4]))
assert (P[3][0]-P3)*(P[4][0]-P4)==Q(-495,128)
for j in range(1,11):
    assert sum(a*Q(1,2)**h for h,a in enumerate(P[j]))!=P[j][0]

assert 2*A2+130<(D-1)**2
for ell in (1,2,4,8,16,32,64,128):
    p=ell*D-2;delta=ell*D+7
    k=1962*ell+3;N=2038*ell+3
    T=2+Q(2*ell*A2+130,(p+1)**2)
    assert p>=20 and p<delta
    assert delta>=8 and k<40*delta and N<432*(delta+1)**2
    assert T<3

ell=64;p=ell*D-2
bound=Q(145152*ell**2*S**2*(p+1)**2,4**ell)
assert bound<Q(1,10**7)
print("Pair-free witness: ell",ell,"p",p,"delta",ell*D+7,
      "N",2038*ell+3,"k",1962*ell+3)
print("Exact upper ratio:",bound,"< 1/10000000")

def pairfree(word):
    values=set(word)
    return not any((n,-e) in values for n,e in values)
assert pairfree([(2*j,-1) for j in range(1,11)]+[(p,-1)])
assert not pairfree([(1,1),(1,-1),(2,1),(24,-1),(28,1),(52,-1)])
assert Q(9,8)**4/4<Q(1,2)

# Obstruction to making the negative-region bound label-uniform.
a=2048;r_label=2*a;s_label=2*a+2
assert a%2==0
assert r_label+1==2*(a+1)-1
assert s_label+1==2*(a+1)+1
assert (a+1)**2>=20*S
assert Q(1)-Q(10*S,(a+1)**2)>=Q(1,2)
assert s_label>=max(20,r_label)
assert s_label<=D+(r_label+s_label)//2
print("Negative-region uniformity obstruction: r=4096, s=4098 PASS")
print("FM-MECH109 ALL CHECKS PASS")