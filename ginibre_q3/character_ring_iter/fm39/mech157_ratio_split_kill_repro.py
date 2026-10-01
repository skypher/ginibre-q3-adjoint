import argparse
from collections import defaultdict
from fractions import Fraction
from math import comb

argparse.ArgumentParser(
    description="FM-MECH157: exact F-star counterexample and ratio-split obstruction."
).parse_args()

def subset_data(ns):
    Z=1<<len(ns);weight=[0]*Z;size=[0]*Z
    for s in range(1,Z):
        bit=s&-s;i=bit.bit_length()-1
        weight[s]=weight[s^bit]+ns[i]
        size[s]=size[s^bit]+1
    return weight,size

def invariants(ns):
    assert all(n>0 and n%2==0 for n in ns)
    weight,size=subset_data(ns);Z=len(weight)
    m=[0]*Z;m[0]=1
    for s in range(1,Z):
        ell=size[s]
        if ell==1:continue
        half=weight[s]//2;j=s;v=0
        while True:
            q=half-weight[j]-size[j]
            if q>=0:
                v+=(-1 if size[j]%2 else 1)*comb(q+ell-2,ell-2)
            if j==0:break
            j=(j-1)&s
        assert v>=0
        m[s]=v
    return m

def transform(m):
    Z=len(m);full=Z-1
    F=[m[s]*m[full^s] for s in range(Z)]
    h=1
    while h<Z:
        for a in range(0,Z,2*h):
            for j in range(a,a+h):
                x,y=F[j],F[j+h]
                F[j],F[j+h]=x+y,x-y
        h*=2
    return F

def fusion_rows(ns):
    rows=[None]*(1<<len(ns));rows[0]=(1,)
    for s in range(1,len(rows)):
        bit=s&-s;i=bit.bit_length()-1;h=ns[i]//2
        old=rows[s^bit];pref=[0]
        for v in old:pref.append(pref[-1]+v)
        row=[]
        for k in range(len(old)+h):
            lo=abs(k-h);hi=min(k+h,len(old)-1)
            row.append(pref[hi+1]-pref[lo] if lo<=hi else 0)
        rows[s]=tuple(row)
    return rows

def coefficient(rows,s,n=0):
    if n%2 or n//2>=len(rows[s]):return 0
    return rows[s][n//2]

def phi_subset(rows,s,negative):
    ans=0;j=s
    while True:
        sign=-1 if (j&negative).bit_count()%2 else 1
        ans+=sign*coefficient(rows,j)*coefficient(rows,s^j)
        if j==0:break
        j=(j-1)&s
    return ans

def mixed(rows,remaining,negative,a,b):
    ans=0;j=remaining
    while True:
        sign=-1 if (j&negative).bit_count()%2 else 1
        ans+=sign*coefficient(rows,remaining^j,a)*coefficient(rows,j,b)
        if j==0:break
        j=(j-1)&remaining
    return ans

ns=tuple(range(40,64,2));L=len(ns);full=(1<<L)-1;negative=5
rows=fusion_rows(ns);m=invariants(ns);F=transform(m)
assert all(rows[s][0]==m[s] for s in range(1<<L))
print("independent invariant identities:",len(m))

signed=tuple(-n if negative>>i&1 else n for i,n in enumerate(ns))
B=signed[:-1];p=signed[-1];W=sum(map(abs,B));delta=(W-p)//2
pair=max(((i,j) for i in range(L-1) for j in range(i+1,L-1)
          if (ns[i]-ns[j])%2==0),
         key=lambda ij:(ns[ij[0]]+ns[ij[1]],ns[ij[1]]))
tp=ns[pair[0]]+ns[pair[1]]
assert len(set(ns))==L and sum(z<0 for z in B)%2==0
assert p>=max(6,max(map(abs,B)))
assert delta>=max(8,max(map(abs,B)))
assert sum(abs(z)>=3 for z in B)>=2
assert (W,p,delta,tp)==(550,62,244,118) and 2*tp<delta
assert pair==(9,10)
parent=F[negative]
assert parent==phi_subset(rows,full,negative)==906415068445728
margins=[]
pairs=[(i,j) for i in range(L) for j in range(i+1,L)]
for i,j in pairs:
    bits=(1<<i)|(1<<j)
    A=mixed(rows,full^bits,negative,ns[i],ns[j])
    margin=(-1 if negative>>j&1 else 1)*A
    assert parent-F[negative^bits]==4*margin
    assert margin<0
    margins.append((margin,ns[i],ns[j]))
assert min(margins)==(-4506106071735,40,44)
assert max(margins)==(-34349665,42,44)
child=phi_subset(rows,full^(1<<9)^(1<<10),negative)
assert child==359292699210
assert parent%2==child%2==0 and child<parent
print("finite residual: W,p,delta,TopPair =",W,p,delta,tp)
print("Phi =",parent,"all negative flip margins =",len(margins))
print("margin minimum/maximum =",min(margins),max(margins))
print("g parent/TopPair child =",parent//2,child//2)

# For n_i=2(t+i), every binomial argument has form
# (|S|-2|J|)t + sum(S)-2sum(J)-|J|.
# Its constant term has absolute value at most 66+12=78.
# Thus the contributing terms are fixed for t>=80.
# Flip differences omit S=empty/full and have degree at most 12-4=8.
T0=1024;mask=6
assert T0>sum(range(L))+L==78

def signed_subsum(m,selected,negative):
    out=0;s=selected
    while True:
        sign=-1 if (s&negative).bit_count()%2 else 1
        out+=sign*m[s]*m[selected^s]
        if s==0:break
        s=(s-1)&selected
    return out

values=[];drops=[];children=[]
for t in range(T0,T0+11):
    mt=invariants(tuple(2*(t+i) for i in range(L)))
    Ft=transform(mt)
    child=signed_subsum(mt,full^(1<<9)^(1<<10),mask)
    assert Ft[mask]%2==child%2==0
    drops.append((Ft[mask]-child)//2)
    children.append(child//2)
    row=[]
    for i,j in pairs:
        diff=Ft[mask^(1<<i)^(1<<j)]-Ft[mask]
        assert diff%4==0
        row.append(diff//4)
    values.append(row)

certificates=[]
for col in range(len(pairs)):
    d=[row[col] for row in values[:9]];c=[]
    while d:
        c.append(d[0])
        d=[d[i+1]-d[i] for i in range(len(d)-1)]
    assert c[0]>0 and all(v>=0 for v in c)
    assert c[7]==c[8]==0
    certificates.append(c)
assert min(c[0] for c in certificates)==5822451256163514
assert min(v for c in certificates for v in c if v>0)==885
assert certificates[0]==[
    5822451256163514,31047417002688,131148758741,
    412024173,856711,885,0,0,0]
print("uniform family: t >= 1024, negative labels at indices 1,2")
print("nonnegative Newton coefficients =",sum(len(c) for c in certificates))
print("minimum constant =",min(c[0] for c in certificates))
print("minimum positive coefficient =",min(v for c in certificates for v in c if v>0))

for t in (T0+9,T0+17,T0+80):
    Ft=transform(invariants(tuple(2*(t+i) for i in range(L))))
    for (i,j),c in zip(pairs,certificates):
        diff=Ft[mask^(1<<i)^(1<<j)]-Ft[mask]
        assert diff==4*sum(v*comb(t-T0,h) for h,v in enumerate(c))
print("additional polynomial identity checks =",3*len(pairs))

# Parent and child have degrees at most 10 and 8, respectively.
for name,data in (("TopPair drop",drops),("TopPair child",children)):
    c=[];d=data.copy()
    while d:
        c.append(d[0]);d=[d[i+1]-d[i] for i in range(len(d)-1)]
    assert c[0]>0 and all(v>=0 for v in c)
    print("uniform",name,"Newton certificate:",len(c),"nonnegative coefficients")

def entry(labels,a=0,b=0):
    labels=sorted(labels,key=abs,reverse=True)
    rem=sum(map(abs,labels));F={(0,0):1}
    for z in labels:
        n=abs(z);sg=1 if z>0 else -1;rem-=n;G=defaultdict(int)
        for (i,j),v in F.items():
            for k in range(abs(i-n),i+n+1,2):
                if abs(k-a)+abs(j-b)<=rem:G[k,j]+=v
            for k in range(abs(j-n),j+n+1,2):
                if abs(i-a)+abs(k-b)<=rem:G[i,k]+=sg*v
        F={ij:v for ij,v in G.items() if v}
    return F.get((a,b),0)

B=[-1]*13+[-2]+[-3]*5+[-4];p=6
C=B.copy();C.remove(-2);C.remove(-4)
g=entry(B,p);gc=entry(C,p)
assert (g,gc)==(700607800,723462700)
assert entry(B+[p])==2*g and entry(C+[p])==2*gc
assert sum(map(abs,B))==34 and Fraction(6,14)==Fraction(3,7)
assert max(map(abs,B))<=14 and sum(z<0 for z in B)%2==0
assert gc>g>0
print("TopPair failure at ratio 3/7: parent/child =",g,gc)
print("PASS")
