from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement
from random import Random
from math import comb
import os, shlex, subprocess

def stamp(s):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
          s, flush=True)

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def charrow(ns):
    r={0:1}
    for n in ns:
        q={}
        for a,v in r.items():
            for b in cg(a,n):
                q[b]=q.get(b,0)+v
        r=q
    return r

@lru_cache(None)
def inv(ns):
    return charrow(ns).get(0,0)

def phi_direct(word):
    d={(0,0):1}
    for z in word:
        n=abs(z)
        eps=1 if z>0 else -1
        q=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n):
                q[c,b]+=v
            for c in cg(b,n):
                q[a,c]+=eps*v
        d={k:v for k,v in q.items() if v}
    return d.get((0,0),0)

def walsh_values(labels):
    L=len(labels)
    full=(1<<L)-1
    m=[inv(tuple(sorted(labels[i] for i in range(L) if s>>i&1)))
       for s in range(1<<L)]
    f=[m[s]*m[full^s] for s in range(1<<L)]
    h=1
    while h<1<<L:
        for st in range(0,1<<L,2*h):
            for j in range(st,st+h):
                x,y=f[j],f[j+h]
                f[j],f[j+h]=x+y,x-y
        h*=2
    return f

def pairfree_words(labels):
    groups=[]
    for n in labels:
        if not groups or groups[-1][0]!=n:
            groups.append([n,1])
        else:
            groups[-1][1]+=1
    for mask in range(1<<len(groups)):
        yield tuple((n if mask>>j&1 else -n)
                    for j,(n,c) in enumerate(groups) for _ in range(c))

# Exact interval certificates for H^(7) and H^(8).
def coefficient_profile(K):
    even=(K==8)
    lens=range(K+1,3*(K+1),2) if even else range(K+1,2*(K+1))
    hs=range(1,3*(K+1),2) if even else range(1,2*(K+1))
    best=[None]*9
    shapes=0
    for l in lens:
        for h in hs:
            for off in range(1-h,l):
                for parity in ((0,) if even else (0,1)):
                    vmin=max(0,K-h+1)
                    u=max(0,vmin-2*off)
                    u+=(parity-u)%2
                    I=range(u,u+2*l,2)
                    J=range(u+2*off,u+2*(off+h),2)
                    b0=len(set(I)&set(J))
                    if not b0:
                        continue
                    shapes+=1
                    for s in range(9):
                        b=sum(abs(i-j)<=2*s<=i+j for i in I for j in J)
                        x=Fraction(b,b0)
                        best[s]=x if best[s] is None else min(best[s],x)
    return shapes,tuple(best)

H7=(Fraction(1),Fraction(2),Fraction(3),Fraction(4),Fraction(5),
    Fraction(43,8),Fraction(41,8),Fraction(9,2),Fraction(7,2))
H8=(Fraction(1),Fraction(2),Fraction(3),Fraction(4),Fraction(5),
    Fraction(6),Fraction(6),Fraction(17,3),Fraction(5))
assert coefficient_profile(7)==(4440,H7)
stamp("H^(7) interval certificate: 4,440 shapes PASS")
assert coefficient_profile(8)==(3393,H8)
stamp("H^(8) interval certificate: 3,393 shapes PASS")

# Layer identity against direct bivariate character multiplication.
def layer_check(word):
    labels=tuple(abs(z) for z in word)
    L=len(labels)
    full=(1<<L)-1
    mm=[inv(tuple(sorted(labels[i] for i in range(L) if s>>i&1)))
        for s in range(1<<L)]
    M=sum((1<<i) for i,z in enumerate(word) if z<0)
    assert M.bit_count()%2==0
    phi=sum((-1 if (s&M).bit_count()%2 else 1)*mm[s]*mm[full^s]
            for s in range(1<<L))
    N=mm[full]
    P=T3p=T3m=T4p=T4m=0
    for s in range(1,full):
        k=s.bit_count()
        if k not in (2,3,4):
            continue
        v=mm[s]*mm[full^s]
        sg=-1 if (s&M).bit_count()%2 else 1
        if k==2:
            P+=sg*v
        elif k==3:
            if sg>0:
                T3p+=v
            else:
                T3m+=v
        else:
            if sg>0:
                T4p+=v
            else:
                T4m+=v
    assert phi==2*(N+P+T3p-T3m+T4p-T4m)
    return phi

rng=Random(106164)
for _ in range(60):
    labels=tuple(rng.randint(1,12) for _ in range(9))
    if sum(labels)%2:
        labels=labels[:-1]+(labels[-1]+1,)
    signs=[1 if rng.randrange(2) else -1 for _ in labels]
    if sum(z<0 for z in signs)%2:
        signs[0]*=-1
    w=tuple(n*e for n,e in zip(labels,signs))
    assert layer_check(w)==phi_direct(w)
stamp("layer identity: 60 subset-sum/direct-bivariate checks PASS")

# Opposite-sign equal-label character identity.
for n in range(1,31):
    assert charrow((n,n))=={2*j:1 for j in range(n+1)}
stamp("pair reduction character identity: n=1..30 PASS")

# Two-odd fusion with an even-label background; j=0 uses the stated convention.
for _ in range(80):
    C=tuple((1 if rng.randrange(2) else -1)*2*rng.randint(1,8)
            for _ in range(4))
    a=2*rng.randint(1,8)+1
    b=2*rng.randint(1,8)+1
    ea=1 if rng.randrange(2) else -1
    eb=1 if rng.randrange(2) else -1
    lhs=phi_direct(C+(ea*a,eb*b))
    rhs=0
    for j in cg(a,b):
        if j==0:
            rhs+=2*phi_direct(C) if ea*eb>0 else 0
        else:
            rhs+=phi_direct(C+((ea*eb)*j,))
    assert lhs==rhs
stamp("two-odd fusion identity: 80 direct checks PASS")

# All nine-factor profiles on labels 1..6, all indexed sign assignments.
stamp("start bounded-label exhaustive sweep")
profiles=signings=0
minimum=None
for P in combinations_with_replacement(range(1,7),9):
    V=walsh_values(P)
    profiles+=1
    signings+=len(V)
    assert min(V)>=0,(P,min(V))
    minimum=min(V) if minimum is not None else min(V)
    if profiles%500==0:
        stamp(f"bounded profiles={profiles}/2002, indexed signings={signings}")
assert (profiles,signings,minimum)==(2002,1025024,0)
stamp("bounded-label sweep PASS: 2,002 profiles / 1,025,024 values")

# 120 random profiles on labels 7..48, all sign assignments per profile.
stamp("start random larger-label sweep")
rng=Random(164106)
random_profiles=random_values=0
for _ in range(120):
    P=tuple(sorted(rng.randint(7,48) for _ in range(9)))
    V=walsh_values(P)
    random_profiles+=1
    random_values+=len(V)
    assert min(V)>=0,(P,min(V))
    if random_profiles%20==0:
        stamp(f"random profiles={random_profiles}/120, signed values={random_values}")
    if random_profiles<=32:
        mask=rng.randrange(512)
        w=tuple(-n if mask>>i&1 else n for i,n in enumerate(P))
        assert phi_direct(w)==V[mask]
assert (random_profiles,random_values)==(120,61440)
stamp("random sweep PASS: 120 profiles / 61,440 values; 32 direct bridges")

def compatible_pairfree(labels):
    for w in pairfree_words(labels):
        if sum(labels)%2==0 and sum(z<0 for z in w)%2==0:
            yield w

base=(1,1,1,1,2,2,8,8,8)
same=sorted((phi_direct(w),w) for w in compatible_pairfree(base))
assert same==[
    (578,(-1,-1,-1,-1,-2,-2,8,8,8)),
    (578,(1,1,1,1,-2,-2,8,8,8)),
    (1226,(-1,-1,-1,-1,2,2,8,8,8)),
    (1226,(1,1,1,1,2,2,8,8,8))]
print("minimum-profile signings:",same,flush=True)

neighbors={}
for i,n in enumerate(base):
    for dn in (-2,2):
        if n+dn<1:
            continue
        q=list(base)
        q[i]+=dn
        q=tuple(sorted(q))
        p=q[-1]
        delta=(sum(q)-2*p)//2
        if (sum(q)%2 or p<6 or delta<8 or q[-2]>delta or
                sum(x>=3 for x in q[:-1])<2):
            continue
        neighbors[q]=sorted(
            (phi_direct(w),w) for w in compatible_pairfree(q))
expected={
    (1,1,1,1,2,4,8,8,8):
        [1030,1030,1130,1130,1570,1570,1694,1694],
    (1,1,1,2,2,3,8,8,8):
        [828,828,1028,1028,1268,1268,1788,1788]}
assert {q:[v for v,w in vv] for q,vv in neighbors.items()}==expected
assert len(neighbors)==2 and all(len(v)==8 for v in neighbors.values())
for q,vv in sorted(neighbors.items()):
    print("admissible +/-2 neighbor:",q,"values:",[v for v,w in vv],
          flush=True)
stamp("minimum and full admissible one-coordinate neighborhood PASS")

# Sixteen direct evaluations sampled from the finite residual box.
rng=Random(164106)
checked=trials=0
stamp("start exact finite-box direct sample")
while checked<16 and trials<200000:
    trials+=1
    small=sorted(rng.randint(1,6) for _ in range(5))
    D=sum(small)
    q=rng.randint(small[-1],2*D)
    r=rng.randint(q,q+D)
    s=rng.randint(r,q+D)
    lo=max(s,2*q-D)
    hi=min(q+D,D+q+r-s)
    if lo>hi:
        continue
    p=rng.randint(lo,hi)
    a=small+[q,r,s,p]
    T=sum(a)
    delta=(T-2*p)//2
    if (T%2 or p<6 or delta<8 or s>delta or
            sum(x>=3 for x in a[:8])<2):
        continue
    if q>=D+(p-q)+1:
        continue
    groups=[]
    for n in a:
        if not groups or groups[-1][0]!=n:
            groups.append([n,1])
        else:
            groups[-1][1]+=1
    for _ in range(100):
        bits=[rng.randrange(2) for _ in groups]
        w=tuple((n if bits[j] else -n)
                for j,(n,c) in enumerate(groups) for _ in range(c))
        if sum(z<0 for z in w)%2==0:
            break
    else:
        continue
    assert phi_direct(w)>=0,(a,w,phi_direct(w))
    checked+=1
    if checked%4==0:
        stamp(f"direct finite-box sample={checked}/16; candidates={trials}")
assert checked==16
stamp("finite-box direct sample PASS: 16 words; maximum label <=23")

# Separated-quartet boundary: pair-partition children can have seven factors.
C=(2,1,1,1,1)
quad=(10,11,12,13)
D=sum(C)
q,r,s,p=quad
assert q==D+(p-q)+1 and p+s<=D+q+r
assert 1 in cg(10,11) and 1 in cg(12,13) and 1+1<=D
assert len(C)+2==7
stamp("quartet boundary: 7-factor child exists; all-seven result covers it")

# Full exact finite-box census. The inclusion-exclusion terms and all
# accumulated values are guarded below 2^55; products use 128-bit temporaries.
CPP=r'''
#include <bits/stdc++.h>
using namespace std; using Z=long long; using Wide=__int128_t;
int a[9],pc[512],sw[512],sp[512]; Z ch[203][8];
Z words=0,signs=0,minimum=LLONG_MAX; array<int,9>witness;
Z inv(int m){
 int k=pc[m],s=sw[m];
 if(!k)return 1;
 if(s%2||k==1)return 0;
 if(k==2){
  int i=__builtin_ctz((unsigned)m);
  int j=__builtin_ctz((unsigned)(m&(m-1)));
  return a[i]==a[j];
 }
 if(k==3){
  int ma=0;
  for(int i=0;i<9;i++)if(m>>i&1)ma=max(ma,a[i]);
  return 2*ma<=s;
 }
 if(k==4){
  int b[4],j=0;
  for(int i=0;i<9;i++)if(m>>i&1)b[j++]=a[i];
  return max(0,(min(b[0]+b[1],b[2]+b[3])-
                max(abs(b[0]-b[1]),abs(b[2]-b[3])))/2+1);
 }
 Z out=0;
 for(int t=m;;t=(t-1)&m){
  int top=s/2-sp[t]+k-2;
  if(top>=k-2){
   assert(top<=202);
   out+=(pc[t]%2?-1:1)*ch[top][k-2];
  }
  assert(abs(out)<(1LL<<55));
  if(!t)break;
 }
 assert(out>=0);
 return out;
}
void progress(){
 time_t now=time(nullptr);
 tm* u=gmtime(&now);
 char b[32];
 strftime(b,sizeof(b),"%Y-%m-%d %H:%M:%S UTC",u);
 cout<<b<<" finite patterns="<<words
     <<" compatible_signings="<<signs<<endl;
}
void word(){
 int T=accumulate(a,a+9,0);
 if(T%2)return;
 int p=a[8],delta=T/2-p;
 if(p<6||delta<8||a[7]>delta)return;
 int big=0;
 for(int i=0;i<8;i++)big+=a[i]>=3;
 if(big<2)return;
 words++;
 for(int m=1;m<512;m++){
  int i=__builtin_ctz((unsigned)m);
  sw[m]=sw[m&(m-1)]+a[i];
  sp[m]=sw[m]+pc[m];
 }
 Z mm[512];
 for(int m=0;m<512;m++)mm[m]=inv(m);
 Z f[512]={};
 f[0]=mm[511];
 for(int m=1;m<511;m++)if(pc[m]>=2&&pc[m]<=4){
  Wide v=Wide(mm[m])*mm[511^m];
  assert(v>=0&&v<(Wide(1)<<55));
  f[m]=Z(v);
 }
 for(int h=1;h<512;h*=2)
  for(int st=0;st<512;st+=2*h)
   for(int j=st;j<st+h;j++){
    Z x=f[j],y=f[j+h];
    assert(Wide(abs(x))+abs(y)<(Wide(1)<<55));
    f[j]=x+y;
    f[j+h]=x-y;
   }
 int groups[9]={},ng=0;
 for(int i=0;i<9;i++){
  if(!i||a[i]!=a[i-1])ng++;
  groups[ng-1]|=1<<i;
 }
 for(int sg=0;sg<(1<<ng);sg++){
  int m=0;
  for(int i=0;i<ng;i++)if(sg>>i&1)m|=groups[i];
  if(pc[m]%2)continue;
  signs++;
  assert(f[m]>=0);
  if(2*f[m]<minimum){
   minimum=2*f[m];
   for(int i=0;i<9;i++)witness[i]=(m>>i&1?-a[i]:a[i]);
  }
 }
 if(words%500000==0)progress();
}
int main(){
 for(int m=0;m<512;m++)pc[m]=__builtin_popcount((unsigned)m);
 for(int n=0;n<=202;n++){
  ch[n][0]=1;
  for(int k=1;k<=7;k++)ch[n][k]=n?ch[n-1][k-1]+ch[n-1][k]:0;
 }
 for(a[0]=1;a[0]<=6;a[0]++)
 for(a[1]=a[0];a[1]<=6;a[1]++)
 for(a[2]=a[1];a[2]<=6;a[2]++)
 for(a[3]=a[2];a[3]<=6;a[3]++)
 for(a[4]=a[3];a[4]<=6;a[4]++){
  int D=a[0]+a[1]+a[2]+a[3]+a[4];
  for(a[5]=a[4];a[5]<=2*D;a[5]++)
  for(a[6]=a[5];a[6]<=a[5]+D;a[6]++)
  for(a[7]=a[6];a[7]<=a[5]+D;a[7]++)
  for(a[8]=max(a[7],2*a[5]-D);
      a[8]<=min(a[5]+D,D+a[5]+a[6]-a[7]);a[8]++){
   assert(a[8]<=90);
   word();
  }
 }
 assert(words==3049184&&signs==160310765&&minimum==578);
 cout<<"finite census "<<words<<" patterns "<<signs
     <<" signings minPhi "<<minimum<<endl;
 cout<<"minimum word";
 for(int x:witness)cout<<" "<<x;
 cout<<endl;
}
'''
obj=os.memfd_create("chk106_obj",0)
exe=os.memfd_create("chk106_exe",0)
subprocess.run(
    ["g++","-O3","-std=c++17","-pipe","-x","c++","-","-c","-o",
     f"/proc/self/fd/{obj}"],
    input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
    ["g++","-###","-fno-use-linker-plugin",f"/proc/self/fd/{obj}",
     "-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines()
          if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
os.close(obj)
os.close(exe)
stamp("FM-CHK106 complete exact verifier PASS")
