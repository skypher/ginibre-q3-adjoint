import argparse,itertools,os,pathlib,re,shlex,subprocess
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb

ap=argparse.ArgumentParser(description="Exact FM-MECH155 certificates; memory only.")
ap.add_argument("--census-log")
args=ap.parse_args()

CPP=r'''
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using Z=long long;
struct Test {vector<int>w; int c; bool exceptional;};
void check(int family) {
 vector<int>L,H; vector<Test>T; int cap,minimum;
 if(family==0) {
  L=H={5,6,7,8,9};cap=22;minimum=4;
  T={{{1,2,3,4,3,2,1},56,false}};
 } else if(family==4) {
  L=H={4,6,8,10};cap=22;minimum=2;
  T={{{1,6,2,1},28,true},{{1,2,4,2,1},28,true},
     {{1,4,3,3,2,1},28,false},{{1,2,3,3,3,2,1},28,false},
     {{1,2,3,4,3,2,1},30,false},{{1,4,3,2,1},30,true},
     {{1,2,3,3,2,1},30,false},{{1,3,5,3,2,1},30,false}};
 } else {
  L={4,6,8,10};H={3,5,7,9,11,13};cap=25;minimum=2;
  T={{{1,4,3,2,1},28,true},{{1,2,3,3,2,1},28,true},
     {{1,3,5,3,2,1},28,false},{{1,2,3,4,3,2,1},28,false}};
 }
 Z count=0,bad=0,exceptions=0;
 #pragma omp parallel for reduction(+:count,bad,exceptions) schedule(dynamic)
 for(int index=0;index<(int)L.size();index++)for(int h:H) {
  int l=L[index];
  for(int x=-l-h+2;x<=cap;x++)
  for(int y=-l+1;y<=cap;y++)
  for(int z=-h+1;z<=cap;z++) {
   int u=x+y,v=x+z,p=y+z;
   if(u<0||v<0||p<minimum)continue;
   if(family && ((u|v|p)&1))continue;
   Z b[7]={};
   for(int i=0;i<l;i++)for(int j=0;j<h;j++)
   for(int s=0;s<=6;s++) {
    int lo=max({-s,s-p,-y-i+j,-z+i-j});
    int hi=min(s,x+i+j);
    b[s]+=max(0,hi-lo+1);
   }
   if(!b[0])continue;
   ++count;
   bool ex=family && l==4 && h==(family==4?4:3)
           && u==0 && v==0 && p==2;
   exceptions+=ex;
   for(auto &t:T) {
    Z a=0;for(int s=0;s<(int)t.w.size();s++)a+=t.w[s]*b[s];
    // Each b[s] <= 13*169; all products here are < 10^6.
    assert(b[0]<=169 && a>=0 && a<1000000);
    if(a<t.c*b[0] && !(ex&&t.exceptional))++bad;
   }
  }
 }
 Z expected=family==0?472650:family==4?79384:170108;
 assert(count==expected && bad==0 && exceptions==(family!=0));
 cout<<"tiles family="<<family<<" cases="<<count
     <<" exceptions="<<exceptions<<" PASS\n"<<flush;
}
int main(int argc,char**argv) {
 if(argc>1 && (string(argv[1])=="-h"||string(argv[1])=="--help")) {
  cout<<"FM155 bounded endpoint verifier; no arguments.\n";return 0;
 }
 check(0);check(4);check(2);
}
'''

obj=os.memfd_create("fm155_object",0)
exe=os.memfd_create("fm155_executable",0)
subprocess.run(["g++","-O3","-std=c++17","-fopenmp","-pipe","-x","c++","-",
                "-c","-o",f"/proc/self/fd/{obj}"],
               input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",
                  f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
                 text=True,capture_output=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
os.close(obj);os.close(exe)

D=(1,2,3,4,3,2,1)
V={"E222":(1,6,2,1),"E224":(1,2,4,2,1),
   "E2aa":(1,4,3,3,2,1),"E2dist":(1,2,3,3,3,2,1),
   "M233":(1,4,3,2,1),"M2dist":(1,2,3,3,2,1),
   "M334":(1,3,5,3,2,1),"D":D}
high1=(1,4,3,3,3,2,1)
high2=(1,2,3,3,3,3,2)
trials=0
for g,h,k in itertools.product(range(8),repeat=3):
    ns=tuple(sorted((g+h,g+k,h+k)))
    if ns[0]<2:continue
    b=[max(0,min(s,g)-max(-s,s-h-k,-h,-k)+1) for s in range(7)]
    aug=[0]*7
    if ns==(2,2,2):key="E222";aug[1]=3
    elif ns==(2,2,4):key="E224";aug[2]=1
    elif ns[0]==2 and ns[1]==ns[2]:
        key="M233" if ns[1]%2 else "E2aa";aug[1]=1
        if ns[1]%2 and ns[1]>=5:
            assert all(x+y>=z for x,y,z in zip(b,aug,high1))
    elif ns[0]==2:
        assert ns[2]==ns[1]+2
        key="M2dist" if ns[1]%2 else "E2dist"
        if ns[1]%2 and ns[1]>=5:
            assert all(x>=z for x,z in zip(b,high2))
    elif ns==(3,3,4):key="M334";aug[2]=1
    else:key="D"
    assert all(b[s]+aug[s]>=v for s,v in enumerate(V[key]))
    if ns[0]>=4:assert all(x>=y for x,y in zip(b,D))
    trials+=1
print("triple slack cases",trials,"PASS")

E=tuple(map(F,(1,F(5,2),3,2,1,F(1,2),0)))
M=tuple(map(F,(1,F(8,3),F(10,3),3,2,1,F(1,3))))
O=tuple(map(F,(1,F(8,3),F(11,3),F(11,3),3,2,1)))
def bands(I,J):
    return [sum(abs(a-b)<=2*s<=a+b for a in I for b in J)
            for s in range(7)]
for r,parities,profile in ((2,(0,1),E),(2,(1,),M),(3,(0,1),O)):
    for d in range(1,7):
        for parity in parities:
            if d<=r:
                t=r+1-d;u=2*t+parity
                pairs=[(range(u-2*t,u+2*d,2),
                        range(u,u+2*(d+t),2))]
            else:
                u=2 if parity==0 else 1
                pairs=[(range(u,u+2*d,2),range(u,u+2*(d+1),2))]
                u=2 if parity==0 else 3
                pairs.append((range(u,u+2*d,2),
                              range(u-2,u+2*d,2)))
            for I,J in pairs:
                b=bands(I,J)
                assert b[0]==d
                assert all(x>=d*y for x,y in zip(b,profile))
    assert all(F(2*s+1)-F(3*s*s-s,14)>=v
               for s,v in enumerate(profile))
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
assert dot(E,O)==30 and dot(M,M)>=30 and dot(M,O)>=30
assert dot(O,O)==49
print("quartet endpoints and large-d formula PASS")

for o in (0,4,6):
    odd=(1<<o)-1
    for sg in range(256):
        m=sg.bit_count()
        if m%2:continue
        counts=[]
        for size in (3,4):
            c=0
            for S in itertools.combinations(range(8),size):
                mask=sum(1<<i for i in S)
                if size==4 and not mask&1:continue
                c+=(mask&odd).bit_count()%2==0 and (mask&sg).bit_count()%2==1
            counts.append(c)
        c3,c4=counts
        if o==0:assert F(c3,56)+F(c4,49)<=1
        else:
            C=28 if o==4 else 30
            assert c3<=C
            if c4:assert c3<=16 and F(c3,C)+F(c4,30)<=1
print("all sign-parity budgets PASS")
b2=(7,18,21,17,10,4,1);b4=(9,24,30,27,19,10,4)
assert dot(high1,b2)>=28*b2[0] and dot(high2,b2)>=28*b2[0]
assert dot(high1,b4)>=30*b4[0]

@lru_cache(None)
def inv(ns):
    n=len(ns);s=sum(ns)
    if not n:return 1
    if s%2 or 2*max(ns)>s:return 0
    if n==1:return int(ns[0]==0)
    out=0
    for mask in range(1<<n):
        a=s//2+n-2-sum(ns[i]+1 for i in range(n) if mask>>i&1)
        if a>=n-2:out+=(-1)**mask.bit_count()*comb(a,n-2)
    return out
def cuts(ns):
    for size in (2,3,4):
        for S in itertools.combinations(range(8),size):
            mask=sum(1<<i for i in S)
            if size==4 and not mask&1:continue
            a=tuple(ns[i] for i in range(8) if mask>>i&1)
            b=tuple(ns[i] for i in range(8) if not mask>>i&1)
            yield mask,size,inv(a)*inv(b)
def stat(word):
    ns=tuple(map(abs,word));sg=sum(1<<i for i,z in enumerate(word) if z<0)
    q=Counter(N=inv(ns),P=0,n3=0,p3=0,n4=0,p4=0,d3=0)
    for mask,size,v in cuts(ns):
        neg=(mask&sg).bit_count()%2
        if size==2:q["P"]+=(-1)**neg*v
        else:q[("n" if neg else "p")+str(size)]+=v
        if size==3 and neg:q["d3"]=max(q["d3"],v)
    q["phi"]=2*(q["N"]+q["P"]+q["p3"]+q["p4"]-q["n3"]-q["n4"])
    return q
def direct(word):
    A={(0,0):1}
    for z in word:
        T=defaultdict(int);n=abs(z);sg=1 if z>0 else -1
        for (a,b),v in A.items():
            for c in range(abs(a-n),a+n+1,2):T[c,b]+=v
            for c in range(abs(b-n),b+n+1,2):T[a,c]+=sg*v
        A={ab:v for ab,v in T.items() if v}
    return A.get((0,0),0)

exceptions=[((2,)*4+(3,)*4,780),
            ((2,)*3+(3,)*4+(4,),764),
            ((2,)*4+(3,)*3+(5,),596),
            ((2,)*2+(3,)*6,1096)]
for ns,minimum in exceptions:
    classes=sorted(set(ns));values=[]
    for mask in range(1<<len(classes)):
        w=tuple(-n if mask>>classes.index(n)&1 else n for n in ns)
        if sum(z<0 for z in w)%2:continue
        v=stat(w)["phi"];assert v==direct(w);values.append(v)
    assert len(values)==4 and min(values)==minimum
    print("exception",ns,"minimum",minimum,"PASS")

for w,phi,lo,hi in [
    ((-1,-2,-3,-4,-4,-5,-5,-6),718,-332,-4),
    ((2,4,-6,8,10,-12,14,16),15660,-2860,-32)]:
    q=stat(w);assert q["phi"]==phi==direct(w)
    diffs=[]
    for i,j in itertools.combinations(range(8),2):
        v=list(w);v[i]*=-1;v[j]*=-1
        diffs.append(phi-stat(tuple(v))["phi"])
    assert (min(diffs),max(diffs))==(lo,hi)
    if abs(w[0])==1:
        assert q["N"]==404 and q["P"]==43 and q["d3"]==20
        assert 28*q["d3"]-q["N"]-q["P"]==113
    print("remaining-sector witness",w,dict(q),"flip range",lo,hi)

if args.census_log:
    lengths=Counter();covered=0;remaining=0
    for line in pathlib.Path(args.census_log).read_text().splitlines():
        m=re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+)",line)
        if not m:continue
        w=tuple(map(int,m[3].split()))+(int(m[2]),)
        lengths[len(w)]+=1
        if len(w)!=8:continue
        ns=tuple(map(abs,w));o=sum(n%2 for n in ns)
        assert stat(w)["phi"]==int(m[4])
        good=(min(ns)>=3 or (min(ns)>=2 and o>0))
        covered+=good;remaining+=not good
    assert lengths=={6:499,7:1597,8:1944,9:1066,10:270,11:52,12:2}
    assert (covered,remaining)==(306,1638)
    print("census eight covered",covered,"remaining",remaining,
          "combined <=7 plus this result",2096+covered,"/5430 PASS")
print("FM-MECH155 PARTIAL THEOREMS PASS")
