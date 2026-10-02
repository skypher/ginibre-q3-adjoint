import argparse,ctypes,os,subprocess
from math import comb
from itertools import combinations
ap=argparse.ArgumentParser(description="FM-STR1: exact intertwiners and locality obstruction; memory only.")
ap.add_argument("--threads",type=int,default=4)
args=ap.parse_args()
assert 1<=args.threads<=32

def parts(n,lo=1):
 if not n:
  yield ()
 for a in range(lo,n+1):
  for z in parts(n-a,a):yield (a,)+z

def polynomial(a,S):
 # Variables: u0,u1,v0,v1,z10,z11,...
 L=len(a);P={(0,)*(4+2*L):1}
 for i,n in enumerate(a):
  e=2 if i in S else 0
  Q={}
  for j in range(n+1):
   ex=[0]*(4+2*L)
   ex[e]=n-j;ex[e+1]=j;ex[4+2*i]=j;ex[5+2*i]=n-j
   Q[tuple(ex)]=(-1)**j*comb(n,j)
  R={}
  for x,c in P.items():
   for y,d in Q.items():
    z=tuple(i+j for i,j in zip(x,y));R[z]=R.get(z,0)+c*d
  P={x:c for x,c in R.items() if c}
 return P

columns=0
for N in range(2,7):
 for a in parts(N):
  for m in range(1,N):
   seen=set()
   for mask in range(1<<len(a)):
    S={i for i in range(len(a)) if mask>>i&1}
    if sum(a[i] for i in S)!=m:continue
    P=polynomial(a,S)
    for flip in (0,1):
     out={}
     for ex,c in P.items():
      for v in range(0,len(ex),2):
       src=v+1-flip;dst=v+flip
       if ex[src]:
        y=list(ex);y[src]-=1;y[dst]+=1;y=tuple(y)
        out[y]=out.get(y,0)+c*ex[src]
     assert not any(out.values()) # both sl2 generators annihilate h_S
    special={}
    for ex,c in P.items():
     if ex[1]==ex[2]==0:
      special[ex[4:]]=special.get(ex[4:],0)+c
    row=tuple(x for i,n in enumerate(a) for x in ((n,0) if i in S else (0,n)))
    assert special=={row:(-1)**m} and row not in seen
    seen.add(row);columns+=1
print("polynomial invariant columns and diagonal minors",columns,"PASS",flush=True)

from itertools import combinations,product
from fractions import Fraction as Q
def eps(v):
 if len(set(v))<3:return 0
 return (-1)**sum(v[i]>v[j] for i in range(3) for j in range(i+1,3))
cols=[(0,)+s for s in combinations(range(1,6),2)]
M=[]
for v in product(range(3),repeat=6):
 row=[]
 for S in cols:
  T=tuple(i for i in range(6) if i not in S)
  row.append(Q(eps(tuple(v[i] for i in S))*eps(tuple(v[i] for i in T))))
 assert row[3]==row[0]-row[1]+row[2]
 assert row[6]==-row[0]-row[4]+row[5]
 assert row[7]==-row[0]+row[1]-row[2]-row[4]+row[5]
 assert row[8]==-row[0]-row[2]+row[5]
 assert row[9]==row[0]-row[1]+row[4]
 if any(row):M.append(row)
row=0;piv=[]
for c in range(10):
 i=next((i for i in range(row,len(M)) if M[i][c]),None)
 if i is None:continue
 M[row],M[i]=M[i],M[row];q=M[row][c];M[row]=[z/q for z in M[row]]
 for i in range(len(M)):
  if i!=row and M[i][c]:
   q=M[i][c];M[i]=[x-q*y for x,y in zip(M[i],M[row])]
 piv.append(c);row+=1
assert row==5 and piv==[0,1,2,4,5]
points=[(0,0,1,1,2,2),(0,0,1,2,1,2),(0,1,0,1,2,2),
        (0,1,0,2,1,2),(0,1,2,0,1,2)]
minor=[]
for v in points:
 minor.append([eps(tuple(v[i] for i in cols[c]))*
               eps(tuple(v[i] for i in range(6) if i not in cols[c]))
               for c in piv])
assert minor==[[0,0,0,0,1],[0,0,0,1,0],[0,0,1,0,0],
               [0,1,0,0,0],[1,0,0,0,1]] # determinant 1
print('six adjoints: ten triple-product tensors have rank',row,'PASS',flush=True)

source=r'''
#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <map>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <atomic>
#include <omp.h>
using namespace std;
using I=boost::multiprecision::cpp_int;
using V=vector<int>;
I binom(int n,int k){
 if(k<0||k>n)return 0;I a=1;
 for(int i=1;i<=k;i++){a*=n-i+1;a/=i;}return a;
}
I cat(int n){return binom(2*n,n)/(n+1);}
I inv(const V&ns){
 vector<I>v(1,I(1));
 for(int n:ns){
  vector<I>w(v.size()+n+2);
  for(int s=0;s<(int)v.size();s++)if(v[s]!=0){w[abs(s-n)]+=v[s];w[s+n+2]-=v[s];}
  for(int s=2;s<(int)w.size();s++)w[s]+=w[s-2];
  w.resize(w.size()-2);v.swap(w);
 }return v[0];
}
struct Graded{
 int d;
 vector<I>v;
 const I& at(int i,int j,int e)const{return v[(i*d+j)*2+e];}
};
Graded table(const V&word){
 int W=0;for(int z:word){assert(z!=0&&abs(z)<=64);W+=abs(z);}
 assert(W<=160&&word.size()<=60);
 int d=W+1,t=0;vector<I>v(2*d*d),w(2*d*d);v[0]=1;
 for(int z:word){
  int n=abs(z);fill(w.begin(),w.end(),I(0));
  for(int i=0;i<=t;i++)for(int j=0;j<=t-i;j++)for(int e=0;e<2;e++){
   const I&c=v[(i*d+j)*2+e];if(c==0)continue;
   for(int q=abs(i-n);q<=i+n;q+=2)w[(q*d+j)*2+e]+=c;
   for(int q=abs(j-n);q<=j+n;q+=2)w[(i*d+q)*2+(e^(z<0))]+=c;
  }t+=n;v.swap(w);
 }return {d,move(v)};
}
I oddhom(const Graded&g){
 I n=0;for(int i=0;i<g.d;i++)for(int j=0;j<g.d;j++)n+=2*g.at(i,j,0)*g.at(i,j,1);return n;
}
void partitions(int n,int lo,V&a,vector<V>&out){
 if(!n){out.push_back(a);return;}
 for(int j=lo;j<=n;j++){a.push_back(j);partitions(n-j,j,a,out);a.pop_back();}
}
extern "C" int run(int threads){
 assert(1<=threads&&threads<=32);omp_set_num_threads(threads);
 struct Case{int n,m;V a;};vector<Case>cases;
 for(int n=1;n<=8;n++)for(int m=n;m<=8;m++){
  V a;vector<V>ps;partitions(n+m,1,a,ps);
  for(auto&p:ps)cases.push_back({n,m,p});
 }
 assert(cases.size()<=10000);
 atomic<int>done{0};atomic<long long>minors{0};
 #pragma omp parallel for schedule(dynamic,4)
 for(int id=0;id<(int)cases.size();id++){
  auto q=cases[id];V w=q.a;w.push_back(-q.n);w.push_back(-q.m);
  auto g=table(w);I k=0;V all=q.a;all.push_back(q.n);all.push_back(q.m);
  assert(q.a.size()<=16);unsigned end=1u<<q.a.size();
  for(unsigned s=0;s<end;s++){
   int sum=0;for(int i=0;i<(int)q.a.size();i++)if(s>>i&1)sum+=q.a[i];
   if(sum==q.m)k++;
  }
  assert(g.at(0,0,1)==2*k);
  // Independent single-SU(2) check of the all-one-colour target.
  I h=inv(all);assert(h>=k);
  assert(g.at(0,0,0)>=2*h);
  assert(g.at(0,0,0)-g.at(0,0,1)>=0);
  assert(k>=0&&k<=65536); // Hence the total counter is below 10000*65536.
  minors+=k.convert_to<long long>();done++;
 }
 cout<<"boundary words "<<done<<"; selected minor columns "<<minors<<" PASS"<<endl;
 int calibration=0;
 for(int n=1;n<=12;n++)for(int m=n;m<=12;m++){
  int N=n+m;V w(N,1);w.push_back(-n);w.push_back(-m);
  auto g=table(w);I value=0;
  for(int j=1;j<=n;j++)value+=2*binom(N,2*j)*cat(j)*binom(N-2*j,n-j);
  assert(g.at(0,0,0)-g.at(0,0,1)==value);
  V all(N,1);all.push_back(n);all.push_back(m);
  assert(inv(all)==binom(N,n));
  if(n==m&&(n==2||n==5||n==10))
   cout<<"boundary n=m="<<n<<" factors="<<N+2<<" dim-="<<g.at(0,0,1)
       <<" dim+="<<g.at(0,0,0)<<" Phi="<<value<<endl;
  calibration++;
 }
 cout<<"positive binomial formula "<<calibration<<" PASS"<<endl;
 int pairchecks=0;
 for(int a=1;a<=8;a++)for(int b=1;b<=8;b++)for(int x:{-1,1})for(int y:{-1,1}){
  if(a==b&&x!=y)continue;
  assert(oddhom(table({x*a,y*b}))==0);pairchecks++;
 }
 cout<<"all pair-free two-slot cases "<<pairchecks<<" PASS"<<endl;
 int local=0;
 for(int n=2;n<=24;n++){
  int h0=(n+1)/2+2;
  for(int t=1;t<=2;t++)for(int s=0;s<=n;s++){
   V w(s,1);w.insert(w.end(),t,-n);auto g=table(w);I h=oddhom(g);
   if(t==1){if(s<n)assert(h==0);else assert(h>0);}
   if(s+t<h0)assert(h==0);
   if(t==2&&s+t==h0){
    assert(h>0);
    int x=n%2?n-1:n,y=(n+1)/2;
    assert(g.at(x,y,0)>0&&g.at(x,y,1)>0);
   }local++;
  }
 }
 cout<<"sharp locality threshold cases "<<local<<" PASS"<<endl;
 auto six=table({-2,-2,-2,-2,-2,-2});
 assert(six.at(0,0,0)==120&&six.at(0,0,1)==20);
 cout<<"six-adjoint transfer obstruction: source=20 image<=10 even=120 Phi=100"<<endl;
 auto channel=table({3,3,3,3,-2,-2});
 assert(channel.at(0,0,0)==92&&channel.at(0,0,1)==12);
 for(int k=0;k<=4;k++)assert(3*k!=2);
 cout<<"maximal-channel obstruction (+3)^4(-2)^2: even=92 odd=12 Phi=80"<<endl;
 int fail=0;
 for(int n=1;n<=8;n++)for(int m=n;m<=8;m++){
  int N=n+m+2;V w(N,1);w.push_back(-n);w.push_back(-m);auto g=table(w);
  I r=(m+1)*binom(N,n)+(n+1)*binom(N,n+2);
  I h=binom(N,n+1)-1;V all(N,1);all.push_back(n);all.push_back(m);
  assert(g.at(0,0,1)==2*r&&inv(all)==h&&r>h);
  assert(g.at(0,0,0)>=g.at(0,0,1));
  if(n==2&&m==2)cout<<"first-excess witness (+1)^6(-2)^2: source="<<2*r
    <<" pure-colour target="<<2*h<<" total even="<<g.at(0,0,0)
    <<" Phi="<<g.at(0,0,0)-g.at(0,0,1)<<endl;
  fail++;
 }
 cout<<"first-excess dimension obstructions "<<fail<<" PASS"<<endl;
 cout<<"FM-STR1 PASS"<<endl;return 0;
}
'''
env=dict(os.environ,TMPDIR="/dev/shm")
assembly=subprocess.run(["g++","-O2","-pipe","-fPIC","-fopenmp","-S",
 "-x","c++","-","-o","-"],input=source.encode(),stdout=subprocess.PIPE,
 env=env,check=True).stdout
obj=os.memfd_create("str1_obj");lib=os.memfd_create("str1_lib")
subprocess.run(["as","-o",f"/proc/self/fd/{obj}"],input=assembly,
 pass_fds=(obj,),env=env,check=True)
def locate(n):
 return subprocess.check_output(["g++","-print-file-name="+n],text=True,env=env).strip()
subprocess.run(["ld","-shared","--eh-frame-hdr",
 "-L"+os.path.dirname(locate("libgcc_s.so")),"-o",f"/proc/self/fd/{lib}",
 locate("crtbeginS.o"),f"/proc/self/fd/{obj}",locate("libstdc++.so"),
 locate("libgcc_s.so"),locate("libgomp.so"),"-lc",locate("crtendS.o")],
 pass_fds=(obj,lib),env=env,check=True)
mod=ctypes.CDLL(f"/proc/self/fd/{lib}");mod.run.argtypes=[ctypes.c_int]
assert mod.run(args.threads)==0
