import argparse, pathlib, re, gzip
ap=argparse.ArgumentParser(description="FM-STR7e exact prefix verifier; memory-only compilation.")
ap.add_argument("--threads",type=int,default=24)
ap.add_argument("--census",default=str(pathlib.Path(__file__).resolve().parent/"sec166_census_w40_noflip.log.gz"))
args=ap.parse_args()
pat=re.compile(r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+) D=")
census={}
_cp=pathlib.Path(args.census)
for line in (gzip.decompress(_cp.read_bytes()).decode().splitlines(True) if _cp.suffix=='.gz' else _cp.open()):
 m=pat.match(line)
 if m and int(m[1])<=40:
  w=tuple(map(int,m[3].split()))+(int(m[2]),)
  assert sum(map(abs,w[:-1]))==int(m[1])
  census[w]=int(m[4])
assert len(census)==5430
jobs=[(0,w,p) for w,p in sorted(census.items(),key=lambda x:(sum(map(abs,x[0][:-1])),x[0]))]
for k in range(1,21):
 W=k*(k+1)//2
 for p in range(k+1,k+11):
  if p<max(k,6) or (W-p)%2 or (W-p)//2<max(k,8):continue
  jobs.append((1,tuple(-n for n in range(1,k+1))+((-1 if k%2 else 1)*p,),-1))
for mask in range(32):
 B=[]
 for j,(n,count) in enumerate([(1,2),(2,1),(3,4),(4,3),(5,4)]):
  B += [(-n if mask&(1<<j) else n)]*count
 sigma=-1 if sum(z<0 for z in B)%2 else 1
 jobs.append((2,tuple(B)+(sigma*8,),-1))
jobs.append((3,(-1,2,3,4,-5,6,7,8),956))
assert len(jobs)==5531
data="".join(str(g)+" "+str(len(w))+" "+str(p)+" "+" ".join(map(str,w))+"\n" for g,w,p in jobs)

CPP=r'''
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;
using W=std::vector<int>;
struct Row{
 int d;std::vector<Z> c;
 Row(int n=0):d(n),c((n+1)*(n+1)){}
 Z get(int u,int v)const{
  if(u<0||v<0||u+v>d)return 0;return c[u*(d+1)+v];
 }
 Z& at(int u,int v){return c[u*(d+1)+v];}
};
void stamp(const std::string&s){
 std::time_t t=std::time(nullptr);std::tm u;gmtime_r(&t,&u);
 std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
int weight(const W&w){int t=0;for(int z:w)t+=std::abs(z);return t;}
Row mul(const Row&f,int z,int cap=-1){
 int n=std::abs(z),e=z>0?1:-1,D=f.d+n;if(cap>=0)D=std::min(D,cap);
 Row g(D);
 for(int axis=0;axis<2;++axis)
 for(int v=0;v<=std::min(f.d,D);++v){
  int last=f.d-v;std::vector<Z> ps(last+3);
  for(int u=0;u<=last;++u)ps[u+2]=ps[u]+(axis?f.get(v,u):f.get(u,v));
  for(int t=0;t+v<=D;++t){
   int lo=std::abs(t-n),hi=std::min(last,t+n);
   if((hi-lo)%2)--hi;if(lo>hi)continue;
   Z q=ps[hi+2]-ps[lo];
   if(axis)g.at(v,t)+=e*q;else g.at(t,v)+=q;
  }
 }return g;
}
Row row(W w,int target=-1){
 std::sort(w.begin(),w.end(),[](int a,int b){return std::abs(a)>std::abs(b);});
 int rem=weight(w);Row f;f.at(0,0)=1;
 for(int z:w){rem-=std::abs(z);f=mul(f,z,target<0?-1:rem+target);}
 return f;
}
W omit(const W&w,int i,int j){
 W C;for(int k=0;k<(int)w.size();++k)if(k!=i&&k!=j)C.push_back(w[k]);return C;
}
std::pair<W,W> split(W C){
 std::sort(C.begin(),C.end(),[](int a,int b){return std::abs(a)>std::abs(b);});
 W A,B;int a=0,b=0;
 for(int z:C)if(a<=b){A.push_back(z);a+=std::abs(z);}else{B.push_back(z);b+=std::abs(z);}
 if(a>b)std::swap(A,B);return {A,B};
}
std::vector<Z> inner_levels(const W&w,int i,int j){
 auto[A,B]=split(omit(w,i,j));B.push_back(w[i]);B.push_back(w[j]);
 Row a=row(A),b=row(B,weight(A));std::vector<Z>L(a.d+1);
 for(int u=0;u<=a.d;++u)for(int v=0;u+v<=a.d;++v)L[u+v]+=a.get(u,v)*b.get(u,v);
 return L;
}
std::vector<Z> terminal_levels(const W&w,int i,int j,Row*save=nullptr){
 int a=std::abs(w[i]),b=std::abs(w[j]);Row C=row(omit(w,i,j),a+b);
 std::vector<Z>L(a+b+1);
 for(int c=std::abs(a-b);c<=a+b;c+=2)L[c]+=2*C.get(c,0);
 L[a+b]+=2*(w[j]>0?1:-1)*C.get(a,b);
 if(save)*save=std::move(C);return L;
}
Z total(const std::vector<Z>&L){Z x=0;for(const Z&z:L)x+=z;return x;}
int bad(const std::vector<Z>&L){Z x=0;for(int t=0;t<(int)L.size();++t){x+=L[t];if(x<0)return t;}return -1;}
std::string wordstr(const W&w){std::ostringstream s;for(int z:w)s<<z<<",";return s.str();}
std::string levstr(const std::vector<Z>&L){std::ostringstream s;for(int t=0;t<(int)L.size();++t)if(L[t]!=0)s<<t<<":"<<L[t]<<",";return s.str();}

struct Job{int group;W w;Z expected;};
struct Stat{
 long long words=0,noflip=0,pairs=0,endpointbad=0,midbad=0,purebad=0,budgetbad=0,allbudgetbad=0;
 long long sp4pairs=0,sp4bad=0,lowerpairs=0,lowerbad=0,lowerallbad=0,tpbad=0;
 long long negpure=0,posmixed=0,prefixes=0;
 void add(const Stat&s){
  words+=s.words;noflip+=s.noflip;pairs+=s.pairs;endpointbad+=s.endpointbad;midbad+=s.midbad;
  purebad+=s.purebad;budgetbad+=s.budgetbad;allbudgetbad+=s.allbudgetbad;sp4pairs+=s.sp4pairs;
  sp4bad+=s.sp4bad;lowerpairs+=s.lowerpairs;lowerbad+=s.lowerbad;lowerallbad+=s.lowerallbad;
  tpbad+=s.tpbad;negpure+=s.negpure;posmixed+=s.posmixed;prefixes+=s.prefixes;
 }
};
std::vector<Z> sp4_levels(const Row&C,int a,int b){
 if(a>b)std::swap(a,b);std::vector<Z>L(a+b-1);
 for(int j=0;j<a;++j)for(int beta=0;beta<a-j;++beta){
  int alpha=a+b-2-2*j-beta,u=alpha+1,v=beta;
  Z z=C.get(u-1,v)+C.get(u+1,v)-C.get(u,v-1)-C.get(u,v+1);
  L[alpha+beta]+=2*z;
 }
 std::vector<Z>P(L.size());Z p=0;
 for(int t=0;t<(int)L.size();++t){p+=L[t];P[t]=p;}
 for(int j=0;j<a;++j){
  Z child=-C.get(a-j,b-j);
  for(int c=b-a;c<=a+b-2*j;c+=2)child+=C.get(c,0);
  assert(P[a+b-2-2*j]==2*child);
 }return L;
}
std::string witness(const W&w,int i,int j,int t,const Z&minimum,
 const std::vector<Z>&P,const std::vector<Z>&M){
 std::ostringstream s;s<<"word="<<wordstr(w)<<" pair="<<w[i]<<","<<w[j]
 <<" first_T="<<t<<" minimum="<<minimum<<" pure="<<levstr(P)<<" mixed="<<levstr(M);return s.str();
}
Stat test(const Job&job,std::string&fail,std::string&bandfail){
 const W&w=job.w;int n=w.size(),p=std::abs(w.back()),Wb=weight(w)-p,mx=0,cores=0,minus=0;
 std::map<int,int> signs;
 for(int k=0;k<n;++k){
  int a=std::abs(w[k]),e=w[k]>0?1:-1;assert(a>0);
  if(signs.count(a))assert(signs[a]==e);signs[a]=e;minus+=w[k]<0;
  if(k<n-1){mx=std::max(mx,a);cores+=a>=3;}
 }
 assert(minus%2==0&&(Wb-p)%2==0&&p>=std::max(6,mx)&&cores>=2);
 assert((Wb-p)/2>=8&&mx<=(Wb-p)/2);
 Z phi=row(w,0).get(0,0);
 if(job.expected>=0)assert(phi==job.expected);assert(phi>=0);
 Stat out;out.words=1;bool nf=true,anybudget=false,anylower=false;int ti=-1,tj=-1;
 for(int i=0;i<n-1;++i)for(int j=i+1;j<n-1;++j)if((std::abs(w[i])+std::abs(w[j]))%2==0)
  if(ti<0||std::make_pair(std::abs(w[i])+std::abs(w[j]),std::max(std::abs(w[i]),std::abs(w[j])))>
           std::make_pair(std::abs(w[ti])+std::abs(w[tj]),std::max(std::abs(w[ti]),std::abs(w[tj])))){ti=i;tj=j;}
 Z best;bool have=false;std::string bestline;
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
  int a=std::abs(w[i]),b=std::abs(w[j]),e=w[j]>0?1:-1;
  Row C;auto E=terminal_levels(w,i,j,&C);assert(total(E)==phi);++out.pairs;
  Z D=e*C.get(a,b);nf&=D<0;
  for(int c=std::abs(a-b);c<=a+b;c+=2)assert(C.get(c,0)>=0);
  out.endpointbad+=bad(E)>=0;
  if(w[i]<0&&w[j]<0){
   auto S=sp4_levels(C,a,b);assert(total(S)==phi);++out.sp4pairs;out.sp4bad+=bad(S)>=0;
  }
  bool eligible=!(a==1&&w[i]<0)&&!(b==1&&w[j]<0);
  if(eligible){
   Z delta=2*(C.get(a+b,0)+e*(C.get(a,b)-C.get(a-1,b-1)));
   assert(phi-delta>=0);++out.lowerpairs;out.lowerbad+=delta<0;anylower|=delta>=0;
   if(i==ti&&j==tj)out.tpbad+=delta<0;
   if(w[i]<0&&w[j]<0&&delta<0&&bandfail.empty()){
    std::ostringstream q;q<<"word="<<wordstr(w)<<" pair="<<w[i]<<","<<w[j]<<" phi="<<phi
     <<" lowered_phi="<<phi-delta<<" top_band="<<delta<<" Sp4_levels="<<levstr(sp4_levels(C,a,b));
    bandfail=q.str();
   }
  }
  auto[A,B]=split(omit(w,i,j));int cap=weight(A);Row x=row(A);
  W V=B;V.push_back(w[i]);V.push_back(w[j]);Row y=row(V,cap);
  V=B;V.push_back(-w[i]);V.push_back(-w[j]);Row z=row(V,cap);
  std::vector<Z>P(cap+1),M(cap+1),L(cap+1);
  Z all=0,pure=0,budget=0,minimum=0;int first=-1;bool mb=false,pb=false,np=false,pm=false;
  for(int t=0;t<=cap;++t){
   for(int u=0;u<=t;++u){
    Z q=x.get(u,t-u)*y.get(u,t-u),f=x.get(u,t-u)*z.get(u,t-u);
    assert((q+f)%2==0);P[t]+=(q+f)/2;M[t]+=(q-f)/2;
   }
   L[t]=P[t]+M[t];all+=L[t];pure+=P[t];budget+=P[t]+(M[t]<0?M[t]:Z(0));
   mb|=all<0;pb|=pure<0;np|=P[t]<0;pm|=M[t]>0;++out.prefixes;
   if(budget<0&&first<0)first=t;if(budget<minimum)minimum=budget;
  }
  assert(all==phi&&total(M)==2*D);
  out.midbad+=mb;out.purebad+=pb;out.negpure+=np;out.posmixed+=pm;out.budgetbad+=first>=0;
  anybudget|=first<0;
  if(!have||minimum>best){have=true;best=minimum;bestline=witness(w,i,j,first,minimum,P,M);}
 }
 out.noflip=nf;if(job.group==0)assert(nf);
 out.allbudgetbad=!anybudget;out.lowerallbad=!anylower;
 if(!anybudget)fail=bestline;
 return out;
}
Row character(int a,int b){
 Row f(a+b);for(int j=0;j<=b;++j)for(int k=0;k<=a-b;++k)f.at(j+k,j+a-b-k)+=1;return f;
}
void module_checks(){
 Row f=character(1,1);for(int j=0;j<4;++j)f=mul(f,-1);assert(f.get(0,0)==-6);
 long long checks=0;
 for(int a=0;a<=4;++a)for(int b=0;b<=a;++b)
 for(int c=0;c<=4;++c)for(int d=0;d<=c;++d)for(int r=0;r<=2;++r){
  Row x=character(a,b),y=character(c,d);
  for(int j=0;j<r;++j)x=mul(x,-1);for(int j=0;j<2-r;++j)y=mul(y,-1);
  Z p=0;for(int t=0;t<=x.d;++t){
   for(int u=0;u<=t;++u)p+=x.get(u,t-u)*y.get(u,t-u);
   assert(p>=0);++checks;
  }
 }
 stamp("two genuine lift prefix checks="+std::to_string(checks)+"; d^4 chi(1,1)=-6");
}
void receipt(){
 W w{-1,2,3,4,-5,6,7,8};int i=4,j=6;auto[A,B]=split(omit(w,i,j));int cap=weight(A);
 Row x=row(A);W V=B;V.push_back(w[i]);V.push_back(w[j]);Row y=row(V,cap);
 V=B;V.push_back(-w[i]);V.push_back(-w[j]);Row z=row(V,cap);
 std::vector<Z>P(cap+1),M(cap+1),stock(cap+1);Z budget=0,fixedMixed=0;
 std::map<std::pair<int,int>,Z>flow;
 std::cout<<"RECEIPT TopPair=(-5,7) A="<<wordstr(A)<<" B="<<wordstr(B)<<std::endl;
 for(int t=0;t<=cap;++t){
  for(int u=0;u<=t;++u){
   Z q=x.get(u,t-u)*y.get(u,t-u),f=x.get(u,t-u)*z.get(u,t-u);
   P[t]+=(q+f)/2;M[t]+=(q-f)/2;
  }
  if(P[t]>0)stock[t]+=P[t];
  Z need=(P[t]<0?-P[t]:Z(0))+(M[t]<0?-M[t]:Z(0));
  for(int s=0;s<=t&&need>0;++s){
   Z take=std::min(stock[s],need);stock[s]-=take;need-=take;flow[{t,s}]+=take;
  }
  assert(need==0);if(M[t]>0)fixedMixed+=M[t];
  budget+=P[t]+(M[t]<0?M[t]:Z(0));assert(budget>=0);
  if(P[t]!=0||M[t]!=0)std::cout<<"LEVEL "<<t<<" pure="<<P[t]<<" mixed="<<M[t]<<" budget="<<budget<<std::endl;
 }
 Z fixedPure=0;for(auto&v:stock)fixedPure+=v;
 assert(fixedPure+fixedMixed==956);
 for(auto&[ts,c]:flow)if(c!=0)std::cout<<"MATCH negative_height="<<ts.first<<" pure_height="<<ts.second<<" count="<<c<<std::endl;
 std::cout<<"FIXED pure="<<fixedPure<<" mixed="<<fixedMixed<<" total=956"<<std::endl;
}
int main(int argc,char**argv){
 int threads=24;
 for(int i=1;i<argc;++i){
  std::string a=argv[i];
  if(a=="-h"||a=="--help"){std::cout<<"FM-STR7e --threads N\n";return 0;}
  if(a=="--threads"&&i+1<argc)threads=std::stoi(argv[++i]);else return 2;
 }
 assert(threads>0);omp_set_num_threads(threads);
 std::vector<Job>jobs;int g,n;Z phi;
 while(std::cin>>g>>n>>phi){W w(n);for(int&z:w)std::cin>>z;jobs.push_back({g,w,phi});}
 assert(jobs.size()==5531);stamp("FM-STR7e begin; cpp_int arithmetic");module_checks();
 std::vector<std::string>fail(jobs.size()),band(jobs.size());Stat sums[4];std::atomic<int>done{0};
 #pragma omp parallel
 {
  Stat local[4];
  #pragma omp for schedule(dynamic,1)
  for(int i=0;i<(int)jobs.size();++i){
   local[jobs[i].group].add(test(jobs[i],fail[i],band[i]));
   int d=++done;if(d%500==0){
    #pragma omp critical
    stamp("processed="+std::to_string(d)+"/"+std::to_string(jobs.size()));
   }
  }
  #pragma omp critical
  for(int i=0;i<4;++i)sums[i].add(local[i]);
 }
 for(int g=0;g<4;++g){
  auto&s=sums[g];
  std::cout<<"GROUP "<<g<<" words="<<s.words<<" noflip="<<s.noflip<<" pairs="<<s.pairs
  <<" endpoint_bad="<<s.endpointbad<<" midpoint_bad="<<s.midbad<<" pure_prefix_bad="<<s.purebad
  <<" pure_target_budget_bad="<<s.budgetbad<<" all_budget_bad_words="<<s.allbudgetbad
  <<" Sp4_minus_pairs="<<s.sp4pairs<<" Sp4_bad="<<s.sp4bad<<" eligible_lower_pairs="<<s.lowerpairs
  <<" lowering_bad="<<s.lowerbad<<" all_lower_bad_words="<<s.lowerallbad<<" TopPair_lower_bad="<<s.tpbad
  <<" negative_pure_level_pairs="<<s.negpure<<" positive_mixed_level_pairs="<<s.posmixed
  <<" prefixes="<<s.prefixes<<std::endl;
 }
 int printed=0;for(auto&s:fail)if(!s.empty()&&printed++<3)std::cout<<"ALL_RESTRICTED_FAIL "<<s<<std::endl;
 for(int g=0;g<4;++g)for(int i=0;i<(int)jobs.size();++i)if(jobs[i].group==g&&!band[i].empty()){
  std::cout<<"NEGATIVE_TOP_BAND group="<<g<<" "<<band[i]<<std::endl;break;
 }
 for(auto&s:sums){assert(s.endpointbad==0&&s.midbad==0&&s.sp4bad==0&&s.purebad==0&&s.budgetbad==0);}
 receipt();
 stamp("FM-STR7e PASS");
}
'''
import os,shlex,subprocess,sys
os.environ["TMPDIR"]="/dev/shm"
obj=os.memfd_create("fmstr7e_obj",0)
exe=os.memfd_create("fmstr7e_exe",0)
subprocess.run(["g++","-Werror=return-type","-O2","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}","--threads",str(args.threads)],input=data,text=True,pass_fds=(exe,),check=True)
