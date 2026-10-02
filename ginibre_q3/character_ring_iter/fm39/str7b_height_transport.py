import os,shlex,subprocess,sys
if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
 print("FM-STR7b exact verifier: --threads N (default 24).")
 sys.exit(0)
os.environ["TMPDIR"]="/dev/shm"
CPP=r'''
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;
using W=std::vector<int>;
struct Row {
 int d; std::vector<Z> c;
 Row(int D=0):d(D),c((D+1)*(D+1)){}
 Z get(int a,int b)const {
  if(a<0||b<0||a+b>d)return 0;
  return c[a*(d+1)+b];
 }
 Z& at(int a,int b){return c[a*(d+1)+b];}
};
Row mul(const Row&f,int z){
 int n=std::abs(z),e=z>0?1:-1; Row g(f.d+n);
 for(int a=0;a<=f.d;++a)for(int b=0;b+a<=f.d;++b){
  Z v=f.get(a,b);if(v==0)continue;
  for(int t=std::abs(a-n);t<=a+n;t+=2)g.at(t,b)+=v;
  for(int t=std::abs(b-n);t<=b+n;t+=2)g.at(a,t)+=e*v;
 }
 return g;
}
Row row(const W&w){Row f;f.at(0,0)=1;for(int z:w)f=mul(f,z);return f;}
W join(W a,const W&b){a.insert(a.end(),b.begin(),b.end());return a;}
Z dot(const Row&a,const Row&b){
 Z z=0;int d=std::min(a.d,b.d);
 for(int i=0;i<=d;++i)for(int j=0;j+i<=d;++j)z+=a.get(i,j)*b.get(i,j);
 return z;
}
Z phi(const W&w){return row(w).get(0,0);}
void stamp(const std::string&s){
 std::time_t t=std::time(nullptr);std::tm u;gmtime_r(&t,&u);
 std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
int key(const std::array<int,4>&a){return a[0]+11*(a[1]+11*(a[2]+11*a[3]));}
std::vector<Z> inv(14641);
int binom[11][11];
void init_inv(){
 binom[0][0]=1;
 for(int n=1;n<=10;++n){binom[n][0]=binom[n][n]=1;
  for(int k=1;k<n;++k)binom[n][k]=binom[n-1][k-1]+binom[n-1][k];}
 for(int a=0;a<=10;++a)for(int b=0;b+a<=10;++b)
 for(int c=0;c+b+a<=10;++c)for(int d=0;d+c+b+a<=10;++d){
  std::array<int,4> v{a,b,c,d};std::vector<Z> f(1,1);
  for(int n=1;n<=4;++n)for(int k=0;k<v[n-1];++k){
   std::vector<Z> g(f.size()+n);
   for(int i=0;i<(int)f.size();++i)for(int j=std::abs(i-n);j<=i+n;j+=2)g[j]+=f[i];
   f.swap(g);
  }
  inv[key(v)]=f[0];
 }
}
using Counts=std::array<int,8>;
Z subset_rec(const Counts&v,const std::array<int,4>&tot,
 int j,std::array<int,4>&s,int sign,long long ways){
 if(j==8){std::array<int,4> t;for(int n=0;n<4;++n)t[n]=tot[n]-s[n];
  return sign*ways*inv[key(s)]*inv[key(t)];}
 Z out=0;int n=j/2;
 for(int k=0;k<=v[j];++k){
  s[n]+=k;out+=subset_rec(v,tot,j+1,s,((j%2)&&(k%2))?-sign:sign,ways*binom[v[j]][k]);s[n]-=k;
 }
 return out;
}
Z subset_phi(const Counts&v){
 std::array<int,4> t{},s{};for(int n=0;n<4;++n)t[n]=v[2*n]+v[2*n+1];
 return subset_rec(v,t,0,s,1,1);
}
struct Job{Counts v;bool pf;int r;};
std::vector<Job> jobs;
void gen(int j,int left,Counts&v){
 if(j==8){
  int r=0,minus=0;bool pf=true;
  for(int n=0;n<4;++n){r+=(v[2*n]%2)+(v[2*n+1]%2);minus+=v[2*n+1];
   if(v[2*n]&&v[2*n+1])pf=false;}
  if((r==2||r==4)&&minus%2==0)jobs.push_back({v,pf,r});
  return;
 }
 for(int k=0;k<=left;++k){v[j]=k;gen(j+1,left-k,v);}
}
struct Stat{
 long long cases=0,even=0,fiber=0,level=0,prefix=0,split=0,both=0;
 void add(const Stat&s){cases+=s.cases;even+=s.even;fiber+=s.fiber;level+=s.level;
  prefix+=s.prefix;split+=s.split;both+=s.both;}
};
Stat check(const Job&job){
 W h,r,full;int wt=0;
 for(int j=0;j<8;++j){
  int z=(j%2?-1:1)*(j/2+1),n=job.v[j];wt+=n*std::abs(z);
  for(int k=0;k<n;++k)full.push_back(z);
  for(int k=0;k<n/2;++k)h.push_back(z);
  if(n%2)r.push_back(z);
 }
 Row a=row(h),b=row(join(h,r));Z value=dot(a,b);
 assert(value==subset_phi(job.v));assert(value>=0);
 if(wt%2)assert(value==0);
 std::vector<Z> levels(a.d+1);bool fiber=false,lev=false,pre=false,split=false;
 for(int u=0;u<=a.d;++u)for(int v=0;v+u<=a.d;++v){
  Z z=a.get(u,v)*b.get(u,v);if(z<0)fiber=true;levels[u+v]+=z;
 }
 Z acc=0;
 for(const Z&z:levels){lev=lev||(z<0);acc+=z;pre=pre||(acc<0);}
 assert(acc==value);
 if(r.size()==2){
  Row c=mul(a,r[0]),d=mul(a,r[1]);assert(dot(c,d)==value);
  for(int u=0;u<=c.d;++u)for(int v=0;u+v<=c.d;++v)
   if(c.get(u,v)*d.get(u,v)<0)split=true;
 }
 return {1,wt%2==0,fiber,lev,pre,split,fiber&&split};
}
void witnesses(){
 W h{1,1,3},r{-2,-4};Row a=row(h),b=row(join(h,r));
 assert(dot(a,b)==136);assert(a.get(2,3)*b.get(2,3)==-3);
 Row c=mul(a,-2),d=mul(a,-4);assert(c.get(2,3)*d.get(2,3)==-2);
 W h2{1,-2};c=row(join(h2,{-2}));d=row(join(h2,{-4}));
 assert(dot(c,d)==20);assert(c.get(1,2)*d.get(1,2)==-1);
 assert(phi({1,3,1,3,-2,-4})==20);
 assert(row({1,3,1,3}).get(2,4)==5);
 assert(phi({1,3})==0 && phi({2,1,3})==2 && phi({2,2,1,3})==4);
 assert(4*phi({1,3})-4*phi({2,1,3})+phi({2,2,1,3})==-4);
 assert(phi({-1,-1,2})==2);
 h={-10,-9,-8,-6,3,12};r={-4,-6};a=row(h);b=row(join(h,r));
 Z top=0;for(int u=0;u<=48;++u)top+=a.get(u,48-u)*b.get(u,48-u);
 assert(top==-3512 && dot(a,b)==Z("2122501142"));
 
 h=W(6,1);r={-2,-4};a=row(h);b=row(join(h,r));
 std::vector<Z> levels(7);
 for(int u=0;u<=6;++u)for(int v=0;v+u<=6;++v)levels[u+v]+=a.get(u,v)*b.get(u,v);
 assert(levels[0]==12600&&levels[2]==97720&&levels[4]==49840&&levels[6]==-4004);
 h=W(8,1);a=row(h);b=row(join(h,r));Z z6=0,z8=0;
 for(int u=0;u<=8;++u){z6+=a.get(u,6-u)*b.get(u,6-u);z8+=a.get(u,8-u)*b.get(u,8-u);}
 assert(z6==-118620&&z8==-344760&&dot(a,b)==17131920);
 W nf{-6,-8,-10,-6,-8,-10,-2,-4};Z best=-1000000;
 assert(phi(nf)==8226);
 for(int i=0;i<8;++i)for(int j=i+1;j<8;++j){
  W rest;for(int k=0;k<8;++k)if(k!=i&&k!=j)rest.push_back(nf[k]);
  Z q=-row(rest).get(-nf[i],-nf[j]);assert(q<0);if(q>best)best=q;
 }
 assert(best==-14);
 stamp("simple level witness: [12600,97720,49840,-4004]; adjacent pair=-463380; no-flip Phi=8226, max D=-14");
 stamp("witnesses: fixed fibers -3/-2; split fiber -1; fusion D=-5; PSD value -4; level 48=-3512, Phi=2122501142");
 for(int k=2;k<=12;++k){
  W child{-2,-2,2};for(int j=2;j<=k;++j){child.push_back(1<<j);child.push_back(1<<j);}
  for(int j=1;j<=k;++j){
   int n=1<<j;auto p=std::find(child.begin(),child.end(),n);assert(p!=child.end());child.erase(p);
   p=std::find(child.begin(),child.end(),-n);assert(p!=child.end());child.erase(p);child.push_back(-2*n);
  }
  W want{-2,-(1<<(k+1))};for(int j=2;j<=k;++j)want.push_back(1<<j);
  std::sort(want.begin(),want.end());std::sort(child.begin(),child.end());assert(want==child);
  assert((int)child.size()==k+1);
 }
 stamp("positive pair-reduction branch: initial |R|=2, final |R|=k+1, k=2..12 PASS");
}

void stress(){
 const int N=6000;std::vector<std::string> fails(N);int badlevel=0,negtotal=0;
 #pragma omp parallel for schedule(dynamic,1) reduction(+:badlevel,negtotal)
 for(int it=0;it<N;++it){
  std::uint64_t state=0x93a37549ULL*(it+1);
  auto rnd=[&](){state^=state<<13;state^=state>>7;state^=state<<17;return state;};
  int maxn=4+rnd()%21,len=1+rnd()%11;std::vector<int> es(maxn+1);
  for(int n=1;n<=maxn;++n)es[n]=rnd()%2?1:-1;
  int a=1+rnd()%maxn,b=1+rnd()%maxn;while(a==b||(a-b)%2)b=1+rnd()%maxn;
  es[b]=es[a];W h,r{a*es[a],b*es[b]};
  for(int k=0;k<len;++k){int n=1+rnd()%maxn;h.push_back(n*es[n]);}
  Row x=row(h),y=row(join(h,r));std::vector<Z> lev(x.d+1);
  for(int u=0;u<=x.d;++u)for(int v=0;v+u<=x.d;++v)lev[u+v]+=x.get(u,v)*y.get(u,v);
  Z acc=0;int f=-1;bool bl=false;
  for(int t=0;t<=x.d;++t){acc+=lev[t];if(acc<0&&f<0)f=t;if(lev[t]<0)bl=true;}
  badlevel+=bl;negtotal+=acc<0;
  if(f>=0){std::ostringstream o;o<<"index="<<it<<" H=";for(int n:h)o<<n<<",";
   o<<" R="<<r[0]<<","<<r[1]<<" first_bad_prefix="<<f<<" levels=";
   for(int t=0;t<=x.d;++t)if(lev[t]!=0)o<<t<<":"<<lev[t]<<",";
   o<<" Phi="<<acc;fails[it]=o.str();}
 }
 int nfail=0;for(auto&s:fails)if(!s.empty()){if(nfail==0)std::cout<<s<<"\n";++nfail;}
 assert(nfail==0&&negtotal==0);
 std::cout<<"stress="<<N<<" negative_level="<<badlevel<<" negative_prefix="<<nfail<<" negative_total="<<negtotal<<"\n";
}

int main(int argc,char**argv){
 int threads=24;
 for(int j=1;j<argc;++j){std::string s=argv[j];
  if(s=="-h"||s=="--help"){std::cout<<"FM-STR7b exact verifier --threads N\n";return 0;}
  if(s=="--threads"&&j+1<argc)threads=std::stoi(argv[++j]);else{std::cerr<<"unknown option\n";return 2;}}
 assert(threads>0);omp_set_num_threads(threads);
 stamp("begin; cpp_int arithmetic; exhaustive signed multisets, labels 1..4, length <=10");
 init_inv();Counts v{};gen(0,10,v);
 std::array<std::array<Stat,2>,2> stats{};
 #pragma omp parallel
 {
  std::array<std::array<Stat,2>,2> local{};
  #pragma omp for schedule(dynamic,16)
  for(int j=0;j<(int)jobs.size();++j){
   Stat s=check(jobs[j]);int k=jobs[j].r==2?0:1;
   local[0][k].add(s);if(jobs[j].pf)local[1][k].add(s);
  }
  #pragma omp critical
  for(int p=0;p<2;++p)for(int k=0;k<2;++k)stats[p][k].add(local[p][k]);
 }
 stamp("census complete");
 for(int p=0;p<2;++p)for(int k=0;k<2;++k){
  Stat s=stats[p][k];
  std::cout<<(p?"pair-free":"all signed")<<" R="<<(k?4:2)<<" cases="<<s.cases
   <<" even_weight="<<s.even<<" negative_fiber="<<s.fiber<<" negative_level="<<s.level
   <<" negative_prefix="<<s.prefix<<" negative_split_fiber="<<s.split
   <<" both_cuts_fail="<<s.both<<"\n";
 }
 assert(stats[1][0].cases==1860&&stats[1][0].even==620&&stats[1][1].cases==280);
 assert(stats[0][0].cases==5940&&stats[0][1].cases==6270);
 assert(stats[0][0].prefix==0&&stats[0][1].prefix==0);
 assert(stats[1][0].fiber==84&&stats[1][0].split==253&&stats[1][0].both==78);
 assert(stats[1][1].fiber==16);
 witnesses();stamp("begin seeded stress");stress();stamp("FM-STR7b PASS");
}
'''
obj=os.memfd_create("fmstr7b_obj",0)
exe=os.memfd_create("fmstr7b_exe",0)
print("compile memory-only C++ verifier",flush=True)
subprocess.run(["g++","-Werror=return-type","-O1","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}",*sys.argv[1:]],pass_fds=(exe,),check=True)
