import os,shlex,subprocess,sys
if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
 print("FM-STR7c exact verifier: --threads N (default 24).")
 sys.exit(0)
import sympy as S
u,v,s=S.symbols("u v s")
lo=u*(v+1)*(s+u+1)-(u+1)*v*(s+v+1)
up=(u+2)*(v+1)*(s+v+1)-(u+1)*(v+2)*(s+u+1)
assert S.expand(lo-(u-v)*(s+(u+1)*(v+1)))==0
assert S.expand((s+u+v+2)*lo+s*up-(u+1)*(u-v)*(v+1)*(u+v+2))==0
print("symbolic boundary identities PASS",flush=True)
os.environ["TMPDIR"]="/dev/shm"
CPP=r'''
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;
using W=std::vector<int>;
struct Row{
 int d;std::vector<Z> c;
 Row(int n=0):d(n),c((n+1)*(n+1)){}
 Z get(int u,int v)const{
  if(u<0||v<0||u+v>d)return 0;
  return c[u*(d+1)+v];
 }
 Z& at(int u,int v){return c[u*(d+1)+v];}
};
Row axis(const Row&f,int n,int a){
 Row g(f.d+n);
 for(int u=0;u<=f.d;++u)for(int v=0;u+v<=f.d;++v){
  Z c=f.get(u,v);if(c==0)continue;int z=a?v:u;
  for(int t=std::abs(z-n);t<=z+n;t+=2)
   if(a)g.at(u,t)+=c;else g.at(t,v)+=c;
 }
 return g;
}
void add(Row&f,const Row&g,int e=1){
 assert(f.d>=g.d);
 for(int u=0;u<=g.d;++u)for(int v=0;u+v<=g.d;++v)f.at(u,v)+=e*g.get(u,v);
}
Row mul(const Row&f,int z){
 Row a=axis(f,std::abs(z),0),b=axis(f,std::abs(z),1);add(a,b,z>0?1:-1);return a;
}
Row row(const W&w){Row q;q.at(0,0)=1;for(int z:w)q=mul(q,z);return q;}
Row cut(const Row&q,int T){
 Row p(q.d);for(int u=0;u<=q.d;++u)for(int v=0;u+v<=q.d;++v)
  if(u+v<=T)p.at(u,v)=q.get(u,v);
 return p;
}
Row kernel(const Row&f,int n){
 Row out(f.d+n-1);
 for(int j=0;j<n;++j)add(out,axis(axis(f,j,0),n-1-j,1));
 return out;
}
Z dot(const Row&a,const Row&b,int T){
 Z z=0;int d=std::min({a.d,b.d,T});
 for(int u=0;u<=d;++u)for(int v=0;u+v<=d;++v)z+=a.get(u,v)*b.get(u,v);
 return z;
}
bool equal(const Row&a,const Row&b){
 int d=std::max(a.d,b.d);
 for(int u=0;u<=d;++u)for(int v=0;u+v<=d;++v)if(a.get(u,v)!=b.get(u,v))return false;
 return true;
}
Z choose(int n,int k){
 if(k<0||k>n)return 0;k=std::min(k,n-k);Z z=1;
 for(int j=1;j<=k;++j){z*=n-k+j;z/=j;}return z;
}
Z ballot(int h,int u,int v){
 if(u<0||v<0||u+v>h||(h-u-v)%2)return 0;
 int s=(h-u-v)/2;
 Z num=(u+1)*(v+1)*choose(h+2,s)*choose(h+2,s+u+1);
 int den=(h+1)*(h+2);assert(num%den==0);return num/den;
}
Z freewalk(int h,int u,int v){
 if((h+u+v)%2)return 0;
 return choose(h,(h+u+v)/2)*choose(h,(h+u-v)/2);
}
void stamp(const std::string&s){
 std::time_t t=std::time(nullptr);std::tm u;gmtime_r(&t,&u);
 std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
Row finite(const W&w,int k){
 Row f(2*k);f.at(0,0)=1;
 for(int z:w){int n=std::abs(z),e=z>0?1:-1;assert(n<=k);Row g(2*k);
  for(int u=0;u<=k;++u)for(int v=0;v<=k;++v){
   Z c=f.get(u,v);if(c==0)continue;
   for(int t=std::abs(u-n);t<=std::min(u+n,2*k-u-n);t+=2)g.at(t,v)+=c;
   for(int t=std::abs(v-n);t<=std::min(v+n,2*k-v-n);t+=2)g.at(u,t)+=e*c;
  }f=std::move(g);
 }return f;
}
W append(W h,int a,int b){h.push_back(a);h.push_back(b);return h;}
void ballot_checks(){
 long long count=0;
 for(int h=0;h<=32;++h){
  Row q=row(W(h,1)),dq=mul(q,-1);
  for(int u=0;u<=h+1;++u)for(int v=0;u+v<=h+1;++v){
   Z z=ballot(h,u,v);assert(q.get(u,v)==z);
   assert(z==freewalk(h,u,v)-freewalk(h,u+2,v)-freewalk(h,u,v+2)+freewalk(h,u+2,v+2));
   assert(dq.get(u,v)*(h+1)*(h+3)==(u-v)*(u+v+2)*ballot(h+1,u,v));
  }
  for(int T=0;T<=h+1;++T){
   Row c=mul(cut(q,T),-1);
   for(int u=0;u<=h+1;++u)for(int v=0;v<u;++v){
    assert(c.get(u,v)>=0);assert(c.get(v,u)==-c.get(u,v));++count;
   }
  }
 }
 stamp("ballot formula, reflection, and boundary checks="+std::to_string(count));
}
void fundamental(){
 struct Job{int h,a,b;};std::vector<Job> jobs;
 for(int h=0;h<=12;++h)for(int a=1;a<=10;++a)for(int b=a+2;b<=10;b+=2)jobs.push_back({h,a,b});
 long long checks=0;
 #pragma omp parallel for schedule(dynamic,1) reduction(+:checks)
 for(int it=0;it<(int)jobs.size();++it){
  auto j=jobs[it];Row q=row(W(j.h,1));
  Row B=kernel(mul(q,-j.a),j.b),g=mul(mul(q,-j.a),-j.b);
  assert(equal(mul(B,-1),g));
  for(int u=0;u<=B.d;++u)for(int v=0;v<u;++v){
   assert(B.get(u,v)>=0);assert(B.get(v,u)==-B.get(u,v));
  }
  for(int T=0;T<=j.h+1;++T){
   Row c=mul(cut(q,T),-1);Z p=dot(q,g,T);
   assert(p>=0&&p==dot(c,B,std::min(c.d,B.d)));++checks;
  }
  for(int e:{-1,1}){
   Row qm=row(W(j.h,-1)),gm=mul(mul(qm,e*j.a),e*j.b);
   int ep=e*(j.a%2?-1:1);Row gp=mul(mul(q,ep*j.a),ep*j.b);
   for(int T=0;T<=j.h+1;++T){assert(dot(qm,gm,T)==dot(q,gp,T));assert(dot(qm,gm,T)>=0);}
  }
 }
 stamp("fundamental positive-character certificate prefixes="+std::to_string(checks));
}
void finite_checks(){
 struct C{int k;W h;int a,b;long long p,g;};
 std::vector<C> cs{
  {3,W(2,1),1,3,28,26},
  {3,W(4,1),1,3,2640,2970},
  {4,{-2},-2,-4,4,2},
  {5,{1,1,1,1,-4},-2,-4,13600,13666}
 };
 for(auto&c:cs){
  W hr=append(c.h,c.a,c.b);Row x=row(c.h),y=row(hr);
  Row xf=finite(c.h,c.k),yf=finite(hr,c.k);
  assert(dot(x,y,c.k)==c.p);assert(dot(xf,yf,2*c.k)==c.g);
 }
 for(int T=1;T<=8;++T){
  Row e(T);e.at(T,0)=1;
  Row AB=cut(axis(cut(axis(e,1,1),T),1,0),T);
  Row BA=cut(axis(cut(axis(e,1,0),T),1,1),T);
  for(int u=0;u<=T;++u)for(int v=0;u+v<=T;++v)
   assert(AB.get(u,v)-BA.get(u,v)==((u==T-1&&v==1)?-1:0));
 }
 int checks=0;
 for(int h=0;h<=5;++h)for(int n=2;n<=6;++n)for(int s:{-1,1}){
  W H(h,1);H.push_back(s*n);int a=2,b=4;W HR=append(H,-a,-b);
  Row x=row(H),y=row(HR);int d=h+n,r=a+b;
  for(int T=0;T<=d;++T){
   int k=std::max({n,a,b,T,(d+r+T+1)/2});
   assert(dot(x,y,T)==dot(finite(H,k),finite(HR,k),T));++checks;
  }
 }
 stamp("finite-level defects PASS; stable projected comparisons="+std::to_string(checks));
}
void obstruction_checks(){
 Row f=row({2}),c=mul(f,-1),B=kernel(mul(f,-3),5);
 assert(c.get(1,0)==1&&c.get(2,1)==-1&&c.get(3,0)==1);
 assert(B.get(1,0)==1&&B.get(2,1)==1&&B.get(3,0)==2);
 assert(dot(c,B,c.d)==4);
 Row q=row(W(9,1)),dq=mul(q,-1);
 B=mul(mul(kernel(mul(q,-2),3),-4),-5);
 assert(dq.get(6,4)==42&&B.get(6,4)==-16722);
 Row g=mul(B,-1);assert(dot(q,g,q.d)==Z("2233171270"));
 stamp("active Sp4 sign defects: one general label -1; four remainders -16722 at coefficient 42");
}
void one_general(){
 struct Job{int h,n,s,a,b,e;};std::vector<Job> jobs;
 for(int h=0;h<=8;++h)for(int n=2;n<=12;++n)for(int s:{-1,1})
 for(int a=1;a<=18;++a)for(int b=a+2;b<=18;b+=2)for(int e:{-1,1}){
  if(h&&a==1&&e==-1)continue;if((a==n||b==n)&&e!=s)continue;
  jobs.push_back({h,n,s,a,b,e});
 }
 assert(jobs.size()==25560);stamp("one-general-label test cases=25560");
 long long prefixes=0;
 #pragma omp parallel for schedule(dynamic,8) reduction(+:prefixes)
 for(int it=0;it<(int)jobs.size();++it){
  auto j=jobs[it];W H(j.h,1);H.push_back(j.s*j.n);
  Row x=row(H),y=mul(mul(x,j.e*j.a),j.e*j.b);Z p=0;
  for(int t=0;t<=x.d;++t){
   for(int u=0;u<=t;++u)p+=x.get(u,t-u)*y.get(u,t-u);
   assert(p>=0);++prefixes;
  }
 }
 stamp("one-general-label prefixes checked="+std::to_string(prefixes)+"; no counterexample");
}
int main(int argc,char**argv){
 int threads=24;
 for(int j=1;j<argc;++j){std::string a=argv[j];
  if(a=="-h"||a=="--help"){std::cout<<"FM-STR7c verifier --threads N\n";return 0;}
  if(a=="--threads"&&j+1<argc)threads=std::stoi(argv[++j]);
  else{std::cerr<<"unknown option\n";return 2;}
 }
 assert(threads>0);omp_set_num_threads(threads);
 stamp("FM-STR7c begin; arbitrary-precision coefficients");
 ballot_checks();fundamental();finite_checks();obstruction_checks();one_general();
 stamp("FM-STR7c PASS");
}
'''
obj=os.memfd_create("fmstr7c_obj",0)
exe=os.memfd_create("fmstr7c_exe",0)
subprocess.run(["g++","-Werror=return-type","-O1","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}",*sys.argv[1:]],pass_fds=(exe,),check=True)
