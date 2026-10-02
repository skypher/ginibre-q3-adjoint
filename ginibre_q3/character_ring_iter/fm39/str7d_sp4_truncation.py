import sys
if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
 print("FM-STR7d exact verifier: --threads N (default 24).");sys.exit(0)
CPP=r'''
#include <algorithm>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;
using Key=std::pair<int,int>;
using Row=std::map<Key,Z>;
using W=std::vector<int>;
Row mul(const Row&f,int z){
 int n=std::abs(z),e=z>0?1:-1;Row g;
 for(auto &[uv,c]:f){auto [u,v]=uv;
  for(int t=std::abs(u-n);t<=u+n;t+=2)g[{t,v}]+=c;
  for(int t=std::abs(v-n);t<=v+n;t+=2)g[{u,t}]+=e*c;
 }
 for(auto it=g.begin();it!=g.end();)if(it->second==0)it=g.erase(it);else++it;
 return g;
}
Row row(const W&w){Row f{{{0,0},1}};for(int z:w)f=mul(f,z);return f;}
Z get(const Row&f,int u,int v){auto it=f.find({u,v});return it==f.end()?Z(0):it->second;}
bool cg(int x,int y,int z){return z>=std::abs(x-y)&&z<=x+y&&(x+y+z)%2==0;}
int overlap(int a,int b,int s,int u){
 if((a+b-s-u)%2)return 0;
 int lo=std::max(std::abs(a-b),std::abs(s-u)),hi=std::min(a+b,s+u);
 return hi<lo?0:(hi-lo)/2+1;
}
std::vector<Z> levels(const Row&f,int a,int b,int e){
 int D=0;for(auto &[uv,c]:f)D=std::max(D,uv.first+uv.second);
 std::vector<Z>L(D+1);
 for(auto &[uv,c]:f){auto[u,v]=uv;Z g=0;
  for(auto &[st,d]:f){auto[s,t]=st;
   int K=(v==t?overlap(a,b,s,u):0)+(u==s?overlap(a,b,t,v):0);
   K+=e*((cg(s,a,u)&&cg(t,b,v))+(cg(s,b,u)&&cg(t,a,v)));
   if(K)g+=d*K;
  }L[u+v]+=c*g;
 }return L;
}
void stamp(const std::string&s){
 std::time_t t=std::time(nullptr);std::tm u;gmtime_r(&t,&u);
 std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
void show(const W&H,int a,int b,int e,const std::vector<Z>&L){
 std::cout<<"H";for(int z:H)std::cout<<" "<<z;
 std::cout<<" R "<<e*a<<" "<<e*b<<" levels";
 for(int t=0;t<(int)L.size();++t)if(L[t]!=0)std::cout<<" "<<t<<":"<<L[t];
 std::cout<<std::endl;
}
struct Dense{
 int d;std::vector<Z> c;Dense(int n=0):d(n),c((n+1)*(n+1)){}
 Z get(int u,int v)const{
  if(u<0||v<0||u+v>d)return 0;return c[u*(d+1)+v];
 }
 Z& at(int u,int v){return c[u*(d+1)+v];}
};
Dense daxis(const Dense&f,int n,int axis){
 Dense g(f.d+n);
 for(int v=0;v<=f.d;++v){
  int last=f.d-v;std::vector<Z> ps(last+3);
  for(int u=0;u<=last;++u)ps[u+2]=ps[u]+(axis?f.get(v,u):f.get(u,v));
  for(int t=0;t+v<=g.d;++t){
   int lo=std::abs(t-n),hi=std::min(last,t+n);
   if((hi-lo)%2) --hi;if(lo>hi)continue;
   Z z=ps[hi+2]-ps[lo];
   if(axis)g.at(v,t)=z;else g.at(t,v)=z;
  }
 }
 return g;
}
Dense dmul(const Dense&f,int n,int e){
 Dense g=daxis(f,n,0),b=daxis(f,n,1);
 for(size_t t=0;t<g.c.size();++t)g.c[t]+=e*b.c[t];return g;
}
Dense fund(int h){
 std::vector<Z> C(h+3);C[0]=1;
 for(int j=1;j<=h+2;++j)C[j]=C[j-1]*(h+3-j)/j;
 Dense f(h);
 for(int u=0;u<=h;++u)for(int v=0;u+v<=h;++v)if((h-u-v)%2==0){
  int j=(h-u-v)/2;Z z=Z(u+1)*(v+1)*C[j]*C[j+u+1];
  f.at(u,v)=z/(Z(h+1)*(h+2));
 }
 return f;
}
std::vector<Z> dlevels(const Dense&f,int a,int b,int e){
 Dense g=dmul(dmul(f,a,e),b,e);std::vector<Z>L(f.d+1);
 for(int u=0;u<=f.d;++u)for(int v=0;u+v<=f.d;++v)L[u+v]+=f.get(u,v)*g.get(u,v);
 return L;
}
Row character(int A,int B){
 Row f;for(int j=0;j<=B;++j)for(int k=0;k<=A-B;++k)f[{j+k,j+A-B-k}]+=1;return f;
}
Row cut(const Row&f,int T){
 Row g;for(auto&[uv,c]:f)if(uv.first+uv.second<=T)g[uv]=c;return g;
}
void add(Row&f,const Row&g,int k=1){
 for(auto&[uv,c]:g)f[uv]+=k*c;
 for(auto it=f.begin();it!=f.end();)if(it->second==0)it=f.erase(it);else++it;
}
Z prefix(const Row&f,const Row&g,int T){
 Z p=0;for(auto&[uv,c]:f)if(uv.first+uv.second<=T)p+=c*get(g,uv.first,uv.second);return p;
}
Row kernel(const Row&f,int n){
 Row g;
 for(auto&[uv,c]:f)for(int j=0;j<n;++j)
  for(int u=std::abs(uv.first-j);u<=uv.first+j;u+=2)
  for(int v=std::abs(uv.second-(n-1-j));v<=uv.second+n-1-j;v+=2)g[{u,v}]+=c;
 return g;
}
void identities(){
 long long checks=0;
 for(int A=0;A<=12;++A)for(int B=0;B<=A;++B){
  Row f=character(A,B),E{{{A+1,B},1},{{B,A+1},-1}};assert(mul(f,-1)==E);
  for(int T=0;T<=2*A+1;++T){
   Row g;if(T>=A-B){int J=std::min(B,(T-A+B)/2);g=character(A-B+J,J);}
   assert(cut(f,T)==g);++checks;
  }
 }
 Row packet=character(3,0);add(packet,character(1,0));
 assert(row({1,2})==packet);
 Row F=mul(character(3,2),-3),G=mul(mul(F,-2),-8);
 assert(prefix(F,G,4)==-2);
 Row q=row(W(5,1)),dq=mul(q,-1);
 assert(get(dq,4,2)==5);
 F=mul(q,-3);G=mul(mul(F,-2),-8);
 assert(prefix(F,G,4)==32180&&prefix(F,G,8)==51936);
 Row V=character(5,0);add(V,character(4,3),3);
 F=mul(V,5);G=mul(mul(F,-2),-4);
 assert(get(F,0,0)==2&&get(G,0,0)==-2&&prefix(F,G,0)==-4);
 Row two=row({1,1});
 assert(get(two,0,0)==2&&get(mul(two,-1),1,0)==1);
 Row C=row({1,3,1,3});
 Z sum=0;for(int c=2;c<=6;c+=2)sum+=get(mul(C,c),0,0);
 assert(get(C,2,4)==5&&sum==30&&get(mul(mul(C,-2),-4),0,0)==20);
 for(int h=0;h<=12;++h)for(int n=h+2;n<=h+10;++n){
  Row A=mul(mul(row(W(h,1)),n),-1);
  assert(get(A,n,h+1)==-1);
 }
 for(int t=0;t<=24;++t){
  Dense q0=fund(t),B=dmul(dmul(dmul(q0,1,-1),1,-1),1,-1);
  Dense M=dmul(fund(t+2),1,-1);Z den=Z(t+1)*(t+2)*(t+3)*(t+4);
  for(int A=0;A<=t+2;++A)for(int b=0;b<=A&&A+b<=t+2;++b){
   Z x=Z(A+b+3)*(A+b+3),y=Z(A-b+1)*(A-b+1);
   Z polynomial=x*y-(3*t+13)*(x+y)+Z(15)*t*t+120*t+241;
   assert(B.get(A+1,b)*den==M.get(A+1,b)*polynomial);
  }
  int A=(t+3)/2,b=(t+2)/2;
  assert(A+b==t+2&&B.get(A+1,b)<0);
 }
 stamp("branching/truncation identities="+std::to_string(checks)+"; exact defects and quartic PASS");
}
void packets(){
 long long checks=0;
 for(int r=0;r<=4;++r)for(int h=r;h<=8;++h){
  W H(h,1);H.insert(H.end(),r,2);Row f=row(H);
  for(int a=1;a<=8;++a)for(int b=a+2;b<=10;b+=2){
   Row B=kernel(mul(f,-a),b),g=mul(mul(f,-a),-b);
   for(auto&[uv,c]:B)if(uv.first>uv.second)assert(c>=0);
   int D=h+2*r;
   for(int T=0;T<=D;++T){
    Row A=mul(cut(f,T),-1);
    for(auto&[uv,c]:A)if(uv.first>uv.second)assert(c>=0);
    assert(prefix(f,g,T)==prefix(A,B,2*D+a+b)&&prefix(f,g,T)>=0);++checks;
   }
  }
 }
 stamp("positive packet certificate prefixes="+std::to_string(checks));
}
void sparse_search(){
 long long one=0,two=0,pref=0;
 #pragma omp parallel for schedule(dynamic,1) reduction(+:one,two,pref)
 for(int it=0;it<5000;++it){
  std::mt19937 rng(769131+it);int k=1+(it%2),h=rng()%11;
  W H(h,1);int n=2+rng()%100;H.push_back((rng()%2?1:-1)*n);
  if(k==2){int m=2+rng()%100,e=rng()%2?1:-1;if(m==n)e=H.back()>0?1:-1;H.push_back(e*m);}
  int a=2+rng()%180,b=2+rng()%180;if((a+b)%2)++b;if(a==b)b+=2;if(a>b)std::swap(a,b);
  int e=rng()%2?1:-1;bool ok=true;
  for(int z:H)if((std::abs(z)==a||std::abs(z)==b)&&(z>0?1:-1)!=e)ok=false;
  if(!ok)continue;
  Row f=row(H);auto L=levels(f,a,b,e);Z p=0;
  for(auto&z:L){p+=z;assert(p>=0);++pref;}
  if(k==1)++one;else++two;
  if(it%1000==0){
   #pragma omp critical
   stamp("mixed-scale index "+std::to_string(it));
  }
 }
 stamp("mixed-scale one="+std::to_string(one)+" two="+std::to_string(two)+" prefixes="+std::to_string(pref));
}
void dense_search(){
 long long cases=0,pref=0;
 #pragma omp parallel for schedule(dynamic,1) reduction(+:cases,pref)
 for(int it=0;it<1000;++it){
  std::mt19937 rng(3413427+it);int h=10+rng()%191,n=2+rng()%81;
  int eta=rng()%2?1:-1,a=2+rng()%110,b=2+rng()%110,e=rng()%2?1:-1;
  if((a+b)%2)++b;if(a==b)b+=2;if(a>b)std::swap(a,b);
  if((a==n||b==n)&&e!=eta)continue;
  Dense f=dmul(fund(h),n,eta);auto L=dlevels(f,a,b,e);Z p=0;
  for(auto&z:L){p+=z;assert(p>=0);++pref;}++cases;
  if(it%200==0){
   #pragma omp critical
   stamp("Gaussian-scale index "+std::to_string(it));
  }
 }
 stamp("Gaussian-scale cases="+std::to_string(cases)+" prefixes="+std::to_string(pref));
}
int main(int argc,char**argv){
 int threads=24;
 for(int j=1;j<argc;++j){
  std::string a=argv[j];
  if(a=="-h"||a=="--help"){std::cout<<"FM-STR7d verifier --threads N (default 24)\n";return 0;}
  if(a=="--threads"&&j+1<argc)threads=std::stoi(argv[++j]);else return 2;
 }
 assert(threads>0);omp_set_num_threads(threads);
 stamp("FM-STR7d begin; cpp_int coefficients");
 identities();packets();sparse_search();dense_search();stamp("FM-STR7d PASS");
}
'''
import os,shlex,subprocess,sys
os.environ["TMPDIR"]="/dev/shm"
obj=os.memfd_create("fmstr7d_obj",0)
exe=os.memfd_create("fmstr7d_exe",0)
subprocess.run(
 ["g++","-Werror=return-type","-O2","-std=c++17","-fopenmp","-pipe",
  "-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],
 input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
 ["g++","-###","-fno-use-linker-plugin","-fopenmp",
  f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
 capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}",*sys.argv[1:]],
               pass_fds=(exe,),check=True)
