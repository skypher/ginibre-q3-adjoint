import argparse
import itertools
import random
import re
import shlex
import subprocess
import os
from datetime import datetime, timezone

ap = argparse.ArgumentParser(description="Independent exact HPP checker; C++ coefficients are arbitrary-precision integers.")
ap.add_argument("--threads", type=int, default=24)
ap.add_argument("--census", default="/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/fx6_52_noflip.log")
args = ap.parse_args()
assert args.threads > 0

def stamp(msg):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), msg, flush=True)

pat = re.compile(r"^NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(-?\d+) D=")
census = []
for line in open(args.census):
    m = pat.match(line)
    if not m or int(m.group(1)) > 52:
        continue
    B = tuple(map(int, m.group(3).split()))
    word = tuple(sorted(B + (int(m.group(2)),)))
    assert sum(abs(z) for z in B) == int(m.group(1))
    census.append((word, int(m.group(4))))
assert len(census) == 75532
rng = random.Random(114)
picked = sorted(rng.sample(range(len(census)), 300))
census_sample = [census[i] for i in picked]
stamp(f"loaded W<=52 rows={len(census)}; seed=114 sample=300")

def box_words(maxlabel, maxlen, pairfree, exactlen=None):
    alphabet = tuple(sorted([s*n for n in range(1, maxlabel+1) for s in (1,-1)]))
    lo = exactlen if exactlen is not None else 2
    hi = exactlen if exactlen is not None else maxlen
    out = []
    for n in range(lo, hi+1):
        for w in itertools.combinations_with_replacement(alphabet, n):
            if sum(z < 0 for z in w) % 2:
                continue
            if sum(map(abs, w)) % 2:
                continue
            if pairfree and any(-z in w for z in w):
                continue
            out.append(w)
    return out

g1 = box_words(4, 8, True)
g2 = box_words(3, 8, False)
g3 = box_words(5, 8, False, 8)
assert (len(g1),len(g2),len(g3)) == (1070,775,6158)
jobs = []
def add(group, word, expected=-1):
    jobs.append((group, tuple(word), expected))
for w in g1: add(1,w)
for w in g2: add(2,w)
for w in g3: add(3,w)
for w,phi in census_sample: add(4,w,phi)
f1 = (1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8)
stamp(f"generated words: pairfree-4={len(g1)}; pair-allowed-3={len(g2)}; label5-length8={len(g3)}")

wire = "".join(str(g)+" "+str(phi)+" "+str(len(w))+" "+" ".join(map(str,w))+"\n" for g,w,phi in jobs)

CPP = r"""
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <random>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>
using Z=boost::multiprecision::cpp_int;
using Key=std::pair<int,int>;
using Tab=std::map<Key,Z>;
using Word=std::vector<int>;
static void stamp(const std::string&s){
 std::time_t t=std::time(nullptr);std::tm u;gmtime_r(&t,&u);
 std::cout<<std::put_time(&u,"%FT%TZ")<<" "<<s<<std::endl;
}
static int mass(const Word&w){int x=0;for(int z:w)x+=std::abs(z);return x;}
static std::vector<int> cg(int a,int b){
 std::vector<int> v;for(int c=std::abs(a-b);c<=a+b;c+=2)v.push_back(c);return v;
}
static Tab multiply_factor(const Tab&f,int signed_n){
 int n=std::abs(signed_n),eps=signed_n>0?1:-1;Tab g;
 for(const auto&[rs,c]:f){
  int r=rs.first,s=rs.second;
  for(int t=std::abs(r-n);t<=r+n;t+=2)g[{t,s}]+=c;
  for(int t=std::abs(s-n);t<=s+n;t+=2)g[{r,t}]+=eps*c;
 }
 for(auto it=g.begin();it!=g.end();)if(it->second==0)it=g.erase(it);else ++it;
 return g;
}
static thread_local std::map<Word,Tab> cache;
static void clear_cache(){cache.clear();}
static const Tab& table(Word w){
 std::sort(w.begin(),w.end());
 auto it=cache.find(w);if(it!=cache.end())return it->second;
 Word order=w;
 std::stable_sort(order.begin(),order.end(),[](int a,int b){
  if(std::abs(a)!=std::abs(b))return std::abs(a)>std::abs(b);
  return a<b;
 });
 Tab f{{{0,0},Z(1)}};
 for(int z:order)f=multiply_factor(f,z);
 return cache.emplace(std::move(w),std::move(f)).first->second;
}
static Word half(const Word&w,unsigned mask){
 Word a;for(int i=0;i<(int)w.size();++i)if(mask>>i&1u)a.push_back(w[i]);return a;
}
static std::pair<Word,Word> canon(Word a,Word b){
 if(b<a)std::swap(a,b);return {std::move(a),std::move(b)};
}
static std::map<int,Z> product_table(const Tab&a,const Tab&b,int lx,int ly){
 std::map<int,Z> out;
 for(const auto&[rs,x]:a){
  auto it=b.find(rs);if(it!=b.end())out[lx*rs.first+ly*rs.second]+=x*it->second;
 }
 for(auto it=out.begin();it!=out.end();)if(it->second==0)it=out.erase(it);else ++it;
 return out;
}
struct Ratio{
 Z n=0,d=1;int T=-1;Word A,B;bool set=false;
};
static Z gcdz(Z a,Z b){if(a<0)a=-a;if(b<0)b=-b;while(b!=0){Z r=a%b;a=b;b=r;}return a;}
static std::string exact_ratio(const Ratio&r){
 if(!r.set)return "undefined";
 Z g=gcdz(r.n,r.d);return (r.n/g).convert_to<std::string>()+"/"+(r.d/g).convert_to<std::string>();
}
struct Metrics{
 long long raw=0,distinct=0,fail=0,fliptests=0,flipfail=0;
 Z minpi=0;bool minpi_set=false;int minpiT=-1;Word minpiA,minpiB;
 Ratio ratio;
 std::map<std::pair<int,int>,long long> dir_fail;
 std::map<std::pair<int,int>,std::tuple<Word,Word,int,Z>> dir_first;
};
static void update_ratio(Metrics&m,const Word&A,const Word&B,int T,const Z&p,const Z&pos){
 if(pos<=0)return;
 if(!m.ratio.set||p*m.ratio.d<m.ratio.n*pos){
  m.ratio.n=p;m.ratio.d=pos;m.ratio.T=T;m.ratio.A=A;m.ratio.B=B;m.ratio.set=true;
 }
}
static void check_split(Metrics&m,const Word&A,const Word&B,bool directions){
 const Tab&fa=table(A);const Tab&fb=table(B);
 auto layer=product_table(fa,fb,1,1);
 int hmax=std::min(mass(A),mass(B));
 Z pref=0,pos=0;
 for(int t=0;t<=hmax;++t){
  auto it=layer.find(t);Z x=(it==layer.end()?Z(0):it->second);
  pref+=x;if(x>0)pos+=x;
  if(!m.minpi_set||pref<m.minpi){
   m.minpi=pref;m.minpiT=t;m.minpiA=A;m.minpiB=B;m.minpi_set=true;
  }
  update_ratio(m,A,B,t,pref,pos);
 }
 if(pref<0||std::any_of(layer.begin(),layer.end(),[](const auto&x){return x.second<0;})){
  // Prefix failure is tested below, at every nonzero height.
 }
 pref=0;
 for(const auto&[t,x]:layer){pref+=x;if(pref<0){++m.fail;break;}}
 if(directions){
  static const std::vector<std::pair<int,int>> dirs={{1,0},{1,-1},{1,2},{1,3},{2,3},{1,10}};
  for(auto d:dirs){
   auto l=product_table(fa,fb,d.first,d.second);Z v=0;bool bad=false;
   for(const auto&[t,x]:l){v+=x;if(v<0){bad=true;++m.dir_fail[d];m.dir_first.emplace(d,std::make_tuple(A,B,t,v));break;}}
  }
 }
 ++m.distinct;
}
static Metrics scan_word(Word w,bool directions){
 std::sort(w.begin(),w.end());int n=w.size();
 Metrics m;
 std::set<std::pair<Word,Word>> seen;
 const unsigned lim=1u<<(n-1);
 for(unsigned mask=1;mask<lim;++mask){
  Word A=half(w,mask),B;
  for(int i=0;i<n;++i)if(!(mask>>i&1u))B.push_back(w[i]);
  ++m.raw;auto ab=canon(A,B);
  if(!seen.insert(ab).second)continue;
  check_split(m,ab.first,ab.second,directions);
 }
 return m;
}
static void merge(Metrics&dst,const Metrics&src){
 dst.raw+=src.raw;dst.distinct+=src.distinct;dst.fail+=src.fail;dst.fliptests+=src.fliptests;dst.flipfail+=src.flipfail;
 if(src.minpi_set&&(!dst.minpi_set||src.minpi<dst.minpi)){
  dst.minpi=src.minpi;dst.minpiT=src.minpiT;dst.minpiA=src.minpiA;dst.minpiB=src.minpiB;dst.minpi_set=true;
 }
 if(src.ratio.set&&(!dst.ratio.set||src.ratio.n*dst.ratio.d<dst.ratio.n*src.ratio.d))dst.ratio=src.ratio;
 for(auto&[d,c]:src.dir_fail)dst.dir_fail[d]+=c;
 for(auto&[d,w]:src.dir_first)dst.dir_first.emplace(d,w);
}
static std::string wstr(const Word&w){std::string s="(";for(int z:w)s+=std::to_string(z)+",";return s+")";}
static void report(const std::string&name,const Metrics&m){
 std::cout<<"RESULT "<<name<<" raw_mask_splits="<<m.raw<<" distinct_multiset_splits="<<m.distinct
 <<" HPP_failures="<<m.fail<<" flip_tests="<<m.fliptests<<" flip_fail="<<m.flipfail<<"\n";
 std::cout<<"  min_prefix="<<m.minpi<<" at T="<<m.minpiT<<" A="<<wstr(m.minpiA)<<" B="<<wstr(m.minpiB)<<"\n";
 std::cout<<"  min_prefix/positive_mass="<<exact_ratio(m.ratio)<<" at T="<<m.ratio.T
 <<" positive_mass="<<m.ratio.d<<" prefix="<<m.ratio.n<<" A="<<wstr(m.ratio.A)<<" B="<<wstr(m.ratio.B)<<"\n";
}
static std::vector<std::vector<Z>> U_monomials(int n){
 std::vector<std::vector<Z>> u(n+1);u[0]={Z(1)};
 if(n>=1)u[1]={Z(0),Z(1)};
 for(int k=1;k<n;++k){
  u[k+1].assign(k+2,Z(0));
  for(int j=0;j<(int)u[k].size();++j)u[k+1][j+1]+=u[k][j];
  for(int j=0;j<(int)u[k-1].size();++j)u[k+1][j]-=u[k-1][j];
 }
 return u;
}
static std::map<Key,Z> poly_multiply(const std::map<Key,Z>&a,const std::map<Key,Z>&b){
 std::map<Key,Z>c;
 for(const auto&[ij,x]:a)for(const auto&[kl,y]:b)c[{ij.first+kl.first,ij.second+kl.second}]+=x*y;
 for(auto it=c.begin();it!=c.end();)if(it->second==0)it=c.erase(it);else ++it;
 return c;
}
static Tab direct_polynomial(const Word&w){
 int degree=mass(w);auto u=U_monomials(degree);
 std::map<Key,Z>p{{{0,0},Z(1)}};
 for(int z:w){
  int n=std::abs(z),eps=z>0?1:-1;std::map<Key,Z>factor;
  for(int k=0;k<(int)u[n].size();++k)if(u[n][k]!=0){
   factor[{k,0}]+=u[n][k];factor[{0,k}]+=eps*u[n][k];
  }
  p=poly_multiply(p,factor);
 }
 std::vector<std::vector<Z>> mono(degree+1);mono[0]={Z(1)};
 for(int i=0;i<degree;++i){
  mono[i+1].assign(i+2,Z(0));
  for(int j=0;j<(int)mono[i].size();++j){
   mono[i+1][j+1]+=mono[i][j];if(j>0)mono[i+1][j-1]+=mono[i][j];
  }
 }
 Tab out;
 for(const auto&[ij,c]:p)for(int r=0;r<(int)mono[ij.first].size();++r)
  for(int s=0;s<(int)mono[ij.second].size();++s)
   out[{r,s}]+=c*mono[ij.first][r]*mono[ij.second][s];
 for(auto it=out.begin();it!=out.end();)if(it->second==0)it=out.erase(it);else ++it;
 return out;
}
static void direct_validation(){
 std::mt19937 gen(1142026);long long checks=0;
 for(int k=0;k<120;++k){
  int n=1+gen()%6;Word w;for(int j=0;j<n;++j){int a=1+gen()%5;w.push_back((gen()&1)?a:-a);}
  const Tab&a=table(w);Tab b=direct_polynomial(w);assert(a==b);++checks;
 }
 stamp("direct x,y polynomial-to-U validation exact PASS words="+std::to_string(checks));
}
static Z pref_at(const std::map<int,Z>&l,int T){Z x=0;for(const auto&[t,v]:l)if(t<=T)x+=v;return x;}
static void f1_and_pair_counterexample(){
 Word f1={1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8};
 clear_cache();Metrics m=scan_word(f1,true);
 assert(m.distinct==599&&m.fail==0);report("F1 all distinct splits",m);
 const Tab&full=table(f1);assert(full.at({0,0})==8150742);
 std::cout<<"F1 distinct_splits=599 Phi=8150742 height_failures=0\n";
 for(auto d:std::vector<std::pair<int,int>>{{1,0},{1,-1},{1,2},{1,3},{2,3},{1,10}}){
  std::cout<<"F1 direction=("<<d.first<<","<<d.second<<") failures="<<m.dir_fail[d];
  auto it=m.dir_first.find(d);
  if(it!=m.dir_first.end()){
   auto[A,B,T,P]=it->second;
   std::cout<<" first A="<<wstr(A)<<" B="<<wstr(B)<<" cutoff="<<T<<" prefix="<<P;
  }
  std::cout<<"\n";
 }
 assert(m.dir_fail[std::make_pair(1,0)]==37&&m.dir_fail[std::make_pair(1,-1)]==13);
 assert(m.dir_fail[std::make_pair(1,2)]==5&&m.dir_fail[std::make_pair(1,3)]==15);
 assert(m.dir_fail[std::make_pair(2,3)]==5&&m.dir_fail[std::make_pair(1,10)]==35);
 Word w10a={3,3,3,5,5,5,5},w10b={-4,-4,-4,-2,1,1,3,8};
 Word w1ma={3,3,3,5,5,5,5,8},w1mb={-4,-4,-4,-2,1,1,3};
 assert(pref_at(product_table(table(w10a),table(w10b),1,0),5)==-607879);
 assert(pref_at(product_table(table(w1ma),table(w1mb),1,-1),1)==-1333203);
 std::cout<<"F1 published witnesses exact: (1,0) T=5 -607879; (1,-1) T=1 -1333203\n";
 Word A={-1,-1},B={1,1,2};
 const Tab&fa=table(A);const Tab&fb=table(B);
 assert(fa.at({0,0})==2&&fa.at({2,0})==1&&fa.at({1,1})==-2&&fa.at({0,2})==1);
 Tab expectedB{{{0,0},2},{{2,0},3},{{0,2},3},{{4,0},1},{{2,2},2},{{1,1},4},{{3,1},2},{{1,3},2},{{0,4},1}};
 assert(fb==expectedB);
 std::map<Key,Z>P;for(auto&[k,x]:fa){auto it=fb.find(k);if(it!=fb.end())P[k]=x*it->second;}
 assert(P.size()==4&&P.at(Key(0,0))==4&&P.at(Key(2,0))==3&&P.at(Key(1,1))==-8&&P.at(Key(0,2))==3);
 auto L10=product_table(fa,fb,1,0),L1m1=product_table(fa,fb,1,-1),L12=product_table(fa,fb,1,2),L11=product_table(fa,fb,1,1);
 assert(pref_at(L10,1)==-1&&pref_at(L12,3)==-1&&pref_at(L1m1,0)==-1);
 assert(pref_at(L11,0)==4&&pref_at(L11,2)==2);
 assert(fa.at({0,0})*fb.at({0,0})==4);
 std::cout<<"PAIR_COUNTEREXAMPLE exact P={(0,0):4,(2,0):3,(1,1):-8,(0,2):3}; "
          <<"height prefixes T=0:4 T=2:2; direction failures (1,0) T=1:-1, (1,2) T=3:-1, (1,-1) T=0:-1\n";
}
int main(int argc,char**argv){
 int threads=24;
 for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--threads"&&i+1<argc)threads=std::stoi(argv[++i]);else if(a=="-h"||a=="--help"){std::cout<<"FM-CHK114 exact HPP screen; --threads N\n";return 0;}else return 2;}
 assert(threads>0);omp_set_num_threads(threads);stamp("FM-CHK114 begin; cpp_int fusion tables");
 direct_validation();
 std::vector<Word> words[5];std::vector<Z> expected[5];
 int g,n;Z phi;
 while(std::cin>>g>>phi>>n){Word w(n);for(int&z:w)std::cin>>z;assert(g>=1&&g<=4);words[g].push_back(w);expected[g].push_back(phi);}
 assert(words[1].size()==1070&&words[2].size()==775&&words[3].size()==6158&&words[4].size()==300);
 Metrics total[5];std::atomic<int> done{0};
 for(int group=1;group<=4;++group){
  long long raw=0,distinct=0;std::vector<Metrics> rowmetrics(words[group].size());
  #pragma omp parallel
  {
   Metrics local;long long localrows=0;
   #pragma omp for schedule(dynamic,1)
   for(int i=0;i<(int)words[group].size();++i){
    clear_cache();Metrics m=scan_word(words[group][i],false);
    const Tab&full=table(words[group][i]);
    if(group==4){
     assert(full.at({0,0})==expected[group][i]);
     std::set<int> labels;for(int z:words[group][i]){assert(!labels.count(-z));labels.insert(z);}
     Z value=full.at({0,0});
     for(int a=0;a<(int)words[group][i].size();++a)for(int b=a+1;b<(int)words[group][i].size();++b){
      Word child=words[group][i];child[a]=-child[a];child[b]=-child[b];
      if(table(child).at({0,0})<=value)++m.flipfail;
      ++m.fliptests;
     }
    }
    rowmetrics[i]=std::move(m);
    int k=++done;long long nrow=++localrows;
    if((group!=4&&i%100==0)||(group==4&&i%20==0)){
     #pragma omp critical
     stamp("group="+std::to_string(group)+" completed_word_index="+std::to_string(i+1)+"/"+std::to_string(words[group].size()));
    }
   }
  }
  for(auto&m:rowmetrics)merge(total[group],m);
  if(group==1)assert(total[group].raw==79842);
  if(group==2)assert(total[group].raw==60553);
  if(group==3)assert(total[group].raw==782066);
  report(group==1?"labels<=4 length<=8 pairfree":group==2?"labels<=3 length<=8 pairs allowed":group==3?"labels<=5 length=8 pairs allowed":"300 seeded W<=52 census rows",total[group]);
  assert(total[group].fail==0);if(group==4)assert(total[group].fliptests>0&&total[group].flipfail==0);
 }
 f1_and_pair_counterexample();
 stamp("FM-CHK114 all exact screens PASS");
}
"""
import os, shlex, subprocess
os.environ["TMPDIR"] = "/dev/shm"
obj = os.memfd_create("fmchk114_obj", 0)
exe = os.memfd_create("fmchk114_exe", 0)
subprocess.run(["g++","-O2","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],
    input=CPP,text=True,pass_fds=(obj,),check=True)
probe=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(line) for line in probe.stderr.splitlines() if "/collect2 " in line)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}","--threads",str(args.threads)],input=wire,text=True,pass_fds=(exe,),check=True)
