import os,shlex,subprocess,sys
if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
 print("FM-STR7 exact memory-only verifier: --threads N (default 24).")
 sys.exit(0)
os.environ["TMPDIR"]="/dev/shm"
CPP = r'''
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>
using namespace std;
using U=uint64_t;
constexpr int Q=33, QQ=Q*Q;
constexpr U LIMIT=100000000ULL;
static_assert(9586981ULL*LIMIT<numeric_limits<U>::max());
void stamp(const string&s){
 time_t t=time(nullptr);tm z;gmtime_r(&t,&z);char b[32];
 strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&z);
 #pragma omp critical(output)
 {cout<<b<<" "<<s<<endl;}
}
struct DP{array<uint32_t,2*QQ> x{};};
DP step(const DP&r,int w,int n,bool minus,bool odd){
 DP v;
 for(int a=0;a<=w;a++)for(int b=(w-a)%2;b<=w-a;b+=2){
  int z=a*Q+b;auto p=r.x[z],q=r.x[QQ+z];if(!(p|q))continue;
  for(int c=abs(a-n);c<=a+n;c+=2){int u=c*Q+b;v.x[u]+=p;v.x[QQ+u]+=q;}
  for(int c=abs(b-n);c<=b+n;c+=2){int u=a*Q+c;v.x[u]+=minus?q:p;v.x[QQ+u]+=minus?p:q;}
 }
 if(odd)for(int a=0;a< Q;a++)v.x[a*Q+a]=v.x[QQ+a*Q+a]=0;
 return v;
}
struct Stats{
 array<U,9> count{},bad{},negative{};
 void merge(const Stats&s){for(int i=0;i<9;i++){count[i]+=s.count[i];bad[i]+=s.bad[i];negative[i]+=s.negative[i];}}
};
void descend(const DP&r,int d,int w,bool parity,int stop,Stats&s){
 if(!parity){s.count[d]++;if(r.x[QQ])s.bad[d]++;s.negative[d]+=r.x[QQ];}
 if(d==stop)return;
 for(int n=1;n<=4;n++)for(int m=0;m<2;m++){
  auto v=step(r,w,n,m,parity^m);
  descend(v,d+1,w+n,parity^m,stop,s);
 }
}
void tails(int threads){
 DP root;root.x[0]=1;Stats ans;descend(root,0,0,false,2,ans);
 atomic<int> done{0};
 #pragma omp parallel num_threads(threads)
 {
  Stats local;
  #pragma omp for schedule(dynamic,1)
  for(int z=0;z<512;z++){
   DP r=root;int w=0;bool par=false;
   for(int j=2;j>=0;j--){int u=(z>>(3*j))&7;int n=u/2+1,m=u%2;par^=m;r=step(r,w,n,m,par);w+=n;}
   descend(r,3,w,par,8,local);
   int k=++done;
   if(k%64==0)stamp("tail prefixes "+to_string(k)+"/512");
  }
  #pragma omp critical(stats)
  ans.merge(local);
 }
 U tot=0,fail=0;
 for(int l=0;l<=8;l++){
  cout<<"TAIL L="<<l<<" profiles="<<ans.count[l]<<" fail="<<ans.bad[l]<<" negative_fixed="<<ans.negative[l]<<endl;
  tot+=ans.count[l];fail+=ans.bad[l];
 }
 assert(tot==9586981 && fail==4675955);
 DP r=root;int w=0;bool par=false;for(int t:vector<int>{0,0,1,1}){par^=t%2;r=step(r,w,t/2+1,t%2,par);w+=t/2+1;}
 assert(r.x[0]==6&&r.x[QQ]==4);
 cout<<"TAIL total="<<tot<<" fail="<<fail<<endl;
}

using Row=map<pair<int,int>,int64_t>;
using Word=vector<int>;
int label(int t){return t/2+1;}
int sg(int t){return t%2?-1:1;}
int parity(const Word&w){int p=0;for(int t:w)p^=t%2;return p;}
Row srow(const Word&w){
 Row r{{{0,0},1}};
 U dim=1;
 for(int t:w){
  int n=label(t);dim*=2*(n+1);assert(dim<=LIMIT);
  Row v;
  for(auto [ab,z]:r){
   auto [a,b]=ab;
   for(int c=abs(a-n);c<=a+n;c+=2)v[{c,b}]+=z;
   for(int c=abs(b-n);c<=b+n;c+=2)v[{a,c}]+=sg(t)*z;
  }
  r.clear();for(auto [ab,z]:v)if(z)r[ab]=z;
 }
 return r;
}
vector<int64_t> ordinary(const Word&w){
 vector<int64_t> r{1};
 for(int t:w){int n=label(t);vector<int64_t> v(r.size()+n);
  for(int a=0;a<(int)r.size();a++)if(r[a])
   for(int c=abs(a-n);c<=a+n;c+=2)v[c]+=r[a];
  r=move(v);
 }return r;
}
int64_t subset_phi(const Word&w){
 int size=1<<w.size(),all=size-1;vector<int64_t> m(size);vector<int> signs(size,1);
 for(int s=0;s<size;s++){Word z;for(int i=0;i<(int)w.size();i++)if(s>>i&1){z.push_back(w[i]);signs[s]*=sg(w[i]);}m[s]=ordinary(z)[0];}
 int64_t ans=0;for(int s=0;s<size;s++)ans+=signs[s]*m[s]*m[all^s];return ans;
}
array<vector<Word>,5> words(){
 array<vector<Word>,5> out;Word w;
 auto gen=[&](auto&&self,int lo,int left)->void{
  if(!left){out[w.size()].push_back(w);return;}
  for(int t=lo;t<8;t++){w.push_back(t);self(self,t,left-1);w.pop_back();}
 };
 for(int k=0;k<=4;k++)gen(gen,0,k);return out;
}
U fact(int k){U x=1;for(int j=2;j<=k;j++)x*=j;return x;}
U orders(const Word&w){U x=fact(w.size());array<int,8> m{};for(int t:w)m[t]++;for(int n:m)x/=fact(n);return x;}
struct Path{int sign=1,mask=0;array<int,5>a{},b{};};
using Fibers=map<pair<int,int>,vector<Path>>;
Fibers paths(const Word&w){
 Fibers f;Path p;
 auto gen=[&](auto&&self,int j)->void{
  if(j==(int)w.size()){f[{p.a[j],p.b[j]}].push_back(p);return;}
  int a=p.a[j],b=p.b[j],n=label(w[j]),old=p.sign,mask=p.mask;
  for(int c=abs(a-n);c<=a+n;c+=2){p.a[j+1]=c;p.b[j+1]=b;self(self,j+1);}
  p.sign*=sg(w[j]);p.mask|=1<<j;
  for(int c=abs(b-n);c<=b+n;c+=2){p.a[j+1]=a;p.b[j+1]=c;self(self,j+1);}
  p.sign=old;p.mask=mask;
 };gen(gen,0);return f;
}
bool cg(int a,int n,int c){return abs(a-n)<=c&&c<=a+n&&((a+n-c)%2==0);}
void valid(const Path&p,const Word&w){
 for(int j=0;j<(int)w.size();j++){
  if(p.mask>>j&1){assert(p.a[j]==p.a[j+1]);assert(cg(p.b[j],label(w[j]),p.b[j+1]));}
  else{assert(p.b[j]==p.b[j+1]);assert(cg(p.a[j],label(w[j]),p.a[j+1]));}
 }
}
U maxflow(const vector<U>&cap,int minusmask){
 int n=cap.size(),s=n,t=n+1,N=n+2;vector<vector<U>>a(N,vector<U>(N));
 for(int i=0;i<n;i++)if(cap[i]&&__builtin_parity(unsigned(i&minusmask))){
  a[s][i]=cap[i];
  for(int j=0;j<n;j++)if(cap[j]&&!__builtin_parity(unsigned(j&minusmask))&&__builtin_popcount(unsigned(i^j))<=2)a[i][j]=LIMIT;
 }else a[i][t]=cap[i];
 U ans=0;
 for(;;){
  vector<int>p(N,-1);queue<int>q;q.push(s);p[s]=s;
  while(!q.empty()&&p[t]<0){int u=q.front();q.pop();for(int v=0;v<N;v++)if(a[u][v]&&p[v]<0){p[v]=u;q.push(v);}}
  if(p[t]<0)return ans;
  U z=LIMIT;for(int v=t;v!=s;v=p[v])z=min(z,a[p[v]][v]);
  for(int v=t;v!=s;v=p[v]){a[p[v]][v]-=z;a[v][p[v]]+=z;}ans+=z;
 }
}
void folds(const array<vector<Word>,5>& ws){
 U words_done=0,atoms=0,fixed=0,ordered=0,local_checks=0,local_bad=0;
 Word first;pair<int,int>first_end;array<U,3> first_flow{};
 for(int h=0;h<=4;h++)for(const Word&w:ws[h]){
  auto fs=paths(w);auto r=srow(w);int64_t phi=0;U total=0;
  Word full=w;full.insert(full.end(),w.rbegin(),w.rend());
  ordered+=orders(full);words_done++;
  int mm=0;for(int j=0;j<h;j++)if(w[j]%2)mm|=1<<j;
  for(const auto& [end,v]:fs){
   int n=v.size();vector<int>plus,minus,partner(n);iota(partner.begin(),partner.end(),0);
   vector<U>cap(1<<h);
   for(int i=0;i<n;i++){valid(v[i],w);cap[v[i].mask]++;(v[i].sign==1?plus:minus).push_back(i);}
   for(int j=0;j<(int)min(plus.size(),minus.size());j++){partner[plus[j]]=minus[j];partner[minus[j]]=plus[j];}
   auto invol=[&](int i,int j){if(partner[i]!=i)return pair{partner[i],j};if(partner[j]!=j)return pair{i,partner[j]};return pair{i,j};};
   int64_t residual=int64_t(plus.size())-int64_t(minus.size());assert(residual==r[end]);
   U here=0;
   for(int i=0;i<n;i++)for(int j=0;j<n;j++){
    auto [a,b]=invol(i,j);assert(a>=0&&a<n&&b>=0&&b<n);
    assert(invol(a,b)==make_pair(i,j));
    if(a==i&&b==j){assert(v[i].sign*v[j].sign==1);here++;}
    else assert(v[i].sign*v[j].sign==-v[a].sign*v[b].sign);
   }
   assert(here==U(residual*residual));fixed+=here;atoms+=U(n)*n;total+=U(n)*n;phi+=residual*residual;
   if(!plus.empty()&&!minus.empty()){
    U f=maxflow(cap,mm);local_checks++;
    if(f<min(plus.size(),minus.size())){
     local_bad++;if(first.empty()){first=w;first_end=end;first_flow={plus.size(),minus.size(),f};}
    }
   }
  }
  assert(phi==subset_phi(full));assert((phi==srow(full)[{0,0}]));assert(total<=LIMIT);
 }
 assert(words_done==495&&local_checks==5160&&local_bad==980);
 assert((first==Word({0,0,3})&&first_end==make_pair(0,0)&&first_flow==array<U,3>{1,1,0}));
 cout<<"FOLD classes="<<words_done<<" ordered_words="<<ordered<<" atoms="<<atoms<<" positive_fixed="<<fixed<<endl;
 cout<<"LOCAL mixed_fibers="<<local_checks<<" deficient="<<local_bad<<" first=(+1,+1,-2), endpoint=(0,0), P=N=1, flow=0"<<endl;
 U proper=0;
 for(int k=3;k<=10;k++){
  vector<int> ns;for(int j=0;j<k-1;j++)ns.push_back(1<<j);ns.push_back((1<<(k-1))-1);
  for(int s=1;s<(1<<k)-1;s++){int sum=0,mx=0;for(int j=0;j<k;j++)if(s>>j&1){sum+=ns[j];mx=max(mx,ns[j]);}assert(2*mx>sum);proper++;}
  assert(accumulate(ns.begin(),ns.end()-1,0)==ns.back());
 }
 cout<<"ARITY proper_subsets="<<proper<<" PASS"<<endl;
}
void cuts(const array<vector<Word>,5>&ws){
 map<Word,Row> rows;for(const auto&v:ws)for(const auto&w:v)rows[w]=srow(w);
 U all=0,fail=0,weighted=0,adaptive_fail=0;Word first_bad;
 array<U,9> expected{0,0,0,0,0,0,23,7,183};
 for(int l=0;l<=8;l++){
  int h=l/2,k=l-h;U cases=0,bad=0;map<Word,bool>coherent;
  for(const Word&a:ws[h])for(const Word&b:ws[k])if(parity(a)==parity(b)){
   bool good=true;int64_t phi=0;
   for(auto [v,c]:rows[a]){auto it=rows[b].find(v);if(it!=rows[b].end()){int64_t z=c*it->second;phi+=z;if(z<0)good=false;}}
   assert(phi>=0);
   cases++;bad+=!good;weighted+=orders(a)*orders(b);
   Word w=a;w.insert(w.end(),b.begin(),b.end());sort(w.begin(),w.end());coherent[w]=coherent[w]||good;
  }
  U no=0;for(auto [w,good]:coherent)if(!good){no++;if(first_bad.empty())first_bad=w;}
  assert(no==expected[l]);
  cout<<"CUT L="<<l<<" half_classes="<<cases<<" bad="<<bad<<" no_balanced_cut="<<no<<endl;
  all+=cases;fail+=bad;adaptive_fail+=no;
 }
 assert(weighted==9586981&&adaptive_fail==213);
 assert(first_bad==Word({0,0,2,5,5,6}));
 cout<<"CUT total="<<all<<" bad="<<fail<<" ordered_words="<<weighted<<endl;
 Word small{0,0,2,5,5,6};int sc=0;
 for(int s=0;s<64;s++)if(__builtin_popcount(unsigned(s))==3){
  Word a,b;for(int j=0;j<6;j++)(s>>j&1?a:b).push_back(small[j]);
  bool bad=false;auto A=srow(a),B=srow(b);for(auto[v,c]:A)if(c*B[v]<0)bad=true;assert(bad);sc++;
 }
 assert(sc==20);cout<<"SMALL first_no_balanced_cut Phi="<<subset_phi(small)<<endl;
 Word w{1,2,4,6,9,10,12,14};
 assert((subset_phi(w)==956&&srow(w)[{0,0}]==956));
 array<U,5> best{};best.fill(LIMIT);array<int,5> count{};
 for(int s=1;s<256;s++){
  int k=__builtin_popcount(unsigned(s));if(k<2||k>4)continue;
  Word a,b;for(int j=0;j<8;j++)(s>>j&1?a:b).push_back(w[j]);
  auto A=srow(a),B=srow(b);int64_t phi=0;U neg=0;
  for(auto[v,c]:A){int64_t z=c*B[v];phi+=z;if(z<0)neg+=U(-z);}
  assert(phi==956&&neg>0);
  if(k==2){int u=label(a[0]),v=label(a[1]);int64_t D=sg(a[1])*B[{u,v}];assert(D<0&&neg==U(-2*D));for(int c=abs(u-v);c<=u+v;c+=2)assert((B[{c,0}]>=0));}
  best[k]=min(best[k],neg);count[k]++;
 }
 assert(count[2]==28&&count[3]==56&&count[4]==70);
 assert(best[2]==4&&best[3]==42&&best[4]==38);
 cout<<"RESIDUAL Phi=956 cut_counts=28,56,70 minimum_negative_mass=4,42,38 PASS"<<endl;
 Word a{1,1,1},b{1,3,3};auto A=srow(a),B=srow(b);
 assert((A[{2,1}]==-3&&B[{2,1}]==1));
 int64_t phi=0;for(auto[v,c]:A)phi+=c*B[v];assert(phi==28);
}

void whole_local(int threads){
 vector<Word> all;Word w;
 auto gen=[&](auto&&self,int lo,int left)->void{
  if(!left){if(!parity(w))all.push_back(w);return;}
  for(int t=lo;t<8;t++){w.push_back(t);self(self,t,left-1);w.pop_back();}
 };
 for(int k=0;k<=8;k++)gen(gen,0,k);
 vector<U>fail(all.size()),deficit(all.size());U checked=0,weighted=0;atomic<int>done{0};
 #pragma omp parallel for num_threads(threads) schedule(dynamic,4) reduction(+:checked,weighted)
 for(int z=0;z<(int)all.size();z++){
  const Word&v=all[z];int L=v.size(),N=1<<L,mask=N-1,mm=0;
  vector<U>m(N),cap(N);for(int i=0;i<L;i++)if(v[i]%2)mm|=1<<i;
  for(int s=0;s<N;s++){Word t;for(int i=0;i<L;i++)if(s>>i&1)t.push_back(v[i]);m[s]=ordinary(t)[0];}
  U p=0,n=0;
  for(int s=0;s<N;s++){cap[s]=m[s]*m[mask^s];(__builtin_parity(unsigned(s&mm))?n:p)+=cap[s];}
  assert(p>=n);
  if(n){checked++;U f=maxflow(cap,mm);if(f<n){fail[z]=1;deficit[z]=n-f;}}
  weighted+=orders(v);
  int d=++done;if(d%1024==0)stamp("whole-word local profiles "+to_string(d)+"/"+to_string(all.size()));
 }
 U bad=accumulate(fail.begin(),fail.end(),U(0)),pfbad=0;int first=-1,pffirst=-1;
 for(int i=0;i<(int)all.size();i++)if(fail[i]){
  if(first<0)first=i;bool pf=true;array<int,4> signs{};
  for(int t:all[i])signs[t/2]|=1<<(t%2);
  for(int b:signs)if(b==3)pf=false;
  if(pf){pfbad++;if(pffirst<0)pffirst=i;}
 }
 assert(weighted==9586981);
 cout<<"WHOLE_LOCAL classes="<<all.size()<<" negative_classes="<<checked<<" deficient="<<bad<<" ordered_words="<<weighted<<endl;
 assert(bad==276&&first>=0&&all[first]==Word({0,0,1,1,2,2})&&deficit[first]==4);
 if(first>=0){cout<<"WHOLE_LOCAL first_types=";for(int t:all[first])cout<<t<<",";cout<<" deficit="<<deficit[first]<<endl;}
 cout<<"WHOLE_LOCAL pair_free_deficient="<<pfbad<<" first_types=";
 if(pffirst>=0)for(int t:all[pffirst])cout<<t<<",";
 cout<<" deficit="<<(pffirst>=0?deficit[pffirst]:0)<<endl;
 for(int scale:vector<int>{1,3,5}){
 Word q{2*(scale-1),2*(scale-1),2*(2*scale-1),2*(2*scale-1),2*(3*scale-1)+1,2*(3*scale-1)+1};vector<U> cap(64);int mm=48;U src=0,dst=0;vector<int> left,nb(64);
 for(int m=0;m<64;m++){Word a,b;for(int i=0;i<6;i++)(m>>i&1?a:b).push_back(q[i]);cap[m]=ordinary(a)[0]*ordinary(b)[0];
  if(__builtin_popcount(unsigned(m&3))==1&&__builtin_popcount(unsigned(m&12))==1&&__builtin_popcount(unsigned(m&48))==1){left.push_back(m);src+=cap[m];}}
 for(int m:left)for(int n=0;n<64;n++)if(cap[n]&&!__builtin_parity(unsigned(n&mm))&&__builtin_popcount(unsigned(m^n))<=2)nb[n]=1;
 for(int n=0;n<64;n++)if(nb[n])dst+=cap[n];assert(src==8&&dst==4);
 cout<<"HALL scale="<<scale<<" demand=8 capacity=4 Phi="<<subset_phi(q)<<endl;
 }
 U tests=0;
 for(int r=3;r<=7;r++){
  int q0=(1<<(r-1))-1,T=2*q0+1;vector<int> ns;
  for(int j=0;j<r-1;j++)ns.push_back(1<<j);ns.push_back(q0);
  for(int j=0;j<r;j++)ns.push_back(T*ns[j]);
  int H=(1<<r)-1,ALL=(1<<(2*r))-1;
  for(int m=0;m<=ALL;m++){
   int sum[2]{},mx[2]{};for(int j=0;j<2*r;j++){int b=(m>>j)&1;sum[b]+=ns[j];mx[b]=max(mx[b],ns[j]);}
   bool possible=true;for(int b=0;b<2;b++)if(sum[b]%2||2*mx[b]>sum[b])possible=false;
   assert(possible==(m==0||m==ALL||m==H||m==(ALL^H)));tests++;
  }
 }
 cout<<"SEPARATED support_tests="<<tests<<" PASS"<<endl;
}

int main(int argc,char**argv){
 int threads=24;
 for(int i=1;i<argc;i++){
  string a=argv[i];
  if(a=="-h"||a=="--help"){cout<<"FM-STR7 exact memory-only verifier: --threads N (default 24).\n";return 0;}
  if(a=="--threads"&&i+1<argc)threads=stoi(argv[++i]);else{cerr<<"bad argument\n";return 2;}
 }
 if(threads<1)return 2;
 auto ws=words();stamp("folded involutions");
 folds(ws);stamp("cut obstructions");cuts(ws);stamp("all ordered tail tests");tails(threads);
 stamp("whole-word bounded-factor surgeries");whole_local(threads);
 cout<<"FM-STR7 PASS"<<endl;
}
'''
obj=os.memfd_create("fmstr7_obj",0)
exe=os.memfd_create("fmstr7_exe",0)
print("compile C++ source",flush=True)
subprocess.run(["g++","-Werror=return-type","-O1","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
print("link C++ verifier",flush=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin","-fopenmp",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}",*sys.argv[1:]],pass_fds=(exe,),check=True)
