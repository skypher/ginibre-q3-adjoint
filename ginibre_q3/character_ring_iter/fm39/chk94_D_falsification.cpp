#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <gmpxx.h>
using namespace std;
using Z = mpz_class;
using B = vector<int>;
using U64 = uint64_t;

struct Profile { string name; B b; int p; };
struct Removal { vector<int> r; int cost; };
struct Result { Z g, best_delta; vector<int> best_R; int tested=0, positive=0; bool have=false; };

static map<string,Z> gp_cache;
static long long gp_calls=0, cache_hits=0;

int weight(const B& b){ int w=0; for(int x:b) w+=abs(x); return w; }
int maxlabel(const B& b){ int m=0; for(int x:b) m=max(m,abs(x)); return m; }
int large_count(const B& b){ int k=0; for(int x:b) if(abs(x)>=3) ++k; return k; }
B canon(B b){
  sort(b.begin(),b.end(),[](int a,int c){
    if(abs(a)!=abs(c)) return abs(a)>abs(c);
    return a<c;
  }); return b;
}
string bkey(B b){
  b=canon(b); string s; for(int x:b){ s+=to_string(x); s.push_back(','); } return s;
}
string gpkey(const B& b,int p){ return to_string(p)+"|"+bkey(b); }
bool pairfree(const B& b){
  map<int,int> sg;
  for(int x:b){ int n=abs(x), e=x>0?1:-1; auto it=sg.find(n);
    if(it!=sg.end()&&it->second!=e) return false; sg[n]=e;
  } return true;
}
bool residual(const B& b,int p){
  int W=weight(b), M=maxlabel(b), k=large_count(b);
  if(!pairfree(b)||p<6||p<M||p>W||((W-p)&1)||k<2) return false;
  int d=(W-p)/2; return d>=8&&M<=d;
}
U64 statekey(int s,int t){ return (U64(uint32_t(s))<<32)|uint32_t(t); }
int sx(U64 k){ return int(k>>32); }
int sy(U64 k){ return int(k&0xffffffffU); }

Z gp_pruned(B b,int p){
  b=canon(b); string ky=gpkey(b,p);
  auto it=gp_cache.find(ky); if(it!=gp_cache.end()){ ++cache_hits; return it->second; }
  ++gp_calls; int W=weight(b);
  if(p<0||p>W||((W-p)&1)){ gp_cache[ky]=0; return 0; }
  unordered_map<U64,Z> cur,nx; cur.reserve(64); cur[statekey(0,0)]=1;
  int rem=W, fi=0; auto last=chrono::steady_clock::now();
  for(int f:b){
    int n=abs(f), e=f>0?1:-1; rem-=n; nx.clear();
    nx.reserve(cur.size()*2+16);
    for(const auto& kv:cur){
      int s=sx(kv.first), t=sy(kv.first); const Z& v=kv.second;
      if(t<=rem){
        for(int s2=abs(s-n);s2<=s+n;s2+=2)
          if(abs(s2-p)<=rem) nx[statekey(s2,t)]+=v;
      }
      if(abs(s-p)<=rem){
        for(int t2=abs(t-n);t2<=t+n;t2+=2)
          if(t2<=rem){ U64 q=statekey(s,t2); if(e>0) nx[q]+=v; else nx[q]-=v; }
      }
    }
    for(auto z=nx.begin();z!=nx.end();){ if(z->second==0) z=nx.erase(z); else ++z; }
    cur.swap(nx); ++fi;
    auto now=chrono::steady_clock::now();
    if(chrono::duration_cast<chrono::seconds>(now-last).count()>=30){
      cerr<<"gp progress factors="<<fi<<"/"<<b.size()<<" remaining="<<rem
          <<" states="<<cur.size()<<" p="<<p<<" W="<<W<<"\n"; cerr.flush(); last=now;
    }
  }
  auto z=cur.find(statekey(p,0)); Z out=z==cur.end()?Z(0):z->second;
  gp_cache.emplace(ky,out); return out;
}
Z gp_direct(B b,int p){
  int W=weight(b), S=W+1; auto ix=[S](int x,int y){return size_t(x)*S+y;};
  vector<Z> cur(size_t(S)*S),nx(size_t(S)*S); cur[ix(0,0)]=1; int done=0;
  for(int f:canon(b)){
    int n=abs(f),e=f>0?1:-1,nd=done+n;
    for(int x=0;x<=nd;++x)for(int y=0;x+y<=nd;++y)nx[ix(x,y)]=0;
    for(int x=0;x<=done;++x)for(int y=0;x+y<=done;++y){
      const Z&v=cur[ix(x,y)];if(v==0)continue;
      for(int q=abs(x-n);q<=x+n;q+=2)nx[ix(q,y)]+=v;
      for(int q=abs(y-n);q<=y+n;q+=2){if(e>0)nx[ix(x,q)]+=v;else nx[ix(x,q)]-=v;}
    }
    cur.swap(nx);done=nd;
  }
  return p<=W?cur[ix(p,0)]:Z(0);
}
void check_pruning(){
  vector<pair<B,int>> cs={
    {{1,-1,3,4},3},{{1,-1,3,4},5},
    {{-1,-1,-1,-1,-2,-3,-3,4},6},
    {{1,1,1,-2,3,3,-5},6},{{-1,2,2,3,-4},4}
  };
  int n=0;
  for(auto &q:cs){ Z a=gp_pruned(q.first,q.second), b=gp_direct(q.first,q.second);
    if(a!=b){ cerr<<"prune check failed B="<<bkey(q.first)<<" p="<<q.second<<"\n"; exit(3); } ++n;
  }
  cout<<"pruned/direct exact checks="<<n<<"\n";
}

vector<Removal> all_even_removals(const B& b){
  map<int,int> c; for(int x:b)++c[x];
  vector<int> cls; for(auto [x,n]:c)cls.push_back(x);
  vector<Removal> out;
  for(int x:cls) if(abs(x)%2==0) out.push_back({{x},abs(x)});
  for(size_t i=0;i<cls.size();++i)for(size_t j=i;j<cls.size();++j){
    int x=cls[i],y=cls[j];
    if((abs(x)&1)!=(abs(y)&1))continue;
    if(i==j&&c[x]<2)continue;
    vector<int> r={x,y}; sort(r.begin(),r.end());
    out.push_back({r,abs(x)+abs(y)});
  }
  sort(out.begin(),out.end(),[](const Removal&a,const Removal&b){
    if(a.cost!=b.cost)return a.cost<b.cost; return a.r<b.r;
  });
  return out;
}
B erase_R(B b,const vector<int>& r){
  for(int x:r){auto it=find(b.begin(),b.end(),x);assert(it!=b.end());b.erase(it);}return b;
}
string rstr(const vector<int>& r){
  string s="(";for(size_t i=0;i<r.size();++i){if(i)s+=",";s+=(r[i]>0?"+":"")+to_string(r[i]);}return s+")";
}
string qstr(const Z& n,const Z& d){
  if(d==0)return "undefined"; mpq_class q(n,d);q.canonicalize();return q.get_num().get_str()+"/"+q.get_den().get_str();
}
Result analyze(const Profile& pr,bool quiet=false){
  assert(residual(pr.b,pr.p));
  Result z; z.g=gp_pruned(pr.b,pr.p);
  auto rs=all_even_removals(pr.b); int total=rs.size();
  for(size_t i=0;i<rs.size();++i){
    Z h=0; int childW=weight(pr.b)-rs[i].cost;
    if(pr.p<=childW) h=gp_pruned(erase_R(pr.b,rs[i].r),pr.p);
    Z d=z.g-h; ++z.tested; if(d>=0)++z.positive;
    if(!z.have||d>z.best_delta){z.have=true;z.best_delta=d;z.best_R=rs[i].r;}
    if(!quiet&&(i+1)%20==0){cout<<"  removal progress "<<i+1<<"/"<<total<<" W="<<weight(pr.b)<<" p="<<pr.p<<"\n";cout.flush();}
  }
  return z;
}
struct GlobalBest { bool set=false; mpq_class score; string name; int p=0; Z g=0,d=0; vector<int> R; };
static GlobalBest global_best;
void report(const Profile&pr,const Result&z){
  cout<<pr.name<<" W="<<weight(pr.b)<<" factors="<<pr.b.size()<<" p="<<pr.p
      <<" delta="<<(weight(pr.b)-pr.p)/2<<" removals="<<z.tested
      <<" monotone-count="<<z.positive<<" g="<<z.g
      <<" max-delta="<<z.best_delta<<" R="<<rstr(z.best_R)
      <<" max-delta/g="<<qstr(z.best_delta,z.g)<<" B="<<bkey(pr.b)<<"\n";cout.flush();
  if(z.best_delta<0)cout<<"D COUNTEREXAMPLE profile="<<pr.name<<" p="<<pr.p<<" max_R="<<rstr(z.best_R)<<" delta="<<z.best_delta<<" B="<<bkey(pr.b)<<"\n";
  if(z.g>0){
    mpq_class q(z.best_delta,z.g);q.canonicalize();
    if(!global_best.set||q<global_best.score){
      global_best.set=true;global_best.score=q;global_best.name=pr.name;
      global_best.p=pr.p;global_best.g=z.g;global_best.d=z.best_delta;global_best.R=z.best_R;
    }
  }
}
void report_global(){
  if(global_best.set)cout<<"smallest exact max_R delta/g="<<global_best.score.get_num().get_str()<<"/"
      <<global_best.score.get_den().get_str()<<" at "<<global_best.name<<" p="<<global_best.p
      <<" g="<<global_best.g<<" max-delta="<<global_best.d<<" R="<<rstr(global_best.R)<<"\n";
  else cout<<"no positive denominator row was fully scored\n";
}
void add(B&b,int n,int e,int c=1){for(int i=0;i<c;++i)b.push_back(e*n);}
vector<Profile> seeds(){
 vector<Profile> v; B b;
 b.clear();add(b,1,-1,2);for(int n=3;n<=15;n+=2)add(b,n,1,2);
 v.push_back({"Rule1p-W128",canon(b),16});
 b.clear();add(b,1,-1);add(b,2,-1,2);for(int n=3;n<=15;n+=2)add(b,n,-1,2);
 v.push_back({"Rule1p-W131",canon(b),15});
 b.clear();add(b,1,-1,4);add(b,2,-1);add(b,3,-1,8);add(b,4,-1);
 v.push_back({"RuleM-W34",canon(b),6});
 b.clear();add(b,1,-1,10);add(b,8,-1,4);for(int n=3;n<=31;n+=2)add(b,n,-1,31/n);
 v.push_back({"RuleW-Minerva-W427",canon(b),31});
 b.clear();for(int n=1;n<=35;n+=2)if(n!=9)add(b,n,1,35/n);add(b,9,-1,4);
 v.push_back({"RuleW-Juno-W540",canon(b),36});
 b.clear();for(int n=1;n<=35;n+=2)if(n!=9)add(b,n,-1,35/n);add(b,9,1,4);
 v.push_back({"RuleW-Juno-reflection-W540",canon(b),36});
 b.clear();int labs[]={1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31};
 int mult[]={62,20,12,8,7,5,4,4,3,2,2,2,2,2,2,2};
 for(int i=0;i<16;++i)add(b,labs[i],labs[i]==9?1:-1,mult[i]);
 v.push_back({"RuleW-Ceres-W869",canon(b),57});
 b.clear();add(b,11,-1,5);add(b,54,-1);add(b,2,-1);
 for(int n=1;n<=53;n+=2)if(n!=11)add(b,n,1,54/n);
 v.push_back({"WM-Juno-W1283",canon(b),55});
 return v;
}
vector<Profile> families(){
 vector<Profile> v;B b;
 add(b,20,-1,25);add(b,25,1,8);add(b,31,-1,4);add(b,39,1,3);
 v.push_back({"W941-heavy-20-class-ladder-25-31-39",canon(b),39});
 b.clear();add(b,1,-1,100);add(b,3,1,250);add(b,5,-1,10);add(b,7,1);
 v.push_back({"W907-dominant-3-small-max-7",canon(b),7});
 b.clear();add(b,1,1,100);add(b,3,-1,80);add(b,4,1,30);add(b,7,-1,20);add(b,12,1,10);
 v.push_back({"W720-mixed-sign-five-classes",canon(b),14});
 return v;
}
vector<Profile> stress_profiles(){
 vector<Profile> v;
 auto put=[&](string name,B b){
   int W=weight(b),M=maxlabel(b),p=W-2*max(8,M);
   b=canon(b);
   if(!residual(b,p)){cerr<<"stress profile invalid "<<name<<" W="<<W<<" p="<<p<<" M="<<M<<"\n";exit(6);}
   v.push_back({name,b,p});
 };
 B b;
 for(int m:{8,16,32,48}){
   b.clear();int labs[]={5,7,9,11,13,17};
   for(int i=0;i<6;++i)add(b,labs[i],(i&1)?-1:1,m);
   add(b,1,-1,m/2);
   put("equal-multiplicity-six-odd-classes-m"+to_string(m),b);
 }
 b.clear();add(b,3,1,100);add(b,47,-1,6);add(b,53,1,2);
 put("near-tied-weight-3x100-vs-47x6",b);
 b.clear();add(b,29,-1,12);add(b,31,1,8);add(b,37,-1,5);add(b,43,1,3);add(b,53,-1,2);
 put("heavy-middle-class-ladder",b);
 b.clear();add(b,1,1,600);add(b,3,-1,30);add(b,5,1,10);add(b,7,-1,2);
 put("dominant-small-label-max7",b);
 b.clear();add(b,1,1,350);add(b,3,-1,80);add(b,4,1,40);add(b,7,-1,20);add(b,12,1,10);add(b,19,-1,5);
 put("mixed-sign-six-class",b);
 b.clear();add(b,2,-1,230);for(int n=3;n<=21;n+=2)add(b,n,-1);
 put("all-minus-odd-run-doubled-even",b);
 b.clear();add(b,1,1,1800);add(b,3,-1,25);add(b,5,1,12);add(b,7,-1,5);add(b,83,1);
 put("one-large-label-many-small",b);
 for(int N:{47,75,107}){
   b.clear();add(b,1,-1,2);for(int n=3;n<=N;n+=2)add(b,n,1);
   put("Juno-minus1-squared-plus-odd-through-"+to_string(N),b);
 }
 return v;
}
int choose_p(const B& b,int target){
 int W=weight(b),M=maxlabel(b),lo=max(6,M),hi=W-2*max(8,M);
 int best=-1,dist=INT32_MAX;
 for(int p=lo;p<=hi;++p)if(((W-p)&1)==0){
   int d=abs(p-target);if(d<dist){best=p;dist=d;}
 }
 return best;
}
B pad_to_target(B b,int target){
  int gap=target-weight(b);if(gap<=0)return canon(b);
  auto sign_for=[&](int n){for(int x:b)if(abs(x)==n)return x>0?1:-1;return 1;};
  while(gap>=49){int e=sign_for(49);add(b,49,e);gap-=49;}
  if(gap>0)add(b,1,sign_for(1),gap);
  return canon(b);
}
uint64_t rng_state=0xd1b54a32d192ed03ULL;
uint64_t rnd(){rng_state^=rng_state<<7;rng_state^=rng_state>>9;return rng_state;}
B mutate(B b,int op){
 map<int,vector<int>> c;for(int x:b)c[abs(x)].push_back(x);
 if(c.empty())return b;auto it=c.begin();advance(it,rnd()%c.size());
 int n=it->first,e=it->second.front()>0?1:-1,m=it->second.size();
 if(op==0){int nn=1+(rnd()%40);auto q=c.find(nn);int sg=e;if(q!=c.end())sg=q->second.front()>0?1:-1;
   for(int i=0;i<1+(rnd()%5);++i)b.push_back(sg*nn);
 }else if(op==1){auto q=find(b.begin(),b.end(),e*n);if(q!=b.end())b.erase(q);}
 else if(op==2){for(int&x:b)if(abs(x)==n)x=-x;}
 else if(op==3){int nn=max(1,n+(int)(rnd()%7)-3);if(nn==n)nn=min(60,n+1);bool clash=false;
   for(int x:b)if(abs(x)==nn&&(x>0?1:-1)!=e)clash=true;
   if(!clash)for(int&x:b)if(abs(x)==n)x=e*nn;
 }else if(op==4){if(m>1){auto q=find(b.begin(),b.end(),e*n);if(q!=b.end())b.erase(q);}else b.push_back(e*n);}
 return canon(b);
}
void run_seed_mode(){
 auto v=seeds(); int expected[]={128,131,34,427,540,540,869,1283};
 for(size_t i=0;i<v.size();++i){auto &p=v[i]; assert(weight(p.b)==expected[i]);
   int W=weight(p.b),M=maxlabel(p.b);
   assert(residual(p.b,p.p));
   cout<<"seed-check "<<p.name<<" W="<<W<<" max="<<M<<" factors="<<p.b.size()<<" p="<<p.p
       <<" delta="<<(W-p.p)/2<<" residual=pass pairfree=pass\n";cout.flush();
   Result z=analyze(p);report(p,z);
   if(z.best_delta<0){cout<<"D COUNTEREXAMPLE "<<p.name<<" p="<<p.p<<" Rmax="<<rstr(z.best_R)<<" delta="<<z.best_delta<<"\n";return;}
 }
}
void run_family_mode(){
 auto v=families();
 for(auto &p:v){
   if(!residual(p.b,p.p)){cerr<<"family residual check failed "<<p.name<<" W="<<weight(p.b)<<" p="<<p.p<<"\n";exit(4);}
   Result z=analyze(p);report(p,z);
   if(z.best_delta<0){cout<<"D COUNTEREXAMPLE "<<p.name<<" p="<<p.p<<" Rmax="<<rstr(z.best_R)<<" delta="<<z.best_delta<<"\n";return;}
 }
}
void run_extra_mode(){
 B b;add(b,50,1,8);add(b,11,1);add(b,10,1,4);add(b,9,1,2);add(b,8,1,2);
 add(b,7,1,2);add(b,6,1,2);add(b,3,1);add(b,2,1,2);add(b,1,-1,2);
 Profile p{"mutated-W97-full-max",canon(b),50};
 if(weight(p.b)!=520||!residual(p.b,p.p))exit(5);
 Result z=analyze(p);report(p,z);
}
bool screen_D(const Profile&pr,bool print=true){
  assert(residual(pr.b,pr.p));Z g=gp_pruned(pr.b,pr.p);auto rs=all_even_removals(pr.b);
  Z best;vector<int>bestR;bool have=false;int done=0;
  for(const auto&r:rs){
    Z h=0;if(pr.p<=weight(pr.b)-r.cost)h=gp_pruned(erase_R(pr.b,r.r),pr.p);
    Z d=g-h;++done;if(!have||d>best){have=true;best=d;bestR=r.r;}
    if(d>=0){
      if(print)cout<<"D witness profile="<<pr.name<<" W="<<weight(pr.b)<<" p="<<pr.p
        <<" R="<<rstr(r.r)<<" g="<<g<<" g_child="<<h<<" delta="<<d
        <<" delta/g="<<qstr(d,g)<<" checked-removals-before-success="<<done<<"/"<<rs.size()<<" B="<<bkey(pr.b)<<"\n";
      return true;
    }
  }
  if(print)cout<<"D COUNTEREXAMPLE profile="<<pr.name<<" W="<<weight(pr.b)<<" p="<<pr.p
    <<" g="<<g<<" max-delta="<<best<<" R="<<rstr(bestR)<<" max-delta/g="<<qstr(best,g)
    <<" removals="<<done<<" B="<<bkey(pr.b)<<"\n";
  return false;
}
void run_stress_mode(){
 auto v=stress_profiles();int pass=0;
 for(auto &pr:v){
   bool ok=screen_D(pr,true);++pass;
   cout<<"structured progress row="<<pass<<"/"<<v.size()<<" W="<<weight(pr.b)<<" p="<<pr.p
       <<" result="<<(ok?"D-held":"D-failed")<<"\n";cout.flush();
   if(!ok)return;
 }
 cout<<"structured summary exact-rows="<<pass<<" gp-calls="<<gp_calls<<" cache-hits="<<cache_hits<<"\n";
}
void run_search(int restarts,int steps,int proposals){
 auto root=seeds();int valid=0;
 int targets[]={500,500,500,500,500,500,500,1283};
 for(int q=0;q<restarts;++q){
   int target=targets[q%8];
   B base=pad_to_target(root[q%root.size()].b,target);
   int p=choose_p(base,root[q%root.size()].p);
   for(int st=0;st<steps;++st)for(int j=0;j<proposals;++j){
     bool accepted=false;
     for(int attempt=0;attempt<100&&!accepted;++attempt){
       B nb=mutate(base,int(rnd()%5));int np=choose_p(nb,p);
       if(weight(nb)<500||weight(nb)>3000||np<0||!residual(nb,np)||bkey(nb)==bkey(base))continue;
       Profile cand{"mutated-"+root[q%root.size()].name+"-target"+to_string(target)+"-"+to_string(q)+"-"+to_string(st)+"-"+to_string(j),canon(nb),np};
       ++valid;
       bool ok=screen_D(cand,true);
       cout<<"mutation progress restart="<<q+1<<"/"<<restarts<<" step="<<st+1<<"/"<<steps
           <<" proposal="<<j+1<<"/"<<proposals<<" residual-mutants="<<valid
           <<" pairfree=pass total-weight="<<weight(cand.b)<<" p="<<cand.p
           <<" result="<<(ok?"D-held":"D-failed")<<"\n";cout.flush();
       if(!ok)return;
       base=canon(nb);p=np;accepted=true;
     }
   }
 }
 cout<<"mutation summary requested-restarts="<<restarts<<" steps="<<steps<<" proposals="<<proposals
     <<" exact-mutants-screened="<<valid<<" gp-calls="<<gp_calls<<" cache-hits="<<cache_hits<<"\n";
}
int main(int argc,char**argv){
 string mode="seeds";int restarts=8,steps=1,proposals=1;
 for(int i=1;i<argc;++i){string a=argv[i];
   if(a=="-h"||a=="--help"){cout<<"usage: fmchk94 [--mode check|seeds|families|stress|search|extra|all] [--restarts N] [--steps N] [--proposals N]\n";return 0;}
   if(a=="--mode"&&i+1<argc)mode=argv[++i];
   else if(a=="--restarts"&&i+1<argc)restarts=atoi(argv[++i]);
   else if(a=="--steps"&&i+1<argc)steps=atoi(argv[++i]);
   else if(a=="--proposals"&&i+1<argc)proposals=atoi(argv[++i]);
   else if(a!="--mode"&&a!="--restarts"&&a!="--steps"&&a!="--proposals"){cerr<<"unknown argument "<<a<<"\n";return 2;}
 }
 check_pruning();
 if(mode=="check")return 0;
 if(mode=="seeds")run_seed_mode();
 else if(mode=="families")run_family_mode();
 else if(mode=="stress")run_stress_mode();
 else if(mode=="search")run_search(restarts,steps,proposals);
 else if(mode=="extra")run_extra_mode();
 else if(mode=="all"){run_seed_mode();run_family_mode();run_search(restarts,steps,proposals);}
 else {cerr<<"unknown mode "<<mode<<"\n";return 2;}
 report_global(); cout<<"global cache calls="<<gp_calls<<" hits="<<cache_hits<<" exact entries="<<gp_cache.size()<<"\n";
 return 0;
}
