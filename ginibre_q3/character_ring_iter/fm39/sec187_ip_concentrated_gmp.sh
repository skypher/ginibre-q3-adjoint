cd /home/yang/q3adjoint
python3 -u - <<'FM187PY'
import os, subprocess, sys
os.environ["TMPDIR"]="/dev/shm"
if any(a in ("-h","--help") for a in sys.argv[1:]):
    print("FM-SEC187 exact GMP interior-prefix search; --help prints this message.")
    raise SystemExit(0)
cpp = r'''
#include <gmpxx.h>
#include <algorithm>
#include <atomic>
#include <ctime>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>
using namespace std;
static string now(){time_t t=time(nullptr);tm q;gmtime_r(&t,&q);char b[32];strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&q);return b;}
static void logx(const string&s){cout<<now()<<" "<<s<<endl;}
static int wt(const vector<int>&w){int s=0;for(int z:w)s+=abs(z);return s;}
static int negs(const vector<int>&w){int s=0;for(int z:w)s+=(z<0);return s;}
static void canon(vector<int>&w){sort(w.begin(),w.end(),[](int a,int b){return abs(a)!=abs(b)?abs(a)>abs(b):a<b;});}
static string wordstr(const vector<int>&w){ostringstream o;o<<"(";for(size_t i=0;i<w.size();++i){if(i)o<<",";o<<w[i];}o<<")";return o.str();}
struct Tab{int cap;vector<mpz_class>v;Tab(int n=0):cap(n),v((size_t)(n+1)*(n+1)){}mpz_class&at(int r,int s){return v[(size_t)r*(cap+1)+s];}const mpz_class&at(int r,int s)const{return v[(size_t)r*(cap+1)+s];}};
static Tab build(vector<int>w){
 canon(w);Tab old(0);old.at(0,0)=1;int S=0;
 for(int z:w){int n=abs(z),ep=z>0?1:-1,T=S+n;Tab out(T);
  for(int s=0;s<=S;++s){vector<mpz_class>p[2];p[0].resize(S+1);p[1].resize(S+1);
   for(int a=0;a<=S;++a){for(int h=0;h<2;++h)p[h][a]=(a?p[h][a-1]:mpz_class(0));p[a&1][a]+=old.at(a,s);}
   for(int r=0;r<=T;++r){int lo=abs(r-n),hi=min(S,r+n),h=((r-n)%2+2)%2;if(lo<=hi){mpz_class q=p[h][hi];if(lo)q-=p[h][lo-1];out.at(r,s)+=q;}}
  }
  for(int r=0;r<=S;++r){vector<mpz_class>p[2];p[0].resize(S+1);p[1].resize(S+1);
   for(int b=0;b<=S;++b){for(int h=0;h<2;++h)p[h][b]=(b?p[h][b-1]:mpz_class(0));p[b&1][b]+=old.at(r,b);}
   for(int s=0;s<=T;++s){int lo=abs(s-n),hi=min(S,s+n),h=((s-n)%2+2)%2;if(lo<=hi){mpz_class q=p[h][hi];if(lo)q-=p[h][lo-1];if(ep<0)q=-q;out.at(r,s)+=q;}}
  }
  old=std::move(out);S=T;
 }
 return old;
}
using Cache=map<vector<int>,shared_ptr<Tab>>;
static shared_ptr<Tab> gettab(vector<int>w,Cache&c,mutex&m){canon(w);{lock_guard<mutex>g(m);auto i=c.find(w);if(i!=c.end())return i->second;}auto p=make_shared<Tab>(build(w));{lock_guard<mutex>g(m);auto i=c.find(w);if(i!=c.end())return i->second;if(c.size()<64)return c.emplace(w,p).first->second;return p;}}
struct Rat{int inf=0;mpz_class n=0,d=1;};
static Rat rat(mpz_class n,mpz_class d){Rat q;if(d==0){q.inf=n>0?1:n<0?-1:0;return q;}if(d<0){d=-d;n=-n;}mpz_class g;mpz_gcd(g.get_mpz_t(),n.get_mpz_t(),d.get_mpz_t());if(g!=0){n/=g;d/=g;}q.n=n;q.d=d;return q;}
static int cmp(const Rat&a,const Rat&b){if(a.inf!=b.inf)return a.inf<b.inf?-1:1;if(a.inf)return 0;mpz_class x=a.n*b.d,y=b.n*a.d;return x<y?-1:x>y?1:0;}
static string rs(const Rat&q){if(q.inf>0)return "+INF";if(q.inf<0)return "-INF";return q.n.get_str()+"/"+q.d.get_str();}
static mpz_class ab(mpz_class x){return x<0?-x:x;}
struct Metric{Rat q;int T=-1;mpz_class prefix=0,den=0;};
struct PairEval{int u=0,v=0;bool pass=true,active=false;int badT=-1;mpz_class badPrefix=0;Metric raw,tail;};
static pair<vector<int>,vector<int>> splitInterior(vector<int>C){canon(C);vector<int>A,B;int wa=0,wb=0;for(int z:C){if(wa<=wb){A.push_back(z);wa+=abs(z);}else{B.push_back(z);wb+=abs(z);}}if(wa>wb)swap(A,B);canon(A);canon(B);return {A,B};}
static PairEval evalPair(const vector<int>&w,int u,int v,Cache&cache,mutex&mx){
 vector<int>C;bool removedU=false,removedV=false;for(int z:w){if(!removedU&&z==u){removedU=true;continue;}if(!removedV&&z==v){removedV=true;continue;}C.push_back(z);}
 auto abside=splitInterior(C);vector<int>B=abside.second;B.push_back(u);B.push_back(v);canon(B);
 auto X=gettab(abside.first,cache,mx),Y=gettab(B,cache,mx);int cap=min(X->cap,Y->cap);
 vector<mpz_class>L(cap+1);for(int r=0;r<=cap;++r)for(int s=0;r+s<=cap;++s){auto&x=X->at(r,s);auto&y=Y->at(r,s);if(x!=0&&y!=0)L[r+s]+=x*y;}
 PairEval e;e.u=u;e.v=v;int parity=wt(abside.first)&1;mpz_class P=0,M=0;bool rawHave=false,tailHave=false;
 for(int T=parity;T<=cap;T+=2){
  P+=L[T];if(L[T]>0)M+=L[T];if(L[T]!=0)e.active=true;
  if(P<0){e.pass=false;if(e.badT<0){e.badT=T;e.badPrefix=P;}}
  if(M>0){Rat q=rat(P,M);if(!rawHave||cmp(q,e.raw.q)<0){e.raw={q,T,P,M};rawHave=true;}}
  else if(P<0){e.raw={rat(-1,0),T,P,0};rawHave=true;}
  mpz_class den=0;for(int j=0;j<5;++j){int h=T-2*j;if(h>=0)den+=ab(L[h]);}
  if(den>0){Rat q=rat(P,den);if(!tailHave||cmp(q,e.tail.q)<0){e.tail={q,T,P,den};tailHave=true;}}
  else if(P<0){e.tail={rat(-1,0),T,P,0};tailHave=true;}
 }
 if(!rawHave)e.raw={rat(1,0),-1,0,0};
 if(!tailHave)e.tail={rat(1,0),-1,0,0};
 return e;
}
struct Result{vector<int>w;string family;int knob=0,W=0,n=0,good=0,total=0;bool fail=false;Metric raw,tail;int ru=0,rv=0,tu=0,tv=0;vector<PairEval>bad;};
static bool valid(const vector<int>&w){return w.size()>=2&&w.size()<=40&&wt(w)<=300&&negs(w)%2==0;}
static Result analyze(vector<int>w,string fam,int knob){
 canon(w);Result r;r.w=w;r.family=fam;r.knob=knob;r.W=wt(w);r.n=w.size();
 vector<int>types;for(int z:w)if(find(types.begin(),types.end(),z)==types.end())types.push_back(z);
 vector<pair<int,int>>jobs;for(size_t i=0;i<types.size();++i)for(size_t j=i;j<types.size();++j){if(i==j&&count(w.begin(),w.end(),types[i])<2)continue;jobs.push_back({types[i],types[j]});}
 r.total=jobs.size();vector<PairEval>out(jobs.size());Cache cache;mutex mx;
 #pragma omp parallel for schedule(dynamic,1)
 for(int k=0;k<(int)jobs.size();++k)out[k]=evalPair(w,jobs[k].first,jobs[k].second,cache,mx);
 bool hr=false,ht=false;
 for(auto&e:out){if(e.pass)++r.good;else r.bad.push_back(e);
  if(!hr||cmp(e.raw.q,r.raw.q)>0){r.raw=e.raw;r.ru=e.u;r.rv=e.v;hr=true;}
  if(!ht||cmp(e.tail.q,r.tail.q)>0){r.tail=e.tail;r.tu=e.u;r.tv=e.v;ht=true;}
 }
 r.fail=(r.good==0);return r;
}
static void printMetric(const char*lab,const Metric&m,int u,int v){cout<<" "<<lab<<"="<<rs(m.q)<<"@T"<<m.T<<" pair=("<<u<<","<<v<<") prefix="<<m.prefix<<" den="<<m.den;}
static vector<int> makeword(int a,int x,int b,int y,int c,int z){vector<int>w;for(int i=0;i<x;++i)w.push_back(a);for(int i=0;i<y;++i)w.push_back(b);for(int i=0;i<z;++i)w.push_back(-c);canon(w);return w;}
static vector<int> make2neg(int a,int x,int b,int y,int c,int z,int d,int q){vector<int>w;for(int i=0;i<x;++i)w.push_back(a);for(int i=0;i<y;++i)w.push_back(b);for(int i=0;i<z;++i)w.push_back(-c);for(int i=0;i<q;++i)w.push_back(-d);canon(w);return w;}
int main(int argc,char**argv){for(int i=1;i<argc;++i)if(string(argv[i])=="-h"||string(argv[i])=="--help"){cout<<"FM-SEC187 exact prefix scan; internal families W<=300, n<=40."<<endl;return 0;}
 omp_set_dynamic(0);omp_set_num_threads(6);set<vector<int>>seen;long long tested=0,oddW=0,invalid=0;int totalFailures=0;
 map<pair<string,int>,Result>group;Result globalRaw,globalTail;bool hasRaw=false,hasTail=false;
 auto eval=[&](const string&fam,int knob,vector<int>w){canon(w);if(!valid(w)){++invalid;return;}if(!seen.insert(w).second)return;if(wt(w)%2){++oddW;return;}
  Result r=analyze(w,fam,knob);++tested;
  if(r.fail){++totalFailures;cout<<"IP_FAILURE family="<<fam<<" knob="<<knob<<" W="<<r.W<<" n="<<r.n<<" word="<<wordstr(w)<<" good="<<r.good<<"/"<<r.total<<endl;for(auto&e:r.bad)cout<<" BAD_PAIR=("<<e.u<<","<<e.v<<") first_negative_T="<<e.badT<<" prefix="<<e.badPrefix<<" raw="<<rs(e.raw.q)<<" tail="<<rs(e.tail.q)<<endl;}
  if(!hasRaw||cmp(r.raw.q,globalRaw.raw.q)<0){globalRaw=r;hasRaw=true;logx("new raw min family="+fam+" knob="+to_string(knob)+" W="+to_string(r.W)+" n="+to_string(r.n)+" score="+rs(r.raw.q));}
  if(!hasTail||cmp(r.tail.q,globalTail.tail.q)<0){globalTail=r;hasTail=true;logx("new tail min family="+fam+" knob="+to_string(knob)+" W="+to_string(r.W)+" n="+to_string(r.n)+" score="+rs(r.tail.q));}
  auto key=make_pair(fam,knob);auto it=group.find(key);if(it==group.end()||cmp(r.raw.q,it->second.raw.q)<0)group[key]=r;
  if(tested%250==0)logx("progress tested="+to_string(tested)+" odd_weight_trivial="+to_string(oddW)+" latest_family="+fam+" knob="+to_string(knob)+" latest_W="+to_string(r.W));
 };
 // The two LP counterexamples, followed by fixed-label one-negative sweeps.
 eval("LP_example_x2_y3_z14",0,makeword(1,2,2,3,3,14));
 eval("LP_example_x2_y7_z10",0,makeword(1,2,2,7,3,10));
 for(int c=3;c<=12;++c)for(int x=1;x<=8;++x)for(int y=1;y<=10;++y)for(int z=2;z<=28;z+=2){vector<int>w=makeword(1,x,2,y,c,z);if((int)w.size()<=40&&wt(w)<=300)eval("one_minus_1_2_c",z,w);}
 // Larger negative label and high-weight edge representatives.
 for(int c:{15,20,25,30})for(int x=1;x<=3;++x)for(int y=1;y<=4;++y)for(int z=2;z<=24;z+=2){auto w=makeword(1,x,2,y,c,z);if((int)w.size()<=40&&wt(w)<=300)eval("one_minus_large_c",z,w);}
 // Vary all three label classes on selected small-label triples.
 vector<tuple<int,int,int>>triples={{1,2,4},{1,2,5},{1,3,4},{1,3,5},{2,3,4},{2,4,7},{3,4,5},{3,5,8}};
 for(auto [a,b,c]:triples)for(int x=1;x<=5;++x)for(int y=1;y<=8;++y)for(int z=2;z<=18;z+=2){auto w=makeword(a,x,b,y,c,z);if((int)w.size()<=40&&wt(w)<=300)eval("varied_one_minus_labels",z,w);}
 // Two repeated minus classes; z+w is even, preserving the consumer parity.
 vector<pair<int,int>>md={{3,4},{3,6},{4,5},{4,7},{5,8},{7,11},{10,11},{14,15},{20,21}};
 for(auto [c,d]:md)for(int x=1;x<=3;++x)for(int y=1;y<=3;++y)for(int z=1;z<=10;++z)for(int q=1;q<=10;++q){if((z+q)%2)continue;auto w=make2neg(1,x,2,y,c,z,d,q);if((int)w.size()<=40&&wt(w)<=300)eval("two_minus_classes",z+q,w);}
 // Add one or two repeated large plus classes above the minus label.
 vector<vector<int>>core={{2,3,14},{2,7,10},{3,4,12},{2,5,18}};
 for(auto v:core)for(int d:{5,7,10,15,20,30})for(int m=1;m<=3;++m){auto w=makeword(1,v[0],2,v[1],3,v[2]);for(int k=0;k<m;++k)w.push_back(d);if((int)w.size()<=40&&wt(w)<=300)eval("one_large_plus_class",m,w);}
 for(auto v:core)for(int d:{5,10,20})for(int e:{12,25,35})if(e>d)for(int m=1;m<=2;++m)for(int q=1;q<=2;++q){auto w=makeword(1,v[0],2,v[1],3,v[2]);for(int k=0;k<m;++k)w.push_back(d);for(int k=0;k<q;++k)w.push_back(e);if((int)w.size()<=40&&wt(w)<=300)eval("two_large_plus_classes",m+q,w);}
 cout<<"SUMMARY tested_even_weight="<<tested<<" odd_weight_trivial="<<oddW<<" invalid="<<invalid<<" unique_seen="<<seen.size()<<" IP_failures="<<totalFailures<<endl;
 if(hasRaw){cout<<"GLOBAL_RAW family="<<globalRaw.family<<" knob="<<globalRaw.knob<<" W="<<globalRaw.W<<" n="<<globalRaw.n<<" word="<<wordstr(globalRaw.w);printMetric("best_prefix_margin",globalRaw.raw,globalRaw.ru,globalRaw.rv);cout<<endl;}
 if(hasTail){cout<<"GLOBAL_TAIL family="<<globalTail.family<<" knob="<<globalTail.knob<<" W="<<globalTail.W<<" n="<<globalTail.n<<" word="<<wordstr(globalTail.w);printMetric("best_tail_prefix_margin",globalTail.tail,globalTail.tu,globalTail.tv);cout<<endl;}
 for(auto&kv:group){auto&r=kv.second;cout<<"GROUP_MIN family="<<r.family<<" knob="<<r.knob<<" W="<<r.W<<" n="<<r.n<<" word="<<wordstr(r.w);printMetric("raw",r.raw,r.ru,r.rv);printMetric("tail",r.tail,r.tu,r.tv);cout<<endl;}
 return totalFailures?3:0;
}
'''
fd=os.memfd_create("fmsec187_gmp",0)
c=subprocess.run(["g++","-std=c++17","-O3","-fopenmp","-x","c++","-","-o",f"/proc/self/fd/{fd}","-lgmpxx","-lgmp"],input=cpp.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(fd,))
if c.returncode:
    print(c.stderr.decode(),flush=True)
    raise SystemExit(c.returncode)
p=subprocess.Popen([f"/proc/self/fd/{fd}"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,pass_fds=(fd,))
for line in p.stdout:
    print(line,end="",flush=True)
rc=p.wait()
print("memfd evaluator exit",rc,flush=True)
FM187PY
