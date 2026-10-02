import os, subprocess, sys
os.environ["TMPDIR"]="/dev/shm"
if any(a in ("-h","--help") for a in sys.argv[1:]):
    print("FM-SEC186 exact GMP layer search; --quick checks fixed seeds only.")
    raise SystemExit(0)
cpp = r'''
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
using namespace std;
static string now(){time_t t=time(nullptr);tm q;gmtime_r(&t,&q);char b[32];strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&q);return b;}
static void logx(const string&s){cout<<now()<<" "<<s<<endl;}
static int weight(const vector<int>&w){int z=0;for(int x:w)z+=abs(x);return z;}
static bool evenminus(const vector<int>&w){int n=0;for(int x:w)n+=(x<0);return n%2==0;}
static bool pf(const vector<int>&w){map<int,int>m;for(int x:w){int s=x>0?1:-1;if(m.count(abs(x))&&m[abs(x)]!=s)return false;m[abs(x)]=s;}return true;}
static void canon(vector<int>&w){sort(w.begin(),w.end(),[](int a,int b){return abs(a)!=abs(b)?abs(a)>abs(b):a<b;});}
static string show(const vector<int>&w){ostringstream o;o<<"(";for(size_t i=0;i<w.size();++i){if(i)o<<",";o<<w[i];}o<<")";return o.str();}
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
static shared_ptr<Tab> table(vector<int>w,Cache&c,mutex&m){canon(w);{lock_guard<mutex>g(m);auto i=c.find(w);if(i!=c.end())return i->second;}auto p=make_shared<Tab>(build(w));{lock_guard<mutex>g(m);auto it=c.find(w);if(it!=c.end())return it->second;if(c.size()<128)return c.emplace(w,p).first->second;return p;}}
struct Rat{int inf=0;mpz_class n=0,d=1;};
static Rat make_rat(mpz_class n,mpz_class d){Rat q;if(d==0){q.inf=n>0?1:n<0?-1:0;return q;}if(d<0){d=-d;n=-n;}mpz_class g;mpz_gcd(g.get_mpz_t(),n.get_mpz_t(),d.get_mpz_t());if(g!=0){n/=g;d/=g;}q.n=n;q.d=d;return q;}
static int compare(const Rat&a,const Rat&b){if(a.inf!=b.inf)return a.inf<b.inf?-1:1;if(a.inf)return 0;mpz_class x=a.n*b.d,y=b.n*a.d;return x<y?-1:x>b.n?1:0;}
static string rs(const Rat&q){if(q.inf==1)return "+INF";if(q.inf==-1)return "-INF";return q.n.get_str()+"/"+q.d.get_str();}
static mpz_class absv(mpz_class x){return x<0?-x:x;}
struct Met{Rat q;int t=-1;mpz_class n=0,d=0;};
static Met minmet(const vector<mpz_class>&L,int kind,const mpz_class&mass=0){Met b;bool have=false;int cap=(int)L.size()-1;for(int t=0;t<=cap;++t)if(L[t]!=0){mpz_class d=0;if(kind==0)d=mass;else if(kind==1){if(t>=2)d=max(d,absv(L[t-2]));if(t+2<=cap)d=max(d,absv(L[t+2]));}else for(int k=-2;k<=2;++k){int u=t+2*k;if(u>=0&&u<=cap)d+=absv(L[u]);}Rat q=make_rat(L[t],d);if(!have||compare(q,b.q)<0){have=true;b={q,t,L[t],d};}}if(!have)b={make_rat(0,0),-1,0,0};return b;}
struct PE{int i=0,j=0;vector<int>A,B;vector<mpz_class>L;Met raw,t2,t5;mpz_class minL=0,mass=0;bool neg=false,active=false;};
static pair<vector<int>,vector<int>> cut(vector<int>C){canon(C);vector<int>A,B;int wa=0,wb=0;for(int z:C){if(wa<=wb){A.push_back(z);wa+=abs(z);}else{B.push_back(z);wb+=abs(z);}}if(wa>wb)swap(A,B);canon(A);canon(B);return {A,B};}
static PE pair_eval(const vector<int>&w,int i,int j,Cache&c,mutex&m){vector<int>C;for(int k=0;k<(int)w.size();++k)if(k!=i&&k!=j)C.push_back(w[k]);auto ab=cut(C);vector<int>B=ab.second;B.push_back(w[i]);B.push_back(w[j]);canon(B);auto X=table(ab.first,c,m),Y=table(B,c,m);int cap=min(X->cap,Y->cap);vector<mpz_class>L(cap+1);for(int r=0;r<=cap;++r)for(int s=0;r+s<=cap;++s){auto&x=X->at(r,s);auto&y=Y->at(r,s);if(x!=0&&y!=0)L[r+s]+=x*y;}PE e;e.i=i;e.j=j;e.A=ab.first;e.B=B;e.L=L;bool first=true;for(auto&x:L)if(x!=0){e.active=true;if(first||x<e.minL){e.minL=x;first=false;}if(x<0)e.neg=true;if(x>0)e.mass+=x;}e.raw=minmet(L,0,e.mass);e.t2=minmet(L,1);e.t5=minmet(L,2);return e;}
struct Sum{vector<int>w;int W=0,good=0,total=0,activePairs=0;bool pairfree=true,fail=false;PE best[3];};
static Sum analyze(vector<int>w,const string&name,bool progress=true){canon(w);Sum s;s.w=w;s.W=weight(w);s.pairfree=pf(w);s.total=(int)w.size()*(int(w.size())-1)/2;if(!s.total)return s;Cache c;mutex m;vector<PE>v(s.total);atomic<int>done(0);mutex pm;
 #pragma omp parallel for schedule(dynamic,1)
 for(int i=0;i<(int)w.size();++i)for(int j=i+1;j<(int)w.size();++j){PE x=pair_eval(w,i,j,c,m);int ix=i*(int)w.size()-i*(i+1)/2+(j-i-1);v[ix]=std::move(x);int d=++done;if(progress&&(d%100==0||d==s.total)){lock_guard<mutex>g(pm);logx("pairs "+name+" "+to_string(d)+"/"+to_string(s.total));}}
 bool h[3]={false,false,false};for(auto&x:v){if(!x.neg)++s.good;if(x.active)++s.activePairs; if(!x.active)continue;Met z[3]={x.raw,x.t2,x.t5};for(int k=0;k<3;++k){if(!h[k]||compare(z[k].q,k==0?s.best[k].raw.q:k==1?s.best[k].t2.q:s.best[k].t5.q)>0){s.best[k]=x;h[k]=true;}}}s.fail=s.good==0;if(s.activePairs==0)for(int k=0;k<3;++k){if(k==0)s.best[k].raw.q=make_rat(1,0);else if(k==1)s.best[k].t2.q=make_rat(1,0);else s.best[k].t5.q=make_rat(1,0);}return s;}
static string metstr(const Met&m){return rs(m.q)+"@t="+to_string(m.t);}
static void report(const Sum&s,const string&name){cout<<"RESULT "<<name<<" W="<<s.W<<" n="<<s.w.size()<<" pairfree="<<s.pairfree<<" good="<<s.good<<"/"<<s.total<<" active_profiles="<<s.activePairs<<" all_pairs_fail="<<s.fail<<endl;for(int k=0;k<3;++k){const PE&p=s.best[k];const Met&m=k==0?p.raw:k==1?p.t2:p.t5;cout<<"  "<<(k==0?"raw":k==1?"tail2":"tail5")<<" best="<<metstr(m)<<" pair=("<<s.w[p.i]<<","<<s.w[p.j]<<") minL="<<p.minL.get_str()<<" positive_mass="<<p.mass.get_str()<<endl;}}
static bool valid(const vector<int>&w){return w.size()>=2&&w.size()<=30&&weight(w)<=300&&evenminus(w);}
static vector<int>mutate(vector<int>w,mt19937_64&r,bool allowPair){int op=r()%5;if(op==0&&!w.empty()){int i=r()%w.size(),sg=w[i]>0?1:-1,n=abs(w[i])+(r()%2?1:-1)*(1+r()%2);if(n<1)n=1;if(n>300)n=300;w[i]=sg*n;}else if(op==1&&w.size()>1){int i=r()%w.size(),j=r()%w.size();while(j==i)j=r()%w.size();w[i]=-w[i];w[j]=-w[j];}else if(op==2&&w.size()<=28){int z=(r()%2?1:-1)*(1+r()%60);w.push_back(z);w.push_back(z);}else if(op==3&&w.size()>=4){int i=r()%w.size(),j=r()%w.size();while(j==i)j=r()%w.size();if((w[i]>0)==(w[j]>0)){if(i>j)swap(i,j);w.erase(w.begin()+j);w.erase(w.begin()+i);}}else if(!w.empty()){int i=r()%w.size();w[i]=(w[i]>0?1:-1)*(1+r()%90);}if(!allowPair&&!pf(w))return {};if(!evenminus(w))return {};canon(w);if(!valid(w))return {};return w;}
static vector<int>random_word(mt19937_64&r,int n,bool wantpf){vector<int>w;for(int i=0;i<n;++i)w.push_back((r()%2?1:-1)*(1+r()%max(2,250/n)));if(!evenminus(w))w[0]=-w[0];if(wantpf&&!pf(w)){for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(abs(w[i])==abs(w[j])&&w[i]*w[j]<0)w[j]=abs(w[j])+1;}canon(w);return w;}
static double approx(const Rat&q){if(q.inf>0)return 1e90;if(q.inf<0)return -1e90;return q.n.get_d()/q.d.get_d();}
int main(int argc,char**argv){bool quick=false;for(int i=1;i<argc;++i){string a=argv[i];if(a=="-h"||a=="--help"){cout<<"FM-SEC186 exact GMP layer search; --quick checks fixed seeds only."<<endl;return 0;}if(a=="--quick")quick=true;else{cerr<<"unknown option "<<a<<endl;return 2;}}omp_set_dynamic(0);omp_set_num_threads(6);
 vector<pair<string,vector<int>>>seeds={
 {"W140_pairbearing",{17,-14,14,14,13,-11,11,11,-10,9,7,4,-3,2}},
 {"W132_pairfree",{21,17,-16,-16,-13,-13,-7,-7,-6,-6,4,-3,-3}},
 {"budget",{-1,-1,-1,-1,-1,-1,-2,-2,-2,-3,-3,-3,-4,-7}},
 {"F1",{1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8}},
 {"SEC181_W60a",{1,-2,-4,-4,-4,1,3,3,3,3,3,3,5,5,5,5,6}},
 {"SEC181_W62b",{1,-2,-4,-4,-4,1,3,3,3,3,3,5,5,5,5,5,6}},
 {"SEC181_W60c",{1,-2,-4,-4,-4,1,1,2,3,3,3,3,3,4,4,5,5,8}},
 {"sign_min",{-1,2,3,4,-5,6,7,8}},
 {"opposite_classes",{1,1,-1,2,2,-2,3,3,-3,4,-4,5,5,-5,6,-6,7,8,-9,10}},
 {"multi_large_pairs",{-30,30,-27,27,-24,24,-21,21,-18,18,-15,15,9,7,5,3}}
 };
 for(int k=6;k<=22;++k){vector<int>w;for(int x=1;x<=k;++x)w.push_back(-x);w.push_back((k%2?-1:1)*(k+3));seeds.push_back({"run"+to_string(k),w});}
 for(int n=1;n<=12;++n){vector<int>w={-n,n}; for(int k=1;k<n;++k){w.push_back(-k);w.push_back(k);} seeds.push_back({"opp24_prefix"+to_string(n),w});}
 {vector<int>w;for(int n=1;n<=14;++n){w.push_back(-n);w.push_back(n);}w.push_back(1);w.push_back(1);seeds.push_back({"opp30",w});}
 map<vector<int>,Sum>done;map<vector<int>,string>names;int serial=0;
 auto eval=[&](string n,vector<int>w){canon(w);if(!valid(w)||done.count(w))return;logx("begin "+n+" W="+to_string(weight(w))+" n="+to_string(w.size()));Sum s=analyze(w,n,true);report(s,n);done.emplace(s.w,s);names[s.w]=n;++serial;if(s.fail)cout<<"ALL_PAIR_FAILURE "<<show(s.w)<<" W="<<s.W<<endl;};
 for(auto&x:seeds)eval(x.first,x.second);
 if(!quick){mt19937_64 rng(186103);for(int k=2;k<=3;++k){vector<int>a={17,-14,14,14,13,-11,11,11,-10,9,7,4,-3,2};for(int&z:a)z*=k;eval("scaled_pairbearing_"+to_string(k),a);vector<int>b={21,17,-16,-16,-13,-13,-7,-7,-6,-6,4,-3,-3};for(int&z:b)z*=k;eval("scaled_pairfree_"+to_string(k),b);}for(int k=0;k<18;++k)eval("random_"+to_string(k),random_word(rng,8+rng()%10,k%2==0));for(int k=0;k<12;++k)eval("random_long_"+to_string(k),random_word(rng,18+rng()%13,k%2==0));
  // Walk each objective from the currently smallest margin, two restarts and five mutation rounds.
  for(int obj=0;obj<3;++obj)for(int restart=0;restart<2;++restart){auto it=done.begin();for(auto q=done.begin();q!=done.end();++q){Rat x=obj==0?q->second.best[0].raw.q:obj==1?q->second.best[1].t2.q:q->second.best[2].t5.q;Rat y=obj==0?it->second.best[0].raw.q:obj==1?it->second.best[1].t2.q:it->second.best[2].t5.q;if(compare(x,y)<0)it=q;}vector<int>cur=it->first;if(restart){it=done.begin();advance(it,rng()%done.size());cur=it->first;}for(int step=0;step<5;++step){vector<int>cand=mutate(cur,rng,true);if(cand.empty()||done.count(cand))cand=mutate(cur,rng,true);if(!cand.empty()&&!done.count(cand))eval("walk"+to_string(obj)+"_"+to_string(restart)+"_"+to_string(step),cand);auto curIt=done.find(cur);if(!cand.empty()&&done.count(cand)){auto c=done.find(cand);Rat x=obj==0?c->second.best[0].raw.q:obj==1?c->second.best[1].t2.q:c->second.best[2].t5.q;Rat y=obj==0?curIt->second.best[0].raw.q:obj==1?curIt->second.best[1].t2.q:curIt->second.best[2].t5.q;if(compare(x,y)<=0)cur=cand;}logx("walk objective="+to_string(obj)+" restart="+to_string(restart)+" step="+to_string(step+1)+"/5 evaluated="+to_string(done.size()));}}
 }
 for(int obj=0;obj<3;++obj){auto it=done.begin();for(auto q=done.begin();q!=done.end();++q){Rat x=obj==0?q->second.best[0].raw.q:obj==1?q->second.best[1].t2.q:q->second.best[2].t5.q;Rat y=obj==0?it->second.best[0].raw.q:obj==1?it->second.best[1].t2.q:it->second.best[2].t5.q;if(compare(x,y)<0)it=q;}cout<<"GLOBAL_EXTREME "<<(obj==0?"raw":obj==1?"tail2":"tail5")<<" score="<<rs(obj==0?it->second.best[0].raw.q:obj==1?it->second.best[1].t2.q:it->second.best[2].t5.q)<<" W="<<it->second.W<<" n="<<it->first.size()<<" pairfree="<<it->second.pairfree<<" name="<<names[it->first]<<" word="<<show(it->first)<<endl;report(it->second,"GLOBAL");}
 int fail=0;for(auto&x:done)if(x.second.fail)++fail;cout<<"TOTAL_LISTS "<<done.size()<<" ALL_PAIR_FAILURES "<<fail<<endl;return fail?3:0;}
'''
fd=os.memfd_create("fmsec186_gmp",0)
c=subprocess.run(["g++","-std=c++17","-O3","-fopenmp","-x","c++","-","-o",f"/proc/self/fd/{fd}","-lgmpxx","-lgmp"],input=cpp.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(fd,))
if c.returncode: print(c.stderr.decode(),flush=True);raise SystemExit(c.returncode)
args=["--quick"] if "--quick" in sys.argv else []
p=subprocess.Popen([f"/proc/self/fd/{fd}",*args],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,pass_fds=(fd,))
for line in p.stdout: print(line,end="",flush=True)
rc=p.wait();print("memfd evaluator exit",rc,flush=True)
