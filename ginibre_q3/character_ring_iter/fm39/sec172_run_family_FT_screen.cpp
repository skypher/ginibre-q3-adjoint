#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <gmpxx.h>

using namespace std;
using Z = mpz_class;
using B = vector<int>;
using U64 = uint64_t;

struct PairIndex { int i, j, priority; };
struct Top { int i=-1, j=-1, weight=-1, largest=-1; };
struct ScanStats {
  long long profiles=0, residual=0, flips=0, noflip=0, top_fail=0;
  long long all_removal_fail=0, channel_identity_checks=0;
};
struct MaxRatio {
  bool set=false;
  mpq_class value;
  string name;
  int p=0, s=0, k=0, step=0, minus_count=0;
  Z parent=0, child=0;
  int ti=-1,tj=-1;
  B b;
};

static long long entry_calls=0;
string stamp();

int W(const B& b){int z=0;for(int x:b)z+=abs(x);return z;}
int Mx(const B& b){int z=0;for(int x:b)z=max(z,abs(x));return z;}
int nlarge(const B& b){int z=0;for(int x:b)if(abs(x)>=3)++z;return z;}
int sigma_of(const B& b){int m=0;for(int x:b)if(x<0)++m;return (m&1)?-1:1;}
bool pairfree(const B& b){
  map<int,int> sign;
  for(int x:b){int n=abs(x),e=(x>0)?1:-1;auto it=sign.find(n);
    if(it!=sign.end()&&it->second!=e)return false;sign[n]=e;
  }
  return true;
}
bool residual(const B& b,int p){
  int w=W(b),mx=Mx(b),d=(w-p)/2;
  return pairfree(b)&&p>=6&&p>=mx&&w>=p&&((w-p)&1)==0&&
         d>=8&&mx<=d&&nlarge(b)>=2;
}
uint64_t key(int s,int t){return (uint64_t(uint32_t(s))<<32)|uint32_t(t);}
int ks(uint64_t q){return int(q>>32);}
int kt(uint64_t q){return int(q&0xffffffffU);}

// Exact coefficient [U_A(X) U_B(Y)] in
// product_f (U_|f|(X)+sign(f)U_|f|(Y)).
Z entry(B lab,int A,int Bt){
  ++entry_calls;
  sort(lab.begin(),lab.end(),[](int a,int b){return abs(a)>abs(b);});
  int rem=0;for(int x:lab)rem+=abs(x);
  unordered_map<U64,Z> cur,nx;cur.reserve(128);cur[key(0,0)]=1;
  static auto hb_last=chrono::steady_clock::now();
  int factor_i=0;
  for(int f:lab){
    int n=abs(f),eps=(f>0)?1:-1;rem-=n;nx.clear();
    nx.reserve(cur.size()*2+32);
    for(const auto& kv:cur){
      int s=ks(kv.first),t=kt(kv.first);const Z& v=kv.second;
      if(abs(t-Bt)<=rem)
        for(int s2=abs(s-n);s2<=s+n;s2+=2)
          if(abs(s2-A)<=rem)nx[key(s2,t)]+=v;
      if(abs(s-A)<=rem)
        for(int t2=abs(t-n);t2<=t+n;t2+=2)
          if(abs(t2-Bt)<=rem){if(eps>0)nx[key(s,t2)]+=v;else nx[key(s,t2)]-=v;}
    }
    for(auto it=nx.begin();it!=nx.end();){if(it->second==0)it=nx.erase(it);else ++it;}
    cur.swap(nx);++factor_i;
    auto hb_now=chrono::steady_clock::now();
    if(chrono::duration_cast<chrono::seconds>(hb_now-hb_last).count()>=15){
      cerr<<"[heartbeat "<<stamp()<<"] entries="<<entry_calls<<" factor="<<factor_i<<"/"<<lab.size()
          <<" A="<<A<<" B="<<Bt<<" rem="<<rem<<" states="<<cur.size()<<"\n";cerr.flush();hb_last=hb_now;
    }
  }
  auto it=cur.find(key(A,Bt));return it==cur.end()?Z(0):it->second;
}
Z gp(const B& b,int p){return entry(b,p,0);}
string showB(const B& b){
  string s="(";for(size_t i=0;i<b.size();++i){if(i)s+=",";s+=(b[i]>0?"+":"")+to_string(b[i]);}return s+")";
}
string qstr(const Z& a,const Z& b){
  if(b==0)return "undefined";
  mpq_class q(a,b);q.canonicalize();
  return q.get_num().get_str()+"/"+q.get_den().get_str();
}
string stamp(){
  time_t t=time(nullptr);char b[16];tm* q=localtime(&t);strftime(b,sizeof(b),"%H:%M:%S",q);return b;
}
Top top_pair(const B& b){
  Top t;
  for(size_t i=0;i<b.size();++i)for(size_t j=i+1;j<b.size();++j){
    int a=abs(b[i]),c=abs(b[j]);if((a+c)&1)continue;
    int w=a+c,m=max(a,c);
    if(w>t.weight||(w==t.weight&&m>t.largest)){t={(int)i,(int)j,w,m};}
  }
  return t;
}
vector<PairIndex> flip_pairs(const B& L,const Top& top,int bsize){
  int pidx=(int)L.size()-1;vector<PairIndex> v;
  for(int i=0;i<(int)L.size();++i)for(int j=i+1;j<(int)L.size();++j){
    int pr=2;
    if(i<bsize&&j<bsize&&((i==top.i&&j==top.j)||(i==top.j&&j==top.i)))pr=0;
    else if(i==pidx||j==pidx)pr=1;
    v.push_back({i,j,pr});
  }
  stable_sort(v.begin(),v.end(),[](const PairIndex&a,const PairIndex&b){
    if(a.priority!=b.priority)return a.priority<b.priority;
    return tie(a.i,a.j)<tie(b.i,b.j);
  });
  return v;
}
struct FlipResult {bool noflip=false;int pairs=0;Z maxD=0;int mi=-1,mj=-1;};
FlipResult scan_flips(const B& b,int p,int sigma){
  B L=b;L.push_back(sigma*p);Top tp=top_pair(b);
  auto prs=flip_pairs(L,tp,(int)b.size());
  FlipResult out;out.noflip=true;out.pairs=(int)prs.size();bool have=false;
  for(const auto& uv:prs){
    B C=L;int u=L[uv.i],v=L[uv.j];
    C.erase(C.begin()+uv.j);C.erase(C.begin()+uv.i);
    Z d=entry(C,abs(u),abs(v));if(v<0)d=-d;
    if(!have||d>out.maxD){have=true;out.maxD=d;out.mi=u;out.mj=v;}
    if(d>=0){out.noflip=false;out.maxD=d;out.mi=u;out.mj=v;return out;}
  }
  return out;
}
struct Removal {vector<int> r;int cost;};
vector<Removal> removals(const B& b){
  map<int,int> count;for(int x:b)++count[x];
  vector<int> cls;for(auto [x,c]:count)cls.push_back(x);
  vector<Removal> r;
  for(int x:cls)if(abs(x)%2==0)r.push_back({{x},abs(x)});
  for(size_t i=0;i<cls.size();++i)for(size_t j=i;j<cls.size();++j){
    int x=cls[i],y=cls[j];if((abs(x)&1)!=(abs(y)&1))continue;
    if(i==j&&count[x]<2)continue;
    r.push_back({{x,y},abs(x)+abs(y)});
  }
  sort(r.begin(),r.end(),[](const Removal&a,const Removal&b){
    if(a.cost!=b.cost)return a.cost<b.cost;return a.r<b.r;
  });
  return r;
}
B minus_r(B b,const vector<int>& r){
  for(int x:r){auto it=find(b.begin(),b.end(),x);if(it==b.end())abort();b.erase(it);}
  return b;
}
ScanStats stats;
MaxRatio max_ratio;
bool stopped=false;
int smax=200,sstep=0,seconds_limit=210;
bool next_numeric_only=false;
bool large_start_probe=false;
bool high_start_probe=false;

void print_full_fallback(const B& b,int p,int sigma,const Top& tp,const Z& parent,const Z& tpchild){
  cout<<"FT_CANDIDATE B="<<showB(b)<<" p="<<p<<" sigma="<<sigma
      <<" W="<<W(b)<<" delta="<<(W(b)-p)/2<<" top="<<showB({b[tp.i],b[tp.j]})
      <<" g="<<parent<<" gTopChild="<<tpchild<<"\n";
  bool dpass=false;Z maxD;vector<int> bestR;bool have=false;
  for(const auto& r:removals(b)){
    B c=minus_r(b,r.r);Z child=(p<=W(c))?gp(c,p):Z(0);Z d=parent-child;
    if(!have||d>maxD){have=true;maxD=d;bestR=r.r;}
    if(d>=0)dpass=true;
  }
  cout<<"ALL_REMOVALS count="<<removals(b).size()<<" D="<<(dpass?"pass":"FAIL")
      <<" maxDelta="<<maxD<<" maxR="<<showB(bestR)<<"\n";
  if(!dpass)++stats.all_removal_fail;
  B L=b;L.push_back(sigma*p);Z phiParent=entry(L,0,0);long long nch=0;
  bool anyChannelNonnegative=false;Z chMin=0,chMax=0;bool first=true;
  for(size_t i=0;i<b.size();++i)for(size_t j=i+1;j<b.size();++j){
    int a=abs(b[i]),c=abs(b[j]);if((a+c)&1)continue;
    B C=L;C.erase(C.begin()+j);C.erase(C.begin()+i);
    int eps=(b[i]>0)==(b[j]>0)?1:-1;Z sum=0;
    for(int h=abs(a-c);h<=a+c;h+=2){
      Z val;
      if(h==0)val=(eps>0)?2*entry(C,0,0):Z(0);
      else {B child=C;child.push_back(eps*h);val=entry(child,0,0);}
      sum+=val;++nch;if(val>=0)anyChannelNonnegative=true;
      if(first||val<chMin)chMin=val;if(first||val>chMax)chMax=val;first=false;
    }
    ++stats.channel_identity_checks;
    if(sum!=phiParent){
      cout<<"CHANNEL_IDENTITY_MISMATCH pair="<<b[i]<<","<<b[j]
          <<" parentPhi="<<phiParent<<" channelSum="<<sum<<"\n";return;
    }
  }
  cout<<"CHANNELS identities=pass terms="<<nch<<" any_nonnegative_term="
      <<(anyChannelNonnegative?"yes":"no")<<" min_term="<<chMin<<" max_term="<<chMax
      <<" parentPhi="<<phiParent<<"\n";
}
void process_profile(const B& b,int p,const string& name,int s,int k,int step,int m){
  ++stats.profiles;
  if(!residual(b,p))return;
  ++stats.residual;
  int sigma=sigma_of(b);
  FlipResult f=scan_flips(b,p,sigma);
  if(!f.noflip){++stats.flips;return;}
  ++stats.noflip;
  Z parent=gp(b,p);Top tp=top_pair(b);
  if(tp.i<0){cerr<<"no TopPair for residual B\n";abort();}
  B childB=b;int x=b[tp.i],y=b[tp.j];
  childB.erase(childB.begin()+tp.j);childB.erase(childB.begin()+tp.i);
  Z child=(p<=W(childB))?gp(childB,p):Z(0);
  cout<<"NOFLIP "<<stamp()<<" name="<<name<<" W="<<W(b)<<" p="<<p
      <<" delta="<<(W(b)-p)/2<<" signsMinus="<<m<<" factors="<<b.size()
      <<" maxD="<<f.maxD<<" at=("<<f.mi<<","<<f.mj<<")"
      <<" B="<<showB(b)<<" TopPair=("<<x<<","<<y<<") g="<<parent
      <<" child="<<child<<" ratio="<<qstr(child,parent)<<"\n";cout.flush();
  if(parent>0){
    mpq_class q(child,parent);q.canonicalize();
    if(!max_ratio.set||q>max_ratio.value){
      max_ratio.set=true;max_ratio.value=q;max_ratio.name=name;max_ratio.p=p;
      max_ratio.s=s;max_ratio.k=k;max_ratio.step=step;max_ratio.minus_count=m;
      max_ratio.parent=parent;max_ratio.child=child;max_ratio.ti=x;max_ratio.tj=y;max_ratio.b=b;
    }
  }else{
    cout<<"NONPOSITIVE_PARENT noflip name="<<name<<" p="<<p<<" g="<<parent<<"\n";
  }
  if(child>parent){
    ++stats.top_fail;
    cout<<"TOPPAIR_FAIL exact name="<<name<<" B="<<showB(b)<<" p="<<p
        <<" g="<<parent<<" child="<<child<<" delta="<<(parent-child)<<"\n";
    print_full_fallback(b,p,sigma,tp,parent,child);
    stopped=true;
  }
}
vector<vector<int>> all_positions(int k,int m){
  vector<vector<int>> out;vector<int> cur;
  function<void(int,int)> rec=[&](int lo,int left){
    if(left==0){out.push_back(cur);return;}
    for(int i=lo;i<=k-left;++i){cur.push_back(i);rec(i+1,left-1);cur.pop_back();}
  };
  rec(0,m);return out;
}
vector<vector<int>> sampled_positions(int k,int m){
  set<vector<int>> q;
  vector<int> a;for(int i=0;i<m;++i)a.push_back(i);q.insert(a);
  a.clear();for(int i=k-m;i<k;++i)a.push_back(i);q.insert(a);
  a.clear();int start=(k-m)/2;for(int i=start;i<start+m;++i)a.push_back(i);q.insert(a);
  a.clear();for(int j=0;j<m;++j)a.push_back((j*(k-1))/(m-1));q.insert(a);
  a.clear();for(int j=0;j<m;++j)a.push_back((j%2==0)?j/2:k-1-j/2);sort(a.begin(),a.end());q.insert(a);
  return vector<vector<int>>(q.begin(),q.end());
}
void progress(const chrono::steady_clock::time_point& t0,int done,int total){
  auto sec=chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-t0).count();
  cout<<"PROGRESS "<<stamp()<<" done="<<done<<" residual="<<stats.residual
      <<" flips="<<stats.flips<<" noflip="<<stats.noflip
      <<" TopPair_fail="<<stats.top_fail<<" elapsed_s="<<sec<<"\n";cout.flush();
}
void run(){
  vector<int> starts=next_numeric_only?vector<int>{1}:(large_start_probe?vector<int>{20,100,200}:(high_start_probe?vector<int>{100,200}:vector<int>{1,2,20,40,100,200}));
  if(sstep>0){vector<int> q;for(int s=1;s<=smax;s+=sstep)q.push_back(s);if(q.empty()||q.back()!=smax)q.push_back(smax);starts=q;}
  else {starts.erase(remove_if(starts.begin(),starts.end(),[&](int x){return x>smax;}),starts.end());}
  int estimated=0;
  for(int k=8;k<=20;++k)for(int st:{1,2})for(int m:{2,4,6})
    estimated+=starts.size()*(k<=14?1:1);
  cout<<"START "<<stamp()<<" starts=";for(int s:starts)cout<<s<<",";cout
      <<" k="<<((large_start_probe||high_start_probe)?"8,12,16,20":"8..20")
      <<" step="<<(next_numeric_only?"2":"1,2")
      <<" p="<<(next_numeric_only?"max+1":((large_start_probe||high_start_probe)?"max,max+step,max+1":"max,max+step"))
      <<" minus=2,4,6 seconds_limit="<<seconds_limit<<"\n";
  auto t0=chrono::steady_clock::now();int done=0,last_report=0;
  vector<int> steps=next_numeric_only?vector<int>{2}:vector<int>{1,2};
  vector<int> lens;if(large_start_probe||high_start_probe)lens={8,12,16,20};else for(int k=8;k<=20;++k)lens.push_back(k);
  for(int s:starts)for(int k:lens)for(int step:steps)for(int m:{2,4,6}) {
    vector<vector<int>> masks=(large_start_probe||high_start_probe)?sampled_positions(k,m):((k<=14)?all_positions(k,m):sampled_positions(k,m));
    for(const auto& neg:masks){
      if(stopped)break;
      auto sec=chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-t0).count();
      if(sec>=seconds_limit){stopped=true;break;}
      B b;for(int i=0;i<k;++i){int sign=binary_search(neg.begin(),neg.end(),i)?-1:1;b.push_back(sign*(s+i*step));}
      int last=s+(k-1)*step;vector<int> ps=next_numeric_only?vector<int>{last+1}:vector<int>{last,last+step};
      if(large_start_probe||high_start_probe){set<int> z{last,last+step,last+1};ps.assign(z.begin(),z.end());}
      for(int p:ps){
        ++done;
        if(residual(b,p)){
          string name="k"+to_string(k)+"-s"+to_string(s)+"-step"+to_string(step)+"-m"+to_string(m)+"-p"+to_string(p);
          process_profile(b,p,name,s,k,step,m);
        }else ++stats.profiles;
        if(done-last_report>=250||sec>=seconds_limit){progress(t0,done,0);last_report=done;}
        if(stopped)break;
      }
    }
  }
  cout<<"SUMMARY "<<stamp()<<" profiles="<<stats.profiles<<" residual="<<stats.residual
      <<" flip_descent="<<stats.flips<<" noflip="<<stats.noflip
      <<" TopPair_fail="<<stats.top_fail<<" D_fail_after_TP="<<stats.all_removal_fail
      <<" channel_pair_checks="<<stats.channel_identity_checks<<" entry_calls="<<entry_calls
      <<" elapsed_s="<<chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-t0).count()
      <<" stopped_at_limit="<<(stopped?"yes":"no")<<"\n";
  if(max_ratio.set)cout<<"MAX_CHILD_PARENT ratio="<<max_ratio.value.get_num()<<"/"<<max_ratio.value.get_den()
      <<" name="<<max_ratio.name<<" p="<<max_ratio.p<<" W="<<W(max_ratio.b)
      <<" s="<<max_ratio.s<<" k="<<max_ratio.k<<" step="<<max_ratio.step
      <<" minus="<<max_ratio.minus_count<<" TopPair=("<<max_ratio.ti<<","<<max_ratio.tj<<")"
      <<" g="<<max_ratio.parent<<" child="<<max_ratio.child<<" B="<<showB(max_ratio.b)<<"\n";
  cout<<"CROSSCHECK requested only if TOPPAIR_FAIL occurs; named independent evaluators were not called otherwise.\n";
}
void check_local(){
  B b;for(int i=0;i<11;++i){int n=40+2*i;b.push_back((i==0||i==2)?-n:n);}
  int p=62;
  if(!residual(b,p))abort();
  FlipResult f=scan_flips(b,p,sigma_of(b));
  Z g=gp(b,p);Top t=top_pair(b);B c=b;int x=b[t.i],y=b[t.j];c.erase(c.begin()+t.j);c.erase(c.begin()+t.i);Z h=gp(c,p);
  if(!f.noflip||f.pairs!=66||f.maxD!=-34349665||g!=Z("453207534222864")||h!=Z("179646349605"))abort();
  cout<<"MECH157_CHECK W="<<W(b)<<" p="<<p<<" flips="<<f.pairs
      <<" noflip="<<(f.noflip?"yes":"no")<<" maxD="<<f.maxD
      <<" g="<<g<<" TopPair=("<<x<<","<<y<<") child="<<h
      <<" ratio="<<qstr(h,g)<<"\n";
}
int main(int argc,char**argv){
  for(int i=1;i<argc;++i){
    string a=argv[i];
    if(a=="--smax"&&i+1<argc)smax=atoi(argv[++i]);
    else if(a=="--sstep"&&i+1<argc)sstep=atoi(argv[++i]);
    else if(a=="--seconds"&&i+1<argc)seconds_limit=atoi(argv[++i]);
    else if(a=="--next-numeric-only")next_numeric_only=true;
    else if(a=="--large-start-probe")large_start_probe=true;
    else if(a=="--high-start-probe")high_start_probe=true;
    else if(a=="--check-only"){check_local();return 0;}
  }
  check_local();
  run();
  return 0;
}
