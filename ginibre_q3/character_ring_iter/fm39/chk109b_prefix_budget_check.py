import os, subprocess, shlex

CPP = r'''#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using I = long long;
using U = unsigned long long;
using W = __uint128_t;

namespace verify {
constexpr int K=12, S=26, Q=1<<20;
constexpr U BQ=1ULL<<40;
int threads=24;

void stamp(){ time_t t=time(nullptr); tm g; gmtime_r(&t,&g); char b[32]; strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&g); cout<<b; }
I tri(I n){ return n<0?0:(n+1)*(n+2)/2; }
I rect(I n,I l,I h){ return tri(n)-tri(n-l)-tri(n-h)+tri(n-l-h); }
I band(I u,I v,I l,I h,I t){
 return rect((v-u+t)/2+h-1,l,h)-rect((v-u-t)/2+h-2,l,h)-rect((t-u-v)/2-1,l,h);
}
struct Prof {
 bool ok=false;
 array<U,S+1> v;
 Prof(){ v.fill(ULLONG_MAX); }
 void take(const vector<I>& a){
  if(a[0]<=0) return;
  ok=true;
  for(int s=0;s<=S;s++){
   assert(a[s]>=0);
   U z=(U)(W(Q)*(U)a[s]/(U)a[0]);
   v[s]=min(v[s],z);
  }
 }
};
struct Pack { array<Prof,3> small; Prof triple; };
vector<vector<int>> tuples, prefixes;
unordered_map<int,int> index_of;
vector<Pack> profiles[2];
Prof H[2], H4[2];

int key_of(const vector<int>& a){ int z=1; for(int n:a) z=z*13+n; return z; }
void enumerate(vector<int>& a,int left,int lo){
 if(left==0){ tuples.push_back(a); if((int)a.size()==7) prefixes.push_back(a); return; }
 for(int n=lo;n<=K;n++){ a.push_back(n); enumerate(a,left-1,n); a.pop_back(); }
}
bool even_tuple(const vector<int>& a){ for(int n:a) if(n&1) return false; return true; }
vector<I> row(const vector<int>& a){
 vector<I> f(1,1);
 for(int n:a){
  vector<I> z(f.size()+n,0);
  for(int j=0;j<(int)f.size();j++) if(f[j])
   for(int t=abs(j-n);t<=j+n;t+=2) z[t]+=f[j];
  f.swap(z);
 }
 return f;
}
void make_H(int mode){
 int L=mode?14:13;
 vector<int> ls,hs;
 for(int l=L+1;l<(mode?3:2)*(L+1);l+=mode?2:1) ls.push_back(l);
 for(int h=1;h<(mode?3:2)*(L+1);h+=mode?2:1) hs.push_back(h);
 uint64_t shapes=0;
 for(int l:ls) for(int h:hs) for(int off=1-h;off<l;off++) for(int par=0;par<(mode?1:2);par++){
  int u=max({0,L-h+1-2*off,-2*off});
  u+=(par-u%2+2)%2;
  int v=u+2*off;
  vector<I> a(S+1);
  for(int s=0;s<=S;s++) a[s]=band(u,v,l,h,2*s);
  if(a[0]>0){ H[mode].take(a); shapes++; }
 }
 assert(H[mode].ok);
 stamp(); cout<<" H-profile mode="<<mode<<" accepted_tiles="<<shapes<<" coords";
 for(U v:H[mode].v) cout<<" "<<v;
 cout<<endl;
}
I rectangle_coefficient(int a,int b,int c,int l,int h){
 return rect(a+h-1,l,h)-rect(-b+h-2,l,h)-rect(c-1,l,h);
}
struct Tile { Prof p; uint64_t shapes=0; };
void make_H4(int mode){
 int L=mode?14:13;
 vector<int> lengths;
 for(int l=L+1;l<(mode?3:2)*(L+1);l+=mode?2:1) lengths.push_back(l);
 vector<Tile> tiles(lengths.size()*lengths.size());
 atomic<int> finished{0};
 #pragma omp parallel for schedule(dynamic)
 for(int job=0;job<(int)tiles.size();job++){
  int l=lengths[job/(int)lengths.size()], h=lengths[job%(int)lengths.size()];
  int R=max({S,l-1,h-1});
  Tile &out=tiles[job];
  for(int a=1-h;a<=l-1+2*S;a++){
   for(int b=max(1-l,-a);b<=h-1+2*S;b++){
    if(mode && (a-b)%2) continue;
    for(int c=-R;c<=min({a,b,l+h-2});c++){
     if(mode && (a-c)%2) continue;
     I den=rectangle_coefficient(a,b,c,l,h);
     if(den==0) continue;
     assert(den>0);
     I pref[2*S+2]={};
     for(int z=-S;z<=S;z++){
      I val=(a+b+2*z<0)?0:rectangle_coefficient(a+z,b+z,c+z,l,h);
      assert(val>=0);
      pref[z+S+1]=pref[z+S]+val;
     }
     vector<I> vals(S+1);
     for(int s=0;s<=S;s++){
      int lo=max(-s,s-a-b);
      vals[s]=pref[s+S+1]-pref[lo+S];
      assert(vals[s]>=0);
     }
     out.p.take(vals); out.shapes++;
    }
   }
  }
  int d=++finished;
  if(d%20==0||d==(int)tiles.size()){
   #pragma omp critical
   { stamp(); cout<<" H4-profile mode="<<mode<<" tile_pairs="<<d<<"/"<<tiles.size()<<endl; }
  }
 }
 H4[mode].ok=true;
 uint64_t total=0;
 for(const Tile& t:tiles){
  assert(t.p.ok); total+=t.shapes;
  for(int s=0;s<=S;s++) H4[mode].v[s]=min(H4[mode].v[s],t.p.v[s]);
 }
 for(int s=0;s<=S;s++) H4[mode].v[s]=max(H4[mode].v[s],H[mode].v[s]);
 stamp(); cout<<" H4-profile mode="<<mode<<" exact_shapes="<<total<<" coords";
 for(U v:H4[mode].v) cout<<" "<<v;
 cout<<endl;
 if(mode==0) assert(total==57920058);
 else assert(total==27446154);
}
void make_small_profile(int mode,int ix){
 const vector<int>& a=tuples[ix];
 if(mode && !even_tuple(a)) return;
 int L=mode?14:13,D=accumulate(a.begin(),a.end(),0),M=D+2*S;
 vector<I> f=row(a);
 vector<vector<I>> g(S+1,vector<I>(M+3,0)), pref(S+1,vector<I>(M+3,0));
 for(int s=0;s<=S;s++){
  for(int t=0;t<=D;t++) if(f[t]){
   g[s][abs(t-2*s)]+=f[t];
   g[s][t+2*s+2]-=f[t];
  }
  for(int t=2;t<=M;t++) g[s][t]+=g[s][t-2];
  pref[s]=g[s];
  for(int t=2;t<=M;t++) pref[s][t]+=pref[s][t-2];
 }
 vector<I> vals(S+1);
 for(int s=0;s<=S;s++) vals[s]=(2*s<=D?f[2*s]:0);
 profiles[mode][ix].small[0].take(vals);
 for(int n=L;n<=D;n+=mode?2:1){
  for(int s=0;s<=S;s++) vals[s]=g[s][n];
  profiles[mode][ix].small[1].take(vals);
 }
 for(int u=D%2;u<=D;u+=2){
  int upper=max(L+1,(D+2*S-u)/2+(mode?2:1));
  for(int len=L+1;len<=upper;len+=mode?2:1){
   int top=min(M,u+2*len-2);
   top-=((top-u)%2+2)%2;
   for(int s=0;s<=S;s++) vals[s]=pref[s][top]-(u>=2?pref[s][u-2]:0);
   profiles[mode][ix].small[2].take(vals);
  }
 }
}
void make_triple_profile(int mode,int ix){
 const vector<int>& a=tuples[ix];
 if(mode && !even_tuple(a)) return;
 if(a.size()>2 || (!a.empty()&&a.back()>4)) return;
 int L=mode?14:13,D=accumulate(a.begin(),a.end(),0),M=D+2*S,C=M+(mode?2:0);
 if(mode && (D&1)) return;
 vector<I> f=row(a);
 Prof& out=profiles[mode][ix].triple;
 for(int g=-D;g<=C;g+=2){
  int hs=max(g,2*L-g);
  for(int h=hs;h<=C;h+=mode?4:2){
   for(int k=h;k<=C;k+=mode?4:2){
    vector<I> pref(M+1,0), vals(S+1,0);
    for(int t=D%2;t<=M;t+=2){
     I lo=max({-t,t-h-k,-h,-k});
     pref[t]=max<I>(0,(min<I>(t,g)-lo)/2+1);
    }
    for(int t=2;t<=M;t++) pref[t]+=pref[t-2];
    for(int u=0;u<=D;u++) if(f[u]){
     for(int s=0;s<=S;s++){
      int lo=abs(2*s-u),hi=2*s+u;
      vals[s]+=f[u]*(pref[hi]-(lo>=2?pref[lo-2]:0));
     }
    }
    out.take(vals);
   }
  }
 }
 if(out.ok) for(int s=0;s<=S;s++) out.v[s]=max(out.v[s],H[mode].v[s]);
}

struct Group {
 array<int,4> g{};
 U ways[5][2][2]{};
};
vector<Group> groups;
U choose_n(int n,int k){
 static U C[5][5]={{1}};
 static bool ready=false;
 if(!ready){ for(int i=0;i<=4;i++){ C[i][0]=C[i][i]=1; for(int j=1;j<i;j++) C[i][j]=C[i-1][j-1]+C[i-1][j]; } ready=true; }
 return C[n][k];
}
void make_groups(){
 for(int a=0;a<=4;a++) for(int b=0;b<=4-a;b++) for(int c=0;c<=4-a-b;c++){
  int d=4-a-b-c;
  Group G; G.g={a,b,c,d};
  for(int r=0;r<=4;r++){
   for(int i=0;i<=min(r,a);i++) for(int j=0;j<=min(r-i,b);j++)
    for(int k=0;k<=min(r-i-j,c);k++){
     int l=r-i-j-k; if(l>d) continue;
     G.ways[r][(k+l)%2][(j+l)%2]+=
      choose_n(a,i)*choose_n(b,j)*choose_n(c,k)*choose_n(d,l);
    }
  }
  groups.push_back(G);
 }
 assert(groups.size()==35);
}
struct Cut { int mask,r,op,k; U cost; };
struct Result { U maxval=0; uint64_t checks=0,fails=0; array<int,4> g{}; int signs=0; };
const Prof& profile_for(int mode,const vector<int>& a,int r){
 int ix=index_of.at(key_of(a));
 if(r==3 && profiles[mode][ix].triple.ok) return profiles[mode][ix].triple;
 if(r>=4) return H4[mode];
 if(r==3) return H[mode];
 return profiles[mode][ix].small[r];
}
Result evaluate(int mode,const vector<int>& a){
 const int h=7,big=4,full=(1<<h)-1;
 vector<Cut> cuts;
 for(int mask=0;mask<=full;mask++){
  vector<int> part,other; int weight=0;
  for(int i=0;i<h;i++){
   if(mask>>i&1){part.push_back(a[i]);weight+=a[i];}
   else other.push_back(a[i]);
  }
  int t=__builtin_popcount((unsigned)mask);
  for(int k:{3,4,5}){
   int r=k-t;
   if(r<0||r>big) continue;
   const Prof& p=profile_for(mode,part,r);
   const Prof& q=profile_for(mode,other,big-r);
   if(!p.ok||!q.ok) continue;
   W dot=0;
   for(int s=0;s<=S;s++){
    assert(p.v[s]<(1ULL<<40)&&q.v[s]<(1ULL<<40));
    dot+=(W)p.v[s]*q.v[s];
   }
   assert(dot>0&&dot<(W(1)<<90));
   W numerator=(W(1)<<80);
   U cost=(U)((numerator+dot-1)/dot);
   cuts.push_back({mask,r,weight%2,k,cost});
  }
 }
 sort(cuts.begin(),cuts.end(),[](const Cut& x,const Cut& y){return x.cost>y.cost;});
 array<vector<I>,128> rows;
 for(int mask=0;mask<=full;mask++){
  vector<int> part; for(int i=0;i<h;i++) if(mask>>i&1) part.push_back(a[i]);
  rows[mask]=row(part);
 }
 int units[128][49]; for(auto& v:units) fill(begin(v),end(v),61);
 for(const Cut& c:cuts) if(c.r==0){
  const vector<I>& complement=rows[full^c.mask];
  I mult=rows[c.mask][0];
  for(int j=0;j<h;j++) for(int i=0;i<j;i++) if(a[i]==a[j]){
   const vector<I>& donor=rows[full^(1<<i)^(1<<j)];
   int u=0;
   for(int t=0;t<(int)complement.size();t++) if(complement[t]){
    if(t>=(int)donor.size()||donor[t]==0){ u=61; break; }
    W num=(W)60*(U)mult*(U)complement[t];
    int need=(int)min<W>(61,(num+(U)donor[t]-1)/(U)donor[t]);
    u=max(u,need);
   }
   if(c.k==3 && (c.mask>>i&1) && (c.mask>>j&1)){
    int third=c.mask^(1<<i)^(1<<j);
    assert(third && !(third&(third-1)));
    int n=a[__builtin_ctz((unsigned)third)];
    assert((n&1)==0);
    U hv=H[mode].v[n/2]; assert(hv>0);
    int need=(int)(((W)60*Q+hv-1)/hv);
    u=min(u,need);
   }
   units[c.mask][i*h+j]=u;
  }
 }
 Result ans;
 for(int sg=0;sg<=full;sg++){
  bool compatible=true;
  for(int i=1;i<h;i++) if(a[i]==a[i-1] && (((sg>>i)^(sg>>(i-1)))&1)) compatible=false;
  if(!compatible) continue;
  vector<pair<int,int>> pairs; vector<int> load;
  for(int j=0;j<h;j++) for(int i=0;i<j;i++) if(a[i]==a[j]){ pairs.push_back({i,j}); load.push_back(0); }
  bool paid[128]={};
  for(const Cut& c:cuts) if(c.r==0 && __builtin_parity((unsigned)(sg&c.mask))){
   int best=-1,bestu=61;
   for(int z=0;z<(int)pairs.size();z++){
    auto [i,j]=pairs[z]; int u=units[c.mask][i*h+j];
    if(u<=60-load[z] && (u<bestu || (u==bestu && (best<0||load[z]<load[best])))){
     best=z; bestu=u;
    }
   }
   if(best>=0){ load[best]+=bestu; paid[c.mask]=true; }
  }
  for(const Group& G:groups){
   int odd=G.g[2]+G.g[3];
   for(int n:a) odd+=n&1;
   if(mode ? odd!=0 : ((odd&1)||odd==0||odd==2)) continue;
   if((__builtin_popcount((unsigned)sg)+G.g[1]+G.g[3])&1) continue;
   W val=0;
   for(const Cut& c:cuts){
    if(c.r==0&&paid[c.mask]) continue;
    int mp=1-__builtin_parity((unsigned)(sg&c.mask));
    val+=(W)c.cost*G.ways[c.r][c.op][mp];
   }
   assert(val<(W(1)<<63));
   U x=(U)val; ans.checks++; ans.fails+=(x>BQ);
   if(x>ans.maxval){ ans.maxval=x; ans.g=G.g; ans.signs=sg; }
  }
 }
 return ans;
}
string show(const vector<int>& a){
 string z="("; for(int i=0;i<(int)a.size();i++){ if(i) z+=","; z+=to_string(a[i]); } return z+")";
}
int run(int argc,char**argv){
 if(argc>1 && (string(argv[1])=="-h"||string(argv[1])=="--help")){
  cout<<"Prefix certificate audit. Usage: verifier [--threads N]. Enumerates every sorted 7-tuple in [1,12]^7; exact integer arithmetic.\n"; return 0;
 }
 for(int i=1;i+1<argc;i++) if(string(argv[i])=="--threads") threads=stoi(argv[i+1]);
 if(threads<1||threads>32) return 2;
 omp_set_num_threads(threads);
 for(int n=0;n<=7;n++){ vector<int> a; enumerate(a,n,1); }
 for(int i=0;i<(int)tuples.size();i++) index_of.emplace(key_of(tuples[i]),i);
 assert(tuples.size()==50388&&prefixes.size()==31824);
 stamp(); cout<<"enumerated sorted prefixes="<<prefixes.size()<<" all_subtuples="<<tuples.size()<<endl;
 for(int m=0;m<2;m++){
  profiles[m].resize(tuples.size());
  make_H(m);
 }
 for(int m=0;m<2;m++) make_H4(m);
 struct Job{int mode,ix;};
 vector<Job> jobs;
 for(int ix=0;ix<(int)tuples.size();ix++){
  jobs.push_back({0,ix});
  if(even_tuple(tuples[ix])) jobs.push_back({1,ix});
 }
 atomic<int> lowdone{0};
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)jobs.size();j++){
  make_small_profile(jobs[j].mode,jobs[j].ix);
  int d=++lowdone;
  if(d%5000==0||d==(int)jobs.size()){
   #pragma omp critical
   { stamp(); cout<<" exact 0/1/2-factor profiles "<<d<<"/"<<jobs.size()<<endl; }
  }
 }
 vector<Job> triple_jobs;
 for(int m=0;m<2;m++) for(int ix=0;ix<(int)tuples.size();ix++){
  const auto& a=tuples[ix];
  if(a.size()<=2&&(a.empty()||a.back()<=4)&&(!m||even_tuple(a))) triple_jobs.push_back({m,ix});
 }
 atomic<int> tdone{0};
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)triple_jobs.size();j++){
  make_triple_profile(triple_jobs[j].mode,triple_jobs[j].ix);
  int d=++tdone;
  if(d%25==0||d==(int)triple_jobs.size()){
   #pragma omp critical
   { stamp(); cout<<" exact 3-large profiles "<<d<<"/"<<triple_jobs.size()<<endl; }
  }
 }
 make_groups();
 U worst_good=0,worst_mode[2]={0,0};
 vector<int> worst_good_prefix; int worst_good_sector=0; Result worst_good_result;
 int good=0,bad=0,even7=0;
 vector<pair<vector<int>,array<Result,2>>> failed;
 uint64_t checks[2]={0,0};
 atomic<int> done=0;
 #pragma omp parallel for schedule(dynamic)
 for(int ix=0;ix<(int)prefixes.size();ix++){
  const vector<int>& a=prefixes[ix];
  array<Result,2> r;
  r[0]=evaluate(0,a);
  if(even_tuple(a)) r[1]=evaluate(1,a);
  #pragma omp critical
  {
   checks[0]+=r[0].checks; worst_mode[0]=max(worst_mode[0],r[0].maxval);
   if(even_tuple(a)){ checks[1]+=r[1].checks; worst_mode[1]=max(worst_mode[1],r[1].maxval); even7++; }
   bool pass=r[0].maxval<=BQ && r[0].fails==0;
   if(even_tuple(a)) pass=pass&&r[1].maxval<=BQ&&r[1].fails==0;
   if(pass){
    good++;
    int sector=(even_tuple(a)&&r[1].maxval>r[0].maxval)?1:0;
    U candidate=sector?r[1].maxval:r[0].maxval;
    if(candidate>worst_good){worst_good=candidate;worst_good_prefix=a;worst_good_sector=sector;worst_good_result=r[sector];}
   } else {bad++; failed.push_back({a,r});}
   int d=++done;
   if(d%1000==0||d==(int)prefixes.size()){
    stamp(); cout<<" prefix budgets "<<d<<"/"<<prefixes.size()<<" good="<<good<<" exceptions="<<bad<<endl;
   }
  }
 }
 assert(good==31802&&bad==22&&even7==792);
 sort(failed.begin(),failed.end(),[](const auto& x,const auto& y){return x.first<y.first;});
 vector<tuple<vector<int>,U,U>> expected={
  {{1,1,1,1,1,1,1},1621337291966ULL,0},
  {{1,1,1,1,1,1,2},2332942233138ULL,0},
  {{1,1,1,1,1,1,3},1136325104908ULL,0},
  {{1,1,1,1,1,2,2},1604581599466ULL,0},
  {{1,1,1,1,1,2,3},1498556251543ULL,0},
  {{1,1,1,1,2,2,2},1735982437898ULL,0},
  {{1,1,1,1,2,2,3},1506479539434ULL,0},
  {{1,1,1,1,2,3,3},1295920927919ULL,0},
  {{1,1,1,2,2,2,2},1462173942544ULL,0},
  {{1,1,1,2,2,2,3},1543637323382ULL,0},
  {{1,1,1,2,2,3,3},1410207480631ULL,0},
  {{1,1,1,2,3,3,3},1127374686154ULL,0},
  {{1,1,2,2,2,2,2},1337247769051ULL,0},
  {{1,1,2,2,2,2,3},1399134694503ULL,0},
  {{1,1,2,2,2,3,3},1319693126067ULL,0},
  {{1,1,2,2,3,3,3},1229508916741ULL,0},
  {{1,1,2,2,3,3,4},1105721880588ULL,0},
  {{2,2,2,2,2,2,2},1866632809993ULL,1764101649107ULL},
  {{2,2,2,2,2,2,4},1578983291282ULL,1493883815905ULL},
  {{2,2,2,2,2,4,4},1538330595366ULL,1467747830904ULL},
  {{2,2,2,2,4,4,4},1307175099713ULL,1256830694533ULL},
  {{2,2,2,4,4,4,4},1176008868237ULL,1148460578721ULL}
 };
 assert(failed.size()==expected.size());
 for(int i=0;i<(int)expected.size();i++){
  const auto& [tuple,m0,m1]=expected[i];
  assert(failed[i].first==tuple);
  assert(failed[i].second[0].maxval==m0);
  assert(failed[i].second[1].maxval==m1);
 }
 assert(worst_good==1093158045286ULL);
 assert(worst_good_prefix==vector<int>({1,1,1,2,2,3,4}));
 assert(worst_good_sector==0&&worst_good_result.signs==88);
 assert((worst_good_result.g==array<int,4>{0,0,3,1}));
 stamp(); cout<<"RESULT prefixes="<<prefixes.size()<<" good="<<good<<" exceptions="<<bad
             <<" all_even_prefixes="<<even7<<" exact_charge_bound="<<BQ
             <<" worst_good_numerator="<<worst_good
             <<" worst_mode0_numerator="<<worst_mode[0]
             <<" worst_mode1_numerator="<<worst_mode[1]
             <<" charge_margin="<<(BQ-worst_good)
             <<" worst_prefix="<<show(worst_good_prefix)
             <<" worst_sector="<<worst_good_sector
             <<" worst_group=("<<worst_good_result.g[0]<<","<<worst_good_result.g[1]<<","<<worst_good_result.g[2]<<","<<worst_good_result.g[3]<<")"
             <<" worst_signmask="<<worst_good_result.signs
             <<" feasible_sign_group_cases_mode0="<<checks[0]
             <<" feasible_sign_group_cases_mode1="<<checks[1]<<endl;
 for(const auto& item:failed){
  const auto& a=item.first; const auto& r=item.second;
  int D=accumulate(a.begin(),a.end(),0);
  cout<<"EXCEPTION prefix="<<show(a)<<" D="<<D<<" q_finite_cutoff="<<2*D
      <<" mode0_max="<<r[0].maxval<<" mode0_failing_cases="<<r[0].fails
      <<" mode0_witness_group=("<<r[0].g[0]<<","<<r[0].g[1]<<","<<r[0].g[2]<<","<<r[0].g[3]<<")"
      <<" mode0_signmask="<<r[0].signs;
  if(even_tuple(a)) cout<<" mode1_max="<<r[1].maxval<<" mode1_failing_cases="<<r[1].fails
      <<" mode1_witness_group=("<<r[1].g[0]<<","<<r[1].g[1]<<","<<r[1].g[2]<<","<<r[1].g[3]<<")"
      <<" mode1_signmask="<<r[1].signs;
  cout<<endl;
 }
 stamp(); cout<<"CHARGE PREFIX SCREEN PASS"<<endl;
 return 0;
}
}
int main(int argc,char**argv){ return verify::run(argc,argv); }
'''

obj=os.memfd_create("fmchk109b_object",0)
exe=os.memfd_create("fmchk109b_executable",0)
subprocess.run(
    ["g++","-O3","-std=c++17","-fopenmp","-pipe","-x","c++","-",
     "-c","-o",f"/proc/self/fd/{obj}"],
    input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
    ["g++","-###","-fno-use-linker-plugin","-fopenmp",
     f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,check=True)
link=next(shlex.split(line) for line in p.stderr.splitlines()
          if "/collect2 " in line)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
for flag in ("-h","--help"):
    subprocess.run([f"/proc/self/fd/{exe}",flag],pass_fds=(exe,),check=True)
subprocess.run([f"/proc/self/fd/{exe}","--threads","24"],
               pass_fds=(exe,),check=True)
