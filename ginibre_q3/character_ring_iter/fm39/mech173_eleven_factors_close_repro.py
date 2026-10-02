CPP = r'''
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
namespace budget{
using Z=unsigned long long;using I=long long;using W=__uint128_t;
constexpr int K=12,S=2*K+2,BASE=K+1;
constexpr Z Q=1ULL<<20,BQ=1ULL<<40;
int threads=8;
void stamp(){time_t now=time(nullptr);tm v;gmtime_r(&now,&v);char b[32];strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&v);cout<<b;}
struct Prof{bool valid=false;array<Z,S+1>v;Prof(){v.fill(ULLONG_MAX);}
 void add(const vector<I>&b){if(b[0]<=0)return;valid=true;for(int s=0;s<=S;s++){assert(b[s]>=0);v[s]=min(v[s],Z(W(Q)*b[s]/b[0]));}}
};
unordered_map<int,int> keys;
vector<vector<int>>smalls;Prof H[2],H2[6];vector<array<Prof,4>> ps[2];
int code(const vector<int>&sm){int q=0;for(int n:sm)q=BASE*q+n;return keys.at(q);}
void rec(vector<int>sm,int h,int lo){if(!h){smalls.push_back(sm);return;}for(int n=lo;n<=K;n++){sm.push_back(n);rec(sm,h-1,n);sm.pop_back();}}
vector<I>row(const vector<int>&sm){vector<I>f(1,1);for(int n:sm){vector<I>g(f.size()+n);for(int j=0;j<(int)f.size();j++)if(f[j])for(int t=abs(j-n);t<=j+n;t+=2)g[t]+=f[j];f.swap(g);}return f;}
I tri(int n){return n>=0?I(n+1)*(n+2)/2:0;}
I rect(int n,int l,int h){return tri(n)-tri(n-l)-tri(n-h)+tri(n-l-h);}
I band(int u,int v,int l,int h,int t){return rect((v-u+t)/2+h-1,l,h)-rect((v-u-t)/2+h-2,l,h)-rect((t-u-v)/2-1,l,h);}
void hprofile(int mode){
 int L=mode?K+2:K+1;vector<int>ls,hs;
 for(int l=L+1;l<(mode?3:2)*(L+1);l+=(mode?2:1))ls.push_back(l);
 for(int h=1;h<(mode?3:2)*(L+1);h+=(mode?2:1))hs.push_back(h);
 I count=0;
 for(int l:ls)for(int h:hs)for(int off=1-h;off<l;off++)for(int par=0;par<(mode?1:2);par++){
  int u=max({0,L-h+1-2*off,-2*off});u+=(par-u%2+2)%2;int v=u+2*off;
  vector<I>b(S+1);for(int s=0;s<=S;s++)b[s]=band(u,v,l,h,2*s);
  if(b[0]){H[mode].add(b);count++;}
 }
 assert(H[mode].valid);stamp();cout<<" H "<<L<<" "<<count;
 for(auto x:H[mode].v)cout<<" "<<x;cout<<endl;
}
Prof H4[2];
I cb(int a,int b,int c,int l,int h){return rect(a+h-1,l,h)-rect(-b+h-2,l,h)-rect(c-1,l,h);}
void fourprofile(int mode){
 int L=mode?K+2:K+1;vector<int>ls;for(int l=L+1;l<(mode?3:2)*(L+1);l+=(mode?2:1))ls.push_back(l);
 vector<Prof> out(ls.size()*ls.size());vector<I> counts(out.size());atomic<int>done{0};
 #pragma omp parallel for schedule(dynamic)
 for(int job=0;job<(int)out.size();job++){
  int l=ls[job/ls.size()],h=ls[job%ls.size()],R=max({S,l-1,h-1});
  for(int a=1-h;a<=l-1+2*S;a++)for(int b=max(1-l,-a);b<=h-1+2*S;b++){
   if(mode&&(a-b)%2)continue;
   for(int c=-R;c<=min({a,b,l+h-2});c++){
    if(mode&&(a-c)%2)continue;
    I den=cb(a,b,c,l,h);if(!den)continue;assert(den>0);
    I pref[2*S+2]={};for(int z=-S;z<=S;z++){I v=a+b+2*z<0?0:cb(a+z,b+z,c+z,l,h);assert(v>=0);pref[z+S+1]=pref[z+S]+v;}
    vector<I>val(S+1);for(int s=0;s<=S;s++){int lo=max(-s,s-a-b);val[s]=pref[s+S+1]-pref[lo+S];}
    out[job].add(val);counts[job]++;
   }
  }
  int d=++done;if(d%20==0||d==(int)out.size()){
   #pragma omp critical
   {stamp();cout<<" four-profile "<<mode<<" tiles "<<d<<"/"<<out.size()<<endl;}
  }
 }
 H4[mode].valid=true;I total=0;
 for(int j=0;j<(int)out.size();j++){assert(out[j].valid);total+=counts[j];for(int s=0;s<=S;s++)H4[mode].v[s]=min(H4[mode].v[s],out[j].v[s]);}
 for(int s=0;s<=S;s++)H4[mode].v[s]=max(H4[mode].v[s],H[mode].v[s]);
 stamp();cout<<" H4 "<<mode<<" shapes "<<total;for(auto v:H4[mode].v)cout<<" "<<v;cout<<endl;
}
void smallprofile(int mode,const vector<int>&sm){
 int L=mode?K+2:K+1,D=accumulate(sm.begin(),sm.end(),0),M=D+2*S,key=code(sm);auto f=row(sm);
 vector<vector<I>>g(S+1,vector<I>(M+3)),pref=g;
 for(int s=0;s<=S;s++){
  for(int t=0;t<=D;t++)if(f[t]){g[s][abs(t-2*s)]+=f[t];g[s][t+2*s+2]-=f[t];}
  for(int t=2;t<=M;t++)g[s][t]+=g[s][t-2];
  pref[s]=g[s];for(int t=2;t<=M;t++)pref[s][t]+=pref[s][t-2];
 }
 vector<I>b(S+1);for(int s=0;s<=S;s++)b[s]=2*s<=D?f[2*s]:0;ps[mode][key][0].add(b);
 for(int n=L;n<=D;n+=(mode?2:1)){for(int s=0;s<=S;s++)b[s]=g[s][n];ps[mode][key][1].add(b);}
 for(int u=D%2;u<=D;u+=2){
  int hi=max(L+1,(D+2*S-u)/2+(mode?2:1));
  for(int l=L+1;l<=hi;l+=(mode?2:1)){
   int top=min(M,u+2*l-2);top-=((top-u)%2+2)%2;
   for(int s=0;s<=S;s++)b[s]=pref[s][top]-(u>=2?pref[s][u-2]:0);
   ps[mode][key][2].add(b);
  }
 }
}
void twosprofile(){
 for(int q=1;q<=5;q++){
  vector<int>sm(q,2);auto f=row(sm);int C=S+q+1,M=2*S+2*q;I cases=0;
  for(int g=-q;g<=C;g++)for(int h=max(g,K+2-g);h<=C;h++)for(int k=h;k<=C;k++){
   if((g-h)%2||(h-k)%2)continue;
   vector<I>a(M+1),b(S+1);
   for(int u=0;2*u<=M;u++)a[2*u]=max(0,min(u,g)-max({-u,u-h-k,-h,-k})+1);
   for(int j=2;j<=M;j++)a[j]+=a[j-2];
   for(int s=0;s<=S;s++)for(int t=0;t<(int)f.size();t++)if(f[t]){
    int lo=abs(2*s-t),hi=2*s+t;b[s]+=f[t]*(a[hi]-(lo>=2?a[lo-2]:0));}
   if(b[0]){H2[q].add(b);cases++;}
  }
  assert(H2[q].valid);
  for(int s=0;s<=S;s++)assert(H2[q].v[s]>=H[1].v[s]);
  stamp();cout<<" H2 "<<q<<" "<<cases;for(auto v:H2[q].v)cout<<" "<<v;cout<<endl;
 }
}
void threeprofile(int mode,const vector<int>&sm){
 int L=mode?K+2:K+1,D=accumulate(sm.begin(),sm.end(),0),M=D+2*S,C=M+(mode?2:0);
 if(D%2&&mode)return;auto f=row(sm);Prof &P=ps[mode][code(sm)][3];
 for(int g=-D;g<=C;g+=2)for(int h=max(g,2*L-g);h<=C;h+=(mode?4:2))
 for(int k=h;k<=C;k+=(mode?4:2)){
  vector<I>a(M+1),b(S+1);
  for(int t=D%2;t<=M;t+=2)a[t]=max(0,(min(t,g)-max({-t,t-h-k,-h,-k}))/2+1);
  for(int t=2;t<=M;t++)a[t]+=a[t-2];
  for(int u=0;u<=D;u++)if(f[u])for(int j=0;j<=S;j++){
   int lo=abs(2*j-u),hi=2*j+u;
   b[j]+=f[u]*(a[hi]-(lo>=2?a[lo-2]:0));
  }
  P.add(b);
 }
 if(P.valid)for(int j=0;j<=S;j++)P.v[j]=max(P.v[j],H[mode].v[j]);
}
void coefficient_bridges(){
 mt19937 rng(173);I checks=0;
 for(int mode=0;mode<2;mode++)for(int z=0;z<240;z++){
  int L=mode?K+2:K+1;vector<int>a;
  for(int j=0;j<4;j++)a.push_back(L+(mode?2*(rng()%10):rng()%20));
  for(int j=0;j<z%4;j++)a.push_back(1+rng()%9);
  auto f=row(a);
  for(int s=0;s<=S;s++){
   vector<I> g(f.size()+2*s+2);
   for(int t=0;t<(int)f.size();t++)if(f[t]){g[abs(t-2*s)]+=f[t];g[t+2*s+2]-=f[t];}
   for(int t=2;t<(int)g.size();t++)g[t]+=g[t-2];
   for(int t=0;t<(int)f.size();t++){assert(W(Q)*g[t]>=W(H4[mode].v[s])*f[t]);checks++;}
  }
 }
 for(int mode=0;mode<2;mode++)for(auto sm:smalls)if(sm.size()<=5&&(sm.empty()||sm.back()<=4)&&(!mode||all_of(sm.begin(),sm.end(),[](int x){return x%2==0;}))){
  auto P=ps[mode][code(sm)][3];if(!P.valid)continue;
  for(int trial=0;trial<4;trial++){
   int L=mode?K+2:K+1;vector<int>a=sm;
   for(int j=0;j<3;j++)a.push_back(L+(mode?2*(rng()%10):rng()%20));
   if(accumulate(a.begin(),a.end(),0)%2)a.back()++;
   auto f=row(a);
   for(int s=0;s<=S;s++){I v=2*s<(int)f.size()?f[2*s]:0;assert(W(Q)*v>=W(P.v[s])*f[0]);checks++;}
  }
 }
 stamp();cout<<" coefficient bridges "<<checks<<" PASS"<<endl;
}
struct Group{array<int,4>g;I ways[6][2][2]{};};
vector<Group>groups[12];I choose[12][12];
void initgroups(){
 for(int n=0;n<=11;n++){choose[n][0]=choose[n][n]=1;for(int k=1;k<n;k++)choose[n][k]=choose[n-1][k-1]+choose[n-1][k];}
 for(int b=4;b<=11;b++)for(int a=0;a<=b;a++)for(int c=0;c<=b-a;c++)for(int d=0;d<=b-a-c;d++){
  Group G;G.g={a,c,d,b-a-c-d};
  for(int r=0;r<=5;r++)for(int i=0;i<=min(r,G.g[0]);i++)for(int j=0;j<=min(r-i,G.g[1]);j++)
  for(int k=0;k<=min(r-i-j,G.g[2]);k++){
   int l=r-i-j-k;if(l>G.g[3])continue;
   G.ways[r][(k+l)%2][(j+l)%2]+=choose[G.g[0]][i]*choose[G.g[1]][j]*choose[G.g[2]][k]*choose[G.g[3]][l];
  }groups[b].push_back(G);
 }
}
struct Cut{int mask,r,op,k;Z cost;};
struct Answer{Z max=0;I checks=0,fails=0;vector<int>sm;array<int,4>g;int sg=0;};
Answer budget(int mode,const vector<int>&sm){
 int h=sm.size(),big=11-h,full=(1<<h)-1;vector<Cut>cuts;
 auto profile=[&](int mask,int r)->const Prof&{
  vector<int>a;for(int j=0;j<h;j++)if(mask>>j&1)a.push_back(sm[j]);
  if(r==3&&ps[mode][code(a)][3].valid)return ps[mode][code(a)][3];
  if(mode&&r==3&&!a.empty()&&all_of(a.begin(),a.end(),[](int n){return n==2;}))return H2[a.size()];
  if(r>=4)return H4[mode];
  if(r>=3)return H[mode];
  return ps[mode][code(a)][r];};
 for(int mask=0;mask<=full;mask++){
  int t=__builtin_popcount((unsigned)mask),w=0;for(int j=0;j<h;j++)if(mask>>j&1)w+=sm[j];
  for(int k:{3,4,5}){
   int r=k-t;if(r<0||r>big)continue;
   auto&p=profile(mask,r);auto&q=profile(full^mask,big-r);if(!p.valid||!q.valid)continue;
   W dot=0;for(int s=0;s<=S;s++){assert(p.v[s]<(1ULL<<40)&&q.v[s]<(1ULL<<40));dot+=W(p.v[s])*q.v[s];}
   assert(dot>0&&dot<(W(1)<<90));
   W num=W(BQ)*Q*Q;Z c=Z((num+dot-1)/dot);cuts.push_back({mask,r,w%2,k,c});
  }
 }
 sort(cuts.begin(),cuts.end(),[](auto a,auto b){return a.cost>b.cost;});
 vector<vector<I>> rr(full+1);
 for(int mask=0;mask<=full;mask++){vector<int>a;for(int j=0;j<h;j++)if(mask>>j&1)a.push_back(sm[j]);rr[mask]=row(a);}
 int units[128][49];fill(&units[0][0],&units[0][0]+128*49,61);
 for(auto c:cuts)if(c.r==0){
  auto &smallrow=rr[full^c.mask];I mul=rr[c.mask][0];
  for(int j=0;j<h;j++)for(int i=0;i<j;i++)if(sm[i]==sm[j]){
   auto &donor=rr[full^(1<<i)^(1<<j)];int unit=0;
   for(int t=0;t<(int)smallrow.size();t++)if(smallrow[t]){
    if(t>=(int)donor.size()||!donor[t]){unit=61;break;}
    unit=max(unit,int(min<I>(61,(60*mul*smallrow[t]+donor[t]-1)/donor[t])));
   }
   if(c.k==3&&(c.mask>>i&1)&&(c.mask>>j&1)){
    int third=c.mask^(1<<i)^(1<<j),t=sm[__builtin_ctz((unsigned)third)];
    assert(t%2==0);unit=min(unit,int((60*Q+H[mode].v[t/2]-1)/H[mode].v[t/2]));
   }
   units[c.mask][i*h+j]=unit;
  }
 }
 Answer out;out.sm=sm;
 for(int sg=0;sg<=full;sg++){
  bool okay=true;for(int i=1;i<h;i++)if(sm[i]==sm[i-1]&&((sg>>i^sg>>(i-1))&1))okay=false;if(!okay)continue;
  vector<pair<int,int>>pairs;vector<int>load;for(int j=0;j<h;j++)for(int i=0;i<j;i++)if(sm[i]==sm[j]){pairs.push_back({i,j});load.push_back(0);}
  bool paid[128]={};
  for(auto c:cuts)if(c.r==0&&__builtin_parity((unsigned)(sg&c.mask))){
   int best=-1,bestunit=61;
   for(int z=0;z<(int)pairs.size();z++){
    auto[i,j]=pairs[z];int unit=units[c.mask][i*h+j];
    if(unit<=60-load[z]&&(unit<bestunit||(unit==bestunit&&(best<0||load[z]<load[best])))){best=z;bestunit=unit;}
   }
   if(best>=0){load[best]+=bestunit;paid[c.mask]=true;}
  }
  for(auto&G:groups[big]){
   int odd=G.g[2]+G.g[3];for(int n:sm)odd+=n%2;
   if(mode?odd!=0:(odd%2||odd==0||odd==2))continue;
   if((__builtin_popcount((unsigned)sg)+G.g[1]+G.g[3])%2)continue;
   Z val=0;
   for(auto c:cuts){if(c.r==0&&paid[c.mask])continue;
    int mp=1-__builtin_parity((unsigned)(sg&c.mask));
    val+=c.cost*G.ways[c.r][c.op][mp];
   }
   assert(val<(1ULL<<55));out.checks++;out.fails+=val>BQ;
   if(val>out.max){out.max=val;out.g=G.g;out.sg=sg;}
  }
 }
 return out;
}
vector<bool> safe;
int run(int argc,char**argv){
 if(argc>1&&(string(argv[1])=="-h"||string(argv[1])=="--help")){cout<<"FM173 exact dyadic cut budgets, cutoff 8; memory only.\n";return 0;}
 for(int i=1;i+1<argc;i++)if(string(argv[i])=="--threads")threads=stoi(argv[i+1]);
 assert(threads>=1&&threads<=32);omp_set_num_threads(threads);for(int h=0;h<=7;h++)rec({},h,1);
 for(int j=0;j<(int)smalls.size();j++){int q=0;for(int n:smalls[j])q=BASE*q+n;keys.emplace(q,j);}
 for(int m=0;m<2;m++){ps[m].resize(smalls.size());hprofile(m);for(int s=1;s<=K/2;s++)assert(H[m].v[s]>=2*Q);}
 fourprofile(0);fourprofile(1);twosprofile();
 struct Job{int mode,idx;};vector<Job>jobs;
 for(int m=0;m<2;m++)for(int j=0;j<(int)smalls.size();j++)
  if(!m||all_of(smalls[j].begin(),smalls[j].end(),[](int n){return n%2==0;}))jobs.push_back({m,j});
 atomic<int>done{0};
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)jobs.size();j++){auto a=jobs[j];smallprofile(a.mode,smalls[a.idx]);int d=++done;
  if(d%5000==0){
   #pragma omp critical
   {stamp();cout<<" small profiles "<<d<<"/"<<jobs.size()<<endl;}
  }}
 stamp();cout<<" small profiles "<<jobs.size()<<"/"<<jobs.size()<<" PASS"<<endl;
 vector<Job> ref;
 for(auto a:jobs)if(smalls[a.idx].size()<=5 && (smalls[a.idx].empty()||smalls[a.idx].back()<=4))ref.push_back(a);
 done=0;
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)ref.size();j++){auto a=ref[j];threeprofile(a.mode,smalls[a.idx]);int d=++done;
  if(d%20==0||d==(int)ref.size()){
   #pragma omp critical
   {stamp();cout<<" exact three profiles "<<d<<"/"<<ref.size()<<endl;}
  }}
 coefficient_bridges();
 initgroups();vector<Answer>results(jobs.size());done=0;
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<(int)jobs.size();j++){auto a=jobs[j];results[j]=budget(a.mode,smalls[a.idx]);int d=++done;
  if(d%5000==0){
   #pragma omp critical
   {stamp();cout<<" budgets "<<d<<"/"<<jobs.size()<<endl;}
  }}
 safe.assign(smalls.size(),true);for(int j=0;j<(int)jobs.size();j++)if(results[j].fails)safe[jobs[j].idx]=false;
 I good=0,bad=0;for(int j=0;j<(int)smalls.size();j++)if(smalls[j].size()==7){if(safe[j])good++;else bad++;}
 stamp();cout<<" PREFIX-CERT good "<<good<<" bad "<<bad<<endl;assert(good==31802&&bad==22);
 for(int m=0;m<2;m++)for(int h=0;h<=7;h++){
  Answer total;for(int j=0;j<(int)jobs.size();j++)if(jobs[j].mode==m&&(int)smalls[jobs[j].idx].size()<=h){
   auto a=results[j];total.checks+=a.checks;total.fails+=a.fails;if(a.max>total.max){I c=total.checks,f=total.fails;total=a;total.checks=c;total.fails=f;}}
  stamp();cout<<" BUDGET "<<m<<" "<<h<<" "<<total.checks<<" "<<total.fails<<" "<<total.max<<"/"<<BQ<<" witness";
  if(h<=6)assert(total.fails==0);
  if(h==6){assert(total.checks==(m?29676:5323784));assert(total.max==(m?1017343725513ULL:1048748276582ULL));}
  for(int n:total.sm)cout<<" "<<n;cout<<" groups";for(int n:total.g)cout<<" "<<n;cout<<" signs "<<total.sg<<endl;
 }
 return 0;
}
}
namespace growth{
using Z=long long;using W=__int128_t;
void stamp(){time_t now=time(nullptr);tm v;gmtime_r(&now,&v);char b[32];strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&v);cout<<b;}
struct A{int D,v,core;Z count=0,sign=0,even=0;};
map<tuple<int,int,int>,A> mp;int K,L;
void rec(vector<int>&b,int todo,int lo){
 if(todo){for(int n=lo;n<=K;n++){b.push_back(n);rec(b,todo-1,n);b.pop_back();}return;}
 int D=accumulate(b.begin(),b.end(),0),v=b.back(),core=0,ng=0,last=0,mul=0;bool even=true;
 for(int n:b){core+=n>=3;if(n!=last){even&=mul%2==0;ng++;last=n;mul=0;}mul++;}even&=mul%2==0;core=min(core,2);
 auto &a=mp[{D,v,core}];a.D=D;a.v=v;a.core=core;a.count++;a.sign+=1LL<<(ng-1);if(even)a.even+=1LL<<(ng-1);
}
void run(int N,int k,int threads){
 L=N;K=k;mp.clear();vector<int>b;rec(b,L-4,1);vector<A> jobs;Z pre=0;
 for(auto [key,a]:mp){jobs.push_back(a);pre+=a.count;}
 W bound=W(pre)*2*(L-4)*K;for(int i=0;i<3;i++)bound*=((L-4)*K+1);bound*=1LL<<L;assert(bound<(W(1)<<63));
 Z words=0,signs=0;atomic<int>done{0};double start=omp_get_wtime();
 #pragma omp parallel for schedule(dynamic) reduction(+:words,signs)
 for(int j=0;j<(int)jobs.size();j++){
  auto a=jobs[j];int D=a.D;
  for(int q=a.v;q<=2*D;q++)for(int r=q;r<=q+D;r++)for(int s=r;s<=q+D;s++){
   if(a.core+(q>=3)+(r>=3)+(s>=3)<2)continue;
   int lo=max({s,2*q-D,6}),hi=min({q+D,D+q+r-s,D+q+r+s-16});
   lo+=(D+q+r+s-lo)&1;if(lo>hi)continue;
   Z w=(hi-lo)/2+1;int ne=(q>a.v)+(r>q)+(s>r);
   words+=a.count*w;signs+=a.sign*(1LL<<ne)*(2*w-(lo==s));
   if(lo==s&&q==r)signs+=a.even*(1LL<<ne);
  }
  int d=++done;if(d%300==0||d==(int)jobs.size()){
   #pragma omp critical
   {stamp();cout<<" count "<<L<<" cutoff "<<K<<" groups "<<d<<"/"<<jobs.size()<<endl;}
  }
 }
 stamp();cout<<" SIZE "<<L<<" "<<K<<" "<<words<<" "<<signs<<" prefixes "<<pre<<" elapsed "<<omp_get_wtime()-start<<endl;
 if(L==11&&K==8)assert(words==358095523&&signs==44895327698LL);
 if(L==12&&K==8)assert(words==1108557860&&signs==173734117438LL);
}
}
namespace census{
using Z=long long;using Wide=__int128_t;
int K=12,threads=24,sample_count=0;int pc[2048];vector<vector<int>>prefixes;
void stamp(){time_t now=time(nullptr);tm v;gmtime_r(&now,&v);char b[32];strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&v);cout<<b;}
void rec(vector<int>a,int h,int lo){if(!h){prefixes.push_back(a);return;}for(int n=lo;n<=K;n++){a.push_back(n);rec(a,h-1,n);a.pop_back();}}
vector<Z>fusion(const vector<int>&ns){
 vector<Z>v(1,1);
 for(int n:ns){vector<Z>w(v.size()+n+2);
  for(int t=0;t<(int)v.size();t++)if(v[t]){w[abs(t-n)]+=v[t];w[t+n+2]-=v[t];}
  for(int t=2;t<(int)w.size();t++)w[t]+=w[t-2];
  w.resize(v.size()+n);v.swap(w);
 }return v;
}
struct Small{int w=0,par=0,jmax=0;vector<Z>f,p0,p1,p2;};
struct Solver{
 array<int,11>a;array<Small,128>sm;Z m[2048],f[2048];int can[128],ways[128]={};vector<int>uniq;
 Solver(const vector<int>&b){
  for(int i=0;i<7;i++)a[i]=b[i];
  for(int s=0;s<128;s++){
   int cm=0;for(int j=0;j<7;){int k=j;while(k<7&&a[k]==a[j])k++;int ct=0;for(int v=j;v<k;v++)ct+=s>>v&1;cm|=((1<<ct)-1)<<j;j=k;}
   can[s]=cm;ways[cm]++;if(cm==s)uniq.push_back(s);
  }
  for(int s:uniq){
   vector<int>ns;for(int j=0;j<7;j++)if(s>>j&1)ns.push_back(a[j]);
   auto&v=sm[s];v.w=accumulate(ns.begin(),ns.end(),0);v.par=v.w%2;v.jmax=(v.w-v.par)/2;v.f=fusion(ns);
   v.p0.resize(v.jmax+1);v.p1.resize(v.jmax+1);v.p2.resize(v.jmax+1);
   for(int j=0;j<=v.jmax;j++){Z z=v.f[v.par+2*j];v.p0[j]=z+(j?v.p0[j-1]:0);v.p1[j]=z*j+(j?v.p1[j-1]:0);v.p2[j]=z*j*j+(j?v.p2[j-1]:0);}
  }
 }
 Z pref(int s,int x){auto&v=sm[s];return x<0?0:v.p0[min(x,v.jmax)];}
 void mults(){
  int bw[16]={},bp[16]={},bm[16]={};
  for(int t=1;t<16;t++){int j=__builtin_ctz((unsigned)t),v=t&(t-1);bw[t]=bw[v]+a[7+j];bp[t]=bp[v]+1;bm[t]=max(bm[v],a[7+j]);}
  for(int t=0;t<16;t++){
   int r=bp[t],par=bw[t]%2;
   pair<int,int>terms[16];int nt=0;
   if(r>=3)for(int u=t;;u=(u-1)&t){int x=(bw[t]-par)/2-bw[u]-bp[u];if(x>=0)terms[nt++]={x,bp[u]%2?-1:1};if(!u)break;}
   for(int s:uniq){
    auto&v=sm[s];int idx=s+(t<<7);m[idx]=0;
    if(v.par!=par||pc[idx]==1)continue;
    if(r&&2*bm[t]>v.w+bw[t])continue;
    if(!r){m[idx]=v.f[0];continue;}
    if(r==1){m[idx]=bw[t]<=v.w?v.f[bw[t]]:0;continue;}
    if(r==2){
     int x=__builtin_ctz((unsigned)t),y=__builtin_ctz((unsigned)(t&(t-1)));
     int lo=abs(a[7+x]-a[7+y]),hi=a[7+x]+a[7+y];
     m[idx]=pref(s,(hi-par)/2)-pref(s,(lo-par)/2-1);continue;
    }
    Z out=0;
    for(int ti=0;ti<nt;ti++){auto[x,sg]=terms[ti];
     int j=min(x,v.jmax);Z val;
     if(r==3)val=Z(x+1)*v.p0[j]-v.p1[j];
     else{
      Wide z=Wide(x+1)*(x+2)*v.p0[j]-Wide(2*x+3)*v.p1[j]+v.p2[j];
      assert(z>=0&&z%2==0&&z<(Wide(1)<<60));val=Z(z/2);
     }
     out+=sg*val;
    }
    assert(out>=0&&out<(1LL<<60));m[idx]=out;
   }
  }
 }
 Z value(int mask){return m[(mask&~127)+can[mask&127]];}
 Z signcount(){
  int ng=0,last=0,mul=0;bool odd=false;
  for(int x:a){if(x!=last){odd|=mul%2;ng++;last=x;mul=0;}mul++;}
  odd|=mul%2;return 1LL<<(ng-(odd?1:0));
 }
 bool check(Z&signs,Z&low){
  mults();Z mid=0;
  for(int t=0;t<16;t++)for(int smask:uniq){int mask=smask+(t<<7);if(pc[mask]>=3&&pc[mask]<=5){
   Wide z=Wide(value(mask))*value(2047^mask)*ways[smask];assert(z>=0&&z<(Wide(1)<<60));mid+=Z(z);
  }}
  assert(mid>=0&&mid<(1LL<<60));signs+=signcount();
  if(value(2047)>=mid){low=min(low,2*(value(2047)-mid));return true;}
  fill(f,f+2048,0);f[0]=value(2047);
  for(int mask=1;mask<2047;mask++)if(pc[mask]>=2&&pc[mask]<=5){
   Wide z=Wide(value(mask))*value(2047^mask);assert(z>=0&&z<(Wide(1)<<60));f[mask]=Z(z);
  }
  for(int h=1;h<2048;h*=2)for(int st=0;st<2048;st+=2*h)for(int j=st;j<st+h;j++){
   Z x=f[j],y=f[j+h];assert(Wide(abs(x))+abs(y)<(Wide(1)<<61));f[j]=x+y;f[j+h]=x-y;
  }
  int gr[11]={},ng=0;
  for(int j=0;j<11;j++){if(!j||a[j]!=a[j-1])ng++;gr[ng-1]|=1<<j;}
  for(int sg=0;sg<(1<<ng);sg++){
   int mask=0;for(int j=0;j<ng;j++)if(sg>>j&1)mask|=gr[j];if(pc[mask]%2)continue;
   if(f[mask]<0){
    #pragma omp critical
    {stamp();cout<<" FAILURE";for(int j=0;j<11;j++)cout<<" "<<(mask>>j&1?-a[j]:a[j]);cout<<" phi "<<2*f[mask]<<endl;}
    abort();
   }low=min(low,2*f[mask]);
  }return false;
 }
};
void bridges(){
 mt19937 rng(169);Z checks=0;
 for(int j=0;j<320;j++){
  vector<int>b(7);for(int&n:b)n=1+rng()%K;sort(b.begin(),b.end());Solver s(b);
  int D=accumulate(b.begin(),b.end(),0);for(int k=7;k<11;k++)s.a[k]=b.back()+rng()%(3*D-b.back()+1);sort(s.a.begin()+7,s.a.end());
  s.mults();
  for(int k=0;k<24;k++){
   int mask=rng()%2048;vector<int>ns;for(int q=0;q<11;q++)if(mask>>q&1)ns.push_back(s.a[q]);
   assert(s.value(mask)==fusion(ns)[0]);checks++;
  }
 }
 stamp();cout<<" moment/fusion bridges "<<checks<<" PASS"<<endl;
}
int run(int argc,char**argv){
 for(int i=1;i<argc;i++){
  string s=argv[i];if(s=="-h"||s=="--help"){cout<<"FM173 finite eleven-factor box: --cutoff 12 --threads 24 [--sample N]. Memory only.\n";return 0;}
  if((s=="--cutoff"||s=="--threads"||s=="--sample")&&i+1<argc){
   int x=stoi(argv[++i]);if(s=="--cutoff")K=x;else if(s=="--threads")threads=x;else sample_count=x;
  }else return 2;
 }
 assert(1<=K&&K==12&&threads>=1&&threads<=32);
 Wide bound=1;for(int i=0;i<7;i++)bound*=K+1;for(int i=0;i<2;i++)bound*=21*K+1;
 assert(2048*bound<(Wide(1)<<60));
 for(int s=0;s<2048;s++)pc[s]=__builtin_popcount((unsigned)s);
 bridges();rec({},7,1);
 vector<int>jobs;
 if(sample_count){for(int j=0;j<min(sample_count,(int)prefixes.size());j++)jobs.push_back(j*((int)prefixes.size()-1)/max(1,sample_count-1));}
 else{jobs.resize(prefixes.size());iota(jobs.begin(),jobs.end(),0);}
 atomic<Z>words{0},signs{0},fast{0};atomic<int>done{0};Z lower=LLONG_MAX;omp_set_num_threads(threads);double started=omp_get_wtime();
 #pragma omp parallel for schedule(dynamic)
 for(int it=0;it<(int)jobs.size();it++){
  int job=jobs[it];auto b=prefixes[job];Solver v(b);int D=accumulate(b.begin(),b.end(),0),base=0;for(int n:b)base+=n>=3;
  Z w=0,sg=0,ff=0,low=LLONG_MAX;
  int qhi=budget::safe[budget::code(b)]?min(K,2*D):2*D;
  for(int q=b.back();q<=qhi;q++)for(int r=q;r<=q+D;r++)for(int s=r;s<=q+D;s++){
   int p0=max({s,2*q-D,6}),hi=min(q+D,D+q+r-s);p0+=(D+q+r+s-p0)&1;
   for(int p=p0;p<=hi;p+=2){
    int delta=(D+q+r+s-p)/2;if(delta<8||s>delta||base+(q>=3)+(r>=3)+(s>=3)<2)continue;
    v.a[7]=q;v.a[8]=r;v.a[9]=s;v.a[10]=p;w++;ff+=v.check(sg,low);
   }
  }
  words+=w;signs+=sg;fast+=ff;int d=++done;
  #pragma omp critical
  {lower=min(lower,low);
   if(d%32==0||d==(int)jobs.size()){stamp();cout<<" finite prefixes "<<d<<"/"<<jobs.size()<<" words "<<words.load()<<" signs "<<signs.load()<<" bulk "<<fast.load()<<" elapsed "<<fixed<<setprecision(1)<<omp_get_wtime()-started<<endl;}
  }
 }
 stamp();cout<<" FINITE "<<K<<" "<<words.load()<<" "<<signs.load()<<" bulk "<<fast.load()<<" lower "<<lower<<" elapsed "<<omp_get_wtime()-started<<" PASS"<<endl;
 if(!sample_count)assert(words==326484586&&signs==47094362366LL);
 return 0;
}
}
int main(int argc,char**argv){
 int th=24;for(int i=1;i<argc;i++){string s=argv[i];
  if(s=="-h"||s=="--help"){cout<<"FM-MECH173 exact eleven-factor verifier; --threads N (default 24), --sample N for timing only. Memory only.\n";return 0;}
  if(s=="--threads"&&i+1<argc)th=stoi(argv[++i]);
  else if(s=="--sample"&&i+1<argc)i++;
  else return 2;
 }
 assert(th>=1&&th<=32);
 budget::run(argc,argv);
 for(int k:{8,10,12})for(int n=11;n<=13;n++)growth::run(n,k,th);
 census::run(argc,argv);
 census::stamp();cout<<" FM-MECH173 PASS"<<endl;
}
'''
import os, shlex, subprocess, sys
obj=os.memfd_create("fm173_obj",0)
exe=os.memfd_create("fm173_exe",0)
subprocess.run(
    ["g++","-Werror=return-type","-O3","-std=c++17","-fopenmp",
     "-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],
    input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
    ["g++","-###","-fno-use-linker-plugin","-fopenmp",
     f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines()
          if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}",*sys.argv[1:]],
               pass_fds=(exe,),check=True)
