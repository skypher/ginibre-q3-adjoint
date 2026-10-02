import argparse,ctypes,os,resource,subprocess
ap=argparse.ArgumentParser(description="Exact FM-STR3 one-generation fusion obstruction; memory only.")
ap.parse_args()
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
CPP=r'''
#include <bits/stdc++.h>
using namespace std;using I=__int128_t;
long long ck(I x){assert(x>=LLONG_MIN&&x<=LLONG_MAX);return (long long)x;}
struct Row{vector<int>L;vector<long long>phi;long long m,cap;};
Row make_row(const vector<int>&L){
 int n=L.size(),W=accumulate(L.begin(),L.end(),0);I bound=1;
 assert(n<=11&&W<=68&&is_sorted(L.begin(),L.end()));
 for(int z:L){assert(z>=0&&z<=23);bound*=2*(z+1);}
 assert(bound<(I(1)<<62));
 vector<int>v,c;
 for(int z:L){if(v.empty()||v.back()!=z){v.push_back(z);c.push_back(1);}else ++c.back();}
 int k=v.size(),states=1,D=W+1;vector<int>stride(k);
 for(int i=0;i<k;i++){stride[i]=states;states*=c[i]+1;}
 vector<long long>mu(states*D);mu[0]=1;
 for(int q=1;q<states;q++){
  int j=0;while((q/stride[j])%(c[j]+1)==0)++j;
  int prev=q-stride[j],a=v[j];
  for(int b=0;b<=W-a;b++)if(mu[prev*D+b])
   for(int t=abs(a-b);t<=a+b;t+=2)
    mu[q*D+t]=ck(I(mu[q*D+t])+mu[prev*D+b]);
 }
 vector<int>idx(n),state(1<<n);
 for(int i=0;i<n;i++)idx[i]=lower_bound(v.begin(),v.end(),L[i])-v.begin();
 vector<long long>f(1<<n);
 for(int s=1;s<(1<<n);s++){
  int j=__builtin_ctz((unsigned)s);state[s]=state[s&(s-1)]+stride[idx[j]];
 }
 for(int s=0;s<(1<<n);s++)f[s]=ck(I(mu[state[s]*D])*mu[(states-1-state[s])*D]);
 for(int step=1;step<(1<<n);step*=2)
  for(int j=0;j<(1<<n);j+=2*step)for(int t=0;t<step;t++){
   long long a=f[j+t],b=f[j+t+step];f[j+t]=ck(I(a)+b);f[j+t+step]=ck(I(a)-b);
  }
 long long mn=LLONG_MAX;
 for(int s=0;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)%2==0)mn=min(mn,f[s]);
 assert(mn>=0);if(n)assert(mn%2==0);
 return {L,f,mu[(states-1)*D],n?mn/2:0};
}
const long long shifts[]={
68353,123406,68188,123571,70966,120793,75055,116704,78342,113417,
80214,111545,81865,109894,84380,107379,84904,106855,91084,100675,
22966,100440,42352,81054,34980,88426,41171,82235,43817,79589,
46139,77267,48055,75351,50614,72792,56594,66812,10205,61166,
78036,16616,43690,59236,72217,30474,34634,58661,67990,25069,
46696,56280,63714,27861,47294,55405,61199,29819,48194,54922,
58824,31946,49456,53946,56411,39153,51772,51474,49360,4369,
27973,53847,63218,10060,50831,58042,12546,33315,49139,54407,
25345,24752,47955,51355,18293,36139,46492,48483,20659,37815,
44990,45943,28464,41373,41137,38433,7812,20434,38547,44803,
49689,12267,22865,31127,37466,42348,45686,14694,31552,36850,
40784,42534,17737,26665,32197,36340,38914,39906,22678,25665,
32807,35913,37185,37511,27358,33126,34636,34341,32561,29737,
5930,16224,23880,35285,38963,41003,8587,18734,34080,36817,
37722,11665,20523,26036,33488,34517,35056,13479,26902,32644,
33207,32375,21964,28174,30165,29460,27252,24270,4434,13505,
19630,29404,32245,33431,33765,7733,15564,21145,25752,28836,
30752,31377,30600,9433,17583,26071,28554,29417,29672,28351,
17666,24436,26670,27586,27038,25503,22965,19895,3638,11696,
17076,21818,27618,28628,28533,27407,5785,13786,18536,27003,
27352,26434,24840,13677,21223,24104,25460,24294,22137,19347,
16172,3389,10451,15285,19414,24527,25105,25180,23928,21802,
11985,18286,21702,23325,23895,23186,21472,18958,16012,12938,
9318,16092,19840,22032,23085,21144,18767,15941,12900,9962};
extern "C" int verify(){
 const vector<int>L={1,2,3,4,5,6,7,8,9,10,13};
 const int n=L.size(),mask=340;const long long t=191759;
 Row root=make_row(L);assert(root.m==222024&&root.phi[mask]==388174);
 assert(root.cap==194087);
 map<vector<int>,int>id;vector<Row>rows;vector<vector<int>>channels;
 vector<pair<int,int>>pairs;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
  pairs.push_back({i,j});vector<int>chans;
  for(int c=abs(L[i]-L[j]);c<=L[i]+L[j];c+=2){
   vector<int>Q;for(int a=0;a<n;a++)if(a!=i&&a!=j)Q.push_back(L[a]);
   Q.push_back(c);sort(Q.begin(),Q.end());
   auto it=id.find(Q);int q;
   if(it==id.end()){q=rows.size();id[Q]=q;rows.push_back(make_row(Q));}else q=it->second;
   chans.push_back(q);
  }channels.push_back(chans);
 }
 assert(rows.size()==240&&size(shifts)==240&&pairs.size()==55);
 long long shorter_checks=0;
 vector<vector<long long>>modified(rows.size());
 for(int q=0;q<(int)rows.size();q++){
  const auto&R=rows[q];assert(0<=shifts[q]&&shifts[q]<=R.cap&&R.cap<=R.m);
  modified[q]=R.phi;
  for(int s=0;s<(int)R.phi.size();s++){
   if(__builtin_popcount((unsigned)s)%2==0){
    modified[q][s]=ck(I(R.phi[s])-2*shifts[q]);assert(modified[q][s]>=0);++shorter_checks;
   }else assert(R.phi[s]==0);
  }
 }
 vector<long long>rp=root.phi;
 for(int s=0;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)%2==0){
  rp[s]=ck(I(rp[s])-2*t);assert(rp[s]>=4656);
 }
 long long maxD=LLONG_MIN;
 for(auto[i,j]:pairs){
  long long a=ck(I(root.phi[mask])-root.phi[mask^(1<<i)^(1<<j)]);
  assert(a%4==0&&a==rp[mask]-rp[mask^(1<<i)^(1<<j)]&&a<0);
  maxD=max(maxD,a/4);
 }
 assert(maxD==-185);
 long long fusion_checks=0;
 for(int q=0;q<(int)pairs.size();q++){
  auto[i,j]=pairs[q];long long a=0,b=0;
  for(int z:channels[q]){a=ck(I(a)+rows[z].m);b=ck(I(b)+rows[z].m-shifts[z]);}
  assert(a==root.m&&b==root.m-t);
  for(int s=0;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)%2==0){
   long long total=0;int ci=0;
   for(int c=abs(L[i]-L[j]);c<=L[i]+L[j];c+=2){
    vector<pair<int,int>>tokens;
    for(int k=0;k<n;k++)if(k!=i&&k!=j)tokens.push_back({L[k],(s>>k)&1});
    tokens.push_back({c,((s>>i)^(s>>j))&1});
    stable_sort(tokens.begin(),tokens.end(),[](auto a,auto b){return a.first<b.first;});
    int sm=0;for(int k=0;k<n-1;k++)sm|=tokens[k].second<<k;
    total=ck(I(total)+modified[channels[q][ci++]][sm]);
   }
   assert(I(rp[s])+rp[s^(1<<i)^(1<<j)]==2*I(total));++fusion_checks;
  }
 }
 pair<int,int>best{-1,-1};int ia=-1,ib=-1;
 for(int i=0;i<n-1;i++)for(int j=i+1;j<n-1;j++)if((L[i]-L[j])%2==0){
  auto key=make_pair(L[i]+L[j],max(L[i],L[j]));
  if(key>best)best=key,ia=i,ib=j;
 }
 assert(L[ia]==8&&L[ib]==10&&((mask>>ia)&1)==0&&((mask>>ib)&1)==0);
 vector<int>child;int cm=0;
 for(int i=0;i<n;i++)if(i!=ia&&i!=ib){cm|=((mask>>i)&1)<<child.size();child.push_back(L[i]);}
 auto cr=make_row(child);assert(cr.phi[cm]==6172);
 assert(rp[mask]/2==2328&&rp[mask]/2-cr.phi[cm]/2==-758);
 long long proper=0;
 for(int sel=0;sel<(1<<n)-1;sel++){
  vector<int>Q;for(int i=0;i<n;i++)if((sel>>i)&1)Q.push_back(L[i]);
  auto R=make_row(Q);auto it=id.find(Q);long long sub=it==id.end()?0:shifts[it->second];
  for(int s=0;s<(int)R.phi.size();s++)if(__builtin_popcount((unsigned)s)%2==0){
   assert(R.phi[s]-2*sub>=0);++proper;
  }
 }
 assert(proper==87550&&shorter_checks==122880&&fusion_checks==56320);
 const auto &Q=rows[0];assert(Q.L==vector<int>({1,3,4,5,6,7,8,9,10,13}));
 long long lower=0;
 for(int c=2;c<=4;c+=2){
  vector<int>A(Q.L.begin()+2,Q.L.end());A.push_back(c);sort(A.begin(),A.end());
  lower=ck(I(lower)+make_row(A).m);
 }
 assert(lower==Q.m&&lower-(Q.m-shifts[0])==68353);
 cout<<"root m: "<<root.m<<" -> "<<root.m-t<<"; g: 194087 -> 2328; child=3086\n"
     <<"all 55 flips unchanged; max D=-185; TopPair drop=-758\n"
     <<"children=240; child-sign checks="<<shorter_checks
     <<"; proper-sign checks="<<proper<<"; signed fusion checks="<<fusion_checks<<"\n"
     <<"second-fusion defect=68353; child invariant "<<Q.m<<" -> "<<Q.m-shifts[0]<<endl;
 return 0;
}
'''
obj=os.memfd_create("fm-str3-object")
lib=os.memfd_create("fm-str3-library")
env={**os.environ,"TMPDIR":"/dev/shm"}
subprocess.run(["g++","-std=c++17","-O2","-pipe","-fPIC","-x","c++","-c","-",
 "-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,env=env,pass_fds=(obj,),check=True)
subprocess.run(["ld","-shared",f"/proc/self/fd/{obj}",
 "/lib/x86_64-linux-gnu/libstdc++.so.6","-lc","-o",f"/proc/self/fd/{lib}"],
 env=env,pass_fds=(obj,lib),check=True)
dll=ctypes.CDLL(f"/proc/self/fd/{lib}")
assert dll.verify()==0
os.close(obj);os.close(lib)
print("ALL EXACT CHECKS PASS")
