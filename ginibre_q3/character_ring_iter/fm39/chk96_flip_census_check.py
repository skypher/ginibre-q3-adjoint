import os, subprocess, sys

src = r'''#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;using i128=__int128_t;using boost::multiprecision::cpp_int;const int L=40;
struct T{i128 x[L+1][L+1];};struct F{int n,e;};thread_local T dp[42];
struct DG{i128 v[L+1][L+3];unsigned s[L+1][L+3];};thread_local DG dx,dy;thread_local unsigned gen=0;
void ad(DG&d,int r,int c,i128 v,unsigned g){if(d.s[r][c]!=g){d.s[r][c]=g;d.v[r][c]=0;}d.v[r][c]+=v;}
bool cg(int a,int n,int t){return t>=abs(a-n)&&t<=a+n&&!((a+n-t)&1);}
void mul(T&o,const T&i,int w,int n,int e){
 unsigned g=++gen;if(!g){memset(dx.s,0,sizeof dx.s);memset(dy.s,0,sizeof dy.s);g=++gen;}int nw=w+n;
 for(int a=0;a<=w;a++)for(int b=0;a+b<=w;b++){i128 z=i.x[a][b];if(!z)continue;int lo=abs(a-n),hi=a+n;ad(dx,b,lo,z,g);ad(dx,b,hi+2,-z,g);lo=abs(b-n);hi=b+n;ad(dy,a,lo,(i128)e*z,g);ad(dy,a,hi+2,-(i128)e*z,g);}
 memset(o.x,0,sizeof(o.x));
 for(int b=0;b<=w;b++){i128 r[2]={0,0};for(int a=0;a<=nw-b;a++){int p=a&1;if(dx.s[b][a]==g)r[p]+=dx.v[b][a];o.x[a][b]+=r[p];}}
 for(int a=0;a<=w;a++){i128 r[2]={0,0};for(int b=0;b<=nw-a;b++){int p=b&1;if(dy.s[a][b]==g)r[p]+=dy.v[a][b];o.x[a][b]+=r[p];}}
}
string dec(i128 x){if(!x)return"0";bool neg=x<0;__uint128_t y=neg?(__uint128_t)(-(x+1))+1:(__uint128_t)x;string s;while(y){s+='0'+y%10;y/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
cpp_int big(i128 x){return cpp_int(dec(x));}
using SM=unordered_map<unsigned long long,i128>;
SM st(const vector<F>&f){int W=0;for(auto a:f)W+=a.n;int base=W+1;SM d;d.reserve(128);d[0]=1;for(auto a:f){SM q;q.reserve(d.size()*max(2,a.n+1));for(auto z:d){int x=z.first/base,y=z.first%base;for(int t=abs(x-a.n);t<=x+a.n;t+=2)q[(unsigned long long)t*base+y]+=z.second;for(int t=abs(y-a.n);t<=y+a.n;t+=2)q[(unsigned long long)x*base+t]+=(i128)a.e*z.second;}d.swap(q);}return d;}
i128 get(const SM&d,int W,int a,int b){auto i=d.find((unsigned long long)a*(W+1)+b);return i==d.end()?0:i->second;}
i128 gs(const vector<F>&f,int p){int W=0;for(auto a:f)W+=a.n;auto d=st(f);return get(d,W,p,0);}
i128 pc(const vector<F>&f,int a,int b){int W=0;for(auto x:f)W+=x.n;auto d=st(f);return get(d,W,a,b);}
i128 app(const T&d,int dw,int p,int e,int a,int b){i128 z=0;for(int r=0;r<=dw;r++)if(cg(r,p,a))z+=d.x[r][b];for(int s=0;s<=dw;s++)if(cg(s,p,b))z+=(i128)e*d.x[a][s];return z;}
bool possible(int w,int last,int lg,int mx){int need=lg>=2?0:(lg==1?max(3,last):6);return max(w+need,3*mx)<=40;}
thread_local vector<F>path(45);thread_local long long cases=0,excluded=0,desc=0,noflip=0,byW[41]={0},failD=0,failTP=0,noTP=0,seen=0,profiles=0;
thread_local vector<F>bestB;thread_local int bestW=0,bestp=0;thread_local i128 bestN=0,bestD=1;thread_local long double bestR=-1;
thread_local uint64_t rs=0x123ab987fedc7654ULL;uint64_t rr(){rs^=rs>>12;rs^=rs<<25;rs^=rs>>27;return rs*2685821657736338717ULL;}
struct Samp{vector<F>b;int p;i128 v;};thread_local vector<Samp>samples;
bool full(const vector<F>&b,int p,int sig){vector<F>all=b;all.push_back({p,sig});for(int i=0;i<(int)all.size();i++)for(int j=i+1;j<(int)all.size();j++){vector<F>c;for(int k=0;k<(int)all.size();k++)if(k!=i&&k!=j)c.push_back(all[k]);if((i128)all[j].e*pc(c,all[i].n,all[j].n)>=0)return true;}return false;}
void nof(const vector<F>&b,int W,int p,i128 parent){
 noflip++;byW[W]++;int bi=-1,bj=-1,bs=-1,bm=-1;
 for(int i=0;i<(int)b.size();i++)for(int j=i+1;j<(int)b.size();j++)if(!((b[i].n+b[j].n)&1)){int s=b[i].n+b[j].n,m=max(b[i].n,b[j].n);if(s>bs||(s==bs&&m>bm)){bi=i;bj=j;bs=s;bm=m;}}
 if(bi<0){noTP++;failTP++;}else{vector<F>c;for(int k=0;k<(int)b.size();k++)if(k!=bi&&k!=bj)c.push_back(b[k]);i128 v=gs(c,p);if(v>parent||parent<=0)failTP++;else if(big(v)*big(bestD)>big(bestN)*big(parent)){bestN=v;bestD=parent;bestW=W;bestp=p;bestB=b;}long double q=(long double)v/(long double)parent;if(q>bestR)bestR=q;}
 bool ok=false;for(int i=0;i<(int)b.size()&&!ok;i++)if(!(b[i].n&1)){vector<F>c;for(int k=0;k<(int)b.size();k++)if(k!=i)c.push_back(b[k]);if(gs(c,p)<=parent)ok=true;}
 for(int i=0;i<(int)b.size()&&!ok;i++)for(int j=i+1;j<(int)b.size()&&!ok;j++)if(!((b[i].n+b[j].n)&1)){vector<F>c;for(int k=0;k<(int)b.size();k++)if(k!=i&&k!=j)c.push_back(b[k]);if(gs(c,p)<=parent)ok=true;}
 if(!ok)failD++;
}
void node(int d,int W,int nm,int bigc,int mx){
 if(!d)return;profiles++;
 if(W>=16&&W<=40&&bigc>=2&&mx<=W/2){int lo=max({6,mx,23-W});if((lo-W)&1)lo++;int hi=min(W-16,W-2*mx);if((hi-W)&1)hi--;int sig=(nm&1)?-1:1;vector<F>b(path.begin(),path.begin()+d);
  for(int p=lo;p<=hi;p+=2){bool bad=false;for(F f:b)if(f.n==p&&f.e==-sig){bad=true;break;}if(bad){excluded++;continue;}cases++;i128 parent=dp[d].x[p][0];seen++;if(samples.size()<160)samples.push_back({b,p,parent});else{uint64_t k=rr()%seen;if(k<samples.size())samples[k]={b,p,parent};}
   bool yes=false;if((i128)b.back().e*dp[d-1].x[p][b.back().n]>=0)yes=true;
   if(!yes&&d>=2){F u=b[d-2],v=b[d-1];i128 t=app(dp[d-2],W-u.n-v.n,p,sig,u.n,v.n);if((i128)v.e*t>=0)yes=true;}
   if(!yes)yes=full(b,p,sig);if(yes)desc++;else nof(b,W,p,parent);
  }
 }
}
void walk(int d,int W,int nm,int bc,int mx){
 if(d&&!possible(W,path[d-1].n,bc,mx))return;
 node(d,W,nm,bc,mx);int last=d?path[d-1].n:0;
 if(d&&W+last<=40&&possible(W+last,last,bc+(last>=3),mx)){F f{last,path[d-1].e};path[d]=f;mul(dp[d+1],dp[d],W,last,f.e);walk(d+1,W+last,nm+(f.e<0),bc+(last>=3),mx);}
 for(int n=last+1;n<=min(13,40-W);n++)for(int e:{1,-1})if(possible(W+n,n,bc+(n>=3),max(mx,n))){path[d]={n,e};mul(dp[d+1],dp[d],W,n,e);walk(d+1,W+n,nm+(e<0),bc+(n>=3),max(mx,n));}
}
struct Root{F a,b;};
long long Gcases=0,Gexcluded=0,Gdesc=0,Gnoflip=0,GbyW[41]={0},GfailD=0,GfailTP=0,GnoTP=0,Gprofiles=0;
i128 GbestN=0,GbestD=1;int GbestW=0,Gbestp=0;vector<F>GbestB;vector<Samp>Gsamples;
int main(){
 vector<Root> roots;
 for(int n=1;n<=13;n++)for(int e:{1,-1}){
  if(n<=13&&possible(2*n,n,2*(n>=3),n))roots.push_back({{n,e},{n,e}});
  for(int n2=n+1;n2<=min(13,40-n);n2++)for(int e2:{1,-1})
   if(possible(n+n2,n2,(n>=3)+(n2>=3),n2))roots.push_back({{n,e},{n2,e2}});
 }
 long long done=0;
 #pragma omp parallel for schedule(dynamic)
 for(int ix=0;ix<(int)roots.size();ix++){
  cases=excluded=desc=noflip=failD=failTP=noTP=seen=profiles=0;fill(byW,byW+41,0);samples.clear();bestB.clear();bestW=bestp=0;bestN=0;bestD=1;bestR=-1;rs=0x123ab987fedc7654ULL+ix;gen=0;memset(dp,0,sizeof dp);memset(dx.s,0,sizeof dx.s);memset(dy.s,0,sizeof dy.s);dp[0].x[0][0]=1;
  Root r=roots[ix];path[0]=r.a;path[1]=r.b;mul(dp[1],dp[0],0,r.a.n,r.a.e);mul(dp[2],dp[1],r.a.n,r.b.n,r.b.e);
  int w=r.a.n+r.b.n,nm=(r.a.e<0)+(r.b.e<0),bc=(r.a.n>=3)+(r.b.n>=3),mx=r.b.n;
  walk(2,w,nm,bc,mx);
  #pragma omp critical
  {
   Gcases+=cases;Gexcluded+=excluded;Gdesc+=desc;Gnoflip+=noflip;GfailD+=failD;GfailTP+=failTP;GnoTP+=noTP;Gprofiles+=profiles;
   for(int w0=0;w0<=40;w0++)GbyW[w0]+=byW[w0];
   if(big(bestN)*big(GbestD)>big(GbestN)*big(bestD)){GbestN=bestN;GbestD=bestD;GbestW=bestW;Gbestp=bestp;GbestB=bestB;}
   for(auto s:samples)if(Gsamples.size()<512)Gsamples.push_back(s);
   done++;if(done%30==0)cerr<<"prefixes "<<done<<"/"<<roots.size()<<" profiles="<<Gprofiles<<" residual="<<Gcases<<" noflip="<<Gnoflip<<endl;
  }
 }
 if(Gcases!=21364508LL||Gexcluded!=994058LL||Gnoflip!=5430LL)throw runtime_error("census totals");
 long long ex[]={2,4,30,14,34,33,82,66,146,128,274,267,416,408,686,629,1102,1109};for(int w=23;w<=40;w++)if(GbyW[w]!=ex[w-23])throw runtime_error("W distribution");
 if(GfailD||GfailTP||GnoTP)throw runtime_error("D or TopPair");
 for(auto s:Gsamples){if(gs(s.b,s.p)!=s.v)throw runtime_error("sparse second evaluator");int m=0;for(auto f:s.b)m+=f.e<0;vector<F>l=s.b;l.push_back({s.p,(m&1)?-1:1});auto t=st(l);int tw=0;for(auto f:l)tw+=f.n;if(get(t,tw,0,0)!=2*s.v)throw runtime_error("Phi=2g");}
 cout<<"profiles="<<Gprofiles<<" residual="<<Gcases<<" excluded="<<Gexcluded<<" noflip="<<Gnoflip<<" desc="<<Gdesc<<" samples="<<Gsamples.size()<<"\nWcounts=";for(int w=23;w<=40;w++)cout<<GbyW[w]<<(w==40?'\n':',');
 cout<<"Dfail="<<GfailD<<" TPfail="<<GfailTP<<" noTP="<<GnoTP<<"\nratio="<<dec(GbestN)<<"/"<<dec(GbestD)<<" "<<setprecision(10)<<(double)((long double)GbestN/(long double)GbestD)<<" W="<<GbestW<<" p="<<Gbestp<<" B=";for(auto f:GbestB)cout<<(f.e>0?'+':'-')<<f.n<<',';cout<<"\n";
}'''

fd = os.memfd_create("fmchk96-census", 0)
os.set_inheritable(fd, True)
env = os.environ.copy()
env["TMPDIR"] = "/dev/shm"
env["OMP_NUM_THREADS"] = "8"
cc = subprocess.run(
    ["g++", "-O3", "-std=c++17", "-fopenmp", "-pipe", "-x", "c++",
     "-o", f"/proc/self/fd/{fd}", "-"],
    input=src, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    pass_fds=(fd,), env=env)
if cc.returncode:
    raise SystemExit(cc.stderr)
run = subprocess.run([f"/proc/self/fd/{fd}"], pass_fds=(fd,), env=env)
raise SystemExit(run.returncode)
