import os,sys,subprocess
from datetime import datetime,timezone
from math import factorial as fac,prod
from fractions import Fraction as Q
def log(*v):
 print(datetime.now(timezone.utc).isoformat(timespec='seconds'),*v,flush=True)
mode=sys.argv[1] if len(sys.argv)>1 else 'all'
if mode in ('-h','--help'):
 print('Usage: python3 -u -B - [all|box|named|family|run9]. Default all includes the box and named runs through k=9.')
 raise SystemExit
assert mode in ('all','box','named','family','run9')
def chars(w):
 d={(0,0):1}
 for v in w:
  n=abs(v);q={}
  for(a,b),z in d.items():
   for c in range(abs(a-n),a+n+1,2):q[c,b]=q.get((c,b),0)+z
   for c in range(abs(b-n),b+n+1,2):q[a,c]=q.get((a,c),0)+(1 if v>0 else -1)*z
  d={k:v for k,v in q.items()if v}
 return d
def hd(w):
 A,B=chars(w[::2]),chars(w[1::2])
 p=[v*B.get(k,0)for k,v in A.items()]
 return sum(max(z,0)for z in p),sum(max(-z,0)for z in p)
for a in range(1,15):
 for b in range(a,15):
  assert hd((-a,)*3+(-b,)*3)[1]==(2 if a%2==b%2==0 and a<b<=2*a else 0)
def triangle(a,b,c):
 return Q(fac(a+b-c)*fac(a+c-b)*fac(b+c-a),fac(a+b+c+1))
for s in range(1,7):
 for t in range(s+1,2*s+1):
  lows=[2*s+2*t,2*s+t,s+2*t,s+3*t]
  highs=[s+3*t,2*s+3*t,3*s+2*t]
  assert max(lows)==min(highs)==s+3*t
  z=s+3*t
  term=Q((-1)**z*fac(z+1),prod(fac(z-x)for x in lows)*prod(fac(x-z)for x in highs))
  r2=prod(triangle(*x)for x in [(s,t,s+t),(s,t,s),(t,t,s),(t,t,s+t)])*term**2
  closed=Q(fac(2*s)*fac(2*t)*fac(s+t)**2*fac(s+3*t+1),
           fac(2*s+2*t+1)*fac(2*s+t+1)*fac(s+2*t+1)*fac(t-s)*fac(2*t-s)*fac(2*s-t))
  assert r2==closed>0

  def nu(a,b,c):
   h=(a+b-c)//2
   return Q(fac(h)*fac(a+b-h+1)*fac(a-h)*fac(b-h),fac(a)*fac(b)*fac(c+1))
  a,b=2*s,2*t
  U=Q(fac(s)*fac(t)*fac(s+t)**2*fac(s+3*t+1),fac(2*s)*fac(2*t)**2*fac(2*s+2*t))
  assert U**2==(a+1)*(b+1)*(a+b+1)**2*nu(a,b,a)*nu(a,b,b)*nu(b,b,a+b)*r2
log('PASS 105 channel classifications and 21 exact minor formulas')
if mode=='all':
 for w,expected in [
  ((1,1,-2,3,3,3,3,-4,-4,-4,5,5,5,5,8),(8437374,286632)),
  ((-1,-1,-2,-3,-3,-3,-3,-4,-4,-4,-5,-5,-5,-5,8),(8437374,286632)),
  ((-40,42,-44,46,48,50,52,54,56,58,60,62),(908860623257408,2445554811680))]:
  log('DIMENSION START',w)
  got=hd(w);assert got==expected
  log('DIMENSION ONLY',w,'H(d0)=',got,'Phi=',got[0]-got[1])

src=r'''
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <omp.h>
using namespace std;
using ll=long long;
constexpr int MOD=1000003;
int plusm(int a,int b){int v=a+b;return v>=MOD?v-MOD:v;}
int mm(int a,int b){return int(ll(a)*b%MOD);}
int power(int a,int n){int r=1;for(;n;n>>=1,a=mm(a,a))if(n&1)r=mm(r,a);return r;}
int inv(int a){assert(a);return power(a,MOD-2);}
int binom(int n,int k){if(k<0||k>n)return 0;int r=1;for(int i=1;i<=k;i++)r=mm(mm(r,n-i+1),inv(i));return r;}
string timestamp(){time_t now=time(nullptr);tm b{};gmtime_r(&now,&b);char s[40];strftime(s,sizeof(s),"%Y-%m-%dT%H:%M:%SZ",&b);return s;}
void logline(const string&s){
 #pragma omp critical(output)
 cout<<timestamp()<<" "<<s<<endl;
}
string wordstr(const vector<int>&w){ostringstream s;s<<"(";for(int i=0;i<(int)w.size();i++){if(i)s<<",";s<<w[i];}s<<")";return s.str();}
vector<int> stride(const vector<int>&ns){vector<int>s(ns.size(),1);for(int i=int(ns.size())-2;i>=0;i--)s[i]=s[i+1]*(ns[i+1]+1);return s;}
using Poly=vector<pair<int,int>>;
using AB=pair<int,int>;
struct Multi{int a;vector<int>path;vector<Poly>states;};
thread_local map<vector<int>,vector<Multi>> CG;
Poly compact(const map<int,int>&p){Poly q;for(auto[v,c]:p)if(c)q.push_back({v,c});return q;}
Poly lowerp(const Poly&p,const vector<int>&ns){
 auto st=stride(ns);map<int,int>q;
 for(auto[v,c]:p)for(int i=0;i<(int)ns.size();i++){
  int h=v/st[i]%(ns[i]+1);if(h<ns[i])q[v+st[i]]=plusm(q[v+st[i]],mm(c,ns[i]-h));
 }
 return compact(q);
}
const vector<Multi>&cg(const vector<int>&ns){
 auto it=CG.find(ns);if(it!=CG.end())return it->second;
 vector<Multi>out;
 if(ns.empty()){out.push_back({0,{},{{{0,1}}}});return CG.emplace(ns,move(out)).first->second;}
 vector<int>old(ns.begin(),ns.end()-1);int n=ns.back();const auto&prev=cg(old);
 for(const auto&v:prev)for(int c=abs(v.a-n);c<=v.a+n;c+=2){
  int t=(v.a+n-c)/2;map<int,int>hv;
  for(int h=0;h<=t;h++){int coef=binom(t,h);if(h%2)coef=MOD-coef;
   for(auto[z,d]:v.states[h])hv[z*(n+1)+t-h]=plusm(hv[z*(n+1)+t-h],mm(coef,d));
  }
  Multi u;u.a=c;u.path=v.path;u.path.push_back(c);u.states.push_back(compact(hv));
  for(int h=0;h<c;h++){auto q=lowerp(u.states.back(),ns);int s=inv(c-h);for(auto&[v,d]:q)d=mm(d,s);u.states.push_back(move(q));}
  assert(lowerp(u.states.back(),ns).empty());out.push_back(move(u));
 }
 ll dim=1,sum=0;for(int n:ns)dim*=n+1;for(auto&m:out)sum+=m.a+1;assert(sum==dim);
 return CG.emplace(ns,move(out)).first->second;
}
struct Copy{int mask,a,b,ix,iy;};
struct Half{
 vector<int>ns,pos,st,metric,weight;int L,dim;
 vector<Copy>copies;map<AB,vector<int>>bytype;
 vector<vector<int>>xl,yl,xpos,ypos;
 map<pair<int,int>,Poly>polycache;
 unordered_map<ll,int>gramcache;
 bool scaled;
 Half(vector<int>n,vector<int>ps,bool scale):ns(n),pos(ps),scaled(scale){
  L=ns.size();st=stride(ns);dim=1;for(int n:ns)dim*=n+1;
  metric.assign(dim,1);
  for(int v=0;v<dim;v++)for(int i=0;i<L;i++)metric[v]=mm(metric[v],inv(binom(ns[i],v/st[i]%(ns[i]+1))));
  weight.assign(1<<L,1);xl.resize(1<<L);yl.resize(1<<L);xpos.resize(1<<L);ypos.resize(1<<L);
  for(int mask=0;mask<(1<<L);mask++){
   for(int i=0;i<L;i++)if(mask>>i&1){xl[mask].push_back(ns[i]);xpos[mask].push_back(i);}else{yl[mask].push_back(ns[i]);ypos[mask].push_back(i);}
   if(mask){int b=__builtin_ctz((unsigned)mask);weight[mask]=mm(weight[mask^(1<<b)],pos[b]+2);}
   const auto&A=cg(xl[mask]);const auto&B=cg(yl[mask]);
   for(int x=0;x<(int)A.size();x++)for(int y=0;y<(int)B.size();y++){
    int id=copies.size();copies.push_back({mask,A[x].a,B[y].a,x,y});bytype[{A[x].a,B[y].a}].push_back(id);
   }
  }
 }
 Poly embed(const Poly&p,const vector<int>&small,const vector<int>&where){
  auto ss=stride(small);Poly q;
  for(auto[v,c]:p){int z=0;for(int i=0;i<(int)small.size();i++)z+=(v/ss[i]%(small[i]+1))*st[where[i]];q.push_back({z,c});}
  return q;
 }
 int nu(int a,int b,int j){
  int t=(a+b-j)/2,r=0;for(int h=0;h<=t;h++)r=plusm(r,mm(mm(binom(t,h),binom(t,h)),inv(mm(binom(a,h),binom(b,t-h)))));
  return r;
 }
 const Poly& coupled(int id,int j){
  pair<int,int>key{id,j};auto it=polycache.find(key);if(it!=polycache.end())return it->second;
  const auto&c=copies[id];assert(j>=abs(c.a-c.b)&&j<=c.a+c.b&&((c.a+c.b-j)%2==0));
  const auto&A=cg(xl[c.mask])[c.ix];const auto&B=cg(yl[c.mask])[c.iy];int t=(c.a+c.b-j)/2;
  map<int,int>ans;
  for(int h=0;h<=t;h++){
   auto p=embed(A.states[h],xl[c.mask],xpos[c.mask]);auto q=embed(B.states[t-h],yl[c.mask],ypos[c.mask]);
   int coef=binom(t,h);if(h%2)coef=MOD-coef;
   for(auto[v,x]:p)for(auto[w,y]:q)ans[v+w]=plusm(ans[v+w],mm(coef,mm(x,y)));
  }
  Poly p=compact(ans);if(scaled){int s=inv(nu(c.a,c.b,j));for(auto&[v,c]:p)c=mm(c,s);}
  return polycache.emplace(key,move(p)).first->second;
 }
 int gram(int u,int v,int j){
  if(u>v)swap(u,v);ll key=(ll(j)*copies.size()+u)*copies.size()+v;
  auto it=gramcache.find(key);if(it!=gramcache.end())return it->second;
  const auto&p=coupled(u,j);const auto&q=coupled(v,j);int a=0,b=0,val=0;
  while(a<(int)p.size()&&b<(int)q.size()){
   if(p[a].first<q[b].first)a++;
   else if(p[a].first>q[b].first)b++;
   else{val=plusm(val,mm(mm(p[a].second,q[b].second),metric[p[a].first]));a++;b++;}
  }
  return gramcache.emplace(key,val).first->second;
 }
 int wt(int u,int v){return weight[(~copies[u].mask)&(~copies[v].mask)&((1<<L)-1)];}
 map<AB,pair<int,vector<int>>>surplus(const vector<int>&w){
  map<AB,pair<int,vector<int>>>out;int M=0;for(int i=0;i<L;i++)if(w[i]<0)M|=1<<i;
  for(auto&[ab,ids]:bytype){
   vector<int>es[2];for(int id:ids)es[__builtin_popcount((unsigned)(M&~copies[id].mask))%2].push_back(id);
   // Generation is mask, x-path, y-path lexicographic order.
   int r=min(es[0].size(),es[1].size());
   for(int e=0;e<2;e++)if((int)es[e].size()>r)out[ab]={e,vector<int>(es[e].begin()+r,es[e].end())};
  }
  return out;
 }
};
map<AB,ll> character(const vector<int>&w){
 map<AB,ll>d{{{0,0},1}};
 for(int v:w){int n=abs(v),sg=v>0?1:-1;map<AB,ll>q;
  for(auto&[ab,z]:d){auto[a,b]=ab;for(int c=abs(a-n);c<=a+n;c+=2)q[{c,b}]+=z;for(int c=abs(b-n);c<=b+n;c+=2)q[{a,c}]+=sg*z;}
  d.clear();for(auto&[ab,z]:q)if(z)d[ab]=z;
 }
 return d;
}
pair<ll,ll> homdims(const vector<int>&w){
 vector<int>a,b;for(int i=0;i<(int)w.size();i++)(i%2?b:a).push_back(w[i]);
 auto A=character(a),B=character(b);ll e=0,o=0;
 for(auto&[ab,v]:A){ll z=v*B[ab];if(z>0)e+=z;else o-=z;}return {e,o};
}
struct Group{AB ab;vector<int>a,b;int offset;};
struct Job{vector<int>w;ll he,ho;int multiplicity;};
struct Result{int rank;int rows;bool full;};

#include <cblas.h>
#include <random>
struct BlockRank {
 int n,r=0;vector<int>perm;vector<double>U;
 BlockRank(int n):n(n),perm(n),U(ll(n)*n){
  iota(perm.begin(),perm.end(),0);
  assert(ll(n)*(MOD-1)*(MOD-1)<(1LL<<53));
 }
 void feed(const vector<vector<int>>&raw){
  int q=raw.size(),tail=n-r;if(!q||!tail)return;
  vector<double>B(ll(q)*n);
  for(int i=0;i<q;i++)for(int j=0;j<n;j++)B[ll(i)*n+j]=raw[i][perm[j]];
  if(r){
   vector<double>Z(ll(q)*tail);
   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,q,tail,r,1,
     B.data(),n,U.data()+r,n,0,Z.data(),tail);
   for(int i=0;i<q;i++)for(int j=0;j<tail;j++){
    double z=Z[ll(i)*tail+j];assert(z>=0&&z<(1LL<<53)&&z==double(ll(z)));
    int x=int(B[ll(i)*n+r+j]),y=ll(z)%MOD;
    B[ll(i)*n+r+j]=x>=y?x-y:x+MOD-y;
   }
   for(int i=0;i<q;i++)fill(B.begin()+ll(i)*n,B.begin()+ll(i)*n+r,0);
  }
  int k=0;
  for(;k<q&&r+k<n;k++){
   int pc=-1,pr=-1;
   for(int c=r+k;c<n&&pc<0;c++)for(int i=k;i<q;i++)if(B[ll(i)*n+c]){pc=c;pr=i;break;}
   if(pc<0)break;
   if(pc!=r+k){
    swap(perm[pc],perm[r+k]);
    for(int i=0;i<r;i++)swap(U[ll(i)*n+pc],U[ll(i)*n+r+k]);
    for(int i=0;i<q;i++)swap(B[ll(i)*n+pc],B[ll(i)*n+r+k]);
   }
   if(pr!=k)for(int c=0;c<n;c++)swap(B[ll(pr)*n+c],B[ll(k)*n+c]);
   int f=inv(int(B[ll(k)*n+r+k]));
   for(int c=r+k;c<n;c++)B[ll(k)*n+c]=mm(int(B[ll(k)*n+c]),f);
   for(int i=0;i<q;i++)if(i!=k&&B[ll(i)*n+r+k]){
    int v=B[ll(i)*n+r+k];B[ll(i)*n+r+k]=0;
    for(int c=r+k+1;c<n;c++){
     int x=B[ll(i)*n+c],y=mm(v,int(B[ll(k)*n+c]));
     B[ll(i)*n+c]=x>=y?x-y:x+MOD-y;
    }
   }
  }
  int rem=n-r-k;
  if(r&&k&&rem){
   vector<double>Z(ll(r)*rem);
   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,r,rem,k,1,
     U.data()+r,n,B.data()+r+k,n,0,Z.data(),rem);
   for(int i=0;i<r;i++)for(int j=0;j<rem;j++){
    double z=Z[ll(i)*rem+j];assert(z>=0&&z<(1LL<<53)&&z==double(ll(z)));
    int x=U[ll(i)*n+r+k+j],y=ll(z)%MOD;
    U[ll(i)*n+r+k+j]=x>=y?x-y:x+MOD-y;
   }
  }
  for(int i=0;i<r;i++)fill(U.begin()+ll(i)*n+r,U.begin()+ll(i)*n+r+k,0);
  for(int i=0;i<k;i++)copy(B.begin()+ll(i)*n,B.begin()+ll(i+1)*n,U.begin()+ll(r+i)*n);
  r+=k;
 }
};
int simple_rank(vector<vector<int>>a){
 int m=a.size(),n=a[0].size(),r=0;
 for(int j=0;j<n&&r<m;j++){
  int p=r;while(p<m&&!a[p][j])p++;if(p==m)continue;swap(a[p],a[r]);
  int v=inv(a[r][j]);for(int t=j;t<n;t++)a[r][t]=mm(a[r][t],v);
  for(int i=r+1;i<m;i++)if(a[i][j]){
   v=a[i][j];for(int t=j;t<n;t++){int z=mm(v,a[r][t]);a[i][t]=a[i][t]>=z?a[i][t]-z:a[i][t]+MOD-z;}
  }r++;
 }return r;
}
void test_block(){
 mt19937 gen(1234);
 for(int t=0;t<30;t++){
  int m=35+t*4,n=20+t*3,k=1+gen()%min(m,n);
  vector<vector<int>>x(m,vector<int>(k)),y(k,vector<int>(n)),a(m,vector<int>(n));
  for(auto&v:x)for(auto&z:v)z=gen()%MOD;
  for(auto&v:y)for(auto&z:v)z=gen()%MOD;
  for(int i=0;i<m;i++)for(int j=0;j<n;j++)for(int z=0;z<k;z++)a[i][j]=plusm(a[i][j],mm(x[i][z],y[z][j]));
  BlockRank R(n);
  for(int i=0;i<m;i+=16)R.feed(vector<vector<int>>(a.begin()+i,a.begin()+min(i+16,m)));
  assert(R.r==simple_rank(a));
 }
 logline("PASS 30 block-rank comparisons against scalar arithmetic");
}
Result evaluate(const Job&job){
 vector<int>na,nb,pa,pb,wa,wb;const auto&w=job.w;
 for(int i=0;i<(int)w.size();i++)if(i%2){nb.push_back(abs(w[i]));pb.push_back(i);wb.push_back(w[i]);}else{na.push_back(abs(w[i]));pa.push_back(i);wa.push_back(w[i]);}
 Half A(na,pa,true),B(nb,pb,false);auto aa=A.surplus(wa),bb=B.surplus(wb);
 vector<Group>src,dst;int ne=0,no=0;
 for(auto&[ab,x]:aa)if(bb.count(ab)){
  auto&y=bb.at(ab);int p=x.first^y.first;auto&vec=p?src:dst;int&ct=p?no:ne;
  vec.push_back({ab,x.second,y.second,ct});ct+=x.second.size()*y.second.size();
 }
 assert(ne==job.he&&no==job.ho);
 vector<tuple<int,int,int>>targets;
 for(int g=0;g<(int)dst.size();g++)for(int ta:dst[g].a)for(int tb:dst[g].b)targets.push_back({g,ta,tb});
 mt19937 gen(314159);shuffle(targets.begin(),targets.end(),gen);
 BlockRank R(no);vector<vector<int>>batch;int rows=0;
 for(auto[g,ta,tb]:targets){
  const auto&target=dst[g];vector<int>row(no);
  for(const auto&source:src){
   int lo=max(abs(target.ab.first-target.ab.second),abs(source.ab.first-source.ab.second));
   int hi=min(target.ab.first+target.ab.second,source.ab.first+source.ab.second);
   for(int j=lo;j<=hi;j+=2){
    vector<int>ga,gb;for(int sa:source.a)ga.push_back(mm(A.gram(ta,sa,j),A.wt(ta,sa)));
    for(int sb:source.b)gb.push_back(mm(B.gram(tb,sb,j),B.wt(tb,sb)));
    for(int ia=0;ia<(int)ga.size();ia++)if(ga[ia]){
     int v=mm(j+1,ga[ia]);int off=source.offset+ia*gb.size();
     for(int ib=0;ib<(int)gb.size();ib++)if(gb[ib])row[off+ib]=plusm(row[off+ib],mm(v,gb[ib]));
    }
   }
  }
  batch.push_back(move(row));rows++;
  if(batch.size()==64||rows==ne){R.feed(batch);batch.clear();}
  if(rows%512==0)logline("FAST "+wordstr(w)+" rows="+to_string(rows)+" rank="+to_string(R.r)+"/"+to_string(no));
  if(R.r==no)return{R.r,rows,true};
 }
 return{R.r,rows,false};
}
int box_main(int argc,char**argv){
 for(int i=1;i<argc;i++)if(string(argv[i])=="-h"||string(argv[i])=="--help"){
  cout<<"Usage: exchange-screen [max-label<=4] [max-length<=10] [threads]; fixed STR8b construction, exact modular ranks.\n";return 0;}
 int maxn=argc>1?stoi(argv[1]):4,maxL=argc>2?stoi(argv[2]):8,threads=argc>3?stoi(argv[3]):4;
 assert(maxn>=1&&maxn<=4&&maxL>=0&&maxL<=10);omp_set_num_threads(threads);
 vector<Job>jobs;ll profiles=0,basepass=0;set<vector<int>>done;
 for(int L=0;L<=maxL;L++){
  vector<int>w(L);function<void(int,int)>gen=[&](int i,int low){
   if(i<L){for(int n=low;n<=maxn;n++){w[i]=n;gen(i+1,n);}return;}
   int present=0;for(int n:w)present|=1<<(n-1);
   for(int M=0;M<(1<<maxn);M++)if((M&present)==M){
    auto signedw=w;int neg=0;for(int&i:signedw)if(M>>(i-1)&1){i=-i;neg++;}
    if(neg%2)continue;profiles++;auto[e,o]=homdims(signedw);
    if(!o){basepass++;continue;}
    auto refl=signedw;for(int&i:refl)if(abs(i)%2)i=-i;
    auto key=min(signedw,refl);
    if(!done.insert(key).second)continue;
    jobs.push_back({key,e,o,key==refl&&key==signedw?1:2});
   }
  };gen(0,1);
 }
 logline("START profiles="+to_string(profiles)+" basepass="+to_string(basepass)+" rank_jobs="+to_string(jobs.size()));
 int passed=0,deficient=0,complete=0; vector<string> defects;
 #pragma omp parallel for schedule(dynamic)
 for(int q=0;q<(int)jobs.size();q++){
  auto result=evaluate(jobs[q]);
  #pragma omp critical(results)
  {
   complete++;
   if(result.full)passed+=jobs[q].multiplicity;else deficient+=jobs[q].multiplicity;
   string msg=string(result.full?"PASS ":"DEFICIENT ")+wordstr(jobs[q].w)+" H0="+to_string(jobs[q].he)+","+to_string(jobs[q].ho)+" rank_mod="+to_string(result.rank)+" rows="+to_string(result.rows)+" job="+to_string(complete)+"/"+to_string(jobs.size()); if(!result.full){defects.push_back(msg);logline(msg);}else if(complete%16==0||jobs[q].ho>=1000)logline(msg);
  }
 }
 for(auto&s:defects)logline("FINAL "+s);
 logline("COMPLETE profiles="+to_string(profiles)+" certified="+to_string(basepass+passed)+" modular_deficient="+to_string(deficient));
 assert(deficient==0&&basepass+passed==profiles);return 0;
}

int entry(Half&A,Half&B,AB beta,int ta,int tb,AB alpha,int sa,int sb){
 int z=0,lo=max(abs(beta.first-beta.second),abs(alpha.first-alpha.second));
 int hi=min(beta.first+beta.second,alpha.first+alpha.second);
 for(int j=lo;j<=hi;j+=2)z=plusm(z,mm(j+1,mm(A.gram(ta,sa,j),B.gram(tb,sb,j))));
 return mm(z,mm(A.wt(ta,sa),B.wt(tb,sb)));
}
void family(){
 for(auto [a,b]:vector<AB>{{2,4},{4,6},{4,8},{6,8}}){
  Half A({a,a,b},{0,2,4},true),B({a,b,b},{1,3,5},false);
  auto aa=A.surplus({-a,-a,-b}),bb=B.surplus({-a,-b,-b});
  AB beta{a+b,a};auto pick=[&](int mask){
   for(int u:aa.at(beta).second)if(A.copies[u].mask==mask)return u;
   assert(false);return -1;
  };
  int tb=bb.at(beta).second.at(0);assert(B.copies[tb].mask==6);
  vector<int>ts{pick(6),pick(5)};vector<AB>as{{a,b},{b,a}};int mat[2][2];
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){
   auto alpha=as[j];assert(aa.at(alpha).second.size()==1&&bb.at(alpha).second.size()==1);
   int sa=aa.at(alpha).second[0],sb=bb.at(alpha).second[0];
   assert(A.copies[sa].mask==(j?6:2)&&B.copies[sb].mask==(j?4:5));
   mat[i][j]=entry(A,B,beta,ts[i],tb,alpha,sa,sb);
  }
  assert(!mat[0][1]&&!mat[1][0]&&mat[0][0]&&mat[1][1]);
  auto h=homdims({-a,-a,-a,-b,-b,-b});assert(h.second==2);
  logline("FAMILY a="+to_string(a)+" b="+to_string(b)+" H0="+to_string(h.first)+",2 rank=2");
 }
}
void named(bool ninth){
 vector<vector<int>>words;
 if(ninth)words.push_back({-1,-2,-3,-4,-5,-6,-7,-8,-9,-9});
 else{
  words.push_back({-1,2,3,4,-5,6,7,8});
  for(int k=3;k<=9;k++){
   vector<int>w;for(int n=1;n<=k;n++)w.push_back(-n);
   int W=k*(k+1)/2,p=max(k,6);
   while((W-p)%2==1||(!(k%2)&&p<=k))p++;
   w.push_back(k%2?-p:p);words.push_back(w);
  }
 }
 for(auto&w:words){
  auto[e,o]=homdims(w);logline("NAMED START "+wordstr(w)+" H0="+to_string(e)+","+to_string(o));
  Result r=o?evaluate({w,e,o,1}):Result{0,0,true};assert(r.full);
  logline("NAMED PASS "+wordstr(w)+" rank="+to_string(r.rank)+" H="+to_string(e-o)+",0");
 }
}
int main(int argc,char**argv){
 string mode=argc>1?argv[1]:"all";
 if(mode=="--help"||mode=="-h"){
  cout<<"Usage: verifier [all|box|named|family|run9]. all = full 4521 box, named runs through k=9, family minors. run9 is a longer extra check.\n";return 0;
 }
 assert(mode=="all"||mode=="box"||mode=="named"||mode=="family"||mode=="run9");
 for(int d=2;ll(d)*d<=MOD;d++)assert(MOD%d);
 omp_set_num_threads(4);openblas_set_num_threads(1);test_block();
 if(mode=="all"||mode=="family")family();
 if(mode=="all"||mode=="box"){
  char a0[]="box",a1[]="4",a2[]="10",a3[]="4";char*av[]={a0,a1,a2,a3};
  assert(box_main(4,av)==0);
 }
 if(mode=="all"||mode=="named"){openblas_set_num_threads(4);named(false);}
 if(mode=="run9"){openblas_set_num_threads(4);named(true);}
 logline("PASS selected mode "+mode);
}
'''

log('compile in memory')
fd=os.memfd_create('str8c-verifier',0)
env=os.environ.copy();env['TMPDIR']='/dev/shm'
try:
 subprocess.run(['g++','-O3','-std=c++17','-fopenmp',
                 '-I/usr/include/x86_64-linux-gnu/openblas-pthread',
                 '-x','c++','-o',f'/proc/self/fd/{fd}','-','-lopenblas'],
                input=src.encode(),pass_fds=(fd,),env=env,check=True)
 subprocess.run([f'/proc/self/fd/{fd}',mode],pass_fds=(fd,),check=True)
finally:os.close(fd)
log('PASS completed verifier')
