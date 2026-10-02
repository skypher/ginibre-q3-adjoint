#include <algorithm>
#include <chrono>
#include <climits>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <omp.h>
using namespace std;
using I=boost::multiprecision::checked_int128_t;
using Q=boost::rational<long long>;
using B=vector<int>;
static long long tick_ms(){return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count();}
static I chooseI(int n,int k){if(k<0||n<k)return 0;k=min(k,n-k);I z=1;for(int i=1;i<=k;i++){z*=n-k+i;z/=i;}return z;}
static long long chooseLL(int n,int k){if(k<0||n<k)return 0;k=min(k,n-k);long long z=1;for(int i=1;i<=k;i++)z=z*(n-k+i)/i;return z;}
static int pop(unsigned x){return __builtin_popcount(x);}
static I absI(I x){return x<0?-x:x;}
static string sI(const I&x){return x.str();}
static string sQ(const Q&q){return to_string(q.numerator())+"/"+to_string(q.denominator());}
static int signpar(const B&v){int k=0;for(int x:v)if(x<0)k++;return k&1;}
static vector<int> cg(int a,int b){vector<int>v;for(int k=abs(a-b);k<=a+b;k+=2)v.push_back(k);return v;}

// Multiplicity table from the Weyl coefficient extraction formula.
static vector<I> inv_table(const B& n){
 int N=n.size(),Z=1<<N,F=Z-1;vector<int>w(Z),sz(Z);
 for(int s=1;s<Z;s++){int bit=s&-s,i=__builtin_ctz((unsigned)bit);w[s]=w[s^bit]+n[i];sz[s]=sz[s^bit]+1;}
 vector<I> m(Z);m[0]=1;
 for(int s=1;s<Z;s++){
  int k=sz[s],W=w[s];if(k==1){m[s]=(n[__builtin_ctz((unsigned)s)]==0?I(1):I(0));continue;}
  if(W&1)continue;
  I v=0;for(int j=s;;j=(j-1)&s){
   int q=W/2-w[j]-sz[j];
   if(q>=0){I z=chooseI(q+k-2,k-2);v+=(pop((unsigned)j)&1)?-z:z;}
   if(j==0)break;
  }
  if(v<0){fprintf(stderr,"negative invariant multiplicity mask=%d\n",s);abort();}
  m[s]=v;
 }
 return m;
}
static I omega(const vector<I>&m,int s){I z=0;for(int a=s;;a=(a-1)&s){z+=m[a]*m[s^a];if(!a)break;}return z;}
struct TopPair{int i=-1,j=-1;};
static TopPair top_pair(const B&n,int upto){
 TopPair t;int bs=-1,bm=-1;
 for(int i=0;i<upto;i++)for(int j=i+1;j<upto;j++)if(((n[i]+n[j])&1)==0){
  int sm=n[i]+n[j],mx=max(n[i],n[j]);
  if(sm>bs||(sm==bs&&mx>bm)){bs=sm;bm=mx;t={i,j};}
 }
 return t;
}
static I family_margin(const B&ns){
 int L=ns.size(),F=(1<<L)-1;auto m=inv_table(ns);auto tp=top_pair(ns,L-1);
 if(tp.i<0)abort();
 int J=F^(1<<tp.i)^(1<<tp.j);
 I P=4*m[F]-omega(m,F);
 for(int s=1;s<F;s++){int k=pop((unsigned)s);if(k==2||k==L-2)P+=2*m[s]*m[F^s];}
 return P-omega(m,J);
}
struct Profile{vector<int>o;int q,low,L;};
static vector<Profile> make_profiles(){
 set<pair<vector<int>,int>>seen;vector<Profile>out;
 for(int L=6;L<=12;L++)for(int q:{1,2})for(int tau=0;tau<=1;tau++)
 for(int j=-1;j<L;j++)for(int eps:{-2,2}){
  if(j==-1&&eps==2)continue;
  vector<int>o;for(int i=0;i<L;i++)o.push_back(tau+q*i+(i==j?eps:0));
  sort(o.begin(),o.end());int sm=accumulate(o.begin(),o.end(),0);if(sm&1)continue;
  int lo=q==1?20:0;while(2*lo+o[0]<1)lo++;
  if(seen.insert({o,lo}).second)out.push_back({o,q,lo,L});
 }
 return out;
}
static I profile_margin(const Profile&p,int t){B n=p.o;for(int&x:n)x+=2*t;return family_margin(n);}
struct FamCounters{long long finite=0,newton=0,profiles=0,q2zero=0;I minfinite=0,mincoef=0;bool init=true;};
extern "C" int verify_families(int threads){
 omp_set_num_threads(threads);auto prof=make_profiles();
 atomic<long long>done{0},finite{0},newton{0},bad{0},q2zero{0};
 atomic<long long>minf{((long long)1<<60)},minc{((long long)1<<60)};
 printf("MECH158 profiles generated=%zu\n",prof.size());fflush(stdout);
 #pragma omp parallel for schedule(dynamic,1)
 for(int z=0;z<(int)prof.size();z++){
  const auto&p=prof[z];bool ok=true;
  long long stab=2*p.L;for(int x:p.o)stab+=abs(x);
  if(stab>=192)ok=false;
  if(p.q==2&&p.low>0)q2zero++;
  for(int t=p.low;t<96;t++){
   I v=profile_margin(p,t);finite++;
   if(v<0)ok=false;
   long long vll;
   if(v>=LLONG_MIN&&v<=LLONG_MAX){vll=v.convert_to<long long>();long long old=minf.load();while(vll<old&&!minf.compare_exchange_weak(old,vll)){}}
  }
  vector<I>v;for(int t=96;t<=96+p.L-2;t++)v.push_back(profile_margin(p,t));
  vector<I>coef;int h=0;
  while(!v.empty()){
   I a=v[0];coef.push_back(a);newton++;
   if(a<0)ok=false;
   if(a>=LLONG_MIN&&a<=LLONG_MAX){long long x=a.convert_to<long long>(),old=minc.load();while(x<old&&!minc.compare_exchange_weak(old,x)){}}
   vector<I>d;for(size_t k=1;k<v.size();k++)d.push_back(v[k]-v[k-1]);v.swap(d);h++;
  }
  if((int)coef.size()!=p.L-1)ok=false;
  for(int u:{17,73,1000}){
   I pred=0,bin=1;
   for(int h=0;h<(int)coef.size();h++){if(h)bin=bin*(u-h+1)/h;pred+=coef[h]*bin;}
   I actual=profile_margin(p,96+u);if(pred!=actual)ok=false;
  }
  if(!ok)bad++;
  long long d=++done;
  if(d%10==0){printf("MECH158 profiles=%lld/%zu q1/q2 certificate stage\n",d,prof.size());fflush(stdout);}
 }
 printf("MECH158 SUMMARY profiles=%lld bad=%lld finite_values=%lld Newton_coefficients=%lld q2_positive_lower_bound=%lld minfinite=%lld minNewton=%lld\n",
 done.load(),bad.load(),finite.load(),newton.load(),q2zero.load(),minf.load(),minc.load());fflush(stdout);
 return bad?1:0;
}

// Sparse two-variable SU(2) character-basis coefficient evaluator.
static uint64_t key2(int a,int b){return (uint64_t(uint32_t(a))<<32)|uint32_t(b);}
static int keya(uint64_t k){return int(k>>32);}
static int keyb(uint64_t k){return int(k&0xffffffffu);}
static I coeff2(B factors,int X,int Y){
 sort(factors.begin(),factors.end(),[](int a,int b){return abs(a)>abs(b);});
 int rem=0;for(int z:factors)rem+=abs(z);
 unordered_map<uint64_t,I>d,nx;d.reserve(256);d[key2(0,0)]=1;
 for(int z:factors){
  int n=abs(z),eps=z<0?-1:1;rem-=n;nx.clear();nx.reserve(d.size()*2+16);
  for(const auto&kv:d){
   int a=keya(kv.first),b=keyb(kv.first);const I&v=kv.second;
   for(int c=abs(a-n);c<=a+n;c+=2)if(abs(c-X)+abs(b-Y)<=rem)nx[key2(c,b)]+=v;
   for(int c=abs(b-n);c<=b+n;c+=2)if(abs(a-X)+abs(c-Y)<=rem)nx[key2(a,c)]+=eps*v;
  }
  for(auto it=nx.begin();it!=nx.end();)if(it->second==0)it=nx.erase(it);else ++it;
  d.swap(nx);
 }
 auto it=d.find(key2(X,Y));return it==d.end()?I(0):it->second;
}
static I gval(const B&B,int p){return coeff2(B,p,0);}
static I phi_direct(const B&L){return coeff2(L,0,0);}
static B apply_signs(const B&n,const vector<int>&neg){
 B r=n;for(size_t i=0;i<r.size();i++)if(neg[i])r[i]=-r[i];return r;
}
static I abs_small(I x){return x<0?-x:x;}
static void random_signs(const B&n,mt19937_64&rng,vector<int>&neg){
 int N=n.size();neg.assign(N,0);map<int,vector<int>>cls;
 for(int i=0;i<N;i++)cls[n[i]].push_back(i);
 for(auto&[v,idx]:cls)if(rng()&1)for(int i:idx)neg[i]=1;
 int parity=0,fix=-1;for(int i=0;i<N;i++){parity^=neg[i];}
 if(parity){for(auto&[v,idx]:cls)if(idx.size()%2){fix=idx[0];break;}if(fix<0)abort();for(int i:cls[n[fix]])neg[i]^=1;}
}

// Random direct checks for Proposition 1 and Theorem 2.
extern "C" int verify_direct(int seed,int count){
 mt19937_64 rng(seed);int bridge=0,family=0,triggers=0;
 for(int z=0;z<count;z++){
  int L=6+(rng()%7),q=(rng()&1)?1:2,t=100+(rng()%151),tau=rng()%2;
  B n;for(int i=0;i<L;i++)n.push_back(2*t+tau+q*i);
  if(rng()&1){int j=rng()%L;n[j]+=(rng()&1)?2:-2;}
  sort(n.begin(),n.end());
  vector<int>neg;random_signs(n,rng,neg);
  int p=n.back();B Bsigned(n.begin(),n.end()-1);vector<int>Bn(neg.begin(),neg.end()-1);
  B Bsign=apply_signs(Bsigned,Bn);
  TopPair tp=top_pair(Bsigned,L-1);B child=Bsign;
  vector<int>rr{tp.i,tp.j};sort(rr.rbegin(),rr.rend());for(int i:rr)child.erase(child.begin()+i);
  I gp=gval(Bsign,p),gc=gval(child,p);
  if(gp<0||gp<gc){printf("MECH158_DIRECT_FAIL family z=%d L=%d t=%d q=%d gp=%s gc=%s\n",z,L,t,q,sI(gp).c_str(),sI(gc).c_str());return 1;}
  // Cross-check 2g against the exact signed subset multiplicity formula.
  vector<I>m=inv_table(n);int F=(1<<L)-1,mask=0;for(int i=0;i<L;i++)if(neg[i])mask|=1<<i;
  I ph=0;for(int s=F;;s=(s-1)&F){I v=m[s]*m[F^s];ph+=(pop((unsigned)(s&mask))&1)?-v:v;if(!s)break;}
  if(ph!=2*gp){printf("MECH158_DIRECT_BRIDGE_FAIL z=%d\n",z);return 1;}family++;
  // Proposition 1 block identities.
  B C;for(int i=0;i<L-1;i++)if(i!=tp.i&&i!=tp.j)C.push_back(Bsign[i]);
  int a=n[tp.i],b=n[tp.j],ea=Bsign[tp.i]>0?1:-1,eb=Bsign[tp.j]>0?1:-1;
  int psign=(signpar(Bsign)+(ea<0)+(eb<0))&1;
  int ep=psign?-1:1;
  auto G=[&](int x,int y){return coeff2(C,x,y);};
  B Cu;for(int x:C)Cu.push_back(abs(x));
  auto H=[&](int x,int y){return coeff2(Cu,x,y);};
  I T=0,X=0,Y=0,Z=0,K=0;
  for(int c:cg(a,b)){for(int s:cg(p,c))T+=G(s,0);Z+=ea*eb*G(p,c);K+=H(p,c);}
  for(int s:cg(p,a)){X+=eb*G(s,b);K+=H(s,b);}
  for(int s:cg(p,b)){Y+=ea*G(s,a);K+=H(s,a);}
  I gp2=gval(Bsign,p);
  if(gp2!=T+X+Y+Z){printf("PROP1 expansion mismatch trial=%d\n",z);return 1;}
  B Lambda=Bsign;Lambda.push_back(ep*p);I phi=phi_direct(Lambda);
  if(phi!=2*gp2){printf("PROP1 phi bridge mismatch trial=%d\n",z);return 1;}
  auto Dflip=[&](vector<int>ix){
   B Fv=Lambda;for(int i:ix)Fv[i]=-Fv[i];I d=phi-phi_direct(Fv);if(d%4!=0)abort();return d/4;
  };
  I Dab=Dflip({tp.i,tp.j}),Dap=Dflip({tp.i,L-1}),Dbp=Dflip({tp.j,L-1});
  if(Dab!=X+Y||Dap!=Y+Z||Dbp!=X+Z){printf("PROP1 flip identity mismatch trial=%d D=%s,%s,%s blocks=%s,%s,%s\n",z,sI(Dab).c_str(),sI(Dap).c_str(),sI(Dbp).c_str(),sI(X+Y).c_str(),sI(Y+Z).c_str(),sI(X+Z).c_str());return 1;}
  I cost=(abs_small(Dab)+abs_small(Dap)+abs_small(Dbp))/2;
  if(cost>abs_small(X)+abs_small(Y)+abs_small(Z)||abs_small(X)+abs_small(Y)+abs_small(Z)>K){printf("PROP1 unsigned bound mismatch trial=%d\n",z);return 1;}
  I surplus=T-G(p,0);
  if(surplus>=K){triggers++;if(gp2<gval(C,p)||surplus<cost){printf("PROP1 sufficient condition mismatch trial=%d\n",z);return 1;}}
  bridge++;
 }
 printf("MECH158 direct tests random=%d family_bridge=%d Prop1 identities=%d sufficient_triggers=%d\n",count,family,bridge,triggers);fflush(stdout);return 0;
}


// Dense, range-add two-coordinate evaluator for bounded direct checks.
static I coeff2_dense(B factors,int X,int Y){
 sort(factors.begin(),factors.end(),[](int a,int b){return abs(a)>abs(b);});
 int total=0;for(int z:factors)total+=abs(z);int D=total+1,old=0,rem=total;
 vector<I>A((size_t)D*D),N((size_t)D*D),d0(D+2),d1(D+2);A[0]=1;
 for(int z:factors){
  int n=abs(z),eps=z<0?-1:1;rem-=n;int now=old+n;fill(N.begin(),N.end(),I(0));
  for(int b=0;b<=old;b++){
   fill(d0.begin(),d0.end(),I(0));fill(d1.begin(),d1.end(),I(0));
   for(int a=0;a<=old-b;a++){const I&v=A[(size_t)a*D+b];if(v==0)continue;int lo=abs(a-n),hi=a+n,par=lo&1,l=lo/2,r=hi/2;
    if(par==0){d0[l]+=v;d0[r+1]-=v;}else{d1[l]+=v;d1[r+1]-=v;}
   }
   I run=0;int idx=0;for(int c=0;c<=now-b;c+=2,idx++){run+=d0[idx];if(abs(c-X)+abs(b-Y)<=rem)N[(size_t)c*D+b]+=run;}
   run=0;idx=0;for(int c=1;c<=now-b;c+=2,idx++){run+=d1[idx];if(abs(c-X)+abs(b-Y)<=rem)N[(size_t)c*D+b]+=run;}
  }
  for(int a=0;a<=old;a++){
   fill(d0.begin(),d0.end(),I(0));fill(d1.begin(),d1.end(),I(0));
   for(int b=0;b<=old-a;b++){const I&v=A[(size_t)a*D+b];if(v==0)continue;int lo=abs(b-n),hi=b+n,par=lo&1,l=lo/2,r=hi/2;I w=eps*v;
    if(par==0){d0[l]+=w;d0[r+1]-=w;}else{d1[l]+=w;d1[r+1]-=w;}
   }
   I run=0;int idx=0;for(int c=0;c<=now-a;c+=2,idx++){run+=d0[idx];if(abs(a-X)+abs(c-Y)<=rem)N[(size_t)a*D+c]+=run;}
   run=0;idx=0;for(int c=1;c<=now-a;c+=2,idx++){run+=d1[idx];if(abs(a-X)+abs(c-Y)<=rem)N[(size_t)a*D+c]+=run;}
  }
  A.swap(N);old=now;
 }
 return (X>=0&&Y>=0&&X<D&&Y<D)?A[(size_t)X*D+Y]:I(0);
}
extern "C" int verify_direct_fast(int seed,int count){
 mt19937_64 rng(seed);auto profiles=make_profiles();vector<Profile>small;
 int dense_sparse=0;
 for(int z=0;z<24;z++){
  B test;int k=3+(rng()%6),sum=0;
  for(int j=0;j<k;j++){int v=1+(rng()%8);if(rng()&1)v=-v;test.push_back(v);sum+=abs(v);}
  int X=rng()%(sum+1),Y=rng()%(sum+1);
  I a=coeff2_dense(test,X,Y),b=coeff2(test,X,Y);
  if(a!=b){printf("DENSE_SPARSE_FAIL z=%d X=%d Y=%d dense=%s sparse=%s\\n",z,X,Y,sI(a).c_str(),sI(b).c_str());return 1;}
  dense_sparse++;
 }
 printf("MECH158 dense-vs-sparse checks=%d\n",dense_sparse);fflush(stdout);
 for(auto&p:profiles)if(p.L<=7)small.push_back(p);
 int family=0,prop=0;
 for(int z=0;z<count;z++){
  const auto&p=small[rng()%small.size()];int t=100+(rng()%51);B n=p.o;for(int&x:n)x+=2*t;
  vector<int>neg;random_signs(n,rng,neg);int L=n.size(),dist=n.back();B Bs;vector<int>bn;
  for(int i=0;i<L-1;i++){Bs.push_back(n[i]);bn.push_back(neg[i]);}Bs=apply_signs(Bs,bn);
  TopPair tp=top_pair(n,L-1);B child=Bs;int i=max(tp.i,tp.j),j=min(tp.i,tp.j);child.erase(child.begin()+i);child.erase(child.begin()+j);
  I gp=coeff2_dense(Bs,dist,0),gc=coeff2_dense(child,dist,0);
  if(gp<0||gp<gc){printf("MECH158_FAST_FAIL z=%d L=%d t=%d q=%d gp=%s gc=%s\n",z,L,t,p.q,sI(gp).c_str(),sI(gc).c_str());return 1;}
  auto m=inv_table(n);int F=(1<<L)-1,mask=0;for(int k=0;k<L;k++)if(neg[k])mask|=1<<k;
  I ph=0;for(int s=F;;s=(s-1)&F){I v=m[s]*m[F^s];ph+=(pop((unsigned)(s&mask))&1)?-v:v;if(!s)break;}
  if(ph!=2*gp){printf("MECH158_FAST_BRIDGE_FAIL z=%d\n",z);return 1;}family++;
  printf("MECH158 direct-family %d/%d L=%d q=%d t=%d g=%s child=%s\n",z+1,count,L,p.q,t,sI(gp).c_str(),sI(gc).c_str());fflush(stdout);
 }
 for(int z=0;z<20;z++){
  int L=6+(rng()%5);B n;for(int i=0;i<L-1;i++)n.push_back(2+(rng()%11));
  int p=*max_element(n.begin(),n.end())+2;n.push_back(p);sort(n.begin(),n.end());
  vector<int>neg;random_signs(n,rng,neg);B Lambda=apply_signs(n,neg);
  int dist=abs(Lambda.back());B Bs(Lambda.begin(),Lambda.end()-1);TopPair tp=top_pair(n,L-1);if(tp.i<0)return 1;
  B C;for(int i=0;i<L-1;i++)if(i!=tp.i&&i!=tp.j)C.push_back(Bs[i]);
  int a=n[tp.i],b=n[tp.j],ea=Bs[tp.i]>0?1:-1,eb=Bs[tp.j]>0?1:-1;
  auto G=[&](int x,int y){return coeff2(C,x,y);};B Cu;for(int x:C)Cu.push_back(abs(x));auto H=[&](int x,int y){return coeff2(Cu,x,y);};
  I T=0,X=0,Y=0,Z=0,K=0;for(int c:cg(a,b)){for(int s:cg(dist,c))T+=G(s,0);Z+=ea*eb*G(dist,c);K+=H(dist,c);}
  for(int s:cg(dist,a)){X+=eb*G(s,b);K+=H(s,b);}for(int s:cg(dist,b)){Y+=ea*G(s,a);K+=H(s,a);}
  I g=coeff2(Bs,dist,0);if(g!=T+X+Y+Z){printf("MECH158_PROP1_EXPANSION_FAIL z=%d\n",z);return 1;}
  I phi=phi_direct(Lambda);if(phi!=2*g)return 1;
  auto D=[&](int i,int j){B v=Lambda;v[i]=-v[i];v[j]=-v[j];I d=phi-phi_direct(v);if(d%4!=0)abort();return d/4;};
  I dab=D(tp.i,tp.j),dap=D(tp.i,L-1),dbp=D(tp.j,L-1);
  if(dab!=X+Y||dap!=Y+Z||dbp!=X+Z)return 1;
  I bound=(abs_small(dab)+abs_small(dap)+abs_small(dbp))/2;
  if(bound>abs_small(X)+abs_small(Y)+abs_small(Z)||abs_small(X)+abs_small(Y)+abs_small(Z)>K)return 1;
  I surplus=T-G(dist,0);if(surplus>=K&&(g<coeff2(C,dist,0)||surplus<bound))return 1;
  prop++;
  printf("MECH158 Prop1 identities %d/20 passed\n",prop);fflush(stdout);
 }
 printf("MECH158 DIRECT PASS family=%d Prop1=%d\n",family,prop);fflush(stdout);return 0;
}

// Exact fusion row for Lemma 2.
static vector<I> fusionrow(const B&ns){
 vector<I>d(1);d[0]=1;int oldmax=0;
 for(int n:ns){
  vector<I>e(oldmax+n+1);for(int a=0;a<=oldmax;a++)if(d[a]!=0)
   for(int c=abs(a-n);c<=a+n;c+=2)e[c]+=d[a];
  oldmax+=n;d.swap(e);
 }
 return d;
}
extern "C" int verify_lemma2(int maxlabel){
 long long lists=0,coeffs=0;
 for(int k:{3,4}){B n(k);
  function<void(int,int)>rec=[&](int pos,int low){
   if(pos==k){auto f=fusionrow(n);I m=f[0];
    for(int r=0;2*r<=n[0];r++){if(f[2*r]<(r+1)*m){printf("LEMMA2_FAIL k=%d labels=",k);for(int x:n)printf("%d,",x);printf(" r=%d mu=%s mu0=%s\n",r,sI(f[2*r]).c_str(),sI(m).c_str());abort();}coeffs++;}
    lists++;return;
   }
   for(int x=low;x<=maxlabel;x++){n[pos]=x;rec(pos+1,x);}
  };rec(0,1);
 }
 printf("MECH159 Lemma2 fusion tests labels=1..%d lists=%lld inequalities=%lld\n",maxlabel,lists,coeffs);fflush(stdout);return 0;
}

// Exact combinatorial negative-partition counts and rational constants.
static long long cLL(int n,int k){return chooseLL(n,k);}
static long long negative_subsets(int L,int q,int s){
 long long z=0;for(int j=1;j<=q&&j<=s;j+=2)z+=cLL(q,j)*cLL(L-q,s-j);return z;
}
extern "C" int verify_constants(){
 map<int,pair<Q,Q>> expected{{7,{Q(160,7),Q(394,8281)}},{8,{Q(127,2),Q(7853,82524)}},{9,{Q(4571,36),Q(58909,666468)}}};
 for(int L:{7,8,9}){
  Q best(0);long long unweighted=0;int qbest=-1;
  for(int q=0;q<=L;q+=2){
   Q v(0);long long uc=0;
   for(int s=3;s<=L/2;s++){
    long long n=negative_subsets(L,q,s);if(2*s==L)n/=2;
    long long den=2LL*L*(L-1);
    v+=Q(n*(den+s*(L-s)),den);uc+=n;
   }
   if(v>best){best=v;qbest=q;}
   unweighted=max(unweighted,uc);
  }
  if(best!=expected[L].first){printf("K_FAIL L=%d found=%s expected=%s\n",L,sQ(best).c_str(),sQ(expected[L].first).c_str());return 1;}
  int u=L==7?12:L==8?22:32;
  int T=(u/2+1)*(u/2+2)/2,H=(1<<(L-3))-L+1;
  Q margin=Q(1)-best/T-Q(u+1+H,(long long)(u+1)*(u+1));
  if(margin!=expected[L].second){printf("RHO_FAIL L=%d\n",L);return 1;}
  if(L==8&&unweighted!=56){printf("NNEG_FAIL L8 got=%lld\n",unweighted);return 1;}
  if(L==9&&unweighted!=112){printf("NNEG_FAIL L9 got=%lld\n",unweighted);return 1;}
  printf("MECH159 constants L=%d q*=%d K=%s T=%d H=%d rho=%s max_unweighted=%lld\n",L,qbest,sQ(best).c_str(),T,H,sQ(margin).c_str(),unweighted);
 }
 Q p8=Q(1)-Q(56,66),p9=Q(1)-Q(112,120);
 if(p8!=Q(5,33)||p9!=Q(1,15))return 1;
 return 0;
}

struct ChildPhi{int U;vector<int>positions;vector<I>phi;};
static void fwht(vector<I>&a){for(size_t h=1;h<a.size();h*=2)for(size_t i=0;i<a.size();i+=2*h)for(size_t j=0;j<h;j++){I x=a[i+j],y=a[i+j+h];a[i+j]=x+y;a[i+j+h]=x-y;}}
static vector<I> phi_for_subset(const vector<I>&m,int U,const vector<int>&positions){
 int k=positions.size(),N=1<<k;vector<I>a(N);
 for(int s=0;s<N;s++){int g=0;for(int j=0;j<k;j++)if(s>>j&1)g|=1<<positions[j];a[s]=m[g]*m[U^g];}
 fwht(a);return a;
}
struct BoxConfig{int L,lo,hi;bool descent,positive;int rhoN,rhoD;string name;};
static void multisets_rec(int pos,int low,const BoxConfig&cfg,B&n,vector<B>&out){
 if(pos==cfg.L){out.push_back(n);return;}
 for(int x=low;x<=cfg.hi;x++){n[pos]=x;multisets_rec(pos+1,x,cfg,n,out);}
}
struct BoxCounts{long long multisets=0,signs=0,positive=0,removals=0;I worstDelta=0,worstC=0;bool set=false;};
static int process_multiset(const B&n,const BoxConfig&cfg,BoxCounts&out){
 int L=n.size(),F=(1<<L)-1;auto m=inv_table(n);
 vector<I>base(1<<L);for(int s=0;s<=F;s++)base[s]=m[s]*m[F^s];fwht(base);
 vector<pair<vector<int>,vector<I>>>childTabs;
 if(cfg.descent){
  for(int i=0;i<L-1;i++)for(int j=i+1;j<L-1;j++)if(((n[i]+n[j])&1)==0){
   int U=F^(1<<i)^(1<<j);vector<int>pos;for(int k=0;k<L;k++)if(U>>k&1)pos.push_back(k);
   childTabs.push_back({pos,phi_for_subset(m,U,pos)});
  }
 }
 vector<vector<int>>classes;
 for(int i=0;i<L;){int j=i+1;while(j<L&&n[j]==n[i])j++;vector<int>x;for(int k=i;k<j;k++)x.push_back(k);classes.push_back(x);i=j;}
 int C=classes.size();BoxCounts local;local.multisets=1;
 for(int cm=0;cm<(1<<C);cm++){
  int negcount=0,neg=0;
  for(int k=0;k<C;k++)if(cm>>k&1){negcount+=classes[k].size();for(int i:classes[k])neg|=1<<i;}
  if(negcount&1)continue;
  local.signs++;I ph=base[neg];
  if(cfg.positive){if(ph<0){printf("FM3_FAIL L=%d labels=",L);for(int x:n)printf("%d,",x);printf(" negmask=%d phi=%s\n",neg,sI(ph).c_str());return 1;}local.positive++;}
  if(!cfg.descent)continue;
  I M=m[F];long long den=2LL*L*(L-1);I dsum=0;
  for(int i=0;i<L;i++)for(int j=i+1;j<L;j++){
   I d=ph-base[neg^(1<<i)^(1<<j)];if(d%4!=0){printf("D nonintegral\n");return 1;}dsum+=d/4;
  }
  for(size_t r=0;r<childTabs.size();r++){
   const auto&positions=childTabs[r].first;const auto&tab=childTabs[r].second;int cn=0,localneg=0,pLocal=-1;
   for(int j=0;j<(int)positions.size();j++){int i=positions[j];if(neg>>i&1){localneg|=1<<j;cn++;}if(i==L-1)pLocal=j;}
   if(pLocal<0)abort();if(cn&1)localneg^=1<<pLocal;
   I cph=tab[localneg];if((ph&1)!=0||(cph&1)!=0) {printf("phi not even\n");return 1;}
   I delta=(ph-cph)/2;
   if(delta*cfg.rhoD<I(cfg.rhoN)*M){printf("DELTA_FAIL L=%d labels=",L);for(int x:n)printf("%d,",x);printf(" Ridx=%zu delta=%s M=%s\n",r,sI(delta).c_str(),sI(M).c_str());return 1;}
   I Cscaled=delta*den+dsum;
   if(Cscaled*cfg.rhoD<I(cfg.rhoN)*M*den){printf("C_FAIL L=%d labels=",L);for(int x:n)printf("%d,",x);printf(" Ridx=%zu delta=%s dsum=%s M=%s\n",r,sI(delta).c_str(),sI(dsum).c_str(),sI(M).c_str());return 1;}
   local.removals++;
   if(!local.set||delta<local.worstDelta)local.worstDelta=delta;
   if(!local.set||Cscaled<local.worstC)local.worstC=Cscaled;
   local.set=true;
  }
 }
 #pragma omp critical
 {
  out.multisets+=local.multisets;out.signs+=local.signs;out.positive+=local.positive;out.removals+=local.removals;
  if(local.set){if(!out.set||local.worstDelta<out.worstDelta)out.worstDelta=local.worstDelta;if(!out.set||local.worstC<out.worstC)out.worstC=local.worstC;out.set=true;}
 }
 return 0;
}
static vector<B> all_multisets(const BoxConfig&cfg){
 vector<B>out;B n(cfg.L);multisets_rec(0,cfg.lo,cfg,n,out);return out;
}
static int run_box(const BoxConfig&cfg,int threads){
 auto lists=all_multisets(cfg);BoxCounts total;atomic<int>done{0},fail{0};
 printf("MECH159 BOX START %s multisets=%zu range=%d..%d L=%d\n",cfg.name.c_str(),lists.size(),cfg.lo,cfg.hi,cfg.L);fflush(stdout);
 #pragma omp parallel for schedule(dynamic,8) num_threads(threads)
 for(int k=0;k<(int)lists.size();k++){
  BoxCounts local;
  if(process_multiset(lists[k],cfg,local)){fail++;continue;}
  int d=++done;
  #pragma omp critical
  {
   total.multisets+=local.multisets;total.signs+=local.signs;total.positive+=local.positive;total.removals+=local.removals;
   if(local.set){if(!total.set||local.worstDelta<total.worstDelta)total.worstDelta=local.worstDelta;if(!total.set||local.worstC<total.worstC)total.worstC=local.worstC;total.set=true;}
   if(d%500==0){printf("MECH159 BOX %s multisets=%d/%zu signs=%lld removals=%lld\n",cfg.name.c_str(),d,lists.size(),total.signs,total.removals);fflush(stdout);}
  }
 }
 printf("MECH159 BOX PASS %s multisets=%lld even_signs=%lld positive=%lld removal_checks=%lld failures=%d\n",
 cfg.name.c_str(),total.multisets,total.signs,total.positive,total.removals,fail.load());fflush(stdout);
 return fail?1:0;
}
extern "C" int verify_boxes(int threads){
 vector<BoxConfig>c{
  {7,12,17,true,false,394,8281,"L7 descent"},
  {8,22,27,true,true,7853,82524,"L8 descent+FM3"},
  {9,32,38,true,true,58909,666468,"L9 descent+FM3"},
  {8,20,25,false,true,0,1,"L8 FM3 threshold"},
  {9,28,34,false,true,0,1,"L9 FM3 threshold"}
 };
 for(const auto&x:c)if(run_box(x,threads))return 1;
 return 0;
}
extern "C" int verify_random_159(int seed,int count,int threads){
 mt19937_64 rng(seed);int done=0,direct=0;
 for(int L:{7,8,9})for(int z=0;z<count;z++){
  int lo=L==7?12:L==8?22:32;B n(L);
  for(int&i:n)i=lo+(rng()%120);sort(n.begin(),n.end());
  BoxConfig cfg{L,lo,lo+120,true,true,L==7?394:L==8?7853:58909,L==7?8281:L==8?82524:666468,"random"};
  BoxCounts row; if(process_multiset(n,cfg,row))return 1;done++;
  if(z<3){
   int F=(1<<L)-1;vector<I>m=inv_table(n),a(1<<L);for(int s=0;s<=F;s++)a[s]=m[s]*m[F^s];fwht(a);
   vector<vector<int>>cls;for(int i=0;i<L;){int j=i+1;while(j<L&&n[j]==n[i])j++;vector<int>q;for(int k=i;k<j;k++)q.push_back(k);cls.push_back(q);i=j;}
   vector<int>neg(L);int bits=0,cnt=0;for(int c=0;c<(int)cls.size();c++)if(rng()&1){bits|=1<<c;cnt+=cls[c].size();for(int i:cls[c])neg[i]=1;}
   if(cnt&1){for(int c=0;c<(int)cls.size();c++)if(cls[c].size()%2){for(int i:cls[c])neg[i]^=1;break;}}
   int nm=0;B Bsign;for(int i=0;i<L;i++){if(neg[i])nm|=1<<i;if(i<L-1)Bsign.push_back(neg[i]?-n[i]:n[i]);}
   int p=n.back();I gp=coeff2(Bsign,p,0);if(a[nm]!=2*gp){printf("RANDOM159 2D bridge failed\n");return 1;}
   auto tp=top_pair(n,L-1);if(tp.i<0)return 1;B child=Bsign;int i=max(tp.i,tp.j),j=min(tp.i,tp.j);child.erase(child.begin()+i);child.erase(child.begin()+j);
   I gc=coeff2(child,p,0);if(gp<0||gp<gc){printf("RANDOM159 direct conclusion fail\n");return 1;}direct++;
  }
 }
 printf("MECH159 random lists=%d direct 2D checks=%d\n",done,direct);fflush(stdout);return 0;
}

extern "C" int run_all(int threads){
 omp_set_num_threads(threads);
 if(verify_families(threads))return 1;
 if(verify_direct_fast(172030,8))return 1;
 if(verify_lemma2(24))return 1;
 if(verify_constants())return 1;
 if(verify_boxes(threads))return 1;
 if(verify_random_159(159030,80,threads))return 1;
 puts("ALL CHECKS PASS");fflush(stdout);return 0;
}
int main(int argc,char**argv){
 int threads=12;bool only159=false,only158direct=false;
 for(int i=1;i<argc;i++){
  string a=argv[i];
  if(a=="--help"||a=="-h"){puts("FM-CHK101 verifier: --threads N --159-only --158-direct-only; exact C++ checks for MECH158/159. No files are written.");return 0;}
  if(a=="--threads"&&i+1<argc)threads=atoi(argv[++i]);
  if(a=="--159-only")only159=true;
  if(a=="--158-direct-only")only158direct=true;
 }
 if(threads<1){fprintf(stderr,"--threads must be positive\n");return 2;}
 if(only158direct){return verify_direct_fast(172030,8);}
 if(only159){
  if(verify_lemma2(24)||verify_constants()||verify_boxes(threads)||verify_random_159(159030,80,threads))return 1;
  puts("MECH159 CHECKS PASS");fflush(stdout);return 0;
 }
 return run_all(threads);
}
