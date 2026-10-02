import os, sys, subprocess
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, prod
from fractions import Fraction as Q
from datetime import datetime, timezone
def log(*a):
 print(datetime.now(timezone.utc).isoformat(timespec='seconds'),*a,flush=True)
if any(a in ('-h','--help') for a in sys.argv[1:]):
 print('Usage: python3 -u -B verifier.py; exact STR8d obstruction and modular triangular-minor audit.')
 raise SystemExit
@lru_cache(None)
def paths(ns):
 if not ns:return ((0,()),)
 return tuple((b,p+(b,)) for a,p in paths(ns[:-1])
              for b in range(abs(a-ns[-1]),a+ns[-1]+1,2))
@lru_cache(None)
def tail(w):
 ns=tuple(map(abs,w)); L=len(w); d={}
 for m in range(1<<L):
  x=tuple(ns[i] for i in range(L) if m>>i&1)
  y=tuple(ns[i] for i in range(L) if not m>>i&1)
  ep=sum(w[i]<0 for i in range(L) if not m>>i&1)%2
  for a,pa in paths(x):
   for b,pb in paths(y):
    d.setdefault((a,b),[[],[]])[ep].append((m,pa,pb))
 out={}
 for ab,v in d.items():
  for z in v:z.sort()
  r=min(map(len,v))
  for ep in (0,1):
   if len(v[ep])>r:out[ab]=(ep,v[ep][r:])
 return out
def info(w):
 A,B=tail(w[::2]),tail(w[1::2]); e=o=c=0
 for ab,(ep,u) in A.items():
  if ab not in B:continue
  eq,v=B[ab]; z=len(u)*len(v)
  if ep==eq:e+=z;continue
  o+=z; ba=ab[::-1]
  if ba not in A or ba not in B:continue
  ax=set(A[ba][1]); by=set(B[ba][1])
  ac=sum((((1<<len(w[::2]))-1)^m,pb,pa) in ax for m,pa,pb in u)
  bc=sum((((1<<len(w[1::2]))-1)^m,pb,pa) in by for m,pa,pb in v)
  c+=ac*bc
 return e,o,c
jobs={}; total=corrected=obstructed=0; first=None
for L in range(11):
 local=0
 for ns in combinations_with_replacement(range(1,5),L):
  present=set(ns)
  for m in range(16):
   if any(m>>(n-1)&1 for n in range(1,5) if n not in present):continue
   w=tuple(-n if m>>(n-1)&1 else n for n in ns)
   if sum(n<0 for n in w)%2:continue
   total+=1; e,o,c=info(w)
   if not o:continue
   corrected+=1
   if first is None:first=w
   if c:obstructed+=1;local+=1;continue
   key=min(w,tuple(-n if abs(n)%2 else n for n in w))
   if key not in jobs:jobs[key]=[e,o,0]
   jobs[key][2]+=1
 log('BOX',L,'complement-obstructed',local)
assert (total,corrected,obstructed)==(4521,531,263)
assert first==(-1,-1,-1,-2,3,4)
assert sum(z[2] for z in jobs.values())==268
log('BOX exact',total,corrected,obstructed,'remaining profiles',268)
expected=[(40,4,4),(196,12,12),(1080,100,46),(7256,822,478),(53872,7664,2988)]
for k,want in zip(range(5,10),expected):
 p=max(k,6); W=k*(k+1)//2
 while (W-p)%2 or (not k%2 and p<=k):p+=1
 w=tuple(-i for i in range(1,k+1))+(p if k%2==0 else -p,)
 assert info(w)==want
 log('RUN',k,'E,O,complement-paired columns',want)
w=(-1,2,3,4,-5,6,7,8)
assert info(w)==(1172,216,174)
log('SIGN MINIMUM',info(w))
w=tuple(-n for n in (1,2,4,5) for _ in range(3))
assert info(w)==(60756,1308,10)
log('FOUR TRIPLED CLASSES',info(w))
@lru_cache(None)
def paths_to(ns,t):
 if not ns:return ((),) if t==0 else ()
 if t<0 or t>sum(ns) or (sum(ns)-t)%2:return ()
 old=ns[:-1]; n=ns[-1]
 return tuple(sorted(p+(t,) for u in range(abs(t-n),min(t+n,sum(old))+1,2)
                     for p in paths_to(old,u)))
def selected(ns,ab):
 vals=[[],[]]; L=len(ns)
 for m in range(1<<L):
  x=tuple(n for i,n in enumerate(ns) if m>>i&1)
  y=tuple(n for i,n in enumerate(ns) if not m>>i&1)
  for px in paths_to(x,ab[0]):
   for py in paths_to(y,ab[1]):vals[(L-m.bit_count())%2].append((m,px,py))
 for z in vals:z.sort()
 r=min(map(len,vals))
 return tuple(map(len,vals)),[z[r:] for z in vals]
for a in range(4,129):
 A=(1,1,2,a,a,a+1); B=(1,2,2,a,a+1,a+1); alpha=(2*a+4,a+1)
 ca,sa=selected(A,alpha); cb,sb=selected(B,alpha)
 assert ca==(4,1) and cb==(7,14 if a==4 else 12)
 u=(46,(1,3,a+3,2*a+4),(1,a+1))
 v=(47,(1,3,5,a+5,2*a+4),(a+1,))
 assert u in sa[0] and v in sb[1]
 ur=(63^u[0],u[2],u[1]); vr=(63^v[0],v[2],v[1])
 _,ta=selected(A,alpha[::-1]); _,tb=selected(B,alpha[::-1])
 assert ur in ta[0] and vr in tb[1]
 paths_to.cache_clear()
log('PASS 125 exact checks of the uniform complementary sources')
def clean(p):return {v:c for v,c in p.items() if c}
def lower(p,ns):
 q={}
 for v,c in p.items():
  for i,n in enumerate(ns):
   if v[i]<n:
    w=list(v);w[i]+=1;w=tuple(w);q[w]=q.get(w,0)+c*(n-v[i])
 return clean(q)
@lru_cache(None)
def cv(ns):
 if not ns:return {():({():Q(1)},)}
 out={}; n=ns[-1]
 for path,states in cv(ns[:-1]).items():
  a=path[-1] if path else 0
  for b in range(abs(a-n),a+n+1,2):
   t=(a+n-b)//2; hv={}
   for h in range(t+1):
    for v,z in states[h].items():
     w=v+(t-h,);hv[w]=hv.get(w,0)+(-1)**h*comb(t,h)*z
   st=[clean(hv)]
   for h in range(b):
    st.append({v:z/Q(b-h) for v,z in lower(st[-1],ns).items()})
   out[path+(b,)]=tuple(st)
 return out
def embed(p,ids,L):
 q={}
 for v,c in p.items():
  w=[0]*L
  for i,h in zip(ids,v):w[i]=h
  q[tuple(w)]=c
 return q
def mul(p,q):
 out={}
 for v,c in p.items():
  for w,d in q.items():
   z=tuple(a+b for a,b in zip(v,w));out[z]=out.get(z,0)+c*d
 return clean(out)
def invpair(ns,I,J,pa,pb,j):
 a=cv(tuple(ns[i] for i in I))[pa];b=cv(tuple(ns[i] for i in J))[pb];out={}
 for h in range(j+1):
  z=mul(embed(a[h],I,len(ns)),embed(b[j-h],J,len(ns)))
  for v,c in z.items():out[v]=out.get(v,0)+(-1)**h*comb(j,h)*c
 return clean(out)
def records(w,par):
 A,B=tail(w[::2]),tail(w[1::2])
 return [(ab,u,v) for ab in sorted(A.keys()&B.keys())
         if (A[ab][0]^B[ab][0])==par for u in A[ab][1] for v in B[ab][1]]
def hv(w,record):
 ns=tuple(map(abs,w));ab,u,v=record
 A=range(0,len(w),2);B=range(1,len(w),2)
 ax=tuple(i for j,i in enumerate(A) if u[0]>>j&1)
 ay=tuple(i for j,i in enumerate(A) if not u[0]>>j&1)
 bx=tuple(i for j,i in enumerate(B) if v[0]>>j&1)
 by=tuple(i for j,i in enumerate(B) if not v[0]>>j&1)
 return sum(1<<i for i in ay+by),mul(invpair(ns,ax,bx,u[1],v[1],ab[0]),
                                      invpair(ns,ay,by,u[2],v[2],ab[1]))
def inner(p,q,ns):
 return sum((c*q.get(v,0)/prod(comb(n,h) for n,h in zip(ns,v))
             for v,c in p.items()),Q(0))
w=(-1,-1,-1,-2,3,4);ns=tuple(map(abs,w))
od=[hv(w,z) for z in records(w,1)]
ev=[hv(w,records(w,0)[i]) for i in (0,9)]
assert [s for s,_ in od]==[50,13] and od[0][1]==od[1][1]
assert [s for s,_ in ev]==[63,0] and ev[0][1]==ev[1][1]
M=[[inner(t,s,ns)*prod(i+2 for i in range(len(w)) if (S&T)>>i&1)
    for S,s in od] for T,t in ev]
assert M==[[630,200],[5,5]] and M[0][0]*M[1][1]-M[0][1]*M[1][0]==2150
log('PASS exact first obstruction and dense minor',M)

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
struct Group{AB ab;vector<int>a,b;int offset;};
struct Job{vector<int>w;ll he,ho;int multiplicity;};
struct MatrixData{vector<vector<int>>m;vector<tuple<int,int,int>>cols,rows;};
MatrixData matrix_of(const vector<int>&w){
 vector<int>na,nb,pa,pb,wa,wb;
 for(int i=0;i<(int)w.size();i++)if(i%2){nb.push_back(abs(w[i]));pb.push_back(i);wb.push_back(w[i]);}else{na.push_back(abs(w[i]));pa.push_back(i);wa.push_back(w[i]);}
 Half A(na,pa,true),B(nb,pb,false);auto aa=A.surplus(wa),bb=B.surplus(wb);
 vector<Group>src,dst;int ne=0,no=0;
 for(auto&[ab,x]:aa)if(bb.count(ab)){
  auto&y=bb.at(ab);int p=x.first^y.first;auto&vec=p?src:dst;int&ct=p?no:ne;
  vec.push_back({ab,x.second,y.second,ct});ct+=x.second.size()*y.second.size();
 }
 MatrixData ans;
 for(int g=0;g<(int)src.size();g++)for(int sa:src[g].a)for(int sb:src[g].b)ans.cols.push_back({g,sa,sb});
 for(int g=0;g<(int)dst.size();g++)for(int ta:dst[g].a)for(int tb:dst[g].b){
  ans.rows.push_back({g,ta,tb});auto&target=dst[g];vector<int>row(no);
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
  ans.m.push_back(move(row));
 }
 return ans;
}
int peel(const vector<vector<int>>&m,bool pairs){
 int E=m.size(),O=m[0].size(),done=0;
 vector<int>left(O,1),deg(E,0);vector<vector<int>>adj(O);set<int>ones;set<pair<int,int>>q;
 for(int r=0;r<E;r++)for(int c=0;c<O;c++)if(m[r][c]){deg[r]++;adj[c].push_back(r);}
 for(int r=0;r<E;r++)if(deg[r]==1)ones.insert(r);
 auto remove=[&](int c){
  assert(left[c]);left[c]=0;done++;
  for(int r:adj[c]){assert(deg[r]>0);deg[r]--;if(deg[r]==1)ones.insert(r);}
 };
 while(done<O){
  while(!ones.empty()){
   int r=*ones.begin();ones.erase(ones.begin());if(deg[r]!=1)continue;
   int c=0;while(!left[c]||!m[r][c])c++;remove(c);
  }
  if(!pairs||done==O)break;
  map<pair<int,int>,int>first;bool found=false;
  for(int r=0;r<E&&!found;r++)if(deg[r]==2){
   int a=0;while(!left[a]||!m[r][a])a++;
   int b=a+1;while(!left[b]||!m[r][b])b++;
   auto key=make_pair(a,b);auto it=first.find(key);
   if(it==first.end()){first[key]=r;continue;}
   int t=it->second;
   if(mm(m[r][a],m[t][b])!=mm(m[r][b],m[t][a])){
    remove(a);remove(b);found=true;
   }
  }
  if(!found)break;
 }
 return done;
}
int main(int argc,char**argv){
 if(argc>1&&(string(argv[1])=="--help"||string(argv[1])=="-h")){
  cout<<"Usage: STR8d modular triangular audit; jobs are read from stdin.\n";return 0;
 }
 for(int d=2;ll(d)*d<=MOD;d++)assert(MOD%d);
 int J;cin>>J;vector<Job>jobs(J);
 for(auto&job:jobs){
  int L;cin>>L>>job.he>>job.ho>>job.multiplicity;job.w.resize(L);
  for(int&v:job.w)cin>>v;
 }
 omp_set_num_threads(4);int complete=0,singles=0,blocks=0;
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<J;j++){
  auto&job=jobs[j];auto M=matrix_of(job.w);
  assert(M.m.size()==size_t(job.he)&&M.cols.size()==size_t(job.ho));
  int a=peel(M.m,false),b=peel(M.m,true);
  #pragma omp critical(results)
  {
   complete++;
   if(a==job.ho)singles+=job.multiplicity;
   if(b==job.ho)blocks+=job.multiplicity;
   if(complete%24==0||complete==J)logline("MODULAR "+to_string(complete)+"/"+to_string(J));
  }
 }
 assert(singles==165&&blocks==236);
 logline("PASS singleton="+to_string(singles)+" block2="+to_string(blocks)+" of 268 profiles");
 vector<int>w;for(int i=1;i<=4;i++)for(int j=0;j<3;j++)w.push_back(-i);
 auto M=matrix_of(w);assert(M.m.size()==16650&&M.cols.size()==224);
 assert(peel(M.m,false)==10&&peel(M.m,true)==29);
 logline("PASS four tripled classes (1,2,3,4): peeled 10 / 29 of 224");
}
'''
log('compile in memory')
fd=os.memfd_create('str8d-verifier',0)
env=dict(os.environ,TMPDIR='/dev/shm')
data=[str(len(jobs))]
for w,(e,o,m) in sorted(jobs.items()):
 data.append(' '.join(map(str,(len(w),e,o,m,*w))))
try:
 subprocess.run(['g++','-O3','-std=c++17','-fopenmp','-x','c++',
                 '-o',f'/proc/self/fd/{fd}','-'],
                input=src.encode(),pass_fds=(fd,),env=env,check=True)
 subprocess.run([f'/proc/self/fd/{fd}'],input=('\n'.join(data)+'\n').encode(),
                pass_fds=(fd,),check=True)
finally:os.close(fd)
log('PASS completed STR8d verifier')
