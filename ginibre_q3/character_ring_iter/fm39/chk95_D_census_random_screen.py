import os, subprocess, sys

cpp = r'''#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using I=__int128_t; using U=__uint128_t;
struct Rem { int x=0,y=0; };
template<int C> struct Profile { array<signed char,C+1> c{}; U code=0; int w=0,mx=0,core=0,nf=0; array<I,C/2+1> g{}; };
template<int C> struct Work { I a[C+1][C+1]{},b[C+1][C+1]{},dx[2][C+3][C+1]{},dy[2][C+3][C+1]{}; };

template<int C> array<I,C/2+1> eval(const Profile<C>&r) {
 static thread_local Work<C> z;
 I(*a)[C+1]=z.a; I(*b)[C+1]=z.b; int s=0; a[0][0]=1;
 for(int n=1;n<=C;++n) for(int k=0;k<abs((int)r.c[n]);++k) {
  int t=s+n,e=r.c[n]>0?1:-1;
  for(int i=0;i<=t;++i) for(int j=0;j<=t;++j) b[i][j]=0;
  for(int q=0;q<2;++q) for(int x=q;x<=t+2;x+=2) for(int j=0;j<=s;++j) z.dx[q][x][j]=0;
  for(int q=0;q<2;++q) for(int y=q;y<=t+2;y+=2) for(int i=0;i<=s;++i) z.dy[q][y][i]=0;
  for(int i=0;i<=s;++i) for(int j=0;j<=s;++j) {
   I v=a[i][j]; if(!v) continue;
   int lo=abs(i-n),hi=i+n,q=(i+n)&1;
   z.dx[q][lo][j]+=v; z.dx[q][hi+2][j]-=v;
   lo=abs(j-n); hi=j+n; q=(j+n)&1;
   z.dy[q][lo][i]+=e*v; z.dy[q][hi+2][i]-=e*v;
  }
  for(int j=0;j<=s;++j) for(int q=0;q<2;++q) {
   I v=0; for(int x=q;x<=t;x+=2) { v+=z.dx[q][x][j]; b[x][j]+=v; }
  }
  for(int i=0;i<=s;++i) for(int q=0;q<2;++q) {
   I v=0; for(int y=q;y<=t;y+=2) { v+=z.dy[q][y][i]; b[i][y]+=v; }
  }
  swap(a,b); s=t;
 }
 array<I,C/2+1> out{};
 for(int p=r.w&1;p<=r.w;p+=2) out[p/2]=a[p][0];
 return out;
}

template<int C> map<pair<int,int>,I> slow(const Profile<C>&r) {
 map<pair<int,int>,I>a,b; a[{0,0}]=1;
 for(int n=1;n<=C;++n) for(int k=0;k<abs((int)r.c[n]);++k) {
  b.clear(); int e=r.c[n]>0?1:-1;
  for(auto [ab,v]:a) {
   auto[x,y]=ab;
   for(int t=abs(x-n);t<=x+n;t+=2) b[{t,y}]+=v;
   for(int t=abs(y-n);t<=y+n;t+=2) b[{x,t}]+=e*v;
  }
  a.swap(b);
 }
 return a;
}

template<int C> Profile<C> child(Profile<C> r,Rem q) {
 auto one=[&](int n) { r.c[n]+=(r.c[n]>0?-1:1); r.w-=n; r.nf--; };
 one(q.x); if(q.y) one(q.y);
 r.mx=r.core=0;
 for(int n=1;n<=C;++n) if(r.c[n]) { r.mx=n; if(n>=3) r.core+=abs((int)r.c[n]); }
 r.g.fill(0); return r;
}

template<int C> vector<Rem> removals(const Profile<C>&r) {
 vector<Rem>v;
 for(int n=2;n<=C;n+=2) if(r.c[n]) v.push_back({n,0});
 for(int n=1;n<=C;++n) if(r.c[n]) for(int m=n;m<=C;++m)
  if(r.c[m]&&(n&1)==(m&1)&&(m!=n||abs((int)r.c[n])>=2)) v.push_back({n,m});
 return v;
}
template<int C> int smallestEven(const Profile<C>&r) {
 for(int n=2;n<=C;n+=2) if(r.c[n]) return n;
 return 0;
}
template<int C> Rem oneClassRule(const Profile<C>&r,int n) {
 if(abs((int)r.c[n])>=2) return {n,n};
 int other=0;
 for(int m=1;m<=C;++m) if(m!=n&&r.c[m]&&(m&1)==(n&1)) other=max(other,m);
 if(other) return {min(n,other),max(n,other)};
 if(!(n&1)) return {n,0};
 int e=smallestEven(r); return e?Rem{e,0}:Rem{};
}
template<int C> Rem ruleW(const Profile<C>&r) {
 int pick=0,bw=-1;
 for(int n=1;n<=C;++n) if(r.c[n]) {
  int q=n*abs((int)r.c[n]);
  if(q>bw||(q==bw&&n>pick)) pick=n,bw=q;
 }
 return oneClassRule(r,pick);
}
template<int C> Rem ruleM(const Profile<C>&r) {
 int n=0; for(int j=1;j<=C;++j) if(r.c[j]) n=j;
 return oneClassRule(r,n);
}
template<int C> Rem ruleF(const Profile<C>&r) {
 int n=0,kmax=0;
 for(int j=1;j<=C;++j) if(r.c[j]) {
  int k=abs((int)r.c[j]);
  if(k>kmax||(k==kmax&&j<n)) n=j,kmax=k;
 }
 if(kmax>=2) return {n,n};
 Rem best{}; int sum=-1;
 for(int a=1;a<=C;++a) if(r.c[a]) for(int b=a+1;b<=C;++b)
  if(r.c[b]&&(a&1)==(b&1)) {
   int z=a+b;
   if(z>sum||(z==sum&&(b>best.y||(b==best.y&&a>best.x)))) best={a,b},sum=z;
  }
 if(sum>=0) return best;
 int e=smallestEven(r); return e?Rem{e,0}:Rem{};
}
template<int C> int lowp(const Profile<C>&r) {
 int p=max(6,r.mx); if((p&1)!=(r.w&1)) ++p; return p;
}
template<int C> int highp(const Profile<C>&r) { return min(r.w-16,r.w-2*r.mx); }
template<int C> bool boundok(const Profile<C>&r,int&bits) {
 U z=1,cap=U(1)<<120;
 for(int n=1;n<=C;++n) for(int j=0;j<abs((int)r.c[n]);++j) {
  U f=U(2*(n+1)); if(z>=(cap+f-1)/f) return false; z*=f;
 }
 bits=0; for(U x=z;x;x>>=1) ++bits; return true;
}
string ustr(U x) { if(!x)return "0"; string s; while(x){s.push_back(char('0'+x%10));x/=10;} reverse(s.begin(),s.end());return s; }
string istr(I x) { return x<0?"-"+ustr(U(-(x+1))+1):ustr(U(x)); }
U gcdU(U a,U b) { while(b){U t=a%b;a=b;b=t;}return a; }
bool lessFrac(U a,U b,U c,U d) {
 bool flip=false;
 for(;;) {
  U x=a/b,y=c/d; if(x!=y) return flip?x>y:x<y;
  U r=a%b,s=c%d;
  if(!r||!s) { if(!r&&!s)return false; bool firstLess=!r; return flip?!firstLess:firstLess; }
  a=b;b=r;c=d;d=s;flip=!flip;
 }
}
template<int C> string showB(const Profile<C>&r) {
 string s="("; bool f=true;
 for(int n=1;n<=C;++n) if(r.c[n]) {
  if(!f)s+=","; f=false;
  s+=(r.c[n]>0?"+":"-")+to_string(n)+"^"+to_string(abs((int)r.c[n]));
 }
 return s+")";
}
template<int C> string keyR(Rem r) { return to_string(r.x)+":"+to_string(r.y); }

template<int C> void selftest() {
 vector<pair<int,int>>v;
 for(int j=0;j<22;++j)v.push_back({1,-1});
 for(int j=0;j<2;++j)v.push_back({2,-1});
 Profile<C>a; for(auto[n,e]:v)a.c[n]+=e,a.w+=n;
 auto f=eval<C>(a); auto s=slow(a);
 if(f[3]!=9557998450LL||s[{6,0}]!=f[3]) throw runtime_error("Rule-M direct witness");
 Profile<C>ac=child(a,{2,2}); auto fc=eval<C>(ac);
 if(fc[3]!=12194599494LL||f[3]-fc[3]>=0) throw runtime_error("Rule-M increase witness");
 v.clear();
 for(int j=0;j<8;++j)v.push_back({1,-1});
 for(int j=0;j<4;++j)v.push_back({2,-1});
 for(int j=0;j<2;++j)v.push_back({3,-1});
 Profile<C>b; for(auto[n,e]:v)b.c[n]+=e,b.w+=n;
 auto fb=eval<C>(b); Profile<C>bc=child(b,{1,3}); auto fbc=eval<C>(bc);
 if(fb[3]!=51874||fbc[3]!=6976) throw runtime_error("SEC164 witness");
 mt19937_64 g(1081);
 for(int t=0;t<100;++t) {
  Profile<C>r; int w=0;
  for(int j=0;j<6&&w<15;++j) {
   int n=1+g()%6; if(w+n>15)break;
   int e=(g()&1)?1:-1;
   if(r.c[n]&&((r.c[n]>0)!=(e>0)))continue;
   r.c[n]+=e;r.w+=n;w+=n;
  }
  auto fast=eval<C>(r); auto direct=slow(r);
  for(int p=0;p<=r.w;++p) {
   I got=((r.w-p)&1)?0:fast[p/2];
   I want=direct.count({p,0})?direct[{p,0}]:0;
   if(got!=want)throw runtime_error("direct sparse comparison");
  }
 }
 if constexpr(C>=120) {
  Profile<C>q; q.c[1]=-2;q.c[2]=2;
  for(int n=6;n<=10;++n)q.c[n]=2;
  q.c[11]=1;q.w=97;q.mx=11;q.nf=15;
  Rem rr=ruleF(q); auto fq=eval<C>(q); auto dq=slow(q);
  auto qc=child(q,rr); auto fqc=eval<C>(qc);
  if(rr.x!=1||rr.y!=1||fq[5]!=4589803663LL||fqc[5]!=6291848748LL||dq[{11,0}]!=fq[5])
   throw runtime_error("frequency-rule exact kill witness");
  cout<<"frequency-rule control: parent="<<istr(fq[5])<<" child="<<istr(fqc[5])
      <<" drop="<<istr(fq[5]-fqc[5])<<"\n";
 }
 cout<<"direct_selftests=PASS (100 small profiles; Rule-M and SEC164 witnesses)\n";
}

struct HashU {
 size_t operator()(U x)const noexcept {
  uint64_t a=(uint64_t)x,b=(uint64_t)(x>>64);
  return size_t(a^(b+0x9e3779b97f4a7c15ULL+(a<<6)+(a>>2)));
 }
};
static vector<Profile<36>> v36;
static array<signed char,37> cur36{};
static array<U,37> place36{};
void gen36(int n,int w,int mx,int core,int nf,U code) {
 if(n>36) {
  Profile<36>r;r.c=cur36;r.code=code;r.w=w;r.mx=mx;r.core=core;r.nf=nf;v36.push_back(r);return;
 }
 int cap=36/n,take=min(cap,(36-w)/n);
 for(int z=-take;z<=take;++z) {
  cur36[n]=z;
  gen36(n+1,w+n*abs(z),z?n:mx,core+(n>=3?abs(z):0),nf+abs(z),code+U(z+cap)*place36[n]);
 }
 cur36[n]=0;
}
void exhaustive() {
 U mult=1;for(int n=1;n<=36;++n){place36[n]=mult;mult*=U(2*(36/n)+1);}
 gen36(1,0,0,0,0,0);
 sort(v36.begin(),v36.end(),[](auto&a,auto&b){return a.w!=b.w?a.w<b.w:a.code<b.code;});
 cout<<"W<=36 profiles="<<v36.size()<<"\n";
 unordered_map<U,size_t,HashU>ix;ix.reserve(v36.size()*1.15);ix.max_load_factor(.72f);
 for(size_t i=0;i<v36.size();++i)ix.emplace(v36[i].code,i);
 #pragma omp parallel for schedule(dynamic,64) num_threads(32)
 for(long long i=0;i<(long long)v36.size();++i)v36[i].g=eval<36>(v36[i]);
 unsigned long long rb=0,pairs=0,comp=0,fail=0,pos=0,zero=0,neg=0;
 bool have=false;U bn=0,bd=1;size_t bi=0;int bp=0;Rem br{};I bpar=0,bchild=0,bdelta=0;
 for(size_t i=0;i<v36.size();++i) {
  const auto&r=v36[i];if(r.core<2)continue;
  int p0=lowp(r),p1=highp(r);if(p0>p1)continue;++rb;
  auto rs=removals(r);
  for(int p=p0;p<=p1;p+=2) {
   ++pairs;I parent=r.g[p/2];
   if(parent>0)++pos;else if(parent<0)++neg;else ++zero;
   bool any=false;I best=0,bc=0;Rem rr{};
   for(Rem q:rs) {
    ++comp;U code=r.code;
    auto one=[&](int n){if(r.c[n]>0)code-=place36[n];else code+=place36[n];};
    if(!q.y)one(q.x);
    else if(q.x==q.y){if(r.c[q.x]>0)code-=2*place36[q.x];else code+=2*place36[q.x];}
    else{one(q.x);one(q.y);}
    auto it=ix.find(code);if(it==ix.end())throw runtime_error("missing child profile");
    const auto&c=v36[it->second];I val=p>c.w?0:c.g[p/2];I d=parent-val;
    if(!any||d>best)any=true,best=d,bc=val,rr=q;
   }
   if(!any||best<0)++fail;
   if(parent>0&&any) {
    U num=best>=0?U(best):U(-(best+1))+1,den=U(parent);
    if(!have||lessFrac(num,den,bn,bd)){have=true;bn=num;bd=den;bi=i;bp=p;br=rr;bpar=parent;bchild=bc;bdelta=best;}
   }
  }
 }
 cout<<"residual_backgrounds="<<rb<<" residual_(B,p)_pairs="<<pairs
     <<" all_removal_comparisons="<<comp<<" D_failures="<<fail
     <<" parent(+/0/-)="<<pos<<"/"<<zero<<"/"<<neg<<"\n";
 if(have) {
  U g=gcdU(bn,bd);
  cout<<"minimum_max_R_margin="<<ustr(bn/g)<<"/"<<ustr(bd/g)
      <<" B="<<showB(v36[bi])<<" p="<<bp<<" R="
      <<(v36[bi].c[br.x]>0?"+":"-")<<br.x;
  if(br.y)cout<<","<<(v36[bi].c[br.y]>0?"+":"-")<<br.y;
  cout<<" parent="<<istr(bpar)<<" child="<<istr(bchild)<<" drop="<<istr(bdelta)<<"\n";
 }
 if(v36.size()!=2298873||pairs!=6590648||fail)throw runtime_error("exhaustive totals");
}

void randomScreen() {
 mt19937_64 g(0xF39A16420261002ULL);
 uniform_int_distribution<int>td(100,400),kd(5,12),coin(0,99),sg(0,1);
 unordered_set<string>seen;seen.reserve(2500);array<int,6>bins{};
 long long pairs=0,wf=0,mf=0,ff=0,df=0;
 int tries=0,wlo=401,whi=0,minlab=400,maxlab=0,maxnf=0,maxbits=0,lowlabels=0;
 while(seen.size()<2000) {
  if(++tries>1000000)throw runtime_error("sampler stalled");
  int target=td(g),k=kd(g),avg=target/k,lo=max(1,avg/2),hi=min(130,max(lo+2,3*avg/2));
  Profile<400>r;vector<int>labs;
  for(int j=0;j<k;++j) {
   bool dup=j>=2&&!labs.empty()&&coin(g)<30;
   if(dup){int n=labs[g()%labs.size()];r.c[n]+=r.c[n]>0?1:-1;continue;}
   int n=0;
   if(j>=2&&coin(g)<25){int z=1+g()%2;if(!r.c[z])n=z;}
   if(!n) {
    int l=max(3,lo);
    for(int t=0;t<500;++t){int z=l+g()%(hi-l+1);if(!r.c[z]){n=z;break;}}
   }
   if(!n)for(int z=3;z<=400;++z)if(!r.c[z]){n=z;break;}
   if(!n)break;r.c[n]=sg(g)?1:-1;labs.push_back(n);
  }
  r.w=0;r.mx=0;r.nf=k;int core=0;
  for(int n=1;n<=400;++n)if(r.c[n]){r.w+=n*abs((int)r.c[n]);r.mx=n;if(n>=3)core+=abs((int)r.c[n]);}
  if(r.w<100||r.w>400||core<2||lowp(r)>highp(r))continue;
  string key;for(int n=1;n<=400;++n)if(r.c[n])key+=to_string(n)+":"+to_string((int)r.c[n])+",";
  if(!seen.insert(key).second)continue;
  int bits=0;if(!boundok(r,bits))throw runtime_error("parent coefficient bound");
  maxbits=max(maxbits,bits);maxnf=max(maxnf,k);wlo=min(wlo,r.w);whi=max(whi,r.w);
  ++bins[min(5,(r.w-100)/50)];bool has=false;for(int n=1;n<=2;++n)has|=r.c[n]!=0;
  if(has)++lowlabels;
  for(int n=1;n<=400;++n)if(r.c[n]){minlab=min(minlab,n);maxlab=max(maxlab,n);}
  r.g=eval<400>(r);auto rs=removals(r);vector<Profile<400>>kids;kids.reserve(rs.size());
  for(Rem q:rs){auto c=child(r,q);int z;if(!boundok(c,z))throw runtime_error("child coefficient bound");c.g=eval<400>(c);kids.push_back(move(c));}
  unordered_map<string,int>ri;for(int j=0;j<(int)rs.size();++j)ri[keyR<400>(rs[j])]=j;
  auto at=[&](Rem q){auto it=ri.find(keyR<400>(q));if(it==ri.end())throw runtime_error("rule chose invalid removal");return it->second;};
  int iw=at(ruleW(r)),im=at(ruleM(r)),ifq=at(ruleF(r));
  for(int p=lowp(r);p<=highp(r);p+=2) {
   ++pairs;I parent=r.g[p/2];
   auto val=[&](int j)->I{return p>kids[j].w?0:kids[j].g[p/2];};
   if(parent-val(iw)<0)++wf;if(parent-val(im)<0)++mf;if(parent-val(ifq)<0)++ff;
   I best=-(I(1)<<126);for(int j=0;j<(int)kids.size();++j)best=max(best,parent-val(j));
   if(best<0)++df;
  }
 }
 cout<<"random_seed=0xF39A16420261002 backgrounds="<<seen.size()<<" attempts="<<tries
     <<" W="<<wlo<<".."<<whi<<" label_range="<<minlab<<".."<<maxlab
     <<" profiles_with_1_or_2="<<lowlabels<<" max_factors="<<maxnf
     <<" max_norm_bound_bits="<<maxbits<<"\nW_bins(100-149,150-199,200-249,250-299,300-349,350-400)=";
 for(int j=0;j<6;++j)cout<<(j?",":"")<<bins[j];
 cout<<"\nrandom_residual_(B,p)_pairs="<<pairs<<"\nRule_W_fail="<<wf<<"/"<<pairs
     <<"\nRule_M_fail="<<mf<<"/"<<pairs<<"\nMost_frequent_pair_fail="<<ff<<"/"<<pairs
     <<"\nD_fail="<<df<<"/"<<pairs<<"\n";
 if(seen.size()!=2000||pairs!=99533||wf||mf||ff||df)throw runtime_error("random screen totals");
}
int main(int argc,char**argv) {
 selftest<36>(); selftest<120>();
 if(argc>1&&string(argv[1])=="random")randomScreen();else exhaustive();
}
'''
fd=os.memfd_create("fm_chk95_verifier",0);os.set_inheritable(fd,True)
env=dict(os.environ);env["TMPDIR"]="/dev/shm";env["OMP_NUM_THREADS"]="32"
b=subprocess.run(["g++","-pipe","-std=c++17","-O3","-fopenmp","-x","c++","-o",f"/proc/self/fd/{fd}","-"],input=cpp.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(fd,),env=env)
if b.returncode:print(b.stderr.decode(),file=sys.stderr);raise SystemExit(b.returncode)
os.fchmod(fd,0o700)
for args in ([],["random"]):
 p=subprocess.run([f"/proc/self/fd/{fd}"]+args,pass_fds=(fd,),env=env)
 if p.returncode:raise SystemExit(p.returncode)
