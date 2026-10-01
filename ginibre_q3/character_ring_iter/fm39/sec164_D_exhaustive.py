import os, subprocess, sys

if any(x in ("-h", "--help") for x in sys.argv[1:]):
    print("Usage: python3 -u - [36|40] < verifier; memory-only exact C++ census.")
    raise SystemExit(0)

source = r"""#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using I=__int128_t; using U=__uint128_t;
constexpr int CAP=40;
struct Rem{int n=0,m=0;};
struct Rec{array<signed char,CAP+1> c{}; U key=0; int W=0,nf=0,mx=0,cores=0; array<I,CAP/2+1> g{};};
struct HashU{size_t operator()(U x)const noexcept{uint64_t a=(uint64_t)x,b=(uint64_t)(x>>64);a^=b+0x9e3779b97f4a7c15ULL+(a<<6)+(a>>2);return (size_t)a;}};
struct Kid{const Rec* r; Rem R;};
struct Rat{bool ok=false;I n=0;U d=1;int W=0,p=0;U key=0;Rem R{};I parent=0,child=0,delta=0;};
struct Wit{bool ok=false;int W=0,p=0;U key=0;Rem R{};I parent=0,child=0,delta=0;};
static int LIM=36; static vector<Rec> V; static array<signed char,CAP+1> cur{}; static array<U,CAP+1> place{}; static unordered_map<U,size_t,HashU> ix;
static U absi(I x){return x<0?U(-(x+1))+1:U(x);}
static string su(U x){if(!x)return "0";string s;while(x){s.push_back(char('0'+x%10));x/=10;}reverse(s.begin(),s.end());return s;}
static string si(I x){return x<0?"-"+su(absi(x)):su(U(x));}
static U gcdU(U a,U b){while(b){U t=a%b;a=b;b=t;}return a;}
static bool lessPos(U a,U b,U c,U d){bool flip=false;for(;;){U qa=a/b,qc=c/d;if(qa!=qc)return flip?qa>qc:qa<qc;U ra=a%b,rc=c%d;if(!ra||!rc){if(!ra&&!rc)return false;bool firstSmaller=!ra;return flip?!firstSmaller:firstSmaller;}a=b;b=ra;c=d;d=rc;flip=!flip;}}
static bool lessRat(I a,U b,I c,U d){if(a<0&&c>=0)return true;if(a>=0&&c<0)return false;if(a<0)return lessPos(absi(c),d,absi(a),b);return lessPos(U(a),b,U(c),d);}
static bool earlier(const Wit&a,const Wit&b){return !b.ok||(a.W!=b.W?a.W<b.W:(a.key!=b.key?a.key<b.key:a.p<b.p));}
static bool earlier(const Rat&a,const Rat&b){if(!b.ok)return true;if(lessRat(a.n,a.d,b.n,b.d))return true;if(lessRat(b.n,b.d,a.n,a.d))return false;return a.W!=b.W?a.W<b.W:(a.key!=b.key?a.key<b.key:a.p<b.p);}
static string showB(const Rec&r){string s="{";bool f=true;for(int n=1;n<=LIM;++n)if(r.c[n]){if(!f)s+=" ";f=false;s+=(r.c[n]>0?"+":"-")+to_string(n)+"^"+to_string(abs((int)r.c[n]));}return s+"}";}
static string showR(Rem R,const Rec&r){auto lab=[&](int n){return string(r.c[n]>0?"+":"-")+to_string(n);};return R.m?("["+lab(R.n)+","+lab(R.m)+"]"):("["+lab(R.n)+"]");}
static void gen(int n,int w,int nf,int mx,int cores,U key){if(n>LIM){Rec r;r.c=cur;r.key=key;r.W=w;r.nf=nf;r.mx=mx;r.cores=cores;V.push_back(r);return;}int M=LIM/n;for(int c=-M;c<=M;++c){int nw=w+n*abs(c);if(nw>LIM)continue;cur[n]=(signed char)c;gen(n+1,nw,nf+abs(c),c?n:mx,cores+(n>=3?abs(c):0),key+U(c+M)*place[n]);}cur[n]=0;}
static inline I sumRange(I pref[2][CAP+2],int par,int lo,int hi){lo=max(lo,0);hi=min(hi,CAP);return lo>hi?0:pref[par][hi+1]-pref[par][lo];}
static array<I,CAP/2+1> calc(const Rec&r){I d[CAP+1][CAP+1]{},q[CAP+1][CAP+1]{};d[0][0]=1;int sm=0;for(int n=1;n<=LIM;++n)for(int rep=0;rep<abs((int)r.c[n]);++rep){int eps=r.c[n]>0?1:-1,next=sm+n;memset(q,0,sizeof(q));I pref[2][CAP+2]{};for(int b=0;b<=sm;++b){pref[0][0]=pref[1][0]=0;for(int x=0;x<=sm;++x){pref[0][x+1]=pref[0][x];pref[1][x+1]=pref[1][x];pref[x&1][x+1]+=d[x][b];}for(int t=0;t<=next;++t){int lo=abs(t-n),hi=min(sm,t+n),par=((t-n)%2+2)%2;if(lo<=hi)q[t][b]+=sumRange(pref,par,lo,hi);}}for(int a=0;a<=sm;++a){pref[0][0]=pref[1][0]=0;for(int y=0;y<=sm;++y){pref[0][y+1]=pref[0][y];pref[1][y+1]=pref[1][y];pref[y&1][y+1]+=d[a][y];}for(int t=0;t<=next;++t){int lo=abs(t-n),hi=min(sm,t+n),par=((t-n)%2+2)%2;if(lo<=hi)q[a][t]+=eps*sumRange(pref,par,lo,hi);}}memcpy(d,q,sizeof(d));sm=next;}array<I,CAP/2+1> out{};for(int p=(r.W&1);p<=r.W;p+=2)out[p/2]=d[p][0];return out;}
static I gp(const Rec&r,int p){return p<0||p>r.W||((r.W-p)&1)?0:r.g[p/2];}
static U childkey(const Rec&r,Rem R){U k=r.key;auto one=[&](int n){if(r.c[n]>0)k-=place[n];else k+=place[n];};if(!R.m)one(R.n);else if(R.n==R.m){if(r.c[R.n]>0)k-=2*place[R.n];else k+=2*place[R.n];}else{one(R.n);one(R.m);}return k;}
static void update(Rat&best,I delta,I parent,I child,const Rec&r,int p,Rem R){if(parent==0)return;Rat a;a.ok=true;a.n=parent<0?-delta:delta;a.d=absi(parent);a.W=r.W;a.p=p;a.key=r.key;a.R=R;a.parent=parent;a.child=child;a.delta=delta;if(earlier(a,best))best=a;}
static void printRatio(const Rat&r){if(!r.ok){cout<<"min_ratio=none\n";return;}U g=gcdU(absi(r.n),r.d);const Rec&b=V[ix.at(r.key)];cout<<"min_ratio="<<si(r.n/(I)g)<<"/"<<su(r.d/g)<<" at W="<<r.W<<" B="<<showB(b)<<" p="<<r.p<<" R="<<showR(r.R,b)<<" parent="<<si(r.parent)<<" child="<<si(r.child)<<" max_delta="<<si(r.delta)<<"\n";}
int main(int argc,char**argv){if(argc>1&&(string(argv[1])=="-h"||string(argv[1])=="--help")){cout<<"Usage: d_census [max_weight<=40]; pair-free residual, all even-weight removals.\n";return 0;}LIM=argc>1?atoi(argv[1]):36;if(LIM<1||LIM>40){cerr<<"max weight must be 1..40\n";return 2;}static_assert(sizeof(I)==16);if((U(1)<<80)>=(U(1)<<127))return 2;U pos=1;for(int n=1;n<=LIM;++n){place[n]=pos;U b=U(2*(LIM/n)+1);if(pos>(~U(0))/b)return 2;pos*=b;}auto start=chrono::steady_clock::now();gen(1,0,0,0,0,0);sort(V.begin(),V.end(),[](const Rec&a,const Rec&b){return a.W!=b.W?a.W<b.W:a.key<b.key;});cout<<"generated_profiles="<<V.size()<<" max_weight="<<LIM<<"\n"<<flush;if((LIM==36&&V.size()!=2298873)||(LIM==40&&V.size()!=6011121)){cerr<<"profile count mismatch\n";return 3;}ix.max_load_factor(.7f);ix.reserve(V.size()*1.2);for(size_t i=0;i<V.size();++i){ix.emplace(V[i].key,i);if((i+1)%500000==0)cout<<"phase=index "<<(i+1)<<"/"<<V.size()<<"\n"<<flush;}
unsigned long long prof=0,resbg=0,respairs=0,checks=0,fails=0,nomove=0,keymiss=0,pospar=0,zeropar=0,negpar=0;Rat minrat;Wit first;
auto low=[&](int w){return lower_bound(V.begin(),V.end(),w,[](const Rec&r,int z){return r.W<z;})-V.begin();};int nt=max(1,omp_get_max_threads());for(int W=0;W<=LIM;++W){size_t lo=low(W),hi=low(W+1);const size_t CH=50000;for(size_t base=lo;base<hi;base+=CH){size_t end=min(hi,base+CH);
#pragma omp parallel for num_threads(nt) schedule(dynamic,16)
for(long long i=(long long)base;i<(long long)end;++i)V[(size_t)i].g=calc(V[(size_t)i]);for(size_t i=base;i<end;++i){const Rec&r=V[i];if(r.cores<2)continue;int p0=max(6,r.mx);if((p0&1)!=(r.W&1))++p0;int pmax=min(r.W-16,r.W-2*r.mx);if(p0>pmax)continue;++resbg;array<Kid,900> kids{};int nk=0;auto add=[&](Rem R){auto it=ix.find(childkey(r,R));if(it==ix.end()){++keymiss;return;}kids[nk++]={&V[it->second],R};};for(int n=2;n<=LIM;n+=2)if(r.c[n])add({n,0});for(int n=1;n<=LIM;++n)if(r.c[n])for(int m=n;m<=LIM;++m)if(r.c[m]&&(n&1)==(m&1)&&(m!=n||abs((int)r.c[n])>=2))add({n,m});if(!nk){++nomove;continue;}for(int p=p0;p<=pmax;p+=2){int delta=(r.W-p)/2;if(delta<8||r.mx>delta)continue;++respairs;I parent=gp(r,p),bc=0;bool have=false;Rem br{};for(int z=0;z<nk;++z){const Rec&ch=*kids[z].r;I cv=p>ch.W?0:gp(ch,p);++checks;if(!have||cv<bc){have=true;bc=cv;br=kids[z].R;}}I d=parent-bc;if(d<0){++fails;Wit w{true,r.W,p,r.key,br,parent,bc,d};if(earlier(w,first))first=w;}if(parent>0)++pospar;else if(parent==0)++zeropar;else ++negpar;update(minrat,d,parent,bc,r,p,br);}}
prof+=end-base;double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();cout<<"progress W="<<W<<" profiles="<<prof<<"/"<<V.size()<<" residual_pairs="<<respairs<<" removal_comparisons="<<checks<<" elapsed_s="<<fixed<<setprecision(1)<<sec<<"\n"<<flush;}
if(W==32||W==36||W==40||W==LIM){cout<<"SUMMARY upto_W="<<W<<" pair_free_profiles="<<prof<<" residual_backgrounds="<<resbg<<" residual_pairs="<<respairs<<" removal_comparisons="<<checks<<" D_failures="<<fails<<" no_removal="<<nomove<<" key_misses="<<keymiss<<" parent(+/0/-)="<<pospar<<"/"<<zeropar<<"/"<<negpar<<"\n";printRatio(minrat);if(first.ok){const Rec&b=V[ix.at(first.key)];cout<<"D_failure W="<<first.W<<" B="<<showB(b)<<" p="<<first.p<<" best_R="<<showR(first.R,b)<<" parent="<<si(first.parent)<<" child="<<si(first.child)<<" max_delta="<<si(first.delta)<<"\n";}else cout<<"no residual pair without a monotone removal so far\n";cout<<flush;}}
return 0;}
"""

limit = int(sys.argv[1]) if len(sys.argv) > 1 else 36
if not 1 <= limit <= 40:
    raise SystemExit("max weight must be in 1..40")
fd = os.memfd_create("fm_sec164_D_census", 0)
os.set_inheritable(fd, True)
env = dict(os.environ)
env["TMPDIR"] = "/dev/shm"
env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "8")
build = subprocess.run(
    ["g++", "-pipe", "-std=c++17", "-O3", "-fopenmp", "-x", "c++",
     "-o", f"/proc/self/fd/{fd}", "-"],
    input=source.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    pass_fds=(fd,), env=env
)
if build.returncode:
    print(build.stderr.decode(), file=sys.stderr)
    raise SystemExit(build.returncode)
os.fchmod(fd, 0o700)
print(f"compiled in memory; max_weight={limit}; threads={env['OMP_NUM_THREADS']}", flush=True)
run = subprocess.run([f"/proc/self/fd/{fd}", str(limit)], pass_fds=(fd,), env=env)
raise SystemExit(run.returncode)