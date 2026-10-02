python3 -u - <<'PY'
import os
import subprocess

src = r'''#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>
#include <random>

using namespace std; using Z=mpz_class;
struct Cand{string fam;vector<int>a;}; static vector<Cand> cs; static set<string> seen;
static int wt(const vector<int>&a){int w=0;for(int x:a)w+=abs(x);return w;}
static string key(vector<int>a){sort(a.begin(),a.end());string s;for(int x:a)s+=to_string(x)+",";return s;}
static void add(string f,vector<int>a){if(a.empty()||a.size()>60||wt(a)>1000)return;int neg=0;for(int x:a)neg+=x<0;if(neg&1)return;if(seen.insert(key(a)).second){sort(a.begin(),a.end(),[](int x,int y){if(abs(x)!=abs(y))return abs(x)<abs(y);return x<y;});cs.push_back({f,a});}}
static void gen(){
 add("LP-seed-14",{1,1,2,2,2});{vector<int>a={1,1,2,2,2};a.insert(a.end(),14,-3);add("LP-seed-14",a);}
 {vector<int>a={1,-2,-4,-4,-4,1,6};for(int i=0;i<6;i++)a.push_back(3);for(int i=0;i<4;i++)a.push_back(5);add("SEC181-segregated",a);}
 add("sign-minimum-1to8",{-1,2,3,4,-5,6,7,8});
 for(int q=1;q<=27;q++){vector<int>a;for(int x:{-1,2,3,4,-5,6,7,8})a.push_back(x*q);add("scaled-sign-minimum",a);}
 for(int q=1;q<=5;q++){vector<int>a;for(int i=1;i<=17;i++)a.push_back(-i*q);a.push_back(-17*q);add("scaled-allminus-1to17",a);}
 {vector<int>a;for(int i=1;i<=41;i++)a.push_back(-i);a.push_back(-43);add("42-factor-no-flip",a);}
 for(int k=2;k<=59;k++){int sig=(k&1)?-1:1;for(int dp:{0,1,3,k,2*k,100}){int p=k+dp;vector<int>a;for(int i=1;i<=k;i++)a.push_back(-i);a.push_back(sig*p);add("negative-run",a);}}
 vector<int>C={3,4,5,6,8,10,16,24,32,48,64,96,128,192,256,384,500};
 vector<pair<int,int>>xy={{1,1},{2,2},{4,2},{2,6},{8,4},{12,8},{20,4},{4,20},{16,16},{30,10},{10,30},{0,12},{12,0},{0,30},{30,0}};
 vector<int>zz={2,4,6,8,10,12,16,20,24,30,40};
 for(int c:C)for(int z:zz)for(auto [x,y]:xy){if(x+y+z>60)continue;vector<int>a(x,1);a.insert(a.end(),y,2);a.insert(a.end(),z,-c);add("one-heavy-minus",a);}
 vector<int>D={3,4,5,8,12,20,32,48,64,96,160,256},mult={1,2,3,4,6,8,10,12,16,20};
 for(int i=0;i<(int)D.size();i++)for(int j=i+1;j<(int)D.size();j++)for(int z:mult)for(int v:mult){if((z+v)&1)continue;int x=(D[i]*z+D[j]*v)/2,y=max(0,(D[i]*z+D[j]*v-x)/2);if(z+v+x+y>60)continue;vector<int>a(x,1);a.insert(a.end(),y,2);a.insert(a.end(),z,-D[i]);a.insert(a.end(),v,-D[j]);add("two-heavy-minus-balanced",a);}
 for(int k=4;k<=60;k+=2){vector<int>a,b;for(int i=1;i<=k;i++){a.push_back((i%3==0||i%5==0)?-i:i);b.push_back(i<=k/2?-i:i);}int n=count_if(a.begin(),a.end(),[](int x){return x<0;});if(n&1)a[0]=-a[0];if((k/2)&1)b[k/2]=-b[k/2];add("opposite-sign-classes",a);add("opposite-sign-half-split",b);}
 mt19937_64 rng(18920261002ULL);
 for(int rep=0;rep<260;rep++){int L=10+rng()%51,cap=max(1,min(90,1000/L));vector<int>a;map<int,int>mults;for(int i=0;i<L;i++){int n=1+rng()%cap;a.push_back(n);mults[n]++;}map<int,int>sg;int neg=0;for(auto [n,m]:mults){int s=(rng()&1)?1:-1;sg[n]=s;neg+=s<0?m:0;}if(neg&1){bool done=false;for(auto it=mults.rbegin();it!=mults.rend();++it)if(it->second&1){sg[it->first]*=-1;done=true;break;}if(!done)continue;}for(int&x:a)x*=sg[x];add("random-pairfree",a);}
 for(int c:{3,4,5,7,10,16,24,32,48,64,96,128})for(int z:{2,4,6,8,10,12,16,20,24,30,40,50}){int target=c*z;for(int x=max(0,target-48);x<=target+48;x+=max(1,target/8)){int y=max(0,(target-x+1)/2);if(x+y+z>60)continue;vector<int>a(x,1);a.insert(a.end(),y,2);a.insert(a.end(),z,-c);add("mass-balanced-one-class",a);}}
 // High-weight balanced mass: choose two positive classes just below the heavy negative class.
 for(int c:{8,12,16,24,32,48,64,96,128})for(int z:{2,4,6,8,10,12,16,20,24,30}){int room=60-z,target=c*z;if(room<=0)continue;int a=max(2,target/room);if(a>=c)a=c-2;int b=a+1;if(b>=c)continue;for(int iy=0;iy<=4;iy++){int y=room*iy/4;int x0=(target-b*y+a/2)/a;for(int dx=-1;dx<=1;dx++){int x=x0+dx;if(x<0||x+y+z>60)continue;vector<int>v(x,a);v.insert(v.end(),y,b);v.insert(v.end(),z,-c);add("mass-balanced-large",v);}}}
}
struct Res{Z phi,l1;};
static Res phi_dp(const vector<int>&a){
 int W=wt(a),d=W+1;vector<Z>cur((size_t)d*d),nx((size_t)d*d);auto ix=[&](int x,int y){return(size_t)x*d+y;};cur[ix(0,0)]=1;int old=0,total=0;vector<Z>p0(d),p1(d);
 for(int val:a){int n=abs(val),eps=val>0?1:-1,newm=old+n;total+=n;for(int x=0;x<=newm;x++)for(int y=0;y<=newm;y++)nx[ix(x,y)]=0;
  for(int y=0;y<=old;y++){Z r0=0,r1=0;for(int i=0;i<=old;i++){if(i&1)r1+=cur[ix(i,y)];else r0+=cur[ix(i,y)];p0[i]=r0;p1[i]=r1;}for(int x=0;x<=newm;x++){if((x+y-total)%2)continue;int lo=abs(x-n),hi=min(old,x+n);if(lo>hi)continue;int par=(x+n)&1;Z s=par?p1[hi]:p0[hi];if(lo)s-=par?p1[lo-1]:p0[lo-1];nx[ix(x,y)]+=s;}}
  for(int x=0;x<=old;x++){Z r0=0,r1=0;for(int j=0;j<=old;j++){if(j&1)r1+=cur[ix(x,j)];else r0+=cur[ix(x,j)];p0[j]=r0;p1[j]=r1;}for(int y=0;y<=newm;y++){if((x+y-total)%2)continue;int lo=abs(y-n),hi=min(old,y+n);if(lo>hi)continue;int par=(y+n)&1;Z s=par?p1[hi]:p0[hi];if(lo)s-=par?p1[lo-1]:p0[lo-1];if(eps>0)nx[ix(x,y)]+=s;else nx[ix(x,y)]-=s;}}
  cur.swap(nx);old=newm;
 }
 Res out;out.phi=cur[ix(0,0)];for(int x=0;x<=W;x++)for(int y=0;y<=W;y++){const Z&z=cur[ix(x,y)];out.l1+=z>=0?z:-z;}return out;
}

static Z gp_pruned(vector<int> lab,int p){
 int total=0;for(int x:lab)total+=abs(x);if(p>total)return 0;
 sort(lab.begin(),lab.end(),[](int a,int b){return abs(a)>abs(b);});
 int d=total+1;auto ix=[&](int s,int t){return(size_t)s*d+t;};
 vector<Z>cur((size_t)d*d),nx((size_t)d*d),p0(d),p1(d);cur[ix(0,0)]=1;
 int processed=0,rem=total,clo=0,chi=0,crem=total;
 int nlo=0,nhi=-1,nrem=-1;
 for(int val:lab){
  int n=abs(val),eps=val>0?1:-1;int oldproc=processed,oldrem=rem,oldlo=clo,oldhi=chi;
  rem-=n;processed+=n;
  int lo=max(0,p-rem),hi=min(processed,p+rem);
  if(nhi>=nlo)for(int s=nlo;s<=nhi;s++)for(int t=0;t<=nrem;t++)nx[ix(s,t)]=0;
  if(lo<=hi){
   int tmax=min(oldrem,rem);
   for(int t=0;t<=tmax;t++){
    Z r0=0,r1=0;for(int s=0;s<=oldproc;s++){if(s&1)r1+=cur[ix(s,t)];else r0+=cur[ix(s,t)];p0[s]=r0;p1[s]=r1;}
    for(int s2=lo;s2<=hi;s2++){if((s2+t-processed)%2)continue;int a=abs(s2-n),b=min(oldproc,s2+n);if(a>b)continue;int par=(s2+n)&1;Z v=par?p1[b]:p0[b];if(a)v-=par?p1[a-1]:p0[a-1];nx[ix(s2,t)]+=v;}
   }
   int slo=max(lo,oldlo),shi=min(hi,oldhi),t2max=min(processed,rem);
   for(int s=slo;s<=shi;s++){
    Z r0=0,r1=0;for(int t=0;t<=oldrem;t++){if(t&1)r1+=cur[ix(s,t)];else r0+=cur[ix(s,t)];p0[t]=r0;p1[t]=r1;}
    for(int t2=0;t2<=t2max;t2++){if((s+t2-processed)%2)continue;int a=abs(t2-n),b=min(oldproc,t2+n);if(a>b)continue;int par=(t2+n)&1;Z v=par?p1[b]:p0[b];if(a)v-=par?p1[a-1]:p0[a-1];if(eps>0)nx[ix(s,t2)]+=v;else nx[ix(s,t2)]-=v;}
   }
  }
  cur.swap(nx);nlo=oldlo;nhi=oldhi;nrem=oldrem;clo=lo;chi=hi;crem=rem;
 }
 if(p<d)return cur[ix(p,0)];return 0;
}
struct GE{Z phi,plus;bool hasplus;};
static GE eval_gp(const vector<int>&a){
 int j=0;for(int i=1;i<(int)a.size();i++)if(abs(a[i])>abs(a[j]))j=i;int p=abs(a[j]);vector<int>b,ap;
 for(int i=0;i<(int)a.size();i++)if(i!=j){b.push_back(a[i]);ap.push_back(abs(a[i]));}
 GE z;z.phi=2*gp_pruned(b,p);z.plus=2*gp_pruned(ap,p);z.hasplus=(z.plus!=0);return z;
}

static Z get(const map<pair<int,int>,Z>&m,int x,int y){auto i=m.find({x,y});return i==m.end()?Z(0):i->second;}
static Z phi_laurent(const vector<int>&a){map<pair<int,int>,Z>c,nx;c[{0,0}]=1;for(int v:a){int n=abs(v),e=v>0?1:-1;nx.clear();for(const auto&kv:c)for(int j=0;j<=n;j++){nx[{kv.first.first+n-2*j,kv.first.second}]+=kv.second;nx[{kv.first.first,kv.first.second+n-2*j}]+=e*kv.second;}c.swap(nx);}return get(c,0,0)-get(c,2,0)-get(c,0,2)+get(c,2,2);}
static bool pairfree(const vector<int>&a){set<int>s;for(int x:a){if(s.count(-x))return false;s.insert(x);}return true;}
static string word(const vector<int>&a){string s="[";for(size_t i=0;i<a.size();i++){if(i)s+=",";s+=(a[i]>0?"+":"")+to_string(a[i]);}return s+"]";}
static string stamp(){time_t t=time(nullptr);struct tm u;gmtime_r(&t,&u);char b[32];strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&u);return b;}
int main(int argc,char**argv){
 long resume=0,limit=-1;int stride=1;bool self=false,quiet=false;
 for(int i=1;i<argc;i++){string q=argv[i];if(q=="-h"||q=="--help"){puts("FM-SEC189 exact C++/GMP search. Options: --limit N, --resume-after N, --stride N, --self-check. W<=1000, factors<=60. Each completed serial is fsynced to a run-local memfd.");return 0;}else if(q=="--limit"&&i+1<argc)limit=atol(argv[++i]);else if(q=="--resume-after"&&i+1<argc)resume=atol(argv[++i]);else if(q=="--stride"&&i+1<argc)stride=max(1,atoi(argv[++i]));else if(q=="--self-check")self=true;else if(q=="--quiet")quiet=true;else{fprintf(stderr,"bad option %s\n",argv[i]);return 2;}}
 gen();fprintf(stderr,"generated=%zu; resume_after=%ld\n",cs.size(),resume);fflush(stderr);
 if(self){vector<vector<int>>t={{1,1,-2,-2},{-1,2,3,4,-5,6,7,8},{1,1,2,2,2,-3,-3},{-1,-2,3,4}};vector<int>lp={1,1,2,2,2};lp.insert(lp.end(),14,-3);t.push_back(lp);vector<int>hi(30,1);hi.insert(hi.end(),30,-32);if(phi_dp(hi).phi!=eval_gp(hi).phi){cerr<<"high-weight pruned/full mismatch\n";return 3;}mt19937_64 rr(18961);for(int h=0;h<100;h++){int L=2+rr()%8;vector<int>a;int neg=0;for(int i=0;i<L;i++){int n=1+rr()%6,s=(rr()&1)?1:-1;a.push_back(s*n);neg+=s<0;}if(neg&1)a[0]=-a[0];t.push_back(a);}for(auto&a:t){Z x=phi_dp(a).phi,y=phi_laurent(a),z=eval_gp(a).phi;if(x!=y||x!=z){cerr<<"fusion/GMP/Laurent mismatch "<<word(a)<<" table="<<x<<" Laurent="<<y<<" gp="<<z<<"\n";return 3;}}cerr<<"fusion/GMP/Laurent cross-check PASS ("<<t.size()<<" words)\n";}
  int cf=memfd_create("fmsec189-checkpoint",MFD_CLOEXEC);if(cf<0){perror("memfd_create");return 4;}
 map<string,pair<mpq_class,int>>best,bestpf;map<pair<string,int>,pair<mpq_class,int>>bins;map<string,int>fcount,fzero;pair<mpq_class,int>global,globalpf;bool have=false,havepf=false,negative=false;int begin=min<long>(resume,cs.size()),end=limit<0?(int)cs.size():min((int)cs.size(),begin+(int)limit),done=0;
 int bad=0,maxW=0,maxL=0;vector<int>maxWord;
 for(int k=begin;k<end;k++){
  const Cand&c=cs[k];bool chosen=(k<40)||(k%stride==0)||(k==279)||(k==1646)||(k==2289)||(wt(c.a)>=980&&c.a.size()>=55);
  if(c.fam=="negative-run"&&((int)c.a.size()<=21||(int)c.a.size()%3==0||(int)c.a.size()==59))chosen=true;
  if((c.fam=="opposite-sign-classes"||c.fam=="opposite-sign-half-split")&&((int)c.a.size()<=20||(int)c.a.size()%4==0||(int)c.a.size()==60))chosen=true;
  if(c.fam=="random-pairfree"&&k%5==0)chosen=true;
  if(c.fam=="mass-balanced-one-class"&&k%10==0)chosen=true;if(c.fam=="mass-balanced-large"&&k%16==0)chosen=true;
  if(!chosen)continue;
  GE r=eval_gp(c.a);bool negcand=(r.phi<0),badlocal=false;Z check;
  if(negcand){check=phi_laurent(c.a);if(check!=r.phi)badlocal=true;}
  mpq_class q;if(r.hasplus)q=mpq_class(r.phi,r.plus);else q=0;q.canonicalize();
  string ps=r.phi.get_str(),ls=r.plus.get_str(),rs=r.hasplus?q.get_str():"undefined",w=word(c.a);
  {
   if(badlocal){bad=1;fprintf(stderr,"Laurent mismatch %s gp=%s Laurent=%s\n",w.c_str(),ps.c_str(),check.get_str().c_str());}
   fcount[c.fam]++;if(r.phi==0)fzero[c.fam]++;if(wt(c.a)>maxW)maxW=wt(c.a),maxL=(int)c.a.size(),maxWord=c.a;
   if(r.hasplus){int wb=min(4,wt(c.a)/200);auto kk=make_pair(c.fam,wb);if(!bins.count(kk)||q<bins[kk].first)bins[kk]={q,k};}
   if(r.hasplus&&(!best.count(c.fam)||q<best[c.fam].first))best[c.fam]={q,k};
   if(r.hasplus&&pairfree(c.a)&&(!bestpf.count(c.fam)||q<bestpf[c.fam].first))bestpf[c.fam]={q,k};
   if(r.hasplus&&(!have||q<global.first)){global={q,k};have=true;}
   if(r.hasplus&&pairfree(c.a)&&(!havepf||q<globalpf.first)){globalpf={q,k};havepf=true;}
   if(negcand){negative=true;printf("NEGATIVE-LAURENT list=%s gp=%s Laurent=%s\n",w.c_str(),ps.c_str(),check.get_str().c_str());}
   char rec[2048];int nr=snprintf(rec,sizeof(rec),"%d|%s|%s|%s|%s\n",k+1,ps.c_str(),ls.c_str(),rs.c_str(),w.c_str());
   if(lseek(cf,0,SEEK_END)<0||write(cf,rec,nr)!=nr||fsync(cf)){bad=1;perror("checkpoint fsync");}
   done++;if(!quiet||done==1||done%25==0){printf("[%s] progress=%d serial=%d/%zu family=%s W=%d L=%zu Phi=%s norm=%s list=%s\n",stamp().c_str(),done,k+1,cs.size(),c.fam.c_str(),wt(c.a),c.a.size(),ps.c_str(),rs.c_str(),w.c_str());fflush(stdout);}
  }
 }
 if(bad)return 6;
 puts("FAMILY COUNTS (evaluated, Phi=0):");for(auto&v:fcount)printf("%s %d %d\n",v.first.c_str(),v.second,fzero[v.first]);puts("FAMILY MINIMA (Phi / all-plus Phi), including pair-free minima separately:");for(auto&v:best){int k=v.second.second;GE r=eval_gp(cs[k].a);printf("%s Phi=%s allplus=%s ratio=%s W=%d L=%zu list=%s pairfree=%s\n",v.first.c_str(),r.phi.get_str().c_str(),r.plus.get_str().c_str(),v.second.first.get_str().c_str(),wt(cs[k].a),cs[k].a.size(),word(cs[k].a).c_str(),pairfree(cs[k].a)?"yes":"no");}
 if(have){int k=global.second;GE r=eval_gp(cs[k].a);printf("GLOBAL MIN family=%s Phi=%s allplus=%s ratio=%s W=%d L=%zu list=%s\n",cs[k].fam.c_str(),r.phi.get_str().c_str(),r.plus.get_str().c_str(),global.first.get_str().c_str(),wt(cs[k].a),cs[k].a.size(),word(cs[k].a).c_str());}
 for(auto&v:bestpf){int k=v.second.second;GE r=eval_gp(cs[k].a);printf("PAIRFREE-MIN family=%s Phi=%s allplus=%s ratio=%s W=%d L=%zu list=%s\n",v.first.c_str(),r.phi.get_str().c_str(),r.plus.get_str().c_str(),v.second.first.get_str().c_str(),wt(cs[k].a),cs[k].a.size(),word(cs[k].a).c_str());}
 if(havepf){int k=globalpf.second;GE r=eval_gp(cs[k].a);printf("GLOBAL PAIRFREE MIN family=%s Phi=%s allplus=%s ratio=%s W=%d L=%zu list=%s\n",cs[k].fam.c_str(),r.phi.get_str().c_str(),r.plus.get_str().c_str(),globalpf.first.get_str().c_str(),wt(cs[k].a),cs[k].a.size(),word(cs[k].a).c_str());}puts("MINIMA BY W BIN (0-199,200-399,400-599,600-799,800-1000):");for(auto&v:bins){int k=v.second.second;GE r=eval_gp(cs[k].a);printf("%s bin=%d-%d ratio=%s W=%d L=%zu\n",v.first.first.c_str(),v.first.second*200,v.first.second==4?1000:v.first.second*200+199,v.second.first.get_str().c_str(),wt(cs[k].a),cs[k].a.size());}
 printf("SAMPLE MAX W=%d L=%d list=%s\n",maxW,maxL,word(maxWord).c_str());printf("SCANNED=%d serials=%d..%d/%zu stride=%d negative=%s\n",done,begin+1,end,cs.size(),stride,negative?"YES":"NO");close(cf);return negative?5:0;
}'''

fd = os.memfd_create("fmsec189-bin", 0)
env = {**os.environ, "TMPDIR": "/dev/shm"}
compile_run = subprocess.run(
    ["g++", "-x", "c++", "-std=c++17", "-O3", "-march=native",
     "-o", f"/proc/self/fd/{fd}", "-", "-lgmpxx", "-lgmp"],
    input=src.encode(), pass_fds=(fd,), env=env)
if compile_run.returncode:
    raise SystemExit(compile_run.returncode)

os.lseek(fd, 0, os.SEEK_SET)
help_run = subprocess.run([f"/proc/self/fd/{fd}", "--help"], pass_fds=(fd,))
if help_run.returncode:
    raise SystemExit(help_run.returncode)

scan = subprocess.run(
    [f"/proc/self/fd/{fd}", "--self-check", "--stride", "40"],
    pass_fds=(fd,), timeout=1800)
raise SystemExit(scan.returncode)
PY
