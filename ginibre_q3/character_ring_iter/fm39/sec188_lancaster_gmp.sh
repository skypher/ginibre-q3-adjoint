cd /home/yang/q3adjoint
python3 -u - <<'FM188PY'
import os, subprocess, sys
os.environ["TMPDIR"]="/dev/shm"
if any(a in ("-h","--help") for a in sys.argv[1:]):
    print("FM-SEC188 exact GMP Lancaster search; --help prints this message.")
    raise SystemExit(0)
cpp = r'''
#include <gmpxx.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;
static string now(){time_t t=time(nullptr);tm q;gmtime_r(&t,&q);char b[32];strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&q);return b;}
static void logx(const string&s){cout<<now()<<" "<<s<<endl;}
static int wt(const vector<int>&w){int s=0;for(int z:w)s+=abs(z);return s;}
static int negs(const vector<int>&w){int s=0;for(int z:w)s+=(z<0);return s;}
static void canon(vector<int>&w){sort(w.begin(),w.end(),[](int a,int b){return abs(a)!=abs(b)?abs(a)>abs(b):a<b;});}
static string wordstr(const vector<int>&w){ostringstream o;o<<"(";for(size_t i=0;i<w.size();++i){if(i)o<<",";o<<w[i];}o<<")";return o.str();}
struct Tab{int cap;vector<mpz_class>v;Tab(int n=0):cap(n),v((size_t)(n+1)*(n+1)){}mpz_class&at(int r,int s){return v[(size_t)r*(cap+1)+s];}};
static Tab build(vector<int>w){
 canon(w);Tab old(0);old.at(0,0)=1;int S=0;
 for(int z:w){int n=abs(z),ep=z>0?1:-1,T=S+n;Tab out(T);
  for(int s=0;s<=S;++s){vector<mpz_class>p[2];p[0].resize(S+1);p[1].resize(S+1);
   for(int a=0;a<=S;++a){for(int h=0;h<2;++h)p[h][a]=(a?p[h][a-1]:mpz_class(0));p[a&1][a]+=old.at(a,s);}
   for(int r=0;r<=T;++r){int lo=abs(r-n),hi=min(S,r+n),h=((r-n)%2+2)%2;if(lo<=hi){mpz_class q=p[h][hi];if(lo)q-=p[h][lo-1];out.at(r,s)+=q;}}
  }
  for(int r=0;r<=S;++r){vector<mpz_class>p[2];p[0].resize(S+1);p[1].resize(S+1);
   for(int b=0;b<=S;++b){for(int h=0;h<2;++h)p[h][b]=(b?p[h][b-1]:mpz_class(0));p[b&1][b]+=old.at(r,b);}
   for(int s=0;s<=T;++s){int lo=abs(s-n),hi=min(S,s+n),h=((s-n)%2+2)%2;if(lo<=hi){mpz_class q=p[h][hi];if(lo)q-=p[h][lo-1];if(ep<0)q=-q;out.at(r,s)+=q;}}
  }
  old=std::move(out);S=T;
 }
 return old;
}
struct Rho{int p,q;};
static void addrho(set<pair<int,int>>&R,int p,int q){int g=gcd(abs(p),q);if(g){p/=g;q/=g;}R.insert({p,q});}
static string rstr(int p,int q){return to_string(p)+"/"+to_string(q);}
static mpz_class evalnum(const vector<mpz_class>&c,int p,int q,mpz_class&den){
 int d=(int)c.size()-1;mpz_class v=c[d],qp=1;
 for(int r=d-1;r>=0;--r){qp*=q;v=v*p+c[r]*qp;}
 den=qp;return v;
}
static string coeffstr(const vector<mpz_class>&c){ostringstream o;o<<"[";for(size_t i=0;i<c.size();++i){if(i)o<<",";o<<c[i];}o<<"]";return o.str();}
static vector<int> word3(int a,int x,int b,int y,int c,int z){vector<int>w;for(int i=0;i<x;++i)w.push_back(a);for(int i=0;i<y;++i)w.push_back(b);for(int i=0;i<z;++i)w.push_back(-c);canon(w);return w;}
static void append(vector<int>&w,int z,int k){for(int i=0;i<k;++i)w.push_back(z);}
struct Best{bool set=false;mpq_class val;int p=0,q=1;string family;vector<int>w;vector<mpz_class>c;};
int main(int argc,char**argv){
 for(int i=1;i<argc;++i)if(string(argv[i])=="-h"||string(argv[i])=="--help"){cout<<"FM-SEC188 exact scan, W<=200 and at most 40 factors."<<endl;return 0;}
 set<pair<int,int>>R;for(int k=-99;k<=99;++k)addrho(R,k,100);for(int k=-19;k<=19;++k)addrho(R,k,20);
 for(int j=2;j<=20;++j){int q=1<<j;addrho(R,q-1,q);addrho(R,-q+1,q);}
 vector<Rho>grid;for(auto [p,q]:R)if(abs(p)<q)grid.push_back({p,q});
 set<vector<int>>seen;map<string,Best>groups;Best global;long long zero=0,oddW=0,invalid=0;int evaluated=0,failures=0;
 auto valid=[&](const vector<int>&w){return w.size()>=2&&w.size()<=40&&wt(w)<=200&&negs(w)%2==0;};
 auto add=[&](string fam,vector<int>w){
  canon(w);if(!valid(w)){++invalid;return;}if(!seen.insert(w).second)return;
  int W=wt(w);if(W&1){++oddW;++zero;return;}
  Tab t=build(w);vector<mpz_class>c;for(int r=0;2*r<=W;++r)c.push_back(t.at(r,r));
  while(c.size()>1&&c.back()==0)c.pop_back();
  mpz_class norm=0;for(auto&a:c)norm+=(a<0?-a:a);
  ++evaluated;if(norm==0){++zero;return;}
  mpz_class dplus=0,dminus=0;for(size_t r=0;r<c.size();++r){dplus+=c[r];dminus+=(r&1)?-c[r]:c[r];}
  if((negs(w)>0&&dplus!=0)||dplus<0||dminus<0){
   cout<<"LC_ENDPOINT_FAILURE family="<<fam<<" W="<<W<<" word="<<wordstr(w)<<" D1="<<dplus<<" Dm1="<<dminus<<" coeff="<<coeffstr(c)<<endl;++failures;return;
  }
  Best local;local.set=true;bool have=false;
  for(auto rh:grid){mpz_class den,num=evalnum(c,rh.p,rh.q,den);
   if(num<0){mpq_class v(num,den);v.canonicalize();cout<<"LC_FAILURE family="<<fam<<" W="<<W<<" n="<<w.size()<<" word="<<wordstr(w)<<" rho="<<rstr(rh.p,rh.q)<<" D="<<v.get_str()<<" normalized=";mpq_class nv(num,den*norm);nv.canonicalize();cout<<nv.get_str()<<" coeff="<<coeffstr(c)<<endl;++failures;return;}
   mpq_class nv(num,den*norm);nv.canonicalize();
   if(!have||nv<local.val){have=true;local.val=nv;local.p=rh.p;local.q=rh.q;}
  }
  local.family=fam;local.w=w;local.c=c;
  if(!have)return;
  auto it=groups.find(fam);if(it==groups.end()||local.val<it->second.val)groups[fam]=local;
  if(!global.set||local.val<global.val)global=local;
  if(evaluated%100==0)logx("progress polynomials="+to_string(evaluated)+" distinct_words="+to_string(seen.size())+" latest="+fam+" W="+to_string(W)+" rho="+rstr(local.p,local.q)+" normalized="+local.val.get_str());
 };
 add("LP1",word3(1,2,2,3,3,14));add("LP2",word3(1,2,2,7,3,10));
 vector<int>sec181;append(sec181,1,2);append(sec181,-2,1);append(sec181,-4,3);append(sec181,3,6);append(sec181,5,4);append(sec181,6,1);add("SEC181_segregated",sec181);
 vector<int>f1;append(f1,1,2);append(f1,-2,1);append(f1,3,4);append(f1,-4,3);append(f1,5,4);append(f1,8,1);add("F1",f1);
 for(int c=2;c<=12;++c)for(int x=1;x<=4;++x)for(int y=0;y<=5;++y)for(int z=2;z<=24;z+=2)add("one_negative_small",word3(1,x,2,y,c,z));
 for(int y:{3,5,7,9,10})for(int z=2;z<=40;z+=2)add("one_negative_LP_axis",word3(1,2,2,y,3,z));
 for(int c:{15,20,30,40,50,75,100})for(int x=1;x<=3;++x)for(int y=0;y<=3;++y)for(int z=2;z<=8;z+=2)add("one_negative_heavy",word3(1,x,2,y,c,z));
 vector<pair<int,int>>md={{2,3},{3,4},{3,6},{4,5},{4,7},{5,8},{7,11},{10,11},{14,15},{20,21}};
 for(auto [c,d]:md)for(int x=1;x<=2;++x)for(int y=0;y<=3;++y)for(int z=1;z<=8;++z)for(int q=1;q<=8;++q)if((z+q)%2==0){vector<int>w;append(w,1,x);append(w,2,y);append(w,-c,z);append(w,-d,q);add("two_negative_classes",w);}
 for(int k=2;k<=20;k+=2){vector<int>w;for(int j=1;j<=k;++j)w.push_back(-j);add("minus_run",w);for(int p=1;p<=30;++p){auto v=w;v.push_back(p);add("minus_run_plus",v);}}
 for(int k=1;k<=19;k+=2){vector<int>w;for(int j=1;j<=k;++j)w.push_back(-j);for(int p=1;p<=20;++p){auto v=w;v.push_back(-p);add("odd_run_extra_minus",v);}}
 for(int n=1;n<=15;++n)for(int a=1;a<=6;++a)for(int b=2;b<=10;b+=2){vector<int>w;append(w,n,a);append(w,-n,b);add("opposite_one_class",w);}
 for(int n=1;n<=10;++n)for(int m=n+1;m<=12;++m)for(int a=1;a<=3;++a)for(int b=1;b<=4;++b)for(int c=1;c<=3;++c)for(int d=1;d<=4;++d)if((b+d)%2==0){vector<int>w;append(w,n,a);append(w,-n,b);append(w,m,c);append(w,-m,d);add("opposite_two_classes",w);}
 vector<vector<int>>cores;cores.push_back(word3(1,2,2,3,3,14));cores.push_back(word3(1,2,2,7,3,10));cores.push_back(sec181);
 for(auto base:cores)for(int d:{5,7,10,15,20,30,50})for(int m=1;m<=3;++m){auto w=base;append(w,d,m);add("large_positive_class",w);}
 for(auto base:cores)for(int d:{5,10,20})for(int e:{12,25,35})if(e>d)for(int m=1;m<=2;++m)for(int q=1;q<=2;++q){auto w=base;append(w,d,m);append(w,e,q);add("two_large_positive_classes",w);}
 mt19937 rng(188);for(int t=0;t<600;++t){int n=6+rng()%20;vector<int>w;int minus=0;for(int i=0;i<n-1;++i){int a=1+rng()%16;int s=(rng()&1)?1:-1;w.push_back(s*a);minus+=(s<0);}int a=1+rng()%16;int s=(minus&1)?-1:1;w.push_back(s*a);add("random_even_minus",w);}
 cout<<"SUMMARY distinct="<<seen.size()<<" evaluated_nonzero="<<evaluated<<" zero_polynomials_or_oddW="<<zero<<" oddW_trivial="<<oddW<<" invalid="<<invalid<<" failures="<<failures<<" rational_grid="<<grid.size()<<endl;
 if(global.set)cout<<"GLOBAL_MIN family="<<global.family<<" W="<<wt(global.w)<<" n="<<global.w.size()<<" word="<<wordstr(global.w)<<" rho="<<rstr(global.p,global.q)<<" normalized="<<global.val.get_str()<<" coeff="<<coeffstr(global.c)<<endl;
 for(auto&kv:groups){auto&b=kv.second;cout<<"FAMILY_MIN "<<kv.first<<" W="<<wt(b.w)<<" n="<<b.w.size()<<" word="<<wordstr(b.w)<<" rho="<<rstr(b.p,b.q)<<" normalized="<<b.val.get_str()<<endl;}
 return failures?2:0;
}
'''
fd=os.memfd_create("fmsec188_gmp",0)
c=subprocess.run(["g++","-std=c++17","-O3","-x","c++","-","-o",f"/proc/self/fd/{fd}","-lgmpxx","-lgmp"],input=cpp.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(fd,))
if c.returncode:
    print(c.stderr.decode(),flush=True)
    raise SystemExit(c.returncode)
p=subprocess.Popen([f"/proc/self/fd/{fd}"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,pass_fds=(fd,))
for line in p.stdout:
    print(line,end="",flush=True)
rc=p.wait()
print("memfd evaluator exit",rc,flush=True)
FM188PY
