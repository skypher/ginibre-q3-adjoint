import os,subprocess,sys
CPP=r"""#include <bits/stdc++.h>
#include <gmpxx.h>
#include <omp.h>
using namespace std; using Z=mpz_class;
static string stamp(){time_t t=time(nullptr);tm*g=gmtime(&t);char b[32];strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",g);return b;}
static string vs(const vector<int>&v){string s="(";for(size_t i=0;i<v.size();i++){if(i)s+=",";s+=to_string(v[i]);}return s+")";}
struct Tab{int w=0;vector<Z>a;};
static Tab mul(const Tab&in,int z){
 int n=abs(z),eps=z>0?1:-1,ow=in.w,nw=ow+n,od=ow+1,nd=nw+1;Tab out;out.w=nw;out.a.resize((size_t)nd*nd);
 vector<Z>p0(od+1),p1(od+1);
 for(int s=0;s<od;s++){p0[0]=0;p1[0]=0;for(int r=0;r<od;r++){const Z&v=in.a[(size_t)r*od+s];p0[r+1]=p0[r]+((r&1)?0:v);p1[r+1]=p1[r]+((r&1)?v:0);}for(int q=0;q<nd;q++){int lo=abs(n-q),hi=min(ow,n+q),par=(n+q)&1;if((lo&1)!=par)++lo;if(lo<=hi){const auto&p=par?p1:p0;out.a[(size_t)q*nd+s]+=p[hi+1]-p[lo];}}}
 for(int r=0;r<od;r++){p0[0]=0;p1[0]=0;for(int s=0;s<od;s++){const Z&v=in.a[(size_t)r*od+s];p0[s+1]=p0[s]+((s&1)?0:v);p1[s+1]=p1[s]+((s&1)?v:0);}for(int q=0;q<nd;q++){int lo=abs(n-q),hi=min(ow,n+q),par=(n+q)&1;if((lo&1)!=par)++lo;if(lo<=hi){const auto&p=par?p1:p0;out.a[(size_t)r*nd+q]+=eps*(p[hi+1]-p[lo]);}}}
 return out;
}
static Tab table(vector<int>w,bool noisy=false,string name=""){
 sort(w.begin(),w.end(),[](int a,int b){if(abs(a)!=abs(b))return abs(a)>abs(b);return a<b;});Tab t;t.a={Z(1)};
 int i=0;for(int z:w){t=mul(t,z);i++;if(noisy)cerr<<stamp()<<" progress "<<name<<" factor="<<i<<"/"<<w.size()<<" degree="<<t.w<<endl;}return t;
}
struct Metric{bool has=false,fail=false;Z num=0,den=1,pi=0,pos=0;int T=0,badT=-1;vector<int>A,B;};
static bool lessM(const Metric&a,const Metric&b){return a.has&&(!b.has||a.num*b.den<b.num*a.den);}
static Metric profile(vector<int>A,vector<int>B,bool noisy=false,string name="cut"){
 Metric m;m.A=A;m.B=B;Tab x=table(A,noisy,name+"/A"),y=table(B,noisy,name+"/B");int total=x.w+y.w,R=min(x.w,y.w),dx=x.w+1,dy=y.w+1;vector<Z>lay(total+1),pos(total+1);
 for(int r=0;r<=R;r++)for(int s=0;s<=R;s++){const Z&u=x.a[(size_t)r*dx+s];const Z&v=y.a[(size_t)r*dy+s];if(u==0||v==0)continue;Z q=u*v;lay[r+s]+=q;if(q>0)pos[r+s]+=q;}
 Z run=0,pm=0;for(int T=0;T<=total;T++){run+=lay[T];pm+=pos[T];if(run<0&&!m.fail){m.fail=true;m.badT=T;m.pi=run;m.pos=pm;}if(pm>0&&(!m.has||run*m.den<m.num*pm)){m.has=true;m.num=run;m.den=pm;m.T=T;}}
 return m;
}
static void printM(string label,const Metric&m){cout<<label<<" ratio="<<(m.has?m.num.get_str()+"/"+m.den.get_str():"none")<<" A="<<vs(m.A)<<" B="<<vs(m.B)<<" T="<<m.T<<" Pi="<<m.num<<" positive="<<m.den<<" fail="<<m.fail;if(m.fail)cout<<" first_bad_T="<<m.badT<<" first_bad_Pi="<<m.pi<<" first_bad_pos="<<m.pos;cout<<"\n";}
static pair<int,int> toppair(const vector<int>&B){int bi=-1,bj=-1,bs=-1,bm=-1;for(int i=0;i<(int)B.size();i++)for(int j=i+1;j<(int)B.size();j++){int a=abs(B[i]),b=abs(B[j]);if((a-b)%2)continue;int s=a+b,m=max(a,b);if(s>bs||(s==bs&&m>bm)){bi=i;bj=j;bs=s;bm=m;}}return {bi,bj};}
static pair<vector<int>,vector<int>> interior(vector<int>C){sort(C.begin(),C.end(),[](int a,int b){if(abs(a)!=abs(b))return abs(a)>abs(b);return a<b;});vector<int>A,B;long wa=0,wb=0;for(int z:C){if(wa<=wb){A.push_back(z);wa+=abs(z);}else{B.push_back(z);wb+=abs(z);}}if(wa>wb)swap(A,B);return{A,B};}
static Metric ipPair(const vector<int>&L,int i,int j,bool noisy=false,string name="IP"){vector<int>C;for(int k=0;k<(int)L.size();k++)if(k!=i&&k!=j)C.push_back(L[k]);auto ab=interior(C);vector<int>Bp=ab.second;Bp.push_back(L[i]);Bp.push_back(L[j]);return profile(ab.first,Bp,noisy,name);}
static int runp(int k){int q=k*(k+1)/2;for(int p=k+1;p<=k+3;p++)if(((q+p)&1)==0)return p;throw runtime_error("p");}
static vector<int>parseLine(const string&l){size_t a=l.find(" p="),b=l.find(" B="),c=l.find("  phi=");if(a==string::npos||b==string::npos||c==string::npos)throw runtime_error("parse census");int p=stoi(l.substr(a+3,b-a-3));vector<int>w;stringstream ss(l.substr(b+3,c-b-3));int z;while(ss>>z)w.push_back(z);w.push_back(p);return w;}
static vector<vector<int>>readRows(string path){ifstream f(path);if(!f)throw runtime_error("census log open");vector<vector<int>>v;string s;while(getline(f,s))v.push_back(parseLine(s));return v;}
static uint64_t nextR(uint64_t&x){x^=x<<13;x^=x>>7;x^=x<<17;return x;}
static void runCensus(string path,int ns,uint64_t seed){
 auto all=readRows(path);if(ns>(int)all.size())throw runtime_error("sample too large");uint64_t x=seed?seed:1;for(int i=0;i<ns;i++){int j=i+nextR(x)%(all.size()-i);swap(all[i],all[j]);}all.resize(ns);cout<<stamp()<<" START IP-census population="<<readRows(path).size()<<" sample="<<ns<<" seed="<<seed<<" rule=TopPair\n";
 vector<Metric>res(ns);atomic<int>done{0};auto st=chrono::steady_clock::now();
 #pragma omp parallel for schedule(dynamic,1)
 for(int q=0;q<ns;q++){auto L=all[q];vector<int>B(L.begin(),L.end()-1);auto ij=toppair(B);res[q]=ipPair(L,ij.first,ij.second);int d=++done;if(d%100==0){
  #pragma omp critical
  {cerr<<stamp()<<" progress IP-census rows="<<d<<"/"<<ns<<" elapsed_s="<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<endl;}
 }}
 Metric best;int fails=0;map<int,int>lc;for(int q=0;q<ns;q++){if(lessM(res[q],best))best=res[q];fails+=res[q].fail;lc[all[q].size()]++;if(res[q].fail)cout<<"IP-CENSUS-FAIL L="<<vs(all[q])<<" A="<<vs(res[q].A)<<" B="<<vs(res[q].B)<<" T="<<res[q].badT<<" Pi="<<res[q].pi<<"\n";}
 cout<<stamp()<<" IP-census summary rows="<<ns<<" fail="<<fails<<" lengths";for(auto[a,b]:lc)cout<<" "<<a<<":"<<b;cout<<"\n";printM("IP-CENSUS-MIN",best);
}
static void runRuns(int maxk){Metric best;int fails=0;for(int k=5;k<=maxk;k++){int p=runp(k),sig=k%2==0?1:-1;vector<int>B;for(int z=1;z<=k;z++)B.push_back(-z);vector<int>L=B;L.push_back(sig*p);auto ij=toppair(B);auto m=ipPair(L,ij.first,ij.second,k==maxk,"run"+to_string(k));if(lessM(m,best))best=m;fails+=m.fail;cout<<stamp()<<" IP-RUN k="<<k<<" signed_p="<<sig*p<<" pair=("<<B[ij.first]<<","<<B[ij.second]<<") fail="<<m.fail<<"\n";printM("IP-RUN-MIN",m);}cout<<stamp()<<" IP-runs summary k=5.."<<maxk<<" failures="<<fails<<"\n";printM("IP-RUN-GLOBAL",best);}
static void run42(){vector<int>B;for(int z=1;z<=41;z++)B.push_back(-z);vector<int>L=B;L.push_back(-43);auto ij=toppair(B);cout<<stamp()<<" IP-42 TopPair=("<<B[ij.first]<<","<<B[ij.second]<<")\n";auto m=ipPair(L,ij.first,ij.second,true,"FM177");cout<<stamp()<<" IP-42 result fail="<<m.fail<<"\n";printM("IP-42-MIN",m);}
static void runWitness(){vector<vector<int>>Ls={{1,-2,-4,-4,-4,1,3,3,3,3,3,3,5,5,5,5,6},{1,-2,-4,-4,-4,1,3,3,3,3,3,5,5,5,5,5,6},{1,-2,-4,-4,-4,1,1,2,3,3,3,3,3,4,4,5,5,8}};for(int q=0;q<3;q++){Metric best;int n=0,fail=0;auto&L=Ls[q];for(int i=0;i<(int)L.size();i++)for(int j=i+1;j<(int)L.size();j++){auto m=ipPair(L,i,j);n++;fail+=m.fail;if(lessM(m,best))best=m;}cout<<stamp()<<" IP-witness id="<<q+1<<" factors="<<L.size()<<" pairs="<<n<<" failing="<<fail<<"\n";printM("IP-WITNESS-MIN",best);}}
int main(int argc,char**argv){if(argc>1&&(!strcmp(argv[1],"-h")||!strcmp(argv[1],"--help"))){cout<<"usage: ip-census LOG [N=2000] [SEED=182] | ip-runs [MAXK=40] | ip-42 | ip-witness\n";return 0;}try{omp_set_num_threads(4);if(argc>1&&!strcmp(argv[1],"ip-census")){runCensus(argv[2],argc>3?atoi(argv[3]):2000,argc>4?strtoull(argv[4],nullptr,10):182);return 0;}if(argc>1&&!strcmp(argv[1],"ip-runs")){runRuns(argc>2?atoi(argv[2]):40);return 0;}if(argc>1&&!strcmp(argv[1],"ip-42")){run42();return 0;}if(argc>1&&!strcmp(argv[1],"ip-witness")){runWitness();return 0;}cerr<<"use --help\n";return 2;}catch(exception&e){cerr<<"ERROR "<<e.what()<<"\n";return 3;}}"""
if any(a in ("-h","--help") for a in sys.argv[1:]):
 print("usage: python3 -u verifier.py ip-census LOG [N=2000] [SEED=182] | ip-runs [MAXK=40] | ip-42 | ip-witness");raise SystemExit(0)
fd=os.memfd_create("fm_sec182_ip",0);os.set_inheritable(fd,True)
env=dict(os.environ,TMPDIR="/dev/shm",OMP_NUM_THREADS="4")
b=subprocess.run(["g++","-std=c++17","-O3","-fopenmp","-x","c++","-o",f"/proc/self/fd/{fd}","-","-lgmpxx","-lgmp"],input=CPP.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(fd,),env=env)
print("compile_rc",b.returncode,b.stderr.decode(),flush=True)
if b.returncode:raise SystemExit(b.returncode)
os.fchmod(fd,0o700)
binary=f"/proc/self/fd/{fd}"
if len(sys.argv)>1:
 raise SystemExit(subprocess.run([binary]+sys.argv[1:],pass_fds=(fd,),env=env).returncode)
log="/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/fx6_52_noflip.log"
for args in (["ip-witness"],["ip-census",log,"2000","182"],["ip-runs","40"],["ip-42"]):
 print("RUN",*args,flush=True)
 rc=subprocess.run([binary]+args,pass_fds=(fd,),env=env).returncode
 if rc:raise SystemExit(rc)
