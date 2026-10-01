// Random search for TopPair failures, exact GMP arithmetic.  Usage: tpsearch_g seed seconds Wmin Wmax maxones
// maxones: generate lists with at most this many factors of label 1 (use 1 to test (T1)); -1 = unrestricted.
#include <gmpxx.h>
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using Z=mpz_class;
static Z gp(vector<int> lab,int p){sort(lab.begin(),lab.end(),[](int a,int b){return abs(a)>abs(b);});long rem=0;for(int x:lab)rem+=abs(x);
  auto key=[](int s,int t){return ((long long)s<<32)|(unsigned)t;};unordered_map<long long,Z> cur,nx;cur[key(0,0)]=1;
  for(int x:lab){int n=abs(x),e=x>0?1:-1;rem-=n;nx.clear();
    for(auto&kv:cur){int s=(int)(kv.first>>32),t=(int)(kv.first&0xffffffff);const Z&v=kv.second;
      if(t<=rem)for(int s2=abs(s-n);s2<=s+n;s2+=2)if(abs(s2-p)<=rem)nx[key(s2,t)]+=v;
      if(abs(s-p)<=rem)for(int t2=abs(t-n);t2<=t+n;t2+=2)if(t2<=rem){if(e>0)nx[key(s,t2)]+=v;else nx[key(s,t2)]-=v;}}
    cur.swap(nx);}
  auto it=cur.find(key(p,0));return it==cur.end()?Z(0):it->second;}
int main(int argc,char**argv){int seed=atoi(argv[1]);double secs=atof(argv[2]);int Wmin=atoi(argv[3]),Wmax=atoi(argv[4]),maxones=atoi(argv[5]);
  double best=0;long tested=0,fails=0;map<int,long> failOnes;auto T0=chrono::steady_clock::now();auto Tl=T0;mutex mu;
  #pragma omp parallel
  {mt19937_64 rng(seed*1000+omp_get_thread_num());
   while(chrono::duration<double>(chrono::steady_clock::now()-T0).count()<secs){
    int W=Wmin+rng()%(Wmax-Wmin+1);int smallmax=2+rng()%3;int nmed=rng()%4;vector<int> lab;int w=0;map<int,int> sg;
    int medmax=max(4,(int)(W*(0.05+0.25*(rng()%1000)/1000.0)));
    for(int k=0;k<nmed;++k){int n=4+rng()%max(1,medmax-3);lab.push_back(n);w+=n;}
    int target=W-(int)(W*(0.1+0.5*(rng()%1000)/1000.0));int ones=0;
    while(w<target){int n=1+rng()%smallmax;if(n==1&&maxones>=0&&ones>=maxones)n=2+rng()%smallmax;if(n==1)++ones;lab.push_back(n);w+=n;}
    for(int& n:lab){if(!sg.count(n))sg[n]=(rng()%2)?1:-1;n*=sg[n];}
    int Wb=0,mx=0,cores=0;for(int z:lab){Wb+=abs(z);mx=max(mx,abs(z));if(abs(z)>=3)++cores;}if(cores<2)continue;
    int pmin=max(6,mx),pmax=Wb-16;if(pmin>pmax)continue;int p=pmin+rng()%(pmax-pmin+1);if((Wb-p)%2)++p;if(p>pmax)continue;
    int delta=(Wb-p)/2;if(delta<8||mx>delta)continue;int neg=0;for(int z:lab)if(z<0)++neg;int sgp=neg%2?-1:1;if(sg.count(p)&&sg[p]*sgp<0)continue;
    int a=0,b=0,tw=-1;for(size_t i=0;i<lab.size();++i)for(size_t j=i+1;j<lab.size();++j)if((abs(lab[i])+abs(lab[j]))%2==0){int ww=abs(lab[i])+abs(lab[j]);if(ww>tw||(ww==tw&&max(abs(lab[i]),abs(lab[j]))>max(abs(a),abs(b)))){tw=ww;a=lab[i];b=lab[j];}}
    if(tw<0)continue;Z g0=gp(lab,p);vector<int> C=lab;C.erase(find(C.begin(),C.end(),a));C.erase(find(C.begin(),C.end(),b));Z g1=(p<=Wb-tw)?gp(C,p):Z(0);
    {lock_guard<mutex> g(mu);++tested;if(g1>g0){++fails;failOnes[ones]++;double r=(double)tw/delta;if(r>best)best=r;
       if(fails<=20||ones<=1){printf("TPFAIL ones=%d ratio=%.3f W=%d p=%d delta=%d wTP=%d N=%zu B=",ones,r,Wb,sgp*p,delta,tw,lab.size()+1);for(int z:lab)printf("%d ",z);printf("\n");fflush(stdout);}}
     auto now=chrono::steady_clock::now();if(chrono::duration<double>(now-Tl).count()>30){Tl=now;time_t tt=time(nullptr);char buf[16];strftime(buf,16,"%H:%M:%S",localtime(&tt));fprintf(stderr,"[hb %s] tested=%ld TopPair_fails=%ld max_fail_ratio=%.3f elapsed=%.0f/%.0fs\n",buf,tested,fails,best,chrono::duration<double>(now-T0).count(),secs);}}
   }}
  printf("DONE tested=%ld fails=%ld max_fail_ratio=%.3f fails_by_ones:",tested,fails,best);for(auto&q:failOnes)printf(" %d:%ld",q.first,q.second);printf("\n");}
