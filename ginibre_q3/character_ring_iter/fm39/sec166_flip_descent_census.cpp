// FM-SEC166 flip-descent census (parallel).  Usage: flipx3 MAXW [heartbeat_seconds]
// For every pair-free residual (B,p) with weight(B) <= MAXW (p >= max(6, max|B|), delta=(W-p)/2 >= 8,
// max|B| <= delta, >= 2 factors >= 3), Lambda = B + (sigma p), sigma = (-1)^(minus count of B) (FM3(Lambda) = 2 g_p(B));
// lists where B holds -sigma p are skipped (they contain a pair):
//  flip descent at factor classes (u,v):  eps_v * A_{n_u n_v}(Lambda - u - v) >= 0  [Phi(Lambda) - Phi(Lambda^{uv}) = 4 * that].
// For lists with no flip descent, tests (D) over all even-weight removals and the rules
// TopPair (two same-parity factors of B with the largest total label; ties -> larger max label) and TopEven (largest even label).
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using I=__int128_t;
static int LIM; static double HB=30;
static string si(I x){if(x==0)return "0";bool neg=x<0;if(neg)x=-x;string s;while(x){s.push_back('0'+(int)(x%10));x/=10;}if(neg)s.push_back('-');reverse(s.begin(),s.end());return s;}
static void table(const vector<int>& L,int S,vector<I>& d){ // d[(a)*(S+1)+b]
  int N=S+1;d.assign(N*N,0);vector<I> q(N*N);d[0]=1;int sm=0;
  for(int z:L){int n=abs(z),e=z>0?1:-1;fill(q.begin(),q.end(),0);
    for(int a=0;a<=sm;++a)for(int b=0;b<=sm-a;++b){I v=d[a*N+b];if(!v)continue;
      for(int c=abs(a-n);c<=a+n;c+=2)q[c*N+b]+=v;
      for(int c=abs(b-n);c<=b+n;c+=2)q[a*N+c]+=e*v;}
    d.swap(q);sm+=n;}}
static vector<array<signed char,41>> P; static array<signed char,41> cur{};
static void gen(int n,int w){if(n>LIM){P.push_back(cur);return;}for(int c=-LIM/n;c<=LIM/n;++c){if(w+n*abs(c)>LIM)continue;cur[n]=c;gen(n+1,w+n*abs(c));}cur[n]=0;}
int main(int argc,char**argv){LIM=atoi(argv[1]);if(argc>2)HB=atof(argv[2]);if(LIM>40)return 2;
  gen(1,0);size_t NP=P.size();
  atomic<unsigned long long> cases{0},nof{0},dfail{0},tpfail{0},tefail{0},both{0},pdone{0};
  auto T0=chrono::steady_clock::now();auto Tlast=T0;mutex mu;
  fprintf(stderr,"profiles=%zu\n",NP);
  #pragma omp parallel for schedule(dynamic,64)
  for(size_t k=0;k<NP;++k){
    const auto& c=P[k];int W=0,mx=0,cores=0,neg=0;vector<int> B;
    for(int n=1;n<=LIM;++n)if(c[n]){int m=abs(c[n]);W+=n*m;mx=n;if(n>=3)cores+=m;if(c[n]<0)neg+=m;for(int t=0;t<m;++t)B.push_back(c[n]>0?n:-n);}
    if(cores>=2){ int sg=(neg%2)?-1:1;
     vector<I> A;
     for(int p=max(6,mx);p<=W;++p){if((W-p)%2)continue;int delta=(W-p)/2;if(delta<8||mx>delta)continue;if(p<=LIM&&c[p]*sg<0)continue;
      ++cases;vector<int> Lam=B;Lam.push_back(sg*p);int T=W+p;
      vector<int> vals;for(int z:Lam)if(find(vals.begin(),vals.end(),z)==vals.end())vals.push_back(z);
      bool desc=false;
      for(size_t i=0;i<vals.size()&&!desc;++i)for(size_t j=i;j<vals.size()&&!desc;++j){
        vector<int> C=Lam;C.erase(find(C.begin(),C.end(),vals[i]));auto it=find(C.begin(),C.end(),vals[j]);if(it==C.end())continue;C.erase(it);
        int a=abs(vals[i]),b=abs(vals[j]);if(a+b>T)continue;table(C,T,A);I v=A[a*(T+1)+b];if((vals[j]>0?v:-v)>=0)desc=true;}
      if(desc)continue;
      ++nof;
      table(B,W,A);I gB=A[p*(W+1)+0];
      int tpA=0,tpB=0,tpw=-1;for(size_t i=0;i<B.size();++i)for(size_t j=i+1;j<B.size();++j)if((abs(B[i])+abs(B[j]))%2==0){int w=abs(B[i])+abs(B[j]);if(w>tpw||(w==tpw&&max(abs(B[i]),abs(B[j]))>max(abs(tpA),abs(tpB)))){tpw=w;tpA=B[i];tpB=B[j];}}
      int te=0;for(int z:B)if(abs(z)%2==0&&abs(z)>abs(te))te=z;
      vector<int> bv;for(int z:B)if(find(bv.begin(),bv.end(),z)==bv.end())bv.push_back(z);
      bool ok=false,tp=false,teok=false;string mono;
      for(size_t i=0;i<bv.size();++i)for(size_t j=i;j<=bv.size();++j){
        vector<int> C=B;C.erase(find(C.begin(),C.end(),bv[i]));int wr=abs(bv[i]);string R=to_string(bv[i]);
        if(j<bv.size()){auto it=find(C.begin(),C.end(),bv[j]);if(it==C.end())continue;C.erase(it);wr+=abs(bv[j]);R+=","+to_string(bv[j]);}
        if(wr%2)continue;int Wc=W-wr;I g=0;if(p<=Wc){table(C,Wc,A);g=A[p*(Wc+1)];}
        if(g<=gB){ok=true;mono+="("+R+")";if(j<bv.size()&&((bv[i]==tpA&&bv[j]==tpB)||(bv[i]==tpB&&bv[j]==tpA)))tp=true;if(j==bv.size()&&bv[i]==te)teok=true;}}
      if(!ok)++dfail;if(!tp)++tpfail;if(!teok)++tefail;if(!tp&&!teok)++both;
      string line="NOFLIP W="+to_string(W)+" p="+to_string(sg*p)+" B=";for(int z:B)line+=to_string(z)+" ";
      line+=" phi="+si(2*gB)+" D="+(ok?"ok":"FAIL")+" TopPair("+to_string(tpA)+","+to_string(tpB)+")="+(tp?"ok":"FAIL")+" TopEven("+to_string(te)+")="+(teok?"ok":"FAIL")+" mono="+mono+"\n";
      {lock_guard<mutex> g(mu);fputs(line.c_str(),stdout);fflush(stdout);}
     }}
    unsigned long long d=++pdone;
    if(d%256==0){lock_guard<mutex> g(mu);auto now=chrono::steady_clock::now();if(chrono::duration<double>(now-Tlast).count()>HB){Tlast=now;time_t tt=time(nullptr);char buf[16];strftime(buf,16,"%H:%M:%S",localtime(&tt));double el=chrono::duration<double>(now-T0).count();
      fprintf(stderr,"[hb %s] profiles %llu/%zu (%.1f%%) residual_cases=%llu no_flip=%llu D_fail=%llu TopPair_fail=%llu TopEven_fail=%llu both_fail=%llu elapsed=%.0fs avg_us_per_case=%.0f\n",buf,d,NP,100.0*d/NP,(unsigned long long)cases,(unsigned long long)nof,(unsigned long long)dfail,(unsigned long long)tpfail,(unsigned long long)tefail,(unsigned long long)both,el,1e6*el*omp_get_num_threads()/max(1ULL,(unsigned long long)cases));fflush(stderr);}}
  }
  double el=chrono::duration<double>(chrono::steady_clock::now()-T0).count();
  printf("SUMMARY W<=%d profiles=%zu residual_cases=%llu no_flip_descent=%llu D_fail=%llu TopPair_fail=%llu TopEven_fail=%llu both_fail=%llu time_s=%.1f\n",LIM,NP,(unsigned long long)cases,(unsigned long long)nof,(unsigned long long)dfail,(unsigned long long)tpfail,(unsigned long long)tefail,(unsigned long long)both,el);}
