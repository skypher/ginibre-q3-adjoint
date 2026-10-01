// Few-factor exhaustive (FT) census via Walsh transforms.
// Usage: walshft N Tmax [heartbeat_s]
// Enumerates every multiset of N positive labels n_1<=...<=n_N with total T <= Tmax, distinguished p = n_N,
// background B = the other N-1 labels (W = T - p): residual iff p >= 6, (W-p) even, delta = (W-p)/2 >= 8,
// max(B) <= delta, at least two labels >= 3 in B.  For every even sign pattern M that is pair-free on the
// whole list (equal labels share a sign), Phi(M) = sum_S (-1)^{|S cap M|} m(S) m(S^c).
// No flip descent at M  <=>  Phi(M) < Phi(M xor {u,v}) for every pair u<v.
// At such M: test TopPair (equal-parity pair of B with largest total label; ties: larger max label) and (D)
// (all even-weight removals of one or two factors from B), comparing Phi(child) with p's sign re-chosen by parity.
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; typedef long long ll; typedef __int128 I;
static string si(I x){if(x==0)return "0";bool n=x<0;if(n)x=-x;string s;while(x){s.push_back('0'+(int)(x%10));x/=10;}if(n)s.push_back('-');reverse(s.begin(),s.end());return s;}
// m(S) for all subsets of labels
static void mults(const vector<int>& L, vector<ll>& m){int N=L.size(),T=0;for(int x:L)T+=x;
  vector<vector<ll>> F(1<<N);F[0].assign(T+1,0);F[0][0]=1;m.assign(1<<N,0);m[0]=1;
  for(int S=1;S<(1<<N);++S){int i=__builtin_ctz(S);const auto& G=F[S&(S-1)];int n=L[i];vector<ll> H(T+1,0);
    for(int a=0;a<=T;++a)if(G[a])for(int c=abs(a-n);c<=a+n&&c<=T;c+=2)H[c]+=G[a];F[S]=move(H);m[S]=F[S][0];}}
static void walsh(vector<I>& v){int n=v.size();for(int h=1;h<n;h<<=1)for(int s=0;s<n;s+=2*h)for(int j=s;j<s+h;++j){I a=v[j],b=v[j+h];v[j]=a+b;v[j+h]=a-b;}}
static map<string,map<int,ll>> ST;
static I phiAt(const vector<int>& L,int M){vector<ll> m;mults(L,m);int N=L.size(),full=(1<<N)-1;I r=0;
  for(int S=0;S<=full;++S){I f=(I)m[S]*m[full^S];if(!f)continue;r+=(__builtin_popcount(S&M)&1)?-f:f;}return r;}
int main(int argc,char**argv){int N=atoi(argv[1]),Tmax=atoi(argv[2]);double HB=argc>3?atof(argv[3]):30;int Wmax=argc>4?atoi(argv[4]):1000000;
  // enumerate multisets (nondecreasing sequences) with sum <= Tmax
  vector<int> cur(N);
  vector<vector<int>> R; size_t nms=0;
  auto keep=[&](const vector<int>& L){int T=0;for(int x:L)T+=x;int p=L[N-1],W=T-p;if(p<6||(W-p)%2||W-p<16||W>Wmax)return false;int d=(W-p)/2;if(L[N-2]>d)return false;int c3=0;for(int i=0;i<N-1;++i)if(L[i]>=3)++c3;return c3>=2;};
  function<void(int,int,int)> rec=[&](int i,int lo,int sum){if(i==N){++nms;if(keep(cur))R.push_back(cur);return;}
    for(int x=lo;sum+x*(N-i)<=Tmax;++x){cur[i]=x;rec(i+1,x,sum+x);}};
  rec(0,1,0);
  // filter residual
  fprintf(stderr,"N=%d Tmax=%d multisets=%zu residual_multisets=%zu\n",N,Tmax,nms,R.size());
  atomic<ll> done{0},pats{0},noflip{0},tpfail{0},dfail{0};mutex mu;auto T0=chrono::steady_clock::now();auto Tl=T0;
  #pragma omp parallel for schedule(dynamic,4)
  for(size_t k=0;k<R.size();++k){const auto& L=R[k];int full=(1<<N)-1;vector<ll> m;mults(L,m);
    vector<I> v(1<<N);for(int S=0;S<=full;++S)v[S]=(I)m[S]*m[full^S];walsh(v); // v[M]=Phi(M)
    // class masks for pair-free: positions with equal labels
    for(int M=0;M<=full;++M){if(__builtin_popcount(M)&1)continue;bool pf=true;for(int i=0;i+1<N&&pf;++i)if(L[i]==L[i+1]&&(((M>>i)^(M>>(i+1)))&1))pf=false;if(!pf)continue;
      ++pats;bool desc=false;for(int u=0;u<N&&!desc;++u)for(int w=u+1;w<N;++w)if(v[M^(1<<u)^(1<<w)]<=v[M]){desc=true;break;}
      if(desc)continue;++noflip;
      {int ones=0,mm=1,run=1;for(int i=0;i<N;++i){if(L[i]==1)++ones;if(i&&L[i]==L[i-1]){++run;mm=max(mm,run);}else run=1;}
       int nminus=__builtin_popcount(M);int nodd=0;for(int x:L)if(x&1)++nodd;
       lock_guard<mutex> g(mu);ST["ones"][ones]++;ST["maxmult"][mm]++;ST["minus"][nminus]++;ST["odd"][nodd]++;}
      int p=L[N-1];vector<int> B(L.begin(),L.end()-1);int Mb=M&((1<<(N-1))-1);
      // TopPair
      int ta=-1,tb=-1,tw=-1;for(int i=0;i<N-1;++i)for(int j=i+1;j<N-1;++j)if((B[i]+B[j])%2==0){int w=B[i]+B[j];if(w>tw||(w==tw&&max(B[i],B[j])>max(ta>=0?B[ta]:0,tb>=0?B[tb]:0))){tw=w;ta=i;tb=j;}}
      auto child=[&](vector<int> rm)->I{vector<int> CL;int CM=0,neg=0;for(int i=0;i<N-1;++i){if(find(rm.begin(),rm.end(),i)!=rm.end())continue;if(Mb>>i&1){CM|=1<<CL.size();++neg;}CL.push_back(B[i]);}
        int Wc=0;for(int x:CL)Wc+=x;if(p>Wc)return 0;if(neg&1)CM|=1<<CL.size();CL.push_back(p);return phiAt(CL,CM);};
      bool tp=true;if(ta>=0)tp=(child({ta,tb})<=v[M]);
      bool dok=false;for(int i=0;i<N-1&&!dok;++i){if(B[i]%2==0&&child({i})<=v[M])dok=true;for(int j=i+1;j<N-1&&!dok;++j)if((B[i]+B[j])%2==0&&child({i,j})<=v[M])dok=true;}
      if(!tp)++tpfail;if(!dok)++dfail;
      if(!tp||!dok){lock_guard<mutex> g(mu);printf("%s N=%d L=",(!dok?"DFAIL":"TPFAIL"),N);for(int i=0;i<N;++i)printf("%d ",(M>>i&1)?-L[i]:L[i]);printf(" phi=%s\n",si(v[M]).c_str());fflush(stdout);}
    }
    ll d=++done;if(d%64==0){lock_guard<mutex> g(mu);auto now=chrono::steady_clock::now();if(chrono::duration<double>(now-Tl).count()>HB){Tl=now;time_t tt=time(nullptr);char b[16];strftime(b,16,"%H:%M:%S",localtime(&tt));double el=chrono::duration<double>(now-T0).count();fprintf(stderr,"[hb %s] multisets %lld/%zu (%.1f%%) patterns=%lld noflip=%lld TopPair_fail=%lld D_fail=%lld elapsed=%.0fs per_multiset_ms=%.2f\n",b,d,R.size(),100.0*d/R.size(),(ll)pats,(ll)noflip,(ll)tpfail,(ll)dfail,el,1e3*el*omp_get_max_threads()/d);fflush(stderr);}}
  }
  for(auto&kv:ST){printf("STAT %s:",kv.first.c_str());for(auto&q:kv.second)printf(" %d:%lld",q.first,q.second);printf("\n");}
  printf("SUMMARY N=%d Tmax=%d residual_multisets=%zu pairfree_even_patterns=%lld noflip=%lld TopPair_fail=%lld D_fail=%lld time_s=%.1f\n",N,Tmax,R.size(),(ll)pats,(ll)noflip,(ll)tpfail,(ll)dfail,chrono::duration<double>(chrono::steady_clock::now()-T0).count());}
