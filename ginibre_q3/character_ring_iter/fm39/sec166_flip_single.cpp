// Flip-descent and TopPair check for one (B,p).  Usage: flip_single p B...
// For each pair of factor classes (u,v) of Lambda = B + (+p): D_uv = eps_v * A_{n_u n_v}(Lambda - u - v)
// (Phi(Lambda) - Phi(Lambda^{uv}) = 4 D_uv).  Flip descent holds iff some D_uv >= 0.
#include <gmpxx.h>
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using Z=mpz_class;
static Z entry(vector<int> lab,int A,int Bt){ // A_{A,Bt}(lab) with pruning
  sort(lab.begin(),lab.end(),[](int a,int b){return abs(a)>abs(b);});
  long rem=0;for(int x:lab)rem+=abs(x);
  auto key=[](int s,int t){return ((long long)s<<32)|(unsigned)t;};
  unordered_map<long long,Z> cur,nx;cur[key(0,0)]=1;
  for(int x:lab){int n=abs(x),e=x>0?1:-1;rem-=n;nx.clear();
    for(auto&kv:cur){int s=(int)(kv.first>>32),t=(int)(kv.first&0xffffffff);const Z&v=kv.second;
      if(abs(t-Bt)<=rem)for(int s2=abs(s-n);s2<=s+n;s2+=2)if(abs(s2-A)<=rem)nx[key(s2,t)]+=v;
      if(abs(s-A)<=rem)for(int t2=abs(t-n);t2<=t+n;t2+=2)if(abs(t2-Bt)<=rem){if(e>0)nx[key(s,t2)]+=v;else nx[key(s,t2)]-=v;}}
    cur.swap(nx);}
  auto it=cur.find(key(A,Bt));return it==cur.end()?Z(0):it->second;}
int main(int argc,char**argv){int p=atoi(argv[1]);vector<int> B;for(int i=2;i<argc;++i)B.push_back(atoi(argv[i]));
  vector<int> Lam=B;Lam.push_back(p);
  vector<int> vals;for(int z:Lam)if(find(vals.begin(),vals.end(),z)==vals.end())vals.push_back(z);
  vector<pair<int,int>> prs;for(size_t i=0;i<vals.size();++i)for(size_t j=i;j<vals.size();++j){
    if(i==j&&count(Lam.begin(),Lam.end(),vals[i])<2)continue;prs.push_back({(int)i,(int)j});}
  printf("classes=%zu pairs=%zu\n",vals.size(),prs.size());fflush(stdout);
  int found=0;vector<string> out(prs.size());
  #pragma omp parallel for schedule(dynamic,1)
  for(size_t k=0;k<prs.size();++k){
    #pragma omp flush(found)
    if(found)continue;
    int a=vals[prs[k].first],b=vals[prs[k].second];vector<int> C=Lam;C.erase(find(C.begin(),C.end(),a));C.erase(find(C.begin(),C.end(),b));
    auto t0=chrono::steady_clock::now();Z v=entry(C,abs(a),abs(b));if(b<0)v=-v;
    {double el=chrono::duration<double>(chrono::steady_clock::now()-t0).count();time_t tt=time(nullptr);char buf[16];strftime(buf,16,"%H:%M:%S",localtime(&tt));
     #pragma omp critical
     {static int done=0;++done;printf("[%s] pair %d/%zu (%+d,%+d) cost %.1fs sign %s\n",buf,done,prs.size(),a,b,el,v>=0?"+":"-");fflush(stdout);}}
    if(v>=0){
      #pragma omp critical
      {found=1;printf("FLIP DESCENT at (%+d,%+d): eps_v*A = %s\n",a,b,v.get_str().c_str());fflush(stdout);}
    }}
  if(!found)printf("NO flip descent among %zu pairs\n",prs.size());
}
