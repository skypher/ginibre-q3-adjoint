// Exact no-flip test over ALL pair-free even sign patterns of a class multiset.
// Usage: manyflip p n1:m1 n2:m2 ...   (labels n_i with multiplicities m_i; distinguished label p, sign by parity)
// For each pattern (one sign per class, even total minus count incl. sigma p), computes every class-pair flip
// D = eps_v A_{n_u n_v}(Lambda - u - v) (pairs within a class need m >= 2; pairs with p included), GMP exact.
// Reports patterns with no flip descent and tests TopPair on them.
#include <gmpxx.h>
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using Z=mpz_class;
static Z entry(vector<int> lab,int A,int Bt){sort(lab.begin(),lab.end(),[](int a,int b){return abs(a)>abs(b);});long rem=0;for(int x:lab)rem+=abs(x);
  auto key=[](int s,int t){return ((long long)s<<32)|(unsigned)t;};unordered_map<long long,Z> cur,nx;cur[key(0,0)]=1;
  for(int x:lab){int n=abs(x),e=x>0?1:-1;rem-=n;nx.clear();for(auto&kv:cur){int s=(int)(kv.first>>32),t=(int)(kv.first&0xffffffff);const Z&v=kv.second;
      if(abs(t-Bt)<=rem)for(int s2=abs(s-n);s2<=s+n;s2+=2)if(abs(s2-A)<=rem)nx[key(s2,t)]+=v;
      if(abs(s-A)<=rem)for(int t2=abs(t-n);t2<=t+n;t2+=2)if(abs(t2-Bt)<=rem){if(e>0)nx[key(s,t2)]+=v;else nx[key(s,t2)]-=v;}}cur.swap(nx);}
  auto it=cur.find(key(A,Bt));return it==cur.end()?Z(0):it->second;}
int main(int argc,char**argv){int p=atoi(argv[1]);vector<int> n,m;for(int i=2;i<argc;++i){int a,b;sscanf(argv[i],"%d:%d",&a,&b);n.push_back(a);m.push_back(b);}
  int k=n.size();int W=0,mx=0,c3=0;for(int i=0;i<k;++i){W+=n[i]*m[i];mx=max(mx,n[i]);if(n[i]>=3)c3+=m[i];}
  int d=(W-p)/2;bool res=(p>=6&&p>=mx&&(W-p)%2==0&&d>=8&&mx<=d&&c3>=2);
  printf("profile W=%d p=%d delta=%d classes=%d factors=%d residual=%s\n",W,p,d,k,accumulate(m.begin(),m.end(),0)+1,res?"yes":"NO");
  if(!res)return 0;
  int noflip=0,pats=0,tpfail=0;
  for(int S=0;S<(1<<k);++S){ // S: classes with minus sign
    int neg=0;for(int i=0;i<k;++i)if(S>>i&1)neg+=m[i];int sp=(neg%2)?-1:1;
    bool haspair=false;for(int i=0;i<k;++i)if(n[i]==p&&((S>>i&1)?-1:1)!=sp)haspair=true;if(haspair)continue;
    ++pats;vector<int> Lam;for(int i=0;i<k;++i)for(int t=0;t<m[i];++t)Lam.push_back((S>>i&1)?-n[i]:n[i]);Lam.push_back(sp*p);
    vector<int> cls;for(int i=0;i<k;++i)cls.push_back((S>>i&1)?-n[i]:n[i]);cls.push_back(sp*p);vector<int> cm=m;cm.push_back(1);
    // merge p into its class if label equal
    vector<pair<int,int>> prs;for(int i=0;i<=k;++i)for(int j=i;j<=k;++j){if(i==j&&cm[i]<2)continue;prs.push_back({i,j});}
    vector<int> pos(prs.size(),0);
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t q=0;q<prs.size();++q){int a=cls[prs[q].first],b=cls[prs[q].second];vector<int> C=Lam;C.erase(find(C.begin(),C.end(),a));C.erase(find(C.begin(),C.end(),b));Z v=entry(C,abs(a),abs(b));if(b<0)v=-v;pos[q]=(v>=0);}
    bool desc=false;for(int x:pos)if(x)desc=true;
    if(!desc){++noflip;
      // TopPair: equal-parity pair of B with the largest total label (ties: larger max label); g_p via entry(.,p,0)
      vector<int> B(Lam.begin(),Lam.end()-1);int ta=-1,tb=-1,tw=-1;
      for(size_t i=0;i<B.size();++i)for(size_t j=i+1;j<B.size();++j)if((abs(B[i])+abs(B[j]))%2==0){int w=abs(B[i])+abs(B[j]);if(w>tw||(w==tw&&max(abs(B[i]),abs(B[j]))>max(ta>=0?abs(B[ta]):0,tb>=0?abs(B[tb]):0))){tw=w;ta=i;tb=j;}}
      Z g0=entry(B,p,0),g1=0;if(ta>=0){vector<int> C=B;int x=B[ta],y=B[tb];C.erase(find(C.begin(),C.end(),x));C.erase(find(C.begin(),C.end(),y));if(p<=W-tw)g1=entry(C,p,0);}
      bool tp=(ta>=0&&g1<=g0);if(!tp)++tpfail;
      printf("NOFLIP pattern minus-classes=");for(int i=0;i<k;++i)if(S>>i&1)printf("%d ",n[i]);printf(" sigma p=%d TopPair(%d,%d) %s g=%s child=%s\n",sp*p,ta>=0?B[ta]:0,tb>=0?B[tb]:0,tp?"ok":"FAIL",g0.get_str().c_str(),g1.get_str().c_str());}
  }
  printf("patterns=%d noflip=%d TopPair_fail=%d\n",pats,noflip,tpfail);}
