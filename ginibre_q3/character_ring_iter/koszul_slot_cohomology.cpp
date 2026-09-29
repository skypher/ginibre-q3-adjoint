// koszul_slot_cohomology: usage ./koszul_slot_cohomology r "k1,k2,..." ["c00,c01,..;c10,.."]  (default C = Vandermonde (i+2)^j).
// H-invariant Koszul cohomology of a = n (x) C^r acting on M = (x)_i Sym^{k_i}(C^4)
// n = Hom(B,A), A=span(e0,e1), B=span(e2,e3), H = SL(A) x SL(B).
// rho(E_e (x) f_j) = sum_i C[j][i] * (slot-i action of E_e).
// Invariant multiplicity via weight inclusion-exclusion: m(0,0)-m(2,0)-m(0,2)+m(2,2).
// Complex splits by s = (A-degree of output) - (cochain degree); computed per block.
// Ranks mod prime by sparse left-looking elimination.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef unsigned long long ull;
static const ll P = 2147483629LL; // prime < 2^31
ll pw(ll a, ll e){ll r=1;a%=P;if(a<0)a+=P;while(e){if(e&1)r=r*a%P;a=a*a%P;e>>=1;}return r;}
ll inv(ll a){return pw(a,P-2);}

struct SparseRank {
  // pivots keyed by leading column
  unordered_map<int, vector<pair<int,ll>>> piv;
  int rank=0;
  void add(vector<pair<int,ll>> row){ // row sorted by col, nonzero vals
    while(!row.empty()){
      int lc=row[0].first;
      auto it=piv.find(lc);
      if(it==piv.end()){
        ll iv=inv(row[0].second);
        for(auto &x:row) x.second=x.second*iv%P;
        piv.emplace(lc,move(row)); rank++; return;
      }
      const auto &pr=it->second; ll f=row[0].second; // pr leading coef is 1
      vector<pair<int,ll>> out; out.reserve(row.size()+pr.size());
      size_t a=0,b=0;
      while(a<row.size()||b<pr.size()){
        if(b>=pr.size()||(a<row.size()&&row[a].first<pr[b].first)){out.push_back(row[a]);a++;}
        else if(a>=row.size()||pr[b].first<row[a].first){ll v=(P-f*pr[b].second%P)%P; if(v) out.push_back({pr[b].first,v}); b++;}
        else {ll v=(row[a].second - f*pr[b].second)%P; if(v<0)v+=P; if(v) out.push_back({row[a].first,v}); a++;b++;}
      }
      row.swap(out);
    }
  }
};

int main(int argc,char**argv){
  // usage: kz r "k1,k2,..." "c00,c01,...;c10,..."
  int r=atoi(argv[1]);
  vector<int> kap; {string s=argv[2]; stringstream ss(s); string t; while(getline(ss,t,',')) if(!t.empty()) kap.push_back(stoi(t));}
  int n=kap.size();
  vector<vector<ll>> C(r, vector<ll>(n,0));
  {string s=argc>3?argv[3]:""; if(s.empty()){ for(int j=0;j<r;j++) for(int i=0;i<n;i++) C[j][i]=pw(i+2,j);} else { stringstream ss(s); string row; int j=0; while(getline(ss,row,';')){ stringstream rs(row); string t; int i=0; while(getline(rs,t,',')){C[j][i++]=stoll(t);} j++;}}}
  // slot monomials
  vector<vector<array<int,4>>> mon(n); vector<map<array<int,4>,int>> mid(n);
  for(int i=0;i<n;i++){int k=kap[i]; for(int a=0;a<=k;a++)for(int b=0;a+b<=k;b++)for(int c=0;a+b+c<=k;c++){array<int,4> m={a,b,c,k-a-b-c}; mid[i][m]=mon[i].size(); mon[i].push_back(m);} }
  vector<ll> radix(n+1,1); for(int i=0;i<n;i++) radix[i+1]=radix[i]*mon[i].size();
  ll Mdim=radix[n];
  // weights of M basis: (wA, wB, adeg)
  vector<int> MwA(Mdim),MwB(Mdim),Mad(Mdim);
  for(ll x=0;x<Mdim;x++){int wa=0,wb=0,ad=0; ll y=x; for(int i=0;i<n;i++){auto &m=mon[i][y%mon[i].size()]; y/=mon[i].size(); wa+=m[0]-m[1]; wb+=m[2]-m[3]; ad+=m[0]+m[1];} MwA[x]=wa;MwB[x]=wb;Mad[x]=ad;}
  // a basis: idx = j*4+e ; E_e: (a,b) pairs
  int NB[4][2]={{0,2},{0,3},{1,2},{1,3}};
  int WA[4]={1,-1,0,0}, WB[4]={0,0,1,-1};
  int dA=4*r;
  vector<int> awA(dA),awB(dA);
  for(int j=0;j<r;j++)for(int e=0;e<4;e++){int a=NB[e][0],b=NB[e][1]; awA[j*4+e]=WA[a]-WA[b]; awB[j*4+e]=WB[a]-WB[b];}
  // index M by (wA,wB,ad)
  map<array<int,3>, vector<ll>> Mby;
  for(ll x=0;x<Mdim;x++) Mby[{MwA[x],MwB[x],Mad[x]}].push_back(x);
  int tws[4][2]={{0,0},{2,0},{0,2},{2,2}};
  int maxad=0; for(int k:kap) maxad+=k;
  vector<vector<ll>> Hres(4, vector<ll>(dA+1,0));
  for(int t=0;t<4;t++){
    int twA=tws[t][0], twB=tws[t][1];
    // for each s block: degree k basis = (mask,x) with |mask|=k, Mad[x]-k = s, wt(x)-wt(mask)=tw
    for(int s=-dA; s<=maxad; s++){
      vector<vector<pair<unsigned,ll>>> basis(dA+1);
      for(unsigned mask=0; mask<(1u<<dA); mask++){
        int k=__builtin_popcount(mask); int ad=s+k; if(ad<0||ad>maxad) continue;
        int mA=0,mB=0; for(int q=0;q<dA;q++) if(mask>>q&1){mA+=awA[q];mB+=awB[q];}
        auto it=Mby.find({twA+mA,twB+mB,ad}); if(it==Mby.end()) continue;
        for(ll x: it->second) basis[k].push_back({mask,x});
      }
      bool any=false; for(auto&b:basis) if(!b.empty()) any=true; if(!any) continue;
      vector<ll> rk(dA+1,0);
      for(int k=0;k<dA;k++){
        if(basis[k].empty()||basis[k+1].empty()) continue;
        // target index
        unordered_map<ull,int> tidx; tidx.reserve(basis[k+1].size()*2);
        for(size_t i=0;i<basis[k+1].size();i++) tidx[((ull)basis[k+1][i].first<<40)|(ull)basis[k+1][i].second]=i;
        SparseRank SR;
        for(auto &bb: basis[k]){
          unsigned mask=bb.first; ll x=bb.second;
          map<int,ll> col;
          for(int q=0;q<dA;q++){ if(mask>>q&1) continue; unsigned T=mask|(1u<<q); int pos=__builtin_popcount(mask&((1u<<q)-1)); ll sg=(pos&1)?-1:1;
            int j=q/4,e=q%4; int a=NB[e][0], b=NB[e][1];
            ll y=x;
            for(int i=0;i<n;i++){ int ms=mon[i].size(); int mi=(y/radix[i])%ms; ll cc=C[j][i]; if(cc==0) continue; auto m=mon[i][mi]; if(m[b]==0) continue; ll co=m[b]; m[b]--; m[a]++; int mi2=mid[i].at(m); ll x2=x+(ll)(mi2-mi)*radix[i];
              auto it=tidx.find(((ull)T<<40)|(ull)x2); if(it==tidx.end()){fprintf(stderr,"missing target\n"); exit(1);}
              ll v=((sg*cc%P)*co)%P; if(v<0)v+=P; col[it->second]=(col[it->second]+v)%P; }
          }
          vector<pair<int,ll>> row; for(auto&z:col) if(z.second) row.push_back(z);
          SR.add(move(row));
        }
        rk[k]=SR.rank;
      }
      for(int k=0;k<=dA;k++){ ll h=(ll)basis[k].size()-rk[k]-(k>0?rk[k-1]:0); Hres[t][k]+=h; }
    }
  }
  vector<ll> invH(dA+1); ll eu=0; bool odd=false;
  for(int k=0;k<=dA;k++){invH[k]=Hres[0][k]-Hres[1][k]-Hres[2][k]+Hres[3][k]; eu+=(k%2?-1:1)*invH[k]; if(k%2&&invH[k]) odd=true;}
  printf("r=%d kappa=(",r); for(int i=0;i<n;i++) printf("%d%s",kap[i],i+1<n?",":""); printf(") H_inv=[");
  for(int k=0;k<=dA;k++) printf("%lld%s",invH[k],k<dA?",":""); printf("] euler=%lld %s\n",eu,odd?"ODD":"");
  fflush(stdout);
}
