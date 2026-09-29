// Built and used by verify_qfm3_r2.py (FM3 note, section "Euler-characteristic mechanism", 2026-09-27).
// Enumerate SSYT with content kappa and <=4 rows; output Kostka-Foulkes K_{lam,kappa}(t) (charge) per shape lam.
// usage: kfc k1,k2,...   (kappa any composition; charge needs partition content -> we sort kappa decreasing)
#include <bits/stdc++.h>
using namespace std;
int N; vector<int> kap;
map<array<int,4>, map<int,long long>> acc;
vector<array<int,4>> chain; // shapes after each letter
// tableau: rows as vectors of letters
int charge_of(const vector<int>& word){
  // Lascoux-Schutzenberger charge for word with partition content
  vector<int> w=word; vector<char> used(w.size(),0); int total=0; int L=w.size(); int remaining=L;
  while(remaining>0){
    // standard subword extraction: rightmost unused 1, then cyclically leftward 2,3,...
    int start=-1; for(int p=L-1;p>=0;p--) if(!used[p]&&w[p]==1){start=p;break;}
    vector<int> sel; sel.push_back(start); used[start]=1; int cur=start; int k=1;
    while(true){ k++; int found=-1; for(int s=1;s<=L;s++){int p=((cur-s)%L+L)%L; if(!used[p]&&w[p]==k){found=p;break;}} if(found<0)break; sel.push_back(found); used[found]=1; cur=found; }
    // charge of standard subword in original order: index(1)=0; if pos(k+1)>pos(k) index++
    int idx=0; for(size_t j=1;j<sel.size();j++){ if(sel[j]>sel[j-1]) idx++; total+=idx; }
    remaining-=sel.size();
  }
  return total;
}
vector<vector<int>> rows(4);
void rec(int letter, array<int,4> sh){
  if(letter==(int)kap.size()){
    vector<int> word; for(int i=3;i>=0;i--) for(int x: rows[i]) word.push_back(x);
    acc[sh][charge_of(word)]++; return;
  }
  int c=kap[letter];
  // add horizontal strip of size c: new row lengths nl[i] with sh[i]<=nl[i]<= (i==0? inf : sh[i-1])
  array<int,4> nl=sh;
  function<void(int,int)> go=[&](int i,int rem){
    if(i==4){ if(rem==0){ for(int r=0;r<4;r++) for(int j=sh[r];j<nl[r];j++) rows[r].push_back(letter+1); rec(letter+1,nl); for(int r=0;r<4;r++) rows[r].resize(sh[r]); } return; }
    int hi = (i==0)? sh[0]+rem : min(sh[i]+rem, sh[i-1]);
    for(int v=sh[i]; v<=hi; v++){ nl[i]=v; go(i+1, rem-(v-sh[i])); }
    nl[i]=sh[i];
  };
  go(0,c);
}
int main(int argc,char**argv){
  string s=argv[1]; stringstream ss(s); string tkn; while(getline(ss,tkn,',')) if(!tkn.empty()) kap.push_back(stoi(tkn));
  sort(kap.rbegin(),kap.rend());
  rec(0,{0,0,0,0});
  for(auto &e: acc){ printf("%d,%d,%d,%d:",e.first[0],e.first[1],e.first[2],e.first[3]); for(auto &f: e.second) printf(" %d:%lld",f.first,f.second); printf("\n"); }
}
