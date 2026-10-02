import argparse,ctypes,os,subprocess
from pathlib import Path
import re
ap=argparse.ArgumentParser(description="FM-MECH158 exact family certificates and optional no-flip census; memory only")
ap.add_argument("--threads",type=int,default=24)
ap.add_argument("--census",type=Path)
args=ap.parse_args()
assert args.threads>=1
source=r'''

#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cassert>
#include <omp.h>
using I=boost::multiprecision::checked_int128_t;
using V=std::vector<I>;
using namespace std;
V mults(const vector<int>&ns){
 int L=ns.size(),Z=1<<L,F=Z-1;
 vector<int>w(Z),sz(Z);
 for(int s=1;s<Z;s++){int b=s&-s,i=__builtin_ctz((unsigned)b);w[s]=w[s^b]+ns[i];sz[s]=sz[s^b]+1;}
 int cap=w[F]/2+L;vector<V> bin(L+1,V(cap+1));
 for(int q=0;q<=cap;q++)bin[0][q]=1;
 for(int d=1;d<=L;d++)for(int q=d;q<=cap;q++)bin[d][q]=bin[d][q-1]+bin[d-1][q-1];
 V out(Z);out[0]=1;
 for(int s=1;s<Z;s++){
  int l=sz[s];if(l<2||w[s]%2)continue;
  I v=0;for(int j=s;;j=(j-1)&s){
   int q=w[s]/2-w[j]-sz[j];if(q>=0)v+=(sz[j]%2?-bin[l-2][q+l-2]:bin[l-2][q+l-2]);
   if(j==0)break;
  }assert(v>=0);out[s]=v;
 }return out;
}
I unsigned_phi(const V&m,int selected){
 I v=0;for(int s=selected;;s=(s-1)&selected){v+=m[s]*m[selected^s];if(!s)break;}return v;
}
pair<int,int> top(const vector<int>&ns){
 pair<int,int>ij{-1,-1};pair<int,int>key{-1,-1};
 for(int i=0;i<(int)ns.size()-1;i++)for(int j=i+1;j<(int)ns.size()-1;j++)
 if((ns[i]-ns[j])%2==0){auto k=make_pair(ns[i]+ns[j],max(ns[i],ns[j]));if(k>key){key=k;ij={i,j};}}
 assert(ij.first>=0);return ij;
}
I margin(const vector<int>&ns){
 auto m=mults(ns);int F=(1<<ns.size())-1;auto[i,j]=top(ns);
 // 2*m(full) minus all proper split terms, minus unsigned child.
 I z=4*m[F]-unsigned_phi(m,F)-unsigned_phi(m,F^(1<<i)^(1<<j));
 for(int s=1;s<F;s++){int h=__builtin_popcount((unsigned)s);if(h==2||h==(int)ns.size()-2)z+=2*m[s]*m[F^s];}
 return z;
}
struct Profile{vector<int>off;int low;};
extern "C" int profiles_audit(int threads){
 omp_set_num_threads(threads);vector<Profile>profiles;set<pair<vector<int>,int>>seen;
 for(int L=6;L<=12;L++)for(int q=1;q<=2;q++)for(int tau=0;tau<=1;tau++)
 for(int j=-1;j<L;j++)for(int eps:{-2,2}){
  if(j==-1&&eps==2)continue;vector<int>off;
  for(int i=0;i<L;i++)off.push_back(tau+q*i+(i==j?eps:0));sort(off.begin(),off.end());
  int sum=0;for(int v:off)sum+=v;if(sum%2)continue;
  int lo=q==1?20:0;while(2*lo+off[0]<1)lo++;
  if(seen.insert({off,lo}).second)profiles.push_back({off,lo});
 }
 atomic<int>done{0},bad{0},good{0};atomic<long long>finite{0},coeffs{0};
 #pragma omp parallel for schedule(dynamic)
 for(int z=0;z<(int)profiles.size();z++){
  auto P=profiles[z];int L=P.off.size(),start=96;V values,certificate;bool ok=true;int first=-1;I smallest=0;
  int bound=2*L;for(int x:P.off)bound+=abs(x);assert(bound<2*start);
  for(int t=P.low;t<=start+L-2;t++){
   vector<int>ns=P.off;for(int&x:ns)x+=2*t;I value=margin(ns);
   if(t<start){finite++;if(value<0){ok=false;if(first<0){first=t;smallest=value;}}}
   else values.push_back(value);
  }
  int h=0;while(!values.empty()){
   coeffs++;certificate.push_back(values[0]);if(values[0]<0){ok=false;if(first<0){first=1000+h;smallest=values[0];}}
   V next;for(int j=1;j<(int)values.size();j++)next.push_back(values[j]-values[j-1]);
   values.swap(next);h++;
  }
  for(int u:{20,73}){
   vector<int>ns=P.off;for(int&x:ns)x+=2*(start+u);I value=0,choose=1;
   for(int h=0;h<(int)certificate.size();h++){
    if(h>0)choose=choose*(u-h+1)/h;value+=certificate[h]*choose;
   }assert(value==margin(ns));
  }
  if(!ok){bad++;
   #pragma omp critical
   {printf("BAD L=%d low=%d first=%d value=%s offsets=",L,P.low,first,smallest.str().c_str());for(int v:P.off)printf("%d,",v);printf("\n");fflush(stdout);}
  }else good++;
  int d=++done;if(d%25==0){
   #pragma omp critical
   {printf("profiles %d/%zu good=%d bad=%d\n",d,profiles.size(),good.load(),bad.load());fflush(stdout);}
  }
 }
 printf("PROFILES %zu good=%d bad=%d finite=%lld Newton=%lld\n",profiles.size(),good.load(),bad.load(),finite.load(),coeffs.load());fflush(stdout);
 return bad;
}


I signed_phi(const V&m,int selected,int neg){
 I ans=0;
 for(int s=selected;;s=(s-1)&selected){I z=m[s]*m[selected^s];ans+=(__builtin_popcount((unsigned)(s&neg))%2?-z:z);if(!s)break;}
 return ans;
}
extern "C" int bridges(){
 bool trapped=false;try{I x=1;for(int k=0;k<130;k++)x*=2;}catch(const std::overflow_error&){trapped=true;}assert(trapped);
 long long count=0;
 for(int sample=0;sample<4;sample++){
  vector<int>ns;for(int i=0;i<12;i++)ns.push_back(sample==0?i+1:sample==1?2*i+2:sample==2?(i==0?2:2*i):40+2*i);
  auto m=mults(ns);int Z=m.size();vector<V>rows(Z);rows[0]={I(1)};
  for(int s=1;s<Z;s++){
   int bit=s&-s,n=ns[__builtin_ctz((unsigned)bit)];auto &old=rows[s^bit];V row(old.size()+n+2);
   for(int k=0;k<(int)old.size();k++)if(old[k]!=0){row[abs(k-n)]+=old[k];row[k+n+2]-=old[k];}
   for(int k=2;k<(int)row.size();k++)row[k]+=row[k-2];
   row.resize(old.size()+n);rows[s]=move(row);
  }
  for(int s=0;s<Z;s++){assert(rows[s][0]==m[s]);count++;}
  if(sample==3){
   int F=Z-1,neg=5;I phi=signed_phi(m,F,neg),sum=0;
   for(auto ij:vector<pair<int,int>>{{9,10},{9,11},{10,11}}){
    I difference=phi-signed_phi(m,F,neg^(1<<ij.first)^(1<<ij.second));
    assert(difference%4==0);I d=difference/4;assert(d<0);sum+=d;
   }
   I child=signed_phi(m,F^(1<<9)^(1<<10),neg)/2;
   I T=(phi-sum)/2;
   assert(phi==I("906415068445728")&&child==I("179646349605"));
   assert(T==I("453312756673043")&&-sum/2==I("105222450179"));
   printf("FM157 triple: surplus=%s cost=%s\n",(T-child).str().c_str(),(-sum/2).str().c_str());
  }
 }
 printf("fusion/inclusion-exclusion identities=%lld; overflow guard PASS\n",count);fflush(stdout);return 0;
}

struct Table{
 int D;vector<I>v;
 I at(int i,int j)const{return i>=0&&j>=0&&i<D&&j<D?v[i*D+j]:I(0);}
};
Table table(const vector<int>&B,bool unsign){
 int W=0;for(int z:B)W+=abs(z);int D=W+3,s=0;
 vector<I>A(D*D),X(D*D),Y(D*D),G(D*D);A[0]=1;
 for(int z:B){int n=abs(z),sg=unsign||z>0?1:-1;
  fill(X.begin(),X.end(),0);fill(Y.begin(),Y.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){
   I v=A[i*D+j];if(v==0)continue;
   X[abs(i-n)*D+j]+=v;X[(i+n+2)*D+j]-=v;
   Y[i*D+abs(j-n)]+=sg*v;Y[i*D+j+n+2]-=sg*v;
  }s+=n;
  fill(G.begin(),G.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){
   if(i>=2)X[i*D+j]+=X[(i-2)*D+j];
   if(j>=2)Y[i*D+j]+=Y[i*D+j-2];
   G[i*D+j]=X[i*D+j]+Y[i*D+j];
  }A.swap(G);
 }return {D,A};
}
struct Row{vector<int>B;int p,W;I phi;};
I ab(I x){return x<0?-x:x;}
extern "C" int census(const char*input,int threads){
 omp_set_num_threads(threads);istringstream in(input);int K;in>>K;
 vector<Row>R(K);for(auto&r:R){int N;string st;in>>r.W>>r.p>>N>>st;r.phi=I(st);
 r.B.resize(N);for(int&z:r.B)in>>z;}
 atomic<int>passed{0},absblocks{0},errors{0},done{0};vector<int>fails;
 #pragma omp parallel for schedule(dynamic,8)
 for(int q=0;q<K;q++){
  auto B=R[q].B;int N=B.size(),p=R[q].p,i0=-1,j0=-1;pair<int,int>best{-1,-1};
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)if((B[i]-B[j])%2==0){
   auto key=make_pair(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])));
   if(key>best){best=key;i0=i;j0=j;}
  }assert(i0>=0);int a=abs(B[i0]),b=abs(B[j0]),ea=B[i0]>0?1:-1,eb=B[j0]>0?1:-1;
  vector<int>C;for(int i=0;i<N;i++)if(i!=i0&&i!=j0)C.push_back(B[i]);
  auto G=table(C,false),U=table(C,true);I T=0,X=0,Y=0,Z=0,cost=0,gc=G.at(p,0);
  for(int c=abs(a-b);c<=a+b;c+=2){
   for(int s=abs(p-c);s<=p+c;s+=2)T+=G.at(s,0);
   Z+=ea*eb*G.at(p,c);cost+=U.at(p,c);
  }
  for(int s=abs(p-a);s<=p+a;s+=2){X+=eb*G.at(s,b);cost+=U.at(s,b);}
  for(int s=abs(p-b);s<=p+b;s+=2){Y+=ea*G.at(s,a);cost+=U.at(s,a);}
  if(2*(T+X+Y+Z)!=R[q].phi||X+Y>=0||X+Z>=0||Y+Z>=0||gc<0)errors++;
  if(T-gc>=cost)passed++;
  else{
   #pragma omp critical
   fails.push_back(q);
  }
  if(T-gc>=ab(X)+ab(Y)+ab(Z))absblocks++;
  int d=++done;
  if(d%5000==0){
   #pragma omp critical
   {printf("census %d/%d certified=%d errors=%d\n",d,K,passed.load(),errors.load());fflush(stdout);}
  }
 }
 sort(fails.begin(),fails.end(),[&](int a,int b){return make_pair(R[a].W,R[a].B)<make_pair(R[b].W,R[b].B);});
 printf("CENSUS rows=%d unsigned_mixed=%d signed_absolute=%d errors=%d remaining=%zu\n",
 K,passed.load(),absblocks.load(),errors.load(),fails.size());
 if(!fails.empty()){auto&r=R[fails[0]];printf("first remaining W=%d p=%d B=",r.W,r.p);for(int x:r.B)printf("%d,",x);printf("\n");}
 fflush(stdout);return errors;
}

'''
assembly=subprocess.run(["g++","-O2","-pipe","-fPIC","-fopenmp","-S","-x","c++","-","-o","-"],input=source.encode(),stdout=subprocess.PIPE,check=True).stdout
obj=os.memfd_create("m158_obj");lib=os.memfd_create("m158_lib")
subprocess.run(["as","-o",f"/proc/self/fd/{obj}"],input=assembly,pass_fds=(obj,),check=True)
def locate(name):
    return subprocess.check_output(["g++","-print-file-name="+name],text=True).strip()
subprocess.run(["ld","-shared","--eh-frame-hdr","-L"+os.path.dirname(locate("libgcc_s.so")),"-o",f"/proc/self/fd/{lib}",locate("crtbeginS.o"),f"/proc/self/fd/{obj}",locate("libstdc++.so"),locate("libgcc_s.so"),locate("libgomp.so"),"-lc",locate("crtendS.o")],pass_fds=(obj,lib),check=True)
module=ctypes.CDLL(f"/proc/self/fd/{lib}")
assert module.bridges()==0
module.profiles_audit.argtypes=[ctypes.c_int]
assert module.profiles_audit(args.threads)==0
if args.census:
    rows=[]
    for W,p,B,phi in re.findall(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+)",args.census.read_text()):
        B=list(map(int,B.split()))
        assert int(W)==sum(map(abs,B))<=48 and len(B)<=15
        rows.append(f"{W} {abs(int(p))} {len(B)} {phi} "+" ".join(map(str,B)))
    assert len(rows)==33487
    data=(str(len(rows))+"\n"+"\n".join(rows)).encode()
    module.census.argtypes=[ctypes.c_char_p,ctypes.c_int]
    assert module.census(data,args.threads)==0
    profiles={}
    for L in range(6,13):
        for q in (1,2):
            for tau in (0,1):
                for j in range(-1,L):
                    for eps in (-2,2):
                        if j==-1 and eps==2: continue
                        off=tuple(sorted(tau+q*i+(eps if i==j else 0) for i in range(L)))
                        if sum(off)%2: continue
                        low=20 if q==1 else 0
                        while 2*low+off[0]<1: low+=1
                        key=tuple(x-off[0] for x in off)
                        profiles.setdefault(key,set()).add((off[0],low))
    covered={}
    for row in rows:
        fields=list(map(int,row.split()))
        ns=sorted([fields[1]]+list(map(abs,fields[4:])))
        key=tuple(x-ns[0] for x in ns)
        if any((ns[0]-base)%2==0 and (ns[0]-base)//2>=low
               for base,low in profiles.get(key,())):
            covered[len(ns)]=covered.get(len(ns),0)+1
    assert covered=={6:6,7:5}
    print("family census coverage",sum(covered.values()),dict(sorted(covered.items())),flush=True)
print("PASS",flush=True)
