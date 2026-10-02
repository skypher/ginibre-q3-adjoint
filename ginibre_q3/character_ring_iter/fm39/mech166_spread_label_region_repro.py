import argparse, ctypes, os, re, subprocess
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from collections import Counter
from itertools import combinations_with_replacement
from random import Random

ap=argparse.ArgumentParser(description="FM-MECH166 exact verifier; memory only.")
ap.add_argument("--threads",type=int,default=24)
ap.add_argument("--skip-census",action="store_true")
ap.add_argument("--census",type=Path,default=Path(
 "/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/"
 "scratchpad/flipdesc/fx3b_48.log"))
args=ap.parse_args()
assert 1<=args.threads<=32

def cg(a,b):return range(abs(a-b),a+b+1,2)
def fusion(ns):
    v={0:1}
    for n in ns:
        w={}
        for s,c in v.items():
            for t in cg(s,n):w[t]=w.get(t,0)+c
        v=w
    return v
checks=0
for ns in combinations_with_replacement(range(1,19),3):
    m=min(ns);v=fusion(ns)
    for r in range(m+1):
        c2=min(2*r+2,m+2);w=fusion(ns+(2*r,))
        for t in range(sum(ns)+1):
            assert 2*w.get(t,0)>=c2*v.get(t,0);checks+=1
assert checks==240711
def triple(a,b,c,t):
    if (a+b+c+t)%2:return 0
    lo=max(abs(a-b),abs(c-t));hi=min(a+b,c+t)
    return max(0,(hi-lo)//2+1)
rng=Random(166)
for _ in range(2000):
    a=rng.randrange(1,41);b=rng.randrange(a,1000001);c=rng.randrange(a,1000001)
    u=abs(a-b)+2*rng.randrange(a+1)
    t=abs(u-c)+2*rng.randrange(min(u,c)+1)
    r=rng.randrange(a+1)
    assert 2*sum(triple(a,b,c,s) for s in cg(t,2*r)) >= \
        min(2*r+2,a+2)*triple(a,b,c,t)
print("profile lemma: 240711 small + 2000 large-label checks PASS")

def capacities(m):
    q=m//2
    if m%2:
        return Q((q+1)*(3*q+5),2),Q((q+1)*(2*q+3)*(8*q+13),12)
    return Q((q+1)*(3*q+2),2),Q((q+1)*(8*q*q+13*q+6),6)
for m in range(1,101):
    c=[Q(min(2*r+2,m+2),2) for r in range(m+1)]
    assert capacities(m)==(sum(c),sum(x*x for x in c))
    aa,bb=capacities(m)
    assert aa>=Q(3*(m+1)**2,8) and bb>=Q((m+1)**3,6)
def nu(L):return 2**(L-2)-(L if L%2==0 else 2*L-2)
def rho(L,m,b=None):
    if b is None:b=m
    aa,bb=capacities(m);J=L-2
    if L==6:return 1-Q(nu(L),1)/bb-(1+Q(3,m+1))/(b+1)
    return 1-Q(nu(L),1)/bb-(1+Q(comb(J,2),1)/aa+
        Q(2**(J-1)-1-J-comb(J,2),1)/bb)/(b+1)
threshold={}
for L in range(6,21):
    m=1
    while rho(L,m)<0:m+=1
    f=1
    while capacities(f)[1]<nu(L):f+=1
    threshold[L]=m
    print("threshold",L,"FM3",f,"removal",m)
def top(B):
    return max(((i,j) for i in range(len(B)) for j in range(i+1,len(B))
                if (B[i]-B[j])%2==0),
               key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                              max(abs(B[ij[0]]),abs(B[ij[1]]))))
def refined(B,p):
    i,j=top(B);L=len(B)+1;m=min(map(abs,B));b=min(abs(B[i]),abs(B[j]))
    full=B+[(-1)**sum(z<0 for z in B)*p]
    child=[z for k,z in enumerate(full) if k not in (i,j)]
    em=sum(z<0 and z%2==0 for z in full)
    ep=sum(z>0 and z%2==0 for z in full)
    om=sum(z<0 and z%2!=0 for z in full)
    op=sum(z>0 and z%2!=0 for z in full)
    if em+om==0 or (em==0 and op==0):neg=0
    else:neg=2**(L-2 if om+op==0 else L-3)-em-em*ep-om*op
    aa,bb=capacities(m)
    if L==6:
        cc=sorted(Counter(map(abs,child)).values())
        pp=3 if cc==[4] else 1 if cc==[2,2] else 0
        return 1-Q(neg,1)/bb-(1+Q(pp,m+1))/(b+1)
    J=len(child);o=sum(abs(z)%2 for z in child);e=J-o
    bulk=2**(J-1 if o==0 else J-2)-1-e-comb(e,2)-comb(o,2)
    equal=sum(comb(n,2) for n in Counter(map(abs,child)).values())
    assert min(neg,bulk,equal)>=0
    return 1-Q(neg,1)/bb-(1+Q(equal,1)/aa+Q(bulk,1)/bb)/(b+1)

# Label profiles: all signs and all same-parity background removals.
profiles=[]
for L in range(6,14):
    m=threshold[L]
    profiles += [[m]*L,[m+2*i for i in range(L)],
                 [m]+sorted(rng.randrange(m,12*m+1) for _ in range(L-1)),
                 [m,m+1]+[100*m+i for i in range(L-2)]]
for ns in profiles:
    ns.sort()
    if sum(ns)%2:ns[-1]+=1

source=r'''
#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cassert>
#include <omp.h>
using I=boost::multiprecision::number<boost::multiprecision::cpp_int_backend<>,boost::multiprecision::et_off>;
using namespace std;
using V=vector<I>;
V mults(const vector<int>&ns){
 int L=ns.size(),Z=1<<L,F=Z-1;vector<int>w(Z),sz(Z),mx(Z);
 for(int s=1;s<Z;s++){
  int bit=s&-s,i=__builtin_ctz((unsigned)bit);
  w[s]=w[s^bit]+ns[i];sz[s]=sz[s^bit]+1;mx[s]=max(mx[s^bit],ns[i]);
 }
 int cap=w[F]/2;vector<V>C(L-1,V(cap+1,I(1)));
 for(int d=1;d<L-1;d++)for(int q=1;q<=cap;q++)C[d][q]=C[d][q-1]+C[d-1][q];
 V out(Z);out[0]=1;
 for(int s=1;s<Z;s++){
  int l=sz[s];if(l<2||w[s]%2||2*mx[s]>w[s])continue;
  I v=0;for(int j=s;;j=(j-1)&s){
   int q=w[s]/2-w[j]-sz[j];
   if(q>=0)v+=(sz[j]%2?-C[l-2][q]:C[l-2][q]);if(!j)break;
  }assert(v>=0);out[s]=v;
 }return out;
}
V walsh(V v){
 for(int h=1;h<(int)v.size();h*=2)
  for(int i=0;i<(int)v.size();i+=2*h)for(int j=i;j<i+h;j++){
   I a=v[j],b=v[j+h];v[j]=a+b;v[j+h]=a-b;
  }return v;
}
V phi(const V&m,const vector<int>&idx){
 int Z=1<<idx.size(),F=0;for(int i:idx)F|=1<<i;V f(Z);
 for(int s=0;s<Z;s++){
  int mask=0;for(int j=0;j<(int)idx.size();j++)if(s>>j&1)mask|=1<<idx[j];
  f[s]=m[mask]*m[F^mask];
 }return walsh(f);
}
I inv_direct(const vector<int>&ns){
 int W=0;for(int n:ns)W+=n;V v(W+3),w(W+3);v[0]=1;int s=0;
 for(int n:ns){
  fill(w.begin(),w.end(),I(0));
  for(int t=0;t<=s;t++)if(v[t]!=0){
   w[abs(t-n)]+=v[t];w[t+n+2]-=v[t];
  }s+=n;
  for(int t=2;t<=s;t++)w[t]+=w[t-2];v.swap(w);
 }return v[0];
}
extern "C" int words(const char*data,int threads){
 omp_set_num_threads(threads);istringstream in(data);int K;in>>K;
 vector<vector<int>>profiles(K);
 for(auto&ns:profiles){int L;in>>L;assert(6<=L&&L<=13);ns.resize(L);
  for(int&n:ns){in>>n;assert(1<=n&&n<=10000);}}
 atomic<long long>tests{0},bridges{0};
 #pragma omp parallel for schedule(dynamic)
 for(int q=0;q<K;q++){
  auto ns=profiles[q];int L=ns.size(),m=ns[0],F=(1<<L)-1,J=L-2;
  auto ms=mults(ns);vector<int>idx(L);for(int i=0;i<L;i++)idx[i]=i;
  auto parent=phi(ms,idx);I A2=0,B4=0;
  for(int r=0;r<=m;r++){I c=min(2*r+2,m+2);A2+=c;B4+=c*c;}
  I neg=(I(1)<<(L-2))-(L%2?2*L-2:L);
  for(int z=0;z<16;z++){
   int mask=(z*137+q*13)&F;vector<int>sub;
   for(int i=0;i<L;i++)if(mask>>i&1)sub.push_back(ns[i]);
   assert(ms[mask]==inv_direct(sub));bridges++;
  }
  for(int a=0;a<L-1;a++)for(int b=a+1;b<L-1;b++)if((ns[a]-ns[b])%2==0){
   vector<int>ci;for(int i=0;i<L;i++)if(i!=a&&i!=b)ci.push_back(i);
   auto child=phi(ms,ci);int v=min(ns[a],ns[b]);
   int childmask=F^(1<<a)^(1<<b);
   assert(ms[F]>=I(v+1)*ms[childmask]);
   I num,den;
   if(L==6){
    den=I(v+1)*(m+1)*B4;
    num=I(v+1)*(m+1)*(B4-4*neg)-I(m+4)*B4;
   }else{
    I pairs=I(J)*(J-1)/2,bulk=(I(1)<<(J-1))-1-J-pairs;
    den=I(v+1)*A2*B4;
    num=I(v+1)*A2*(B4-4*neg)-A2*B4-2*pairs*B4-4*bulk*A2;
   }assert(num>=0);
   for(int mask=0;mask<=F;mask++){
    if(__builtin_popcount((unsigned)mask)%2)continue;
    bool ok=true;for(int i=1;i<L;i++)if(ns[i]==ns[i-1]&&
       ((mask>>i&1)!=(mask>>(i-1)&1)))ok=false;
    if(!ok)continue;
    int cm=0;for(int j=0;j<J-1;j++)cm|=((mask>>ci[j])&1)<<j;
    cm|=(__builtin_popcount((unsigned)cm)%2)<<(J-1);
    assert(B4*parent[mask]>=2*(B4-4*neg)*ms[F]);
    assert(den*(parent[mask]-child[cm])>=2*num*ms[F]);tests++;
   }
  }
 }
 printf("caller checks=%lld independent fusion bridges=%lld profiles=%d PASS\n",
        tests.load(),bridges.load(),K);fflush(stdout);return 0;
}
struct Table{
 int D;V v;I at(int i,int j)const{return i>=0&&j>=0&&i<D&&j<D?v[i*D+j]:I(0);}
};
Table table(const vector<int>&B){
 int W=0;for(int z:B)W+=abs(z);int D=W+3,s=0;
 V A(D*D),X(D*D),Y(D*D),G(D*D);A[0]=1;
 for(int z:B){
  int n=abs(z),sg=z>0?1:-1;
  fill(X.begin(),X.end(),I(0));fill(Y.begin(),Y.end(),I(0));
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++)if(A[i*D+j]!=0){
   I v=A[i*D+j];X[abs(i-n)*D+j]+=v;X[(i+n+2)*D+j]-=v;
   Y[i*D+abs(j-n)]+=sg*v;Y[i*D+j+n+2]-=sg*v;
  }s+=n;fill(G.begin(),G.end(),I(0));
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){
   if(i>=2)X[i*D+j]+=X[(i-2)*D+j];
   if(j>=2)Y[i*D+j]+=Y[i*D+j-2];
   G[i*D+j]=X[i*D+j]+Y[i*D+j];
  }A.swap(G);
 }return {D,A};
}
extern "C" int census(const char*data,int threads){
 omp_set_num_threads(threads);istringstream in(data);int K;in>>K;
 struct Row{int p,a,b;vector<int>B;I value,num,den;};vector<Row>rows(K);
 for(auto&r:rows){
  int n;string f,u,v;in>>r.p>>n>>r.a>>r.b>>f>>u>>v;
  assert(5<=n&&n<=15&&r.p<=48);r.value=I(f);r.num=I(u);r.den=I(v);
  r.B.resize(n);for(int&z:r.B)in>>z;
 }
 atomic<int>done{0};
 #pragma omp parallel for schedule(dynamic,8)
 for(int t=0;t<K;t++){
  auto&r=rows[t];vector<int>C,ns;
  int W=0;for(int k=0;k<(int)r.B.size();k++){
   int z=r.B[k];W+=abs(z);ns.push_back(abs(z));
   if(k!=r.a&&k!=r.b)C.push_back(z);
  }assert(W<=48);ns.push_back(r.p);
  I parent=table(r.B).at(r.p,0),child=table(C).at(r.p,0),dim=inv_direct(ns);
  assert(2*parent==r.value);
  assert(r.den*(parent-child)>=r.num*dim&&r.num>=0);
  int d=++done;if(d%1000==0){
   #pragma omp critical
   {printf("census callers %d/%d\n",d,K);fflush(stdout);}
  }
 }
 printf("census quantitative callers=%d PASS\n",K);fflush(stdout);return 0;
}
'''
assembly=subprocess.run(["g++","-O2","-pipe","-fPIC","-fopenmp","-S",
 "-x","c++","-","-o","-"],input=source.encode(),stdout=subprocess.PIPE,check=True).stdout
obj=os.memfd_create("fm166_obj");lib=os.memfd_create("fm166_lib")
subprocess.run(["as","-o",f"/proc/self/fd/{obj}"],input=assembly,pass_fds=(obj,),check=True)
def locate(n):
    return subprocess.check_output(["g++","-print-file-name="+n],text=True).strip()
subprocess.run(["ld","-shared","--eh-frame-hdr","-L"+os.path.dirname(locate("libgcc_s.so")),
 "-o",f"/proc/self/fd/{lib}",locate("crtbeginS.o"),f"/proc/self/fd/{obj}",
 locate("libstdc++.so"),locate("libgcc_s.so"),locate("libgomp.so"),
 "-lc",locate("crtendS.o")],pass_fds=(obj,lib),check=True)
module=ctypes.CDLL(f"/proc/self/fd/{lib}")
module.words.argtypes=[ctypes.c_char_p,ctypes.c_int]
data=str(len(profiles))+"\n"+"\n".join(
 str(len(ns))+" "+" ".join(map(str,ns)) for ns in profiles)
assert module.words(data.encode(),args.threads)==0

if not args.skip_census:
    counts=[Counter(),Counter(),Counter()];total=0;rows=[]
    for W,p,B,phi in re.findall(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+)",
                                args.census.read_text()):
        B=list(map(int,B.split()));p=abs(int(p));L=len(B)+1
        assert sum(map(abs,B))==int(W)<=48 and L>=6
        total+=1;m=min(map(abs,B));i,j=top(B);b=min(abs(B[i]),abs(B[j]))
        vals=[rho(L,m),rho(L,m,b),refined(B,p)]
        for q,v in enumerate(vals):
            if v>=0:counts[q][L]+=1
        v=vals[-1]
        if v>=0:
            rows.append(f"{p} {len(B)} {i} {j} {phi} {v.numerator} {v.denominator} "+
                        " ".join(map(str,B)))
    assert total==33487
    assert [sum(c.values()) for c in counts]==[1607,2161,4391]
    for q,c in enumerate(counts):
        print("census criterion",q,"covered",sum(c.values()),"by L",dict(sorted(c.items())))
    module.census.argtypes=[ctypes.c_char_p,ctypes.c_int]
    assert module.census((str(len(rows))+"\n"+"\n".join(rows)).encode(),args.threads)==0
print("FM-MECH166 PASS")
