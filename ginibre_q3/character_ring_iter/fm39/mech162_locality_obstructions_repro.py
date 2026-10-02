import argparse, ctypes, os, re, subprocess
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
import sympy as sp

ap=argparse.ArgumentParser(description="FM-MECH162: exact tests; no files written.")
ap.add_argument("--threads",type=int,default=24)
ap.add_argument("--skip-census",action="store_true")
ap.add_argument("--census",type=Path,default=Path(
 "/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/"
 "scratchpad/flipdesc/fx3b_48.log"))
args=ap.parse_args()
assert 1<=args.threads<=32

def table(B):
    out={(0,0):1}
    for z in B:
        n=abs(z);sg=1 if z>0 else -1;new={}
        for (i,j),v in out.items():
            for k in range(abs(i-n),i+n+1,2):
                new[k,j]=new.get((k,j),0)+v
            for k in range(abs(j-n),j+n+1,2):
                new[i,k]=new.get((i,k),0)+sg*v
        out={k:v for k,v in new.items() if v}
    return out

def top(B):
    return max(((i,j) for i in range(len(B)) for j in range(i+1,len(B))
        if (B[i]-B[j])%2==0),
        key=lambda ij:(abs(B[ij[0]])+abs(B[ij[1]]),
                       max(abs(B[ij[0]]),abs(B[ij[1]]))))
def blocks(B,p):
    i,j=top(B);a,b=abs(B[i]),abs(B[j])
    ea,eb=(1 if B[i]>0 else -1),(1 if B[j]>0 else -1)
    C=[z for k,z in enumerate(B) if k not in (i,j)]
    G=table(C);T=X=Y=Z=0
    for c in range(abs(a-b),a+b+1,2):
        T+=sum(G.get((s,0),0) for s in range(abs(p-c),p+c+1,2))
        Z+=ea*eb*G.get((p,c),0)
    X=eb*sum(G.get((s,b),0) for s in range(abs(p-a),p+a+1,2))
    Y=ea*sum(G.get((s,a),0) for s in range(abs(p-b),p+b+1,2))
    return T,X,Y,Z,G.get((p,0),0)

# An actual failure of the three-pair implication.
B=[-1]*13+[-2]+[-3]*5+[-4]
T,X,Y,Z,h=blocks(B,6)
assert (T,X,Y,Z,h)==(8595054780,-5602622750,-7494310590,
                     5202486360,723462700)
assert max(X+Y,X+Z,Y+Z)<0 and T+X+Y+Z-h==-22854900
print("actual triple-only obstruction: drop=-22854900")

@lru_cache(None)
def invariant(ns):
    v={0:1}
    for n in ns:
        w={}
        for i,c in v.items():
            for j in range(abs(i-n),i+n+1,2):w[j]=w.get(j,0)+c
        v=w
    return v.get(0,0)

# Full subset expansion, and the origin-coordinate obstruction.
B=(-1,-3,-3,-4,-5,-5,-6);p=7
ns=tuple(map(abs,B))+(p,);L=len(ns);F=(1<<L)-1
ms=[invariant(tuple(sorted(ns[i] for i in range(L) if s>>i&1)))
    for s in range(1<<L)]
w=[]
for s in range(1<<(L-1)):
    neg=sum(B[i]<0 for i in range(L-1) if s>>i&1)
    w.append((-1)**neg*ms[s]*ms[F^s])
T,X,Y,Z,h=blocks(B,p)
assert (T,X,Y,Z,h)==(741,-98,-164,71,31)
assert [sum(v for s,v in enumerate(w)
    if (2*((s>>3)&1)+((s>>6)&1))==q) for q in range(4)]==[T,X,Y,Z]
cuts=[sum(v for s,v in enumerate(w)
          if ((s>>i)&1)!=((s>>j)&1)) for i,j in combinations(range(L),2)]
assert max(cuts)==-3 and sum(w)==550 and ms[F]==624
assert max(abs(v) for v in w[1:])==46
def phi(selected,neg):
    out=0;s=selected
    while True:
        out+=(-1 if (s&neg).bit_count()%2 else 1)*ms[s]*ms[selected^s]
        if not s:return out
        s=(s-1)&selected
proper=0
for selected in range(F):
    neg=selected
    while True:
        if neg.bit_count()%2==0:
            assert phi(selected,neg)>=0;proper+=1
        if not neg:break
        neg=(neg-1)&selected
assert proper==3153
assert min(phi(F,q) for q in range(1<<L) if q.bit_count()%2==0)==1060
# Only m(full) changes: all cuts and proper-list values remain fixed.
ms[F]=104
assert min(phi(F,q) for q in range(1<<L) if q.bit_count()%2==0)==20
assert phi(F,F)//2==30<h
assert T-520-h-abs(X)-abs(Y)-abs(Z)==-143
print("origin test: 28 negative cuts; 3153 proper sign checks; "
      "full minimum=20; child=31 > parent=30")

# Exact joint moments; direct Catalan expansion supplies an independent bridge.
@lru_cache(None)
def moment(m,j):
    return Q(2*factorial(2*m)*factorial(2*m+1)*factorial(2*j)*factorial(2*j+1),
      factorial(m)**2*factorial(j)**2*factorial(m+j+1)*factorial(m+j+2))
def catmoment(n):
    return 0 if n%2 else comb(n,n//2)//(n//2+1)
for m in range(7):
    for j in range(7):
        direct=sum(comb(2*m,i)*comb(2*j,k)*(-1)**k*
          catmoment(2*m+2*j-i-k)*catmoment(i+k)
          for i in range(2*m+1) for k in range(2*j+1))
        assert moment(m,j)==direct

def multiply(P,R):
    out={}
    for (i,j),v in P.items():
        for (k,l),u in R.items():
            out[i+k,j+l]=out.get((i+k,j+l),Q(0))+v*u
    return {key:v for key,v in out.items() if v}
def factor(n,eps):
    out={}
    for h0 in range(n//2+1):
        d=n-2*h0;c=Q((-1)**h0*comb(n-h0,h0),2**d)
        for j in range(d+1):
            v=c*comb(d,j)*(1+eps*(-1)**j)
            if v:out[d-j,j]=out.get((d-j,j),Q(0))+v
    return out
def word(ns):
    out={(0,0):Q(1)}
    for z in ns:out=multiply(out,factor(abs(z),1 if z>0 else -1))
    assert all(i%2==j%2==0 for i,j in out)
    return {(i//2,j//2):v for (i,j),v in out.items()}
def val(P,A):
    return sum(v*moment(A+i,j) for (i,j),v in P.items())

cores=[-2,-4,-6,-8,10]
profiles=[cores,[-2,-4,10]]
for i,j in combinations(range(5),2):
    ns0=cores.copy();ns0[i]*=-1;ns0[j]*=-1;profiles.append(ns0)
for ns0 in profiles:
    P=word(ns0)
    for A in range(4):
        assert val(P,A)==table([1]*(2*A)+ns0).get((0,0),0)
print("moment bridges=49; word bridges=48")

A=sp.Symbol("A")
@lru_cache(None)
def basis(D,i,j):
    return sp.Poly(16**i*sp.rf(A+sp.Rational(1,2),i)*
      sp.rf(A+sp.Rational(3,2),i)*
      sp.Rational(factorial(2*j)*factorial(2*j+1),factorial(j)**2)*
      sp.rf(A+2+i+j,D-i-j)*sp.rf(A+3+i+j,D-i-j),A)
def numerator(P,D):
    out=sp.Poly(0,A)
    for (i,j),v in P.items():
        out+=sp.Rational(v.numerator,v.denominator)*basis(D,i,j)
    return out
P=word(cores);D=15;base=numerator(P,D)
certs=[(base,1024),(numerator(word(profiles[1]),D)-base,1024)]
certs += [(numerator(word(ns0),D)-base,1024) for ns0 in profiles[2:]]
# The +1,+1 flip: suffix index A-1, polynomial (s^2-d^2)P.
ones={}
for (i,j),v in P.items():
    ones[i+1,j]=ones.get((i+1,j),Q(0))+v
    ones[i,j+1]=ones.get((i,j+1),Q(0))-v
certs.append((numerator(ones,D+1),1023))
coeffs=0
for poly,shift in certs:
    cc=poly.shift(shift).all_coeffs()
    assert min(cc)>0
    coeffs+=len(cc)
assert val(P,1024)<val(word(profiles[1]),1024)
print("uniform A>=1024 certificates",len(certs),"positive coefficients",coeffs)
if args.skip_census:
    print("PASS (census skipped)")
    raise SystemExit

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
I ab(I x){return x<0?-x:x;}
struct Table{int D;vector<I>v;I at(int i,int j)const{return i>=0&&j>=0&&i<D&&j<D?v[i*D+j]:I(0);}};
Table table(const vector<int>&B){
 int W=0;for(int z:B)W+=abs(z);int D=W+3,s=0;vector<I>A(D*D),X(D*D),Y(D*D),G(D*D);A[0]=1;
 for(int z:B){int n=abs(z),sg=z>0?1:-1;fill(X.begin(),X.end(),0);fill(Y.begin(),Y.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){I v=A[i*D+j];if(v==0)continue;
   X[abs(i-n)*D+j]+=v;X[(i+n+2)*D+j]-=v;Y[i*D+abs(j-n)]+=sg*v;Y[i*D+j+n+2]-=sg*v;
  }s+=n;fill(G.begin(),G.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=s-i;j++){
   if(i>=2)X[i*D+j]+=X[(i-2)*D+j];if(j>=2)Y[i*D+j]+=Y[i*D+j-2];G[i*D+j]=X[i*D+j]+Y[i*D+j];
  }A.swap(G);
 }return {D,A};
}
struct Blocks{I T,X,Y,Z,h,g,S;int ia,ib;vector<int>C;};
Blocks blocks(const vector<int>&B,int p){
 int N=B.size(),ia=-1,ib=-1;pair<int,int>best{-1,-1};
 for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)if((B[i]-B[j])%2==0){
  auto key=make_pair(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])));if(key>best){best=key;ia=i;ib=j;}}
 assert(ia>=0);vector<int>C;for(int i=0;i<N;i++)if(i!=ia&&i!=ib)C.push_back(B[i]);
 int a=abs(B[ia]),b=abs(B[ib]),ea=B[ia]>0?1:-1,eb=B[ib]>0?1:-1;
 auto A=table(C);I T=0,X=0,Y=0,Z=0,h=A.at(p,0);
 for(int c=abs(a-b);c<=a+b;c+=2){for(int s=abs(p-c);s<=p+c;s+=2)T+=A.at(s,0);Z+=ea*eb*A.at(p,c);}
 for(int s=abs(p-a);s<=p+a;s+=2)X+=eb*A.at(s,b);
 for(int s=abs(p-b);s<=p+b;s+=2)Y+=ea*A.at(s,a);
 return {T,X,Y,Z,h,T+X+Y+Z,T-h-ab(X)-ab(Y)-ab(Z),ia,ib,C};
}
I flip(const vector<int>&L,int i,int j){
 vector<int>C;for(int k=0;k<(int)L.size();k++)if(k!=i&&k!=j)C.push_back(L[k]);
 auto A=table(C);return (L[j]>0?1:-1)*A.at(abs(L[i]),abs(L[j]));
}

struct Row{int W,p;vector<int>B;I phi;Blocks V;};
extern "C" int audit(const char*data,int threads){
 omp_set_num_threads(threads);istringstream in(data);int K;in>>K;
 assert(K==33487);vector<Row>rs(K);atomic<int>done{0};
 for(auto&r:rs){int n;string ph;in>>r.W>>r.p>>n>>ph;
  assert(r.W<=48&&n<=15&&r.p<=48);r.phi=I(ph);
  r.B.resize(n);for(int&z:r.B)in>>z;}
 #pragma omp parallel for schedule(dynamic,8)
 for(int k=0;k<K;k++){auto&r=rs[k];r.V=blocks(r.B,r.p);auto&v=r.V;
  assert(2*v.g==r.phi&&v.S>=0&&v.X+v.Y<0&&v.X+v.Z<0&&v.Y+v.Z<0);
  int d=++done;if(d%5000==0){
   #pragma omp critical
   {printf("census %d/%d\n",d,K);fflush(stdout);}
  }
 }
 sort(rs.begin(),rs.end(),[](const Row&a,const Row&b){
  I x=a.V.S*b.V.g,y=b.V.S*a.V.g;return x==y?a.W<b.W:x<y;});
 assert(rs[0].V.S==377&&rs[0].V.g==550);
 vector<vector<int>>seen;int printed=0;
 for(auto&r:rs){
  vector<int>key;for(int z:r.B)key.push_back(abs(z));key.push_back(r.p);
  if(find(seen.begin(),seen.end(),key)!=seen.end())continue;
  seen.push_back(key);auto&v=r.V;vector<int>L=r.B;int sig=1;
  for(int z:r.B)if(z<0)sig=-sig;L.push_back(sig*r.p);
  bool first=true;I near;pair<int,int>at;
  for(int i=0;i<(int)L.size();i++)for(int j=i+1;j<(int)L.size();j++){
   I z=flip(L,i,j);assert(z<0);
   if(first||z>near){first=false;near=z;at={L[i],L[j]};}
  }
  printf("EXTREME W=%d p=%d S/g=%s/%s nearest=(%d,%d):%s\n",
    r.W,sig*r.p,v.S.str().c_str(),v.g.str().c_str(),
    at.first,at.second,near.str().c_str());fflush(stdout);
  if(++printed==3)break;
 }
 puts("CENSUS 33487/33487 signed bounds PASS");fflush(stdout);return 0;
}

'''

assembly=subprocess.run(["g++","-O2","-pipe","-fPIC","-fopenmp","-S",
 "-x","c++","-","-o","-"],input=source.encode(),stdout=subprocess.PIPE,check=True).stdout
obj=os.memfd_create("fm162_obj");lib=os.memfd_create("fm162_lib")
subprocess.run(["as","-o",f"/proc/self/fd/{obj}"],input=assembly,pass_fds=(obj,),check=True)
def locate(n):
    return subprocess.check_output(["g++","-print-file-name="+n],text=True).strip()
subprocess.run(["ld","-shared","--eh-frame-hdr","-L"+os.path.dirname(locate("libgcc_s.so")),
 "-o",f"/proc/self/fd/{lib}",locate("crtbeginS.o"),f"/proc/self/fd/{obj}",
 locate("libstdc++.so"),locate("libgcc_s.so"),locate("libgomp.so"),
 "-lc",locate("crtendS.o")],pass_fds=(obj,lib),check=True)
module=ctypes.CDLL(f"/proc/self/fd/{lib}")
module.audit.argtypes=[ctypes.c_char_p,ctypes.c_int]
rows=[]
for W,p,B,phi in re.findall(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)\s+phi=(\d+)",
                            args.census.read_text()):
    B=list(map(int,B.split()))
    assert int(W)==sum(map(abs,B))<=48 and len(B)<=15
    rows.append(f"{W} {abs(int(p))} {len(B)} {phi} "+" ".join(map(str,B)))
assert len(rows)==33487
assert module.audit((str(len(rows))+"\n"+"\n".join(rows)).encode(),args.threads)==0
print("PASS")
