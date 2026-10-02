import argparse,ctypes,os,pathlib,re,subprocess
from collections import defaultdict
from fractions import Fraction as Q
from math import comb
from random import Random

ap=argparse.ArgumentParser(description="Read-only exact FM-MECH156 verifier.")
ap.add_argument("--no-logs",action="store_true")
ap.add_argument("--random",type=int,default=0,help="Seeded random attempts, at most 20000.")
ap.add_argument("--threads",type=int,default=12)
ap.add_argument("--logs",default="/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc")
args=ap.parse_args()
assert 0<=args.random<=20000 and 1<=args.threads<=32

def kernel(B,deform=False):
 M={(0,0,0):1}
 for z in B:
  n=abs(z);eps=1 if z>0 else -1;F=defaultdict(int)
  for j in range(n+1):F[j,j,0]+=1
  if deform:
   for r in range(n//2+1):
    k=n-2*r
    for j in range(k+1):
     F[j+r,k-j+r,k]+=eps*(-1)**r*comb(n-r,r)*comb(k,j)
  else:
   for j in range(n+1):F[j,n-j,0]+=eps
  O=defaultdict(int)
  for (i,j,k),v in M.items():
   for (s,t,l),c in F.items():O[i+s,j+t,k+l]+=v*c
  M={key:v for key,v in O.items() if v}
 return M

rng=Random(156)
profiles=[(),(1,),(-1,),(-2,),(-3,2),(-3,-3,1,1),(-1,2,3,-4)]
profiles += [tuple(rng.choice((-1,1))*rng.randrange(1,5)
 for i in range(rng.randrange(1,5))) for j in range(60)]
checks=0
for B in profiles:
 T=sum(map(abs,B));M=kernel(B);V=kernel(B,True)
 def H(i,j):
  i,j=max(i,j),min(i,j)
  return M.get((i,j,0),0)-M.get((i+1,j-1,0),0)
 def moment(k,i,j):
  i,j=max(i,j),min(i,j);t=i-j
  z=Q(k+t,2);f=Q(k-t,2);c=1/(z+1);v=0
  for s in range(j+1):
   v+=c*H(i+s,j-s);c*=Q(f-s,z+s+2)
  return v
 for k in range(5):
  direct=defaultdict(Q)
  for (i,j,r),v in V.items():direct[i,j]+=Q(2*v,k+r+2)
  for i in range(T+1):
   for j in range(T+1):
    assert direct[i,j]==moment(k,i,j);checks+=1
 K=[[Q(H(i,j),abs(i-j)+1) for j in range(T+1)]
    for i in range(T+1)]
 for k in range(T+1):
  assert K[k][k]>0
  for i in range(k+1,T+1):
   for j in range(i,T+1):
    K[i][j]=K[j][i]=K[i][j]-K[i][k]*K[k][j]/K[k][k]
print("exact moment identities",checks,"Gram profiles",len(profiles),"PASS")

rows=[]
if not args.no_logs:
 for name in ("fx5_40.log","tpg_one.log","tpg_all.log"):
  for line in (pathlib.Path(args.logs)/name).read_text().splitlines():
   if not line.startswith(("NOFLIP ","TPFAIL ")):continue
   B=tuple(map(int,re.search(r"\bB=(.*?)(?:  |$)",line)[1].split()))
   p=abs(int(re.search(r"\bp=(-?\d+)",line)[1]))
   tag=(0 if line.startswith("NOFLIP ") else 1) if name=="fx5_40.log" else 2
   rows.append((tag,p,B))
 assert [sum(t==j for t,p,B in rows) for j in range(3)]==[5430,126,64]
rows += [(3,8,(1,)*24+(-3,-3)),
         (3,6,(1,)*6+(-2,)+(3,)*6+(-4,)),
         (3,26,(1,)*14+(-25,-25,-12))]
rng=Random(1562026)
for it in range(args.random):
 sg={n:rng.choice((-1,1)) for n in range(1,41)}
 if it%2:
  ns=[rng.randrange(1,41) for j in range(rng.randrange(5,16))]
 else:
  ns=[n for n in range(1,4) for j in range(rng.randrange(0,32))]
  ns += [rng.randrange(4,25) for j in range(rng.randrange(1,7))]
 W=sum(ns)
 if not ns or W>220:continue
 m=max(ns);lo=max(8,m);hi=(W-max(6,m))//2
 if lo>hi:continue
 d=rng.randrange(lo,hi+1);p=W-2*d
 B=tuple(sg[n]*n for n in ns)
 sig=(-1)**sum(z<0 for z in B)
 if -sig*p in B:continue
 rows.append((4,p,B))
for tag,p,B in rows:
 W=sum(map(abs,B));d=(W-p)//2;sig=(-1)**sum(z<0 for z in B)
 assert all(0<abs(z)<=220 for z in B) and len(B)<=220 and W<=220
 assert (W-p)%2==0 and p>=max(6,max(map(abs,B))) and d>=8
 assert max(map(abs,B))<=d and -sig*p not in B
 assert all(-z not in B for z in B)

CPP=r'''
#include <bits/stdc++.h>
#include <gmpxx.h>
#include <omp.h>
using namespace std;
using I=mpz_class;
using Q=mpq_class;
int weight(const vector<int>&B){int w=0;for(int z:B)w+=abs(z);return w;}
vector<I> table(const vector<int>&B,bool torus){
 int w=weight(B),D=w+1,s=0;vector<I>A(D*D),O(D*D);A[0]=1;
 for(int z:B){
  int n=abs(z),e=z>0?1:-1;fill(O.begin(),O.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=(torus?s:s-i);j++){
   const I&v=A[i*D+j];if(v==0)continue;
   if(torus)for(int r=0;r<=n;r++){
    O[(i+r)*D+j+r]+=v;O[(i+r)*D+j+n-r]+=e*v;
   }else{
    for(int r=abs(i-n);r<=i+n;r+=2)O[r*D+j]+=v;
    for(int r=abs(j-n);r<=j+n;r+=2)O[i*D+r]+=e*v;
   }
  }A.swap(O);s+=n;
 }return A;
}
struct Record{int tag,p;vector<int>B;};
extern "C" int audit(const char*text,int threads){
 omp_set_num_threads(threads);istringstream in(text);int n;in>>n;
 if(!in||n<0||n>50000)return 1;
 vector<Record>rs(n);
 for(auto&r:rs){
  int k;in>>r.tag>>r.p>>k;
  if(!in||k<2||k>220||r.p<0||r.p>220)return 2;
  r.B.resize(k);
  for(int&z:r.B){in>>z;if(!in||z==0||z<-220||z>220)return 3;}
 }
 atomic<int>bad{0},done{0},nf{0},cf{0},gf{0},target{0},cert{0};
 atomic<int>rnd{0},rhigh{0},rfail{0};
 #pragma omp parallel for schedule(dynamic,1)
 for(int rr=0;rr<n;rr++){
  auto [tag,p,B]=rs[rr];int W=weight(B),N=B.size();
  if(W>220||p>W||(W-p)%2){++bad;continue;}
  int ia=-1,ib=-1,w=-1,ma=-1;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)
   if((B[i]-B[j])%2==0){
    auto key=make_pair(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])));
    if(key>make_pair(w,ma)){ia=i;ib=j;w=key.first;ma=key.second;}
   }
  if(ia<0){++bad;continue;}
  int za=B[ia],zb=B[ib];if(abs(za)<abs(zb))swap(za,zb);
  int a=abs(za),b=abs(zb),ea=za>0?1:-1,eb=zb>0?1:-1;
  auto C=B;C.erase(C.begin()+ib);C.erase(C.begin()+ia);
  int T=weight(C),D=T+1,d=(W-p)/2,A=d-a,B0=d-b,L=d-a-b-1;
  int h=d-(a+b)/2;
  if(p<a||A<0){++bad;continue;}
  auto F=table(C,false),FB=table(B,false),M=table(C,true);
  auto G=[&](int i,int j)->I{
   return min(i,j)<0||i+j>T?I(0):F[i*D+j];};
  auto m=[&](int i,int j)->I{
   return min(i,j)<0||max(i,j)>T?I(0):M[i*D+j];};
  auto H=[&](int i,int j)->I{
   if(min(i,j)<0)return I(0);
   if(i<j)swap(i,j);return I(m(i,j)-m(i+1,j-1));};
  auto P=[&](int i)->I{return H(i,i);};
  auto J=[&](int k,int i,int j)->Q{
   if(min(i,j)<0)return Q(0);
   if(i<j)swap(i,j);
   int t=i-j;Q z(k+t,2),f(k-t,2);
   z.canonicalize();f.canonicalize();Q c=1/(z+1),v=0;
   for(int s=0;s<=j;s++){v+=c*H(i+s,j-s);c*=(f-s)/(z+s+2);}
   return v;
  };
  I S=0,X=0,Y=0,Z=0;
  for(int c=a-b;c<=a+b;c+=2){
   for(int s=abs(p-c);s<=p+c;s+=2)S+=G(s,0);
   Z+=ea*eb*G(p,c);
  }
  for(int s=p-a;s<=p+a;s+=2)X+=eb*G(s,b);
  for(int s=p-b;s<=p+b;s+=2)Y+=ea*G(s,a);
  I S2=0;
  for(int j=B0;j<=d;j++)S2+=P(j);
  for(int j=L;j<A;j++)S2-=P(j);
  I X2=eb*(H(d,B0)-H(A-1,L));
  I Y2=ea*(H(d,A)-H(B0-1,L));
  I Z2=ea*eb*(m(B0,A)-m(d+1,L)-m(B0-1,A-1)+m(d,L-1));
  I parent=FB[p*(W+1)],child=G(p,0);
  if(S!=S2||X!=X2||Y!=Y2||Z!=Z2||
     parent!=S+X+Y+Z||child!=P(h)-P(h-1))++bad;
  Q u=(b+1)*(b+1)*J(2*b,B0,B0)+(a+1)*(a+1)*J(2*a,A,A)
       +2*ea*eb*(a+1)*(b+1)*J(a+b,A,B0);
  Q v=0;
  if(L>=0)v=(b+1)*(b+1)*J(2*b,A-1,A-1)
       +(a+1)*(a+1)*J(2*a,B0-1,B0-1)
       +2*ea*eb*(a+1)*(b+1)*J(a+b,A-1,B0-1);
  I q=S+Z-child;
  Q x=Q(P(d))*u,y=Q(P(L))*v,z=Q(q*q)-x-y;
  bool pass=q>=0&&z>=0&&z*z>=4*x*y;
  I star1=eb*H(d,B0)+ea*H(d,A);
  I star2=eb*H(A-1,L)+ea*H(B0-1,L);
  if(x<0||y<0||Q(star1*star1)>x||Q(star2*star2)>y)++bad;
  if(u<0||v<0||Q((X+Y)*(X+Y))>
       (P(d)+P(L))*(u+v))++bad;
  if(pass){++cert;if(parent<child)++bad;}
  if(2*w>=d){++target;if(parent<child)++bad;}
  if(tag==0){++nf;if(!pass||parent<child)++bad;}
  if(tag==1){++cf;if(parent>=child)++bad;}
  if(tag==2){++gf;if(parent>=child)++bad;}
  if(tag==3){
   if(p==8&&(parent!=I("442880519250")||child!=I("67588703325")||
       S!=I("1840899488025")||X!=I("-897733207800")||
       Y!=X||Z!=I("397447446825")||!pass))++bad;
   if(p==6&&(parent!=460357||child!=391286||S!=4677859||
       X!=-3988179||Y!=-2904678||Z!=2675355||pass||
       q!=6961928||x!=Q(I("437959030418208"),7)||
       y!=Q(I("19417224291546"),35)))++bad;
   if(p==26){
    auto lower=C;lower.push_back(-24);lower.push_back(-24);
    auto V=table(lower,false);
    if(parent!=126032636||child!=1||
       V[p*(weight(lower)+1)]!=126032684||!pass)++bad;
   }
  }
  if(tag==4){++rnd;if(w>=d){++rhigh;if(!pass)++rfail;}}
  int finished=++done;
  if(finished%2000==0){
   #pragma omp critical
   {cout<<"audited "<<finished<<"/"<<n<<endl;}
  }
 }
 cout<<"records="<<n<<" noflip="<<nf<<" census_fail="<<cf
     <<" GMP_fail="<<gf<<" ratio_region="<<target<<" certificate="<<cert
     <<" random="<<rnd<<" random_w_ge_delta="<<rhigh
     <<" random_certificate_misses_there="<<rfail<<" errors="<<bad<<endl;
 return bad.load();
}
'''

data=str(len(rows))+"\n"+"\n".join(
 f"{tag} {p} {len(B)} "+" ".join(map(str,B)) for tag,p,B in rows)
obj=os.memfd_create("fm156-object");lib=os.memfd_create("fm156-library")
subprocess.run(["g++","-std=c++17","-O2","-pipe","-fPIC","-fopenmp",
 "-x","c++","-c","-","-o",f"/proc/self/fd/{obj}"],
 input=CPP,text=True,pass_fds=(obj,),check=True)
subprocess.run(["ld","-shared",f"/proc/self/fd/{obj}",
 "/lib/x86_64-linux-gnu/libgmpxx.so.4","/lib/x86_64-linux-gnu/libgmp.so.10",
 "/lib/x86_64-linux-gnu/libstdc++.so.6","/lib/x86_64-linux-gnu/libgomp.so.1",
 "-lc","-o",f"/proc/self/fd/{lib}"],pass_fds=(obj,lib),check=True)
dll=ctypes.CDLL(f"/proc/self/fd/{lib}")
dll.audit.argtypes=[ctypes.c_char_p,ctypes.c_int]
assert dll.audit(data.encode(),args.threads)==0
os.close(obj);os.close(lib)
print("ALL CHECKS PASS")
