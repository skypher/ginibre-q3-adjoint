import argparse,ctypes,os,pathlib,re,subprocess
from collections import defaultdict
from fractions import Fraction as Q
from math import comb

ap=argparse.ArgumentParser(description="Exact read-only FM-MECH161 checks.")
ap.add_argument("--no-census",action="store_true")
ap.add_argument("--threads",type=int,default=12)
ap.add_argument("--logs",default="/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc")
args=ap.parse_args()
assert 1<=args.threads<=32

def kernel(B,deform):
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
   for (a,b,c),x in F.items():O[i+a,j+b,k+c]+=v*x
  M={key:v for key,v in O.items() if v}
 return M

checks=0
for B in ((1,-2,3),(-1,2,-3,4),(1,1,-2,3,-4),(-2,-4,5)):
 T=sum(map(abs,B));M=kernel(B,False);V=kernel(B,True)
 def H(i,j):
  i,j=max(i,j),min(i,j)
  return M.get((i,j,0),0)-M.get((i+1,j-1,0),0)
 for k in range(13):
  direct=defaultdict(Q)
  for (i,j,r),v in V.items():direct[i,j]+=Q(2*v,k+r+2)
  for i in range(T+1):
   for j in range(i+1):
    alpha=Q(k-i+j,2);beta=Q(k+i-j,2)
    c=1/(beta+1);v=0
    for t in range(j+1):
     v+=c*H(i+t,j-t);c*=Q(alpha-t,beta+t+2)
    assert v==direct[i,j];checks+=1
print("direct moment checks",checks,"PASS")

rows=[]
base=pathlib.Path(args.logs)
if not args.no_census:
 for line in (base/"fx3b_48.log").read_text().splitlines():
  if not line.startswith("NOFLIP "):continue
  B=tuple(map(int,re.search(r"\bB=(.*?)(?:  |$)",line)[1].split()))
  p=abs(int(re.search(r"\bp=(-?\d+)",line)[1]))
  rows.append((0,p,B))
 assert len(rows)==33487
 for name in ("fx5_40.log","tpg_one.log","tpg_all.log"):
  for line in (base/name).read_text().splitlines():
   if not line.startswith("TPFAIL "):continue
   B=tuple(map(int,re.search(r"\bB=(.*?)(?:  |$)",line)[1].split()))
   p=abs(int(re.search(r"\bp=(-?\d+)",line)[1]))
   rows.append((1,p,B))
 assert sum(tag==1 for tag,p,B in rows)==190
B=(1,1,-2)+(3,)*4+(-4,)*3+(5,)*4
for p,C in ((8,B),(6,B+(-2,)),(8,B+(8,))):
 for D in (C,tuple(-z if z%2 else z for z in C)):
  rows.append((2,p,D))
rows.append((2,62,(-40,42,-44,46,48,50,52,54,56,58,60)))
rows.append((3,6,(-1,)*13+(-2,)+(-3,)*5+(-4,)))
for tag,p,B in rows:
 W=sum(map(abs,B));d=(W-p)//2;sig=(-1)**sum(z<0 for z in B)
 assert len(B)<=1200 and all(0<abs(z)<=1200 for z in B) and W<=1200
 assert p>=max(6,max(map(abs,B))) and (W-p)%2==0
 assert d>=max(8,max(map(abs,B))) and sum(abs(z)>=3 for z in B)>=2
 assert -sig*p not in B and all(-z not in B for z in B)

CPP=r'''
#include <bits/stdc++.h>
#include <gmpxx.h>
#include <omp.h>
using namespace std;using I=mpz_class;using Q=mpq_class;
int weight(const vector<int>&B){int w=0;for(int z:B)w+=abs(z);return w;}
vector<I> table(const vector<int>&B,bool torus){
 int W=weight(B),D=W+1,s=0;vector<I>F(D*D),O(D*D);F[0]=1;
 for(int z:B){
  int n=abs(z),e=z>0?1:-1;fill(O.begin(),O.end(),0);
  for(int i=0;i<=s;i++)for(int j=0;j<=(torus?s:s-i);j++){
   const I&v=F[i*D+j];if(v==0)continue;
   if(torus)for(int r=0;r<=n;r++){
    O[(i+r)*D+j+r]+=v;O[(i+r)*D+j+n-r]+=e*v;
   }else{
    for(int r=abs(i-n);r<=i+n;r+=2)O[r*D+j]+=v;
    for(int r=abs(j-n);r<=j+n;r+=2)O[i*D+r]+=e*v;
   }
  }F.swap(O);s+=n;
 }return F;
}
I sqrt_up(const Q&v){
 I a=v.get_num(),b=v.get_den(),t=a/b,r;
 mpz_sqrt(r.get_mpz_t(),t.get_mpz_t());
 if(r*r*b<a)++r;return r;
}
struct Row{int tag,p;vector<int>B;};
extern "C" int audit(const char*data,int threads){
 omp_set_num_threads(threads);istringstream in(data);int count;in>>count;
 if(!in||count<0||count>50000)return 1;
 vector<Row>rows(count);
 for(auto&r:rows){
  int n;in>>r.tag>>r.p>>n;
  if(!in||n<2||n>1200||r.p<0||r.p>1200)return 2;
  r.B.resize(n);
  for(int&z:r.B){in>>z;if(!in||z==0||z<-1200||z>1200)return 3;}
 }
 atomic<int>bad{0},done{0},oldpass{0},newpass{0},noflip{0},controls{0};
 #pragma omp parallel for schedule(dynamic,1)
 for(int idx=0;idx<count;idx++){
  auto [tag,p,B]=rows[idx];int W=weight(B),N=B.size();
  if(W>1200||p>W||(W-p)%2){++bad;continue;}
  int ia=-1,ib=-1,w=-1,mx=-1;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)if((B[i]-B[j])%2==0){
   auto key=make_pair(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])));
   if(key>make_pair(w,mx))ia=i,ib=j,w=key.first,mx=key.second;
  }
  if(ia<0){++bad;continue;}
  int za=B[ia],zb=B[ib];if(abs(za)<abs(zb))swap(za,zb);
  int a=abs(za),b=abs(zb),ea=za>0?1:-1,eb=zb>0?1:-1;
  auto C=B;C.erase(C.begin()+ib);C.erase(C.begin()+ia);
  int T=weight(C),D=T+1,d=(W-p)/2,r=d-a,s=d-b,l=d-a-b-1,h=d-w/2;
  if(p<a||r<0){++bad;continue;}
  auto F=table(C,false),FB=table(B,false),M=table(C,true);
  auto g=[&](int i,int j)->I{
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
   Q alpha(k-i+j,2),beta(k+i-j,2);
   alpha.canonicalize();beta.canonicalize();Q c=1/(beta+1),v=0;
   for(int t=0;t<=j;t++){v+=c*H(i+t,j-t);c*=(alpha-t)/(beta+t+2);}
   return v;
  };
  I S=0,X=0,Y=0,Z=0;
  for(int c=a-b;c<=a+b;c+=2){
   for(int t=abs(p-c);t<=p+c;t+=2)S+=g(t,0);
   Z+=ea*eb*g(p,c);
  }
  for(int t=p-a;t<=p+a;t+=2)X+=eb*g(t,b);
  for(int t=p-b;t<=p+b;t+=2)Y+=ea*g(t,a);
  I S2=0;for(int j=s;j<=d;j++)S2+=P(j);
  for(int j=l;j<r;j++)S2-=P(j);
  I upper=eb*H(d,s)+ea*H(d,r);
  I wrap=eb*H(r-1,l)+ea*H(s-1,l);
  I Z2=ea*eb*(m(s,r)-m(d+1,l)-m(s-1,r-1)+m(d,l-1));
  I parent=FB[p*(W+1)],child=g(p,0),q=S+Z-child,A=q-wrap;
  if(S!=S2||Z!=Z2||X+Y!=upper-wrap||
     parent!=S+X+Y+Z||child!=P(h)-P(h-1))++bad;
  Q E0=0,best=0;int bestk=0;
  for(int k:set<int>{0,b-b%2,b+b%2}){
   Q u=(b+1)*(b+1)*J(2*b-k,s,s)+(a+1)*(a+1)*J(2*a-k,r,r)
       +2*ea*eb*(a+1)*(b+1)*J(a+b-k,r,s);
   Q e=J(k,d,d)*u;
   if(u<0||e<0||Q(upper*upper)>e)++bad;
   if(k==0)E0=best=e;
   else if(e<best)best=e,bestk=k;
  }
  Q v=0;
  if(l>=0)v=(b+1)*(b+1)*J(2*b,r-1,r-1)
       +(a+1)*(a+1)*J(2*a,s-1,s-1)
       +2*ea*eb*(a+1)*(b+1)*J(a+b,r-1,s-1);
  Q Eminus=Q(P(l))*v,disc=Q(q*q)-E0-Eminus;
  if(v<0||Eminus<0||Q(wrap*wrap)>Eminus)++bad;
  bool old=q>=0&&disc>=0&&disc*disc>=4*E0*Eminus;
  bool improved=A>=0&&Q(A*A)>=best;
  if(old&&!improved)++bad;
  if(improved&&(parent<child||parent-child<A-sqrt_up(best)))++bad;
  if(tag==0){
   ++noflip;if(old)++oldpass;if(improved)++newpass;
   if(!improved||parent<child)++bad;
  }
  if(tag==1){++controls;if(improved||parent>=child)++bad;}
  if(tag==2){
   if(!improved)++bad;
   if(W==48&&(parent!=4075371||child!=160741||old||
      best!=Q(I("77923386713832"),7)||bestk!=6||A!=4459446))++bad;
   if(W==50&&(parent!=10625415||child!=424915||old||
      best!=Q(I("1857470360352348"),25)||bestk!=4||A!=11444760))++bad;
   if(W==56&&(parent!=27576641||child!=961895||old||
      best!=Q(I("113034438518476594"),495)||bestk!=4||A!=28000199))++bad;
   if(W==550&&(parent!=I("453207534222864")||child!=I("179646349605")||
      !old||bestk!=58||A!=I("453122062854658")))++bad;
   #pragma omp critical
   {cout<<"special W="<<W<<" p="<<p<<" old="<<old<<" new="<<improved
      <<" k="<<bestk<<" lower="<<A-sqrt_up(best)<<endl;}
  }
  if(tag==3){
   if(parent!=700607800||child!=723462700||improved||
      X+Y!=I("-13096933340")||Y+Z!=I("-400136390")||
      X+Z!=I("-2291824230"))++bad;
   auto L=B;L.erase(find(L.begin(),L.end(),-1));
   L.erase(find(L.begin(),L.end(),-1));L.push_back(6);
   auto V=table(L,false);int dim=weight(L)+1;
   if(-V[dim+1]!=345935230)++bad;
   #pragma omp critical
   {cout<<"three-flip counterexample PASS"<<endl;}
  }
  int completed=++done;
  if(completed%4000==0){
   #pragma omp critical
   {cout<<"audited "<<completed<<"/"<<count<<endl;}
  }
 }
 cout<<"census="<<noflip<<" old="<<oldpass<<" improved="<<newpass
     <<" failure_controls="<<controls<<" errors="<<bad<<endl;
 return bad.load();
}
'''

data=str(len(rows))+"\n"+"\n".join(
 f"{tag} {p} {len(B)} "+" ".join(map(str,B)) for tag,p,B in rows)
obj=os.memfd_create("fm161-object");lib=os.memfd_create("fm161-library")
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
