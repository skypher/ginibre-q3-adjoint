import argparse,ctypes,os,pathlib,re,resource,subprocess
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
from math import comb

ap=argparse.ArgumentParser(description="Exact, memory-only FM-MECH165 verifier.")
ap.add_argument("--no-census",action="store_true")
ap.add_argument("--threads",type=int,default=12)
ap.add_argument("--log",default="/tmp/claude-1006/-home-yang-q3adjoint/"
 "2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/fx3b_48.log")
args=ap.parse_args()
assert 1<=args.threads<=32
corners=[]
for m in (6,7):
 for x in range(m+1):
  for y in range(m-x+1):
   c=comb(x,2)
   for e in (-1,1):
    P1=m+c
    P2=comb(m+1,2)-x+c*(m-2+e*y)+comb(y,2)+2*comb(x,4)
    for t in (-1,1):corners.append(P2-2*P1-2*t*(e*y+c))
assert len(corners)==256 and min(corners)==0
print("256 base corners PASS")
rows=[]
if not args.no_census:
 for line in pathlib.Path(args.log).read_text().splitlines():
  if line.startswith("NOFLIP "):
   B=tuple(map(int,re.search(r"\bB=(.*?)(?:  |$)",line)[1].split()))
   p=abs(int(re.search(r"\bp=(-?\d+)",line)[1]))
   rows.append((0,p,B))
 assert len(rows)==33487
A=(1,1,-2)+(3,)*4+(-4,)*3+(5,)*4
rows += [
 (1,6,A+(-2,)),(1,8,A+(8,)),
 (1,62,(-40,42,-44,46,48,50,52,54,56,58,60)),
 (2,7,(1,1)+(3,)*4+(-4,)*4+(5,)*5),
 (3,7,(1,1,-2)+(3,)*5+(-4,)*3+(5,)*6),
 (4,6,(-1,)*13+(-2,)+(-3,)*5+(-4,)),
 (5,8,A),(6,13,(1,2,-3,4,4,6,8,-11)),
 (7,6,(1,-2,3,3,-4,-4,5,5,5))]
for tag,p,B in rows:
 W=sum(map(abs,B));d=(W-p)//2;sig=(-1)**sum(z<0 for z in B)
 assert 0<p<=W<=1200 and len(B)<=1200 and (W-p)%2==0
 assert p>=max(6,max(map(abs,B))) and d>=max(8,max(map(abs,B)))
 assert all(z and -z not in B for z in B) and -sig*p not in B

CPP=r'''
#include <bits/stdc++.h>
#include <gmpxx.h>
#include <omp.h>
using namespace std;using I=mpz_class;using Q=mpq_class;
Q rat(I x,I y){Q q(x,y);q.canonicalize();return q;}
int wt(const vector<int>&v){int t=0;for(int z:v)t+=abs(z);return t;}
I choose(int n,int k){if(k>n)return 0;I v=1;for(int j=1;j<=k;j++)v=v*(n-j+1)/j;return v;}
vector<I> table(const vector<int>&B,bool torus){
 int W=wt(B),D=W+1,t=0;vector<I>F(D*D),O(D*D);F[0]=1;
 for(int z:B){int n=abs(z),e=z>0?1:-1;fill(O.begin(),O.end(),0);
  for(int i=0;i<=t;i++)for(int j=0;j<=(torus?t:t-i);j++){
   const I&v=F[i*D+j];if(v==0)continue;
   if(torus)for(int h=0;h<=n;h++){O[(i+h)*D+j+h]+=v;O[(i+h)*D+j+n-h]+=e*v;}
   else{for(int h=abs(i-n);h<=i+n;h+=2)O[h*D+j]+=v;
        for(int h=abs(j-n);h<=j+n;h+=2)O[i*D+h]+=e*v;}
  }F.swap(O);t+=n;
 }return F;
}
struct Kernel{
 int T,D;vector<I>M;
 Kernel(const vector<int>&B):T(wt(B)),D(T+1),M(table(B,true)){}
 I at(int i,int j)const{return min(i,j)<0||max(i,j)>T?I(0):M[i*D+j];}
 I H(int i,int j)const{if(min(i,j)<0)return 0;if(i<j)swap(i,j);return I(at(i,j)-at(i+1,j-1));}
 I P(int j)const{return H(j,j);}
 Q J(int k,int i,int j)const{
  if(min(i,j)<0)return 0;if(i<j)swap(i,j);
  Q al(k-i+j,2),be(k+i-j,2);al.canonicalize();be.canonicalize();
  Q c=1/(be+1),v=0;
  for(int h=0;h<=j;h++){v+=c*H(i+h,j-h);c*=(al-h)/(be+h+2);}return v;
 }
};
struct Result{
 int a,b,d,r,s,m,k;I Q0,upper,S,X,Y,Z,parent,child,pd,p1,p2;Q best,E0;
};
Result evaluate(const vector<int>&B,int p){
 int W=wt(B),n=B.size(),ia=-1,ib=-1,w=-1,mx=-1;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if((B[i]-B[j])%2==0){
  auto key=make_pair(abs(B[i])+abs(B[j]),max(abs(B[i]),abs(B[j])));
  if(key>make_pair(w,mx))ia=i,ib=j,w=key.first,mx=key.second;
 }
 assert(ia>=0);int za=B[ia],zb=B[ib];if(abs(za)<abs(zb))swap(za,zb);
 int a=abs(za),b=abs(zb),ea=za>0?1:-1,eb=zb>0?1:-1,d=(W-p)/2;
 auto C=B;C.erase(C.begin()+ib);C.erase(C.begin()+ia);Kernel K(C);
 int r=d-a,s=d-b,l=d-a-b-1,h=d-w/2,m=C.size();
 Result R{};R.a=a;R.b=b;R.d=d;R.r=r;R.s=s;R.m=m;
 for(int j=s;j<=d;j++)R.S+=K.P(j);
 for(int j=l;j<r;j++)R.S-=K.P(j);
 R.X=eb*(K.H(d,s)-K.H(r-1,l));R.Y=ea*(K.H(d,r)-K.H(s-1,l));
 R.Z=ea*eb*(K.at(s,r)-K.at(d+1,l)-K.at(s-1,r-1)+K.at(d,l-1));
 R.upper=eb*K.H(d,s)+ea*K.H(d,r);
 R.child=K.P(h)-K.P(h-1);R.parent=R.S+R.X+R.Y+R.Z;
 R.Q0=R.parent-R.child-R.upper;R.pd=K.P(d);R.p1=K.P(1);R.p2=K.P(2);
 for(int k=0;k<=2*b;k++){
  Q v=(b+1)*(b+1)*K.J(2*b-k,s,s)+(a+1)*(a+1)*K.J(2*a-k,r,r)
      +2*ea*eb*(a+1)*(b+1)*K.J(a+b-k,r,s);
  Q e=K.J(k,d,d)*v;assert(e>=Q(R.upper*R.upper));
  if(k==0)R.best=R.E0=e;else if(e<R.best)R.best=e,R.k=k;
 }
 int x=0,y=0,e1=1,y2=0,y3=0,y4=0;
 for(int z:C){int e=z>0?1:-1;
  if(abs(z)==1)x++,e1=e;if(abs(z)==2)y++,y2+=e;
  if(abs(z)==3)y3+=e;if(abs(z)==4)y4+=e;
 }
 I c=choose(x,2),B2=y2+c,P2=choose(m+1,2)-x+c*(m-2+y2)+choose(y,2)+2*choose(x,4);
 assert(R.p1==m+c&&R.p2==P2);
 I D2=choose(m+1,2)-x;
 I v0=D2-y2*(m-1)+y4+choose(y,2);
 I v1=2*y2*(m-1)-6*y4+2*c*(m-2)-4*e1*x*y3-4*choose(y,2)-2*y2*c;
 I v2=6*y4+6*e1*x*y3+6*choose(y,2)+6*y2*c+6*choose(x,4);
 assert(v0>=0&&v0+v1+v2>=0);
 if(m>=3)assert(v0+v1/2>=0);
 for(int k=0;k<=3;k++)
  assert(K.J(k,2,2)==rat(2*v0,k+2)+rat(2*v1,k+4)+rat(2*v2,k+6));
 if(m>=6){
  assert(P2>=2*(R.p1+abs(B2)));
  if(w>=2*d-2){
   assert(R.Q0>0&&Q(R.Q0*R.Q0)>=R.E0);
   if(r==0&&s==0){
    assert(R.Q0>=R.pd+1+(d-1)*m&&R.E0==Q(4*(d+1)*R.pd));
   }else if(r==1&&s==1){
    assert(R.Q0>=R.pd+(d-2)*P2&&R.E0<=Q(2*d*R.pd*P2));
   }else{
    assert(r==0&&s==2);
    assert(Q(R.Q0)>=Q(R.pd)+rat(2*d-5,2)*P2+1);
    Q u=(rat(3*(d-1)*(d-1),d+1)+rat(d-1,2))*P2+d+1;
    assert(R.E0<=R.pd*u);
   }
  }
 }
 return R;
}
struct Row{int tag,p;vector<int>B;};
extern "C" int audit(const char*txt,int threads){
 omp_set_num_threads(threads);istringstream in(txt);int nr;in>>nr;
 if(nr<0||nr>50000)return 1;vector<Row>rs(nr);
 for(auto&r:rs){int n;in>>r.tag>>r.p>>n;if(n<2||n>1200)return 2;
  r.B.resize(n);for(int&z:r.B){in>>z;if(!in||z==0||abs(z)>1200)return 3;}
  if(wt(r.B)>1200||r.p<0||r.p>wt(r.B))return 4;
 }
 atomic<int>done{0},census{0},active{0},band{0};Q minrho=1;
 #pragma omp parallel for schedule(dynamic,1)
 for(int z=0;z<nr;z++){
  auto [tag,p,B]=rs[z];auto R=evaluate(B,p);int W=wt(B);
  bool cb=R.Q0>=0&&Q(R.Q0*R.Q0)>=R.best;
  if(tag==0){
   assert(cb);++census;if(B.size()>=8){++active;if(R.a+R.b>=2*R.d-2)++band;}
   Q rho=1-R.best/Q(R.Q0*R.Q0);
   #pragma omp critical
   {if(rho<minrho)minrho=rho;}
  }else{
   auto F=table(B,false);assert(F[p*(W+1)]==R.parent);
   if(tag==2||tag==3||tag==4)assert(!cb);else assert(cb);
   if(W<=70){
    auto L=B;int sig=1;for(int b:B)if(b<0)sig=-sig;L.push_back(sig*p);
    map<pair<int,int>,I>ds;I ma,star,opp;bool first=true,fs=true,fo=true;
    for(int i=0;i<(int)L.size();i++)for(int j=i+1;j<(int)L.size();j++){
     auto key=minmax(L[i],L[j]);I v;
     if(ds.count(key))v=ds[key];
     else{
      auto C=L;C.erase(C.begin()+j);C.erase(C.begin()+i);
      int T=wt(C),u=abs(L[i]),v0=abs(L[j]),q=(T-u-v0)/2;
      auto G=table(C,false);v=u+v0>T?I(0):I((L[j]>0?1:-1)*G[u*(T+1)+v0]);
      Kernel K(C);Q mm=(L[j]>0?1:-1)*(v0+1)*
       (K.J(v0,q+v0,q)-K.J(v0,q+v0-1,q-1));
      assert(Q(v)==mm);ds[key]=v;
     }
     if(first||v>ma)ma=v,first=false;
     if(j==(int)L.size()-1&&(fs||v>star))star=v,fs=false;
     if((L[i]-L[j])%2&&(fo||v>opp))opp=v,fo=false;
    }
    if(tag==1||tag>=5)assert(ma<0);
    if(tag==2)assert(star<0&&R.X+R.Y<0&&ds[make_pair(-4,1)]==2781040&&ma==2781040&&R.parent==38181497&&R.child==1432432);
    if(tag==3)assert(star<0&&opp<0&&R.X+R.Y<0&&ds[make_pair(3,5)]==13932121&&ma==13932121&&R.parent==424148816&&R.child==14744015);
    if(tag==4)assert(R.parent==700607800&&R.child==723462700&&R.X+R.Y<0&&R.X+R.Z<0&&R.Y+R.Z<0&&ma==345935230);
    if(tag==5)assert(ma==-20406&&ds[make_pair(-4,5)]==ma);
    if(tag==6)assert(ma==-1&&ds[make_pair(-3,1)]==ma&&ds[make_pair(-11,13)]==ma);
    if(tag==7)assert(ma==-11&&ds[make_pair(-4,5)]==ma);
    #pragma omp critical
    {cout<<"W="<<W<<" p="<<p<<" CB="<<cb<<" k="<<R.k<<" Q="<<R.Q0
     <<" parent="<<R.parent<<" child="<<R.child<<" maxD="<<ma
     <<" maxStar="<<star<<" maxOpp="<<opp<<" D_TP="<<R.X+R.Y<<endl;}
   }else{
    assert(W==550&&R.parent==I("453207534222864")&&R.child==I("179646349605"));
    #pragma omp critical
    {cout<<"W=550 CB PASS k="<<R.k<<endl;}
   }
  }
  int n=++done;if(n%6000==0){
   #pragma omp critical
   {cout<<"audited "<<n<<"/"<<nr<<endl;}}
 }
 if(census){
  assert(census==33487&&active==12436&&band==1024);
  assert(minrho==Q(I("93926567757613893746579518725304147657"),
                   I("208726237512666607346314528948441436025")));
  cout<<"census="<<census<<" active="<<active<<" band="<<band
      <<" minimum_margin="<<minrho<<endl;
 }
 return 0;
}
'''
obj=os.memfd_create("fm165-object");lib=os.memfd_create("fm165-library")
subprocess.run(["g++","-std=c++17","-O2","-pipe","-fPIC","-fopenmp",
 "-x","c++","-c","-","-o",f"/proc/self/fd/{obj}"],
 input=CPP,text=True,pass_fds=(obj,),check=True)
subprocess.run(["ld","-shared",f"/proc/self/fd/{obj}",
 "/lib/x86_64-linux-gnu/libgmpxx.so.4","/lib/x86_64-linux-gnu/libgmp.so.10",
 "/lib/x86_64-linux-gnu/libstdc++.so.6","/lib/x86_64-linux-gnu/libgomp.so.1",
 "-lc","-o",f"/proc/self/fd/{lib}"],pass_fds=(obj,lib),check=True)
dll=ctypes.CDLL(f"/proc/self/fd/{lib}")
dll.audit.argtypes=[ctypes.c_char_p,ctypes.c_int]
data=str(len(rows))+"\n"+"\n".join(
 f"{t} {p} {len(B)} "+" ".join(map(str,B)) for t,p,B in rows)
assert dll.audit(data.encode(),args.threads)==0
os.close(obj);os.close(lib)
print("ALL CHECKS PASS")
