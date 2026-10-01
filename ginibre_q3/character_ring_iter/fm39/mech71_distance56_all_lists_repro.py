import sympy as S, time
from itertools import combinations_with_replacement
from math import comb, factorial
from functools import lru_cache
t,b,k,x=S.symbols('t b k x')
Q=S.Rational
def choose(v,j):
 if j<0:return S.Integer(0)
 return S.prod(v-i for i in range(j))*Q(1,factorial(j))
A=[S.Integer(1),x];P=[S.Integer(1),x]
for j in range(1,12):
 A.append(S.expand((x*A[j]-(t-j+1)*A[j-1])/(j+1)))
 P.append(S.expand((x*P[j]-(t-2*b-j+1)*P[j-1])/(j+1)))
@lru_cache(None)
def moment(j,word):
 row={0:1}
 for n in (1,)*j+word:
  new={}
  for i,v in row.items():
   for q in range(abs(i-n),i+n+1,2):new[q]=new.get(q,0)+v
  row=new
 return row.get(0,0)
@lru_cache(None)
def rempoly(r,R,j):
 return S.expand(sum(choose(k-r+h-2,h)*choose(t-j,R-h)
                     for h in range(R+1)))
@lru_cache(None)
def kernel(d,word):
 w=sum(word);r=len(word);out=0
 for j in range(w%2,2*d-w+1,2):
  D=d-(w+j)//2
  for v in range(D//2+1):
   for u in range(D-2*v+1):
    m=moment(j+2*u,word)
    if m:
     out+=A[j]*choose(b,u)*choose(b-u,v)*m*rempoly(r,D-u-2*v,j)
 return S.expand(out)
def diagonal(poly,d):
 rem=poly;out={}
 for j in range(d,-1,-1):
  out[j]=S.factor(rem.coeff(x,2*j)*factorial(j)**2)
  rem=S.expand(rem-out[j]*P[j]**2)
 assert rem==0
 return out

import resource
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
import os,ctypes,subprocess,glob
def mem(name,data=b''):
 fd=os.memfd_create(name,0)
 if data:os.write(fd,data)
 return fd
def build(src):
 start=time.monotonic()
 sf=mem('fm71.cpp',src.encode())
 cc=subprocess.run(
  ['g++','-std=c++17','-O3','-fPIC','-pthread','-S','-x','c++',
   '-o','-',f'/proc/self/fd/{sf}'],
  pass_fds=(sf,),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 if cc.returncode:raise RuntimeError(cc.stderr.decode())
 af=mem('fm71.s',cc.stdout);of=mem('fm71.o')
 rr=subprocess.run(
  ['as','--64','-o',f'/proc/self/fd/{of}',f'/proc/self/fd/{af}'],
  pass_fds=(of,af),capture_output=True)
 if rr.returncode:raise RuntimeError(rr.stderr.decode())
 so=mem('fm71.so')
 gcc=sorted(glob.glob('/usr/lib/gcc/x86_64-linux-gnu/*/libgcc.a'))[-1]
 gb=os.path.dirname(gcc)
 rr=subprocess.run(
  ['ld.gold','-shared','-o',f'/proc/self/fd/{so}',gb+'/crtbeginS.o',
   f'/proc/self/fd/{of}',gcc,'/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
   '-lpthread','-lc','-lm',gb+'/crtendS.o'],
  pass_fds=(of,so),capture_output=True)
 if rr.returncode:raise RuntimeError(rr.stderr.decode())
 lib=ctypes.CDLL(f'/proc/self/fd/{so}')
 print('memory compilation seconds',round(time.monotonic()-start,2),flush=True)
 return lib

CPP_HEADER=r'''
#include <boost/multiprecision/cpp_int.hpp>
#include <pthread.h>
#include <atomic>
#include <stdio.h>
#include <stdint.h>
#include <algorithm>
#include <stdlib.h>
using I=__int128;
using boost::multiprecision::cpp_int;
static I absI(I a){return a<0?-a:a;}
static I binom(int n,int r){
 if(r<0||r>n)return 0;
 I a=1;for(int i=1;i<=r;i++)a=a*(n-i+1)/i;return a;
}
static cpp_int big(I x){
 bool neg=x<0;unsigned __int128 u=neg?-x:x;
 cpp_int a=(uint64_t)(u>>64);a<<=64;a+=(uint64_t)u;
 return neg?-a:a;
}
'''
CPP_SUFFIX=r'''
static void evaluate(int d,int t,int b,int k,I*v){
 if(d==5)eval5(t,b,k,v);else eval6(t,b,k,v);
}
static I field(I*v,int j){return j<0?0:v[j];}
static bool covers(int d,int t,int b,int k,I*v){
 I D[7]={},O[6]={},E[7]={};
 int nr=d==5?11:22;
 for(int n=0;n<2*d;n++){
  I u=0;
  for(int r=1;r<=2*d/3;r++){
   I mx=0;
   for(int a=0;a<nr;a++){
    int rr=d==5?wr5[a]:wr6[a];if(rr!=r)continue;
    int id=d==5?ker5[a][n]:ker6[a][n];
    mx=std::max(mx,absI(field(v,id)));
   }
   u+=binom(k,r)*mx;
  }
  if(n%2)O[n/2]=u;else E[n/2]=u;
 }
 for(int j=0;j<=d;j++){
  I mx=0;
  for(int i=0;i<d-3;i++){
   int id=d==5?dep5[i][j]:dep6[i][j];
   mx=std::max(mx,absI(field(v,id)));
  }
  E[j]+=k*mx;
  D[j]=field(v,d==5?aa5[j]:aa6[j])-E[j];
  if(D[j]<=0)return false;
 }
 cpp_int prev=1,cur=big(2*D[d]);
 for(int j=d-1;j>=0;j--){
  cpp_int off=big(O[j]);
  cpp_int nxt=big(2*D[j])*cur-off*off*prev;
  if(nxt<=0)return false;
  prev=cur;cur=nxt;
 }
 return true;
}
static void profile(int d,int t,int b,int k,int x,
                    I*v,I&base,I*cc,I*dd){
 I p[7]={1,x};int scale=d==5?720:5040;
 for(int j=1;j<d;j++){
  I num=x*p[j]-(t-2*b-j+1)*p[j-1];
  if(num%(j+1))exit(1);
  p[j+1]=num/(j+1);
 }
 base=0;
 for(int j=0;j<=d;j++)
  base+=field(v,d==5?aa5[j]:aa6[j])*p[j]*p[j];
 if(base%scale)exit(1);base/=scale;
 for(int i=0;i<d-3;i++){
  I z=0;
  for(int j=0;j<=d;j++)
   z+=field(v,d==5?dep5[i][j]:dep6[i][j])*p[j]*p[j];
  if(z%scale)exit(1);dd[i]=z/scale;
 }
 for(int a=0;a<(d==5?11:22);a++){
  I z=0;
  for(int n=0;n<=2*d;n++)
   z+=field(v,d==5?ker5[a][n]:ker6[a][n])*p[n/2]*p[(n+1)/2];
  if(z%scale)exit(1);cc[a]=z/scale;
 }
}
static bool scalar(int d,int k,I base,I*cc,I*dd){
 I dep=0,single=0,dpos=0,dneg=0,off=0,triple=0;
 for(int i=0;i<d-3;i++)dep=std::max(dep,dd[i]);
 for(int a=0;a<(d==5?11:22);a++){
  int r=d==5?wr5[a]:wr6[a];
  const int*w=d==5?wn5[a]:wn6[a];I z=cc[a];
  if(r==1)single=std::max(single,absI(z));
  else if(r==2){
   if(w[0]==w[1]){
    dpos=std::max(dpos,z);dneg=std::max(dneg,-z);
   }else off=std::max(off,absI(z));
  }else if(r==3)triple=std::max(triple,absI(z));
  else if(z!=4)exit(1);
 }
 I lb=2*base-2*k*dep-2*k*single-k*dpos
      -I(k)*k*std::max(dneg,off)-2*binom(k,3)*triple;
 if(d==6&&k>=4)lb-=2*I(k)*k;
 return lb>=0;
}
static I floordiv(I a,I b){
 I q=a/b,r=a%b;if(r<0)q--;return q;
}
static int clip(I q,int J){return q<0?0:q>J?J:(int)q;}
static I ec(int c,int s,int r){
 if(r==0)return 1;if(r==1)return s;
 I a=1,b=s;
 for(int j=1;j<r;j++){
  I num=s*b-(c-j+1)*a;if(num%(j+1))exit(1);
  I z=num/(j+1);a=b;b=z;
 }
 return b;
}
static I fullvalue(int d,int k,I base,I*cc,I*dd,
                   const int*counts,const int*signs){
 I ans=base;
 for(int i=0;i<d-3;i++)ans-=counts[i]*dd[i];
 for(int a=0;a<(d==5?11:22);a++){
  int mult[4]={};int r=d==5?wr5[a]:wr6[a];
  const int*w=d==5?wn5[a]:wn6[a];
  for(int j=0;j<r;j++)mult[w[j]-3]++;
  I term=cc[a];
  for(int i=0;i<d-2;i++)term*=ec(counts[i],signs[i],mult[i]);
  ans+=term;
 }
 return ans;
}
static I getcc(int d,I*cc,int key){
 const int*ids=d==5?id5:id6;
 for(int j=0;j<(d==5?11:22);j++)if(ids[j]==key)return cc[j];
 return 0;
}
struct Counts{long long brute,quad,points;I min;};
static void check(I val,int d,int t,int b,int k,int x,Counts&out){
 if(val<0){
  printf("NEGATIVE d=%d t=%d b=%d k=%d x=%d\n",d,t,b,k,x);
  fflush(stdout);exit(1);
 }
 out.min=std::min(out.min,val);
}
static void brute(int d,int t,int b,int k,int x,
                  I base,I*cc,I*dd,Counts&out){
 int ns=d-2,cs[4]={},ss[4]={};
 auto recur=[&](auto&&self,int a,int left,int weight)->void{
  if(a==ns){
   if(t+2*b+weight+(d+1)*left<2*d)return;
   I val=fullvalue(d,k,base,cc,dd,cs,ss);
   check(6*val,d,t,b,k,x,out);out.brute++;return;
  }
  for(int c=0;c<=left;c++){
   cs[a]=c;
   for(int s=-c;s<=c;s+=2){
    ss[a]=s;self(self,a+1,left-c,weight+(a+3)*c);
   }
  }
 };
 recur(recur,0,k,0);
}
static void finite6(int t,int b,int k,int x,
                    I base,I*cc,I*dd,Counts&out){
 if(t+2*b+3*k<12){brute(6,t,b,k,x,base,cc,dd,out);return;}
 I c3=getcc(6,cc,3),c4=getcc(6,cc,4),c5=getcc(6,cc,5),
   c6=getcc(6,cc,6),c33=getcc(6,cc,33),c34=getcc(6,cc,34),
   c35=getcc(6,cc,35),c36=getcc(6,cc,36),c44=getcc(6,cc,44),
   c45=getcc(6,cc,45),c46=getcc(6,cc,46),c55=getcc(6,cc,55),
   c56=getcc(6,cc,56),c66=getcc(6,cc,66),
   c333=getcc(6,cc,333),c334=getcc(6,cc,334),
   c335=getcc(6,cc,335),c336=getcc(6,cc,336),
   c344=getcc(6,cc,344),c345=getcc(6,cc,345),
   c444=getcc(6,cc,444),c3333=getcc(6,cc,3333);
 if(c66!=1||c336!=1||c345!=1||c444!=1||c3333!=4)exit(1);
 for(int y=-k;y<=k;y++)
 for(int z=-(k-abs(y));z<=k-abs(y);z++){
  int rem=k-abs(y)-abs(z);
  for(int w=-rem;w<=rem;w++){
   int cy=abs(y),fz=abs(z),gw=abs(w),M=rem-abs(w);
   I qf=2*dd[1]+c44+y*c344+z*c444;
   I qg=2*dd[2]+c55;
   I qm=std::max(I(1),std::max(qf,qg));
   int which=qf==qm?1:qg==qm?2:3;
   I H=2*dd[0]+c33+y*c333+z*c334+w*c335+2*I(y)*y+2;
   I dv=2*(c6+y*c36+z*c46+w*c56)+I(y)*y;
   I C0=6*base+6*(y*c3+z*c4+w*c5)
       +3*I(y)*y*c33+6*I(y)*z*c34+6*I(y)*w*c35
       +3*I(z)*z*c44+6*I(z)*w*c45+3*I(w)*w*c55
       +(I(y)*y*y+2*y)*c333+3*I(y)*y*z*c334
       +3*I(y)*y*w*c335+3*I(y)*z*z*c344
       +6*I(y)*z*w*c345+(I(z)*z*z+2*z)*c444
       +I(y)*y*y*y+8*I(y)*y;
   for(int q=0;q<2;q++){
    if(q>M)continue;
    int r=(M-q)%2,J=(M-q)/2;
    for(int sg:{-1,1}){
     I li=12*cy-6*H-6*sg*q+6*qm;
     I lj=-6*sg*cy+6*dv*sg+12*q-6+6*qm;
     I least=(I(1)<<120);
     auto test=[&](I ii,I jj){
      if(ii<0||jj<0||ii+jj>J)return;
      int i=(int)ii,j=(int)jj,u=q+2*j,c=cy+2*i,v=sg*u;
      I val=C0+3*I(c)*c-3*H*c-3*I(v)*c
             +3*dv*v+3*I(u)*u-3*qf*fz-3*qg*gw
             -3*u-3*qm*(M-2*i-u-r);
      least=std::min(least,val);out.points++;
      int cs[4]={c,fz,gw,u},ss[4]={y,z,w,v};
      cs[which]+=M-2*i-u-r;
      if(val!=6*fullvalue(6,k,base,cc,dd,cs,ss)){
       printf("QP_BRIDGE_FAIL\n");fflush(stdout);exit(1);
      }
     };
     test(0,clip(floordiv(12-lj,24),J));
     test(clip(floordiv(12-li,24),J),0);
     I aa=12*(2+sg),lin=li-lj-aa*J;
     int edge=clip(floordiv(aa-lin,2*aa),J);
     test(edge,J-edge);
     I vi=floordiv(-2*li-sg*lj,36);
     I vj=floordiv(-sg*li-2*lj,36);
     for(int di=-1;di<=2;di++)
      for(int dj=-1;dj<=2;dj++)test(vi+di,vj+dj);
     check(least,6,t,b,k,x,out);out.quad++;
    }
   }
  }
 }
}
struct Task{int d,tid,T,B,K;};
struct Result{long long good,bad,scalars,profiles;};
static Result res[4];
static Counts cnt[4];
static std::atomic<int> nextt(0);
static void* work(void*arg){
 Task*a=(Task*)arg;Result&r=res[a->tid];
 for(;;){
  int t=nextt.fetch_add(1);if(t>=a->T)break;
  for(int b=0;b<a->B;b++)for(int k=0;k<a->K;k++){
   if(t+2*b+(a->d+1)*k<2*a->d)continue;
   I v[80];evaluate(a->d,t,b,k,v);
   if(covers(a->d,t,b,k,v)){r.good++;continue;}
   r.bad++;
   for(int x=t%2;x<=t;x+=2){
    I base,cc[22],dd[3];
    profile(a->d,t,b,k,x,v,base,cc,dd);
    if(scalar(a->d,k,base,cc,dd)){r.scalars++;continue;}
    r.profiles++;
    if(a->d==5)brute(5,t,b,k,x,base,cc,dd,cnt[a->tid]);
    else finite6(t,b,k,x,base,cc,dd,cnt[a->tid]);
   }
  }
 }
 return nullptr;
}
extern "C" void run(int d){
 for(auto&r:res)r=Result{};
 for(auto&r:cnt){r=Counts{};r.min=(I(1)<<120);}
 nextt=0;
 int T=d==5?59:129,B=d==5?32:73,K=d==5?56:152;
 pthread_t th[4];Task a[4];
 for(int i=0;i<4;i++){
  a[i]={d,i,T,B,K};pthread_create(&th[i],nullptr,work,&a[i]);
 }
 for(int i=0;i<4;i++)pthread_join(th[i],nullptr);
 Result z{};Counts zz{};zz.min=(I(1)<<120);
 for(auto&r:res){
  z.good+=r.good;z.bad+=r.bad;
  z.scalars+=r.scalars;z.profiles+=r.profiles;
 }
 for(auto&r:cnt){
  zz.brute+=r.brute;zz.quad+=r.quad;
  zz.points+=r.points;zz.min=std::min(zz.min,r.min);
 }
 printf("FINITE d=%d brute=%lld quadratics=%lld points=%lld minimum6F=%lld\n",
        d,zz.brute,zz.quad,zz.points,(long long)zz.min);
 printf("GRID d=%d matrix=%lld scalar=%lld remaining=%lld\n",
        d,z.good,z.scalars,z.profiles);
 fflush(stdout);
}
extern "C" const char* value(int d,int t,int b,int k,int x,
                            const int*cs,const int*ss){
 I f[80],base,cc[22],dd[3];
 evaluate(d,t,b,k,f);profile(d,t,b,k,x,f,base,cc,dd);
 I n=2*fullvalue(d,k,base,cc,dd,cs,ss);
 static thread_local char out[128];
 int i=0;bool neg=n<0;if(neg)n=-n;
 do{out[i++]='0'+n%10;n/=10;}while(n);
 if(neg)out[i++]='-';
 std::reverse(out,out+i);out[i]=0;return out;
}
'''

from sympy.printing.c import C99CodePrinter
from math import lcm
import argparse
argparse.ArgumentParser(
 description="FM-MECH71 exact distance-5/6 proof").parse_args()
begun=time.monotonic()
bases={d:kernel(d,()) for d in range(7)}
diags={d:diagonal(bases[d],d) for d in range(7)}
cores={}
for d in (5,6):
 cores[d]={}
 for r in range(1,2*d//3+1):
  for word in combinations_with_replacement(range(3,d+1),r):
   if sum(word)<=2*d:cores[d][word]=kernel(d,word)
def balanced(pol):
 rem=pol;out={}
 for n in range(int(S.degree(rem,x)),-1,-1):
  v=S.expand(rem).coeff(x,n)*factorial(n//2)*factorial((n+1)//2)
  if v:
   out[n]=S.factor(v)
   rem=S.expand(rem-v*P[n//2]*P[(n+1)//2])
 assert rem==0
 return out
bal={d:{w:balanced(p) for w,p in cores[d].items()} for d in (5,6)}
def upper(polys):
 ds=[S.Poly(p,t,b,k).as_dict() for p in polys]
 mons=set().union(*(q.keys() for q in ds))
 return sum(max(abs(q.get(m,0)) for q in ds)
            *t**m[0]*b**m[1]*k**m[2] for m in mons)
weights={5:[Q(32,5),Q(16,5),Q(16,5),2,1],
         6:[Q(48,5),Q(24,5),Q(24,5),4,Q(16,7),1]}
tails={5:(59,32,56),6:(129,73,152)}
for d in (5,6):
 E=[S.Integer(0)]*(d+1);O=[S.Integer(0)]*d
 for n in range(2*d):
  bound=0
  for r in range(1,2*d//3+1):
   qs=[qs.get(n,0) for w,qs in bal[d].items() if len(w)==r]
   if qs:bound+=k**r*upper(qs)/factorial(r)
  if n%2:O[n//2]=S.expand(bound)
  else:E[n//2]=S.expand(bound+k*upper(
       [diags[d-i-1].get(n//2,0) for i in range(3,d)]))
 R=[S.expand(diags[d][j]-E[j]
       -(weights[d][j]*O[j]/2 if j<d else 0)
       -(O[j-1]/(2*weights[d][j-1]) if j else 0))
       for j in range(d+1)]
 for var,H in zip((t,b,k),tails[d]):
  assert all(min(S.Poly(p.subs(var,var+H),t,b,k).coeffs())>=0
             for p in R)
 print("polynomial tails",d,tails[d],"PASS",flush=True)

for d in range(7):
 D=[]
 for i in range(d+1):
  vv=0
  for v in range((d-i)//2+1):
   for u in range(d-i-2*v+1):
    vv+=choose(b,u)*choose(b-u,v)*Q(comb(2*(i+u),i+u),i+u+1) \
        *rempoly(0,d-i-u-2*v,2*i)
  D.append(S.expand(vv))
 hs=[S.expand(sum((-1)**(i-r)*choose(b,i-r)*D[i]
                 for i in range(r,d+1))) for r in range(d+1)]
 for j in range(d+1):
  val=sum((-1)**(r-j)*(choose(t-2*b-j-r,r-j)
           +choose(t-2*b-j-r-1,r-j-1))*hs[r]/comb(2*r,r)
          for r in range(j,d+1))
  assert S.expand(val-diags[d][j])==0
print("general diagonal formula through d=6 PASS",flush=True)
for d in (5,6):
 assert diags[d][d-2].subs({t:2*d-2,b:1,k:1}) \
        == -Q(6*(d-4),d*(d+1))
def gi(n,r):
 if r<0:return 0
 a=1
 for j in range(r):a=a*(n-j)//(j+1)
 return a
for M in range(-12,25):
 for j in range(13):
  for ss in range(j,13):
   assert sum((-1)**(r-j)*(gi(M-j-r,r-j)+gi(M-j-r-1,r-j-1))
              *gi(M-2*r,ss-r) for r in range(j,ss+1))==(j==ss)
print("inverse matrix through index 12 PASS",flush=True)

class Printer(C99CodePrinter):
 def _print_Pow(self,ex):
  assert ex.exp.is_Integer and ex.exp>=0
  return '('+'*'.join(['('+self._print(ex.base)+')']*int(ex.exp))+')' \
         if ex.exp else '1'
pr=Printer()
def arr(a):
 return '{'+','.join(arr(v) if isinstance(v,list) else str(v)
                     for v in a)+'}'
source=CPP_HEADER
for d in (5,6):
 fields=[]
 def add(p):
  p=S.expand(p)
  if p==0:return -1
  fields.append(p);return len(fields)-1
 aa=[add(diags[d][j]) for j in range(d+1)]
 deps=[[add(diags[d-i-1].get(j,0)) for j in range(d+1)]
       for i in range(3,d)]
 words=list(cores[d])
 kr=[[add(bal[d][w].get(n,0)) for n in range(2*d+1)] for w in words]
 den=lcm(*(int(v.q) for p in fields for v in S.Poly(p,t,b,k).coeffs()))
 assert den==(720 if d==5 else 5040)
 expr=[S.expand(den*p) for p in fields]
 assert all(v.q==1 for p in expr for v in S.Poly(p,t,b,k).coeffs())
 rr,ee=S.cse(expr)
 code='\n'.join(['I '+str(v)+'='+pr.doprint(q)+';' for v,q in rr]
      +['v['+str(i)+']='+pr.doprint(q)+';' for i,q in enumerate(ee)])
 source+=f'\nstatic const int aa{d}[{d+1}]={arr(aa)};\n'
 source+=f'static const int dep{d}[{d-3}][{d+1}]={arr(deps)};\n'
 source+=f'static const int ker{d}[{len(words)}][{2*d+1}]={arr(kr)};\n'
 source+=f'static const int wr{d}[{len(words)}]={arr([len(w) for w in words])};\n'
 source+=f'static const int wn{d}[{len(words)}][4]={arr([list(w)+[0]*(4-len(w)) for w in words])};\n'
 source+=f'static const int id{d}[{len(words)}]={arr([int("".join(map(str,w))) for w in words])};\n'
 source+=f'static void eval{d}(I t,I b,I k,I*v){{\n{code}\n}}\n'
source+=CPP_SUFFIX
lib=build(source)
lib.run.argtypes=(ctypes.c_int,)
for d in (5,6):lib.run(d)

from collections import defaultdict
from random import Random
lib.value.argtypes=(ctypes.c_int,)*5+(ctypes.POINTER(ctypes.c_int),)*2
lib.value.restype=ctypes.c_char_p
def cg(word):
 out={(0,0):1}
 for tok in word:
  n=abs(tok);sg=1 if tok>0 else -1;nxt=defaultdict(int)
  for (i,j),v in out.items():
   for q in range(abs(i-n),i+n+1,2):nxt[q,j]+=v
   for q in range(abs(j-n),j+n+1,2):nxt[i,q]+=sg*v
  out={ij:v for ij,v in nxt.items() if v}
 return out.get((0,0),0)
def check_word(word):
 rest=list(word)
 ix=max(range(len(rest)),key=lambda i:abs(rest[i]))
 n=abs(rest.pop(ix))
 d=(sum(map(abs,rest))-n)//2
 assert d in (5,6)
 t=sum(abs(v)==1 for v in rest)+2*rest.count(-2)
 b=rest.count(2);x=rest.count(1)-rest.count(-1)
 k=sum(abs(v)>=3 for v in rest)
 cs=(ctypes.c_int*4)(*[rest.count(i)+rest.count(-i) for i in range(3,7)])
 ss=(ctypes.c_int*4)(*[rest.count(i)-rest.count(-i) for i in range(3,7)])
 ans=int(lib.value(d,t,b,k,x,cs,ss))
 actual=cg(word)
 assert ans==actual,(word,ans,actual)
 assert ans>=0
 return ans
count=0
for L in range(1,9):
 for word in combinations_with_replacement((-1,1,-2,2,-3,3,-4,4),L):
  W=sum(map(abs,word))
  if W%2 or sum(v<0 for v in word)%2:continue
  d=W//2-max(map(abs,word))
  if d in (5,6):check_word(word);count+=1
assert count==1101
print('independent small-list bridges',count,flush=True)
rng=Random(7106);large=0
while large<80:
 d=rng.choice((5,6))
 rest=[rng.choice((-1,1))*rng.randrange(1,19)
       for _ in range(rng.randrange(3,9))]
 n=sum(map(abs,rest))-2*d
 if n<max(map(abs,rest)):continue
 word=rest+[(-1)**sum(v<0 for v in rest)*n]
 check_word(word);large+=1
print('independent large-label bridges',large,flush=True)
for d,expected in ((5,136),(6,400)):
 word=[-1]*(d-1)+[1]*(d-1)+[2,d+1,(-1)**(d-1)*(d+1)]
 assert check_word(word)==expected
 print('negative diagonal coefficient, positive word',d,expected,flush=True)
print('ALL CHECKS PASS; seconds',round(time.monotonic()-begun,2),flush=True)
