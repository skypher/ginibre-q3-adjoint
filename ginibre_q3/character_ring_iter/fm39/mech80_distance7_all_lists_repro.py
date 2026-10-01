import argparse
argparse.ArgumentParser(description='FM-MECH80: exact repair and distance-seven certificate; memory only').parse_args()
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
for j in range(1,15):
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

started=time.monotonic()
def repair6():
 import sympy as S
 from itertools import combinations_with_replacement
 from math import factorial
 c,f,g,h,y,z,w,v=S.symbols('c f g h y z w v')
 base,d0,d1,d2=S.symbols('B d0 d1 d2')
 words=[ww for r in range(1,5) for ww in combinations_with_replacement(range(3,7),r) if sum(ww)<=12]
 C={ww:S.Symbol('C'+''.join(map(str,ww))) for ww in words}
 for ww,val in [((6,6),1),((3,3,6),1),((3,4,5),1),((4,4,4),1),((3,3,3,3),4)]:C[ww]=S.Integer(val)
 def ec(n,s,r):
  a,b=S.Integer(1),s
  if r==0:return a
  for j in range(1,r):a,b=b,S.expand((s*b-(n-j+1)*a)/(j+1))
  return b
 F=base-c*d0-f*d1-g*d2
 for ww,coef in C.items():
  term=coef
  for m,co,si in zip(range(3,7),(c,f,g,h),(y,z,w,v)):
   term*=ec(co,si,ww.count(m))
  F+=term
 qf=2*d1+C[4,4]+y*C[3,4,4]+z
 qg=2*d2+C[5,5]
 H=2*d0+C[3,3]+y*C[3,3,3]+z*C[3,3,4]+w*C[3,3,5]+2*y*y+2
 dv=2*(C[6,]+y*C[3,6]+z*C[4,6]+w*C[5,6])+y*y
 C0=(6*base+6*(y*C[3,]+z*C[4,]+w*C[5,])
  +3*y*y*C[3,3]+6*y*z*C[3,4]+6*y*w*C[3,5]
  +3*z*z*C[4,4]+6*z*w*C[4,5]+3*w*w*C[5,5]
  +(y**3+2*y)*C[3,3,3]+3*y*y*z*C[3,3,4]
  +3*y*y*w*C[3,3,5]+3*y*z*z*C[3,4,4]
  +6*y*z*w+(z**3+2*z)+y**4+8*y*y)
 canonical=C0+3*c*c-3*H*c-3*v*c+3*dv*v+3*v*v-3*qf*f-3*qg*g-3*h
 assert S.expand(6*F-canonical)==0
 print("full symbolic count identity: PASS")
 c0,f0,g0,i,j,q,r,AF,AG,AH,AL,Qm=S.symbols('c0 f0 g0 i j q r AF AG AH AL Qm')
 u=q+2*j;cc=c0+2*i
 M=2*i+u+r+2*(AF+AG+AH+AL)
 for eps in (-1,1):
  vv=eps*u
  full=canonical.subs({c:cc,f:f0+2*AF,g:g0+2*AG,h:u+2*AH,v:vv})
  lower=(C0+3*cc*cc-3*H*cc-3*vv*cc+3*dv*vv+3*u*u
        -3*qf*f0-3*qg*g0-3*u-3*Qm*(M-2*i-u-r))
  gap=6*((Qm-qf)*AF+(Qm-qg)*AG+(Qm-1)*AH+Qm*AL)
  assert S.expand(full-lower-gap)==0
 print("all-allocation gap identities, both signs: PASS")
repair6()

def balanced(poly):
 rem=poly;out={}
 for n in range(int(S.degree(rem,x)),-1,-1):
  v=S.expand(rem).coeff(x,n)*factorial(n//2)*factorial((n+1)//2)
  if v:
   out[n]=S.factor(v)
   rem=S.expand(rem-v*P[n//2]*P[(n+1)//2])
 assert rem==0
 return out
words=[w for r in range(1,5)
       for w in combinations_with_replacement(range(3,8),r) if sum(w)<=14]
base=[kernel(d,()) for d in range(8)]
diags=[diagonal(base[d],d) for d in range(8)]
raw=[kernel(7,w) for w in words]
bal=[balanced(p) for p in raw]
D={'words':words}
def upper(ps):
 ds=[S.Poly(p,t,b,k).as_dict() for p in ps]
 mons=set().union(*(v.keys() for v in ds))
 return sum(max(abs(v.get(m,0)) for v in ds)*t**m[0]*b**m[1]*k**m[2]
            for m in mons)
E=[S.Integer(0)]*8;O=[S.Integer(0)]*7
for n0 in range(14):
 bound=sum(k**r/factorial(r)*upper([v.get(n0,0)
              for w,v in zip(words,bal) if len(w)==r]) for r in range(1,5))
 if n0%2:O[n0//2]=S.expand(bound)
 else:E[n0//2]=S.expand(bound+k*upper([diags[d].get(n0//2,0) for d in range(4)]))
O[1]+=k*k/3;O[0]+=k*k*(t+2*b+1)/3;E[0]+=k**3/2
weights=list(map(S.Rational,['64/5','32/5','32/5','6','32/7','5/2','1']))
R=[S.expand(diags[7][j]-E[j]-(weights[j]*O[j]/2 if j<7 else 0)
            -(O[j-1]/(2*weights[j-1]) if j else 0)) for j in range(8)]
for v,H in zip((t,b,k),(249,137,333)):
 for r in R:
  assert all(a>=0 for a in S.Poly(r.subs(v,v+H),t,b,k).coeffs())
print('distance-7 shifted certificates: 24 PASS',flush=True)
assert kernel(3,(3,))==P[3] and kernel(3,(3,3))==1
cs=S.symbols('c f g h n');ys=S.symbols('y z w v u')
c,f,g,h,n=cs;y,z,w,v,u=ys
Cc={word:S.Symbol('C'+''.join(map(str,word))) for word in words}
special={(7,7):1,(3,3,7):x,(3,4,7):1,
         (3,3,3,4):7*x,(3,3,3,5):3,(3,3,4,4):4}
rr=dict(zip(words,raw))
for word,val in special.items():
 assert S.expand(rr[word]-val)==0
 Cc[word]=S.sympify(val)
assert S.expand(rr[3,3,3,3]-(26*b+8*k-5*t+13*x*x-40)/2)==0
B,D3,D2,D1,D0,p3=S.symbols('B D3 D2 D1 D0 p3')
def ec(co,si,r):
 aa,bb=S.Integer(1),si
 if r==0:return aa
 for j in range(1,r):
  aa,bb=bb,S.expand((si*bb-(co-j+1)*aa)/(j+1))
 return bb
F=B-c*D3-f*D2-g*D1-h*D0
for word,cv in Cc.items():
 F+=cv*S.prod(ec(co,si,word.count(m))
              for m,co,si in zip(range(3,8),cs,ys))
F=S.expand(24*(F-(c-1)*y*p3-(c-2)*ec(c,y,2)))
cp=S.Poly(F,*cs,u)
ce=[cp.coeff_monomial(mon) for mon in (c*c,c,f,g,h,n,u,S.Integer(1))]
assert S.expand(F-(ce[0]*c*c+24*c*f-12*x*c*u+ce[1]*c+ce[2]*f+
                  ce[3]*g+ce[4]*h+ce[5]*n+12*u*u+ce[6]*u+ce[7]))==0
V=S.Poly(F.subs({g:0,h:0,n:0}),w,v,u)
HP=S.Poly(V.coeff_monomial(1),c,f)
LP=[S.Poly(V.coeff_monomial(si),c,f) for si in (w,v,u)]
assert HP.coeff_monomial(c*f)==24
se=[HP.coeff_monomial(mon) for mon in (c*c,c,f,S.Integer(1))]
se +=[ll.coeff_monomial(mon) for ll in LP for mon in (c,f,S.Integer(1))]
se +=[S.diff(F,co) for co in (g,h,n)]
K=2*b+k-3;RR=(K*(K+1)-6*b)/2
checks={(3,5,5):2*x,(4,5,5):1,(3,5,6):1,(5,7):P[2],(6,7):x,
        (5,6):x*(P[2]+K),(6,6):x*x+K,(5,5):P[2]**2+K*x*x+RR}
for word,rhs in checks.items():assert S.expand(rr[word]-rhs)==0
G=S.hessian(F,(w,v,u))/24
want=S.Matrix([[P[2]**2+K*x*x+RR+2*x*y+z,x*(P[2]+K)+y,P[2]],
               [x*(P[2]+K)+y,x*x+K,x],[P[2],x,1]])
rs={S.Symbol('C'+''.join(map(str,word))):rv for word,rv in rr.items()}
assert all(S.expand(val)==0 for val in G.subs(rs)-want)
print('canonical allocation and Schur identities: PASS',flush=True)

from sympy.printing.c import C99CodePrinter
from math import lcm
class Printer(C99CodePrinter):
 def _print_Pow(self,ex):
  assert ex.exp.is_Integer and ex.exp>=0
  return '('+'*'.join(['('+self._print(ex.base)+')']*int(ex.exp))+')' if ex.exp else '1'
pr=Printer()
def arr(a):
 return '{'+','.join(arr(v) if isinstance(v,(list,tuple)) else str(v) for v in a)+'}'
fields=[]
def add(p):
 p=S.expand(p)
 if p==0:return -1
 fields.append(p);return len(fields)-1
aa=[add(diags[7][j]) for j in range(8)]
deps=[[add(diags[d].get(j,0)) for j in range(8)] for d in (3,2,1,0)]
kr=[[add(q.get(n,0)) for n in range(15)] for q in bal]
den=lcm(*(int(v.q) for p in fields for v in S.Poly(p,t,b,k).coeffs()))
def emit(name,args,expr,out='out'):
 rr,ee=S.cse(expr,symbols=S.numbered_symbols('a'))
 body=['I '+str(v)+'='+pr.doprint(q)+';' for v,q in rr]
 body +=[out+'['+str(i)+']='+pr.doprint(q)+';' for i,q in enumerate(ee)]
 return 'static void '+name+'('+args+'){\n'+'\n'.join(body)+'\n}\n'
code=f'\nstatic const int SCALE={den},NF={len(fields)},NW={len(words)};\n'
code+=f'static const int aa[8]={arr(aa)},dep[4][8]={arr(deps)};\n'
code+=f'static const int ker[NW][15]={arr(kr)};\n'
code+=f'static const int wr[NW]={arr([len(w) for w in words])};\n'
code+=f'static const int wn[NW][4]={arr([list(w)+[0]*(4-len(w)) for w in words])};\n'
code+=f'static const int ids[NW]={arr([int("".join(map(str,w))) for w in words])};\n'
code+=emit('eval','I t,I b,I k,I*v',[S.expand(den*p) for p in fields],'v')
subs={B:S.Symbol('base'),D3:S.Symbol('dd[0]'),D2:S.Symbol('dd[1]'),
      D1:S.Symbol('dd[2]'),D0:S.Symbol('dd[3]')}
subs.update({S.Symbol('C'+''.join(map(str,ww))):S.Symbol('cc['+str(j)+']')
             for j,ww in enumerate(words)})
canon=emit('canonical_coeff','I x,I y,I z,I w,I v,I base,I*cc,I*dd,I p3,I*out',
           [S.expand(a.subs(subs)) for a in ce])
schurcode=emit('schur_coeff','I x,I y,I z,I base,I*cc,I*dd,I p3,I*out',
               [S.expand(a.subs(subs)) for a in se])

import os,ctypes,subprocess,glob,resource
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
def mem(name,data=b''):
 fd=os.memfd_create(name,0)
 if data:os.write(fd,data)
 return fd
def build(src):
 st=time.monotonic()
 sf=mem('fm80.cpp',src.encode())
 p=subprocess.run(['g++','-std=c++17','-O3','-ftrapv','-fPIC','-pthread',
  '-S','-x','c++','-o','-',f'/proc/self/fd/{sf}'],
  pass_fds=(sf,),capture_output=True)
 if p.returncode:raise RuntimeError(p.stderr.decode())
 af=mem('fm80.s',p.stdout);of=mem('fm80.o')
 p=subprocess.run(['as','--64','-o',f'/proc/self/fd/{of}',f'/proc/self/fd/{af}'],
                  pass_fds=(of,af),capture_output=True)
 if p.returncode:raise RuntimeError(p.stderr.decode())
 so=mem('fm80.so')
 gcc=sorted(glob.glob('/usr/lib/gcc/x86_64-linux-gnu/*/libgcc.a'))[-1]
 gd=os.path.dirname(gcc)
 p=subprocess.run(['ld.gold','-shared','-o',f'/proc/self/fd/{so}',
  gd+'/crtbeginS.o',f'/proc/self/fd/{of}',gcc,
  '/usr/lib/x86_64-linux-gnu/libstdc++.so.6','-lpthread','-lc','-lm',
  gd+'/crtendS.o'],pass_fds=(of,so),capture_output=True)
 if p.returncode:raise RuntimeError(p.stderr.decode())
 lib=ctypes.CDLL(f'/proc/self/fd/{so}')
 print('checked memory compilation seconds',round(time.monotonic()-st,2),flush=True)
 return lib

HEADER=r'''
#include <boost/multiprecision/cpp_int.hpp>
#include <pthread.h>
#include <atomic>
#include <stdio.h>
#include <stdint.h>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <array>
#include <iostream>
#include <chrono>
using I=__int128;
using boost::multiprecision::cpp_int;
static I ab(I x){return x<0?-x:x;}
static I binom(int n,int r){if(r<0||r>n)return 0;I v=1;for(int j=1;j<=r;j++)v=v*(n-j+1)/j;return v;}
static cpp_int big(I x){bool neg=x<0;unsigned __int128 u=neg?-x:x;cpp_int a=(uint64_t)(u>>64);a<<=64;a+=(uint64_t)u;return neg?-a:a;}
static I fd(I a,I b){I q=a/b,r=a%b;if(r<0)q--;return q;}
static int clip(I a,int J){return a<0?0:a>J?J:(int)a;}
static I minquad(I a,I b,I c,int J){
 I ans=std::min(c,a*J*J+b*J+c);
 if(a>0){int j=clip(fd(a-b,2*a),J);ans=std::min(ans,a*j*j+b*j+c);}
 return ans;
}
'''

GRID=r'''
static I field(I*v,int j){return j<0?0:v[j];}
static bool covers(int t,int b,int k,I*v){
 I D[8]={},O[7]={},E[8]={};
 for(int n=0;n<14;n++){
  I z=0;
  for(int r=1;r<=4;r++){
   I mx=0;
   for(int a=0;a<NW;a++)if(wr[a]==r)mx=std::max(mx,ab(field(v,ker[a][n])));
   z+=binom(k,r)*mx;
  }
  if(n%2)O[n/2]=z;else E[n/2]=z;
 }
 if(k>=2){O[1]+=I(k)*(k-1)*(SCALE/3);O[0]+=I(k)*(k-1)*ab(t-2*b-1)*(SCALE/3);}
 if(k>=3)E[0]+=binom(k,2)*(k-2)*SCALE;
 for(int j=0;j<8;j++){
  I mx=0;for(int h=0;h<4;h++)mx=std::max(mx,ab(field(v,dep[h][j])));
  E[j]+=k*mx;D[j]=field(v,aa[j])-E[j];if(D[j]<=0)return false;
 }
 cpp_int prev=1,cur=big(2*D[7]);
 for(int j=6;j>=0;j--){
  cpp_int off=big(O[j]),nxt=big(2*D[j])*cur-off*off*prev;
  if(nxt<=0)return false;prev=cur;cur=nxt;
 }
 return true;
}
static void profile(int t,int b,int k,int x,I*f,I&base,I*cc,I*dd,I&p3){
 I P[8]={1,x};
 for(int j=1;j<7;j++){
  I v=I(x)*P[j]-(t-2*b-j+1)*P[j-1];
  if(v%(j+1))std::exit(3);P[j+1]=v/(j+1);
 }
 p3=P[3];base=0;
 for(int j=0;j<8;j++)base+=field(f,aa[j])*P[j]*P[j];
 if(base%SCALE)std::exit(3);base/=SCALE;
 for(int h=0;h<4;h++){
  I z=0;for(int j=0;j<8;j++)z+=field(f,dep[h][j])*P[j]*P[j];
  if(z%SCALE)std::exit(3);dd[h]=z/SCALE;
 }
 for(int a=0;a<NW;a++){
  I z=0;for(int n=0;n<15;n++)z+=field(f,ker[a][n])*P[n/2]*P[(n+1)/2];
  if(z%SCALE)std::exit(3);cc[a]=z/SCALE;
 }
}
static I get(I*cc,int key){for(int i=0;i<NW;i++)if(ids[i]==key)return cc[i];return 0;}
static I quartic(int k,int x,I p3,I A){
 I ans=0;
 for(int c=0;c<=k;c++)for(int y=c%2;y<=c;y+=2){
  I e2=(I(y)*y-c)/2,e3=(I(y)*y*y-(3*c-2)*I(y))/6;
  I e4=(I(y)*y*y*y+(8-6*c)*I(y)*y+3*I(c)*c-6*c)/24;
  int K=k-c;I ae3=ab(e3),coef=3*ae3-2*e2;
  I C=A*e4-(c-2)*e2-ab((c-1)*I(y)*p3)-3*ae3*K;
  I qA=2*e2,qB=-7*ab(I(x)*e3);
  if(coef>=0)ans=std::min(ans,minquad(qA,qB+coef,C,K));
  else for(int q=0;q<=1&&q<=K;q++){
   int f=K-(K-q)%2,J=(K-q)/2;
   ans=std::min(ans,minquad(4*qA,4*qA*q+2*qB,qA*q*q+qB*q+C+coef*f,J));
  }
 }
 return ans;
}
static int scalar(int k,int x,I base,I*cc,I*dd,I p3){
 I depmax=0,single=0,dpos=0,dneg=0,off=0,triple=0,four=0;
 for(int i=0;i<4;i++)depmax=std::max(depmax,dd[i]);
 for(int a=0;a<NW;a++){
  I z=cc[a];int r=wr[a];
  if(r==1)single=std::max(single,ab(z));
  else if(r==2){
   if(wn[a][0]==wn[a][1]){dpos=std::max(dpos,z);dneg=std::max(dneg,-z);}
   else off=std::max(off,ab(z));
  }else if(r==3)triple=std::max(triple,ab(z));
  else four=std::max(four,ab(z));
 }
 I common=2*base-2*k*depmax-2*k*single-k*dpos-I(k)*k*std::max(dneg,off)-2*binom(k,3)*triple;
 I cheap=common-2*binom(k,4)*four;
 if(k>=2)cheap-=2*I(k)*(k-1)*ab(p3);
 if(k>=3)cheap-=2*binom(k,2)*(k-2);
 if(cheap>=0)return 1;
 I low=quartic(k,x,p3,get(cc,3333));
 return common+2*low>=0?2:0;
}
static bool schur(int t,int b,int k,int x,I base,I*cc,I*dd,I p3){
 I K=2*b+k-3;if(K<=0)return false;
 I R=(K*(K+1)-6*b)/2;
 if(std::min(K*(R-k),K*R-I(k)*k)<=0)return false;
 I P2=(I(x)*x-t+2*b)/2;
 for(int y=-k;y<=k;y++)for(int z=-(k-abs(y));z<=k-abs(y);z++){
  I f[16];schur_coeff(x,y,z,base,cc,dd,p3,f);
  I pen=std::min(I(0),std::min(f[13],std::min(f[14],f[15])));
  I A=K*x*x+R+2*I(x)*y+z,B=K*x+y,det=K*(R+z)-I(y)*y;
  if(det<=0)std::exit(10);
  int c0=abs(y),f0=abs(z),J=(k-c0-f0)/2;
  auto val=[&](int i,int branch)->I{
   I c=c0+2*i,ff=branch?k-c:f0;
   I H=f[0]*c*c+24*c*ff+f[1]*c+f[2]*ff+f[3]+pen*(k-c-ff);
   I lw=f[4]*c+f[5]*ff+f[6],lv=f[7]*c+f[8]*ff+f[9],lu=f[10]*c+f[11]*ff+f[12];
   I aw=lw-P2*lu,av=lv-x*lu;
   return 48*det*H-det*lu*lu-K*aw*aw+2*B*aw*av-A*av*av;
  };
  for(int branch=0;branch<2;branch++){
   I v0=val(0,branch),v1=val(1,branch),v2=val(2,branch);
   I aa=v2-2*v1+v0;if(aa%2)std::exit(10);aa/=2;
   I bb=v1-v0-aa;
   if(minquad(aa,bb,v0,J)<0)return false;
  }
 }
 return true;
}
struct Result{long long matrix=0,failed=0,cheap=0,quart=0,sch=0;std::vector<std::array<int,4>> rows;};
static Result result[4];
static std::atomic<int> nextt(0),done(0);
static void* work(void*arg){
 int id=(intptr_t)arg;auto&r=result[id];
 for(;;){
  int t=nextt.fetch_add(1);if(t>=249)break;
  for(int b=0;b<137;b++)for(int k=0;k<333;k++){
   if(t+2*b+8*k<14)continue;
   I f[NF];eval(t,b,k,f);
   if(covers(t,b,k,f)){r.matrix++;continue;}r.failed++;
   for(int x=t%2;x<=t;x+=2){
    I base,cc[NW],dd[4],p3;profile(t,b,k,x,f,base,cc,dd,p3);
    int s=scalar(k,x,base,cc,dd,p3);
    if(s==1)r.cheap++;else if(s==2)r.quart++;
    else if(schur(t,b,k,x,base,cc,dd,p3))r.sch++;
    else r.rows.push_back({t,b,k,x});
   }
  }
  int nd=done.fetch_add(1)+1;
  if(nd%16==0||nd==249){printf("GRID progress completed_t=%d/249\n",nd);fflush(stdout);}
 }
 return nullptr;
}
static std::vector<std::array<int,4>> final_rows;
extern "C" long long rowcount(){return final_rows.size();}
extern "C" const int* rowdata(){return reinterpret_cast<const int*>(final_rows.data());}
extern "C" void run(){
 pthread_t th[4];for(int i=0;i<4;i++)pthread_create(&th[i],nullptr,work,(void*)(intptr_t)i);
 for(int i=0;i<4;i++)pthread_join(th[i],nullptr);
 long long ma=0,fa=0,ch=0,qu=0,sh=0;std::vector<std::array<int,4>> rows;
 for(auto&r:result){ma+=r.matrix;fa+=r.failed;ch+=r.cheap;qu+=r.quart;sh+=r.sch;rows.insert(rows.end(),r.rows.begin(),r.rows.end());}
 std::sort(rows.begin(),rows.end());
 printf("GRID matrix=%lld failed=%lld cheap=%lld quartic=%lld schur=%lld remaining=%zu\n",ma,fa,ch,qu,sh,rows.size());fflush(stdout);
 if(ma!=11099400||fa!=260161||ch!=2112211||qu!=2934||sh!=39035||rows.size()!=580)std::exit(23);
 final_rows=rows;int mt=0,mb=0,mk=0,mx=0,mn=0;
 unsigned long long total=0;
 for(auto&a:rows){
  mt=std::max(mt,a[0]);mb=std::max(mb,a[1]);mk=std::max(mk,a[2]);mx=std::max(mx,a[3]);mn=std::max(mn,a[0]+2*a[1]+a[2]);
  int k=a[2];
  for(int s=0;s<=k;s++){
   I ns=s==0?1:0;
   for(int r=1;r<=4&&r<=s;r++)ns+=binom(4,r)*binom(s-1,r-1)*(1<<r);
   total+=(unsigned long long)(ns*10*(1+(s<k)));
  }
 }
 if(mn!=33)std::exit(23);
 printf("RESIDUAL max_N=%d\n",mn);
 printf("RESIDUAL max_t=%d max_b=%d max_k=%d max_x=%d vertex_quadratics=%llu\n",mt,mb,mk,mx,total);fflush(stdout);
}
'''

FINISH=r'''
static I ec(int c,int y,int r){
 if(r==0)return 1;I a=1,b=y;
 for(int j=1;j<r;j++){I num=I(y)*b-(c-j+1)*a;if(num%(j+1))std::exit(3);I z=num/(j+1);a=b;b=z;}
 return b;
}
static I fullvalue(int k,I base,I*cc,I*dd,I p3,const int*c,const int*s){
 I v=base;for(int j=0;j<4;j++)v-=c[j]*dd[j];
 for(int a=0;a<NW;a++){
  int mult[5]={};for(int j=0;j<wr[a];j++)mult[wn[a][j]-3]++;
  I z=cc[a];for(int j=0;j<5;j++)z*=ec(c[j],s[j],mult[j]);v+=z;
 }
 v-=(c[0]-1)*I(s[0])*p3+(c[0]-2)*ec(c[0],s[0],2);
 return v;
}
struct FStats{unsigned long long brute=0,quads=0,candidates=0;I least=I(1)<<120;};
static FStats fs[4];
static I tri_min(I A,I B,I D,I E,I F,int J,FStats&st){
 if(B%48)std::exit(5);I best=I(1)<<120;
 for(int r=0;r<=1&&r<=J;r++){
  I a=fd(48-B*r-E,96),b=-B/48;
  for(int reg=0;reg<3;reg++){
   I lo=0,hi=(J-r)/2,aa,bb;
   auto le=[&](I c,I d){
    if(c>0)hi=std::min(hi,fd(d,c));
    else if(c<0)lo=std::max(lo,-fd(d,-c));
    else if(d<0)hi=-1;
   };
   if(reg==0){le(b,-a);aa=0;bb=0;}
   else if(reg==1){le(-(b+2),a-J+r);aa=J-r;bb=-2;}
   else{le(-b,a);le(b+2,J-r-a);aa=a;bb=b;}
   if(lo>hi)continue;
   auto val=[&](I hh)->I{
    I i=r+2*hh,j=aa+bb*hh;
    return A*i*i+B*i*j+48*j*j+D*i+E*j+F;
   };
   I v0=val(0),v1=val(1),v2=val(2);
   I qa=v2-2*v1+v0;if(qa%2)std::exit(5);qa/=2;
   I qb=v1-v0-qa;
   auto test=[&](I hh){
    I i=r+2*hh,j=aa+bb*hh;
    if(i<0||j<0||i+j>J)std::exit(5);
    best=std::min(best,val(hh));st.candidates++;
   };
   test(lo);test(hi);
   if(qa>0){I hh=fd(qa-qb,2*qa);hh=std::max(lo,std::min(hi,hh));test(hh);}
  }
 }
 return best;
}
static void finite(int t,int b,int k,int x,FStats&st){
 I fields[NF],base,cc[NW],dd[4],p3;
 eval(t,b,k,fields);profile(t,b,k,x,fields,base,cc,dd,p3);
 if(t+2*b+3*k<14){
  int c[5]={},s[5]={};
  auto rec=[&](auto&&self,int j,int left,int wt)->void{
   if(j==5){
    if(t+2*b+wt+8*left<14)return;
    I val=24*fullvalue(k,base,cc,dd,p3,c,s);
    if(val<0){printf("NEG_SMALL %d %d %d %d\n",t,b,k,x);fflush(stdout);std::exit(8);}
    st.least=std::min(st.least,val);st.brute++;return;
   }
   for(int co=0;co<=left;co++){c[j]=co;for(int si=-co;si<=co;si+=2){
    s[j]=si;self(self,j+1,left-co,wt+(j+3)*co);
   }}
  };rec(rec,0,k,0);return;
 }
 for(int y=-k;y<=k;y++)
 for(int z=-(k-abs(y));z<=k-abs(y);z++)
 for(int w=-(k-abs(y)-abs(z));w<=k-abs(y)-abs(z);w++)
 for(int v=-(k-abs(y)-abs(z)-abs(w));v<=k-abs(y)-abs(z)-abs(w);v++){
  int c0=abs(y),f0=abs(z),g0=abs(w),h0=abs(v),M=k-c0-f0-g0-h0;
  I co[8];canonical_coeff(x,y,z,w,v,base,cc,dd,p3,co);
  for(int q=0;q<=1&&q<=M;q++)for(int eps:{-1,1}){
   int J=(M-q)/2;
   for(int which=0;which<5;which++){
    auto val=[&](int i,int j)->I{
     I c=c0+2*i,f=f0,g=g0,h=h0,n=q+2*j,u=eps*n,extra=2*(J-i-j);
     if(which==1)f+=extra;
     if(which==2)g+=extra;
     if(which==3)h+=extra;
     if(which==4)n+=extra;
     return co[0]*c*c+24*c*f-12*I(x)*c*u+co[1]*c+co[2]*f+co[3]*g+co[4]*h+co[5]*n+12*u*u+co[6]*u+co[7];
    };
    I F=val(0,0),V10=val(1,0),V01=val(0,1);
    I A=val(2,0)-2*V10+F,C=val(0,2)-2*V01+F;
    if(A%2||C!=96)std::exit(5);A/=2;
    I B=val(1,1)-V10-V01+F,D=V10-F-A,E=V01-F-48;
    if(B!=-48*eps*x-(which==1?96:0))std::exit(5);
    I low=tri_min(A,B,D,E,F,J,st);st.quads++;
    if(low<0){
     printf("NEG_Q t=%d b=%d k=%d x=%d signs=%d,%d,%d,%d q=%d eps=%d which=%d\n",t,b,k,x,y,z,w,v,q,eps,which);
     fflush(stdout);std::exit(8);
    }
    st.least=std::min(st.least,low);
   }
  }
 }
}
static const int* inputrows;
static long long ninput;
static std::atomic<long long> nextrow(0),donerow(0);
static void* fwork(void*arg){
 int id=(intptr_t)arg;
 for(;;){
  long long j=nextrow.fetch_add(1);if(j>=ninput)break;
  const int*r=inputrows+4*j;finite(r[0],r[1],r[2],r[3],fs[id]);
  long long nd=donerow.fetch_add(1)+1;
  if(nd%40==0||nd==ninput){printf("FINITE progress rows=%lld/%lld\n",nd,ninput);fflush(stdout);}
 }
 return nullptr;
}
extern "C" void finish(const int*rows,long long n){
 inputrows=rows;ninput=n;
 pthread_t th[4];for(int i=0;i<4;i++)pthread_create(&th[i],nullptr,fwork,(void*)(intptr_t)i);
 for(int i=0;i<4;i++)pthread_join(th[i],nullptr);
 FStats out;
 for(auto&s:fs){out.brute+=s.brute;out.quads+=s.quads;out.candidates+=s.candidates;out.least=std::min(out.least,s.least);}
 if(out.brute!=8157||out.quads!=11135300||out.candidates!=45853270||out.least!=0)std::exit(24);
 printf("FINITE brute=%llu quadratics=%llu candidates=%llu minimum24F=%lld\n",out.brute,out.quads,out.candidates,(long long)out.least);
 fflush(stdout);
}
'''

BRIDGE=r'''
extern "C" const char* val7(int t,int b,int k,int x,const int*c,const int*s){
 static std::string text;
 I fields[NF],base,cc[NW],dd[4],p3;
 eval(t,b,k,fields);profile(t,b,k,x,fields,base,cc,dd,p3);
 text=big(2*fullvalue(k,base,cc,dd,p3,c,s)).convert_to<std::string>();
 return text.c_str();
}
'''

COUNT=r'''
#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
using boost::multiprecision::cpp_int;
using I=__int128;
static cpp_int bn(int n,int r){
 if(r<0||n<r)return 0;cpp_int v=1;
 for(int j=1;j<=r;j++){v*=n-j+1;v/=j;}return v;
}
static cpp_int seed(int n,bool plus){
 cpp_int ans=0;
 for(int h=0;h<=n/2;h++)
  ans+=(plus?cpp_int(1):cpp_int(h+1))*bn(n-2*h+(plus?6:11),plus?6:11);
 return ans;
}
static std::vector<cpp_int> seq(int max,bool plus){
 std::vector<cpp_int> out(max+1);int deg=plus?7:13;
 for(int p=0;p<2;p++){
  std::vector<cpp_int> v(deg+1),dif;
  for(int j=0;j<=deg;j++)v[j]=seed(2*j+p,plus);
  for(int d=0;d<=deg;d++){
   dif.push_back(v[0]);
   for(int j=0;j<deg-d;j++)v[j]=v[j+1]-v[j];
  }
  for(int n=p;n<=max;n+=2){
   out[n]=dif[0];
   for(int j=0;j<deg;j++)dif[j]+=dif[j+1];
  }
 }
 for(int n=0;n<=std::min(max,100);n++)if(out[n]!=seed(n,plus))std::exit(2);
 return out;
}
extern "C" void count(){
 const int C=1687040,A=24320;const long long den=14745600;
 std::vector<int> arg(C,-1);int q=0,mx=0,last=-1;
 for(int n=0;n<C;n++){
  if(n>=A){
   I num=I(n)*(3*n-72960)*(3*n-72960);
   while(I(q)*q*den<=num)q++;
  }
  if(q<=n){arg[n]=n-q;mx=std::max(mx,n-q);last=n;}
 }
 auto G=seq(mx,false),F=seq(mx,true);
 cpp_int all=0,plus=0;
 for(int n=0;n<C;n++)if(arg[n]>=0){
  if(n<A)all+=G[n]+(n?G[n-1]:cpp_int(0));
  else all+=2*G[arg[n]];
  plus+=F[arg[n]];
 }
 int lab[14]={1,1,2,3,3,4,4,5,5,6,6,7,7,8};
 bool neg[14]={false,true,false,false,true,false,true,false,true,false,true,false,true,false};
 long long bad=0,small=0;
 auto rec=[&](auto&&self,int i,int W,int mm,bool mn)->void{
  if(i==14){small++;if(mn&&W-14<std::max(5,mm))bad++;return;}
  for(int c=0;W+c*lab[i]<=21;c++)
   self(self,i+1,W+c*lab[i],c?std::max(mm,lab[i]):mm,mn||(c&&neg[i]));
 };
 rec(rec,0,0,0,false);
 cpp_int two=0;
 for(int n=15;n<A;n++){
  long long w=n-1,h=w/2;
  two+=(h+1)*(3*(w-h)+1);
 }
 std::cout<<"max_argument "<<mx<<" last_N "<<last<<"\n";
 std::cout<<"Gaussian_failed "<<all<<"\nall_plus "<<plus
          <<"\nsmall_inadmissible "<<bad<<" small_enumerated "<<small
          <<"\ntwo_core "<<two<<"\nEXACT_RESIDUAL "<<all-plus-bad-two<<"\n";
 if(all-plus-bad-two!=cpp_int("786505985094237249570124700669525047622199528819431922103292935154"))std::exit(25);
 std::cout.flush();
}
'''

CHECKS=r'''
extern "C" void checks(){
 unsigned long long seed=8007;
 auto rnd=[&](int m)->int{seed=seed*6364136223846793005ULL+1;return (seed>>32)%m;};
 for(int it=0;it<2500;it++){
  I A=rnd(30000)-10000,B=48*(rnd(61)-30),D=rnd(200001)-100000;
  I E=rnd(200001)-100000,F=rnd(200001)-100000;int J=rnd(51);
  FStats st;I z=I(1)<<120;
  for(int i=0;i<=J;i++)for(int j=0;j<=J-i;j++)
   z=std::min(z,A*i*i+B*i*j+48*j*j+D*i+E*j+F);
  if(z!=tri_min(A,B,D,E,F,J,st))std::exit(21);
 }
 for(int it=0;it<160;it++){
  int k=rnd(11),x=rnd(21)-10;I p3=rnd(81)-40,A=rnd(131)-50,z0=0;
  for(int c=0;c<=k;c++)for(int y=-c;y<=c;y+=2)
  for(int f=0;f<=k-c;f++)for(int z=-f;z<=f;z+=2){
   I v=A*ec(c,y,4)+7*x*ec(c,y,3)*z-3*ab(ec(c,y,3))*(k-c-f)
       +4*ec(c,y,2)*ec(f,z,2)-(c-1)*I(y)*p3-(c-2)*ec(c,y,2);
   z0=std::min(z0,v);
  }
  if(z0!=quartic(k,x,p3,A))std::exit(22);
 }
 printf("quadratic minimizer: 2500; quartic minimizer: 160 exact checks PASS\n");fflush(stdout);
}
'''

lib=build(HEADER+code+schurcode+GRID+canon+FINISH+BRIDGE+COUNT+CHECKS)
lib.checks()

from random import Random
@lru_cache(None)
def inv(ns):
 if not ns:return 1
 if sum(ns)%2 or 2*max(ns)>sum(ns):return 0
 return moment(0,ns)
def even_word(word):
 ans=0
 for mask in range(1<<len(word)):
  left=[];right=[];sg=1
  for j,(n,e) in enumerate(word):
   if mask>>j&1:left.append(n);sg*=e
   else:right.append(n)
  ans+=sg*inv(tuple(sorted(left)))*inv(tuple(sorted(right)))
 return ans
def gf7(rest):
 row={(0,0):1}
 for nn,ee in rest:
  block={(2*q,0):1 for q in range(min(nn,7)+1)}
  if nn<=14:
   for q in range(nn//2+1):
    pos=(nn,nn-2*q)
    block[pos]=block.get(pos,0)+ee*(-1)**q*comb(nn-q,q)
  new={}
  for (i,j),a in row.items():
   for (ii,jj),bb in block.items():
    if i+ii<=14:new[i+ii,j+jj]=new.get((i+ii,j+jj),0)+a*bb
  row={m:a for m,a in new.items() if a}
 return 2*sum(a*comb(j,j//2)//(j//2+1)*((i==14)-(i==12))
              for (i,j),a in row.items() if j%2==0)
def formula(rest):
 tt=bb=kk=xx=0;cs=[0]*5;ss=[0]*5
 for nn,ee in rest:
  if nn==1:tt+=1;xx+=ee
  elif nn==2:
   if ee==1:bb+=1
   else:tt+=2
  else:
   kk+=1
   if nn<=7:cs[nn-3]+=1;ss[nn-3]+=ee
 return int(lib.val7(tt,bb,kk,xx,(ctypes.c_int*5)(*cs),(ctypes.c_int*5)(*ss)))
lib.val7.argtypes=[ctypes.c_int]*4+[ctypes.POINTER(ctypes.c_int)]*2
lib.val7.restype=ctypes.c_char_p
rng=Random(8007)
for it in range(100):
 while True:
  rest=[(rng.choice([1,1,2,2,3,3,3,4,5,6,7,8,19,53]),rng.choice([-1,1]))
        for j in range(rng.randrange(5,11))]
  pp=sum(n for n,e in rest)-14
  if pp>=max(n for n,e in rest):break
 ep=(-1)**sum(e<0 for n,e in rest)
 assert formula(rest)==gf7(rest)==even_word(rest+[(pp,ep)])
for it in range(100):
 rest=[(rng.randrange(1,10),rng.choice([-1,1])) for j in range(rng.randrange(5,45))]
 assert formula(rest)==gf7(rest)
print('distance-7 direct character bridges: 100; generating-function bridges: 200 PASS',flush=True)

lib.run()
lib.rowdata.restype=ctypes.POINTER(ctypes.c_int)
lib.rowcount.restype=ctypes.c_longlong
lib.finish.argtypes=[ctypes.POINTER(ctypes.c_int),ctypes.c_longlong]
assert lib.rowcount()==580
lib.finish(lib.rowdata(),lib.rowcount())
lib.count()
print('ALL CHECKS PASS; seconds',round(time.monotonic()-started,2),flush=True)