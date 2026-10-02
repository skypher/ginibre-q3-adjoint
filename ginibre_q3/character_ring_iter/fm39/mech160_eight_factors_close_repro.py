import argparse,itertools,os,re,shlex,subprocess
from fractions import Fraction as F
from functools import lru_cache
ap=argparse.ArgumentParser(description="FM-MECH160 exact eight-factor certificate; memory only.")
ap.add_argument("--census-log")
args=ap.parse_args()
P={"D":(1,2,3,4,3,2,1),"K6":(1,2,3,4,5,5,4),
   "G":(1,2,3,3,2,1),"F":(1,2,3,3,3,2,1),"J4":(1,2,2,2,2,1)}
checks=0
for g,h,k in itertools.product(range(9),repeat=3):
    ns=tuple(sorted((g+h,g+k,h+k)))
    if ns[0]==0:continue
    b=[max(0,min(s,g)-max(-s,s-h-k,-h,-k)+1) for s in range(7)]
    names=[]
    if ns[0]>=3 and ns!=(3,3,4):names.append("D")
    if ns[0]>=6:names.append("K6")
    if ns[0]==2 and ns[-1]>=5:names.append("G")
    if ns[0]==2 and ns[-1]>=6 and all(n%2==0 for n in ns):names.append("F")
    if ns[0]==1 and ns[-1]>=5:names.append("J4")
    assert all(x>=y for x,y in zip(b,(1,2,1)))
    if ns[0]>=2:assert all(x>=y for x,y in zip(b,(1,2,2,1)))
    for name in names:
        assert all(x>=y for x,y in zip(b,P[name]));checks+=1
print("triple profiles",checks,"PASS")
O=tuple(map(F,(1,F(8,3),F(11,3),F(11,3),3,2,1)))
def band(I,J):
    return [sum(abs(a-b)<=2*s<=a+b for a in I for b in J) for s in range(7)]
for d in range(1,7):
    for par in (0,1):
        if d<=3:
            t=4-d;u=2*t+par
            pairs=[(range(u-2*t,u+2*d,2),range(u,u+2*(d+t),2))]
        else:
            u=2 if par==0 else 1
            pairs=[(range(u,u+2*d,2),range(u,u+2*(d+1),2))]
            u=2 if par==0 else 3
            pairs.append((range(u,u+2*d,2),range(u-2,u+2*d,2)))
        for I,J in pairs:assert all(x>=d*y for x,y in zip(band(I,J),O))
assert all(F(2*s+1)-F(3*s*s-s,14)>=v for s,v in enumerate(O))
def row(ns):
    a={0:1}
    for n in ns:
        q={}
        for j,v in a.items():
            for k in range(abs(n-j),n+j+1,2):q[k]=q.get(k,0)+v
        a=q
    return [a.get(2*s,0) for s in range(7)]
def norm(I,J):
    K=set(I)&set(J)
    if not K or 0 in K:return None
    return tuple(F(x,len(K)) for x in band(I,J))
@lru_cache(None)
def profiles(sm,L):
    if not sm:return (O,)
    H=L+1;out=set()
    if len(sm)==1:
        l=sm[0]+1
        for off in range(1-H,l):
            for par in (0,1):
                u=max(L-sm[0],-2*off,0);u+=(par-u)%2
                if u==0 and off==0:u+=2
                z=norm(range(u,u+2*l,2),range(u+2*off,u+2*(off+H),2))
                if z:out.add(z)
    elif len(sm)==2:
        a,b=sm;I=range(abs(a-b),a+b+1,2)
        for v in range((a+b)%2,a+b+1,2):
            z=norm(I,range(v,v+2*H,2))
            if z:out.add(z)
    else:
        for n in range(L,sum(sm)+1):
            b=row(sm+(n,))
            if b[0]:out.add(tuple(F(x,b[0]) for x in b))
    return tuple(out)
for mode,L,alphabet,first,expected,count in [
    ("fund",5,range(1,5),1,F(122,3),1916),
    ("even",6,(2,4),2,F(176,3),946)]:
    best=None;trials=0
    for k in range(1,4):
        for sm in itertools.combinations_with_replacement(alphabet,k):
            if sm[0]!=first:continue
            for mask in range(1<<k):
                a=tuple(sm[i] for i in range(k) if mask>>i&1)
                b=tuple(sm[i] for i in range(k) if not mask>>i&1)
                for x,y in itertools.product(profiles(a,L),profiles(b,L)):
                    v=sum(z*t for z,t in zip(x,y));trials+=1
                    best=v if best is None else min(best,v)
    assert best==expected and trials==count
    print("quartets",mode,trials,"minimum",best,"PASS")
for o in (0,4,6):
    odd=(1<<o)-1
    for sg in range(256):
        if sg.bit_count()%2:continue
        counts=[]
        for size in (3,4):
            c=0
            for S in itertools.combinations(range(8),size):
                m=sum(1<<i for i in S)
                if size==4 and not m&1:continue
                c+=(m&odd).bit_count()%2==0 and (m&sg).bit_count()%2==1
            counts.append(c)
        c3,c4=counts
        for weak in (0,1):
            if weak>c3:continue
            if o:
                assert F(weak,10)+F(c3-weak,33)+F(c4,40)<=1
            else:
                assert F(weak,20)+F(c3-weak,63)+F(c4,56)<=1
print("sign budgets PASS")

CPP=r'''
#include <bits/stdc++.h>
using namespace std;using Z=long long;
struct Test{vector<int>v;int c;};
void tiles(int id){
 vector<int>F={6,7,8,9,10,11},E={7,9,11,13,15,17,19};
 vector<int>L=id<3?F:E,H;vector<Test>ts;int pm;
 if(id==0){H={2,3};pm=1;ts={{{1,2,3,4,3,2,1},33}};}
 if(id==1){H={2,3};pm=5;ts={{{1,2,2,2,2,1},33},{{1,2,3,3,2,1},33}};}
 if(id==2){H=F;pm=5;ts={{{1,2,1},10}};}
 if(id==3){H={3,5,7};pm=6;ts={{{1,2,3,3,3,2,1},63},{{1,2,3,4,3,2,1},63}};}
 if(id==4){H={3,5,7};pm=2;ts={{{1,2,3,4,5,5,4},63}};}
 if(id==5){H=E;pm=6;ts={{{1,2,2,1},20}};}
 Z count=0;
 for(int l:L)for(int h:H)
 for(int x=-l-h+2;x<=41;x++)
 for(int y=-l+1;y<=41;y++)
 for(int z=-h+1;z<=41;z++){
  int u=x+y,v=x+z,p=y+z;
  if(u<0||v<0||p<pm||(id>=3&&((u|v|p)&1)))continue;
  Z b[7]={};
  for(int i=0;i<l;i++)for(int j=0;j<h;j++)for(int s=0;s<7;s++)
   b[s]+=max(0,min(s,x+i+j)-max({-s,s-p,-y-i+j,-z+i-j})+1);
  if(!b[0])continue;count++;
  for(auto t:ts){
   Z val=0;for(int s=0;s<(int)t.v.size();s++)val+=t.v[s]*b[s];
   assert(val>=t.c*b[0]);
   // l,h <=19: b[s]<=13*19^2; every scalar is <10^6.
   assert(val<1000000);
  }
 }
 cout<<"tiles "<<id<<" "<<count<<" PASS"<<endl;
}
int a[8],pc[256],sw[256],sp[256];Z ch[111][7];
Z words=0,signs=0,minimum=LLONG_MAX;array<int,8>witness;
Z inv(int m){
 int k=pc[m],s=sw[m];if(!k)return 1;if(s%2)return 0;if(k==1)return 0;
 if(k==2){int x=__builtin_ctz((unsigned)m),y=__builtin_ctz((unsigned)(m&(m-1)));return a[x]==a[y];}
 if(k==3){int ma=0;for(int i=0;i<8;i++)if(m>>i&1)ma=max(ma,a[i]);return 2*ma<=s;}
 if(k==4){
  int b[4],j=0;for(int i=0;i<8;i++)if(m>>i&1)b[j++]=a[i];
  return max(0,(min(b[0]+b[1],b[2]+b[3])-max(abs(b[0]-b[1]),abs(b[2]-b[3])))/2+1);
 }
 Z out=0;
 for(int t=m;;t=(t-1)&m){
  int top=s/2-sp[t]+k-2;
  if(top>=k-2){assert(top<=110);out+=(pc[t]%2?-1:1)*ch[top][k-2];}
  if(!t)break;
 }
 assert(out>=0);return out;
}
void word(){
 int odd=0;for(int x:a)odd+=x%2;
 if(!((a[0]==1&&(odd==4||odd==6))||(a[0]==2&&odd==0)))return;
 int T=accumulate(a,a+8,0),p=a[7],d=T/2-p;
 if(T%2||p<6||d<8||a[6]>d)return;
 int big=0;for(int i=0;i<7;i++)big+=a[i]>=3;if(big<2)return;
 words++;
 for(int m=1;m<256;m++){int i=__builtin_ctz((unsigned)m);sw[m]=sw[m&(m-1)]+a[i];sp[m]=sw[m]+pc[m];}
 Z mm[256];for(int m=0;m<256;m++)mm[m]=inv(m);
 Z f[256]={};f[0]=mm[255];
 for(int m=1;m<255;m++)if(pc[m]==2||pc[m]==3||(pc[m]==4&&(m&1)))
  f[m]=mm[m]*mm[255^m];
 for(int h=1;h<256;h*=2)for(int st=0;st<256;st+=2*h)for(int j=st;j<st+h;j++){
  Z x=f[j],y=f[j+h];f[j]=x+y;f[j+h]=x-y;
  assert(abs(f[j])<(1LL<<54)&&abs(f[j+h])<(1LL<<54));
 }
 int groups[8]={},ng=0;
 for(int i=0;i<8;i++){if(!i||a[i]!=a[i-1])ng++;groups[ng-1]|=1<<i;}
 for(int sg=0;sg<(1<<ng);sg++){
  int m=0;for(int i=0;i<ng;i++)if(sg>>i&1)m|=groups[i];
  if(pc[m]%2)continue;signs++;
  assert(f[m]>=0);
  if(f[m]<minimum){minimum=f[m];for(int i=0;i<8;i++)witness[i]=(m>>i&1?-a[i]:a[i]);}
 }
}
void finite(){
 for(int m=0;m<256;m++)pc[m]=__builtin_popcount((unsigned)m);
 for(int n=0;n<=110;n++){ch[n][0]=1;for(int k=1;k<=6;k++)ch[n][k]=(n?ch[n-1][k-1]+ch[n-1][k]:0);}
 for(a[0]=1;a[0]<=2;a[0]++)for(a[1]=a[0];a[1]<=4;a[1]++)
 for(a[2]=a[1];a[2]<=4;a[2]++)for(a[3]=a[2];a[3]<=4;a[3]++){
  int D=a[0]+a[1]+a[2]+a[3];
  for(a[4]=a[3];a[4]<=2*D;a[4]++)for(a[5]=a[4];a[5]<=a[4]+D;a[5]++)
  for(a[6]=a[5];a[6]<=a[4]+D;a[6]++)
  for(a[7]=max(a[6],2*a[4]-D);a[7]<=min(a[4]+D,D+a[4]+a[5]-a[6]);a[7]++){
   assert(a[7]<=48);word();
  }
 }
 assert(words==14903&&signs==517426&&2*minimum==270);
 cout<<"finite "<<words<<" words "<<signs<<" signs minPhi "<<2*minimum<<" PASS"<<endl;
 cout<<"minimum witness";for(int x:witness)cout<<" "<<x;cout<<endl;
}
int main(){for(int i=0;i<6;i++)tiles(i);finite();}
'''

obj=os.memfd_create("fm160_obj",0);exe=os.memfd_create("fm160_exe",0)
subprocess.run(["g++","-O3","-std=c++17","-pipe","-x","c++","-","-c",
                "-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin",f"/proc/self/fd/{obj}",
                  "-o",f"/proc/self/fd/{exe}"],text=True,capture_output=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
os.close(obj);os.close(exe)
if args.census_log:
    counts={}
    for line in open(args.census_log):
        m=re.search(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?)  phi=(\d+)",line)
        if m:
            k=len(m[3].split())+1;counts[k]=counts.get(k,0)+1
    assert sum(counts.values())==5430 and counts.get(8)==1944
    assert sum(v for k,v in counts.items() if k<=8)==4040
    print("census",counts,"settled <=8: 4040/5430 PASS")
print("FM-MECH160 PASS")
