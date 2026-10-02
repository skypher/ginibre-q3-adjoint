import os,shlex,subprocess
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement,product
from math import comb
S=8
def cg(a,b):return range(abs(a-b),a+b+1,2)
@lru_cache(None)
def row(ns):
    a={0:1}
    for n in ns:
        b={}
        for j,v in a.items():
            for k in cg(j,n):b[k]=b.get(k,0)+v
        a=b
    return a
def coefficient_profile(L):
    even=L==8
    ls=range(L+1,3*(L+1),2) if even else range(L+1,2*(L+1))
    hs=range(1,3*(L+1),2) if even else range(1,2*(L+1))
    best=[None]*(S+1);count=0
    for l in ls:
        for h in hs:
            for off in range(1-h,l):
                for par in ((0,) if even else (0,1)):
                    vmin=max(0,L-h+1)
                    u=max(0,vmin-2*off);u+=(par-u)%2
                    I=range(u,u+2*l,2);J=range(u+2*off,u+2*(off+h),2)
                    d=len(set(I)&set(J))
                    if not d:continue
                    count+=1
                    for s in range(S+1):
                        v=F(sum(abs(a-b)<=2*s<=a+b for a in I for b in J),d)
                        best[s]=v if best[s] is None else min(best[s],v)
    expected=(1,2,3,4,5,F(43,8),F(41,8),F(9,2),F(7,2)) if L==7 else (1,2,3,4,5,6,6,F(17,3),5)
    assert tuple(best)==expected
    assert count==(4440 if L==7 else 3393)
    print("coefficient profile",L,count,"PASS",flush=True)
    return tuple(best)
@lru_cache(None)
def compositions(n,k=4):
    if k==1:return ((n,),)
    return tuple((j,)+t for j in range(n+1) for t in compositions(n-j,k-1))
def budgets(L,H):
    even=L==8
    @lru_cache(None)
    def profile(sm,r):
        if r>=3:return H
        f=row(sm);D=sum(sm)
        if r==0:
            d=f.get(0,0)
            return tuple(F(f.get(2*s,0),d) for s in range(S+1)) if d else None
        rows=[]
        if r==1:
            for n in range(L,D+1,2 if even else 1):
                d=f.get(n,0)
                if d:rows.append(tuple(F(sum(f.get(t,0) for t in cg(n,2*s)),d) for s in range(S+1)))
        else:
            for u in range(D%2,D+1,2):
                hi=max(L+1,(D+2*S-u)//2+(2 if even else 1))
                for l in range(L+1,hi+1,2 if even else 1):
                    I=range(u,u+2*l,2);d=sum(f.get(t,0) for t in I)
                    if d:rows.append(tuple(F(sum(v for t,v in f.items() for n in I if abs(n-t)<=2*s<=n+t),d) for s in range(S+1)))
        return tuple(min(v[s] for v in rows) for s in range(S+1)) if rows else None
    @lru_cache(None)
    def cost(sm,big,mask,r):
        a=tuple(x for j,x in enumerate(sm) if mask>>j&1)
        b=tuple(x for j,x in enumerate(sm) if not mask>>j&1)
        p,q=profile(a,r),profile(b,big-r)
        return sum(x*y for x,y in zip(p,q)) if p and q else None
    assert all(H[t]>=2 for t in (1,2,3))
    checks=0;best=F(0)
    for h in range(5):
        big=9-h
        alphabet=(2,4,6) if even else range(1,7)
        for sm in combinations_with_replacement(alphabet,h):
            for groups in compositions(big):
                # Types: even+, even-, odd+, odd-.
                odd=sum(n%2 for n in sm)+groups[2]+groups[3]
                if (odd!=0 if even else odd%2 or odd in (0,2)):continue
                for sg in range(1<<h):
                    if any(sm[i]==sm[j] and (((sg>>i)^(sg>>j))&1) for i in range(h) for j in range(i)):continue
                    if (sg.bit_count()+groups[1]+groups[3])%2:continue
                    budget=F(0)
                    for mask in range(1<<h):
                        t=mask.bit_count();a=tuple(sm[j] for j in range(h) if mask>>j&1)
                        for k in (3,4):
                            r=k-t
                            if r<0 or r>big:continue
                            # Positive equal-pair terms pay these triples.
                            if k==3 and r==0 and len(set(a))<3:continue
                            c=cost(sm,big,mask,r)
                            if not c:continue
                            ways=0
                            for z in compositions(r):
                                if any(z[i]>groups[i] for i in range(4)):continue
                                if (sum(a)+z[2]+z[3])%2:continue
                                if ((mask&sg).bit_count()+z[1]+z[3])%2!=1:continue
                                ways+=comb(groups[0],z[0])*comb(groups[1],z[1])*comb(groups[2],z[2])*comb(groups[3],z[3])
                            budget+=F(ways,c)
                    assert budget<=1,(L,sm,groups,sg,budget)
                    checks+=1;best=max(best,budget)
    expected=(435,F(3163355402463,3869812496225)) if even else (16747,F(32920,47411))
    assert (checks,best)==expected
    print("budgets",L,checks,best,"PASS",flush=True)
for L in (7,8):budgets(L,coefficient_profile(L))

# Exact receipt showing why positive middle terms remain in the finite check.
w=(1,2,2,2,-3,-3,4,5,6)
mi=[row(tuple(abs(w[i]) for i in range(9) if m>>i&1)).get(0,0) for m in range(512)]
N=mi[511];P=pos=neg=0
for m in range(1,511):
    if m.bit_count() not in (2,3,4):continue
    v=mi[m]*mi[511^m]
    minus=sum(w[i]<0 for i in range(9) if m>>i&1)%2
    if m.bit_count()==2:
        assert not (v and minus)
        P+=v
    elif minus:neg+=v
    else:pos+=v
assert (N,P,neg,pos,2*(N+P+pos-neg))==(604,276,902,567,1090)
print("middle-layer receipt PASS",flush=True)

# Shifted-triple lemma: a small exact test, separate from its counting proof.
checks=0
for ns in combinations_with_replacement(range(1,19),3):
    f=row(ns)
    for s in range(1,min(ns)//2+1):
        for r,v in f.items():
            assert sum(f.get(t,0) for t in cg(r,2*s))>=(s+1)*v
            checks+=1
assert checks==47230
for length,expected in ((9,112),(10,246)):
    best=F(0)
    for m in range(0,length+1,2):
        total=F(0)
        for k in range(3,length//2+1):
            q=sum(comb(m,j)*comb(length-m,k-j) for j in range(1,min(m,k)+1,2) if 0<=k-j<=length-m)
            total+=F(q,2 if 2*k==length else 1)
        best=max(best,total)
    assert best==expected
assert sum(k*k for k in range(1,10))==285
assert 21051+7860==28911 and 33487-28911==4576
assert comb(202,7)*512<2**55 and 512*7**5*91**4<2**55
print("shifted triples, cut counts, integer bounds PASS",flush=True)

CPP=r'''
#include <bits/stdc++.h>
using namespace std;using Z=long long;using Wide=__int128_t;
int a[9],pc[512],sw[512],sp[512];Z ch[203][8];
Z words=0,signs=0,minimum=LLONG_MAX;array<int,9>witness;
Z inv(int m){
 int k=pc[m],s=sw[m];if(!k)return 1;if(s%2||k==1)return 0;
 if(k==2){int i=__builtin_ctz((unsigned)m),j=__builtin_ctz((unsigned)(m&(m-1)));return a[i]==a[j];}
 if(k==3){int ma=0;for(int i=0;i<9;i++)if(m>>i&1)ma=max(ma,a[i]);return 2*ma<=s;}
 if(k==4){int b[4],j=0;for(int i=0;i<9;i++)if(m>>i&1)b[j++]=a[i];
  return max(0,(min(b[0]+b[1],b[2]+b[3])-max(abs(b[0]-b[1]),abs(b[2]-b[3])))/2+1);}
 Z out=0;
 for(int t=m;;t=(t-1)&m){
  int top=s/2-sp[t]+k-2;
  if(top>=k-2){assert(top<=202);out+=(pc[t]%2?-1:1)*ch[top][k-2];}
  assert(abs(out)<(1LL<<55));if(!t)break;
 }
 assert(out>=0);return out;
}
void word(){
 int T=accumulate(a,a+9,0);if(T%2)return;
 int p=a[8],delta=T/2-p;
 if(p<6||delta<8||a[7]>delta)return;
 int big=0;for(int i=0;i<8;i++)big+=a[i]>=3;if(big<2)return;
 words++;
 for(int m=1;m<512;m++){int i=__builtin_ctz((unsigned)m);sw[m]=sw[m&(m-1)]+a[i];sp[m]=sw[m]+pc[m];}
 Z mm[512];for(int m=0;m<512;m++)mm[m]=inv(m);
 Z f[512]={};f[0]=mm[511];
 for(int m=1;m<511;m++)if(pc[m]>=2&&pc[m]<=4){
  Wide v=Wide(mm[m])*mm[511^m];assert(v>=0&&v<(Wide(1)<<55));f[m]=Z(v);
 }
 for(int h=1;h<512;h*=2)for(int st=0;st<512;st+=2*h)for(int j=st;j<st+h;j++){
  Z x=f[j],y=f[j+h];assert(Wide(abs(x))+abs(y)<(Wide(1)<<55));
  f[j]=x+y;f[j+h]=x-y;
 }
 int groups[9]={},ng=0;
 for(int i=0;i<9;i++){if(!i||a[i]!=a[i-1])ng++;groups[ng-1]|=1<<i;}
 for(int sg=0;sg<(1<<ng);sg++){
  int m=0;for(int i=0;i<ng;i++)if(sg>>i&1)m|=groups[i];
  if(pc[m]%2)continue;signs++;assert(f[m]>=0);
  if(2*f[m]<minimum){minimum=2*f[m];for(int i=0;i<9;i++)witness[i]=(m>>i&1?-a[i]:a[i]);}
 }
 if(words%500000==0)cout<<"finite progress "<<words<<" "<<signs<<endl;
}
int main(){
 for(int m=0;m<512;m++)pc[m]=__builtin_popcount((unsigned)m);
 for(int n=0;n<=202;n++){ch[n][0]=1;for(int k=1;k<=7;k++)ch[n][k]=n?ch[n-1][k-1]+ch[n-1][k]:0;}
 for(a[0]=1;a[0]<=6;a[0]++)for(a[1]=a[0];a[1]<=6;a[1]++)
 for(a[2]=a[1];a[2]<=6;a[2]++)for(a[3]=a[2];a[3]<=6;a[3]++)for(a[4]=a[3];a[4]<=6;a[4]++){
  int D=a[0]+a[1]+a[2]+a[3]+a[4];
  for(a[5]=a[4];a[5]<=2*D;a[5]++)for(a[6]=a[5];a[6]<=a[5]+D;a[6]++)
  for(a[7]=a[6];a[7]<=a[5]+D;a[7]++)
  for(a[8]=max(a[7],2*a[5]-D);a[8]<=min(a[5]+D,D+a[5]+a[6]-a[7]);a[8]++){
   assert(a[8]<=90);word();
  }
 }
 assert(words==3049184&&signs==160310765&&minimum==578);
 cout<<"finite "<<words<<" "<<signs<<" minPhi "<<minimum<<" PASS"<<endl;
 cout<<"minimum witness";for(int x:witness)cout<<" "<<x;cout<<endl;
}
'''
obj=os.memfd_create("fm164_obj",0);exe=os.memfd_create("fm164_exe",0)
subprocess.run(["g++","-O3","-std=c++17","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=CPP,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(["g++","-###","-fno-use-linker-plugin",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True,check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x);link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
os.close(obj);os.close(exe)
print("FM-MECH164 PASS")
