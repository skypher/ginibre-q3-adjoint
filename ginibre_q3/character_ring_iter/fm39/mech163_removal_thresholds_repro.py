import argparse, itertools, random
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
ap=argparse.ArgumentParser(description="FM-MECH163 exact verifier; memory only.")
ap.add_argument("--census",action="store_true")
args=ap.parse_args()
def ch(n,k):return comb(n,k) if 0<=k<=n else 0
@lru_cache(None)
def profile(u,s,v,neg=False):
    ns=[u+(u%2)]*(s-v)+[u+(1-u%2)]*v
    if s==2:return (1,)*(max(ns)+1)
    if s==3:
        m=min(ns);S=sum(ns)//2
        return tuple(1+min(m,r,S-r) for r in range(S+1))
    ns.sort(reverse=True);hi,third=ns[0],ns[2]
    if s==4 and neg and min(ns)>=2:
        hi=max(hi,min(ns)+2)
        base=(Q(1),Q(8,3),Q(10,3),Q(3))
        return tuple(max(r+1 if 2*r<=third else 1,base[r] if r<4 else Q(0)) for r in range(max(hi,3)+1))
    return tuple(r+1 if 2*r<=third else 1 for r in range(hi+1))
@lru_cache(None)
def cap(u,s,v,t,w,neg=False):return sum(a*b for a,b in zip(profile(u,s,v,neg),profile(u,t,w,neg)))
@lru_cache(None)
def ec(L,o,q,j,s,v):
    return sum(ch(j,i)*ch(o-j,v-i)*ch(q-j,h)*ch(L-o-q+j,s-v-h)
               for i in range(v+1) for h in range(s-v+1) if (i+h)%2)
def coef(u,L):
    worst=Q(0);wd=None
    for o in range(0,L+1,2):
        alphas=[]
        for rp in (0,1):
            if (o if rp else L-o)<2:continue
            oc=o-2*rp;N=L-2;P=Q(1)
            for s in range(2,N//2+1):
                for v in range(0,s+1,2):
                    if v>oc or s-v>N-oc:continue
                    P+=Q(ch(oc,v)*ch(N-oc,s-v),
                         cap(u,s,v,N-s,oc-v)*(2 if 2*s==N else 1))
            alpha=P/(u+(rp-u)%2+1)
            alphas.append((alpha,rp))
        alpha,rp=max(alphas)
        for q in range(0,L+1,2):
            for j in range(max(0,o+q-L),min(o,q)+1):
                beta=Q(0);D=Q(0)
                for s in range(3,L//2+1):
                    for v in range(0,s+1,2):
                        if v>o or s-v>L-o:continue
                        E=ec(L,o,q,j,s,v)
                        if not E:continue
                        z=Q(E,cap(u,s,v,L-s,o-v,True)*(2 if 2*s==L else 1))
                        beta+=z;D+=z*s*(L-s)
                total=alpha+beta+Q(D,4096)
                if total>worst:worst=total;wd=(o,q,j,rp,alpha,beta,D)
    return 1-worst,wd

thresholds={7:3,8:7,9:11,10:15,11:21}
expected={
7:Q(106131,1713920),
8:Q(135652459,1597736448),
9:Q(640107145,3846941696),
10:Q(10856088501515,89337731741696),
11:Q(27696112712089,139929209853440)}
for L,u in thresholds.items():
    rho,_=coef(u,L)
    assert rho==expected[L] and rho>Q(1,20)
    print("threshold",L,u,rho)

def cg(a,b):return range(abs(a-b),a+b+1,2)
def fusion(ns):
    f={0:1}
    for n in ns:
        g=Counter()
        for a,v in f.items():
            for b in cg(a,n):g[b]+=v
        f=g
    return f
checks=0
for r in range(1,7):
    for abc in itertools.combinations_with_replacement(range(2*r,2*r+7),3):
        f=fusion(abc)
        for t in range(6*r+10):
            assert sum(f.get(s,0) for s in cg(t,2*r))>=(r+1)*f.get(t,0)
            checks+=1
rng=random.Random(163)
for _ in range(1000):
    r=rng.randrange(1,9)
    ns=[rng.randrange(2*r,2*r+10) for _ in range(3)]
    ns += [rng.randrange(0,2*r) for _ in range(rng.randrange(9))]
    f=fusion(ns)
    assert f.get(2*r,0)>=(r+1)*f.get(0,0)
    checks+=1
assert checks==16624
f=fusion((1,10,11))
assert (f[0],f[4])==(1,2)

# Negative-quartet profiles from FM-MECH153.
def qprofile(ns):
    if min(ns)>=2:return (Q(1),Q(8,3),Q(10,3),Q(3))
    if sum(n%2 for n in ns)==2:return (Q(1),Q(5,2),Q(5,2),Q(3,2))
    return (Q(1),Q(5,2),Q(3),Q(1))
for ns in itertools.combinations_with_replacement(range(1,9),4):
    f=fusion(ns);d=f.get(0,0)
    if not d or not any(v%2 for v in Counter(ns).values()):continue
    assert max(ns)>=min(ns)+2
    assert all(f.get(2*r,0)>=d*v for r,v in enumerate(qprofile(ns)))

def fwt(v):
    v=v[:];h=1
    while h<len(v):
        for i in range(0,len(v),2*h):
            for j in range(i,i+h):
                a,b=v[j],v[j+h];v[j],v[j+h]=a+b,a-b
        h*=2
    return v
def table(ns):
    N=1<<len(ns);w=[0]*N;k=[0]*N;hi=[0]*N;m=[0]*N;m[0]=1
    for s in range(1,N):
        bit=s&-s;i=bit.bit_length()-1;t=s^bit
        w[s]=w[t]+ns[i];k[s]=k[t]+1;hi[s]=max(hi[t],ns[i])
        if w[s]%2 or 2*hi[s]>w[s] or k[s]==1:continue
        t=s
        while True:
            v=w[s]//2-w[t]-k[t]
            if v>=0:m[s]+=(-1)**k[t]*comb(v+k[s]-2,k[s]-2)
            if not t:break
            t=(t-1)&s
    return m,fwt([m[s]*m[N-1-s] for s in range(N)])
def child_table(ns,m,idx):
    masks=[sum(1<<idx[j] for j in range(len(idx)) if s>>j&1)
           for s in range(1<<len(idx))]
    full=masks[-1]
    return fwt([m[t]*m[full^t] for t in masks])
comparisons=0
for L,u in thresholds.items():
    v=list(range(u,u+L));v[-1]+=sum(v)%2
    for ns in (tuple(v),(u,)*(L-1)+(u+(L*u)%2,)):
        m,values=table(ns);M=m[-1]
        children=[]
        for i,j in itertools.combinations(range(L-1),2):
            if (ns[i]+ns[j])%2:continue
            idx=[h for h in range(L) if h not in (i,j)]
            children.append((idx,child_table(ns,m,idx)))
        for sg in range(1<<L):
            if sg.bit_count()%2:continue
            if any(ns[i]==ns[j] and ((sg>>i^sg>>j)&1)
                   for i,j in itertools.combinations(range(L),2)):continue
            diff=sum(values[sg]-values[sg^(1<<i)^(1<<j)]
                     for i,j in itertools.combinations(range(L),2))
            assert diff%4==0;D=diff//4
            for idx,cv in children:
                cs=sum(((sg>>i)&1)<<j for j,i in enumerate(idx[:-1]))
                cs|=(cs.bit_count()%2)<<(L-3)
                d=values[sg]-cv[cs]
                assert d%2==0;drop=d//2
                assert 20*drop>=M
                assert 20*(4096*drop+D)>=4096*M
                comparisons+=1

def local_profile(ns,negative):
    ns=sorted(ns,reverse=True);h=ns[0]
    if len(ns)==3:
        S=sum(ns)//2
        return [Q(max(0,1+min(ns[-1],2*r,S-r,S-h+r)))
                for r in range(S+1)]
    third=ns[2] if len(ns)>=3 else 0
    f=[Q(r+1 if 2*r<=third else 1) for r in range(h+1)]
    if negative and len(ns)==4:
        for r,x in enumerate(qprofile(ns)):f[r]=max(f[r],x)
    return f
def budget(word):
    ns=tuple(map(abs,word));m,val=table(ns);N=len(ns);full=(1<<N)-1
    sg=sum(1<<i for i,x in enumerate(word) if x<0)
    costs=[Q(0),Q(0)]
    for s in range(1,1<<(N-1)):
        if not m[s]*m[full^s]:continue
        neg=(sg&s).bit_count()%2
        a=[ns[i] for i in range(N) if s>>i&1]
        b=[ns[i] for i in range(N) if not s>>i&1]
        K=sum(x*y for x,y in zip(local_profile(a,neg),local_profile(b,neg)))
        costs[neg]+=1/K
    return m,val,sg,costs
word=(-1,-2,-3,-3,-4,-5,-5,-7)
m,val,sg,costs=budget(word)
cm,cv,cs,ccosts=budget((-1,-2,-3,-3,-4,-7))
beta=costs[1];alpha=(1+ccosts[0])/6
ds=[(val[sg]-val[sg^(1<<i)^(1<<j)])//4
    for i,j in itertools.combinations(range(8),2)]
assert (m[-1],val[sg]//2,cm[-1],cv[cs]//2)==(343,288,15,14)
assert (beta,alpha,1-beta-alpha)==(Q(294743,218025),Q(1,5),-Q(120323,218025))
assert (min(ds),max(ds),sum(ds))==(-181,-2,-863)
print("PASS",checks,"low-channel controls;",comparisons,"removal/mixed checks")

# General factor-count bound, checked here through L=64.
def ceilcube(n):
    a,b=0,1
    while b**3<n:b*=2
    while a+1<b:
        c=(a+b)//2
        if c**3<n:a=c
        else:b=c
    return b
for L in range(7,65):
    u=ceilcube(24*2**L);h=u//2
    H=(h+1)*(h+2)*(2*h+3)//6;T=(h+1)*(h+2)//2
    assert u>=2*L
    charge=Q(5,4)*Q(2**(L-2),H)
    charge+=(1+Q(comb(L-2,2),T)+Q(2**(L-3),H))/(u+1)
    assert charge<Q(1,2)
assert ceilcube(24*2**12)==47

if args.census:
    source=r"""#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <array>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
using namespace std;
using boost::multiprecision::cpp_rational;
struct St {int w=0,n=0,lo=0,odd=0;array<int,3> z={0,0,0};};
St add(St x,int a){x.lo=x.n?min(x.lo,a):a;x.odd+=a%2;x.w+=a;++x.n;for(int j=0;j<3;j++)if(a>x.z[j])swap(a,x.z[j]);return x;}
bool ok(const St&x){return !(x.w%2)&&2*x.z[0]<=x.w;}
int edge(const St&a){return a.n==3?a.w/2:a.z[0];}
int chan(const St&a,int r,bool neg){if(a.n==3){int S=a.w/2;return 6*max(0,1+min({a.z[2],2*r,S-r,S-a.z[0]+r}));}
int v=r<=a.z[0]?(2*r<=a.z[2]?6*(r+1):6):0;
if(neg&&a.n==4&&r<4){
static const int hi[4]={6,16,20,18},mix[4]={6,15,15,9},od[4]={6,15,18,6};
v=max(v,a.lo>=2?hi[r]:a.odd==2?mix[r]:od[r]);}
return v;}
int cap(const St&a,const St&b,bool neg){int R=min(edge(a),edge(b)),v=0;for(int r=0;r<=R;r++)v+=chan(a,r,neg)*chan(b,r,neg);return v;}
cpp_rational cost(const map<int,long long>&d){cpp_rational x=0;for(auto [k,v]:d)x+=cpp_rational(36*v)/k;return x;}
struct Row{int p;vector<int>B;};
int main(int argc,char**argv){
 if(argc>1&&(string(argv[1])=="-h"||string(argv[1])=="--help")){cout<<"Read signed p and background rows on stdin; exact TopPair capacity audit.\n";return 0;}
 vector<Row>rows;string s;
 while(getline(cin,s)){istringstream in(s);Row r;in>>r.p;int x;while(in>>x)r.B.push_back(x);if(!r.B.empty())rows.push_back(r);}
 array<array<array<long long,3>,20>,2> counts{};
 #pragma omp parallel for schedule(dynamic)
 for(int ri=0;ri<(int)rows.size();ri++){
 const auto&r=rows[ri];assert(r.B.size()<16);int n=r.B.size(),L=n+1,full=(1<<L)-1,bg=(1<<n)-1,sm=0;vector<int>ns=r.B;ns.push_back(r.p);
 int W=0;for(int x:r.B)W+=abs(x);
 assert(L<=16&&W<=48&&abs(r.p)<=48&&abs(r.p)<=W);
 map<int,int>signs;for(int x:ns){assert(x&&abs(x)<=48);int e=x>0?1:-1;
 assert(!signs.count(abs(x))||signs[abs(x)]==e);signs[abs(x)]=e;}
 for(int i=0;i<L;i++){if(ns[i]<0)sm|=1<<i;ns[i]=abs(ns[i]);}
 assert(!__builtin_parity((unsigned)sm)&&!( (W+abs(r.p))%2));
 vector<St>a(1<<L);for(int m=1;m<=full;m++){int i=__builtin_ctz((unsigned)m);a[m]=add(a[m^(1<<i)],ns[i]);}
 map<int,long long>neg,negw;
 for(int m=1;m<=bg;m++)if(__builtin_parity((unsigned)(m&sm))&&ok(a[m])&&ok(a[full^m])){
 int k=cap(a[m],a[full^m],true);neg[k]++;negw[k]+=(long long)a[m].n*a[full^m].n;}
 auto beta=cost(neg);int I=-1,J=-1;pair<int,int>key={-1,-1};
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if((ns[i]+ns[j])%2==0){
 auto cur=make_pair(ns[i]+ns[j],max(ns[i],ns[j]));if(cur>key){key=cur;I=i;J=j;}}
 bool hit=false,parent=beta<1,zero=false;
 if(I>=0){int child=full^(1<<I)^(1<<J),cbg=child&bg;cpp_rational alpha=0;
 if(ok(a[child])){map<int,long long>pos;
 for(int m=cbg;m;m=(m-1)&cbg)if(!__builtin_parity((unsigned)(m&sm))&&ok(a[m])&&ok(a[child^m]))pos[cap(a[m],a[child^m],false)]++;
 alpha=(1+cost(pos))/(min(ns[I],ns[J])+1);
 }else zero=true;
 hit=1-beta-alpha>=0;
 }
 #pragma omp critical
 {for(int h=0;h<2;h++)if(W<=(h?48:40)){
 counts[h][L][0]++;counts[h][L][1]+=hit;
 int u=L==7?3:L==8?7:L==9?11:L==10?15:L==11?21:99;
 counts[h][L][2]+=*min_element(ns.begin(),ns.end())>=u;}}
 }
 for(int h=0;h<2;h++){long long total=0,hit=0,uni=0;
 cout<<"W <= "<<(h?48:40)<<"\n";
 for(int L=0;L<20;L++)if(counts[h][L][0]){
 cout<<L;for(auto v:counts[h][L])cout<<" "<<v;cout<<"\n";
 total+=counts[h][L][0];hit+=counts[h][L][1];uni+=counts[h][L][2];}
 assert(total==(h?33487:5430)&&hit==(h?23739:4134)&&uni==(h?2742:348));
 cout<<"TOTAL "<<total<<" "<<hit<<" "<<uni<<"\n";
 }
}"""
    import os,pathlib,re,subprocess
    obj=os.memfd_create("fm163-object");exe=os.memfd_create("fm163-executable")
    subprocess.run(["g++","-std=c++17","-O3","-fopenmp","-pipe","-x","c++","-c",
        "-o",f"/proc/self/fd/{obj}","-"],input=source.encode(),pass_fds=(obj,),check=True)
    lib="/usr/lib/gcc/x86_64-linux-gnu/13/"
    subprocess.run(["ld","-m","elf_x86_64","-dynamic-linker","/lib64/ld-linux-x86-64.so.2",
        "-o",f"/proc/self/fd/{exe}","/usr/lib/x86_64-linux-gnu/crt1.o",
        "/usr/lib/x86_64-linux-gnu/crti.o",lib+"crtbegin.o",f"/proc/self/fd/{obj}",
        "-L"+lib,"-lstdc++","-lm","-lgomp","-lgcc_s","-lgcc","-lc","-lgcc_s","-lgcc",
        lib+"crtend.o","/usr/lib/x86_64-linux-gnu/crtn.o"],pass_fds=(obj,exe),check=True)
    path=pathlib.Path("/tmp/claude-1006/-home-yang-q3adjoint/2d613d32-0be2-46ab-91e8-7dc4d59487ae/scratchpad/flipdesc/fx3b_48.log")
    data=[]
    for line in path.read_text().splitlines():
        m=re.match(r"NOFLIP W=(\d+) p=(-?\d+) B=(.*?) phi=(\d+)",line)
        if m:data.append(m[2]+" "+m[3])
    assert len(data)==33487
    subprocess.run([f"/proc/self/fd/{exe}"],input=("\n".join(data)+"\n").encode(),
        pass_fds=(exe,),env=dict(os.environ,OMP_NUM_THREADS="12"),check=True)
print("FM-MECH163 PASS")
