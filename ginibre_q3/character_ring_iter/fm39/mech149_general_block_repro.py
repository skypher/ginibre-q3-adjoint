import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
from random import Random

ap=argparse.ArgumentParser(description="FM-MECH149 exact block verifier.")
ap.add_argument("--census",action="store_true",help="also run the exact C++ census")
ap.add_argument("--limit",type=int,default=30,help="census maximum W, 22..30")
args=ap.parse_args()

def cg(a,b):
    return range(abs(a-b),a+b+1,2)

@lru_cache(None)
def row(word):
    A={(0,0):1}
    for z in word:
        T=defaultdict(int); n=abs(z); sg=1 if z>0 else -1
        for (a,b),v in A.items():
            for c in cg(a,n): T[c,b]+=v
            for c in cg(b,n): T[a,c]+=sg*v
        A={k:v for k,v in T.items() if v}
    return A

def phi(word):
    return row(tuple(word)).get((0,0),0)

def data(C,I):
    D=sum(map(abs,C)); Q=row(I); A=row(C)
    even=all(z%2==0 for z in C)
    sigma=(-1)**sum(z<0 for z in C)
    h=[Q.get((c,0),0) for c in range(D+1)]
    assert h[0]%2==0
    h[0]//=2
    children=[]
    for a in range(1,D//2+1):
        for b in range(a,D-a+1):
            if even and ((a|b)&1): continue
            t=Q.get((a,b),0)
            if not t: continue
            beta=1 if t>0 else -1
            w=abs(t)
            if a==b:
                assert w%2==0
                w//=2
            children.append((w,(sigma*beta*a,beta*b)))
            for c in cg(a,b): h[c]-=w
    R=sum(v*A.get((c,0),0) for c,v in enumerate(h))
    return h,R,children

rng=Random(149)
for trial in range(36):
    C=tuple(rng.choice((-1,1))*rng.randint(1,3)
            for _ in range(rng.randrange(1,4)))
    r=rng.choice((4,6,8))
    I=[rng.choice((-1,1))*rng.randint(1,4) for _ in range(r-1)]
    sg=(-1)**(sum(z<0 for z in C)+sum(z<0 for z in I))
    n=rng.randint(1,4)
    if (sum(map(abs,C))+sum(map(abs,I))+n)%2: n+=1
    I=tuple(I+[sg*n])
    h,R,children=data(C,I)
    rhs=2*R+sum(w*phi(C+E) for w,E in children)
    assert phi(C+I)==rhs
print("36 general block identities PASS")

def fusion(ns):
    A={0:1}
    for n in ns:
        T=defaultdict(int)
        for a,v in A.items():
            for c in cg(a,n): T[c]+=v
        A=dict(T)
    return A

for D in range(11):
    for c in range(D%2,D+1,2):
        pairs=[(a,b) for a in range(D+1) for b in range(D+1-a)
               if c in cg(a,b)]
        s1=sum(b+1 for a,b in pairs)
        s2=sum((a+1)*(b+1) for a,b in pairs)
        assert s1==F((c+1)*(D-c+2)*(D+c+4),8)
        assert s2==F((c+1)*(D+3)*(D-c+2)*(D+c+4),24)
for trial in range(24):
    D=rng.randrange(4); r=rng.choice((4,6,8))
    t=[0]+[rng.randrange(3) for _ in range(r-1)]
    M=max(t); q=M+2*D+1+rng.randrange(3)
    ns=tuple(q+x for x in t); B=(q-M-2*D+1)//2
    mu=fusion(ns)
    for c in range(sum(ns)%2,D+1,2):
        assert mu.get(c,0)>=(c+1)*B**(r-3)
    for s in range(3,r+1):
        mu=fusion(ns[:s])
        for a in range(D+1):
            assert mu.get(a,0)<=(a+1)*(q+M+1)**(s-3)
# Even a full block of weight at most 44 cannot meet the large-block
# parameter bound: D=M=0 is the most favorable case at fixed r,q.
for r in range(6,45,2):
    for q in range(1,44//r+1):
        b=(q+1)//2; X=q+1
        assert b**(r-3)<X**(r-6)*(comb(r,2)*X+
               2**(r-1)-1-r-comb(r,2))
print("uniform-bound ingredients PASS")

C=(-1,-1); I=(-3,)*8
h,R,_=data(C,I)
assert h[0]==-364 and h[2]==4752 and R==4024
print("no size-independent quartet cutoff:",h)

# J(m,h)=E[(x+y)^(2m)(x-y)^(2h)], using Catalan moments.
def J(m,h):
    T=2*(m+h); ans=0
    for j in range(0,T+1,2):
        co=sum((-1)**k*comb(2*h,k)*comb(2*m,j-k)
               for k in range(max(0,j-2*m),min(2*h,j)+1))
        ans+=co*comb(j,j//2)//(j//2+1)*(
             comb(T-j,(T-j)//2)//((T-j)//2+1))
    return ans

P=[9072000,-5290560,5700756,-1619690,339830,
   -44895,7293,-60,120,5,1]
T=[1814400,12009600,7555236,5645710,1911390,
   753525,204993,37140,5160,425,21]
ev=lambda A,r:sum(a*r**i for i,a in enumerate(A))
shift=[sum(P[i]*comb(i,j)*2**(i-j) for i in range(j,11))
       for j in range(11)]
assert min(shift)>0 and ev(P,1)>0
assert min(T)>0
for r in range(1,12):
    I=(1,)*(2*r-1)+(3,3,-6)
    Q=row(I); base=Q.get((1,0),0)+Q.get((3,0),0)
    cross=Q.get((1,2),0)
    den=(r+4)*(r+9)
    for j in range(5,9): den*=(r+j)**2
    assert (base+2*cross)*den==28*J(r+1,1)*ev(P,r)
    assert (base-2*cross)*den==28*J(r+1,1)*ev(T,r)
    h,R,children=data((1,-2),I)
    assert R==base-2*abs(cross)>0
    assert phi((1,-2)+I)==2*R+sum(
           w*phi((1,-2)+E) for w,E in children)
assert phi((1,-2,1,-2))==6
assert phi((1,-2,-1,2))==2
print("actual-coefficient family, degree-10 identities PASS")

word=(-1,)*15+(-3,-4,-6)
classes=sorted(Counter(word).items())
tested=0; best=None
for ks in product(*(range(m+1) for n,m in classes)):
    C=tuple(n for (n,m),k in zip(classes,ks) for _ in range(k))
    I=tuple(n for (n,m),k in zip(classes,ks) for _ in range(m-k))
    if not C or len(I)<4 or len(I)%2: continue
    h,R,_=data(C,I)
    assert min(h)<0 and R<0
    assert all(row(C).get((c,0),0)>=0
               for c in range(sum(map(abs,C))+1))
    tested+=1
    if best is None or R>best[0]: best=(R,C)
assert tested==55 and best==(-6673282,(-1,-1))
assert phi(word)==57523716
Q=row((-1,)*13+(-3,-4,-6))
assert (Q[0,0],Q[2,0],Q[1,1])==(4657456,9930346,-14174056)
assert 2*(Q[0,0]+Q[2,0]-Q[1,1])==phi(word)
assert row((-2,-2,-3,-3))[1,1]==4
print("55-block obstruction PASS; mixed-sign rescue PASS")
print("FM-MECH149 PROOF VERIFIER PASS")

if args.census:
    import os,shlex,subprocess
    src=r'''#include <bits/stdc++.h>
    using namespace std;
    using Z=__int128_t;
    const int K=46;
    int LIM=24, SMAX=32, MX=8;
    uint64_t place[K];
    struct Bound{int reg=-1,even=-1;};
    struct Off{int a,b;Z w;};
    struct Entry{Bound z;vector<Z> axis;vector<Off> off;};
    unordered_map<uint64_t,Entry> cacheB;
    unsigned long long oracle_calls=0,raw_queries=0;
    struct Counts{long long all=0,pair=0,odd=0,sep=0,block=0,weighted=0,left=0,paritygain=0;};
    Counts stats[K];
    map<int,long long> oddleft,lenleft,onesleft;
    long long minusleft=0,plusleft=0,refplusleft=0;
    vector<int> ns;
    string zs(Z x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;string s;while(x){s+=char('0'+x%10);x/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
    Entry& oracle(array<int,K> cnt,int weight){
     int first=0;for(int n=1;n<=MX;n+=2)if(cnt[n]){first=n;break;}
     if(first&&cnt[first]<0)for(int n=1;n<=MX;n+=2)cnt[n]=-cnt[n];
     uint64_t key=0;for(int n=1;n<=MX;n++)if(cnt[n])key+=place[n]*(2*abs(cnt[n])-(cnt[n]>0));
     auto it=cacheB.find(key);if(it!=cacheB.end())return it->second;
     ++oracle_calls;
     static Z A[K][K],B[K][K];
     memset(A,0,sizeof(A));A[0][0]=1;int s=0;
     for(int n=1;n<=MX;n++)for(int rep=0;rep<abs(cnt[n]);rep++){
      int sg=cnt[n]>0?1:-1,t=s+n;memset(B,0,sizeof(B));
      for(int b=0;b<=s;b++){
       Z pref[2][K+1]{};for(int a=0;a<=s;a++){pref[0][a+1]=pref[0][a];pref[1][a+1]=pref[1][a];pref[a&1][a+1]+=A[a][b];}
       for(int c=0;c<=t;c++){int lo=abs(c-n),hi=min(s,c+n);if(lo<=hi)B[c][b]+=pref[(c+n)&1][hi+1]-pref[(c+n)&1][lo];}
      }
      for(int a=0;a<=s;a++){
       Z pref[2][K+1]{};for(int b=0;b<=s;b++){pref[0][b+1]=pref[0][b];pref[1][b+1]=pref[1][b];pref[b&1][b+1]+=A[a][b];}
       for(int c=0;c<=t;c++){int lo=abs(c-n),hi=min(s,c+n);if(lo<=hi)B[a][c]+=sg*(pref[(c+n)&1][hi+1]-pref[(c+n)&1][lo]);}
      }
      memcpy(A,B,sizeof(A));s=t;
     }
     assert(s==weight);
     Entry data;data.axis.resize(weight+1);for(int c=0;c<=weight;c++)data.axis[c]=A[c][0];
     for(int a=1;a<=min(weight,SMAX-weight)/2;a++)for(int b=a;b<=min(weight,SMAX-weight)-a;b++){
      Z z=A[a][b];if(z<0)z=-z;if(a==b){assert(z%2==0);z/=2;}
      if(z)data.off.push_back({a,b,z});
     }
     Bound ans;
     for(int even=0;even<2;even++){
      if(even&&(weight&1))continue;
      Z h[K]{};h[0]=A[0][0]/2;assert(A[0][0]%2==0);
      for(int c=1;c<=SMAX-weight;c++)h[c]=A[c][0];
      int result=-1;
      for(int d=weight&1;d<=SMAX-weight;d+=2){
       for(int a=1;a<=d/2;a++){
        int b=d-a;if(even&&((a|b)&1))continue;
        Z z=A[a][b];if(z<0)z=-z;if(a==b){assert(z%2==0);z/=2;}
        if(!z)continue;
        for(int c=abs(a-b);c<=a+b;c+=2)h[c]-=z;
       }
       bool ok=true;for(int c=weight&1;c<=d;c+=2)if(h[c]<0){ok=false;break;}
       if(!ok)break;result=d;
      }
      if(even)ans.even=result;else ans.reg=result;
     }
     data.z=ans;return cacheB.emplace(key,move(data)).first->second;
    }
    struct Choice{array<unsigned char,K> c{};int D=0,nf=0;bool even=true;};
    vector<Choice> choices;
    array<int,K> fc;
    void makeChoices(int pos,int D,int nf,int fullnf,Choice &cur,vector<int>&labs){
     if(pos==(int)labs.size()){
      int r=fullnf-nf;
      if(nf&&r>=4&&r%2==0){cur.D=D;cur.nf=nf;cur.even=true;for(int n:labs)if(n%2&&cur.c[n])cur.even=false;choices.push_back(cur);}
      return;
     }
     int n=labs[pos];for(int c=0;c<=fc[n];c++){cur.c[n]=c;makeChoices(pos+1,D+c*n,nf+c,fullnf,cur,labs);}cur.c[n]=0;
    }
    void process(int W){
     if(ns.empty())return;
     int cores=0;for(int n:ns)cores+=n>=3;if(cores<2)return;
     int mx=ns.back(),p0=max(6,mx);if((p0-W)&1)p0++;
     int pmax=min(W-16,W-2*mx);
     if(p0>pmax)return;
     array<int,K> bc{};vector<int>bl;for(int n:ns)bc[n]++;for(int n=1;n<=MX;n++)if(bc[n])bl.push_back(n);
     for(int p=p0;p<=pmax;p+=2){
      vector<int>full=ns;full.push_back(p);sort(full.begin(),full.end());
      int nf=full.size(),S=W+p,odd=0;for(int n:full)odd+=n%2;
      bool sep=false;
      if(nf>=4){int q=full[nf-4],M=full.back()-q,D=S-accumulate(full.end()-4,full.end(),0);sep=q>=D+M+1;}
      fc=bc;fc[p]++;vector<int>labs=bl;if(!bc[p])labs.push_back(p);
      choices.clear();
      if(!sep&&odd!=2){Choice c;makeChoices(0,0,0,nf,c,labs);sort(choices.begin(),choices.end(),[](auto&a,auto&b){return a.D!=b.D?a.D<b.D:a.nf<b.nf;});}
      for(int mask=0;mask<(1<<(int)bl.size());mask++){
       Counts &st=stats[W];st.all++;static unsigned long long progress=0;if(++progress%50000==0)cout<<"PROGRESS "<<progress<<" cache "<<cacheB.size()<<"\n"<<flush;array<int,K>sg{};int sigma=1;
       for(int i=0;i<(int)bl.size();i++){int n=bl[i];sg[n]=(mask>>i&1)?-1:1;if((bc[n]&1)&&sg[n]<0)sigma=-sigma;}
       if(bc[p]&&sg[p]!=sigma){st.pair++;continue;}sg[p]=sigma;
       if(odd==2){st.odd++;continue;}
       if(sep){st.sep++;continue;}
       bool found=false,pg=false;
       for(auto &C:choices){
        array<int,K>cnt{};for(int n:labs)cnt[n]=sg[n]*(fc[n]-C.c[n]);
        ++raw_queries;Bound z=oracle(cnt,S-C.D).z;
        if(z.reg>=C.D){found=true;break;}
        if(C.even&&z.even>=C.D){found=true;pg=true;break;}
       }
       if(found){st.block++;st.paritygain+=pg;}
       else{
        bool weighted=false;
        for(auto &C:choices){
         array<int,K>cnt{},ctx{};for(int n:labs){cnt[n]=sg[n]*(fc[n]-C.c[n]);ctx[n]=sg[n]*C.c[n];}
         Entry &e=oracle(cnt,S-C.D);
         Z h[K]{};h[0]=e.axis[0]/2;
         for(int c=1;c<=C.D&&c<(int)e.axis.size();c++)h[c]=e.axis[c];
         for(auto o:e.off)if(o.a+o.b<=C.D&&(!C.even||((o.a|o.b)%2==0)))
          for(int c=abs(o.a-o.b);c<=o.a+o.b;c+=2)h[c]-=o.w;
         Entry &ctxe=oracle(ctx,C.D);
         Z source=0;for(int c=0;c<=C.D;c++)source+=h[c]*ctxe.axis[c];
         if(source>=0){weighted=true;
          static int shown=0;if(shown++<3){
           cout<<"WEIGHTED W="<<W<<" p="<<p<<" C=";for(int n:labs)if(C.c[n])cout<<sg[n]*n<<"^"<<(int)C.c[n]<<",";
           cout<<" I=";for(int n:labs)if(fc[n]-C.c[n])cout<<sg[n]*n<<"^"<<fc[n]-C.c[n]<<",";
           cout<<" source="<<zs(source)<<"\n"<<flush;
          }
          break;
         }
        }
        if(weighted)st.weighted++;
        else{
         st.left++;oddleft[odd]++;lenleft[nf]++;onesleft[fc[1]]++;
         bool am=true,ap=true,rp=true;for(int n:labs){am&=sg[n]<0;ap&=sg[n]>0;rp&=sg[n]==(n%2?-1:1);}
         minusleft+=am;plusleft+=ap;refplusleft+=rp;
         if(nf%2){cout<<"ODDLEN_LEFT W="<<W<<" p="<<p<<" word=";for(int n:labs)cout<<sg[n]*n<<"^"<<fc[n]<<",";cout<<"\n"<<flush;}
         static int minW=100;static string example;
         bool allminus=true;for(int n:labs)if(sg[n]>0)allminus=false;
         if(allminus&&W<minW){minW=W;cout<<"ALLMINUS_LEFT W="<<W<<" p="<<p<<" word=";for(int n:labs)cout<<-n<<"^"<<fc[n]<<",";cout<<"\n"<<flush;}
        }
       }
      }
     }
    }
    void gen(int last,int W){
     process(W);
     for(int n=last;n<=LIM-W;n++){ns.push_back(n);gen(n,W+n);ns.pop_back();}
    }
    int main(int argc,char**argv){
     if(argc>1&&string(argv[1])=="--help"){cout<<"usage: census [limit <= 30]\n";return 0;}
     LIM=argc>1?atoi(argv[1]):24;assert(LIM<=30&&LIM>=22);MX=LIM-16;SMAX=2*LIM-16;
     uint64_t pos=1;for(int n=1;n<=MX;n++){place[n]=pos;uint64_t b=2*(SMAX/n)+1;assert(pos<=UINT64_MAX/b);pos*=b;}
     auto start=chrono::steady_clock::now();cacheB.reserve(1000000);gen(1,0);
     Counts tot;for(int w=0;w<=LIM;w++){auto&a=stats[w];tot.all+=a.all;tot.pair+=a.pair;tot.odd+=a.odd;tot.sep+=a.sep;tot.block+=a.block;tot.weighted+=a.weighted;tot.left+=a.left;tot.paritygain+=a.paritygain;if(w>=22)cout<<"CUM "<<w<<" "<<tot.all<<" "<<tot.pair<<" "<<tot.odd<<" "<<tot.sep<<" "<<tot.block<<" "<<tot.weighted<<" "<<tot.left<<" paritygain "<<tot.paritygain<<"\n";}
     cout<<"CACHE "<<cacheB.size()<<" queries "<<raw_queries<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
     if(LIM==30){assert(tot.all==777264&&tot.pair==58313&&tot.odd==62652&&tot.sep==902&&tot.block==518062&&tot.weighted==60914&&tot.left==76421);}
     cout<<"SIGNS allminus "<<minusleft<<" allplus "<<plusleft<<" reflectedplus "<<refplusleft<<"\n";
     cout<<"ONES";for(auto[k,v]:onesleft)cout<<" "<<k<<":"<<v;cout<<"\n";
     cout<<"ODD";for(auto [k,v]:oddleft)cout<<" "<<k<<":"<<v;cout<<"\nLEN";for(auto[k,v]:lenleft)cout<<" "<<k<<":"<<v;cout<<"\n";
    }
    '''
    obj=os.memfd_create("fm149_object",0);exe=os.memfd_create("fm149_exec",0)
    subprocess.run(["g++","-O3","-std=c++17","-pipe","-x","c++","-","-c","-o",f"/proc/self/fd/{obj}"],input=src,text=True,pass_fds=(obj,),check=True)
    p=subprocess.run(["g++","-###","-fno-use-linker-plugin",f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],capture_output=True,text=True)
    link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x);link[0]="/usr/bin/ld"
    subprocess.run(link,pass_fds=(obj,exe),check=True)
    subprocess.run([f"/proc/self/fd/{exe}",str(args.limit)],pass_fds=(exe,),check=True)
