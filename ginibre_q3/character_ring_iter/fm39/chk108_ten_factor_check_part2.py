import os, subprocess, shlex

src = r'''
#include <bits/stdc++.h>
using namespace std;
using I=long long;
using W=__int128_t;
constexpr int L=10,N=1<<L;

I chooseN(int n,int k){
    if(n<0||k<0||k>n)return 0;
    W x=1;
    for(int j=1;j<=k;j++)x=x*(n-j+1)/j;
    return (I)x;
}

I m0(const array<int,L>&a,int s){
    int r=__builtin_popcount((unsigned)s);
    if(!r)return 1;
    int w=0;
    for(int j=0;j<L;j++)if(s>>j&1)w+=a[j];
    if(r==1||(w&1))return 0;
    W z=0;
    for(int v=s;;v=(v-1)&s){
        int t=w/2+r-2;
        for(int j=0;j<L;j++)if(v>>j&1)t-=a[j]+1;
        I c=chooseN(t,r-2);
        z+=__builtin_parity((unsigned)v)?-c:c;
        if(!v)break;
    }
    return (I)z;
}

array<I,N> table(const array<int,L>&a){
    array<I,N>m{},f{};
    for(int s=0;s<N;s++)m[s]=m0(a,s);
    for(int s=0;s<N;s++)f[s]=m[s]*m[N-1-s];
    for(int h=1;h<N;h*=2)
        for(int b=0;b<N;b+=2*h)
            for(int j=b;j<b+h;j++){
                I x=f[j],y=f[j+h];
                f[j]=x+y;
                f[j+h]=x-y;
            }
    return f;
}

I directPhi(const array<int,L>&a,int signs){
    map<pair<int,int>,W>d,z;
    d[{0,0}]=1;
    for(int k=0;k<L;k++){
        z.clear();
        int e=(signs>>k&1)?-1:1;
        for(auto [ab,v]:d){
            auto [x,y]=ab;
            for(int t=abs(x-a[k]);t<=x+a[k];t+=2)z[{t,y}]+=v;
            for(int t=abs(y-a[k]);t<=y+a[k];t+=2)z[{x,t}]+=e*v;
        }
        d.swap(z);
    }
    return (I)d[{0,0}];
}

bool inBox(const array<int,L>&a){
    int D=accumulate(a.begin(),a.begin()+6,0);
    int q=a[6],r=a[7],s=a[8],p=a[9],v=D+q+r+s-p,large=0;
    for(int j=0;j<6;j++)large+=a[j]>=3;
    return a[5]<=8&&q>=a[5]&&p>=max({s,2*q-D,6})
        &&p<=min(q+D,D+q+r-s)&&v%2==0&&v/2>=8&&s<=v/2
        &&large+(q>=3)+(r>=3)+(s>=3)>=2;
}

struct Stat{
    long long n=0,bad=0;
    I lo=LLONG_MAX,hi=LLONG_MIN;
    array<int,L>w{};
    int mask=0;
};

void add(Stat&z,const array<int,L>&a,const array<I,N>&v,int m){
    I x=v[m];
    z.n++;
    z.bad+=x<0;
    if(x<z.lo){z.lo=x;z.w=a;z.mask=m;}
    z.hi=max(z.hi,x);
}

int main(){
    Stat small;
    array<int,L>a{};
    long long profiles=0;
    function<void(int,int)>go=[&](int j,int lo){
        if(j==L){
            profiles++;
            auto v=table(a);
            for(int s=0;s<N;s++)add(small,a,v,s);
            return;
        }
        for(int n=lo;n<=5;n++){a[j]=n;go(j+1,n);}
    };
    go(0,1);
    cout<<"SMALL profiles="<<profiles<<" signings="<<small.n
        <<" negatives="<<small.bad<<" min="<<small.lo
        <<" max="<<small.hi<<"\n";

    mt19937_64 rng(0x169108);
    vector<array<int,6>>pre;
    array<int,6>b6{};
    function<void(int,int)>g6=[&](int j,int lo){
        if(j==6){pre.push_back(b6);return;}
        for(int n=lo;n<=8;n++){b6[j]=n;g6(j+1,n);}
    };
    g6(0,1);
    Stat box;
    set<array<int,L>>seen;
    long long tries=0;
    while(box.n<1000000){
        tries++;
        auto sm=pre[rng()%pre.size()];
        array<int,L>b{};
        copy(sm.begin(),sm.end(),b.begin());
        int D=accumulate(sm.begin(),sm.end(),0);
        int q=sm[5]+rng()%(2*D-sm[5]+1);
        int r=q+rng()%(D+1);
        int s=r+rng()%(q+D-r+1);
        int lo=max({s,2*q-D,6}),hi=min(q+D,D+q+r-s);
        if(lo>hi)continue;
        lo+=(D+q+r+s-lo)&1;
        if(lo>hi)continue;
        int p=lo+2*(rng()%((hi-lo)/2+1));
        b[6]=q;b[7]=r;b[8]=s;b[9]=p;
        if(!inBox(b)||!seen.insert(b).second)continue;
        auto v=table(b);
        int groups=0,last=-1;
        array<int,L>g{};
        for(int j=0;j<L;j++){
            if(j==0||b[j]!=last){g[groups++]=1<<j;last=b[j];}
            else g[groups-1]|=1<<j;
        }
        for(int u=0;u<(1<<groups)&&box.n<1000000;u++){
            int m=0;
            for(int j=0;j<groups;j++)if(u>>j&1)m|=g[j];
            if(__builtin_parity((unsigned)m))continue;
            add(box,b,v,m);
        }
    }
    cout<<"BOX profiles="<<seen.size()<<" attempts="<<tries
        <<" signings="<<box.n<<" negatives="<<box.bad
        <<" min="<<box.lo<<" max="<<box.hi<<" minmask="<<box.mask<<" word=";
    for(int x:box.w)cout<<x<<",";
    cout<<"\n";

    mt19937 rng2(971);
    long long bridge=0;
    for(int it=0;it<80;it++){
        array<int,L>x{};
        for(int&n:x)n=1+rng2()%12;
        sort(x.begin(),x.end());
        auto v=table(x);
        for(int j=0;j<10;j++){
            int m=rng2()%N;
            assert(v[m]==directPhi(x,m));
            bridge++;
        }
    }
    cout<<"DIRECT_FUSION_BRIDGE="<<bridge<<" PASS\n";
}
'''

obj=os.memfd_create("chk108_core_obj",0)
exe=os.memfd_create("chk108_core_exe",0)
subprocess.run(
    ["g++","-O3","-std=c++17","-pipe","-x","c++","-","-c",
     "-o",f"/proc/self/fd/{obj}"],
    input=src,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
    ["g++","-###","-fno-use-linker-plugin",f"/proc/self/fd/{obj}",
     "-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,pass_fds=(obj,exe),check=True)
link=next(shlex.split(x) for x in p.stderr.splitlines() if "/collect2 " in x)
link[0]="/usr/bin/ld"
subprocess.run(link,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
