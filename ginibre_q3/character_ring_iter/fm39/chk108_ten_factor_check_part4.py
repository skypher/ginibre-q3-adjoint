import os,subprocess,shlex
src=r'''
#include <bits/stdc++.h>
using namespace std;using I=long long;using Z=unsigned long long;using W=__uint128_t;
constexpr int S=18;constexpr Z Q=1ULL<<20,BQ=1ULL<<40;using V=array<Z,S+1>;
I tri(int n){return n>=0?I(n+1)*(n+2)/2:0;}
I rect(int n,int l,int h){return tri(n)-tri(n-l)-tri(n-h)+tri(n-l-h);}
I band(int u,int v,int l,int h,int t){
    return rect((v-u+t)/2+h-1,l,h)-rect((v-u-t)/2+h-2,l,h)
        -rect((t-u-v)/2-1,l,h);
}
V catalog(int L,bool even){
    V best;best.fill(ULLONG_MAX);int stop=(even?3:2)*(L+1);long long count=0;
    for(int l=L+1;l<stop;l+=even?2:1)
    for(int h=1;h<stop;h+=even?2:1)
    for(int d=1-h;d<l;d++)
    for(int p=0;p<(even?1:2);p++){
        int u=max({0,L-h+1-2*d,-2*d});u+=(p-u%2+2)%2;int v=u+2*d;
        I b0=band(u,v,l,h,0);if(!b0)continue;count++;
        for(int s=0;s<=S;s++)best[s]=min(best[s],Z((W)Q*band(u,v,l,h,2*s)/b0));
    }
    cerr<<"CAT L="<<L<<" shapes="<<count<<"\n";return best;
}
vector<I> row(const vector<int>&ns){
    vector<I>f(1,1);
    for(int n:ns){vector<I>g(f.size()+n);
        for(int j=0;j<(int)f.size();j++)if(f[j])
            for(int t=abs(j-n);t<=j+n;t+=2)g[t]+=f[j];
        f.swap(g);
    }
    return f;
}
array<V,5> H2;V H[2];
void h2init(){
    for(int q=1;q<=4;q++){
        int C=S+q+1;auto f=row(vector<int>(q,2));V best;best.fill(ULLONG_MAX);
        long long cnt=0;int maxj=2*(S+q);
        for(int g=-q;g<=C;g++)for(int h=max(g,10-g);h<=C;h++)
        for(int k=h;k<=C;k++){
            if((g-h)%2||(h-k)%2)continue;
            vector<I>a(maxj+1);
            for(int u=0;2*u<=maxj;u++)
                a[2*u]=max(0,min(u,g)-max({-u,u-h-k,-h,-k})+1);
            V b{};
            for(int s=0;s<=S;s++)for(int t=0;t<(int)f.size();t++)if(f[t])
                for(int j=abs(2*s-t);j<=2*s+t&&j<=maxj;j+=2)b[s]+=f[t]*a[j];
            if(!b[0])continue;cnt++;
            for(int s=0;s<=S;s++)best[s]=min(best[s],Z((W)Q*b[s]/b[0]));
        }
        H2[q]=best;
        for(int s=0;s<=S;s++)assert(best[s]>=H[1][s]);
        cerr<<"H2 q="<<q<<" rows="<<cnt<<"\n";
    }
}
struct Key{vector<int>s;int m,r;bool operator<(Key const&o)const{return tie(m,r,s)<tie(o.m,o.r,o.s);}};
map<Key,pair<bool,V>>cache;
V makeprof(const vector<int>&sm,int mode,int r,bool&ok){
    V best;best.fill(ULLONG_MAX);int L=mode?10:9;auto f=row(sm);int D=accumulate(sm.begin(),sm.end(),0);
    if(r==0){
        if(!f[0]){ok=false;return best;}
        for(int s=0;s<=S;s++)best[s]=Z((W)Q*(2*s<(int)f.size()?f[2*s]:0)/f[0]);
        ok=true;return best;
    }
    if(r==1){
        for(int n=L;n<=D;n+=mode?2:1){
            vector<int>x=sm;x.push_back(n);auto full=row(x);if(!full[0])continue;
            for(int s=0;s<=S;s++)best[s]=min(best[s],Z((W)Q*(2*s<(int)full.size()?full[2*s]:0)/full[0]));
        }
    }else if(r==2){
        int M=D+2*S;
        for(int u=D%2;u<=D;u+=2){
            int hi=max(L+1,(M-u)/2+(mode?2:1));
            for(int ell=L+1;ell<=hi;ell+=mode?2:1){
                int top=min(M,u+2*ell-2);top-=((top-u)%2+2)%2;V b{};
                for(int t=0;t<(int)f.size();t++)if(f[t]&&(t-u)%2==0){
                    if(u<=t&&t<=top)b[0]+=f[t];
                    for(int s=1;s<=S;s++){
                        int lo=max(u,abs(t-2*s)),up=min(top,t+2*s);
                        if(lo>up)continue;
                        lo+=((u-lo)%2+2)%2;
                        if(lo<=up)b[s]+=f[t]*((up-lo)/2+1);
                    }
                }
                if(!b[0])continue;
                for(int s=0;s<=S;s++)best[s]=min(best[s],Z((W)Q*b[s]/b[0]));
            }
        }
    }else{ok=false;return best;}
    ok=best[0]!=ULLONG_MAX;return best;
}
V getprof(const vector<int>&sm,int mode,int r,bool&ok){
    if(r>=3){
        ok=true;
        if(mode&&r==3&&!sm.empty()&&all_of(sm.begin(),sm.end(),[](int n){return n==2;}))
            return H2[sm.size()];
        return H[mode];
    }
    Key k{sm,mode,r};auto it=cache.find(k);
    if(it!=cache.end()){ok=it->second.first;return it->second.second;}
    V v=makeprof(sm,mode,r,ok);cache[k]={ok,v};return v;
}
struct Cut{int mask,r,op,k;Z cost;};
struct Res{long long checks=0,fails=0;Z max=0;vector<int>sm,g;int sg=0;};
long long Csmall(int n,int k){
    if(k<0||k>n)return 0;long long z=1;
    for(int j=1;j<=k;j++)z=z*(n-j+1)/j;return z;
}
struct Group{array<int,4>g;long long ways[6][2][2]{};};
vector<Group> groups(int b){
    vector<Group>out;
    for(int a=0;a<=b;a++)for(int c=0;c<=b-a;c++)for(int d=0;d<=b-a-c;d++){
        Group G;G.g={a,c,d,b-a-c-d};
        for(int i=0;i<=G.g[0];i++)for(int j=0;j<=G.g[1];j++)
        for(int k=0;k<=G.g[2];k++)for(int l=0;l<=G.g[3];l++){
            int r=i+j+k+l;if(r<=5)
                G.ways[r][(k+l)&1][(j+l)&1]+=Csmall(G.g[0],i)*Csmall(G.g[1],j)*Csmall(G.g[2],k)*Csmall(G.g[3],l);
        }
        out.push_back(G);
    }
    return out;
}
Res budget(const vector<int>&sm,int mode,const vector<Group>&gs){
    int h=sm.size(),big=10-h,full=(1<<h)-1;vector<Cut>cuts;
    for(int mask=0;mask<=full;mask++){
        vector<int>a,b;int w=0;
        for(int j=0;j<h;j++)if(mask>>j&1){a.push_back(sm[j]);w+=sm[j];}else b.push_back(sm[j]);
        int t=a.size();
        for(int k:{3,4,5}){
            int r=k-t;if(r<0||r>big)continue;
            bool oka,okb;V x=getprof(a,mode,r,oka),y=getprof(b,mode,big-r,okb);
            if(!oka||!okb)continue;
            W dot=0;for(int s=0;s<=S;s++)dot+=(W)x[s]*y[s];if(k==5)dot*=2;
            W num=(W)BQ*Q*Q;Z cost=Z((num+dot-1)/dot);
            cuts.push_back({mask,r,w%2,k,cost});
        }
    }
    sort(cuts.begin(),cuts.end(),[](auto a,auto b){return a.cost>b.cost;});
    Res R;R.sm=sm;
    for(int sg=0;sg<=full;sg++){
        bool good=true;
        for(int i=1;i<h;i++)if(sm[i]==sm[i-1]&&(((sg>>i)^(sg>>(i-1)))&1))good=false;
        if(!good)continue;
        vector<pair<int,int>>pairs;
        for(int i=0;i<h;i++)for(int j=i+1;j<h;j++)if(sm[i]==sm[j])pairs.push_back({i,j});
        vector<int>load(pairs.size());bool paid[32]={};
        for(auto c:cuts)if(c.k==3&&c.r==0&&__builtin_parity((unsigned)(sg&c.mask))){
            int zbest=-1;
            for(int z=0;z<(int)pairs.size();z++){
                auto[i,j]=pairs[z];
                if((c.mask>>i&1)&&(c.mask>>j&1)&&load[z]<2&&(zbest<0||load[z]<load[zbest]))zbest=z;
            }
            if(zbest>=0){load[zbest]++;paid[c.mask]=true;}
        }
        for(auto&G:gs){
            int odd=G.g[2]+G.g[3];for(int n:sm)odd+=n&1;
            if(mode?odd!=0:((odd&1)||odd==0||odd==2))continue;
            if((__builtin_popcount((unsigned)sg)+G.g[1]+G.g[3])&1)continue;
            Z val=0;
            for(auto c:cuts){
                if(c.k==3&&c.r==0&&paid[c.mask])continue;
                int mp=1-__builtin_parity((unsigned)(sg&c.mask));
                val+=c.cost*G.ways[c.r][c.op][mp];
            }
            R.checks++;if(val>BQ)R.fails++;
            if(val>R.max){R.max=val;R.g.assign(G.g.begin(),G.g.end());R.sg=sg;}
        }
    }
    return R;
}
void rec(vector<int>&a,int lo,int h,int mode,long long&checks,long long&fails,Z&maxv){
    if(!h){
        auto gs=groups(10-a.size());Res r=budget(a,mode,gs);
        checks+=r.checks;fails+=r.fails;maxv=max(maxv,r.max);return;
    }
    for(int n=lo;n<=8;n+=mode?2:1){
        if(mode&&n<2)n=2;
        a.push_back(n);rec(a,n,h-1,mode,checks,fails,maxv);a.pop_back();
    }
}
int main(){
    H[0]=catalog(9,false);H[1]=catalog(10,true);h2init();
    for(int mode=0;mode<2;mode++){
        long long checks=0,fails=0;Z maxv=0;
        for(int h=0;h<=5;h++){vector<int>a;rec(a,mode?2:1,h,mode,checks,fails,maxv);}
        cout<<"BUDGET mode="<<mode<<" checks="<<checks<<" fails="<<fails
            <<" max="<<maxv<<"/"<<BQ<<"\n";
    }
}
'''

obj=os.memfd_create("chk108_budget_obj",0)
exe=os.memfd_create("chk108_budget_exe",0)
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
