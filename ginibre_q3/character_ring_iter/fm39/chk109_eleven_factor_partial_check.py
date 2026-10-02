import os, subprocess, shlex

src = r'''#include <bits/stdc++.h>
#include <omp.h>
using namespace std; using Z=__int128_t; using U=__uint128_t;
const int S=26,Q=1<<20;
string stamp(){time_t t=time(nullptr);tm v;gmtime_r(&t,&v);char b[32];strftime(b,sizeof b,"%Y-%m-%dT%H:%M:%SZ",&v);return b;}
string dec(Z x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;string s;while(x){s.push_back(char('0'+x%10));x/=10;}if(neg)s.push_back('-');reverse(s.begin(),s.end());return s;}
long long tri(int n){return n>=0?1LL*(n+1)*(n+2)/2:0;}
long long rect(int n,int l,int h){return tri(n)-tri(n-l)-tri(n-h)+tri(n-l-h);}
long long B(int a,int b,int c,int l,int h){return rect(a+h-1,l,h)-rect(-b+h-2,l,h)-rect(c-1,l,h);}
long long band(int u,int v,int l,int h,int t){return rect((v-u+t)/2+h-1,l,h)-rect((v-u-t)/2+h-2,l,h)-rect((t-u-v)/2-1,l,h);}
struct Prof{array<long long,S+1>x;bool ok=false;Prof(){x.fill(LLONG_MAX);}void add(const vector<long long>&v){if(v[0]<=0)return;ok=true;for(int s=0;s<=S;s++){if(v[s]<0)abort();x[s]=min(x[s],(long long)(U(Q)*v[s]/v[0]));}}};
vector<int> lens(int mode){vector<int>v;int lo=mode?15:14,hi=mode?43:27,st=mode?2:1;for(int n=lo;n<=hi;n+=st)v.push_back(n);return v;}
Prof H[2],H4[2];
void makeH(int mode){int L=mode?14:13;auto ls=lens(mode);vector<int>hs;for(int n=1;n<=2*(L+1)-1;n+=mode?2:1)hs.push_back(n);long long count=0;for(int l:ls)for(int h:hs)for(int off=1-h;off<l;off++)for(int parity=0;parity<(mode?1:2);parity++){int u=max({0,L-h+1-2*off,-2*off});u+=(parity-u%2+2)%2;int v=u+2*off;vector<long long>w(S+1);for(int s=0;s<=S;s++)w[s]=band(u,v,l,h,2*s);if(w[0]){H[mode].add(w);count++;}}cout<<stamp()<<" H-mode="<<mode<<" shapes="<<count<<" dyadic=";for(auto n:H[mode].x)cout<<n<<",";cout<<"\n";}
struct Tile{int l,h;};
void makeH4(int mode,int threads){auto ls=lens(mode);vector<Tile>jobs;for(int l:ls)for(int h:ls)jobs.push_back({l,h});vector<Prof>per(jobs.size());vector<long long>counts(jobs.size());atomic<int>done{0};
#pragma omp parallel for num_threads(threads) schedule(dynamic)
for(int job=0;job<(int)jobs.size();job++){int l=jobs[job].l,h=jobs[job].h,R=max({S,l-1,h-1});for(int a=1-h;a<=l-1+2*S;a++)for(int b=max(1-l,-a);b<=h-1+2*S;b++){if(mode&&(a-b)%2)continue;for(int c=-R;c<=min({a,b,l+h-2});c++){if(mode&&(a-c)%2)continue;long long den=B(a,b,c,l,h);if(den<=0)continue;long long pref[2*S+2]={};for(int z=-S;z<=S;z++){long long v=a+b+2*z<0?0:B(a+z,b+z,c+z,l,h);if(v<0)abort();pref[z+S+1]=pref[z+S]+v;}vector<long long>num(S+1);for(int s=0;s<=S;s++){int lo=max(-s,s-a-b);num[s]=pref[s+S+1]-pref[lo+S];}per[job].add(num);counts[job]++;}}
int n=++done;if(n%25==0||n==(int)jobs.size()){
#pragma omp critical
{cout<<stamp()<<" H4 mode="<<mode<<" tiles="<<n<<"/"<<jobs.size()<<"\n";}}
}
H4[mode].ok=true;long long all=0;for(int j=0;j<(int)jobs.size();j++){if(!per[j].ok)abort();all+=counts[j];for(int s=0;s<=S;s++)H4[mode].x[s]=min(H4[mode].x[s],per[j].x[s]);}for(int s=0;s<=S;s++)H4[mode].x[s]=max(H4[mode].x[s],H[mode].x[s]);cout<<stamp()<<" H4 mode="<<mode<<" shapes="<<all<<" dyadic=";for(auto n:H4[mode].x)cout<<n<<",";cout<<"\n";}
vector<Z> row(vector<int>a){vector<Z>v(1,1);for(int n:a){vector<Z>w(v.size()+n);for(int t=0;t<(int)v.size();t++)if(v[t])for(int d=abs(t-n);d<=t+n;d+=2)w[d]+=v[t];v.swap(w);}return v;}
vector<Z> times(const vector<Z>&v,int n){vector<Z>w(v.size()+n);for(int t=0;t<(int)v.size();t++)if(v[t])for(int d=abs(t-n);d<=t+n;d+=2)w[d]+=v[t];return w;}
void check_profiles(mt19937_64&rng){long long checks=0;for(int mode=0;mode<2;mode++)for(int it=0;it<80;it++){vector<int>a;for(int j=0;j<4;j++)a.push_back(mode?14+2*(rng()%34):13+rng()%68);vector<int>back;int nb=rng()%4;for(int j=0;j<nb;j++)back.push_back(1+rng()%12);vector<int>all=a;all.insert(all.end(),back.begin(),back.end());auto f=row(all);for(int s=0;s<=S;s++){auto g=times(f,2*s);for(int t=0;t<(int)f.size();t++){if(U(Q)*g[t]<Z(H4[mode].x[s])*f[t]){cerr<<"H4 coefficient failure\n";exit(7);}checks++;}}}cout<<stamp()<<" direct coefficient tests="<<checks<<" PASS\n";}
vector<Z> conv(vector<Z>p,int n){int T=p.size()-1;vector<Z>q(T+n+1);Z w=0;for(int k=0;k<=T+n;k++){if(k<=T)w+=p[k];if(k-n-1>=0&&k-n-1<=T)w-=p[k-n-1];q[k]=w;}return q;}
struct Eval{vector<int>a;vector<Z>m,f;vector<int>sum;Eval(vector<int>x):a(x){int N=1<<a.size(),all=N-1;sum.assign(N,0);vector<vector<Z>>p(N);m.assign(N,0);p[0]={1};m[0]=1;for(int s=1;s<N;s++){int b=__builtin_ctz((unsigned)s),q=s&(s-1);sum[s]=sum[q]+a[b];p[s]=conv(p[q],a[b]);if(!(sum[s]&1)){int k=sum[s]/2;m[s]=p[s][k]-(k?p[s][k-1]:0);}if(m[s]<0)abort();}f.resize(N);for(int s=0;s<N;s++)f[s]=m[s]*m[all^s];for(int h=1;h<N;h<<=1)for(int st=0;st<N;st+=2*h)for(int j=st;j<st+h;j++){Z x=f[j],y=f[j+h];f[j]=x+y;f[j+h]=x-y;}}};
struct Profile{array<int,11>a;};int mod2(int n){return (n%2+2)%2;}
Profile randprof(mt19937_64&rng){for(;;){Profile p;for(int i=0;i<7;i++)p.a[i]=1+rng()%12;sort(p.a.begin(),p.a.begin()+7);int D=accumulate(p.a.begin(),p.a.begin()+7,0),loq=max(p.a[6],6);int q=loq+rng()%(2*D-loq+1);int r=q+rng()%(D+1);int s=r+rng()%(q+D-r+1);int lo=max({s,6,2*q-D}),hi=min({q+D,D+q+r-s,D+q+r+s-16}),par=mod2(D+q+r+s);if(mod2(lo)!=par)lo++;if(lo>hi)continue;int cores=0;for(int i=0;i<7;i++)cores+=p.a[i]>=3;cores+=(q>=3)+(r>=3)+(s>=3);if(cores<2)continue;int n=(hi-lo)/2+1;p.a[7]=q;p.a[8]=r;p.a[9]=s;p.a[10]=lo+2*(rng()%n);return p;}}
string pkey(const Profile&p){string x;for(int n:p.a){x+=to_string(n);x+=',';}return x;}
bool validbox(const Profile&p){int D=accumulate(p.a.begin(),p.a.begin()+7,0),q=p.a[7],r=p.a[8],s=p.a[9],z=p.a[10],cores=0;for(int i=0;i<7;i++)cores+=p.a[i]>=3;cores+=(q>=3)+(r>=3)+(s>=3);int d=(D+q+r+s-z)/2;return is_sorted(p.a.begin(),p.a.end())&&q<=2*D&&z<=3*D&&z>=max({s,6,2*q-D})&&z<=min({q+D,D+q+r-s,D+q+r+s-16})&&((D+q+r+s+z)%2==0)&&d>=8&&s<=d&&cores>=2;}
void reductions(mt19937_64&rng){long long pairs=0,odds=0,pay=0;for(int it=0;it<80;it++){vector<int>c;int len=1+rng()%8;for(int i=0;i<len;i++)c.push_back(1+rng()%9);int cs=0;for(int i=0;i<len;i++)if(rng()%2)cs|=1<<i;int n=1+rng()%8;vector<int>x=c;x.push_back(n);x.push_back(n);Eval e(x);Z lhs=e.f[cs|(1<<(len+1))],rhs=0;for(int k=1;k<=n;k++){vector<int>v=c;v.push_back(2*k);Eval q(v);rhs+=q.f[cs|(1<<len)];}if(lhs!=rhs)exit(8);pairs++;int a=1+2*(rng()%5),b=1+2*(rng()%5);vector<int>ec;for(int i=0;i<len;i++)ec.push_back(2*(1+rng()%5));int es=0;for(int i=0;i<len;i++)if(rng()%2)es|=1<<i;int ea=rng()%2,eb=rng()%2;vector<int>top=ec;top.push_back(a);top.push_back(b);Eval u(top);Z lval=u.f[es|(ea<<len)|(eb<<(len+1))],rval=0;for(int j=abs(a-b);j<=a+b;j+=2){vector<int>child=ec;child.push_back(j);Eval d(child);rval+=d.f[es|((ea^eb)<<len)];}if(lval!=rval)exit(9);odds++;}
for(int it=0;it<100;it++){vector<int>a(11);for(int&n:a)n=1+rng()%12;sort(a.begin(),a.end());if(accumulate(a.begin(),a.end(),0)%2)a.back()++;Eval e(a);int sg=0;for(int i=0;i<11;){int j=i;while(j<11&&a[j]==a[i])j++;if(rng()%2)for(int k=i;k<j;k++)sg|=1<<k;i=j;}if(__builtin_popcount((unsigned)sg)%2){for(int i=0;i<11;){int j=i;while(j<11&&a[j]==a[i])j++;if((j-i)%2){for(int k=i;k<j;k++)sg^=1<<k;break;}i=j;}}
Z rhs=2*e.m[2047];for(int k=1;k<=5;k++)for(int m=0;m<2048;m++)if(__builtin_popcount((unsigned)m)==k){Z v=e.m[m]*e.m[2047^m];if(k==2)rhs+=2*v;else rhs+=2*((__builtin_parity((unsigned)(sg&m))?-1:1)*v);}if(rhs!=e.f[sg])exit(10);pay++;}cout<<stamp()<<" pair="<<pairs<<" two-odd="<<odds<<" payment="<<pay<<" PASS\n";}
void bridges(mt19937_64&rng){long long n=0;for(int it=0;it<160;it++){Profile p=randprof(rng);vector<int>a(p.a.begin(),p.a.end());Eval e(a);for(int t=0;t<64;t++){int mask=rng()%2048;vector<int>x;for(int j=0;j<11;j++)if(mask>>j&1)x.push_back(a[j]);if(row(x)[0]!=e.m[mask])exit(11);n++;}}cout<<stamp()<<" fusion bridges="<<n<<" PASS\n";}
void small_labels(){long long profiles=0,signings=0;Z minv=Z(1)<<120;array<int,11>a{};function<void(int,int)>rec=[&](int i,int lo){if(i==11){Eval e(vector<int>(a.begin(),a.end()));for(Z x:e.f){if(x<0)exit(12);minv=min(minv,x);signings++;}profiles++;return;}for(int n=lo;n<=4;n++){a[i]=n;rec(i+1,n);}};rec(0,1);if(profiles!=364||signings!=745472)exit(12);cout<<stamp()<<" labels1..4 profiles="<<profiles<<" all-signings="<<signings<<" minimum="<<dec(minv)<<" PASS\n";}
void sample_box(mt19937_64&rng){long long vals=0;Z mn=Z(1)<<120;set<string>seen;while(seen.size()<512){Profile p=randprof(rng);if(!validbox(p))exit(13);if(!seen.insert(pkey(p)).second)continue;Eval e(vector<int>(p.a.begin(),p.a.end()));for(Z x:e.f){if(x<0){cerr<<"sample negative ";for(int n:p.a)cerr<<n<<",";cerr<<" value="<<dec(x)<<"\n";exit(13);}mn=min(mn,x);vals++;}}cout<<stamp()<<" residual profiles="<<seen.size()<<" signing-values="<<vals<<" minimum="<<dec(mn)<<" PASS\n";Profile b[3];for(int i=0;i<7;i++){b[0].a[i]=12;b[1].a[i]=12;b[2].a[i]=3;}b[0].a[7]=12;b[0].a[8]=12;b[0].a[9]=12;b[0].a[10]=96;b[1].a[7]=168;b[1].a[8]=252;b[1].a[9]=252;b[1].a[10]=252;b[2].a[7]=42;b[2].a[8]=63;b[2].a[9]=63;b[2].a[10]=63;long long n=0;Z bm=Z(1)<<120;for(auto&p:b){if(!validbox(p))exit(14);Eval e(vector<int>(p.a.begin(),p.a.end()));Z x=*min_element(e.f.begin(),e.f.end());if(x<0)exit(14);n+=e.f.size();bm=min(bm,x);cout<<stamp()<<" boundary="<<pkey(p)<<" signs="<<e.f.size()<<" min="<<dec(x)<<"\n";}cout<<stamp()<<" boundary-signings="<<n<<" minimum="<<dec(bm)<<" PASS\n";}
int main(){mt19937_64 rng(0x173109);U product=1;for(int i=0;i<7;i++)product*=13;for(int i=0;i<4;i++)product*=253;if(U(2048)*product>=(U(1)<<127))return 2;cout<<stamp()<<" dimension-bound="<<dec((Z)product)<<" Walsh-bound="<<dec((Z)(U(2048)*product))<<"\n";
long long bchecks=0;for(int it=0;it<30000;it++){int l=1+rng()%43,h=1+rng()%43,a=1-h+rng()%(l+h+4),b=1-l+rng()%(l+h+4),c=-50+rng()%100;if(a+b<0||c>min({a,b,l+h-2}))continue;long long actual=0;for(int x=0;x<l;x++)for(int y=0;y<h;y++)if(x-y<=a&&y-x<=b&&x+y>=c)actual++;if(B(a,b,c,l,h)!=actual)return 3;bchecks++;}
long long bandchecks=0;for(int it=0;it<1000;it++){int l=14+rng()%30,h=14+rng()%30,u=rng()%60,v=rng()%60;if((u-v)%2)v++;int t=2*(rng()%40);long long actual=0;for(int x=0;x<l;x++)for(int y=0;y<h;y++){int aa=u+2*x,bb=v+2*y;if(abs(aa-bb)<=t&&t<=aa+bb)actual++;}if(band(u,v,l,h,t)!=actual)return 4;bandchecks++;}cout<<stamp()<<" tile checks B="<<bchecks<<" band="<<bandchecks<<" PASS\n";
for(int mode=0;mode<2;mode++){makeH(mode);for(int s=1;s<=6;s++)if(H[mode].x[s]<2*Q)return 5;}for(int mode=0;mode<2;mode++)makeH4(mode,12);check_profiles(rng);reductions(rng);bridges(rng);small_labels();sample_box(rng);cout<<stamp()<<" independent screens PASS\n";}
'''

obj=os.memfd_create("fm_chk109_obj",0)
exe=os.memfd_create("fm_chk109_exe",0)
subprocess.run(
    ["g++","-O3","-std=c++17","-fopenmp","-pipe","-x","c++","-","-c",
     "-o",f"/proc/self/fd/{obj}"],
    input=src,text=True,pass_fds=(obj,),check=True)
p=subprocess.run(
    ["g++","-###","-fno-use-linker-plugin","-fopenmp",
     f"/proc/self/fd/{obj}","-o",f"/proc/self/fd/{exe}"],
    capture_output=True,text=True,pass_fds=(obj,exe),check=True)
cmd=next(shlex.split(line) for line in p.stderr.splitlines()
         if "/collect2 " in line)
cmd[0]="/usr/bin/ld"
subprocess.run(cmd,pass_fds=(obj,exe),check=True)
subprocess.run([f"/proc/self/fd/{exe}"],pass_fds=(exe,),check=True)
