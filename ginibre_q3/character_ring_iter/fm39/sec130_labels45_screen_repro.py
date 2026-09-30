import argparse, os, subprocess, sys
from datetime import datetime, timezone
from fractions import Fraction
from math import comb, factorial

argparse.ArgumentParser(
    description="Read-only exact FM3 census and Catalan bridge verifier"
).parse_args()

def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

cpp = r"""
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;
using i128=__int128_t;
using u128=__uint128_t;

struct Type { int n, eps; };
struct Rat { cpp_int num, den; };
struct Best { bool set=false; cpp_int num=0, den=1; string word; long long ties=0; };
struct Sample { int degree; string value; vector<int> counts; };
struct Stats {
    long long leaves=0, words=0, zeros=0, oddzero=0, evenzero=0, negatives=0;
    string zero_word, neg_word;
    Best minpos;
    unsigned long long sampled=0;
    vector<Sample> samples;
    vector<string> even_zero_words;
};

string utc() {
    time_t t=time(nullptr); tm z=*gmtime(&t); char b[32];
    strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%SZ",&z); return b;
}
cpp_int absint(cpp_int a) { return a<0?-a:a; }
cpp_int gcdint(cpp_int a,cpp_int b) {
    a=absint(a); b=absint(b);
    while(b!=0) { cpp_int r=a%b; a=b; b=r; }
    return a;
}
cpp_int from128(i128 v) {
    bool neg=v<0; u128 a=neg?(u128)(-v):(u128)v;
    cpp_int r=(uint64_t)(a>>64); r<<=64; r+=(uint64_t)a;
    return neg?-r:r;
}
string str128(i128 x) { return from128(x).str(); }
string fracstr(cpp_int a,cpp_int b) {
    cpp_int g=gcdint(a,b); if(g!=0) { a/=g; b/=g; }
    return a.str()+"/"+b.str();
}
string wordstr(const vector<Type>& ts,const vector<int>& cs) {
    stringstream s; s<<"("; bool first=true;
    for(size_t i=0;i<ts.size();++i) for(int j=0;j<cs[i];++j) {
        if(!first) s<<","; first=false;
        s<<(ts[i].eps>0?"+":"-")<<ts[i].n;
    }
    s<<")"; return s.str();
}

vector<cpp_int> facts;
vector<vector<Rat>> mus;
void init_mu(int D) {
    facts.resize(2*D+5); facts[0]=1;
    for(size_t i=1;i<facts.size();++i) facts[i]=facts[i-1]*i;
    int M=D/2+2; mus.assign(M,vector<Rat>(M));
    for(int m=0;m<M;++m) for(int k=0;k<M;++k) {
        cpp_int n=2*facts[2*m]*facts[2*m+1]*facts[2*k]*facts[2*k+1];
        cpp_int d=facts[m]*facts[m]*facts[k]*facts[k]
                  *facts[m+k+1]*facts[m+k+2];
        cpp_int g=gcdint(n,d); mus[m][k]={n/g,d/g};
    }
}
vector<i128> mult(const vector<i128>& f,int deg,int D,const Type& t) {
    int side=D+1; vector<i128> out((size_t)side*side,0);
    for(int a=0;a<=deg;++a) for(int b=0;b<=deg-a;++b) {
        i128 v=f[(size_t)a*side+b]; if(!v) continue;
        for(int c=abs(a-t.n);c<=a+t.n;c+=2) out[(size_t)c*side+b]+=v;
        for(int c=abs(b-t.n);c<=b+t.n;c+=2)
            out[(size_t)a*side+c]+=(i128)t.eps*v;
    }
    return out;
}
void offer(Best& b,i128 value,const Rat& mu,const string& w) {
    cpp_int n=from128(value)*mu.den,d=mu.num;
    if(!b.set||n*b.den<b.num*d) {
        cpp_int g=gcdint(n,d); b.num=n/g; b.den=d/g;
        b.word=w; b.set=true; b.ties=1;
    } else if(n*b.den==b.num*d) b.ties++;
}

struct Run {
    int D,K,special,target,sample_cap;
    string name;
    vector<Type> ts;
    Stats st;
    vector<int> counts;
    mt19937_64 rng;
    chrono::steady_clock::time_point last;

    Run(int DD,int KK,int sp,int tg,int scap,string nm,unsigned seed)
      :D(DD),K(KK),special(sp),target(tg),sample_cap(scap),name(nm),rng(seed) {
        if(sp<0) {
            for(int n=1;n<=KK;++n) { ts.push_back({n,1}); ts.push_back({n,-1}); }
        } else {
            int upto=max(3,sp);
            for(int n=1;n<=upto;++n) if(n<=3||n==sp) {
                ts.push_back({n,1}); ts.push_back({n,-1});
            }
        }
        counts.assign(ts.size(),0); last=chrono::steady_clock::now();
    }
    void progress(int at) {
        auto now=chrono::steady_clock::now();
        if(now-last>=chrono::seconds(15)) {
            last=now;
            cerr<<"PROGRESS "<<utc()<<" phase="<<name<<" leaves="<<st.leaves
                <<" eligible_words="<<st.words<<" type_index="<<at
                <<" negatives="<<st.negatives<<"\n"<<flush;
        }
    }
    void leaf(const vector<i128>& f,int degree,int minus,int A,int E,int special_count) {
        st.leaves++;
        if(special>=0&&special_count!=target) return;
        if(degree==0&&special>=0) return;
        if(minus%2) return;
        st.words++; string w=wordstr(ts,counts); i128 val=0;
        if((degree&1)==0) val=f[0];
        if(degree&1) { st.zeros++; st.oddzero++; }
        else if(val<0) { st.negatives++; if(st.neg_word.empty()) st.neg_word=w; }
        else if(val==0) {
            st.zeros++; st.evenzero++;
            if(st.zero_word.empty()) st.zero_word=w;
            if(sample_cap) st.even_zero_words.push_back(w);
        } else {
            int m=A/2,k=E/2;
            if((A&1)||(E&1)||m>=(int)mus.size()||k>=(int)mus.size()) exit(2);
            offer(st.minpos,val,mus[m][k],w);
        }
        if(sample_cap) {
            st.sampled++; Sample x{degree,str128(val),counts};
            if((int)st.samples.size()<sample_cap) st.samples.push_back(std::move(x));
            else {
                uniform_int_distribution<unsigned long long> d(0,st.sampled-1);
                auto j=d(rng);
                if(j<(unsigned long long)sample_cap) st.samples[(size_t)j]=std::move(x);
            }
        }
        if(st.words%250000==0) progress((int)ts.size());
    }
    void dfs(int ix,int deg,int minus,int A,int E,int special_count,
             const vector<i128>& f) {
        if(ix==(int)ts.size()) { leaf(f,deg,minus,A,E,special_count); return; }
        const Type t=ts[ix]; int cap=(D-deg)/t.n;
        if(special==t.n) cap=min(cap,target-special_count);
        vector<i128> cur=f;
        for(int c=0;c<=cap;++c) {
            counts[ix]=c;
            int da=t.eps>0?(t.n%2?c:0):(t.n%2?0:c);
            int de=t.eps<0?c:0;
            dfs(ix+1,deg+c*t.n,minus+(t.eps<0?c:0),A+da,E+de,
                special_count+(special==t.n?c:0),cur);
            if(c<cap) cur=mult(cur,deg+c*t.n,D,t);
        }
    }
    void run() {
        vector<i128> f((size_t)(D+1)*(D+1),0); f[0]=1;
        cerr<<"STAGE "<<utc()<<" phase="<<name<<" types="<<ts.size()
            <<" degree_cap="<<D<<"\n"<<flush;
        dfs(0,0,0,0,0,0,f);
        cout<<"RESULT\t"<<name<<"\tdegree_cap="<<D<<"\teligible="<<st.words
            <<"\tzeros="<<st.zeros<<"\todd_zeros="<<st.oddzero
            <<"\teven_zeros="<<st.evenzero<<"\tnegatives="<<st.negatives
            <<"\tzero_word="<<(st.zero_word.empty()?"none":st.zero_word)
            <<"\tnegative_word="<<(st.neg_word.empty()?"none":st.neg_word)
            <<"\tmin_positive="<<(st.minpos.set?fracstr(st.minpos.num,st.minpos.den):"none")
            <<"\tmin_word="<<(st.minpos.set?st.minpos.word:"none")
            <<"\tmin_ties_seen="<<st.minpos.ties<<"\n";
        if(sample_cap) for(auto& w:st.even_zero_words) cout<<"EZERO\t"<<name<<"\t"<<w<<"\n";
        for(const auto& s:st.samples) {
            cout<<"SAMPLE\t"<<name<<"\t"<<s.degree<<"\t"<<s.value<<"\t";
            for(size_t j=0;j<s.counts.size();++j) {
                if(j) cout<<","; cout<<s.counts[j];
            }
            cout<<"\n";
        }
        cout.flush();
    }
};
int main() {
    init_mu(40);
    Run(40,4,-1,0,256,"sector4",124).run();
    Run(30,5,-1,0,256,"sector5",125).run();
    for(int n=3;n<=12;++n) for(int t=1;t<=2;++t)
        Run(40,n,n,t,0,"tight_n"+to_string(n)+"_count"+to_string(t),
            1000+100*n+t).run();
}
"""

print("BUILD",stamp(),"C++ screen in anonymous memory",flush=True)
asm=subprocess.run(
    ["g++","-std=c++20","-O2","-S","-pipe","-x","c++","-o","-","-"],
    input=cpp.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if asm.returncode:
    print(asm.stderr.decode(),file=sys.stderr); raise SystemExit(asm.returncode)
obj=os.memfd_create("fm130_object",0); exe=os.memfd_create("fm130_executable",0)
r=subprocess.run(["as","--64","-o",f"/proc/self/fd/{obj}","-"],
                 input=asm.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(obj,))
if r.returncode:
    print(r.stderr.decode(),file=sys.stderr); raise SystemExit(r.returncode)
def gp(name): return subprocess.check_output(["g++","-print-file-name="+name],text=True).strip()
paths={x:gp(x) for x in ("crt1.o","crti.o","crtbegin.o","crtend.o","crtn.o")}
cmd=["ld","-o",f"/proc/self/fd/{exe}","-dynamic-linker","/lib64/ld-linux-x86-64.so.2",
     paths["crt1.o"],paths["crti.o"],paths["crtbegin.o"],f"/proc/self/fd/{obj}",
     "-L/usr/lib/gcc/x86_64-linux-gnu/13","-L/usr/lib/x86_64-linux-gnu",
     "-lstdc++","-lm","-lgcc_s","-lgcc","-lc",paths["crtend.o"],paths["crtn.o"]]
r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,pass_fds=(obj,exe))
if r.returncode:
    print(r.stderr.decode(),file=sys.stderr); raise SystemExit(r.returncode)
print("RUN",stamp(),"sector censuses and tightness sweep",flush=True)
p=subprocess.Popen([f"/proc/self/fd/{exe}"],stdout=subprocess.PIPE,stderr=None,pass_fds=(exe,))
out,_=p.communicate()
if p.returncode: raise SystemExit(p.returncode)

samples=[]
for line in out.decode().splitlines():
    if line.startswith("SAMPLE\t"): samples.append(line.split("\t"))
    else: print(line)

# Generate the requested s^A d^E G(Z,P) factor catalog.
import sympy as sp
x,y,s,d,Z,P=sp.symbols("x y s d Z P")
def Usp(n,z):
    return sum((-1)**j*comb(n-j,j)*z**(n-2*j) for j in range(n//2+1))
def symmetric_to_sp(expr):
    terms=dict(sp.Poly(sp.expand(expr),x,y).terms()); out=sp.Integer(0)
    maxr=max((abs(a-b) for a,b in terms),default=0); power=[sp.Integer(2),s]
    for r in range(2,maxr+1): power.append(sp.expand(s*power[-1]-P*power[-2]))
    done=set()
    for (a,b),c in terms.items():
        if (a,b) in done: continue
        done.add((a,b)); done.add((b,a))
        if a==b: out+=c*P**a
        else:
            assert terms.get((b,a),0)==c
            lo,hi=min(a,b),max(a,b); out+=c*P**lo*power[hi-lo]
    return sp.expand(out)
def factor_catalog(n,eps):
    raw=sp.expand(Usp(n,x)+eps*Usp(n,y)); A=E=0
    if eps>0:
        if n%2: A=1; raw=sp.cancel(raw/(x+y))
    else:
        E=1; raw=sp.cancel(raw/(x-y))
        if n%2==0: A=1; raw=sp.cancel(raw/(x+y))
    assert sp.denom(raw)==1
    sym=symmetric_to_sp(raw); out=sp.Integer(0)
    for (r,t),c in sp.Poly(sym,s,P).terms():
        assert r%2==0
        out+=c*(Z+2+2*P)**(r//2)*P**t
    G=sp.Poly(sp.expand(out),Z,P)
    sub={Z:(s**2+d**2)/2-2,P:(s**2-d**2)/4}
    assert sp.expand(s**A*d**E*G.as_expr().subs(sub)
                     -(Usp(n,(s+d)/2)+eps*Usp(n,(s-d)/2)))==0
    return A,E,G
catalog={(n,e):factor_catalog(n,e) for n in range(1,13) for e in (1,-1)}
assert catalog[2,1][2].as_expr()==Z
assert catalog[3,1][2].as_expr()==Z-P and catalog[3,-1][2].as_expr()==Z+P
assert catalog[4,1][2].as_expr()==Z**2+Z-2*P**2
assert catalog[4,-1][2].as_expr()==Z-1
assert catalog[5,1][2].as_expr()==Z**2-Z*P-P**2+2*P-1
assert catalog[5,-1][2].as_expr()==Z**2+Z*P-P**2-2*P-1
print("FACTOR_CATALOG",stamp(),"24 factors through label 12 verified")

def catalan(n): return comb(2*n,n)//(n+1)
def moment_mu(m,k):
    return Fraction(2*factorial(2*m)*factorial(2*m+1)*factorial(2*k)*factorial(2*k+1),
                    factorial(m)**2*factorial(k)**2*factorial(m+k+1)*factorial(m+k+2))
def poly_mul(A,B):
    C={}
    for (i,j),u in A.items():
        for (k,l),v in B.items():
            key=(i+k,j+l); C[key]=C.get(key,0)+u*v
    return {key:v for key,v in C.items() if v}
def Ucoef(n): return {n-2*j:(-1)**j*comb(n-j,j) for j in range(n//2+1)}
def direct_xy(word):
    poly={(0,0):1}
    for n,eps in word:
        f={}
        for r,c in Ucoef(n).items():
            f[(r,0)]=f.get((r,0),0)+c
            f[(0,r)]=f.get((0,r),0)+eps*c
        poly=poly_mul(poly,f)
    return sum(c*catalan(i//2)*catalan(j//2)
               for (i,j),c in poly.items() if i%2==j%2==0)
def scaled_sd_factor(n,eps):
    f={}
    for j in range(n//2+1):
        r=n-2*j; base=(-1)**j*comb(n-j,j)*(1<<(2*j))
        for e in range(r+1):
            c=base*comb(r,e)*(1+eps*((-1)**e))
            if c:
                key=(r-e,e); f[key]=f.get(key,0)+c
    return f
def moment_sd(word):
    poly={(0,0):1}; degree=0
    for n,eps in word:
        poly=poly_mul(poly,scaled_sd_factor(n,eps)); degree+=n
    out=Fraction(0)
    for (a,e),c in poly.items():
        if a%2==e%2==0: out+=c*moment_mu(a//2,e//2)
    return out/(1<<degree)

checked=0
for row in samples:
    _,phase,degree_s,value_s,cs=row
    K=4 if phase=="sector4" else 5
    counts=[int(z) for z in cs.split(",")]; word=[]
    for n in range(1,K+1):
        word.extend([(n,1)]*counts[2*(n-1)])
        word.extend([(n,-1)]*counts[2*(n-1)+1])
    assert direct_xy(word)==moment_sd(word)==int(value_s)
    checked+=1
assert checked==512
print("CATALAN_BRIDGES",stamp(),"512 exact word checks passed")
os.close(obj); os.close(exe)
