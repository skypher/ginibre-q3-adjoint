"""FM-MECH57: identities, counterfamily, outer propagation, large census."""
import argparse, ctypes, os, subprocess
from fractions import Fraction as F
import sympy as S

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--large", action="store_true")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--max-n", type=int, default=80)
args = ap.parse_args()

def data(a, e):
    N = a+e
    sig, de = N+2, a-e
    C = sig*sig-de*de
    c = [1]
    for k in range(N):
        z, rem = divmod(
            de*c[k]-(N-k+1)*(c[k-1] if k else 0), k+1)
        assert rem == 0
        c.append(z)
    v = lambda k: c[k] if 0 <= k <= N else 0
    B = [v(k-1)+v(k+1) for k in range(N+2)]
    D = [v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H = [sig*(v(k)**2+v(k-1)**2)-2*de*v(k)*v(k-1)
         for k in range(N+2)]
    return N, sig, de, C, v, B, D, H

# Fixed polynomial certificate: all m >= 28.
m, z = S.symbols("m z", nonnegative=True)
N = 2*m+3
sig, de, j = N+2, 2*m-3, m+2
C = sig**2-de**2
cc = {-1: S.Integer(-1), 0: S.Integer(1)}
for h in range(5):
    k = j+h
    cc[h+1] = S.factor(
        (de*cc[h]-(N-k+1)*cc[h-1])/(k+1))
B = lambda h: cc[h-1]+cc[h+1]
D = lambda h: S.factor(cc[h]**2-cc[h-1]*cc[h+1])
H = lambda h: S.factor(
    sig*(cc[h]**2+cc[h-1]**2)-2*de*cc[h]*cc[h-1])
ratio = S.factor(
    C*(sig+9)**2*(D(0)-D(4))**2/(324*H(0)*H(4)))
num, den = S.fraction(S.factor(1-ratio))
for p in (num, den):
    pc = S.Poly(S.expand(p.subs(m, 28+z)), z)
    assert pc.degree() == 16
    assert len(pc.all_coeffs()) == 17
    assert min(pc.all_coeffs()) == 81
for h in range(1, 5):
    W = S.factor(cc[0]*B(h)-B(0)*cc[h])
    dot = S.factor(cc[0]*cc[h]+B(0)*B(h))
    for p in (W, dot):
        for q in S.fraction(p):
            assert min(S.Poly(
                S.expand(q.subs(m, 28+z)), z).all_coeffs()) > 0
print("COUNTERFAMILY m >= 28: degree 16, 17 coefficients, "
      "minimum 81; arc checks PASS")

N, sig, de, C, v, B, D, H = data(3, 56)
j, i = 30, 34
X = 2*i-N
R = F(C*(sig+X)**2*(D[j]-D[i])**2,
      4*X*X*H[j]*H[i])
assert R == F(498126628810859275, 517544213794592688)
W = v(j)*B[i]-B[j]*v(i)
assert D[j]-D[i]-abs(W) == 3477468916788060894219093966
print("UNFOLDED WITNESS", (3, 56, j, i), "R", R)
print("E SLACK", D[j]-D[i]-abs(W))

# Exact increment and the outer-region comparison.
sig, X, C, H0, loss, T, delta, Hj = S.symbols(
    "sig X C H0 loss T delta Hj")
t = X/(sig+X)
tn = (X+2)/(sig+X+2)
L = T*T-4*t*t*Hj*H0/C
Ln = (T+delta)**2-4*tn*tn*Hj*(H0-loss)/C
inc = (2*T*delta+delta**2
       +4*Hj*(tn*tn*loss-(tn*tn-t*t)*H0)/C)
assert S.cancel(Ln-L-inc) == 0
assert S.expand(
    X*(X+2)**2-2*(X*(sig+X+2)+sig)
    -X*(X**2-2*sig)-2*(X**2-sig)) == 0
print("SYMBOLIC INCREMENT AND OUTER COMPARISON PASS")

rows = outer = prop = 0
for N in range(8, args.max_n+1):
    for e in range(3, (N-2)//2+1):
        a = N-e
        N, sig, de, C, v, B, D, H = data(a, e)
        rows += 1
        for i in range(N//2+1, N+1):
            X = 2*i-N
            if X*X < C:
                continue
            ell = sig*v(i-1)-de*v(i)
            assert ell*ell >= C*v(i)**2
            assert C >= 2*sig
            t, tn = F(X, sig+X), F(X+2, sig+X+2)
            assert t*t*H[i] >= tn*tn*H[i+1]
            outer += 1
            for j in range(N//2+1, i):
                T, delta = D[j]-D[i], D[i]-D[i+1]
                inc = (2*T*delta+delta**2
                       +F(4, C)*H[j]*
                       (t*t*H[i]-tn*tn*H[i+1]))
                assert inc >= 0
                prop += 1
print("OUTER CHECKS: ROWS", rows, "STEPS", outer,
      "PROPAGATIONS", prop)
print("PASS")

CPP = r'''
extern "C" { void* __dso_handle = nullptr; }
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <omp.h>
using Z=boost::multiprecision::cpp_int;
using R=boost::multiprecision::cpp_dec_float_50;

struct Result {
    long long rows=0,pairs=0,fail=0,efail=0;
    Z num=0,den=1;
    int a=0,e=0,j=0,i=0;
    void offer(const Z& p,const Z& q,int aa,int ee,int jj,int ii){
        if(num==0 || p*den<num*q){
            num=p;den=q;a=aa;e=ee;j=jj;i=ii;
        }
    }
    void add(const Result& b){
        rows+=b.rows;pairs+=b.pairs;
        fail+=b.fail;efail+=b.efail;
        if(b.pairs)offer(b.num,b.den,b.a,b.e,b.j,b.i);
    }
};

Result row(int a,int e){
    int N=a+e,de=a-e,sig=N+2,C=sig*sig-de*de;
    std::vector<Z> c(N+3),B(N+2),D(N+2),H(N+2);
    auto v=[&](int k)->Z{return k<0||k>N?Z(0):c[k];};
    c[0]=1;
    for(int k=0;k<N;k++){
        Z z=de*c[k]-(N-k+1)*v(k-1);
        if(z%(k+1)!=0)throw std::runtime_error("division");
        c[k+1]=z/(k+1);
    }
    for(int k=0;k<=N+1;k++){
        B[k]=v(k-1)+v(k+1);
        D[k]=v(k)*v(k)-v(k-1)*v(k+1);
        H[k]=sig*(v(k)*v(k)+v(k-1)*v(k-1))
             -2*de*v(k)*v(k-1);
    }
    for(int k=N/2+1;k<=N;k++){
        Z ell=sig*v(k-1)-de*v(k);
        if(D[k]<D[k+1] || D[k]<0 ||
           Z(k+1)*(k+1)*(H[k]-H[k+1])
             !=Z(2*k-N)*ell*ell)
            throw std::runtime_error("energy/OL");
    }
    Result out;out.rows=1;
    for(int j=N/2+1;j<=N-4;j++){
        bool crossed=false;
        for(int i=j+1;i<=N;i++){
            Z W=c[j]*B[i]-B[j]*c[i];
            if(W<0 || (W==0 && c[j]*c[i]+B[j]*B[i]<0))
                crossed=true;
            if(i-j<4 || !crossed)continue;
            out.pairs++;
            Z drop=D[j]-D[i];
            int X=2*i-N;
            Z num=Z(C)*(sig+X)*(sig+X)*drop*drop;
            Z den=Z(4)*X*X*H[j]*H[i];
            if(den<=0)throw std::runtime_error("denominator");
            out.fail+=(num<den);
            out.efail+=(drop<abs(W));
            out.offer(num,den,a,e,j,i);
        }
    }
    return out;
}
Z gcdz(Z a,Z b){
    while(b!=0){Z r=a%b;a=b;b=r;}
    return a;
}
void phase(const std::vector<std::pair<int,int>>& jobs,
           const char* name,int threads){
    Result all;int done=0;
    #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
    for(int l=0;l<(int)jobs.size();l++){
        auto [a,e]=jobs[l];
        Result z=row(a,e);
        #pragma omp critical
        {
            all.add(z);done++;
            if(z.fail)
                std::cout<<"FAIL_ROW "<<a<<" "<<e
                         <<" COUNT "<<z.fail<<"\n"<<std::flush;
            if(done%1000==0 ||
               (jobs.size()<1000 && done%32==0))
                std::cout<<"PROGRESS "<<name<<" "<<done
                         <<"/"<<jobs.size()<<" LONG "<<all.pairs
                         <<" JFAIL "<<all.fail<<"\n"<<std::flush;
        }
    }
    std::cout<<"RESULT "<<name<<" ROWS "<<all.rows
             <<" LONG "<<all.pairs<<" JFAIL "<<all.fail
             <<" EFAIL "<<all.efail<<"\n";
    if(all.pairs){
        Z g=gcdz(all.num,all.den);
        std::cout<<"MIN "<<all.a<<" "<<all.e<<" "
                 <<all.j<<" "<<all.i<<" RATIO "
                 <<all.num/g<<"/"<<all.den/g<<"\n";
        R ratio=R(all.num.convert_to<std::string>())
               /R(all.den.convert_to<std::string>());
        std::cout<<std::setprecision(25)
                 <<"DECIMAL "<<ratio<<"\n";
    }
    std::cout<<std::flush;
}
extern "C" int scan(int bound,int count,int rbound,int threads){
    std::vector<std::pair<int,int>> grid,randoms;
    for(int a=5;a<=bound;a++)
        for(int e=3;e<=a-2;e++)grid.push_back({a,e});
    phase(grid,"GRID",threads);
    std::mt19937 rng(20261001);
    for(int l=0;l<count;l++){
        int e=3+rng()%(rbound-4);
        int a=e+2+rng()%(rbound-e-1);
        randoms.push_back({a,e});
    }
    phase(randoms,"RANDOM",threads);
    return 0;
}
extern "C" int lopsided(int threads){
    std::vector<std::pair<int,int>> jobs;
    for(int a:{250,500,1000,2000,4000})
        for(int e=3;e<=8;e++)jobs.push_back({a,e});
    phase(jobs,"LOPSIDED",threads);
    return 0;
}
'''

if args.large:
    obj = os.memfd_create("mech57-object")
    lib = os.memfd_create("mech57-library")
    asm = subprocess.run(
        ["g++", "-O3", "-fPIC", "-fopenmp", "-pipe",
         "-x", "c++", "-S", "-o", "-", "-"],
        input=CPP.encode(), stdout=subprocess.PIPE,
        check=True).stdout
    subprocess.run(
        ["as", "--64", "-o", f"/proc/self/fd/{obj}"],
        input=asm, pass_fds=(obj,), check=True)
    subprocess.run(
        ["ld", "-shared", "-o", f"/proc/self/fd/{lib}",
         f"/proc/self/fd/{obj}"],
        pass_fds=(obj, lib), check=True)
    ctypes.CDLL("libstdc++.so.6", mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    dll = ctypes.CDLL(f"/proc/self/fd/{lib}")
    dll.scan.argtypes = [ctypes.c_int]*4
    dll.lopsided.argtypes = [ctypes.c_int]
    assert dll.scan(200, 256, 1000, args.threads) == 0
    assert dll.lopsided(args.threads) == 0

