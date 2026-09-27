// Randomized exact stress test of Conjecture FM3:
//   for F a product of generators S_p (p>=1), H_q (q>=2) and m>=0,
//   column 0 of F * D_1^m has nonnegative coefficients.
// Also reports the smallest ratio signed/positive-part over nonzero cases,
// where positive part = column 0 of F^+ * S_1^m, F^+ = even-b part of F
// (twisted form: column0(F D_1^m) = column0(F^eps S_1^m)).
#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
using boost::multiprecision::cpp_int;
struct Mat { int D; std::vector<cpp_int> v; Mat(int d): D(d), v((size_t)d*d) {} cpp_int& at(int a,int b){return v[(size_t)a*D+b];} const cpp_int& at(int a,int b) const {return v[(size_t)a*D+b];} };
static int LEVEL=-1;
static inline bool ok(int u,int p,int c){ if(LEVEL<0) return true; if(u>LEVEL||p>LEVEL||c>LEVEL) return false; return c<=2*LEVEL-u-p; }
static void mul_side(const Mat&F,int p,int side,long s,Mat&o){
  int D=F.D;
  for(int a=0;a<D;++a)for(int b=0;b<D;++b){const cpp_int&f=F.at(a,b); if(f==0)continue; int u=side?b:a;
    for(int c=(u>p?u-p:p-u);c<=u+p;c+=2){ if(!ok(u,p,c))continue; if(c>=D){std::fprintf(stderr,"DIM\n");std::exit(3);} if(side) o.at(a,c)+=s*f; else o.at(c,b)+=s*f; }}
}
static Mat mulS(const Mat&F,int p){Mat o(F.D);mul_side(F,p,0,1,o);mul_side(F,p,1,1,o);return o;}
static Mat mulH(const Mat&F,int q){Mat o(F.D);for(int a=0;a<=q-1;++a){int b=q-1-a; if(LEVEL>=0&&(a>LEVEL||b>LEVEL))continue; Mat t(F.D);mul_side(F,a,0,1,t);mul_side(t,b,1,1,o);}return o;}
static Mat mulD1(const Mat&F,long sy){Mat o(F.D);mul_side(F,1,0,1,o);mul_side(F,1,1,sy,o);return o;}
int main(int argc,char**argv){
  if(argc<2||!std::strcmp(argv[1],"-h")||!std::strcmp(argv[1],"--help")){
    std::printf("usage: stress_su2_fm3_random SAMPLES NMAX LMAX MMAX [--level k] [--threads t] [--seed s] [--hprob x] [--control]\n"
      "  --control: adjoin X1=V_1xV_1+1 to half the samples (negative control; failures expected).\n"
      "  random products of 1..NMAX generators (H_q with prob hprob, else S_p; labels<=LMAX),\n"
      "  checks column 0 of F*D_1^m >=0 for all m<=MMAX; prints failures, exact zeros summary,\n"
      "  and the minimal signed/positive ratio with its witness.\n"); return 0;}
  long SAMPLES=atol(argv[1]); int NMAX=atoi(argv[2]), LMAX=atoi(argv[3]), MMAX=atoi(argv[4]);
  int threads=16; unsigned long seed=1; double hprob=0.5; bool control=false;
  for(int i=5;i<argc;++i){ if(!strcmp(argv[i],"--level"))LEVEL=atoi(argv[++i]); else if(!strcmp(argv[i],"--threads"))threads=atoi(argv[++i]); else if(!strcmp(argv[i],"--seed"))seed=strtoul(argv[++i],0,10); else if(!strcmp(argv[i],"--hprob"))hprob=atof(argv[++i]); else if(!strcmp(argv[i],"--control"))control=true; }
  omp_set_num_threads(threads);
  auto T0=std::chrono::steady_clock::now();
  long done=0, fails=0, checks=0; double best=1e300; std::string bestw;
  #pragma omp parallel for schedule(dynamic,4) reduction(+:checks)
  for(long s=0;s<SAMPLES;++s){
    std::mt19937_64 rng(seed*1000003ULL+ (unsigned long)s);
    int n=1+(int)(rng()%NMAX); int Lcap = (LEVEL>=0 && LEVEL<LMAX)?LEVEL:LMAX;
    std::string word; int tot=0; std::vector<std::pair<char,int>> g;
    for(int i=0;i<n;++i){ bool H = (std::uniform_real_distribution<double>(0,1)(rng)<hprob) && Lcap>=2;
      int lab = H? 2+(int)(rng()%(Lcap-1)) : 1+(int)(rng()%Lcap); g.push_back({H?'H':'S',lab}); tot+=lab; }
    int D=(LEVEL>=0)?LEVEL+1:tot+MMAX+2;
    Mat F(D); F.at(0,0)=1;
    if(control && (rng()%2)){ g.push_back({'X',1}); }
    for(auto&[k,l]:g){ if(k=='X'){ Mat t(D); mul_side(F,1,0,1,t); Mat o(D); mul_side(t,1,1,1,o); for(size_t z=0;z<o.v.size();++z) o.v[z]+=F.v[z]; F=o; }
      else { F = (k=='S')?mulS(F,l):mulH(F,l); }
      word+=" "; word+=k; word+=std::to_string(l); }
    Mat Fp(D); for(int a=0;a<D;++a)for(int b=0;b<D;b+=2)Fp.at(a,b)=F.at(a,b);
    Mat R=F, P=Fp;
    for(int m=0;m<=MMAX;++m){
      if(m){R=mulD1(R,-1);P=mulD1(P,+1);}  // P tracks F^+ * S_1^m restricted suitably
      for(int a=0;a<D;++a){ const cpp_int& val=R.at(a,0); ++checks;
        if(val<0){
          #pragma omp critical
          { ++fails; if(fails<=20){std::printf("FAIL m=%d F=[%s ] target=%d value=%s\n",m,word.c_str(),a,val.str().c_str());std::fflush(stdout);} }
        } else if(val>0){
          // positive part of the twisted expansion: column0(F^+ S_1^m) where F^+ has even b
          // (computed through P, which uses +1 on y-steps on the even-b part of F)
          const cpp_int& pos=P.at(a,0);
          if(pos>0){ double r = val.convert_to<double>()/pos.convert_to<double>();
            if(r<best){
              #pragma omp critical
              { if(r<best){best=r; bestw="m="+std::to_string(m)+" target="+std::to_string(a)+" F=["+word+" ] value="+val.str()+" pos="+pos.str();} } } }
        }
      }
    }
    #pragma omp atomic
    ++done;
    if(done%2000==0){
      #pragma omp critical
      { double el=std::chrono::duration<double>(std::chrono::steady_clock::now()-T0).count();
        std::printf("progress samples=%ld/%ld fails=%ld best_ratio=%.3e elapsed=%.0fs\n",done,SAMPLES,fails,best,el); std::fflush(stdout);} }
  }
  std::printf("done samples=%ld checks=%ld fails=%ld level=%d\nbest_ratio=%.6e witness: %s\n",SAMPLES,checks,fails,LEVEL,best,bestw.c_str());
  return fails?1:0;
}
