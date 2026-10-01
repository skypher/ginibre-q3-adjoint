import os, subprocess, ctypes, glob, sys, time

SRC = r"""
#include <cstdio>
#include <cstdlib>
#include <cstring>
typedef __int128 I;
const int Z=61;
static I P[Z][Z][Z],pref[Z][Z],H[Z];
static int wt[7]={1,1,2,3,3,4,4};
static int sg[7]={1,-1,1,1,-1,1,-1};
static int counts[7],limit,mode;
static unsigned long long profiles,coeffs,pairs,triples;

static I get(int depth,int D,int a,int j){
 if(a<0||a>2*D||j<0)return 0;
 if(a>D)a=2*D-a;
 return j>a?0:P[depth][a][j];
}
static I cumulative(int D,int i){return i<0?0:H[i>D?D:i];}
static I range(int D,int a,int lo,int hi){
 if(a<0||a>2*D)return 0;
 if(a>D)a=2*D-a;
 if(lo>a)return 0;
 if(hi>a)hi=a;
 return pref[a][hi]-(lo?pref[a][lo-1]:0);
}
static void fail(int D,int d,int m,int n){
 printf("FAIL D=%d delta=%d m=%d n=%d counts=",D,d,m,n);
 for(int j=0;j<7;j++)printf("%d,",counts[j]);
 putchar('\n');fflush(stdout);abort();
}
static void visit(int depth,int D,int start){
 profiles++;
 I previous=0;
 for(int d=0;2*d<=D;d++){
  I v=P[depth][2*d][0];
  if(v<previous)fail(D,d,0,0);
  previous=v;coeffs++;
 }
 for(int d=8;2*d<=D;d++)for(int n=5;n<=d;n++){
  I R=P[depth][2*d][0];
  if(d-n-1>=0)R-=P[depth][2*(d-n-1)][0];
  I C=P[depth][2*d-n][n]-P[depth][2*d-n-2][n];
  if(R<C||R< -C)fail(D,d,0,n);
  pairs++;
 }
 if(mode){
  for(int i=0;i<=D;i++){
   H[i]=get(depth,D,2*i,0)+(i?H[i-1]:0);
   I accum=0;
   for(int j=0;j<=D;j++){
    accum+=P[depth][i][j];pref[i][j]=accum;
   }
  }
  for(int d=8;d<D;d++){
   int low=2*d-D;if(low<5)low=5;
   for(int m=low;m<=d+1;m++)for(int n=m;n<=d+1;n++){
    I R=cumulative(D,d)-cumulative(D,d-m-1)
       -cumulative(D,d-n-1)+cumulative(D,d-m-n-2);
    I A=get(depth,D,2*d-n,n)-get(depth,D,2*d-n-2*m-2,n);
    I B=get(depth,D,2*d-m,m)-get(depth,D,2*d-m-2*n-2,m);
    I C=range(D,2*d-n-m,n-m,n+m)
       -range(D,2*d-n-m-2,n-m,n+m);
    I X=A+B,Y=A-B;
    if(X<0)X=-X;if(Y<0)Y=-Y;
    if(R+C<X||R-C<Y)fail(D,d,m,n);
    triples++;
   }
  }
 }
 if(profiles%500000==0){
  printf("progress %llu\n",profiles);fflush(stdout);
 }
 for(int token=start;token<7;token++){
  int n=wt[token],dn=D+n;
  if(dn>limit)continue;
  memset(P[depth+1],0,sizeof(P[depth+1]));
  for(int q=0;q<=n;q++)for(int a=2*q;a<=dn;a++){
   int aa=a-2*q;
   if(aa>2*D)continue;
   if(aa>D)aa=2*D-aa;
   for(int j=aa&1;j<=aa;j+=2)
    P[depth+1][a][j]+=P[depth][aa][j];
  }
  for(int a=n;a<=dn;a++){
   int aa=a-n;
   if(aa>2*D)continue;
   if(aa>D)aa=2*D-aa;
   for(int j=aa&1;j<=aa;j+=2){
    I v=sg[token]*P[depth][aa][j];
    if(!v)continue;
    int lo=j-n;if(lo<0)lo=-lo;
    for(int k=lo;k<=j+n;k+=2)P[depth+1][a][k]+=v;
   }
  }
  counts[token]++;visit(depth+1,dn,token);counts[token]--;
 }
}
extern "C" void run(int cap,int md,unsigned long long* out){
 limit=cap;mode=md;profiles=coeffs=pairs=triples=0;
 memset(P,0,sizeof(P));memset(counts,0,sizeof(counts));
 P[0][0][0]=1;visit(0,0,0);
 out[0]=profiles;out[1]=coeffs;out[2]=pairs;out[3]=triples;
}
"""
def memfd(name,data=b""):
    f=os.memfd_create(name,0)
    if data: os.write(f,data)
    return f

def call(args,fds):
    p=subprocess.run(args,pass_fds=fds,capture_output=True)
    assert p.returncode==0,p.stderr.decode()
    return p.stdout

s=memfd("fm87.cpp",SRC.encode())
assembly=call(["g++","-std=c++17","-O3","-fPIC","-S","-x","c++",
               "-o","-",f"/proc/self/fd/{s}"],(s,))
a=memfd("fm87.s",assembly);o=memfd("fm87.o")
call(["as","--64","-o",f"/proc/self/fd/{o}",f"/proc/self/fd/{a}"],(a,o))
so=memfd("fm87.so")
gcc=sorted(glob.glob("/usr/lib/gcc/x86_64-linux-gnu/*/libgcc.a"))[-1]
call(["ld.gold","-shared","-o",f"/proc/self/fd/{so}",
      f"/proc/self/fd/{o}",gcc],(o,so))
lib=ctypes.CDLL(f"/proc/self/fd/{so}")
lib.run.argtypes=(ctypes.c_int,ctypes.c_int,
                  ctypes.POINTER(ctypes.c_ulonglong))

runs=((24,0),(24,1)) if "--quick" in sys.argv else ((60,0),(48,1))
for cap,mode in runs:
    dp=[1]+[0]*cap
    for w in (1,1,2,3,3,4,4):
        for j in range(w,cap+1):dp[j]+=dp[j-w]
    expected=[
        sum(dp),
        sum(dp[D]*(D//2+1) for D in range(cap+1)),
        sum(dp[D]*sum(d-4 for d in range(8,D//2+1))
            for D in range(cap+1)),
        mode*sum(dp[D]*sum(
            sum(d+2-m for m in range(max(5,2*d-D),d+2))
            for d in range(8,D)) for D in range(cap+1))
    ]
    out=(ctypes.c_ulonglong*4)()
    start=time.monotonic()
    lib.run(cap,mode,out)
    assert list(out)==expected,(list(out),expected)
    print("PASS",cap,mode,list(out),"seconds",time.monotonic()-start)