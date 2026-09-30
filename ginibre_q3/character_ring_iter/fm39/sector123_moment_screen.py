# {1,2,3} sector: E[s^A d^E Z^b (Z+P)^al (Z-P)^ga] for A,E even, al<=E, ga<=A ; also test without constraints
import sys
from math import comb
from fractions import Fraction as Q
def cat(j): return comb(2*j,j)//(j+1)
def bmul(P,R):
    out={}
    for (i,j),a in P.items():
        for (k,l),b in R.items(): out[(i+k,j+l)]=out.get((i+k,j+l),0)+a*b
    return out
def bpow(P,n):
    R={(0,0):1}
    for _ in range(n): R=bmul(R,P)
    return R
s={(1,0):1,(0,1):1}; d={(1,0):1,(0,1):-1}; Z={(2,0):1,(0,2):1,(0,0):-2}
ZP={(2,0):1,(0,2):1,(1,1):1,(0,0):-2}; ZM={(2,0):1,(0,2):1,(1,1):-1,(0,0):-2}
def E(F): return sum(c*cat(i//2)*cat(j//2) for (i,j),c in F.items() if i%2==0 and j%2==0)
K=int(sys.argv[1]); inside=0; neg_in=0; outside=0; neg_out=0; ex=[]
cache={}
for A in range(0,K+1,2):
  for Ee in range(0,K+1,2):
    for b in range(0,K//2+1):
      for al in range(0,K+1):
        for ga in range(0,K+1):
          if A+Ee+2*b+2*al+2*ga>2*K: continue
          F=bmul(bmul(bmul(bpow(s,A),bpow(d,Ee)),bpow(Z,b)),bmul(bpow(ZP,al),bpow(ZM,ga)))
          v=E(F)
          if al<=Ee and ga<=A:
              inside+=1; neg_in+= v<0
              if v<0: print("NEG inside",A,Ee,b,al,ga,v)
          else:
              outside+=1
              if v<0: neg_out+=1; ex.append((A,Ee,b,al,ga,v))
print("inside",inside,"neg",neg_in,"| outside constraints",outside,"neg",neg_out,ex[:8])
