# Main agent: independent (u,v)-expansion positivity check of random labels<=4 consumer words.  Usage: SEED COUNT
import random, sys
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache
@lru_cache(None)
def mu(m,k): return Q(2*factorial(2*m)*factorial(2*m+1)*factorial(2*k)*factorial(2*k+1), factorial(m)**2*factorial(k)**2*factorial(m+k+1)*factorial(m+k+2))
def pmul(a,b):
    r={}
    for (i,j),c in a.items():
        for (k,l),d in b.items(): r[(i+k,j+l)]=r.get((i+k,j+l),0)+c*d
    return r
BASE={'Z':{(1,0):Q(1,2),(0,1):Q(1,2),(0,0):Q(-2)},'H':{(1,0):Q(3,4),(0,1):Q(1,4),(0,0):Q(-2)},'G':{(1,0):Q(1,4),(0,1):Q(3,4),(0,0):Q(-2)},
      'Y':{(1,0):Q(1,2),(0,1):Q(1,2),(0,0):Q(-3)}}
# S4 = Z^2+Z-2P^2, P=(u-v)/4
Z=BASE['Z']; P={(1,0):Q(1,4),(0,1):Q(-1,4)}
S4={}
for part,c in ((pmul(Z,Z),1),(Z,1),(pmul(P,P),-2)):
    for k,v in part.items(): S4[k]=S4.get(k,0)+c*v
BASE['S']=S4
@lru_cache(None)
def ppow(name,n):
    if n==0: return {(0,0):Q(1)}
    h=ppow(name,n//2); r=pmul(h,h)
    return pmul(r,BASE[name]) if n%2 else r
random.seed(int(sys.argv[1])); cnt=0; mn=None
for _ in range(int(sys.argv[2])):
    # consumer word: h2^al h3^m h1^a S2^b S3^ga S4^n at level r, al+m<=2r ; A=a+ga+m, E=2r
    r=random.randint(1,8); al=random.randint(0,2*r); m=random.randint(0,2*r-al); ga=random.randint(0,6); n=random.randint(0,5); b=random.randint(0,10)
    a=random.randint(0,8)
    A=a+ga+m
    if A%2: a+=1; A+=1
    E=2*r
    Pp=pmul(pmul(pmul(ppow('Z',b),ppow('H',al)),pmul(ppow('G',ga),ppow('Y',m))),ppow('S',n))
    v=sum(c*mu(A//2+i,E//2+j) for (i,j),c in Pp.items())
    assert v>=0,(r,al,m,a,b,ga,n,v)
    cnt+=1
    rat=v/mu(A//2,E//2)
    if mn is None or rat<mn[0]: mn=(rat,(r,al,m,a,b,ga,n))
print("labels<=4 consumer words checked",cnt,"all >= 0; min value/mu",float(mn[0]),mn[1])
