
from fractions import Fraction
from math import comb
from functools import lru_cache
from collections import Counter

@lru_cache(None)
def row(a,e):
    N=a+e
    return tuple(sum((-1)**t*comb(e,t)*comb(a,k-t)
                     for t in range(e+1) if 0<=k-t<=a)
                 for k in range(N+1))

def E_from_initial(N,d,j,x,y):
    c={j:Fraction(x),j+1:Fraction(y)}
    for k in (j+1,j+2,j+3):
        c[k+1]=(d*c[k]-(N-k+1)*c[k-1])/(k+1)
    c[j-1]=(d*c[j]-(j+1)*c[j+1])/(N-j+1)
    g=lambda k:c.get(k,Fraction(0))
    D=lambda k:g(k)**2-g(k-1)*g(k+1)
    B=lambda k:g(k-1)+g(k+1)
    return D(j)-D(j+3)+B(j+3)*g(j)-g(j+3)*B(j)

def ABC(N,d,j):
    A=E_from_initial(N,d,j,1,0)
    C=E_from_initial(N,d,j,0,1)
    B=E_from_initial(N,d,j,1,1)-A-C
    return A,B,C

def E_direct(a,e,j):
    N=a+e
    c=row(a,e)
    g=lambda k:c[k] if 0<=k<=N else 0
    D=lambda k:g(k)**2-g(k-1)*g(k+1)
    B=lambda k:g(k-1)+g(k+1)
    return D(j)-D(j+3)+B(j+3)*g(j)-g(j+3)*B(j)

full=Counter()
minimum=None
identity_bad=[]
for a in range(41):
  for e in range(41):
    N=a+e
    c=row(a,e)
    for j in range((N+1)//2,N):
      val=E_direct(a,e,j)
      full['pairs']+=1
      if val<0:
        full['negative']+=1
      if minimum is None or val<minimum[0]:
        minimum=(val,a,e,j)
      if E_from_initial(N,a-e,j,c[j],
                        c[j+1] if j+1<=N else 0)!=val:
        identity_bad.append((a,e,j))

counts=Counter()
cert_bad=[]
denom_bad=[]
nondef_not_outer=[]
nondef=[]
for a in range(41):
  for e in range(3,a+1):
    d=a-e
    if d<2:
      continue
    N=a+e
    c=row(a,e)
    for j in range((N+1)//2,N-1):
      counts['open_folded']+=1
      A,B,C=ABC(N,d,j)
      disc=B*B-4*A*C
      scale=(j+2)**8*(j+3)**6*(j+4)**2*(N-j+1)**2
      if any((z*scale).denominator!=1 for z in (A,B,C,disc)):
        denom_bad.append((a,e,j))
      definite=(A>0 and disc<0)
      counts['definite' if definite else 'nondefinite']+=1
      if not definite:
        nondef.append((a,e,j,A,B,C,disc))
      n=N-j
      U=j+2
      v=n-1
      outer=(d*d>=4*U*v)
      if outer:
        counts['outer']+=1
        lam=C/Fraction(U*v)
        r1=B+lam*d*n
        endpoint=A+B*Fraction(n,d)
        if min(lam,r1,endpoint)<0:
          cert_bad.append((a,e,j,lam,r1,endpoint))
        else:
          counts['outer_certificate_pass']+=1
      if not definite and not outer:
        nondef_not_outer.append((a,e,j))

first=min(nondef,key=lambda z:(z[0]+z[1],z[0],z[1],z[2]))
a,e,j,A,B,C,disc=first
N=a+e
d=a-e
n=N-j
U=j+2
v=n-1
c=row(a,e)
lam=C/Fraction(U*v)
r1=B+lam*d*n
endpoint=A+B*Fraction(n,d)

print('direct grid a,e<=40:',dict(full),'minimum (value,a,e,j)=',minimum)
print('binary identity mismatches=',len(identity_bad),
      'common-denominator failures=',len(denom_bad))
print('folded open grid:',dict(counts),
      'outer-certificate failures=',len(cert_bad),
      'nondefinite-not-outer=',len(nondef_not_outer))
print('first nondefinite (a,e,j)=',(a,e,j),'N,d,n,U,v=',(N,d,n,U,v))
print('A,B,C,disc=',A,B,C,disc)
print('c_j,c_(j+1),ratio=',c[j],c[j+1],Fraction(c[j+1],c[j]))
print('lambda,r1,A+B*n/d,actual E=',lam,r1,endpoint,E_direct(a,e,j))
print('boundary checks (a,e,E,D1 formula):')
for aa,ee in ((0,2),(4,6),(19,3),(40,40)):
    NN=aa+ee
    jj=NN-1
    dd=aa-ee
    print((aa,ee),E_direct(aa,ee,jj),Fraction(dd*dd+NN,2))

def delta_cross(a1,e1,a2,e2,k):
    f=row(a1,e1)
    g=row(a2,e2)
    H=lambda t: Fraction(f[t]*g[t])-Fraction(
        f[t-1]*g[t+1]+g[t-1]*f[t+1],2)
    return H(k)-H(k+1)

print('CS witness deltaF,deltaG,cross=',
      delta_cross(6,6,6,6,7),
      delta_cross(4,8,4,8,7),
      delta_cross(6,6,4,8,7))
print('q2-plus at (4,6,5)=',E_direct(4,6,5))
