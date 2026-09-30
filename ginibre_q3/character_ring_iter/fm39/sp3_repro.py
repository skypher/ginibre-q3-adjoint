
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

# Original-word evaluator using two semicircle variables and exact Catalan moments.
def add(A,B,s=1):
    C=Counter(A)
    for k,v in B.items(): C[k]+=s*v
    return {k:v for k,v in C.items() if v}
def mul(A,B):
    C=Counter()
    for (i,j),v in A.items():
        for (k,l),w in B.items(): C[i+k,j+l]+=v*w
    return {k:v for k,v in C.items() if v}
def power(A,n):
    R={(0,0):1}
    for _ in range(n): R=mul(R,A)
    return R
X={(1,0):1}; Y={(0,1):1}
@lru_cache(None)
def U(n,var):
    V=X if var==0 else Y
    if n==0: return {(0,0):1}
    if n==1: return V
    return add(mul(V,U(n-1,var)),U(n-2,var),-1)
@lru_cache(None)
def h(k):
    R={}
    for j in range(k+1): R=add(R,mul(U(j,0),U(k-j,1)))
    return R
@lru_cache(None)
def shat(k): return add(U(k,0),U(k,1))
@lru_cache(None)
def cat(n): return 0 if n%2 else comb(2*(n//2),n//2)//(n//2+1)
def moment(F): return sum(v*cat(i)*cat(j) for (i,j),v in F.items())
def direct(r,a,k,L):
    fs=([h(L[0]),h(L[1]),shat(L[2])] if k=='A' else
        [h(L[0]),shat(L[1]),shat(L[2])] if k=='B' else
        [shat(L[0]),shat(L[1]),shat(L[2])])
    F={(0,0):1}
    for f in fs: F=mul(F,f)
    K=mul(power(add(X,Y),a),power(add(X,Y,-1),2*r))
    v=moment(mul(F,K)); assert v%2==0
    return v//2

# Closed row formula using W_c(p,q)=B_i c_j-c_i B_j and SU(2) fusion.
def cg(p,q): return range(abs(p-q),p+q+1,2)
def fuse(labels):
    out=Counter({0:1})
    for ell in labels:
        z=Counter()
        for p,m in out.items():
            for t in cg(p,ell): z[t]+=m
        out=z
    return out
def row(a,e):
    c=[1]
    for sign,n in ((1,a),(-1,e)):
        for _ in range(n):
            z=[0]*(len(c)+1)
            for k,v in enumerate(c): z[k]+=v; z[k+1]+=sign*v
            c=z
    return c
def W(c,p,q,e):
    N=len(c)-1
    if p+q>N or (N+p+q)%2: return 0
    if p<q: return (-1)**e*W(c,q,p,e)
    C=lambda k:c[k] if 0<=k<=N else 0
    j=(N+p-q)//2; i=(N+p+q)//2+1
    return (C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
def Wprod(c,L,R,e):
    return sum(m*n*W(c,p,q,e)
               for p,m in fuse(L).items() for q,n in fuse(R).items())
def closed(r,a,k,L):
    m={'A':2,'B':1,'C':0}[k]; e=2*r-m; c=row(a,e); eps=e%2
    if k=='A':
        u,v,p=L; A,B=u+1,v+1
        return (Wprod(c,(A,B,p),(),eps)+Wprod(c,(A,B),(p,),eps)
                -Wprod(c,(A,p),(B,),eps)-Wprod(c,(B,p),(A,),eps))
    if k=='B':
        u,p,q=L; A=u+1
        return (Wprod(c,(A,p,q),(),eps)+Wprod(c,(A,),(p,q),eps)
                +Wprod(c,(A,p),(q,),eps)+Wprod(c,(A,q),(p,),eps))
    p,q,s=L
    return (Wprod(c,(p,q,s),(),eps)+Wprod(c,(p,q),(s,),eps)
            +Wprod(c,(p,s),(q,),eps)+Wprod(c,(q,s),(p,),eps))

def g0_window(r,a,k,L):
    m={'A':2,'B':1,'C':0}[k]; e=2*r-m; N=a+e
    signed=([(L[0]+1,-1),(L[1]+1,-1),(L[2],1)] if k=='A' else
            [(L[0]+1,-1),(L[1],1),(L[2],1)] if k=='B' else
            [(L[0],1),(L[1],1),(L[2],1)])
    ((Xv,sx),(Yv,sy),(Zv,sz))=sorted(signed,key=lambda t:(-t[0],t[1]))
    if (N+Xv+Yv+Zv)%2 or Xv+Yv-Zv<=N: return None
    n=N+Xv-Yv-Zv
    if n%2: return None
    x=n//2
    if x<0 or x+Zv>N: return None
    sigma=((1 if sx==sy==-1 else -1) if k=='A' else
           (-1 if sz==-1 else 1) if k=='B' else 1)
    c=row(a,e); at=lambda j:c[j] if 0<=j<=N else 0
    D=lambda j:at(j)**2-at(j-1)*at(j+1)
    S=sum(D(j) for j in range(x,x+Zv+1))
    T=at(x)*at(x+Zv)-at(x-1)*at(x+Zv+1)
    return Xv,Yv,Zv,x,sigma,S,T
def uniform_band(r,a,k,Z):
    e=2*r-{'A':2,'B':1,'C':0}[k]; d=a-e; N=a+e
    # Z=2 uses rho<=1/2; Z>=3 uses the broader energy band.
    return d*d <= (N+1 if Z==2 else 2*(N+1))

# Requested boundary tests: smallest labels, equal labels, and suffix a=0.
print('boundary: kind r a labels direct closed')
for r in (1,2):
  for k in 'ABC':
    for a in (0,1):
      L=(2,2,2); d=direct(r,a,k,L); f=closed(r,a,k,L); assert d==f
      print(k,r,a,L,d,f)
for r in (1,2):
  for q in (2,3,4):
    for k in 'ABC':
      L=(q,q,q); d=direct(r,0,k,L); f=closed(r,0,k,L); assert d==f
      print(k,r,0,L,d,f)

# Independent formula check on an exhaustive small box.
checks=0; labels=range(2,7)
for r in range(1,5):
  for a in range(7):
    for uv in combinations_with_replacement(labels,2):
      for p in labels:
        L=(uv[0],uv[1],p)
        assert direct(r,a,'A',L)==closed(r,a,'A',L); checks+=1
    for pq in combinations_with_replacement(labels,2):
      for u in labels:
        L=(u,pq[0],pq[1])
        assert direct(r,a,'B',L)==closed(r,a,'B',L); checks+=1
    for L in combinations_with_replacement(labels,3):
      assert direct(r,a,'C',L)==closed(r,a,'C',L); checks+=1
print('independent Catalan-vs-row formula checks:',checks)

# Check the G0 endpoint form, then screen its complement to the uniform band.
g0checks=Counter()
sc={k:Counter() for k in 'ABC'}
first={k:None for k in 'ABC'}
minimum={k:None for k in 'ABC'}
for r in range(1,7):
  for a in range(17):
    for k in 'ABC':
      if k=='A':
        words=[(uv[0],uv[1],p)
               for uv in combinations_with_replacement(range(2,9),2)
               for p in range(2,9)]
      elif k=='B':
        words=[(u,pq[0],pq[1])
               for pq in combinations_with_replacement(range(2,9),2)
               for u in range(2,9)]
      else:
        words=list(combinations_with_replacement(range(2,9),3))
      for L in words:
        v=closed(r,a,k,L); C=sc[k]; C['total']+=1
        g=g0_window(r,a,k,L)
        covered=bool(g) and uniform_band(r,a,k,g[2])
        C['g0_windows']+=bool(g); C['g0_band']+=covered
        C['remainder']+=not covered
        if g:
          Xv,Yv,Zv,x,sig,S,T=g; g0checks[k]+=1
          assert v==S+sig*T,(k,r,a,L,g,v)
        if not covered:
          C['remainder_negative']+=v<0
          C['remainder_zero']+=v==0
          if v>0 and first[k] is None: first[k]=(r,a,L,v)
          if minimum[k] is None or v<minimum[k][3]:
            minimum[k]=(r,a,L,v)
print('valid G0 endpoint checks:',dict(g0checks))
print('screen box: r=1..6, a=0..16, labels=2..8; covered=valid G0 window and uniform band')
for k in 'ABC':
  print(k,dict(sc[k]),
        'first positive remainder',first[k],
        'minimum remainder',minimum[k])
