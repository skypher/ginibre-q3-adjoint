#!/usr/bin/env python3
"""Exact all-exponent certificate for the rank-five SU(2)_9 orbit ring."""
import ast,base64,gzip,itertools,sys
from fractions import Fraction as Q
from pathlib import Path

class I:
 def __init__(s,a,b=None):s.a,s.b=Q(a),Q(a if b is None else b);assert s.a<=s.b
 def __add__(s,o):o=ii(o);return I(s.a+o.a,s.b+o.b)
 __radd__=__add__
 def __neg__(s):return I(-s.b,-s.a)
 def __sub__(s,o):return s+-ii(o)
 def __rsub__(s,o):return ii(o)-s
 def __mul__(s,o):
  o=ii(o);v=(s.a*o.a,s.a*o.b,s.b*o.a,s.b*o.b);return I(min(v),max(v))
 __rmul__=__mul__
 def inv(s):assert not(s.a<=0<=s.b);return I(1/s.b,1/s.a)
 def __truediv__(s,o):return s*ii(o).inv()
 def __pow__(s,n):
  if n<0:return s.inv()**-n
  if n==0:return I(1)
  if s.a>=0:return I(s.a**n,s.b**n)
  if s.b<=0:return I(s.a**n,s.b**n) if n&1 else I(s.b**n,s.a**n)
  return I(s.a**n,s.b**n) if n&1 else I(0,max((-s.a)**n,s.b**n))
 def ab(s):return s if s.a>=0 else (-s if s.b<=0 else I(0,max(-s.a,s.b)))
def ii(x):return x if isinstance(x,I) else I(x)

def trim(p):
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def der(p):return [Q(i)*p[i] for i in range(1,len(p))] or [Q(0)]
def pdm(a,b):
 a,b=trim(a[:]),trim(b[:]);q=[Q(0)]*max(1,len(a)-len(b)+1)
 while a!=[0] and len(a)>=len(b):
  d=len(a)-len(b);c=a[-1]/b[-1];q[d]=c
  for i,x in enumerate(b):a[d+i]-=c*x
  trim(a)
 return trim(q),trim(a)
def sturm(p):
 z=[trim(p[:]),trim(der(p))]
 while z[-1]!=[0]:
  _,r=pdm(z[-2],z[-1]);r=trim([-x for x in r])
  if r==[0]:break
  z.append(r)
 return z
def pe(p,x):
 y=Q(0)
 for c in reversed(p):y=y*x+c
 return y
def var(z,x):
 s=[1 if pe(p,x)>0 else -1 for p in z if pe(p,x)]
 return sum(a!=b for a,b in zip(s,s[1:]))
def dec(t):
 n=t.startswith('-');t=t.lstrip('-');a,b=t.split('.');x=Q(int(a+b),10**len(b));return -x if n else x
F=[Q(1),Q(3),Q(-3),Q(-4),Q(1),Q(1)]
C=('1.68250706566236233772362329784','0.830830026003772851058548298459','-0.284629676546570280887585337233','-1.30972146789057012811385014493','-1.91898594722899477978073611413')
def roots():
 z=sturm(F);r=[];e=Q(1,10**29)
 for t in C:
  x=dec(t);v=I(x-e,x+e);assert var(z,v.a)-var(z,v.b)==1 and pe(F,v.a)*pe(F,v.b)<0;r.append(v)
 assert var(z,Q(-3))-var(z,Q(3))==5
 return r
def chars(x):
 x2=x*x;x3=x2*x;return[x+1,x2+x-1,x3+x2-2*x-1,x2*x2+x3-3*x2-2*x+1]
def signs(code):
 s=[]
 for _ in range(4):d=code%3;code//=3;s.append(0 if d==0 else(1 if d==1 else -1))
 return tuple(s)

def chamber(sc,pc,R):
 sg=signs(sc);ac=tuple(i for i,x in enumerate(sg) if x)
 if not ac or not any(x<0 for x in sg) or pc>=1<<len(ac):return None
 pw=[0]*4;mp=0
 for j,l in enumerate(ac):
  pw[l]=1 if pc>>j&1 else 2
  if sg[l]<0 and pw[l]&1:mp^=1
 if mp:return None
 V=[chars(x) for x in R];W=[(2-x)/11 for x in R];T=[]
 for i in range(5):
  for j in range(i+1,5):
   c=2*W[i]*W[j];si=1;la=[I(1)]*4
   for l in ac:
    b=V[i][l]+sg[l]*V[j][l];assert b.a>0 or b.b<0
    if b.b<0 and pw[l]&1:si=-si
    c*=b.ab()**pw[l];la[l]=b**2
   T.append((si,c,tuple(la)))
 ne=tuple(i for i,t in enumerate(T) if t[0]<0);po=tuple(i for i,t in enumerate(T) if t[0]>0)
 return sg,ac,tuple(pw),tuple(T),ne,po
def mag(t,f):
 x=t[1]
 for i,n in enumerate(f):x*=t[2][i]**n
 return x
def base(t,c):
 x=I(1)
 for i in c:x*=t[2][i]
 return x
def ge(a,b):return a.a>=b.b

def hall(ch,f,coords):
 T,ne,po=ch[3:];M=[mag(t,f) for t in T];E=[]
 for n in ne:E.append([k for k,p in enumerate(po) if all(ge(base(T[p],c),base(T[n],c)) for c in coords)])
 for mask in range(1,1<<len(ne)):
  d=I(0);nb=set()
  for k,n in enumerate(ne):
   if mask>>k&1:d+=M[n];nb.update(E[k])
  c=I(0)
  for k in nb:c+=M[po[k]]
  if not ge(c,d):return False
 return True
def braid(ch,f,v):
 for p in itertools.permutations(v):
  c=[];z=[]
  for x in p:z.append(x);c.append(z[:])
  if not hall(ch,f,c):return False
 return True
def outs(a,b):return range(abs(a-b),min(a+b,9-a-b)+1)
def corner(ch,f):
 st={(0,0):1}
 for l in ch[1]:
  for _ in range(ch[2][l]+2*f[l]):
   ns={}
   for (a,b),v in st.items():
    for x in outs(a,l+1):ns[x,b]=ns.get((x,b),0)+v
    for x in outs(b,l+1):ns[a,x]=ns.get((a,x),0)+ch[0][l]*v
   st={k:v for k,v in ns.items() if v}
 return st.get((0,0),0)
def residual(ch,rc):
 f=[0]*4;v=[]
 for j,l in enumerate(ch[1]):
  if rc>>j&1:f[l]=1;v.append(l)
 return f,v
def sep(ch):return ch[0]==(-1,0,0,-1) and ch[2]==(1,0,0,1)
def frac(x):
 return I(Q(1,4) if abs(x-.25)<1e-9 else Q(1,2) if abs(x-.5)<1e-9 else Q(3,4))

def endpoint(ch,f,cert):
 hard,lo,hi,fr,A,B,C,easy,V,F0=cert
 if tuple(f)!=tuple(F0):return False
 T,ne,po=ch[3:];M=[mag(t,f) for t in T];V=tuple(V);mid=[x for x in V if x not in(0,3)];fr=frac(fr)
 al=T[lo][2][0]/T[hard][2][0];bl=T[lo][2][3]/T[hard][2][3];ah=T[hi][2][0]/T[hard][2][0];bh=T[hi][2][3]/T[hard][2][3]
 if not(al.a>1 and bl.b<1 and ah.b<1 and bh.a>1):return False
 if any(not ge(T[lo][2][x],T[hard][2][x]) or not ge(T[hi][2][x],T[hard][2][x]) for x in mid):return False
 cap=[]
 for p in po:cap.append(I(0) if p==hi else((I(1)-fr)*M[p] if p==lo else M[p]))
 E=[]
 for n in easy:E.append([k for k,p in enumerate(po) if cap[k].b>0 and all(ge(T[p][2][x],T[n][2][x]) for x in V)])
 for mask in range(1,1<<len(easy)):
  d=I(0);nb=set()
  for k,n in enumerate(easy):
   if mask>>k&1:d+=M[n];nb.update(E[k])
  c=I(0)
  for k in nb:c+=cap[k]
  if not ge(c,d):return False
 lc=fr*M[lo]/M[hard];hc=M[hi]/M[hard]
 return all(x.a>=1 for x in(al**B*bl**A,ah**B*bh**A,lc**B*bl**C,hc**B*bh**(C+1)))
def rec(sc,pc,ch,cert,f=None,v=None):
 m,d=cert
 if m=='leaf':return corner(ch,f or[0]*4)>=0
 if m=='direct':V,F0=d;return hall(ch,f or list(F0),[[x] for x in(v or V)])
 if m=='braid':V,F0=d;return braid(ch,f or list(F0),v or V)
 if m=='shift':
  V,F0,s,tail,faces=d
  if f is not None and tuple(f)!=tuple(F0):return False
  ff=list(F0)
  for x in V:ff[x]+=s
  if tail[0]!='braid' or not braid(ch,ff,V):return False
  for x,o,c in faces:
   f2=list(F0);f2[x]+=o
   if not rec(sc,pc,ch,c,f2,[y for y in V if y!=x]):return False
  return True
 raise AssertionError(m)

def main():
 R=roots();path=Path(__file__).resolve().parents[1]/'certificates'/'su2_odd_orbit_rank5_ledger.gz.b64'
 L=ast.literal_eval(gzip.decompress(base64.b64decode(path.read_bytes())).decode());D={};
 for sc,pc,rc,m,c in L:
  k=(sc,pc,rc);assert k not in D;D[k]=(m,c)
 E=set();chs=pts=ptr=0
 for sc in range(81):
  for pc in range(16):
   ch=chamber(sc,pc,R)
   if ch is None:continue
   chs+=1;n=1<<len(ch[1])
   if not ch[4]:pts+=1;ptr+=n
   else:
    for rc in range(n):E.add((sc,pc,rc))
 assert set(D)==E
 S={x:0 for x in('direct','endpoint','leaf','shift','separated')};mn=None;zl=0
 for k in sorted(E):
  sc,pc,rc=k;ch=chamber(sc,pc,R);f,v=residual(ch,rc);m,c=D[k]
  if m=='separated':assert sep(ch)
  elif m=='endpoint':assert endpoint(ch,f,c)
  elif m=='leaf':
   z=corner(ch,f);assert z>=0;mn=z if mn is None else min(mn,z);zl+=z==0
  else:assert rec(sc,pc,ch,c,f,v)
  S[m]+=1
 assert S=={'direct':1828,'endpoint':332,'leaf':205,'shift':39,'separated':4}
 assert(chs,pts,ptr)==(272,66,560)
 print('SU2_ODD_ORBIT_RANK5_EXACT PASS chambers=272 pointwise_chambers=66 pointwise_regimes=560 certified_regimes=2408 direct=1828 endpoint=332 leaf=205 shift=39 separated=4 minimum_leaf=%s zero_leaves=%s'%(mn,zl))
if __name__=='__main__':main()
