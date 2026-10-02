cd /home/yang/q3adjoint
python3 -u - <<'FM188LAURENT'
from collections import defaultdict
from fractions import Fraction
from decimal import Decimal, getcontext
from math import comb
from datetime import datetime, timezone

def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
def wt(w):
    return sum(abs(z) for z in w)
def laurent_coeff(word):
    poly={(0,0):1}
    for z in word:
        n=abs(z);eps=1 if z>0 else -1
        terms=[(k,0,1) for k in range(-n,n+1,2)]
        terms += [(0,k,eps) for k in range(-n,n+1,2)]
        nxt=defaultdict(int)
        for (i,j),a in poly.items():
            for x,y,v in terms:nxt[(i+x,j+y)]+=a*v
        poly={k:v for k,v in nxt.items() if v}
    W=wt(word);c=[]
    for r in range(W//2+1):
        c.append(poly.get((r,r),0)-poly.get((r+2,r),0)
                 -poly.get((r,r+2),0)+poly.get((r+2,r+2),0))
    while len(c)>1 and c[-1]==0:c.pop()
    return c
def bernstein(c,lo,hi):
    d=len(c)-1;width=hi-lo
    a=[Fraction(0) for _ in range(d+1)]
    for k,ck in enumerate(c):
        for j in range(k+1):
            a[j]+=ck*comb(k,j)*lo**(k-j)*width**j
    return [sum(a[j]*Fraction(comb(i,j),comb(d,j))
                for j in range(i+1)) for i in range(d+1)]
def certify(c,lo=Fraction(-1),hi=Fraction(1),depth=0,maxdepth=18):
    b=bernstein(c,lo,hi)
    if min(b)>=0:return True,depth
    if b[0]<0 or b[-1]<0:return False,depth
    if depth>=maxdepth:return None,depth
    mid=(lo+hi)/2
    left,l=certify(c,lo,mid,depth+1,maxdepth)
    right,r=certify(c,mid,hi,depth+1,maxdepth)
    if left is False or right is False:return False,max(l,r)
    if left is None or right is None:return None,max(l,r)
    return True,max(l,r)

w=[50]+[-3]*14+[2]*3+[1]*2
c=laurent_coeff(w)
rho=Fraction(1048575,1048576)
D=sum(Fraction(a)*rho**i for i,a in enumerate(c))
norm=sum(abs(a) for a in c)
ratio=D/norm
getcontext().prec=30
print(stamp(),"word",w,"W",wt(w),"degree",len(c)-1)
print("coeff",c)
print("D1",sum(c),"Dminus1",sum((-1)**i*a for i,a in enumerate(c)))
print("rho",rho,"D",D,"L1",norm,"normalized",ratio)
print("normalized_decimal",Decimal(ratio.numerator)/Decimal(ratio.denominator))
print("Bernstein_subdivision",certify([Fraction(a) for a in c]))
assert ratio==Fraction(17228238974574007225045386509025693319241962915074523713429035210036191924664664215786865301496201368911921575247890907929817543575130322217653499982084856271533073915947954355267938503622523844911123622424927189256370084049,20358663536539079098020075932140034400666691422405138541431257379036669971373786327370769197498221304424298785433509789678555760428153966828370697391927292127475678772166038872300041749618505169568112946090788302468805688408057658362414533497546419537685273307604430118554816680139986334990690769331814400)
FM188LAURENT
