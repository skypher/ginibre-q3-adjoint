# Check: for P=(1+z)^a(1-z)^e prod(1+t z+z^2), the quadratic W_c(p,q) equals E[K(x,y) U_p(x) U_q(y)],
# K = (x+y)^a (x-y)^e prod (x^2 + t x y + y^2 + t^2 - 4), x,y independent semicircle on [-2,2].
import random
from fractions import Fraction as Fr
from math import comb
from collections import defaultdict
exec(open('wformula_check.py').read().split("n=0")[0])
exec(open('recip_test.py').read().split("def tools")[0].split("import random")[1].replace("from fractions import Fraction as Fr\nfrom math import comb\n",""))
def pmul(f,g):
    o=defaultdict(Fr)
    for (i,j),x in f.items():
        for (k,l),y in g.items(): o[i+k,j+l]+=x*y
    return dict(o)
def Kpoly(a,e,ts):
    K={(0,0):Fr(1)}
    for _ in range(a): K=pmul(K,{(1,0):1,(0,1):1})
    for _ in range(e): K=pmul(K,{(1,0):1,(0,1):-1})
    for t in ts: K=pmul(K,{(2,0):1,(1,1):t,(0,2):1,(0,0):t*t-4})
    return K
def Elin(K,p,q): return sum(v*_mom(i,p)*_mom(j,q) for (i,j),v in K.items())
random.seed(3); n=0
for _ in range(200):
    a=random.randint(0,4); e=random.randint(0,4); ts=[Fr(random.randint(-30,30),7) for _ in range(random.randint(1,2))]
    c=mkrow(a,e,ts); N=len(c)-1; K=Kpoly(a,e,ts)
    for p in range(0,N+2):
        for q in range(0,N+2):
            if (N+p+q)%2: continue
            assert Wq(c,p,q,e%2)==Elin(K,p,q),(a,e,ts,p,q,Wq(c,p,q,e%2),Elin(K,p,q)); n+=1
print('identity W_c(p,q) = E[K U_p U_q] holds in',n,'cases (any real t, not only |t|>=2)')
