import argparse, ast
from pathlib import Path
from fractions import Fraction as Q
from math import comb, factorial as fac
from functools import lru_cache

argparse.ArgumentParser(description="ADV-1b exact SU(2) midpoint test").parse_args()
path=Path("ginibre_q3/character_ring_iter/fm39/mech54_genfun_obstructions_repro.py")
lib={}
exec(compile(ast.Module(body=[n for n in ast.parse(path.read_text()).body
    if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))],
    type_ignores=[]),str(path),"exec"),lib)

@lru_cache(None)
def beta(m,g):
    return Q(2*fac(2*m)*fac(2*g+2),
             4**(m+g+1)*fac(m)*fac(g+1)*fac(m+g+1))

def R(e,a,b,p,l,j):
    g=(e+l)//2+j
    return sum(Q((-1)**k*fac(p-k)*2**(p-l-2*k),
                 fac(l)*fac(k)*fac(p-l-2*k))
               *comb(b-j,h)*4**h*(-1)**(b-j-h)
               *beta((a+p-l-2*k+2*h)//2,g)
               for k in range((p-l)//2+1) for h in range(b-j+1))

def I(e,l,j):
    return sum(Q((-1)**k*fac(2*l-2*k),
                 2**l*fac(k)*fac(l-k)*fac(l-2*k))
               *Q(comb(j,h)*3**h*(-1)**(j-h),
                   2**j*(e+l-2*k+2*h+1))
               for k in range(l//2+1) for h in range(j+1))

assert R(6,8,3,8,2,1)==Q(-77,131072)
assert R(6,10,3,8,2,1)==Q(261,262144)
assert R(3,11,4,8,1,1)==Q(-21,32768)
assert R(3,13,4,8,1,1)==Q(741,262144)

for e,a,b,d in ((3,11,3,6),(6,8,3,6),(3,11,4,7),(6,8,4,7)):
    p=e+a+2*b-2*d
    groups=[Q(0)]*(b+1)
    negative=0
    for j in range(b+1):
        for l in range(e%2,min(p,e+2*j)+1,2):
            A=Q(4**l*fac(l)**2*(2*l+1)*fac(p-l),fac(p+l+1))
            term=comb(b,j)*8**j*A*I(e,l,j)*R(e,a,b,p,l,j)*R(e,a+2,b,p,l,j)
            groups[j]+=term
            negative+=term<0
    value=4**(e+a+1)*Q(2,3)**b*sum(groups)
    assert value==lib["profile"](e,a,b)[p]
    assert min(groups)>=0
    print((e,a,b,d),"F =",value,"negative channels =",negative)
