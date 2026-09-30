import itertools
from fractions import Fraction as Fr
from math import comb
src=open("qs_family_screen.py").read().split("import random")[0]
exec(src)
def bernstein(coeffs):  # coeffs[k] of s^k, degree d -> b_j with sum b_j C(d,j) s^j (1-s)^(d-j)
    d=len(coeffs)-1
    return [sum(Fr(comb(j,k),comb(d,k))*coeffs[k] for k in range(j+1)) for j in range(d+1)]
lists=[]
for L in range(2,7):
    for lab in itertools.combinations_with_replacement(range(1,4),L):
        if sum(lab)%2==0 and sum(lab)<=12: lists.append(lab)
tot=0; negB=0; nonmono=0; ex=[]
for lab in lists:
    L=len(lab)
    for T in range(1<<L):
        if bin(T).count('1')%2: continue
        p=Fqs(lab,T)
        c={}
        for (i,j),v in p.items():
            if i==0: c[j]=c.get(j,0)+v
        if not c: continue
        d=max(c); co=[c.get(k,0) for k in range(d+1)]
        tot+=1
        b=bernstein(co)
        if min(b)<0:
            negB+=1
            if len(ex)<6: ex.append((lab,T,co,[str(x) for x in b]))
        # monotone decreasing on [0,1]? derivative coefficients
        vals=[sum(co[k]*Fr(t,20)**k for k in range(d+1)) for t in range(21)]
        if any(vals[i+1]>vals[i] for i in range(20)): nonmono+=1
print("profiles",tot,"negative Bernstein",negB,"non-decreasing somewhere",nonmono)
for e in ex: print(e)
