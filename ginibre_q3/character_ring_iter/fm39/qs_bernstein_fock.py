import random, time, sys
from fractions import Fraction as Fr
from math import comb
from qs_fock import F
def interp(xs,ys):
    # exact Lagrange -> monomial coefficients
    n=len(xs); coef=[Fr(0)]*n
    for i in range(n):
        num=[Fr(1)]; den=Fr(1)
        for j in range(n):
            if j==i: continue
            num=[ (num[k-1] if k>0 else 0) - xs[j]*(num[k] if k<len(num) else 0) for k in range(len(num)+1)]
            den*=xs[i]-xs[j]
        for k in range(n): coef[k]+=ys[i]*num[k]/den
    while len(coef)>1 and coef[-1]==0: coef.pop()
    return coef
def bernstein(co):
    d=len(co)-1
    return [sum(Fr(comb(j,k),comb(d,k))*co[k] for k in range(j+1)) for j in range(d+1)]
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 5)
t0=time.time(); tot=0; neg=0; done=0
lists=[]
while len(lists)<60:
    L=rng.randint(5,8); lab=tuple(sorted(rng.randint(1,4) for _ in range(L)))
    if sum(lab)%2 or sum(lab)>16 or lab in lists: continue
    lists.append(lab)
for lab in lists:
    D=sum(lab); dmax=(D//4)**2+1
    xs=[Fr(k,dmax) for k in range(dmax+1)]
    L=len(lab)
    for T in range(1<<L):
        if bin(T).count('1')%2: continue
        ys=[F(lab,T,Fr(0),x) for x in xs]
        co=interp(xs,ys)
        b=bernstein(co); tot+=1
        if min(b)<0:
            neg+=1; print("NEG",lab,T,[str(c) for c in co],flush=True)
    done+=1
    print(time.strftime("%H:%M:%S"),"lists",done,"/",len(lists),"profiles",tot,"negBern",neg,"last",lab,"elapsed %.0fs"%(time.time()-t0),flush=True)
print("DONE profiles",tot,"negative Bernstein",neg)
