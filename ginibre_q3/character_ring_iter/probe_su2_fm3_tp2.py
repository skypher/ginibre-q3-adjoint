#!/usr/bin/env python3
"""Coefficientwise TP_k test for exp(s A_p), A_p[j,k]=N_{j p}^k (SU(2) fusion with V_p).
usage: probe_su2_fm3_tp2.py P ORDER MAXIDX [KMINOR]"""
import sys, itertools
from fractions import Fraction
if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); sys.exit(0)
P=int(sys.argv[1]); ORD=int(sys.argv[2]); M=int(sys.argv[3]); KM=int(sys.argv[4]) if len(sys.argv)>4 else 2
N=M+ORD*P+2
def A(j,k): return 1 if (abs(j-P)<=k<=j+P and (j+k+P)%2==0) else 0
Amat=[[A(j,k) for k in range(N)] for j in range(N)]
# E[j][k] = list of coefficients (power series in s) of exp(sA)
# compute A^n restricted exactly: entries of A^n with j,k<M need indices < M+n*P
powers=[[[Fraction(int(j==k)) for k in range(N)] for j in range(N)]]
for n in range(1,ORD+1):
    prev=powers[-1]
    new=[[sum(prev[j][l]*Amat[l][k] for l in range(max(0,k-P),min(N,k+P+1))) for k in range(N)] for j in range(N)]
    powers.append(new)
fact=[1]
for n in range(1,ORD+1): fact.append(fact[-1]*n)
def ser(j,k): return [powers[n][j][k]/fact[n] for n in range(ORD+1)]
def smul(a,b): 
    c=[Fraction(0)]*(ORD+1)
    for i,x in enumerate(a):
        if x==0: continue
        for jj,y in enumerate(b[:ORD+1-i]): c[i+jj]+=x*y
    return c
def sadd(a,b,sg=1): return [x+sg*y for x,y in zip(a,b)]
def sdet(rows,cols):
    if len(rows)==1: return ser(rows[0],cols[0])
    tot=[Fraction(0)]*(ORD+1)
    for c in range(len(cols)):
        sub=sdet(rows[1:],cols[:c]+cols[c+1:])
        tot=sadd(tot,smul(ser(rows[0],cols[c]),sub),(-1)**c)
    return tot
bad=0; worst=None; count=0
for rows in itertools.combinations(range(M),KM):
    for cols in itertools.combinations(range(M),KM):
        d=sdet(list(rows),list(cols)); count+=1
        negs=[(n,v) for n,v in enumerate(d) if v<0]
        if negs:
            bad+=1
            if bad<=8: print("NEG minor rows=%s cols=%s first negative coeff s^%d = %s"%(rows,cols,negs[0][0],negs[0][1]))
print("p=%d order<=%d idx<%d k=%d minors=%d with negative coefficient: %d"%(P,ORD,M,KM,count,bad))
