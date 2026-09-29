#!/usr/bin/env python3
"""Check Conjecture qFM3_2 (FM3 note, 'Euler-characteristic mechanism and the graded r = 2 conjecture').

Phi_kappa(t) = sum_{l(lambda)<=4} K_{lambda,kappa}(t) Ehat_lambda(t) must lie in N[t], where K is the
charge Kostka-Foulkes polynomial (computed by kf_charge_sp4.cpp) and Ehat is the closed form below.
At t=1, Phi_kappa(1) = L(kappa) = 4I(kappa+{2}) + 8I(kappa) - 3I(kappa+{1,1}).

usage: verify_qfm3_r2.py Nmax        (checks all partitions of even size 2..Nmax; prints per-size summary)
"""
import os, subprocess, sys, tempfile
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
def build():
    exe=os.path.join(tempfile.gettempdir(),'kf_charge_sp4_%d'%os.getuid())
    subprocess.run(['g++','-O2','-o',exe,os.path.join(HERE,'kf_charge_sp4.cpp')],check=True)
    return exe
def Ehat(l):
    a,b,c=l[0]-l[1],l[1]-l[2],l[2]-l[3]
    if a==0 and c==0: return {0:1,1:1,2:1,3:1,4:1} if b==0 else {0:1,4:1}
    if (a,c) in [(2,0),(0,2)]: return {2:1}
    if (a,c)==(1,1): return {1:-1,3:-1}
    return {}
def phi(exe,kap):
    out=subprocess.run([exe,",".join(map(str,kap))],capture_output=True,text=True,check=True).stdout
    P=defaultdict(int)
    for line in out.strip().split('\n'):
        if not line: continue
        sh,rest=line.split(':',1); l=tuple(map(int,sh.split(',')))
        E=Ehat(l)
        for tok in rest.split():
            d,c=map(int,tok.split(':'))
            for e,v in E.items(): P[d+e]+=c*v
    return {k:v for k,v in P.items() if v}
def partitions(n,maxpart=None):
    if maxpart is None: maxpart=n
    if n==0: yield (); return
    for k in range(min(n,maxpart),0,-1):
        for rest in partitions(n-k,k): yield (k,)+rest
def main():
    if len(sys.argv)<2 or sys.argv[1] in ('-h','--help'): print(__doc__); return
    Nmax=int(sys.argv[1]); exe=build()
    # boundary values
    assert phi(exe,(2,))=={2:1}
    assert phi(exe,(1,1))=={0:1,3:1,4:1}
    assert phi(exe,(2,1,1))=={4:1,5:1}
    total=0; bad=0
    for N in range(2,Nmax+1,2):
        cnt=0; neg=[]
        for kap in partitions(N):
            P=phi(exe,kap); cnt+=1
            if any(v<0 for v in P.values()): neg.append(kap)
        total+=cnt; bad+=len(neg)
        print(f'N={N}: {cnt} partitions, negative coefficients in {len(neg)} {neg[:3]}',flush=True)
    print('TOTAL',total,'partitions; failures',bad); sys.exit(1 if bad else 0)
if __name__=='__main__': main()
