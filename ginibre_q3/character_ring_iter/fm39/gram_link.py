# (1) link the Astra R_t factorization to direct phi_r values (Catalan moments), both parities.
import sympy as S
from fractions import Fraction as F
from math import comb, factorial
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
src=open('gram_repro.py').read().split("def cx(")[0]
G={}
exec(src,G)
Rpoly=G['Rpoly']; Nsym=G['N']; Tsym=G['T']
from functools import lru_cache
Rpoly=lru_cache(None)(Rpoly)
def formula(t,NN,kk):
    R=Rpoly(t); pp,nn=2*kk-NN,NN-kk
    fall=factorial(NN)//factorial(NN-t)
    return F(comb(NN,kk)**2*(pp+1)*(NN-t+1)*(NN-t+2)*int(R.eval({Nsym:NN,Tsym:pp*(pp+2)})),(nn+1)*(kk+1)**2*(kk+2)*fall**2)
Rc={}
bad=0; cnt=0
for t in range(4,11):
    for NN in range(t,t+14):
        a_=NN-t
        if t%2==0:   # hat S_q, t=2r, k=(N+q)/2
            r_=t//2
            for q in range(2,NN+1):
                if (NN+q)%2: continue
                kk=(NN+q)//2; direct=phi(r_,shat(q),a_); f=formula(t,NN,kk); cnt+=1
                if direct!=f: bad+=1; print('MISMATCH hatS',t,NN,q,direct,f) if bad<4 else None
        else:        # Sym^s = h_s, t=2r-1, k=(N+s+1)/2
            r_=(t+1)//2
            for s_ in range(2,NN):
                if (NN+s_+1)%2: continue
                kk=(NN+s_+1)//2
                if kk>NN: continue
                direct=phi(r_,h(s_),a_); f=formula(t,NN,kk); cnt+=1
                if direct!=f: bad+=1; print('MISMATCH Sym',t,NN,s_,direct,f) if bad<4 else None
    print(f't={t}: cumulative checks {cnt}, mismatches {bad}',flush=True)
