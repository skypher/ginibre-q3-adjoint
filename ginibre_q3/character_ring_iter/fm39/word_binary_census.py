# Binary-form certificates for arbitrary consumer words: phi(word) as a quadratic form in (c_x0, c_(x0+1)) via the recurrence.
# PSD form (A >= 0, C >= 0, disc <= 0) certifies phi >= 0.  Exact rationals, heartbeat.
import sys, time, itertools
from fractions import Fraction as Fr
from collections import Counter
exec(open('general_row_words.py').read().split("random.seed(")[0])
def rec_row(a,e,x0,p0,p1):
    N=a+e; d=a-e; c={x0:Fr(p0), x0+1:Fr(p1)}
    for k in range(x0+1,N): c[k+1]=(d*c[k]-(N-k+1)*c[k-1])/(k+1)
    for k in range(x0,0,-1): c[k-1]=(d*c[k]-(k+1)*c[k+1])/(N-k+1)
    return [c[k] for k in range(N+1)]
def word_form(a,e,Ls,Ps):
    N=a+e; x0=max(0,min(N-1,N//2))
    vals=[phi_row(rec_row(a,e,x0,p0,p1),Ls,Ps,e%2) for (p0,p1) in ((1,0),(0,1),(1,1))]
    A,C=vals[0],vals[1]; B=vals[2]-A-C
    return A,B,C
pattern=sys.argv[1]; R2=int(sys.argv[2]); A2=int(sys.argv[3]); LMAX=int(sys.argv[4])
nh={'hhS':2,'hSS':1,'SSS':0}[pattern]; ns=3-nh
st=Counter(); bad=[]; t0=time.time(); last=t0
for r in range(1,R2+1):
    e=2*r-nh
    if e<0: continue
    for a in range(0,A2+1):
        N=a+e
        if N<1: continue
        for hl in itertools.combinations_with_replacement(range(2,LMAX+1),nh):
            for sl in itertools.combinations_with_replacement(range(2,LMAX+1),ns):
                Ls=[u+1 for u in hl]; Ps=list(sl)
                if (N+sum(Ls)+sum(Ps))%2: continue
                A,B,C=word_form(a,e,Ls,Ps); st['words']+=1
                disc=B*B-4*A*C
                if disc<0 and A>0: st['definite']+=1
                elif A>=0 and C>=0 and disc<=0: st['psd']+=1
                else:
                    st['not psd']+=1
                    if len(bad)<8: bad.append((r,a,hl,sl))
                if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat',pattern,'r=%d a=%d'%(r,a),dict(st),flush=True)
print(time.strftime('%H:%M:%S'),'done',pattern,dict(st),'elapsed',round(time.time()-t0,1)); print('not-psd examples:',bad)
