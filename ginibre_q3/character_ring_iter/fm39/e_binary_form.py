# (E)_(q, sign) at (j, j+q+1) as a binary quadratic form in (c_j, c_(j+1)) via the three-term recurrence (exact rationals).
# Definite (disc < 0, leading > 0) => (E) holds for any ratio: a ratio-free proof, as in step 1 of Theorem OL.
from fractions import Fraction as Fr
from collections import Counter
import sys
def window_vals(d,N,j,q,p0,p1):
    c={j:Fr(p0), j+1:Fr(p1)}
    for k in range(j+1, j+q+2):
        c[k+1]=(d*c[k]-(N-k+1)*c[k-1])/(k+1)
    for k in (j, j-1):
        c[k-1]=(d*c[k]-(k+1)*c[k+1])/(N-k+1)
    return c
def Evalue(d,N,j,q,sg,p0,p1):
    c=window_vals(d,N,j,q,p0,p1); i=j+q+1
    C=lambda k: c[k]
    D=lambda k: C(k)**2-C(k-1)*C(k+1); B=lambda k: C(k-1)+C(k+1)
    W=B(i)*C(j)-C(i)*B(j)
    return D(j)-D(i)-sg*W
def form(d,N,j,q,sg):
    f10=Evalue(d,N,j,q,sg,1,0); f01=Evalue(d,N,j,q,sg,0,1); f11=Evalue(d,N,j,q,sg,1,1)
    A=f10; C=f01; B=f11-A-C
    return A,B,C
qmax=int(sys.argv[1]) if len(sys.argv)>1 else 6
st=Counter()
for a in range(3,41):
    for e in range(3,41):
        if abs(a-e)<2: continue
        N=a+e; d=a-e
        for j in range((N+1)//2, N+1):
            for q in range(1,qmax+1):
                i=j+q+1
                if i>N+1: continue
                for sg in (1,-1):
                    A,B,C=form(d,N,j,q,sg)
                    key=(q,'-W' if sg==1 else '+W')
                    st[key+('pairs',)]+=1
                    if B*B-4*A*C<0 and A>0: st[key+('definite',)]+=1
for q in range(1,qmax+1):
    for s in ('-W','+W'):
        n=st[(q,s,'pairs')]; dfn=st[(q,s,'definite')]
        print('q=%d %s: open-region pairs %6d  definite %6d  (%.1f%%)'%(q,s,n,dfn,100*dfn/max(n,1)))
