# FM-STR12e (main agent): Hermite (q = 1) analogue of HPP.  Expected output: 71,508 splits, no negative Phi,
# no failure for the height direction (1,1); failures for (1,0), (1,-1), (1,2).  Pairs allowed, labels <= 4, length 3..7.
# Hermite (q=1) analogue: F_X = prod (He_n(x) + eps He_n(y)), x,y iid N(0,1); coefficients in He_r(x)He_s(y);
# Pi_T = sum_{r+s<=T} r! s! f_A(r,s) f_B(r,s).  Test every split, several directions.
import itertools
from math import comb, factorial
from collections import defaultdict
def hmul(m,n):  # He_m He_n = sum_k C(m,k)C(n,k)k! He_{m+n-2k}
    return [(m+n-2*k, comb(m,k)*comb(n,k)*factorial(k)) for k in range(min(m,n)+1)]
def table(word):
    d={(0,0):1}
    for z in word:
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for c,w in hmul(a,n): q[c,b]+=v*w
            for c,w in hmul(b,n): q[a,c]+=e*v*w
        d={k:v for k,v in q.items() if v}
    return d
labs=[s*n for n in range(1,5) for s in (1,-1)]
dirs=[(1,1),(1,0),(1,-1),(1,2)]
bad=defaultdict(int); n=0; negphi=0
for N in range(3,8):
    for w in itertools.combinations_with_replacement(sorted(labs),N):
        if sum(1 for z in w if z<0)%2 or sum(map(abs,w))%2: continue
        for m in range(1,1<<(N-1)):
            A=[w[i] for i in range(N) if m>>i&1]; B=[w[i] for i in range(N) if not m>>i&1]
            fa=table(A); fb=table(B); P={k:v*fb[k]*factorial(k[0])*factorial(k[1]) for k,v in fa.items() if k in fb}
            n+=1
            if sum(P.values())<0: negphi+=1
            for d in dirs:
                lay=defaultdict(int)
                for (r,s),v in P.items(): lay[d[0]*r+d[1]*s]+=v
                run=0
                for t in sorted(lay):
                    run+=lay[t]
                    if run<0: bad[d]+=1; break
print("Hermite splits",n,"negative Phi",negphi,"prefix failures by direction",dict(bad))
