# Independent checks of FM-MECH6 (one-label monotonicity D_k >= D_(k+1) for all t).
import time, math
from fractions import Fraction as F
from math import comb, isqrt
t0=time.time()
# 1) coverage (Lemma 5): folded N>=2t, t>=4, 1<=p<(N-2)/3, p==N mod 2: F<=0 or d^2>=4(n-1)(k+2) or N>=3ts^2
bad=0; pts=0
for t in range(4,61):
    for N in range(2*t,3001):
        d=N-2*t
        for p in range(1+((N-1)%2==0 and 0 or 0),N):
            if (N-p)%2 or not (3*p<N-2): continue
            s=p+1; n=(N-p)//2; k=(N+p)//2; pts+=1
            Fv=s*s*d*d-4*p*(p+2)*(n+1)*(k+2)
            if Fv<=0 or d*d>=4*(n-1)*(k+2) or N>=3*t*s*s: continue
            bad+=1
            if bad<5: print('UNCOVERED',t,N,p)
print(f'1) coverage: {pts} middle-range points (t<=60, N<=3000), uncovered {bad} ({time.time()-t0:.0f}s)',flush=True)
# 2) reciprocal-root bounds (Lemma 3), exact, t<=80, N in [2t, 2t+400]
def lowcoef(t,N):  # coefficients of x^0..x^3 of monic Krawtchouk K_t(x;N)
    o=[1,0,0,0]; c=[0,1,0,0]
    if t==0: return o
    for j in range(1,t):
        aj=j*(N-j+1); nw=[(c[i-1] if i else 0)-aj*o[i] for i in range(4)]; o,c=c,nw
    return c
bad2=0; cnt2=0
for t in range(4,81):
    for N in range(2*t,2*t+401):
        L=lowcoef(t,N)
        if t%2==0: act=-F(L[2],L[0]); bnd=F(t,2*(N-t+2))
        else: act=-F(L[3],L[1]); bnd=F(t-1,4*(N-t+3))
        cnt2+=1
        if not (0<act<=bnd): bad2+=1
print(f'2) reciprocal-root bounds: {cnt2} (t,N), failures {bad2} ({time.time()-t0:.0f}s)',flush=True)
# 3) theorem + intermediate lemma claims, exact, t<=30, N<=400
def cvec(a,t):
    # coefficients of (1+z)^a (1-z)^t
    N=a+t; c=[0]*(N+1)
    for i in range(t+1):
        ci=(-1)**i*comb(t,i)
        for j in range(a+1): c[i+j]+=ci*comb(a,j)
    return c
bad3=0; cnt3=0; cls={'band':0,'outer':0,'central':0}
for t in range(4,31):
    for N in range(2*t,401):
        a=N-t; c=cvec(a,t); cc=lambda j: c[j] if 0<=j<=N else 0
        d=N-2*t
        for p in range(1,N):
            if (N-p)%2 or not (3*p<N-2): continue
            s=p+1; n=(N-p)//2; k=(N+p)//2
            Dk=cc(k)**2-cc(k-1)*cc(k+1); Dk1=cc(k+1)**2-cc(k)*cc(k+2)
            cnt3+=1
            if Dk<Dk1: bad3+=1; print('THEOREM FAILS',t,N,p)
            Fv=s*s*d*d-4*p*(p+2)*(n+1)*(k+2)
            if Fv<=0: cls['band']+=1; continue
            r=F(cc(k+1),cc(k))
            lam=(n-1)*(k+2)
            if d*d>=4*lam:
                cls['outer']+=1
                # r>0 and lam r^2 - n d r + n^2 >= 0 with r <= smaller root (r <= n d/(2 lam))
                if not (r>0 and lam*r*r-n*d*r+n*n>=0 and r<=F(n*d,2*lam)): bad3+=1; print('OUTER claim fails',t,N,p)
            else:
                cls['central']+=1
                assert N>=3*t*s*s
                if t%2==0:
                    if not (0<r<=1): bad3+=1; print('CENTRAL even claim fails',t,N,p)
                else:
                    U=F(s+1,s-1); b=F(4*t+2,3)
                    if not (r>=U*(1-b*F(s,N))): bad3+=1; print('CENTRAL odd claim fails',t,N,p)
print(f'3) theorem + lemma claims at {cnt3} points (t<=30, N<=400): failures {bad3}; classes {cls} ({time.time()-t0:.0f}s)')
