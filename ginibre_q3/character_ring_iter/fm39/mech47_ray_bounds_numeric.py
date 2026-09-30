# Main agent: numerical (scipy quad) check of the FM-MECH47 ray lower bounds (4), (5) on random parameters.
import random, math
from scipy.integrate import quad
th=4*math.sqrt(2)/9
def I(q,N,b,al,ga):
    beta=4*(1+q*q); lp=4*(1+q*q+q); lm=4*(1+q*q-q); B=beta-2; Bp=lp-2; Bm=lm-2
    f=lambda z: z**N*math.sqrt((1-z)*(1-q*q*z))*((beta*z-2)/B)**b*((lp*z-2)/Bp)**al*((lm*z-2)/Bm)**ga
    return quad(f,0,1,limit=400,points=[2/3,2/5,2/9])[0]
def bound4(N,L,b):
    K=L+N/3+1
    return 1/(9*(K+b)*(K+b+1)) - 4*th**(L-2)*(2/3)**(N-L/2)*(1/(4*(b+1)*(b+2))+(32/225)*(5/9)**b)
def bound5(N,L,b):
    K=L+N/3+1
    t1=math.factorial(N)/2**(N+1)
    for j in range(1,N+2): t1/=(b+j)
    t2=((2/3)**(N+1)-(2/5)**(N+1))*(5/9)**b/(N+1)
    return 1/(9*(K+b)*(K+b+1)) - 2**L*(t1+t2)
random.seed(3); viol=0; n=0
for _ in range(3000):
    N=random.randint(1,30); L=random.randint(2,min(24,2*N) if N>=1 else 2); al=random.randint(0,L); ga=L-al; b=random.randint(0,50); q=random.uniform(-1,1)
    if N<L/2: continue
    v=I(q,N,b,al,ga); n+=1
    for bd in (bound4(N,L,b) if N>=L/2 else -1e9, bound5(N,L,b)):
        if v < bd - 1e-12*max(1,abs(bd)): viol+=1; print("VIOL",q,N,L,al,b,v,bd)
print("ray samples",n,"bound violations",viol)
