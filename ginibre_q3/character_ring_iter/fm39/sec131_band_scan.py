# Main agent: the band R minus OL-outer for the q=2 plus sign; distribution and central-condition margin.
import math
def polys(n,k,D):
    j=n+k
    H=(n+1)*(j+2)**2*(j+3)**2*(j+4)
    alpha=D*n*(j*j+6*k+7*n+9)+(j+3)**2*(n+1)*(j*j+6*j+(n-1)**2+7)
    b2=k**3+2*k*k*n+8*k*k+k*n*n+11*k*n+21*k+n*n+13*n+18
    b0=(k**5+4*k**4*n+14*k**4+7*k**3*n*n+41*k**3*n+73*k**3+6*k*k*n**3+47*k*k*n*n+143*k*k*n+176*k*k
        +2*k*n**4+22*k*n**3+90*k*n*n+194*k*n+192*k+2*n**4+12*n**3+42*n*n+80*n+72)
    beta=D*b2+b0
    eta=k**4+3*k**3*n+9*k**3+3*k*k*n*n+19*k*k*n+27*k*k+k*n**3+12*k*n*n+34*k*n+29*k+2*n**3+8*n*n+14*n+8
    gamma=-D*D*(n+1)+D*eta+(j+2)**2*(j+4)*(k+4)*(j*j+4*j+n*n+2)
    return D*beta*beta-4*(k+2)*alpha*gamma
# for each (n,k), find smallest integer d>=0 (parity matching N) with Delta>=0, compare to outer threshold
worst=[]
for n in range(1,200,7):
    for k in range(0,200,7):
        dmin=None
        for d in range(0,2*n+k+1):
            if polys(n,k,d*d)>=0: dmin=d; break
        dout=2*math.sqrt(max(0,(n-1)*(n+k+2)))
        if dmin is not None:
            worst.append((dmin/dout if dout>0 else 99, n,k,dmin,round(dout,2)))
worst.sort(); print(worst[:12])
print("--- band R minus outer: min of N/(3t(k+1)^2)")
res=[]
for n in range(1,301,3):
    for k in range(0,301,3):
        dout2=4*(n-1)*(n+k+2)
        for d in range(0,2*n+k+1):
            if (2*n+k-d)%2: continue
            if d*d>=dout2: break
            if polys(n,k,d*d)>=0:
                t=(2*n+k-d)//2
                if t==0: continue
                res.append(((2*n+k)/(3*t*(k+1)**2),n,k,d,t))
res.sort(); print(len(res),res[:10])
from collections import Counter
print("t distribution in band:",sorted(Counter(r[4] for r in res).items()))
# full (not sampled) scan of small n,k for band points with t>=3
big=[]
for n in range(1,160):
    for k in range(0,160):
        dout2=4*(n-1)*(n+k+2)
        for d in range(0,2*n+k+1):
            if (2*n+k-d)%2: continue
            if d*d>=dout2: break
            if polys(n,k,d*d)>=0:
                t=(2*n+k-d)//2
                if t>=3: big.append((n,k,d,t,(2*n+k)>=3*t*(k+1)**2))
print("band points with t>=3 (n,k<160):",len(big),big[:10], "all central:",all(b[4] for b in big))
