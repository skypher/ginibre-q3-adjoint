# Symbolic proof check of the multiplier identity for (E):
#   D_j - D_i - W_ij = sum_(i0 < q/2) delta_(j+q)[(z^i0 - z^(q-i0)) c]
#   D_j - D_i + W_ij = sum_(i0 < q/2) delta_(j+q)[(z^i0 + z^(q-i0)) c] + [q even] 2 delta_(j+q)[z^(q/2) c]
import sympy as sp
L=24; c=sp.symbols('c0:%d'%L)
C_=lambda k: c[k] if 0<=k<L else 0
def D(d,k):
    g=lambda t: d(t)
    return g(k)**2-g(k-1)*g(k+1)
def dl(d,k): return D(d,k)-D(d,k+1)
ok=True
for q in range(1,8):
    j=8; i=j+q+1
    B=lambda k: C_(k-1)+C_(k+1); W=B(i)*C_(j)-C_(i)*B(j)
    base=D(C_,j)-D(C_,i)
    for sg in (1,-1):
        tot=0
        for i0 in range(0,(q+1)//2):
            if 2*i0>=q: break
            m=q-2*i0
            d=lambda t,i0=i0,m=m: C_(t-i0)-sg*C_(t-i0-m)
            tot+=dl(d,j+q)
        if q%2==0 and sg==-1:
            d=lambda t: C_(t-q//2)
            tot+=2*dl(d,j+q)
        ok&= sp.expand(base-sg*W-tot)==0
print('multiplier identity for (E), q = 1..7, both signs, symbolic:',ok)
