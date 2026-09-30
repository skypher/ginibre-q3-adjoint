# Closed form of (HT) in the coordinates m_k = sigma c_k - d B_k/2, n_k = X_k A_k/2 (sigma = N+2, d = a-e, X_k = 2k-N,
# A_k = c_(k-1) - c_(k+1), B_k = c_(k-1) + c_(k+1), V = (a+1)(e+1), sigma^2 - d^2 = 4V).  Claims, checked exactly:
#  (1) 4V D_k = m_k^2 + (4V - X_k^2) A_k^2/4                      (D_k is the diagonal metric diag(1, 4V/X_k^2 - 1) in (m,n))
#  (2) 2V (D_k - D_(k+1)) = m_k n_(k+1) - n_k m_(k+1)             (fan area of the curve p_k = (m_k, n_k))
#  (3) 2V W_ij = m_j n_i - n_j m_i  so (E) <=> 2V(D_j - D_i) >= |m_j n_i - n_j m_i|
#  (4) Q_t = (sigma f - d g)^2/(4V) + (d f - sigma g)^2 (1/t - 1/(4V)), so min over t of the metric bound is the
#      permanent: (HT) <=> 2V(D_j - D_i) >= |m_j n_i| + |n_j m_i|.  Compared with the t-optimization of e_sweep_split.
import sys
from math import comb
from fractions import Fraction
from collections import Counter
def row(a,e):
    N=a+e; c=[0]*(N+1)
    for u in range(e+1):
        s=(-1)**u*comb(e,u)
        for k in range(a+1): c[k+u]+=s*comb(a,k)
    return c
st=Counter()
for a in range(0,31):
    for e in range(0,31):
        c=row(a,e); N=a+e; V=(a+1)*(e+1); d=a-e; sg=N+2; C_=lambda k: c[k] if 0<=k<=N else 0
        assert sg*sg-d*d==4*V
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); B=lambda k: C_(k-1)+C_(k+1); A=lambda k: C_(k-1)-C_(k+1)
        m2=lambda k: 2*sg*C_(k)-d*B(k)          # 2 m_k (integer)
        n2=lambda k: (2*k-N)*A(k)               # 2 n_k (integer)
        for k in range(0,N+2):
            assert 16*V*D(k)==m2(k)**2+(4*V-(2*k-N)**2)*A(k)**2; st['(1)']+=1
            assert 8*V*(D(k)-D(k+1))==m2(k)*n2(k+1)-n2(k)*m2(k+1); st['(2)']+=1
        for j in range((N+1)//2,N+1):
            for i in range(j+1,N+2):
                W=C_(j)*B(i)-B(j)*C_(i)
                assert 8*V*W==(m2(j)*n2(i)-n2(j)*m2(i)); st["(3)"]+=1
                E=D(j)-D(i)
                perm = 8*V*E >= abs(m2(j)*n2(i))+abs(n2(j)*m2(i))
                # t-optimization (exact quadratic, as FM-SEC48)
                x=2*j-N; y=2*i-N
                P=lambda k: 4*D(k)-A(k)**2; Q=lambda k: (2*k-N)**2*A(k)**2
                L=4*E*E+P(j)*P(i); b=16*V*E*E-P(j)*Q(i)-P(i)*Q(j); c0=Q(j)*Q(i)
                opt=(L>0 and 0<b<8*L*V and b*b>=4*L*c0) or (L==0 and b==0 and c0==0)
                if E>=0 and (perm!=opt): st['perm != t-opt']+=1
                st['pairs']+=1
print(dict(st))
