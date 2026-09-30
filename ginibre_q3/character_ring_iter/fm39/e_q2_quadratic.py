# (E)_(q=2,+) at (j, j+3) as a binary quadratic form in (c_j, c_(j+1)) via the recurrence (k+1)c_(k+1) = d c_k - (N-k+1) c_(k-1).
import sympy as sp
d,N,j=sp.symbols('d N j'); p,q=sp.symbols('p q')   # p = c_j, q = c_(j+1)
c={j:p, j+1:q}
# forward: c_(k+1) = (d c_k - (N-k+1) c_(k-1))/(k+1); backward: c_(k-1) = (d c_k - (k+1) c_(k+1))/(N-k+1)
for k in [j+1,j+2,j+3]:   # c_(k+1)
    c[k+1]=sp.simplify((d*c[k]-(N-k+1)*c[k-1])/(k+1))
for k in [j, j-1]:        # c_(k-1)
    c[k-1]=sp.simplify((d*c[k]-(k+1)*c[k+1])/(N-k+1))
C=lambda k: c[k]
D=lambda k: C(k)**2-C(k-1)*C(k+1)
B=lambda k: C(k-1)+C(k+1)
i=j+3
W=B(i)*C(j)-C(i)*B(j)
for name,expr in [('(E)_(2,+)',D(j)-D(i)+W),('(E)_(2,-)',D(j)-D(i)-W)]:
    e=sp.together(sp.expand(expr)); num,den=sp.fraction(e)
    P=sp.Poly(sp.expand(num),p,q)
    A_=P.coeff_monomial(p**2); B_=P.coeff_monomial(p*q); C2=P.coeff_monomial(q**2)
    disc=sp.factor(B_**2-4*A_*C2)
    print(name,': denominator',sp.factor(den))
    print('   coeff p^2:',sp.factor(A_)); print('   coeff q^2:',sp.factor(C2)); print('   discriminant:',disc)
# --- coverage of the ratio-free (definite) region over the consumer range, and how it meets the residual
from collections import Counter
e=sp.together(sp.expand(D(j)-D(i)+W)); num,den=sp.fraction(e); P=sp.Poly(sp.expand(num),p,q)
A_=sp.lambdify((d,N,j),P.coeff_monomial(p**2)); B_=sp.lambdify((d,N,j),P.coeff_monomial(p*q)); C2=sp.lambdify((d,N,j),P.coeff_monomial(q**2))
st=Counter()
exec(open('split3.py').read())
for a in range(0,41):
    for ee in range(0,41):
        NN=a+ee; dd=a-ee
        for jj in range((NN+1)//2, NN-1):
            A1,B1,C1=A_(dd,NN,jj),B_(dd,NN,jj),C2(dd,NN,jj)
            st['pairs']+=1
            if B1*B1-4*A1*C1<0 and A1>0: st['definite (ratio-free)']+=1
print('(E)_(2,+) consumer pairs, a,e <= 40:',dict(st))
