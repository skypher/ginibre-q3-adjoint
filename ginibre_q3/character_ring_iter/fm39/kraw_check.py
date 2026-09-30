# Independent check of FM-MECH2 mechanism 2: phi_3(hatS_p h_1^a) >= 0 for all p>=2, a>=0.
import sympy as sp
from fractions import Fraction as Fr
from math import comb
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
N_,v,u,z=sp.symbols('N v u z')
rows={0:(1,[1575,-3150,4095,-1740,345,-30,1]),1:(16,[20475,-44730,58710,-25800,5235,-470,16]),
      2:(96,[96075,-232830,312525,-142620,29752,-2760,96]),3:(256,[193725,-527670,732030,-348432,75056,-7200,256]),
      4:(128,[70875,-220050,318420,-158784,35456,-3520,128])}
b={i:sp.Rational(1,mlt)*sum(c*u**j for j,c in enumerate(vec)) for i,(mlt,vec) in rows.items()}
# 1) Sturm / real-root positivity of each b_i on u>=0
for i in range(5):
    P=sp.Poly(b[i]*rows[i][0],u)
    nroots=len([r for r in sp.real_roots(P) if r>=0])
    print(f'b_{i}: roots in [0,oo): {nroots}; b_{i}(0)={b[i].subs(u,0)}; leading coeff sign {sp.sign(P.LC())}')
# 2) reconstruct R6(N,v) = N^6 * sum_i b_i(v/N) C(4,i) z^i (1-z)^(4-i), z = 8/N
R6=sp.expand(sp.together(N_**6*sum(b[i].subs(u,v/N_)*sp.binomial(4,i)*(8/N_)**i*(1-8/N_)**(4-i) for i in range(5))))
R6=sp.simplify(R6); print('R6 is polynomial:',sp.Poly(sp.expand(R6),N_,v).is_polynomial if hasattr(sp.Poly(sp.expand(R6),N_,v),'is_polynomial') else True, ' total degree', sp.Poly(sp.expand(R6),N_,v).total_degree())
R6f=sp.lambdify((N_,v),sp.expand(R6),'sympy')
def claimed(N,p):
    if (N-p)%2 or p>N: return Fr(0)
    k=(N+p)//2; n=(N-p)//2
    num=Fr(comb(N,k))**2*(p+1)*(N-5)*(N-4)*Fr(int(sp.Integer(R6f(N,p*(p+2)))) if True else 0)
    den=Fr((n+1)*(k+1)**2*(k+2))*Fr(sp.prod([N-i for i in range(6)]))**2
    return num/den
bad=0; checked=0; neg=0
for N in range(6,61):
    a=N-6
    for p in range(2,N+1):
        if (N-p)%2: continue
        direct=phi(3,shat(p),a)
        if direct<0: neg+=1
        if N>=8:
            c=claimed(N,p); checked+=1
            if c!=direct: bad+=1; 
            if bad and bad<=3 and c!=direct: print('MISMATCH',N,p,direct,c)
print(f'direct phi_3(hatS_p h_1^a): negatives {neg}; identity checked on {checked} (N,p) with N in 8..60: mismatches {bad}')
print('N=6 values p=2,4,6:',[phi(3,shat(p),0) for p in (2,4,6)],' N=7 p=3,5,7:',[phi(3,shat(p),1) for p in (3,5,7)])
