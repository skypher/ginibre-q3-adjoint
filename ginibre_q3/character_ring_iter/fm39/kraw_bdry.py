import sympy as sp
N,v,u=sp.symbols('N v u')
rows={0:(1,[1575,-3150,4095,-1740,345,-30,1]),1:(16,[20475,-44730,58710,-25800,5235,-470,16]),
      2:(96,[96075,-232830,312525,-142620,29752,-2760,96]),3:(256,[193725,-527670,732030,-348432,75056,-7200,256]),
      4:(128,[70875,-220050,318420,-158784,35456,-3520,128])}
b={i:sp.Rational(1,m)*sum(c*u**j for j,c in enumerate(vec)) for i,(m,vec) in rows.items()}
R6=sp.expand(sp.cancel(N**6*sum(b[i].subs(u,v/N)*sp.binomial(4,i)*(sp.Integer(8)/N)**i*(1-sp.Integer(8)/N)**(4-i) for i in range(5))))
def Cpoly(top,m):  # binomial(top, m) as polynomial in N for fixed integer m
    if m<0: return sp.Integer(0)
    return sp.expand(sp.prod([top-i for i in range(m)])/sp.factorial(m))
for n in range(0,8):
    k=N-n; a=N-6
    def c(e):  # c_{k+e} = sum_i (-1)^i C(6,i) C(a, k+e-i) = C(a, a-(k+e-i)) = C(N-6, n-6-e+i)
        return sp.expand(sum((-1)**i*sp.binomial(6,i)*Cpoly(a,n-6-e+i) for i in range(7)))
    D=lambda e: c(e)**2-c(e-1)*c(e+1)
    lhs=sp.expand(D(0)-D(1))
    p=N-2*n
    CNn=Cpoly(N,n)
    rhs=CNn**2*(p+1)*(N-5)*(N-4)*R6.subs(v,p*(p+2))/((n+1)*(k+1)**2*(k+2)*sp.prod([N-i for i in range(6)])**2)
    print('n=',n,' lhs-rhs =',sp.simplify(sp.together(lhs-rhs)))
