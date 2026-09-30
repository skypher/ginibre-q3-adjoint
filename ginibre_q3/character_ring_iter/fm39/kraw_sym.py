# Symbolic check (generic interior case) of phi_3(hatS_p h_1^a) = D_k - D_{k+1} = claimed factorized form.
import sympy as sp
a,k,v,Ns=sp.symbols('a k v N')
exec(open('kraw_check.py').read().split("# 1) Sturm")[0].split("exec(open('n0check.py')")[0]) if False else None
rows={0:(1,[1575,-3150,4095,-1740,345,-30,1]),1:(16,[20475,-44730,58710,-25800,5235,-470,16]),
      2:(96,[96075,-232830,312525,-142620,29752,-2760,96]),3:(256,[193725,-527670,732030,-348432,75056,-7200,256]),
      4:(128,[70875,-220050,318420,-158784,35456,-3520,128])}
u=sp.symbols('u')
b={i:sp.Rational(1,m)*sum(c*u**j for j,c in enumerate(vec)) for i,(m,vec) in rows.items()}
R6=sp.expand(sp.cancel(Ns**6*sum(b[i].subs(u,v/Ns)*sp.binomial(4,i)*(8/Ns)**i*(1-8/Ns)**(4-i) for i in range(5))))
def ratio(e):  # C(a,k+e)/C(a,k) as a rational function
    if e==0: return sp.Integer(1)
    if e>0: return sp.prod([(a-k-i) for i in range(e)])/sp.prod([(k+1+i) for i in range(e)])
    e=-e;   return sp.prod([(k-i) for i in range(e)])/sp.prod([(a-k+1+i) for i in range(e)])
def c(j):  # c_j / C(a,k), j = k + offset
    return sum((-1)**i*sp.binomial(6,i)*ratio(j-i) for i in range(7))
def D(off): return c(off)**2-c(off-1)*c(off+1)
lhs=sp.together(D(0)-D(1))
N=a+6; p=2*k-N; n=N-k
CNk_over_Cak=sp.prod([(a+i) for i in range(1,7)])/sp.prod([(a-k+i) for i in range(1,7)])
rhs=CNk_over_Cak**2*(p+1)*(N-5)*(N-4)*R6.subs({Ns:N,v:p*(p+2)})/((n+1)*(k+1)**2*(k+2)*sp.prod([N-i for i in range(6)])**2)
print("floats in R6,lhs,rhs:",bool(R6.atoms(sp.Float)),bool(lhs.atoms(sp.Float)),bool(rhs.atoms(sp.Float)))
R6=sp.nsimplify(R6,rational=True) if R6.atoms(sp.Float) else R6
rhs=CNk_over_Cak**2*(p+1)*(N-5)*(N-4)*R6.subs({Ns:N,v:p*(p+2)})/((n+1)*(k+1)**2*(k+2)*sp.prod([N-i for i in range(6)])**2)
diff=sp.cancel(sp.together(lhs-rhs))
print('generic rational identity lhs - rhs =',diff)
