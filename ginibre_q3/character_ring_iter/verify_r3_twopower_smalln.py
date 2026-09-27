import sympy as sp, sys
d=sp.symbols('d', nonnegative=True)
def Fsym(n):
    a=d+2*n
    def rho(k):
        if k>=0:
            r=sp.Integer(1)
            for i in range(k): r*= (n+d-i)/sp.Integer(1)/(n+i+1)
            return r
        h=-k; r=sp.Integer(1)
        for i in range(h): r*= sp.Integer(n-i)/(n+d+1+i)
        return r
    u={j: rho((j-1)//2)-rho(-(j+1)//2) for j in (1,3,5,7)}
    k7=u[5]-5*u[3]+10*u[1]; k5=-u[7]+14*u[3]-35*u[1]; k3=5*u[7]-14*u[5]+35*u[1]; k1=-10*u[7]+35*u[5]-35*u[3]
    al,be,ga,de=(k7,-5*k7+k5,6*k7-3*k5+k3,-k7+k5-k3+k1)
    def B(top,x):
        if x<0: return sp.Integer(0)
        return sp.prod([top-i for i in range(x)])/sp.factorial(x)
    return [sp.together(al*B(a+6,x)+be*B(a+4,x-1)+ga*B(a+2,x-2)+de*B(a,x-3)) for x in range(n+3)]
for n in range(int(sys.argv[1]),int(sys.argv[2])+1):
    F=Fsym(n); end=F[-1]; ok=True; notes=[]
    for x in range(0,n+2):
        num,den=sp.fraction(sp.factor(end-F[x]))
        Pn=sp.Poly(sp.expand(num),d); Pd=sp.Poly(sp.expand(den),d)
        rn=Pn.count_roots(0,None) if Pn.degree()>0 else 0
        rd=Pd.count_roots(0,None) if Pd.degree()>0 else 0
        v0=(num/den).subs(d,0)
        if not (rn==0 and rd==0 and v0>=0):
            # allow roots only if nonnegativity holds: report
            ok=False; notes.append((x,rn,rd,v0))
    print(f"n={n}: F(n+2)-F(x) >= 0 for all x<=n+1 and all d>=0 via Sturm: {ok}", notes[:4], flush=True)
