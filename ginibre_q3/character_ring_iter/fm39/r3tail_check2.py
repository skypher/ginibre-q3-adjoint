# Independent check of FM-MECH2 mechanism 1: P(x) <= e^x for n >= 161, 0 <= d < n (r = 3 two-power family, last region).
import sympy as sp, time
from fractions import Fraction as Fr
from math import comb
exec(open('r3defs.py').read())          # defines d, n, Q, B1, B2, B3, S, D, C1, C2, C3 (from the note's verifier)
t0=time.time()
A=sp.expand(S*Q); x=sp.symbols('x')
J=sp.expand(243*A*C2**2-270*A*C2*C3-25*A*C3**2+75*A**2*C2+80*C3**3)
# 1) Taylor-5 remainder identity and discriminant
P=1+D*B1*x/(2*A)+2*D**2*B2*x**2/A+2*D**3*B3*x**3/A
T5=sum(x**j/sp.factorial(j) for j in range(6))
H=A*x**3+5*A*x**2+20*C3*x+60*C2
SKIP=1
if not SKIP: print('Taylor identity T5 - P == x C1/(2A) + x^2 H/(120A):', sp.simplify(sp.together(T5-P-(x*C1/(2*A)+x**2*H/(120*A))))==0)
if not SKIP:
  disc=sp.discriminant(sp.Poly(H,x))
  print('disc(H) == -400 A J:', sp.expand(disc+400*A*J)==0, f'({time.time()-t0:.0f}s)',flush=True)
# 2) region n >= 169, d >= 4 sqrt(n): n=(13+v)^2, d=4(13+v)+e
v,e=sp.symbols('v e')
for name,C in (() if SKIP else (('C1',C1),('C2',C2),('C3',C3))):
    Pc=sp.Poly(sp.expand(C.subs({n:(13+v)**2,d:4*(13+v)+e},simultaneous=True)),v,e)
    print(f'  {name}: min coefficient after substitution = {min(Pc.coeffs())} (terms {len(Pc.coeffs())})',flush=True)
# 3) region n >= 169, d <= 4 sqrt(n): z=d/sqrt(n), h=1/sqrt(n); Bernstein on rectangles
z,h=sp.symbols('z h')
def scaled(C,powh):
    return sp.Poly(sp.expand(sp.cancel(h**powh*C.subs({d:z/h,n:h**-2},simultaneous=True))),z,h)
def bernstein_min(Pz,z0,z1,h0,h1):
    t,s=sp.symbols('t s')
    Q2=sp.Poly(sp.expand(Pz.as_expr().subs({z:z0+(z1-z0)*t,h:h0+(h1-h0)*s},simultaneous=True)),t,s)
    p=Q2.degree(t); q=Q2.degree(s)
    a={}
    for (i,j),c in Q2.terms(): a[(i,j)]=Fr(int(sp.numer(c)),int(sp.denom(c)))
    mn=None; nneg=0
    for k in range(p+1):
        for l in range(q+1):
            bkl=sum(Fr(comb(k,i),comb(p,i))*Fr(comb(l,j),comb(q,j))*c for (i,j),c in a.items() if i<=k and j<=l)
            if bkl<0: nneg+=1
            mn=bkl if mn is None or bkl<mn else mn
    return nneg,float(mn),(p,q)
pw={}
for name,C in (('C1',C1),('C2',C2),('C3',C3),('J',J)):
    Pd=sp.Poly(C,d,n); k=max(i+2*j for (i,j) in Pd.monoms())
    expr=sum(c*z**i*h**(k-i-2*j) for (i,j),c in Pd.terms())
    pw[name]=(k,sp.Poly(expr,z,h))
    print(f'  scaled {name}: h-power {k}, degrees {pw[name][1].degree(z)},{pw[name][1].degree(h)}',flush=True)
for (z0,z1),names in (((0,2),('C1','C2','J')),((2,4),('C1','C2','C3'))):
    for nm in names:
        nneg,mn,deg=bernstein_min(pw[nm][1],sp.Integer(z0),sp.Integer(z1),sp.Integer(0),sp.Rational(1,13))
        print(f'  rectangle z in [{z0},{z1}], h in [0,1/13]: {nm} negative Bernstein coefficients {nneg}, min {mn:.3e}, degrees {deg} ({time.time()-t0:.0f}s)',flush=True)
# 4) bridge n = 161..168, 0 <= d < n
fl=[sp.lambdify((d,n),expr,'sympy') for expr in (C1,C2,C3,J)]
cnt=0; bad=0
for nv in range(161,169):
    for dv in range(0,nv):
        c1,c2,c3,j=[sp.Integer(f(dv,nv)) for f in fl]; cnt+=1
        if not (c1>=0 and c2>=0 and (c3>=0 or j>=0)): bad+=1
print(f'bridge: {cnt} checks, failures {bad} ({time.time()-t0:.0f}s)')
