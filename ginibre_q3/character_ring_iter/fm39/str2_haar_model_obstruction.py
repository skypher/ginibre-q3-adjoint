from math import comb
import sympy as s

r,a,b = s.symbols('r a b')
La = 16*(r+1)-(8+4*r)*a-4*r*b
Lb = -8*r-(6+2*r)*b
Gaa, Gab, Gbb = 16*a-4*a*a+8*b*b, 16*b-2*a*b, 4*a-2*b*b
def L(f):
    return s.expand(s.diff(f,a)*La+s.diff(f,b)*Lb
        +s.diff(f,a,2)*Gaa+2*s.diff(f,a,b)*Gab+s.diff(f,b,2)*Gbb)
M = {(0,0):s.Integer(1)}
for degree in range(1,4):
    inds = [(i,degree-i) for i in range(degree+1)]
    unknown = s.symbols('v0:'+str(len(inds)))
    current = dict(M); current.update(zip(inds,unknown))
    equations = [sum(c*current[ij] for ij,c in
        s.Poly(L(a**i*b**j),a,b).terms()) for i,j in inds]
    solution = s.solve(equations,unknown)
    M.update({ij:s.factor(solution[v]) for ij,v in zip(inds,unknown)})
def E(f):
    return s.factor(sum(c*M[ij] for ij,c in s.Poly(s.expand(f),a,b).terms()))
A, B, S2, H3 = a+2*b, (a-2)**2, a-2, a+b-2
p = r**4+5*r**3+8*r**2+6*r-5
cov = -1152*(2*r+3)*p/((r+2)**2*(r+3)**3*(r+4)**2*(r+5))
cov1 = -96*(r-1)*(2*r+3)/((r+2)**2*(r+3)**2*(r+4))
assert s.factor(E(A*B)-E(A)*E(B)-cov) == 0
assert s.factor(E(A*S2)-E(A)*E(S2)-cov1) == 0
assert s.factor(E(H3)-2*r*(r-1)/((r+2)*(r+3))) == 0
assert s.factor(E(A)-12/((r+2)*(r+3))) == 0
u = s.symbols('u')
assert s.expand(p.subs(r,u+1)) == u**4+9*u**3+29*u**2+41*u+15
assert [E(f).subs(r,1) for f in (A,B,A*B)] == [1,3,2]
assert [E(f).subs(r,2) for f in (A,S2,A*S2)] == [s.Rational(3,5),s.Rational(9,5),s.Rational(4,5)]

# Independent Catalan integration; all arithmetic is exact.
x,y = s.symbols('x y')
def moment(n):
    return 0 if n%2 else comb(n,n//2)//(n//2+1)
def J(f):
    return sum(c*moment(i)*moment(j) for (i,j),c in
               s.Poly(s.expand(f),x,y).terms())
def xy(f):
    return s.expand(f.subs({a:x*x+y*y,b:x*y},simultaneous=True))
for level in range(1,9):
    weight = (x-y)**(2*level); Z = J(weight)
    for f in (A,B,A*B,S2,A*S2,H3):
        assert J(weight*xy(f))/Z == E(f).subs(r,level)
    num,den = 2*level*(level-1),(level+2)*(level+3)
    if level >= 2: assert s.denom(s.Rational(num,den)) > 1
    assert 24 % s.igcd(num,den) == 0

# Character and measure-domination controls.
psi20,psi11 = xy(H3),x*y+1
assert s.expand(xy(S2)-psi20+psi11-1) == 0
assert J((x-y)**2*psi11**2)/2 == 1
assert J((x-y)**2*xy(S2)*psi11)/2 == -1
assert J((x-y)**2*xy(B)*psi11)/2 == -3
density = 1-xy(S2)/6
assert J(density) == 1
assert J(density*xy(S2)) == -s.Rational(1,3)

# Laurent-polynomial checks of the angular identities.
z,w = s.symbols('z w',nonzero=True)
X,Y = z*w+1/(z*w),z/w+w/z
Ux,Uy = [s.Integer(1),X],[s.Integer(1),Y]
for n in range(2,7):
    Ux.append(s.expand(X*Ux[-1]-Ux[-2]))
    Uy.append(s.expand(Y*Uy[-1]-Uy[-2]))
def C(v,n): return (v**n+v**(-n))/2
def T(v,n): return (v**n-v**(-n))/(2*s.I)
for n in range(1,7):
    plus = (2 if n%2==0 else 0)+4*sum(C(z,m)*C(w,m) for m in range(n,0,-2))
    minus = -4*sum(T(z,m)*T(w,m) for m in range(n,0,-2))
    assert s.expand(Ux[n]+Uy[n]-plus) == 0
    assert s.expand(Ux[n]-Uy[n]-minus) == 0
def CT(f): return s.expand(f).coeff(z,0).coeff(w,0)
fz,fw = (2-z*z-z**(-2))/4,(2-w*w-w**(-2))/4
angular_weight = 4*(fz-fw)**2
assert CT(angular_weight) == 1
assert CT(2*angular_weight) == 2
assert CT(4*C(z,2)*C(w,2)*angular_weight) == -2

# Exact inertia of the difference-power kernel.
for level in range(1,7):
    size = 2*level+1
    Q = s.zeros(size)
    for i in range(size): Q[i,2*level-i] = (-1)**i*comb(2*level,i)
    columns = []
    for i in range(level):
        for sign in (1,-1):
            v = s.zeros(size,1); v[i],v[2*level-i] = 1,sign
            columns.append(v)
    v = s.zeros(size,1); v[level] = 1; columns.append(v)
    P = s.Matrix.hstack(*columns); diagonal = P.T*Q*P
    assert diagonal.is_diagonal()
    signs = [s.sign(diagonal[i,i]) for i in range(size)]
    assert signs.count(1) == level+(level%2==0)
    assert signs.count(-1) == level+(level%2==1)
print('PASS: symbolic covariance, independent moments, character controls, angular identities, inertia')
print('Cov at r=1:',cov.subs(r,1))
print('Cov(S1^2,S2) at r=2:',cov1.subs(r,2))
