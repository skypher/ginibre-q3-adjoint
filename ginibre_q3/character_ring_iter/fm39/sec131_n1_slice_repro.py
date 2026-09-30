import sympy as sp

k = sp.Symbol("k", nonnegative=True)
D = sp.Symbol("D", nonnegative=True)
r = sp.Symbol("r")
d = sp.Symbol("d", nonnegative=True)

n = 1
j = k + 1
alpha = D*n*(j**2 + 6*k + 7*n + 9) + (j+3)**2*(n+1)*(
    j**2 + 6*j + n**2 - 2*n + 8
)
beta2 = (k**3 + 2*k**2*n + 8*k**2 + k*n**2 + 11*k*n + 21*k
         + n**2 + 13*n + 18)
beta0 = (k**5 + 4*k**4*n + 14*k**4 + 7*k**3*n**2 + 41*k**3*n
         + 73*k**3 + 6*k**2*n**3 + 47*k**2*n**2 + 143*k**2*n
         + 176*k**2 + 2*k*n**4 + 22*k*n**3 + 90*k*n**2
         + 194*k*n + 192*k + 2*n**4 + 12*n**3 + 42*n**2
         + 80*n + 72)
beta = D*beta2 + beta0
eta = (k**4 + 3*k**3*n + 9*k**3 + 3*k**2*n**2 + 19*k**2*n
       + 27*k**2 + k*n**3 + 12*k*n**2 + 34*k*n + 29*k
       + 2*n**3 + 8*n**2 + 14*n + 8)
gamma = (-D**2*(n+1) + D*eta
         + (j+2)**2*(j+4)*(k+4)*(j**2 + 4*j + n**2 + 2))
delta = sp.factor(D*beta**2 - 4*(k+2)*alpha*gamma)
T = 3*k + 16

assert sp.factor(delta.subs(D, 0)) == (
    -8*(k+2)**2*(k+3)**2*(k+4)**4*(k+5)*(k**2 + 8*k + 14)
)
assert sp.factor(delta.subs(D, T)) == (
    -k*(k+3)**2*(k+4)**2*(k+5)**2
    *(5*k**4 + 86*k**3 + 549*k**2 + 1568*k + 1728)
)
assert sp.factor(sp.Poly(delta, D).coeff_monomial(D**2)) == (
    2*(k+3)**2*(k+4)**3*(k**3 + 8*k**2 + 17*k + 14)
)
assert sp.factor(sp.Poly(delta, D).coeff_monomial(D**3)) == (
    (k+3)**4*(k+4)**2
)

omega0 = (k+3)*D - 2*k**2 - 20*k - 48
omega1 = -d*((k+3)*D - k**2 - 15*k - 48)
omega2 = k*(k+5)*D
assert sp.factor(
    (omega0 + omega1*r + omega2*r**2
     - (1-d*r)*(omega0-omega2*r/d)).subs(D, d**2)
) == 0
assert sp.factor(
    omega0 - omega2/D - (k+3)*(D-3*k-16)
) == 0

print("n=1 exact identities: PASS")

