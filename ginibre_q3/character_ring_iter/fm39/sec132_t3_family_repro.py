import sympy as sp

n, k, t = sp.symbols("n k t", integer=True, nonnegative=True)
d = 2*n + k - 2*t
D = d**2
j = n + k

alpha = (
    D*n*(j**2 + 6*k + 7*n + 9)
    + (j+3)**2*(n+1)*(j**2 + 6*j + (n-1)**2 + 7)
)
beta2 = (
    k**3 + 2*k**2*n + 8*k**2 + k*n**2 + 11*k*n + 21*k
    + n**2 + 13*n + 18
)
beta0 = (
    k**5 + 4*k**4*n + 14*k**4 + 7*k**3*n**2 + 41*k**3*n
    + 73*k**3 + 6*k**2*n**3 + 47*k**2*n**2 + 143*k**2*n
    + 176*k**2 + 2*k*n**4 + 22*k*n**3 + 90*k*n**2
    + 194*k*n + 192*k + 2*n**4 + 12*n**3 + 42*n**2
    + 80*n + 72
)
beta = D*beta2 + beta0
eta = (
    k**4 + 3*k**3*n + 9*k**3 + 3*k**2*n**2 + 19*k**2*n
    + 27*k**2 + k*n**3 + 12*k*n**2 + 34*k*n + 29*k
    + 2*n**3 + 8*n**2 + 14*n + 8
)
gamma = (
    -D**2*(n+1) + D*eta
    + (j+2)**2*(j+4)*(k+4)*(j**2 + 4*j + n**2 + 2)
)
delta = sp.expand(D*beta**2 - 4*(k+2)*alpha*gamma)

delta_t3 = sp.Poly(sp.expand(delta.subs(t, 3)), n)
assert delta_t3.degree() == 10
assert sp.factor(delta_t3.coeff_monomial(n**10)) == 144

dout = 4*(n-1)*(n+k+2)
band_gap_t3 = sp.expand((dout - D).subs(t, 3))
central_margin_t3 = sp.expand(
    (2*n + k - 3*t*(k+1)**2).subs(t, 3)
)
assert sp.expand(band_gap_t3 - (28*n + 8*k - k**2 - 44)) == 0
assert sp.expand(
    central_margin_t3 - (2*n + k - 9*(k+1)**2)
) == 0

witness = {n: 1283, k: 9, t: 3}
values = (
    d.subs(witness),
    (2*n+k).subs(witness),
    delta.subs(witness),
    band_gap_t3.subs(witness),
    central_margin_t3.subs(witness),
)
assert tuple(map(int, values)) == (
    2569,
    2575,
    830827222272722822584897152000,
    35871,
    1675,
)

T2 = n**2 + n*k + n - k - 2
assert sp.factor(dout - T2 - 3*(n-1)*(n+k+2)) == 0

print("t=3 leading coefficient:", delta_t3.coeff_monomial(n**10))
print("witness (d, N, Delta, Dout-D, central margin):", values)
print("outer omega2 identity: PASS")

