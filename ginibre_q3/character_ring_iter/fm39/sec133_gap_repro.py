import sympy as sp

n, k, t = sp.symbols("n k t", integer=True, nonnegative=True)
d = 2*n + k - 2*t
D = d**2
j = n + k

alpha = D*n*(j**2 + 6*k + 7*n + 9) + (j+3)**2*(n+1)*(
    j**2 + 6*j + (n-1)**2 + 7
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
dout = 4*(n-1)*(n+k+2)

assert sp.expand(
    dout-D-(8*n*t+4*n+4*k*t-4*k-k**2-4*t**2-8)
) == 0

delta3 = sp.Poly(sp.expand(delta.subs(t, 3)), n)
assert delta3.degree() == 10
assert delta3.coeff_monomial(n**10) == 144

q = sp.Symbol("q", integer=True, nonnegative=True)
delta_q = sp.Poly(
    sp.expand(delta.subs({n: 16*q**2, k: q, t: 3})), q
)
assert delta_q.degree() == 20
assert delta_q.coeff_monomial(q**20) == 69423851372544

band_q = sp.expand((dout-D).subs({n: 16*q**2, k: q, t: 3}))
central_q = sp.expand(
    (2*n+k-3*t*(k+1)**2).subs({n: 16*q**2, k: q, t: 3})
)
assert band_q == 447*q**2 + 8*q - 44
assert central_q == 23*q**2 - 17*q - 9

q9 = {n: 1296, k: 9, t: 3}
assert int(delta.subs(q9)) == 21094807986113839754859332527104
assert int((dout-D).subs(q9)) == 36235
assert int((2*n+k-3*t*(k+1)**2).subs(q9)) == 1701

T2 = n**2 + n*k + n - k - 2
assert sp.factor(dout-T2-3*(n-1)*(n+k+2)) == 0

def delta_int(nv, kv, tv):
    dv = 2*nv + kv - 2*tv
    Dv = dv*dv
    jv = nv + kv
    av = Dv*nv*(jv*jv+6*kv+7*nv+9) + (jv+3)**2*(nv+1)*(
        jv*jv+6*jv+(nv-1)**2+7
    )
    b2v = (
        kv**3 + 2*kv**2*nv + 8*kv**2 + kv*nv**2
        + 11*kv*nv + 21*kv + nv**2 + 13*nv + 18
    )
    b0v = (
        kv**5 + 4*kv**4*nv + 14*kv**4 + 7*kv**3*nv**2
        + 41*kv**3*nv + 73*kv**3 + 6*kv**2*nv**3
        + 47*kv**2*nv**2 + 143*kv**2*nv + 176*kv**2
        + 2*kv*nv**4 + 22*kv*nv**3 + 90*kv*nv**2
        + 194*kv*nv + 192*kv + 2*nv**4 + 12*nv**3
        + 42*nv**2 + 80*nv + 72
    )
    betav = Dv*b2v + b0v
    etav = (
        kv**4 + 3*kv**3*nv + 9*kv**3 + 3*kv**2*nv**2
        + 19*kv**2*nv + 27*kv**2 + kv*nv**3
        + 12*kv*nv**2 + 34*kv*nv + 29*kv
        + 2*nv**3 + 8*nv**2 + 14*nv + 8
    )
    gv = (
        -Dv*Dv*(nv+1) + Dv*etav
        + (jv+2)**2*(jv+4)*(kv+4)*(jv*jv+4*jv+nv*nv+2)
    )
    return Dv*betav*betav - 4*(kv+2)*av*gv

count = 0
for tv in range(3, 8):
    for kv in range(16):
        low = kv*kv + 4*kv + 4*tv*tv + 8 - 4*kv*tv
        nmin = max(1, low//(8*tv+4) + 1)
        high = 3*tv*(kv+1)**2 - kv
        nmax = (high-1)//2
        for nv in range(nmin, nmax+1):
            dv = 2*nv + kv - 2*tv
            if dv < 0:
                continue
            count += 1
            assert delta_int(nv, kv, tv) < 0
assert count == 55583

print("D_out-D identity: PASS")
print("t=3 leading coefficient in n:", delta3.coeff_monomial(n**10))
print("scaled-family leading coefficient in q:", delta_q.coeff_monomial(q**20))
print("q=9 witness exact checks: PASS")
print("outer omega2 identity: PASS")
print("finite gap screen t=3..7, k=0..15:", count, "PASS")

