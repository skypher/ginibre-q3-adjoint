# FM-SEC71 (luna_max_mars) printed verifier: exact weighted-homogeneous positivity certificate for Q(p,s), p >= 0, s >= 5,
# and the exact small-s forms; hence phi_2(h_(p+2) h_2^2 h_1^(p+2s+2)) > 0 for all p, s >= 0.
import sympy as sp

p, s, S, T, z, Y = sp.symbols('p s S T z Y')
N = p + 2*s + 3

def c_low(j):  # c_(s+j) / binom(N,s-5)
    ratio = sp.prod((p+s+9-m)/(s-5+m) for m in range(1, j+6))
    return sp.factor((p+3-2*j)*ratio/N)

def c_high(j):  # c_(p+s+j) / binom(N,s-5)
    return -c_low(3-j)

def P_of(f):
    D = lambda k: f(k)**2 - f(k-1)*f(k+1)
    return sum(D(k) for k in range(4)) - f(0)*f(3) + f(-1)*f(4)

P_hi = P_of(c_high)
f_low = lambda k: c_low(k-4)
P_lo = P_of(f_low)
U = (c_high(0)+c_low(0), c_high(-1)+c_low(-1))
V = (c_high(4)+c_low(-4), c_high(3)+c_low(-5))
F = sp.factor(P_hi-P_lo-(U[0]*V[1]-U[1]*V[0]))

den = (s**2*(s-4)**2*(s-3)**2*(s-2)**2*(s-1)**2
       *(s+1)**2*(s+2)**2*(s+3)**2*(s+4)*N)
outer = (p+3)*(p+4)*(p+5)*(p+s+8)*(p+2*s+4)
Q = sp.cancel(F*den/outer)
assert sp.denom(Q) == 1
Qpoly = sp.Poly(Q, p, s)
assert (Qpoly.degree(p), Qpoly.degree(s)) == (12, 10)

# Group monomials p^i S^j by weight w=i+2j.
Qs = sp.Poly(sp.expand(Q.subs(s, 5+S)), p, S)
groups = {}
for (i, j), coeff in zip(Qs.monoms(), Qs.coeffs()):
    w = i + 2*j
    eps = w % 2
    assert (i-eps) % 2 == 0
    degree = (i-eps)//2
    groups.setdefault(w, {})[degree] = groups.setdefault(w, {}).get(degree, 0) + coeff
H = {w: sp.Poly.from_dict({(k,):v for k,v in d.items()}, z)
     for w,d in groups.items()}
assert set(H) == set(range(21))

# Verify Q(TY,5+Y^2) = sum_w Y^w T^(w mod 2) H_w(T^2).
lhs = sp.expand(Q.subs({p:T*Y, s:5+Y**2}))
rhs = sum(Y**w*T**(w%2)*H[w].as_expr().subs(z,T**2)
          for w in range(21))
assert sp.expand(lhs-rhs) == 0

positive_weights = (7, 5, 4, 3, 2, 1, 0)
for w in positive_weights:
    assert H[w].as_expr() != 0
    assert all(a >= 0 for a in H[w].all_coeffs())

quadratics = {
    20:(100, 1, -12, 60),
    19:(200, 2, -19, 102),
    17:(40, 467, -4886, 35453),
    15:(8, 47566, -520493, 5475344),
    13:(2, 2226170, -23097548, 394537840),
    11:(2, 16436592, -131531842, 4568259080),
     9:(2, 77562495, -223433420, 35258510164),
}
for w, (m,a,b,c) in quadratics.items():
    rem = sp.Poly(sp.expand(H[w].as_expr()-m*(a*z**2+b*z+c)), z)
    assert all(v >= 0 for v in rem.all_coeffs())
    assert a > 0 and c > 0 and b*b-4*a*c < 0

splits = {
    18:(20, 9), 16:(4, 6), 14:(4,10), 12:(1,12),
    10:(1,20), 8:(1,4), 6:(1,2),
}
negative_powers = {
    18:{1}, 16:{1,2}, 14:{2}, 12:{2}, 10:{2}, 8:{2}, 6:{2}
}
margins = {}
for w, (m,R) in splits.items():
    hp = sp.Poly(H[w].as_expr()/m, z)
    co = {k:int(hp.nth(k)) for k in range(hp.degree()+1)}
    assert {k for k,v in co.items() if v < 0} == negative_powers[w]
    if w == 18:
        low, high = co[0]-2393*R, 36*R*R-2393
    elif w == 16:
        low = co[0]-17982*R*R-87684*R
        high = 8342*R*R-17982*R-87684
    elif w == 14:
        low, high = co[0]-552412*R*R, 162493*R-552412
    elif w == 12:
        low, high = co[0]-26521268*R*R, 6813496*R-26521268
    elif w == 10:
        low, high = co[0]+co[1]*R+co[2]*R*R, 40545008*R-159555812
    elif w == 8:
        low, high = co[0]-447604864*R*R, 129856458*R-447604864
    else:
        low, high = co[0]-254798128*R*R, 174855155*R-254798128
    assert low > 0 and high > 0
    margins[w] = (low, high)

# Requested boundaries.
Q5 = sp.Poly(sp.expand(Q.subs(s,5)), p)
assert all(v > 0 for v in Q5.all_coeffs())
Q6 = sp.Poly(sp.expand(Q.subs(s,6)), p)
assert {k for k in range(Q6.degree()+1) if Q6.nth(k) < 0} == {4}
assert Q6.nth(3) > -Q6.nth(4) and Q6.nth(5) > -Q6.nth(4)
expected_p0 = 48*(s+2)*(s+3)**2*(s+4)*(
    125*s**6+1725*s**5+10415*s**4+34101*s**3
    +64178*s**2+64092*s+26460)
assert sp.expand(Q.subs(p,0)-expected_p0) == 0
assert H[20].as_expr() == 100*(z**2-12*z+60)
assert sp.expand(
    100*S**8*(p**2-6*S)**2+2400*S**10
    -(100*p**4*S**8-1200*p**2*S**9+6000*S**10)
) == 0

# Generate exact mirror forms for s=0,...,4.
def binom_row_low(k, ss):
    if k < 0:
        return sp.Integer(0)
    A = p+2*ss+2
    return sp.binomial(A,k)-sp.binomial(A,k-1)

def phi_fixed_s(ss):
    hi = lambda j: -binom_row_low(ss+3-j, ss)
    lo = lambda k: binom_row_low(k, ss)
    D = lambda f,k: f(k)**2-f(k-1)*f(k+1)
    ph = sum(D(hi,k) for k in range(4))-hi(0)*hi(3)+hi(-1)*hi(4)
    f = lambda k: lo(ss-4+k)
    pl = sum(D(f,k) for k in range(4))-f(0)*f(3)+f(-1)*f(4)
    u = (hi(0)+lo(ss), hi(-1)+lo(ss-1))
    v = (hi(4)+lo(ss-4), hi(3)+lo(ss-5))
    return sp.factor(sp.expand_func(ph-pl-(u[0]*v[1]-u[1]*v[0])))  # main agent: closing parenthesis restored

small = {ss:phi_fixed_s(ss) for ss in range(5)}
assert sp.expand(small[0]-(p**6+3*p**5+7*p**4+9*p**3
                             +280*p**2+996*p+1296)/144) == 0
assert sp.expand(small[1]-(p+3)*(p+4)*(
    p**6+9*p**5+35*p**4-37*p**3+348*p**2+4972*p+13680)/2880) == 0
assert sp.expand(small[2]-(p+3)*(p+4)*(p+5)*(
    p**7+23*p**6+217*p**5+725*p**4-866*p**3+7892*p**2
    +175608*p+580320)/86400) == 0
assert sp.expand(small[3]-(p+3)*(p+4)*(p+5)*(
    p**9+48*p**8+1002*p**7+11100*p**6+61209*p**5+94692*p**4
    +60908*p**3+9434880*p**2+77729760*p+187608960)/3628800) == 0
assert sp.expand(small[4]-(p+3)*(p+4)*(p+5)*(
    p**11+79*p**10+2806*p**9+57594*p**8+725697*p**7
    +5410671*p**6+19830152*p**5+24831896*p**4
    +504685952*p**3+7937334960*p**2+43803934272*p
    +86588974080)/203212800) == 0

assert all(v > 0 for v in sp.Poly(small[0],p).all_coeffs())
assert all(v > 0 for v in sp.Poly(
    small[3]/((p+3)*(p+4)*(p+5)),p).all_coeffs())
assert all(v > 0 for v in sp.Poly(
    small[4]/((p+3)*(p+4)*(p+5)),p).all_coeffs())
I1 = p**6+9*p**5+35*p**4-37*p**3+348*p**2+4972*p+13680
I2 = p**7+23*p**6+217*p**5+725*p**4-866*p**3+7892*p**2+175608*p+580320
assert I1.subs(p,0)>0 and I1.subs(p,1)>0 and 35*2-37>0
assert I2.subs(p,0)>0 and I2.subs(p,1)>0 and 725*2-866>0

print("split margins:", margins)
for ss in range(5):
    print("phi(s=%d) =" % ss, small[ss])
