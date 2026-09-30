"""Exact symbolic and finite verifier for FM-SEC129. Writes no files."""
from argparse import ArgumentParser
from collections import Counter
from math import comb
import sympy as sp


def closed(n, kappa, d):
    j = n+kappa
    H = (n+1)*(j+2)**2*(j+3)**2*(j+4)
    alpha = (d*d*n*(j*j+6*kappa+7*n+9)
             + (j+3)**2*(n+1)*(j*j+6*j+n*n-2*n+8))
    beta2 = (kappa**3+2*kappa**2*n+8*kappa**2
             +kappa*n**2+11*kappa*n+21*kappa+n**2+13*n+18)
    beta0 = (kappa**5+4*kappa**4*n+14*kappa**4
             +7*kappa**3*n**2+41*kappa**3*n+73*kappa**3
             +6*kappa**2*n**3+47*kappa**2*n**2
             +143*kappa**2*n+176*kappa**2
             +2*kappa*n**4+22*kappa*n**3+90*kappa*n**2
             +194*kappa*n+192*kappa+2*n**4+12*n**3
             +42*n**2+80*n+72)
    beta = d*d*beta2+beta0
    eta = (kappa**4+3*kappa**3*n+9*kappa**3
           +3*kappa**2*n**2+19*kappa**2*n+27*kappa**2
           +kappa*n**3+12*kappa*n**2+34*kappa*n+29*kappa
           +2*n**3+8*n**2+14*n+8)
    gamma = (-d**4*(n+1)+d*d*eta
             +(j+2)**2*(j+4)*(kappa+4)*(j*j+4*j+n*n+2))
    discnum = d*d*beta*beta-4*(kappa+2)*alpha*gamma
    return H, alpha, beta, gamma, discnum


def omega_coeffs(n, kappa, d):
    j = n+kappa
    den = (n+1)*(j+2)*(j+3)*(j+4)
    q0 = n*(d*d*(kappa+3)
            -(kappa*kappa*n+kappa*kappa+kappa*n*n
              +10*kappa*n+9*kappa+6*n*n+24*n+18))
    q1 = -d*(d*d*(kappa+3)
             -(kappa*kappa*n+kappa*n*n+10*kappa*n
               +4*kappa+9*n*n+27*n+12))
    q2 = kappa*(j+4)*(d*d-(n*n+n*kappa+n-kappa-2))
    return den, q0, q1, q2


def symbolic_checks():
    N,d,j,x,y = sp.symbols('N d j x y')
    n,kappa = sp.symbols('n kappa')
    c = {j:x, j+1:y}
    for q in (j+1,j+2,j+3):
        c[q+1] = (d*c[q]-(N-q+1)*c[q-1])/(q+1)
    c[j-1] = (d*c[j]-(j+1)*c[j+1])/(N-j+1)
    C = lambda q: c.get(q,sp.Integer(0))
    D = lambda q: C(q)**2-C(q-1)*C(q+1)
    B = lambda q: C(q-1)+C(q+1)
    E = sp.together(D(j)-D(j+3)+B(j+3)*C(j)-C(j+3)*B(j))
    W = sp.together(B(j+3)*C(j)-C(j+3)*B(j))
    sub = {N:2*n+kappa,j:n+kappa}
    H,alpha,beta,gamma,_ = closed(n,kappa,d)

    num,den = sp.fraction(E)
    pe = sp.Poly(num,x,y)
    qA = pe.coeff_monomial(x**2)/den
    qB = pe.coeff_monomial(x*y)/den
    qC = pe.coeff_monomial(y**2)/den
    assert sp.cancel(qA.subs(sub)-(kappa+2)*alpha/H) == 0
    assert sp.cancel(qB.subs(sub)+d*beta/H) == 0
    assert sp.cancel(qC.subs(sub)-gamma/H) == 0
    disc = sp.cancel(qB*qB-4*qA*qC)
    discnum = d*d*beta*beta-4*(kappa+2)*alpha*gamma
    assert sp.cancel(disc.subs(sub)-discnum/H**2) == 0

    nw,dw = sp.fraction(W)
    pw = sp.Poly(nw,x,y)
    wA = pw.coeff_monomial(x**2)/dw
    wB = pw.coeff_monomial(x*y)/dw
    wC = pw.coeff_monomial(y**2)/dw
    denW,q0,q1,q2 = omega_coeffs(n,kappa,d)
    assert sp.cancel(wA.subs(sub)-q0/denW) == 0
    assert sp.cancel(wB.subs(sub)-q1/denW) == 0
    assert sp.cancel(wC.subs(sub)-q2/denW) == 0
    return 'PASS'


def row(a,e):
    N=a+e
    return tuple(sum((-1)**h*comb(e,h)*comb(a,q-h)
                     for h in range(e+1) if 0 <= q-h <= a)
                 for q in range(N+1))


def direct_slack_W(c,N,j):
    v=lambda q: c[q] if 0 <= q <= N else 0
    D=lambda q: v(q)**2-v(q-1)*v(q+1)
    B=lambda q: v(q-1)+v(q+1)
    W=B(j+3)*v(j)-v(j+3)*B(j)
    return D(j)-D(j+3)+W,W


def sign_int(q):
    return (q > 0)-(q < 0)


def sign_a_plus_b_sqrt(a,b,s):
    assert s >= 0
    if b == 0 or s == 0:
        return sign_int(a)
    if a == 0:
        return sign_int(b)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    if a > 0:
        return sign_int(a*a-b*b*s)
    return sign_int(b*b*s-a*a)


def omega_nonnegative_on_outer_interval(n,kappa,d,j,q0,q1,q2):
    if q0 < 0:
        return False
    if n == 1:
        assert d > 0
        if q0*d*d+q1*d+q2 < 0:
            return False
        derivative = q1*d+2*q2
        if q2 > 0 and q1 < 0 and derivative >= 0:
            return 4*q0*q2-q1*q1 >= 0
        return True
    L=(n-1)*(j+2)
    S=d*d-4*L
    assert S >= 0
    h0=L*q0-n*n*q2
    h1=L*q1+d*n*q2
    qR = sign_a_plus_b_sqrt(d*h0+2*n*h1,h0,S)
    if qR < 0:
        return False
    derivativeR = sign_a_plus_b_sqrt(d*q1+4*n*q2,q1,S)
    if q2 > 0 and q1 < 0 and derivativeR >= 0:
        return 4*q0*q2-q1*q1 >= 0
    return True


def finite_screen(amax):
    st=Counter()
    outside=[]; badW=[]; badOmega=[]
    for a in range(amax+1):
        for e in range(amax+1):
            N=a+e
            if N < 2:
                continue
            af,ef=max(a,e),min(a,e)
            d=af-ef
            c=row(af,ef)
            for j in range((N+1)//2,N):
                n=N-j
                kappa=2*j-N
                H,alpha,beta,gamma,discnum=closed(n,kappa,d)
                assert alpha > 0 and H > 0
                val,W=direct_slack_W(c,N,j)
                assert H*val == ((kappa+2)*alpha*c[j]**2
                                 -d*beta*c[j]*c[j+1]
                                 +gamma*c[j+1]**2)
                denW,q0,q1,q2=omega_coeffs(n,kappa,d)
                assert denW*W == (q0*c[j]**2
                                  +q1*c[j]*c[j+1]
                                  +q2*c[j+1]**2)
                st['rows']+=1
                if discnum < 0:
                    continue
                st['R']+=1
                outer=d*d >= 4*(n-1)*(j+2)
                t=ef
                central=N >= 3*t*(kappa+1)**2
                st['R_outer']+=outer
                st['R_central']+=central
                if not outer and central:
                    st['R_central_only']+=1
                if not outer and not central:
                    st['R_outside_both']+=1
                    if len(outside)<5:
                        outside.append((a,e,j,n,kappa,d))
                st['R_W_negative']+=W<0
                if W<0 and len(badW)<5:
                    badW.append((a,e,j,W))
                if outer:
                    passed=omega_nonnegative_on_outer_interval(
                        n,kappa,d,j,q0,q1,q2)
                    st['R_outer_interval_pass']+=passed
                    if not passed and len(badOmega)<5:
                        badOmega.append((a,e,j,n,kappa,d,q0,q1,q2))
    print('symbolic recurrence identities:',symbolic_checks())
    print('finite box 0<=a,e<=',amax,dict(st))
    print('R outside outer and central (first examples):',outside)
    print('R with W<0 (first examples):',badW)
    print('outer interval failures (first examples):',badOmega)


if __name__ == '__main__':
    parser=ArgumentParser(
        description='Exact q=2 plus PD-region verifier; writes no files.')
    parser.add_argument('--amax',type=int,default=70,
                        help='scan 0 <= a,e <= amax (default: 70)')
    finite_screen(parser.parse_args().amax)

