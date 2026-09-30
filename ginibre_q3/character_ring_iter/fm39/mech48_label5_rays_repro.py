import argparse
from fractions import Fraction as Q
from math import comb
import sympy as sp

argparse.ArgumentParser(
    description="FM-MECH48: exact label-5 ray certificates"
).parse_args()

q, z = sp.symbols("q z")
A = 1+q+q*q+q**3+q**4
B = 1+q+q*q
F = 16*A*z*z-16*B*z+3
D = sp.expand(F.subs(z, 1))
rho = Q(6, 7)

assert sp.expand(
    D-sp.Rational(21, 16)
    -(4*q+3)**2*(16*q*q-8*q+3)/16
) == 0

def transforms(deg, lo, hi):
    return [
        [
            sum(
                Q(comb(i, t))*lo**(i-t)*(hi-lo)**t
                *Q(comb(k, t), comb(deg, t))
                for t in range(min(i, k)+1)
            )
            for i in range(deg+1)
        ]
        for k in range(deg+1)
    ]

def certify(expr, nq=16, nz=16):
    poly = sp.Poly(sp.expand(expr), q, z)
    dq, dz = poly.degree(q), poly.degree(z)
    coeff = [
        [Q(poly.coeff_monomial(q**i*z**j))
         for j in range(dz+1)]
        for i in range(dq+1)
    ]
    tz = [
        transforms(dz, rho*j/nz, rho*(j+1)/nz)
        for j in range(nz)
    ]
    minimum = None
    bad = 0
    for iq in range(nq):
        tq = transforms(
            dq, Q(-1)+Q(2*iq, nq),
            Q(-1)+Q(2*(iq+1), nq))
        part = [
            [sum(tq[k][i]*coeff[i][j]
                 for i in range(dq+1))
             for j in range(dz+1)]
            for k in range(dq+1)
        ]
        for T in tz:
            vals = [
                sum(part[k][j]*T[l][j]
                    for j in range(dz+1))
                for k in range(dq+1)
                for l in range(dz+1)
            ]
            mm = min(vals)
            minimum = mm if minimum is None else min(minimum, mm)
            bad += mm < 0
    return minimum, bad

certs = {
    "root cutoff": F.subs(z, sp.Rational(6, 7)),
    "vertex cutoff": 2*sp.Rational(6, 7)*16*A-16*B,
    "weighted suppression 3/4":
        sp.Rational(9, 16)*D**2-z*F**2,
    "raw lower -3": 3*D+F,
    "raw upper 3": 3*D-F,
}
for name, expr in certs.items():
    result = certify(expr)
    print(name, result)
    assert result[1] == 0

constant = Q(36369, 200)
def cutoff_bound(H):
    return (constant*Q(4, 5)**(H-2)
            *(Q(15, 7)*H+1)*(Q(15, 7)*H+2))

assert cutoff_bound(70) > 1
assert cutoff_bound(71) < 1
assert Q(4, 5)*Q(72, 71)**2 < 1
print("High-count cutoff:", 71)
print("New core count:", comb(74, 4)-comb(73, 3))
print("PASS")
