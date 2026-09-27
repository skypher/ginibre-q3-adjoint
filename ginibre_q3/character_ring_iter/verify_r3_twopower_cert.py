import argparse
from fractions import Fraction
import sympy as sp

argparse.ArgumentParser(
    description="Exact r=3 normalized-cubic certificate."
).parse_args()

d, n, u, v = sp.symbols("d n u v")

Q = (d**8-4*d**7-24*d**6*n-50*d**6+72*d**5*n+164*d**5
     +216*d**4*n**2+1128*d**4*n+1381*d**4
     +144*d**3*n**2+792*d**3*n+992*d**3
     -480*d**2*n**3-2664*d**2*n**2-5088*d**2*n-3348*d**2
     +480*d*n**3+2304*d*n**2+3120*d*n+864*d
     +720*n**4+5760*n**3+16560*n**2+20160*n+8640)

B1 = (d**8+4*d**7-24*d**6*n-82*d**6-120*d**5*n-300*d**5
      +216*d**4*n**2+1448*d**4*n+1861*d**4
      +1584*d**3*n**2+7032*d**3*n+7304*d**3
      -480*d**2*n**3-3624*d**2*n**2-10208*d**2*n-9556*d**2
      -3360*d*n**3-26016*d*n**2-64368*d*n-51072*d
      +720*n**4+1920*n**3-9360*n**2-36480*n-31680)

B2 = (d**6-3*d**5-24*d**4*n-63*d**4+24*d**3*n+3*d**3
      +180*d**2*n**2+716*d**2*n+686*d**2
      -180*d*n**2-716*d*n-624*d
      -480*n**3-3120*n**2-6480*n-4320)

B3 = (d**4-6*d**3-20*d**2*n-49*d**2+20*d*n+54*d
      +60*n**2+300*n+360)

S = d+n+4
D = d+2*n+7
C1 = 2*S*Q-D*B1
C2 = S*Q-4*D**2*B2
C3 = S*Q-12*D**3*B3

# Exact coefficient checks on the two cones.
cone_a = {d: 3*(4+v)+u, n: 4+v}
cone_b = {d: 16+v+u, n: 16+v}
for name, poly in (("Q", Q), ("C1", C1), ("C2", C2), ("C3", C3)):
    for region, subs in (("d>=3n,n>=4", cone_a),
                         ("d>=n,n>=16", cone_b)):
        coeffs = [int(c) for c in
                  sp.Poly(sp.expand(poly.subs(subs)), u, v).coeffs()]
        assert min(coeffs) > 0
        print(name, region, len(coeffs), min(coeffs), max(coeffs),
              flush=True)

polys = [sp.Poly(poly, d, n) for poly in (Q, B1, B2, B3)]
terms = [[(int(c), i, j) for (i, j), c in poly.terms()]
         for poly in polys]

def evaluate(poly_terms, dv, nv):
    return sum(c * dv**i * nv**j for c, i, j in poly_terms)

total = 0
for nlo, nhi, d_factor, K in ((4, 15, 3, 512), (16, 160, 1, 128)):
    for nv in range(nlo, nhi+1):
        count = 0
        for dv in range(d_factor*nv):
            q, b1, b2, b3 = (evaluate(ts, dv, nv) for ts in terms)
            assert q > 0
            den = (dv+nv+4)*q
            big_d = dv+2*nv+7
            base = K*big_d
            for y in range(1, nv+3):
                w = y*(dv+y+2)
                P = Fraction(den+b1*w+8*b2*w*w+16*b3*w*w*w, den)
                assert P <= Fraction(base+2*w, base)**K, (dv, nv, y)
                count += 1
        total += count
        print("screen n=%d d=%d y-cases=%d total=%d" %
              (nv, d_factor*nv, count, total), flush=True)

assert total == 1406802
print("PASS", total, flush=True)
