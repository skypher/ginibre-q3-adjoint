"""FM-SEC123 (luna_max_venus): law-level necessary conditions for positive Walsh realizations; quantile-encoding failures at U_4; B checks for (2^a,4^b), (4^M)."""
import sympy as sp
from fractions import Fraction as F
from math import factorial

th, theta, t = sp.symbols("th theta t", real=True)
u4 = sp.chebyshevu(4, sp.cos(th))  # U_n(2 cos(th)) in this project's convention
density_product = sp.simplify(4*sp.sin(th)**2*u4/sp.pi)
assert sp.trigsimp(sp.expand_trig(
    density_product - 2*(sp.cos(4*th)-sp.cos(6*th))/sp.pi
)) == 0
coef4 = 2*sp.integrate(
    2*(sp.cos(4*th)-sp.cos(6*th))/sp.pi, (th, 0, theta))
assert sp.simplify(
    coef4 - (sp.sin(4*theta)/sp.pi - 2*sp.sin(6*theta)/(3*sp.pi))
) == 0

for k in range(8):
    integrand = 2*(sp.cos((2*k+1)*th)-sp.cos((2*k+3)*th))/sp.pi
    val = sp.simplify(sp.integrate(integrand, (th, 0, sp.pi/2)))
    expected = sp.Rational(8*(-1)**k*(k+1),
                           (2*k+1)*(2*k+3))/sp.pi
    assert sp.simplify(val-expected) == 0

P = -1 + 3*sp.cos(t) + 4*sp.cos(t)**2 - 6*sp.cos(t)**3
assert sp.trigsimp(
    sp.Rational(2,3)*sp.sin(3*t)-sp.Rational(1,2)*sp.sin(4*t)
    - sp.Rational(2,3)*sp.sin(t)*P
) == 0

def atan_partial(x, n):
    return sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(n))

def atan_bounds(x, n):
    p, q = atan_partial(x,n), atan_partial(x,n+1)
    return min(p,q), max(p,q)

a5 = atan_bounds(F(1,5), 12)
a239 = atan_bounds(F(1,239), 12)
pi_lo = 16*a5[0] - 4*a239[1]
pi_hi = 16*a5[1] - 4*a239[0]
assert pi_lo > F(312,100) and pi_hi < 4  # Machin's formula and alternating bounds

x = F(23,10)
sin_poly = sum((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(9))
sin_err = x**18/factorial(18)
cos_poly = sum((-1)**k*x**(2*k)/factorial(2*k) for k in range(9))
cos_err = x**17/factorial(17)
assert sin_poly-sin_err > F(74,100)
assert cos_poly+cos_err < F(-66,100)

p = lambda c: -1+3*c+4*c*c-6*c**3
dp = lambda c: 3+8*c-18*c*c
assert p(F(-33,50)) == F(30461,62500)
assert dp(F(-33,50)) == F(-12651,1250)
print("Odd absolute moments k=0..7: PASS")
print("U4 coefficient identity: PASS")
print("pi, sin, cos rational bounds: PASS")
print("U6 first-sign coefficient factor: positive")


# ---- block ----
from itertools import combinations
from functools import lru_cache
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

@lru_cache(None)
def inv(labels):
    row = {0: 1}
    for n in labels:
        out = {}
        for p,v in row.items():
            for q in range(abs(p-n), p+n+1, 2):
                out[q] = out.get(q, 0) + v
        row = out
    return row.get(0, 0)

def subspaces(d):
    for k in range(d+1):
        for piv in combinations(range(d), k):
            free = [(i,j) for i,p in enumerate(piv)
                    for j in range(p+1,d) if j not in piv]
            for mask in range(1 << len(free)):
                basis = [1 << p for p in piv]
                for h,(i,j) in enumerate(free):
                    if mask >> h & 1:
                        basis[i] |= 1 << j
                H = [0]
                for b in basis:
                    H += [x ^ b for x in H]
                yield H

def test(labels):
    L = len(labels)
    G = 1 << (L-1)
    labs = tuple(sorted(set(labels)))
    total = tuple(labels.count(z) for z in labs)

    def profile(x):
        c = tuple(sum(1 for i,n in enumerate(labels)
                      if n == lab and (x >> i) & 1)
                  for lab in labs)
        return min(c, tuple(total[j]-c[j] for j in range(len(labs))))

    profiles = sorted({profile(x) for x in range(G)})
    ix = {p:i for i,p in enumerate(profiles)}
    qprof = [ix[profile(x)] for x in range(G)]
    sizes = [qprof.count(i) for i in range(len(profiles))]
    target = []
    for p in profiles:
        left = tuple(labs[j] for j,n in enumerate(p) for _ in range(n))
        right = tuple(labs[j] for j,n in enumerate(p)
                      for _ in range(total[j]-n))
        target.append(inv(left)*inv(right)*sizes[len(target)])
    if not any(target):
        return len(profiles), 0, 0, Q(0)

    columns = set()
    for H in subspaces(L-1):
        hist = [0]*len(profiles)
        for x in H:
            hist[qprof[x]] += 1
        columns.add(tuple(hist))
    columns = list(columns)

    lp = linprog(
        [1+(j*65537 % 100003)/100003 for j in range(len(columns))],
        A_eq=np.asarray(columns, dtype=float).T,
        b_eq=np.asarray(target, dtype=float),
        bounds=(0,None), method="highs")
    assert lp.success, (labels, lp.message)
    active = [j for j,x in enumerate(lp.x) if x > 1e-8]
    mat = Matrix([[columns[j][i] for j in active]
                  for i in range(len(target))])
    sol, params = mat.gauss_jordan_solve(Matrix(target))
    assert not params.rows
    coeff = [Q(int(x.p),int(x.q)) for x in sol]
    assert min(coeff) >= 0
    got = [sum(coeff[h]*columns[j][i] for h,j in enumerate(active))
           for i in range(len(target))]
    assert got == [Q(x) for x in target]
    return len(profiles), len(columns), len(active), min(coeff)

rows = []
for L in range(1,7):
    for a in range(L+1):
        labels = (2,)*a + (4,)*(L-a)
        rows.append((labels, test(labels)))
for L in (7,8,9):
    labels = (4,)*L
    rows.append((labels, test(labels)))

print("exact B checks", len(rows), "lists")
for labels, data in rows:
    print(labels, "profiles/columns/terms/mincoef =", data)
print("All rational subgroup re-expansions passed.")
