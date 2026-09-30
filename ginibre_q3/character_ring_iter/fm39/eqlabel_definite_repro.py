
from fractions import Fraction as F
from math import comb, prod
import sympy as sp

def qadd(p, q):
    return tuple(p[i] + q[i] for i in range(3))

def qsub(p, q):
    return tuple(p[i] - q[i] for i in range(3))

def qmul(p, q):
    return (p[0]*q[0], p[0]*q[1] + p[1]*q[0], p[1]*q[1])

def form_shift(m, q):
    # Coefficient pair for c_(alpha+h), relative to (c_alpha,c_(alpha+1)).
    M, d = m + q, 2*m - 2*q + 1
    v = {0: (F(1), F(0)), 1: (F(0), F(1))}
    v[-1] = tuple(
        (F(d)*v[0][t] - F(M-1)*v[1][t]) / F(M+2)
        for t in (0, 1)
    )
    for h in range(1, 8):
        v[h+1] = tuple(
            (F(d)*v[h][t] - F(M+2-h)*v[h-1][t]) / F(M+h-1)
            for t in (0, 1)
        )

    def D(h):
        return qsub(qmul(v[h], v[h]), qmul(v[h-1], v[h+1]))

    def P(x):
        z = (F(0), F(0), F(0))
        for h in range(x, x+4):
            z = qadd(z, D(h))
        return qadd(qsub(z, qmul(v[x], v[x+3])), qmul(v[x-1], v[x+4]))

    def AB(x):
        return (
            tuple(v[x][t] - v[x+3][t] for t in (0, 1)),
            tuple(v[x-1][t] - v[x+4][t] for t in (0, 1)),
        )

    A0, B0 = AB(0)
    A4, B4 = AB(4)
    wedge = qsub(qmul(A0, B4), qmul(B0, A4))
    return qsub(qsub(P(0), P(4)), wedge)

def form_clipped(a, e, u, v, w):
    # Source-style evaluator with zero extension outside 0 <= k <= N.
    A, B, C = u+1, v+1, w+1
    N, d = a+e, a-e
    alpha = (N+A-B-C)//2
    gamma = (N+A+B-C)//2
    tau = gamma + 1
    x0 = max(0, min(N-1, N//2))
    lo = min(alpha-A-2, -2)
    hi = max(tau+C+3, N+3)

    def value(p0, p1):
        c = {x0: F(p0), x0+1: F(p1)}
        for k in range(x0+1, hi):
            c[k+1] = (d*c[k] - (N-k+1)*c[k-1]) / (k+1)
        for k in range(x0, lo, -1):
            c[k-1] = (d*c[k] - (k+1)*c[k+1]) / (N-k+1)

        cc = lambda k: c.get(k, F(0)) if 0 <= k <= N else F(0)
        D = lambda k: cc(k)**2 - cc(k-1)*cc(k+1)

        def P(x):
            return (sum((D(k) for k in range(x, x+C+1)), F(0))
                    - cc(x)*cc(x+C) + cc(x-1)*cc(x+C+1))

        U = lambda x: cc(x)-cc(x+C)
        V = lambda x: cc(x-1)-cc(x+C+1)
        return P(alpha)-P(tau)-(U(alpha)*V(tau)-V(alpha)*U(tau))

    xx, yy = value(1, 0), value(0, 1)
    xy = value(1, 1) - xx - yy
    Q = (xx, xy, yy)
    det = xx*yy - xy*xy/F(4)
    return alpha, gamma, Q, det

def detQ(Q):
    return Q[0]*Q[2] - Q[1]*Q[1]/4

def delta_line(q, sign):
    M = 2*q if sign == 1 else 2*q-1
    return prod(M+h-1 for h in range(1, 8))

def derive_line(sign):
    t = sp.symbols('t')
    R = []
    for i, label in enumerate(("Qxx", "Qxy", "Qyy")):
        points = []
        for q in range(4, 19):
            m = q if sign == 1 else q-1
            z = form_shift(m, q)[i] * F(delta_line(q, sign)**2)
            points.append((sp.Integer(q),
                           sp.Rational(z.numerator, z.denominator)))
        p = sp.factor(sp.interpolate(points, t))
        assert sp.degree(p, t) <= 14
        R.append(p)
        print(("a=e+1" if sign == 1 else "a=e-1"),
              label, "over Delta^2 =", p)
    det_num = sp.factor(R[0]*R[2] - R[1]**2/4)
    print(("a=e+1" if sign == 1 else "a=e-1"),
          "det numerator over Delta^4 =", det_num)

derive_line(1)
derive_line(-1)

print("Support-clipped small q cases:")
for sign in (1, -1):
    for q in (1, 2, 3):
        a = 2*q if sign == 1 else 2*q-2
        e = 2*q-1
        alpha, gamma, Q, det = form_clipped(a, e, 2, 2, 2)
        print(sign, q, a, e, alpha, gamma, Q, det)

print("Boundary tests:")
cases = [
    ("r2 gamma=N+1", 0, 1, 2, 2, 2),
    ("r2 gamma=N",   2, 1, 2, 2, 2),
    ("r3 gamma=N",   0, 3, 2, 2, 2),
    ("r2 a=e+2",     3, 1, 3, 2, 2),
    ("r3 a=e-2",     1, 3, 3, 2, 2),
    ("r3 a=e+2",     5, 3, 3, 2, 2),
]
for label, a, e, u, v, w in cases:
    alpha, gamma, Q, det = form_clipped(a, e, u, v, w)
    assert (a+u+v+w) % 2 == 0
    print(label, (a, e, u, v, w), "N=", a+e,
          "alpha=", alpha, "gamma=", gamma, "Q=", Q, "det=", det)

print("Fixed r=2, equal labels, low m:")
for m in range(6):
    alpha, gamma, Q, det = form_clipped(2*m, 1, 2, 2, 2)
    assert Q[0] > 0 and Q[2] > 0 and det > 0
    print(m, "a=", 2*m, "gamma=", gamma, "Q=", Q, "det=", det)

for m in range(6, 136):
    Q = form_shift(m, 1)
    assert Q[0] > 0 and Q[2] > 0 and detQ(Q) > 0
print("Exact scan m=6..135: all forms are positive definite.")

# Exact interpolation for the determinant on r=2, labels (2,2,2).
t = sp.symbols('t')
Delta_m = lambda m: prod(m+j for j in range(1, 8))
R = []
for i in range(3):
    points = []
    for m in range(7, 22):
        z = form_shift(m, 1)[i] * F(Delta_m(m)**2)
        points.append((sp.Integer(m),
                       sp.Rational(z.numerator, z.denominator)))
    p = sp.factor(sp.interpolate(points, t))
    assert sp.degree(p, t) <= 14
    R.append(p)
print("r2 equal-label Q numerators over Delta_m^2:", *R, sep="\n")
print("r2 equal-label determinant numerator over Delta_m^4 =",
      sp.factor(R[0]*R[2] - R[1]**2/4))

for m in (135, 136):
    Q = form_shift(m, 1)
    det = detQ(Q)
    c = lambda k: ((comb(2*m, k) if 0 <= k <= 2*m else 0)
                   - (comb(2*m, k-1) if 0 <= k-1 <= 2*m else 0))
    X, Y = F(c(m-1)), F(c(m))
    actual_phi = Q[0]*X*X + Q[1]*X*Y + Q[2]*Y*Y
    print("m=", m, "a=", 2*m, "N=", 2*m+1, "gamma=", m+2,
          "Q=", Q, "det=", det, "actual_phi_positive=", actual_phi > 0)
