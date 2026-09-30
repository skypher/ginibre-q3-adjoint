
from math import comb
from fractions import Fraction
from collections import Counter
import sympy as sp

a, k = sp.symbols('a k', integer=True, nonnegative=True)
Q = lambda j: (a - 2*j + 1)**2 + a + 1

c0 = (a + 1 - 2*k) / (a - k + 1)
cminus = k*(a - 2*k + 3) / ((a - k + 1)*(a - k + 2))
cplus = (a - 2*k - 1) / (k + 1)
Dnorm = sp.factor(c0**2 - cminus*cplus)
Dformula = (a + 1)*Q(k) / ((k + 1)*(a - k + 1)**2*(a - k + 2))
assert sp.factor(Dnorm - Dformula) == 0

Dn = lambda j: Q(j) / ((j + 1)*(a - j + 1)**2*(a - j + 2))
ratio = sp.factor(
    (k/(a-k+1)*(a-k)/(k+1))**2 * Dn(k-1)*Dn(k+1) / Dn(k)**2
)
logconc_factor = (
    2*(a+2)*(a-2*k+1)**2*((a-2*k+1)**2+3*a+5)
    / ((k+2)*(a-k+3)*Q(k)**2)
)
assert sp.factor((1-ratio) - logconc_factor) == 0
print('symbolic e=1 Turan formula and log-concavity factorization: PASS')

def row(a, e):
    N = a + e
    c = [
        sum((-1)**i * comb(e, i) * comb(a, m-i)
            for i in range(e+1) if 0 <= m-i <= a)
        for m in range(N+1)
    ]
    cp = [0] + c + [0]
    D = [cp[m+1]**2 - cp[m]*cp[m+2] for m in range(N+1)]
    return c, D, N

windows = 0
for aa in range(81):
    c, D, N = row(aa, 1)
    assert D[N] == 1
    for j in range(aa+1):
        b = comb(aa, j)
        q = (aa-2*j+1)**2 + aa + 1
        assert Fraction(
            b*b*(aa+1)*q,
            (j+1)*(aa-j+1)**2*(aa-j+2)
        ) == D[j]
    assert all(D[j]**2 >= D[j-1]*D[j+1] for j in range(1, N))
    for x in range(N+1):
        for y in range(x, N+1):
            S = sum(D[x:y+1])
            assert S*S >= (y-x+1)**2 * D[x] * D[y]
            windows += 1
print('exact e=1 grid: a=0..80; all checks pass; windows=', windows)

stats = Counter()
misses = []
for r in range(2, 11):
    e = 2*r - 3
    for aa in range(60):
        c, D, N = row(aa, e)
        for C in range(3, N+1):
            for x in range(N-C+1):
                y = x + C
                if 2*x < N-C:
                    continue
                A, B = D[x], D[y]
                if A <= 0 or B <= 0 or not (A > B and A < C*C*B):
                    continue
                stats['medium'] += 1
                if x <= N//2:
                    stats['crossing'] += 1
                    h, delta = N-2*x, 2*x+C-N
                    if (h+1)**2*A >= delta**2*B:
                        stats['plateau_covers'] += 1
                    else:
                        stats['plateau_misses'] += 1
                        if len(misses) < 5:
                            misses.append((r, aa, x, C, N, A, B,
                                           h, delta, Fraction(A, B)))
                else:
                    stats['right_only'] += 1
print('exact medium-window screen r=2..10,a=0..59:', dict(stats))
print('first plateau-condition misses:', misses)

c, D, N = row(7, 1)
x, C = 4, 3
A, B = D[x], D[x+C]
S = sum(D[x:x+C+1])
h, delta = N-2*x, 2*x+C-N
print('e=1 plateau-miss example: c=', c, 'Dwindow=', D[x:x+C+1],
      'A,B,S=', A, B, S)
print('h,delta,lower bound=', h, delta, (h+1)*A+delta*B,
      'target-square margin=', S*S-(C+1)**2*A*B)

c, D, N = row(33, 3)
x, C = 17, 9
h, delta = N-2*x, 2*x+C-N
A, B = D[x], D[x+C]
print('user residual example (r,a,x,C)=(3,33,17,9): N,h,delta=',
      N, h, delta, 'D endpoints=', A, B,
      'ratio^2=', Fraction(A, B),
      'plateau test difference=', (h+1)**2*A-delta**2*B)

c, D, N = row(11, 3)
x, C = 8, 3
j = 10
S = sum(D[x:x+C+1])
print('e=3,a=11 consumer window x=8,C=3: N=', N, 'coefficients=', c)
print('Dwindow=', D[x:x+C+1],
      'local log-concavity difference=', D[j]**2-D[j-1]*D[j+1])
print('endpoint ratio squared=', Fraction(D[x], D[x+C]), 'S=', S,
      'target-square margin=',
      S*S-(C+1)**2*D[x]*D[x+C])
