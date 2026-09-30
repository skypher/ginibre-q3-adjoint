
from math import comb
from collections import Counter
import sympy as sp

# Symbolic identities for delta log-concavity at e=0 and e=1.
a, k, h, t = sp.symbols('a k h t')

g0 = lambda x: (a+1)*(a+2)*(2*x-a+1) / (
    (x+1)**2*(x+2)*(a-x+1)
)
ratio0 = sp.factor(
    (k/(a-k+1))**2 * ((a-k)/(k+1))**2
    * g0(k-1)*g0(k+1)/g0(k)**2
)
factor0 = 4*(a**2 + 2*a*h**2 + 2*a + 5*h**2 + 1) / (
    h**2*(a-h+5)*(a+h+5)
)
assert sp.factor((1-ratio0).subs(k, (a+h-1)/2)-factor0) == 0

H0 = lambda x: (a+1)/((x+1)*(a-x+1))
assert sp.factor(
    H0(k)*(1-((a-k)/(k+1))**2*H0(k+1)/H0(k))-g0(k)
) == 0

g1 = lambda x: (
    -(a+1)*(a+2)*(a-2*x)*(a-2*x-1)*(a-2*x+1)
    / ((x+1)**2*(x+2)*(a-x+1)**2*(a-x+2))
)
ratio1 = sp.factor(
    (k/(a-k+1))**2 * ((a-k)/(k+1))**2
    * g1(k-1)*g1(k+1)/g1(k)**2
)
P = lambda aa, hh: (
    3*aa**2*hh**2 + 6*aa**2*hh - 6*aa**2
    + 2*aa*hh**4 + 8*aa*hh**3 + 22*aa*hh**2
    + 28*aa*hh - 24*aa
    + 5*hh**4 + 20*hh**3 + 43*hh**2 + 46*hh - 18
)
factor1 = 4*P(a, h) / (
    h*(h+1)**2*(h+2)*(a-h+5)*(a+h+7)
)
assert sp.factor((1-ratio1).subs(k, (a+h+1)/2)-factor1) == 0

H1 = lambda x: (
    (a+1)*((a-2*x+1)**2+a+1)
    / ((x+1)*(a-x+1)**2*(a-x+2))
)
assert sp.factor(
    H1(k)*(1-((a-k)/(k+1))**2*H1(k+1)/H1(k))-g1(k)
) == 0

p_form = (
    (2*a+5)*t**4 + (16*a+40)*t**3
    + (3*a**2+58*a+133)*t**2
    + (12*a**2+104*a+212)*t
    + 3*a**2+36*a+96
)
assert sp.expand(P(a, 1+t)-p_form) == 0
assert sp.factor(
    (a*(a-1)/2)**2
    - ((a-3)*(a-2)*(a-1)*(a+2)/12)
) == (a-1)*(a**3+2*a-6)/6
print("symbolic e=0/e=1 delta identities and LC factors: PASS")


def row(a, e):
    N = a + e
    return [
        sum((-1)**u * comb(e, u) * comb(a, k-u)
            for u in range(e+1) if 0 <= k-u <= a)
        for k in range(N+1)
    ]


def at(c, k):
    return c[k] if 0 <= k < len(c) else 0


def ld_region(L, s, x, y):
    # Exact test of L >= s*(sqrt(x)+sqrt(y)) for nonnegative integers x,y.
    if L < 0:
        return False
    z = L*L - s*s*(x+y)
    return z >= 0 and z*z >= 4*s**4*x*y


def values(c, j, i):
    C = lambda k: at(c, k)
    D = lambda k: C(k)**2 - C(k-1)*C(k+1)
    delta = lambda k: D(k)-D(k+1)
    B = lambda k: C(k-1)+C(k+1)

    s = i-j
    L = D(j)-D(i)
    W = B(i)*C(j)-C(i)*B(j)
    u, v = delta(j), delta(i-1)
    g = u*v
    x, y = D(j)*D(i-1), D(j+1)*D(i)

    i_ok = g >= 0 and W*W <= s*s*g
    ii_ok = L >= 0 and L*L >= s*s*g
    E_ok = L >= abs(W)
    ld_ok = ld_region(L, s, x, y)
    return L, W, s, u, v, g, x, y, i_ok, ii_ok, E_ok, ld_ok


counts = Counter()
open_counts = Counter()
odd_open_counts = Counter()
failures = []
odd_failures = []
open_residual = []
odd_open_residual = []
open_i_failures = []
span_failures = Counter()
odd_open_spans = Counter()
q1_count = q1_failures = 0
center_windows = center_failures = 0
odd_center_windows = odd_center_failures = 0
support_edge_ii_failures = []
support_edge_E_failures = 0
min_odd_residual_slack = None

four = {(4, 6, 5, 8), (4, 7, 6, 9),
        (4, 8, 6, 9), (4, 9, 7, 10)}
four_data = []

for A in range(41):
    for e in range(41):
        N = A+e
        c = row(A, e)
        open_region = min(A, e) >= 3 and abs(A-e) >= 2

        for j in range((N+1)//2, N+1):
            for i in range(j+2, N+2):
                L,W,s,u,v,g,x,y,i_ok,ii_ok,E_ok,ld_ok = values(c,j,i)
                assert u >= 0 and v >= 0

                counts["pairs"] += 1
                counts["i pass"] += i_ok
                counts["ii pass"] += ii_ok
                counts["both pass"] += i_ok and ii_ok
                counts["E fail"] += not E_ok

                if not ii_ok:
                    failures.append((A,e,j,i))
                    span_failures[s] += 1
                    if e % 2:
                        odd_failures.append((A,e,j,i))

                if s == 2:
                    q1_count += 1
                    q1_failures += not ii_ok

                if N % 2 == 0 and j == N//2:
                    center_windows += 1
                    center_failures += not ii_ok
                    if e % 2:
                        odd_center_windows += 1
                        odd_center_failures += not ii_ok
                        assert u == 0 and ii_ok

                if i == N+1:
                    if not ii_ok:
                        support_edge_ii_failures.append((A,e,j,i))
                    support_edge_E_failures += not E_ok

                if open_region:
                    open_counts["pairs"] += 1
                    open_counts["i pass"] += i_ok
                    open_counts["ii pass"] += ii_ok
                    open_counts["both pass"] += i_ok and ii_ok
                    if not i_ok:
                        open_i_failures.append((A,e,j,i))
                    if not (i_ok and ii_ok):
                        open_residual.append((A,e,j,i))

                    if e % 2:
                        odd_open_counts["pairs"] += 1
                        odd_open_counts["i pass"] += i_ok
                        odd_open_counts["ii pass"] += ii_ok
                        odd_open_counts["both pass"] += i_ok and ii_ok

                        if not (i_ok and ii_ok):
                            odd_open_residual.append((A,e,j,i))
                            odd_open_spans[s] += 1
                            assert i_ok and not ld_ok
                            slack = L-abs(W)
                            candidate = (slack,A,e,j,i,L,W,u,v)
                            if (min_odd_residual_slack is None
                                    or slack < min_odd_residual_slack[0]):
                                min_odd_residual_slack = candidate

                if (A,e,j,i) in four:
                    four_data.append(
                        ((A,e,j,i),L,W,s,u,v,g,L*L-s*s*g,E_ok)
                    )

assert counts["E fail"] == 0
assert q1_failures == 0
assert odd_center_failures == 0
assert len(open_i_failures) == 6
assert support_edge_E_failures == 0

print("all q>=1 counts:", dict(counts))
print("all-parity ii-failure spans:", dict(span_failures))
print("open all-parity:", dict(open_counts),
      "union residual:", len(open_residual))
print("open odd-e:", dict(odd_open_counts),
      "residual spans:", dict(odd_open_spans),
      "residual:", len(odd_open_residual))
print("odd-e ii failures before the known E strips:", len(odd_failures))
print("q=1 checks/failures:", q1_count, q1_failures)
print("centre-start checks/failures:", center_windows, center_failures)
print("odd-e centre-start checks/failures:",
      odd_center_windows, odd_center_failures)
print("open i-psi failures:", open_i_failures)
print("support-edge ii failures:", support_edge_ii_failures)
print("support-edge E failures:", support_edge_E_failures)
print("least E slack on odd-open residual:", min_odd_residual_slack)
print("four listed witnesses:", four_data)

# These are the full finite sets, generated by the exact integer tests above.
print("F40 =", failures)
print("odd-open consumer residual =", odd_open_residual)
