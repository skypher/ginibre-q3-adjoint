"""FM-SEC124 (luna_max_mars): explicit joint-moment expansion of M_(b+1) - Q_b; nonnegative H blocks; alternating outer sum; 945 checks."""
import argparse
from fractions import Fraction as Q
from math import comb, factorial

parser = argparse.ArgumentParser(
    description='Exact explicit-sum and inequality census for FM-SEC124.'
)
parser.parse_args()

def mu(m, k):
    return Q(
        2 * factorial(2*m) * factorial(2*m+1)
        * factorial(2*k) * factorial(2*k+1),
        factorial(m)**2 * factorial(k)**2
        * factorial(m+k+1) * factorial(m+k+2)
    )

def M(A, E, b):
    a, e = A//2, E//2
    out = Q(0)
    for u in range(b+1):
        outer = Q(comb(b,u) * (-2)**(b-u), 2**u)
        for v in range(u+1):
            out += outer * comb(u,v) * mu(a+v, e+u-v)
    return out

def Qmoment(A, E, b):
    a, e = A//2, E//2
    out = Q(0)
    for u in range(b+1):
        outer = Q(comb(b,u) * (-2)**(b-u), 2**u)
        for v in range(u+1):
            m, k = a+v, e+u-v
            out += outer * comb(u,v) * (
                mu(m+1,k) - mu(m,k+1)
            ) / 4
    return out

def H(m, k):
    return (mu(m+1,k) + 3*mu(m,k+1))/4 - 2*mu(m,k)

def Hclosed(m, k):
    F = m*m - m + 5*k*k + 7*k - 2*m*k
    return Q(2*F, (m+k+2)*(m+k+3)) * mu(m,k)

def explicit_sum(A, E, b):
    a, e = A//2, E//2
    out = Q(0)
    for u in range(b+1):
        outer = Q(comb(b,u) * (-2)**(b-u), 2**u)
        for v in range(u+1):
            out += outer * comb(u,v) * H(a+v, e+u-v)
    return out

for m in range(31):
    for k in range(31):
        assert H(m,k) == Hclosed(m,k)
        assert H(m,k) >= 0

count = 0
minimum = None
first_failure = None
for E in range(2,21,2):
    for A in range(E+2,21,2):
        prev = Q(0)
        for b in range(21):
            Mb = M(A,E,b)
            Qb = Qmoment(A,E,b)
            assert (A+E+2*b+6)*Qb == 4*(A-E)*Mb + 12*b*prev
            diff = M(A,E,b+1) - Qb
            assert explicit_sum(A,E,b) == diff
            count += 1
            if diff <= 0 and first_failure is None:
                first_failure = (A,E,b,diff)
            if minimum is None or diff < minimum[0]:
                minimum = (diff,(A,E,b))
            prev = Qb

assert count == 945
assert first_failure is None

# First termwise sign failure in the H-block expansion.
A, E, b = 4, 2, 1
a, e = A//2, E//2
term_u0 = -2*H(a,e)
term_u1 = (H(a,e+1) + H(a+1,e))/2
assert (term_u0,term_u1,term_u0+term_u1) == (
    Q(-8), Q(16), Q(8)
)

# The moment-only comparison for the surrogate sequence M_j = 2^j.
q = Q(0)
for j in range(5):
    q = (8*2**j + 12*j*q)/Q(12+2*j)
assert q == Q(3872,105)
assert Q(32)-q == Q(-512,105)

actual_b4 = M(4,2,5) - Qmoment(4,2,4)
assert actual_b4 == 100

print('exact target cases:', count, 'failures:', first_failure,
      'minimum:', minimum)
print('H-block mixed-sign expansion:', term_u0, term_u1,
      'total:', term_u0+term_u1)
print('surrogate margin at (4,2,4):', Q(32)-q)
print('actual semicircle margin at (4,2,4):', actual_b4)
print('H closed form and nonnegativity for 0<=m,k<=30: PASS')