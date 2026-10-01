"""FM-MECH86: corrected energy and a uniform one-label band."""
import argparse
from functools import lru_cache
from fractions import Fraction as Q
from math import comb
import sympy as sp

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--box", type=int, default=16)
ap.add_argument("--max-b", type=int, default=10)
ap.add_argument("--family", type=int, default=30)
args = ap.parse_args()

def choose(n, j):
    return comb(n, j) if 0 <= j <= n else 0

@lru_cache(None)
def row(a, e):
    N = a + e
    c = [1]
    for k in range(N):
        z, rem = divmod(
            (a-e)*c[k] - (N-k+1)*(c[k-1] if k else 0),
            k+1)
        assert rem == 0
        c.append(z)
    return tuple(c)

def data(a, e):
    N = a + e
    cc = row(a, e)

    def c(k):
        return cc[k] if 0 <= k <= N else 0

    @lru_cache(None)
    def v(h, k):
        if h < 0:
            return 0
        return sum(comb(h, l)*c(k+h-2*l) for l in range(h+1))

    @lru_cache(None)
    def norm(B, k):
        return sum(comb(B, h)*2**(B-h)*v(h, k)**2
                   for h in range(B+1))

    @lru_cache(None)
    def T(B, k):
        if B < 0:
            return 0
        return sum(
            comb(B, h)*2**(B-h)*
            (v(h, k)**2-v(h, k-1)*v(h, k+1))
            for h in range(B+1))

    def F(B, k):
        return T(B, k)-T(B, k+1)

    def W(B, j, i):
        return sum(
            comb(B, h)*2**(B-h)*
            (v(h, j)*v(h+1, i)-v(h+1, j)*v(h, i))
            for h in range(B+1))

    def H(B, k):
        S = N+2*B+2
        return sum(
            comb(B, h)*2**(B-h)*(
                S*(v(h, k)**2+v(h, k-1)**2)
                -2*v(h, k-1)*(
                    (a-e)*v(h, k)+(B-h)*v(h+1, k)
                    +4*h*v(h-1, k)))
            for h in range(B+1))

    def energy(B, k):
        return H(B, k)+4*B*T(B-1, k)

    def sos(B, k):
        base = sum(
            comb(B, h)*2**(B-h)*(
                (e+1)*(v(h, k)+v(h, k-1))**2
                +(a+1)*(v(h, k)-v(h, k-1))**2)
            for h in range(B+1))
        twice = 0
        for h in range(B):
            weight = comb(B, h)*2**(B-h)*(B-h)
            plus = (2*(v(h, k)+v(h, k-1))
                    -v(h+1, k)-v(h+1, k-1))
            minus = (2*(v(h, k)-v(h, k-1))
                     +v(h+1, k)-v(h+1, k-1))
            twice += weight*(plus*plus+minus*minus)
        assert twice % 2 == 0
        return base+twice//2

    def dot_A(B, k, l):
        return sum(
            comb(B, h)*2**(B-h)*v(h, k)*(
                (a-e)*v(h, l)+(B-h)*v(h+1, l)
                +4*h*v(h-1, l))
            for h in range(B+1))

    return N, v, norm, T, F, W, H, energy, sos, dot_A

def band(a, e, B, p):
    S = a+e+2*B+2
    if not 0 < p < S:
        return False
    R = p*(p+2)*(S-p)*(S+p+2)
    A = abs(a-e)*(p+1)
    C = B*(S+3*p+4)
    Z = 2*R-2*A*A-C*C
    return Z >= 0 and Z*Z >= 8*A*A*C*C

# Polynomial identities valid for symbolic parameters.
p, z = sp.symbols("p z")
S = 2*p+z
R = p*(p+2)*(S-p)*(S+p+2)
lhs = 7056*R-((41*p+21)*S+60*p*p+80*p)**2
rhs = (
    p*p*(1004*p*p+21800*p+13340)
    +z*p*(16580*p*p+54592*p+23100)
    +z*z*(5375*p*p+12390*p-441))
assert sp.expand(lhs-rhs) == 0

x, y, u, w = sp.symbols("x y u w")
assert sp.expand(
    ((2*(x+y)-(u+w))**2+(2*(x-y)+(u-w))**2)/2
    -(4*x*x+4*y*y+u*u+w*w-4*y*u-4*x*w)) == 0
print("SYMBOLIC CERTIFICATES PASS", flush=True)

rec = energies = sos_checks = local = insert = 0
nband = rectangle = nontrivial = 0
for a in range(args.box+1):
    for e in range(a+1):
        N, v, norm, T, F, W, H, energy, sos, dot_A = data(a, e)
        for B in range(args.max_b+1):
            for k in range((N+1)//2, N+B+2):
                S = N+2*B+2
                p = 2*k-N
                t = k+B+1
                u = N+B-k+1
                for h in range(B+1):
                    assert t*v(h, k+1) == (
                        (a-e)*v(h, k)-u*v(h, k-1)
                        +(B-h)*v(h+1, k)+4*h*v(h-1, k))
                    rec += 1

                square = p*sum(
                    comb(B, h)*2**(B-h)*
                    (v(h, k-1)-v(h, k+1))**2
                    for h in range(B+1))
                assert energy(B, k)-energy(B, k+1) == square
                prev = F(B-1, k) if B else 0
                assert H(B, k)-H(B, k+1) == square-4*B*prev
                energies += 1

                assert energy(B, k) == sos(B, k)
                assert energy(B, k) >= (
                    2*(e+1)*(norm(B, k)+norm(B, k-1)))
                sos_checks += 1

                if B:
                    assert F(B, k) == (
                        T(B-1, k-1)-T(B-1, k+2)
                        +W(B-1, k-1, k+2))
                    insert += 1

                if 0 < p < S:
                    left = dot_A(B, k, k+1)
                    right = dot_A(B, k+1, k)
                    qform = (
                        u*(p+2)*norm(B, k)
                        +(t+1)*p*norm(B, k+1)
                        +u*left-(t+1)*right)
                    assert u*(t+1)*F(B, k) == qform
                    assert left-right == -2*B*prev
                    assert 2*u*(t+1)*F(B, k) == (
                        2*u*(p+2)*norm(B, k)
                        +2*(t+1)*p*norm(B, k+1)
                        -(p+1)*(left+right)
                        -2*B*(S+1)*prev)
                    local += 1

                    rect = (
                        B >= 1 and p >= 3*B and 2*p <= S
                        and 4*(a-e) <= S)
                    if rect:
                        assert band(a, e, B, p)
                        rectangle += 1
                    if band(a, e, B, p):
                        assert F(B, k) >= 0
                        nband += 1
                        if (a >= 2 and e >= 2 and B >= 3
                                and p >= 11 and N+2*B-p >= 16):
                            nontrivial += 1

print("RECURRENCES", rec, "ENERGIES", energies,
      "SOS", sos_checks, flush=True)
print("LOCAL FORMS", local, "INSERTIONS", insert, flush=True)
print("BAND", nband, "RECTANGLE", rectangle,
      "B>=3, p>=11, distance>=8", nontrivial, flush=True)

# Independent Catalan moments followed by ballot conversion.
@lru_cache(None)
def cat(j):
    return comb(2*j, j)//(j+1)

def direct(a, e, B):
    N = a+e
    degree = N+2*B
    mon = [0]*(degree+1)
    cc = row(a, e)
    for j in range(0, N+1, 2):
        for l in range(B+1):
            for h in range(l+1):
                mon[N-j+2*(l-h)] += (
                    cc[j]*comb(B, l)*(-2)**(B-l)
                    *comb(l, h)*cat(j//2+h))
    return [
        sum(mon[d]*(
            choose(d, (d-p)//2)-choose(d, (d-p)//2-1))
            for d in range(p, degree+1, 2))
        for p in range(degree+1)]

bridges = 0
for a in range(7):
    for e in range(7):
        N, v, norm, T, F, W, H, energy, sos, dot_A = data(a, e)
        for B in range(5):
            for p, val in enumerate(direct(a, e, B)):
                if (N+p) % 2:
                    assert val == 0
                else:
                    assert val == F(B, (N+p)//2)
                bridges += 1
print("CATALAN BRIDGES", bridges, flush=True)

# Bounded vetting of the proved uniform family.
for B in range(1, args.family+1):
    a = e = 4*B
    p = 4*B
    N, v, norm, T, F, W, H, energy, sos, dot_A = data(a, e)
    assert 3*B <= p and 2*p <= N+2*B+2
    assert band(a, e, B, p)
    assert F(B, (N+p)//2) > 0
print("UNIFORM FAMILY", args.family, "PASS", flush=True)

# Exact obstruction to positivity on arbitrary vector states.
N, B, p, h = 100, 32, 12, 16
k = (N+p)//2
t = k+B+1
u = N+B-k+1
weight = lambda h: comb(B, h)*2**(B-h)
av = weight(h)*Q(p+2, t+1)
bv = weight(h+1)*Q(p, u)
cv = weight(h)*(B-h)*(Q(1, t+1)-Q(2, u))
tau = -cv/(2*bv)
bad = av+cv*tau+bv*tau*tau
assert tau == Q(1751, 1080)
assert bad == Q(-15241796664229888, 10395)
assert Q(B-h+1, u) == Q(17, 77) != 1
print("UNRESTRICTED-VECTOR CONTROL", tau, bad,
      "compatibility 17/77 != 1", flush=True)

N, v, norm, T, F, W, H, energy, sos, dot_A = data(50, 50)
actual = F(32, 56)
assert actual == 17857708380871947518719985545060893888
print("ACTUAL ROW AT (50,50,32,12)", actual, flush=True)
print("PASS", flush=True)