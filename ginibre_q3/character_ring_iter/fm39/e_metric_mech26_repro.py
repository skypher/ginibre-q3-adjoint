
import argparse
from fractions import Fraction
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--max", type=int, default=40)
ap.add_argument("--q1-max", type=int, default=50)
args = ap.parse_args()

def row(a, e):
    N, d = a+e, a-e
    c, prev = [1], 0
    for k in range(N):
        z = d*c[-1]-(N-k+1)*prev
        assert z % (k+1) == 0
        prev, nxt = c[-1], z//(k+1)
        c.append(nxt)
    return c

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def D(c, k):
    return at(c, k)**2-at(c, k-1)*at(c, k+1)

def B(c, k):
    return at(c, k-1)+at(c, k+1)

def A(c, k):
    return at(c, k-1)-at(c, k+1)

def delta(c, k):
    return D(c, k)-D(c, k+1)

def ld(E, s, p, q):
    z = E*E-s*s*(p+q)
    return z >= 0 and z*z >= 4*s**4*p*q

counts, left = Counter(), []
for a in range(args.max+1):
    for e in range(args.max+1):
        c = row(a, e)
        N, V = a+e, (a+1)*(e+1)
        ds = [D(c, k) for k in range(N+2)]
        bs = [B(c, k) for k in range(N+2)]
        aa = [A(c, k) for k in range(N+2)]
        for j in range((N+1)//2, N+1):
            for i in range(j+1, N+2):
                s, x, y = i-j, 2*j-N, 2*i-N
                P, Q = ds[j], ds[i]
                E = P-Q
                W = at(c, j)*bs[i]-bs[j]*at(c, i)
                assert E >= abs(W)
                counts["pairs"] += 1
                old = (
                    abs(a-e) <= 1 or min(a, e) <= 2 or s == 1
                    or s == 2 and abs(a-e) in (2, 3)
                    or (j, i) == (N-1, N+1)
                    or ld(E, s, P*ds[i-1], ds[j+1]*Q)
                )
                mi = mj = False
                if y*y < 4*V:
                    rhs = Q*(4*y*y*P-(y*y-x*x)*aa[j]**2)
                    assert 0 <= rhs
                    assert (4*V-y*y)*W*W <= rhs
                    mi = (4*V-y*y)*E*E >= rhs
                if 0 < x*x < 4*V:
                    rhs = P*(4*x*x*Q+(y*y-x*x)*aa[i]**2)
                    assert (4*V-x*x)*W*W <= rhs
                    mj = (4*V-x*x)*E*E >= rhs
                if not old:
                    counts["old_residual"] += 1
                    if mi or mj:
                        counts["new_metric"] += 1
                    else:
                        left.append((a, e, j, i))

assert all(i-j == 2 for a, e, j, i in left)
print(dict(counts))
print("left for q=1:", left)

q1 = 0
for a in range(args.q1_max+1):
    for e in range(args.q1_max+1):
        c = row(a, e)
        cp, cm = row(a+1, e), row(a, e+1)
        N = a+e
        for j in range((N+1)//2, N):
            i = j+2
            E = D(c, j)-D(c, i)
            W = at(c, j)*B(c, i)-B(c, j)*at(c, i)
            assert E+W == delta(cp, j+1) >= 0
            assert E-W == delta(cm, j+1) >= 0
            q1 += 1
print("q=1 shifted-OL checks:", q1)

def half(v):
    return 0 if v[1] > 0 or v[1] == 0 and v[0] > 0 else 1

def lt(v, w):
    if half(v) != half(w):
        return half(v) < half(w)
    return v[0]*w[1]-v[1]*w[0] > 0

def short_arc(c, j, i):
    v = [(at(c, k), B(c, k)) for k in range(j, i+1)]
    if any(delta(c, k) <= 0 for k in range(j, i)):
        return False
    W = v[0][0]*v[-1][1]-v[0][1]*v[-1][0]
    turns = sum(lt(q, p) for p, q in zip(v, v[1:]))
    return W > 0 and turns == int(lt(v[-1], v[0]))

def K(l, X, n):
    p, q = 1, X
    if l == 0:
        return p
    for h in range(1, l):
        p, q = q, X*q-h*(n-h+1)*p
    return q

for a, e, j, i in [(296, 4, 169, 173), (297, 3, 164, 168)]:
    c = row(a, e)
    N, V = a+e, (a+1)*(e+1)
    x, y = 2*j-N, 2*i-N
    P, Q = D(c, j), D(c, i)
    E = P-Q
    W = at(c, j)*B(c, i)-B(c, j)*at(c, i)

    # R(t)=F_j(t)F_i(t)-4t(4V-t)E^2.
    U0, U1 = 4*P-A(c, j)**2, 4*Q-A(c, i)**2
    V0, V1 = x*x*A(c, j)**2, y*y*A(c, i)**2
    alpha = U0*U1+4*E*E
    beta = U0*V1+U1*V0-16*V*E*E
    gamma = V0*V1
    tstar = Fraction(-beta, 2*alpha)
    minimum = Fraction(4*alpha*gamma-beta*beta, 4*alpha)

    assert alpha > 0 and 0 < tstar < 4*V and minimum > 0
    assert alpha*tstar*tstar+beta*tstar+gamma == minimum
    assert K(e, x, N+2)*K(e, y, N+2) < 0
    assert short_arc(c, j, i) and E >= abs(W)
    print("norm test fails for every t; short arc proves E:",
          (a, e, j, i),
          "root endpoint values:", K(e, x, N+2), K(e, y, N+2))
