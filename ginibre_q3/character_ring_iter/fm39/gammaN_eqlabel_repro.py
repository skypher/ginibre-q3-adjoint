
from math import comb

def row(a, e):
    c = [0] * (a + e + 1)
    for h in range(e + 1):
        for k in range(a + 1):
            c[h + k] += (-1)**h * comb(e, h) * comb(a, k)
    return c

def at(c, k):
    return c[k] if 0 <= k < len(c) else 0

def g(c, k):
    return (at(c, k), at(c, k - 1))

def add(u, v):
    return (u[0] + v[0], u[1] + v[1])

def wedge(u, v):
    return u[0] * v[1] - u[1] * v[0]

def D(c, k):
    return at(c, k)**2 - at(c, k - 1) * at(c, k + 1)

def P(c, x, C):
    T = at(c, x) * at(c, x + C) - at(c, x - 1) * at(c, x + C + 1)
    return sum(D(c, k) for k in range(x, x + C + 1)) - T

def upoly(n):
    # U_0=1, U_1=x, U_(n+1)=x U_n-U_(n-1)
    if n == 0:
        return [1]
    if n == 1:
        return [0, 1]
    u0, u1 = [1], [0, 1]
    for _ in range(2, n + 1):
        u2 = [0] * (len(u1) + 1)
        for k, v in enumerate(u1):
            u2[k + 1] += v
        for k, v in enumerate(u0):
            u2[k] -= v
        while u2 and u2[-1] == 0:
            u2.pop()
        u0, u1 = u1, u2
    return u1

def bimul(f, h):
    out = {}
    for (i, j), x in f.items():
        for (k, l), y in h.items():
            out[i + k, j + l] = out.get((i + k, j + l), 0) + x * y
    return {ij: v for ij, v in out.items() if v}

def phi_definition(e, a, A, B, C):
    # Expand 1/2 E[(x-y)^e (x+y)^a prod(U_L(x)-U_L(y)).
    f = {(e-k, k): (-1)**k * comb(e, k) for k in range(e + 1)}
    f = bimul(f, {(a-k, k): comb(a, k) for k in range(a + 1)})
    for n in (A, B, C):
        u = upoly(n)
        diff = {}
        for k, v in enumerate(u):
            diff[k, 0] = diff.get((k, 0), 0) + v
            diff[0, k] = diff.get((0, k), 0) - v
        f = bimul(f, {ij: v for ij, v in diff.items() if v})

    def semicircle_moment(m):
        if m % 2:
            return 0
        t = m // 2
        return comb(2*t, t) // (t + 1)

    total = sum(v * semicircle_moment(i) * semicircle_moment(j)
                for (i, j), v in f.items())
    assert total % 2 == 0
    return total // 2

def low_R(c, a, e, p, q, s):
    N, d = a + e, a - e
    if s == 0:
        return at(c, p) + at(c, q)
    if s == 1:
        return d * (at(c, p + 1) + at(c, q + 1)) \
               - at(c, p) - at(c, q) - 1
    return at(c, 2) * (at(c, p + 2) + at(c, q + 2)) \
           - d * (at(c, p + 1) + at(c, q + 1)) - 1 - (d*d + N)//2

def report(e, a, A, B, C):
    N = a + e
    gamma = (N + A + B - C) // 2
    s = N - gamma
    alpha = (N + A - B - C) // 2
    tau = gamma + 1
    p, q = A - C, B - C
    c = row(a, e)

    U = add(g(c, alpha), g(c, tau - A - 1))
    V = add(g(c, tau), g(c, alpha - A - 1))
    X = wedge(U, V)
    drop = P(c, alpha, C) - P(c, alpha - A - 1, C)
    value = phi_definition(e, a, A, B, C)
    assert value == drop - X

    print((e, a, A, B, C),
          "s,alpha,p,q,P(alpha),P(mirror),drop,X,phi =",
          (s, alpha, p, q, P(c, alpha, C),
           P(c, alpha - A - 1, C), drop, X, value))
    if s <= 2:
        R = low_R(c, a, e, p, q, s)
        assert value == P(c, p + s, C) + R
        print("  R_s =", R)

report(3, 1, 4, 3, 3)
report(3, 5, 6, 4, 4)
for s in range(3):
    e, C, p, q = 3, 3, 0, 0
    N = C + p + q + 2*s
    a = N - e
    report(e, a, C + p, C + q, C)

# The combined-vector H-area certificate at the stated case.
e, a, C, p, q, s = 1, 23, 3, 1, 0, 10
N = C + p + q + 2*s
alpha, tau, A = p + s, N - s + 1, C + p
c = row(a, e)
U = add(g(c, alpha), g(c, q + s))
V = add(g(c, tau), g(c, s - C - 1))
drop = P(c, alpha, C) - P(c, s - C - 1, C)
def H(v):
    x, y = v
    return ((e + 1)*(x + y)**2 + (a + 1)*(x - y)**2) // 2
detH = (a + 1) * (e + 1)
margin = drop**2 * detH - H(U) * H(V)
assert phi_definition(e, a, C+p, C+q, C) == drop - wedge(U, V)
print("metric case:", drop, wedge(U, V), drop - wedge(U, V))
print("U,V,H(U),H(V),detH,margin:",
      U, V, H(U), H(V), detH, margin)
