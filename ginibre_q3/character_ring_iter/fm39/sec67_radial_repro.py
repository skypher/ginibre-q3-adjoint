# FM-SEC67 (luna_max_uranus) printed code: sharp radial cutoffs for three cores; parity images of the 31 residual words.
import sympy as sp
from pathlib import Path

root = Path("/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39")
me = {}
exec((root / "mech28_eval.py").read_text(), me)

s, d, R, n, xi = sp.symbols("s d R n xi", integer=True)

def rising(x, k):
    out = sp.Integer(1)
    for j in range(k):
        out *= x + j
    return sp.expand(out)

def core_poly(hs, ss):
    f = me["word"](hs, ss)
    expr = sum(v * ((s+d)/2)**i * ((s-d)/2)**j
               for (i, j), v in f.items())
    return sp.Poly(sp.expand(expr), s, d)

def qrad(hs, ss, Rv, nv):
    poly = core_poly(hs, ss)
    k0 = min((i+1)//2 for (i, j), _ in poly.terms())
    assert k0 == 2
    m0 = nv + k0
    total = 0
    for (i, j), b in poly.terms():
        assert i % 2 == 1 and j % 2 == 0
        alpha = (i+1)//2 - k0
        jj = j//2
        ratio = (
            rising(2*m0+1, 2*alpha)
            * rising(2*m0+2, 2*alpha)
            / rising(m0+1, alpha)**2
        )
        ratio *= (
            rising(2*Rv+1, 2*jj)
            * rising(2*Rv+2, 2*jj)
            / rising(Rv+1, jj)**2
        )
        ratio /= (
            rising(m0+Rv+2, alpha+jj)
            * rising(m0+Rv+3, alpha+jj)
        )
        total += b * ratio
    return sp.factor(total)

families = [
    ("hatS3^3", (), (3,3,3), 1),
    ("h3^3", (3,3,3), (), 0),
    ("h5hatS3^2", (5,), (3,3), 0),
]
for name, hs, ss, shift in families:
    Q = qrad(hs, ss, R, n)
    expr = sp.factor(Q.subs(n, R + shift + xi))
    num, den = sp.fraction(expr)
    num_poly = sp.Poly(sp.expand(num), xi)
    den_poly = sp.Poly(sp.expand(den), R, xi)
    assert all(c >= 0 for c in den_poly.coeffs())  # main agent: coeffs() for the multivariate denominator
    for ((k,), coeff) in num_poly.terms():
        assert all(c > 0 for c in sp.Poly(coeff, R).all_coeffs())
    print(name, "denominator =", sp.factor(den))
    print(name, "numerator coefficients =", [
        sp.factor(num_poly.nth(k)) for k in range(num_poly.degree()+1)
    ])

for name, ss in [
    ("hatS5^3", (5,5,5)),
    ("hatS5hatS3^2", (5,3,3)),
]:
    Q = qrad((), ss, sp.Integer(1), n)
    num, den = sp.fraction(sp.factor(Q))
    content, primitive = sp.Poly(sp.expand(num), n).primitive()
    assert all(c > 0 for c in primitive.all_coeffs())
    print(name, "content =", content)
    print(name, "primitive numerator =", primitive.as_expr())
    print(name, "denominator =", sp.factor(den))

T = sp.symbols("T")
for name, hs, ss in [
    ("hatS3^3", (), (3,3,3)),
    ("h3^3", (3,3,3), ()),
    ("h5hatS3^2", (5,), (3,3)),
    ("hatS5^3", (), (5,5,5)),
    ("hatS5hatS3^2", (), (5,3,3)),
]:
    jet = sp.Poly(sp.expand(core_poly(hs, ss).as_expr().subs(s, 4-T)), T, d)
    terms = {(i,j): c for (i,j), c in jet.terms() if i+j <= 2}
    print("jet", name, terms)

e = sp.symbols("e", positive=True)
beta_mean = 4*(2*R+3)/(e+2*R+4)
beta_second = 16*(2*R+3)*(2*R+4)/((e+2*R+4)*(e+2*R+5))
print("Beta bounds:", beta_mean, beta_second)


# ---- part 2 ----
from pathlib import Path
from collections import defaultdict

root = Path("/home/yang/q3adjoint/ginibre_q3/character_ring_iter/fm39")
me = {}
exec((root / "mech28_eval.py").read_text(), me)
direct = {}
exec((root / "split3.py").read_text(), direct)

failures = [
    (17,4,3,0,0,9), (17,31,4,0,0,22),
    (19,4,3,0,0,10), (19,16,3,0,0,16),
    (19,31,4,0,0,23), (19,33,4,0,0,24), (19,35,4,0,0,25),
    (21,4,3,0,0,11), (21,5,3,3,0,10), (21,6,3,0,0,12),
    (21,18,3,0,0,18), (21,33,4,0,0,25), (21,35,4,0,0,26),
    (21,37,4,0,0,27), (21,39,4,0,0,28),
    (23,4,3,0,0,12), (23,5,3,3,0,11), (23,6,3,0,0,13),
    (23,18,3,0,0,19), (23,20,3,0,0,20), (23,35,4,0,0,27),
    (23,37,4,0,0,28), (23,39,4,0,0,29),
    (25,4,3,0,0,13), (25,5,3,3,0,12), (25,6,3,0,0,14),
    (25,8,3,0,0,15), (25,20,3,0,0,21), (25,22,3,0,0,22),
    (25,37,4,0,0,29), (25,39,4,0,0,30),
]

groups = defaultdict(list)
for e, a, C, p, q, suffix in failures:
    ks = (C+p-1, C+q-1, C-1)
    odd_count = sum(k % 2 for k in ks)
    assert (a + odd_count) % 2 == 0
    R0 = (a + odd_count)//2
    hs = tuple(k for k in ks if k % 2)
    hats = tuple(k+1 for k in ks if k % 2 == 0)
    transformed = me["phi"](
        R0,
        me["mul"](me["word"](hs, hats),
                  me["power"](me["h"](1), e))
    )
    original = direct["phi3"]((e+3)//2, a, *ks)
    assert transformed == original and original > 0
    groups[(hs, hats, R0)].append((e, int(original)))

assert len(failures) == 31
for key, rows in groups.items():
    print(key, rows)

for a in range(2, 9):
    vals = []
    for r in range(5, 21):
        e = 2*r - 3
        val = direct["phi3"](r, a, 2, 2, 2)
        if a % 2:
            assert val == 0
        else:
            image = me["phi"](
                a//2,
                me["mul"](me["word"]((), (3,3,3)),
                          me["power"](me["h"](1), e))
            )
            assert image == val and val > 0
        vals.append(int(val))
    print("a, min, max, zeros =", a, min(vals), max(vals),
          sum(v == 0 for v in vals))

low_a2 = [
    (11,2,3,0,0,5), (13,2,3,0,0,6),
    (21,2,3,2,0,9), (25,2,5,0,0,11),
]
for e, a, C, p, q, suffix in low_a2:
    ks = (C+p-1, C+q-1, C-1)
    R0 = (a + sum(k % 2 for k in ks))//2
    hs = tuple(k for k in ks if k % 2)
    hats = tuple(k+1 for k in ks if k % 2 == 0)
    val = me["phi"](
        R0, me["mul"](me["word"](hs, hats),
                      me["power"](me["h"](1), e))
    )
    assert val == direct["phi3"]((e+3)//2, a, *ks)
    print("low-a case:", (e,a,C,p,q,suffix), "R =", R0,
          "core =", (hs,hats), "value =", int(val))
