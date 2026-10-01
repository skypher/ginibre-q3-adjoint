import argparse, ast
from fractions import Fraction as Q
from pathlib import Path

argparse.ArgumentParser(description='Exact mixed-Gram certificate verifier').parse_args()
path = Path('ginibre_q3/character_ring_iter/fm39/mech54_genfun_obstructions_repro.py')
ns = {}
nodes = [n for n in ast.parse(path.read_text()).body
          if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), ns)
fpoly, craw, profile = ns['fpoly'], ns['craw'], ns['profile']

def mul(P, R):
    out = [Q(0)] * (len(P) + len(R) - 1)
    for i, x in enumerate(P):
        for j, y in enumerate(R):
            out[i + j] += x*y
    return out

def even(P):
    return [P[2*i] if 2*i < len(P) else Q(0)
            for i in range((len(P) + 1)//2)]

def square(P):
    return even(mul(P, P))

def shifted_square(P, R, theta):
    v = list(P) + [Q(0)]*max(0, len(R)-len(P))
    for i, x in enumerate(R):
        v[i] += theta*x
    return square(v)

def add_scaled(dst, src, scale):
    for i, x in enumerate(src):
        dst[i] += scale*x

T, b, d = 4, 8, 6
M = T - 2*b
c = craw(M, d)
F = fpoly(T, b, d)
atoms = [
    (4, Q(91873, 26040), None),
    (5, Q(83, 42), None),
    (6, Q(1, 7), None),
    (0, Q(575104, 155), Q(-1, 8)),
    (1, Q(1551776, 3255), Q(-1, 8)),
    (1, Q(993424, 3255), Q(-1, 4)),
    (2, Q(7306, 93), Q(-1, 4)),
]
reconstructed = [Q(0)]*(d + 1)
diag = [Q(0)]*(d + 1)
off = [Q(0)]*(d - 1)
for j, weight, theta in atoms:
    if theta is None:
        atom = square(c[j])
        diag[j] += weight
    else:
        atom = shifted_square(c[j], c[j+2], theta)
        diag[j] += weight
        diag[j+2] += weight*theta*theta
        off[j] += weight*theta
    assert weight > 0
    add_scaled(reconstructed, atom, weight)
assert reconstructed == F
assert ns['peval'](F, 0) == profile(2, 2, 8)[8] == 2280

def pivots(parity):
    indices = list(range(parity, d + 1, 2))
    ans = [diag[indices[0]]]
    for j in indices[1:]:
        ans.append(diag[j] - off[j-2]*off[j-2]/ans[-1])
    return ans

even_pivots = pivots(0)
odd_pivots = pivots(1)
assert even_pivots == [Q(575104,155), Q(7306,93), Q(91873,26040), Q(1,7)]
assert odd_pivots == [Q(24240,31), Q(3010881877,1035578250), Q(83,42)]
assert all(x > 0 for x in even_pivots + odd_pivots)
print('T,b,d,M=', T, b, d, M)
print('F_coefficients=', list(map(str, F)))
print('Gram_diagonal=', list(map(str, diag)))
print('Gram_offdiagonal=', list(map(str, off)))
print('even_Schur_pivots=', list(map(str, even_pivots)))
print('odd_Schur_pivots=', list(map(str, odd_pivots)))
print('exact_sum_of_squares_and_consumer_bridge=PASS')

