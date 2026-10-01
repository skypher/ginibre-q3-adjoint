import argparse, ast
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path

argparse.ArgumentParser(description='Exact FM-SEC136 grid verifier').parse_args()
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
            out[i + j] += x * y
    return out

def even(P):
    return [P[2*i] if 2*i < len(P) else Q(0)
            for i in range((len(P) + 1)//2)]

def coeff(P, i):
    return P[i] if i < len(P) else Q(0)

def gram_family(T, b, d, beta):
    M = T - 2*b
    f = fpoly(T, b, d)
    c = craw(M, d)
    sq = [even(mul(c[j], c[j])) for j in range(d + 1)]
    cross = [even(mul(c[j], c[j+2])) for j in range(d - 1)]
    A = [Q(0)] * (d + 1)
    for m in range(d, -1, -1):
        rhs = f[m]
        rhs -= sum(A[k] * coeff(sq[k], m) for k in range(m + 1, d + 1))
        rhs -= 2 * sum(beta[j] * coeff(cross[j], m)
                       for j in range(d - 1))
        A[m] = rhs / coeff(sq[m], m)
    reconstructed = [
        sum(A[j] * coeff(sq[j], m) for j in range(d + 1))
        + 2*sum(beta[j] * coeff(cross[j], m) for j in range(d - 1))
        for m in range(d + 1)
    ]
    assert reconstructed == f
    return A

def diagonal_member(T, b, d):
    f = fpoly(T, b, d)
    c = craw(T - 2*b, d)
    sq = [even(mul(c[j], c[j])) for j in range(d + 1)]
    A = [Q(0)] * (d + 1)
    for m in range(d, -1, -1):
        A[m] = (f[m] - sum(A[k] * coeff(sq[k], m)
                            for k in range(m + 1, d + 1))) / coeff(sq[m], m)
    assert [sum(A[j] * coeff(sq[j], m) for j in range(d + 1))
            for m in range(d + 1)] == f
    return A

def evalpoly(P, u):
    value = Q(0)
    for x in reversed(P):
        value = value*u + x
    return value

rows = negative_diagonal = sign_checks = sign_failures = 0
first = None
print('grid_start_utc=', datetime.now(timezone.utc).isoformat(timespec='seconds'),
      flush=True)
for T in range(4, 61):
    for b in range(3, 21):
        for d in range(6, 21):
            if T + 2*b - 2*d < 7:
                continue
            A = diagonal_member(T, b, d)
            rows += 1
            if any(v < 0 for v in A):
                negative_diagonal += 1
                if first is None:
                    first = (T, b, d, T - 2*b, tuple(map(str, A)))
            F = fpoly(T, b, d)
            for xi in range(61):
                sign_checks += 1
                if evalpoly(F, xi*xi) < 0:
                    sign_failures += 1
            if rows % 500 == 0:
                print('grid_progress rows=', rows,
                      'negative_diag=', negative_diagonal,
                      'last=', (T, b, d),
                      'utc=', datetime.now(timezone.utc).isoformat(timespec='seconds'),
                      flush=True)

assert rows == 12951 and negative_diagonal == 7660 and sign_checks == 790011
assert sign_failures == 0
assert first == (4, 8, 6, -12,
    ('56134/35', '-363/20', '-3396/35', '-411/140',
     '886/105', '83/42', '1/7'))
A0 = gram_family(4, 8, 6, [Q(0)]*5)
A1 = gram_family(4, 8, 6, [Q(1)] + [Q(0)]*4)
assert A0 != A1
assert evalpoly(fpoly(4, 8, 6), 0) == profile(2, 2, 8)[8] == 2280
print('rows=', rows, 'negative_diagonal_rows=', negative_diagonal)
print('first_diagonal_failure=', first)
print('integer_xi_checks=', sign_checks, 'negative_values=', sign_failures)
print('nonunique_family_check=PASS; consumer_bridge=2280')
print('grid_end_utc=', datetime.now(timezone.utc).isoformat(timespec='seconds'))

