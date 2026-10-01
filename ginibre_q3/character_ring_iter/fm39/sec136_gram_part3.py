import argparse, ast
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path

ap = argparse.ArgumentParser(description='Exact finite-slope Gram-menu verifier')
ap.add_argument('--full-sweep', action='store_true',
                help='run the 595-row census and expanded-slope check')
args = ap.parse_args()
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
            out[i+j] += x*y
    return out

def even(P):
    return [P[2*i] if 2*i < len(P) else Q(0)
            for i in range((len(P)+1)//2)]

def coeff(P, i):
    return P[i] if i < len(P) else Q(0)

def square(P):
    return even(mul(P, P))

def shifted_square(P, R, theta):
    v = list(P) + [Q(0)]*max(0, len(R)-len(P))
    for i, x in enumerate(R):
        v[i] += theta*x
    return square(v)

def menu_columns(T, b, d, exponent_limit):
    c = craw(T-2*b, d)
    cols = [square(v) for v in c]
    for j in range(d-1):
        for exponent in range(-exponent_limit, exponent_limit+1):
            q = Q(2**exponent) if exponent >= 0 else Q(1, 2**(-exponent))
            cols.append(shifted_square(c[j], c[j+2], q))
            cols.append(shifted_square(c[j], c[j+2], -q))
    return cols

def dot_poly(y, P):
    return sum(y[i]*coeff(P, i) for i in range(len(y)))

T, b, d, E = 4, 12, 10, 10
F = fpoly(T, b, d)
expected_F = list(map(Q, [
    '624514', '302858271017/436590', '5654853187772/35083125',
    '76503792694039/7072758000', '2353995177979/8840947500',
    '3839849737/1347192000', '16591876897/1131641280000',
    '34421071/905313024000', '42857/862202880000',
    '443/14485008384000', '1/144850083840000']))
assert F == expected_F and all(x > 0 for x in F)
assert ns['peval'](F, 0) == profile(2, 2, 12)[8] == 624514
columns = menu_columns(T, b, d, E)
y = list(map(Q, [
    '45315/2356009092926060468371456',
    '5035/294501136615757558546432',
    '-65455/147250568307878779273216',
    '669655/73625284153939389636608',
    '55091039/36812642076969694818304',
    '-1037158585/2300790129810605926144',
    '227112384503/2300790129810605926144',
    '-749217563575/35949845778290717596',
    '513660797165717/143799383113162870384',
    '2577589027235735/8987461444572679399', '-1']))
atom_pairings = [dot_poly(y, col) for col in columns]
dual_value = dot_poly(y, F)
assert min(atom_pairings) == 0
assert dual_value == Q(-1667447749697, 1133829375970666600403763200) < 0
c = craw(T-2*b, d)
a00 = dot_poly(y, square(c[0]))
b02 = dot_poly(y, even(mul(c[0], c[2])))
a22 = dot_poly(y, square(c[2]))
local_det = a00*a22 - b02*b02
assert local_det == Q(
    -3067498225,
    1387694711487569557826861382409066114635198889984
) < 0
print('first_finite_menu_miss=', (T,b,d), 'M=', T-2*b,
      'output_label=', T+2*b-2*d)
print('F_coefficients=', list(map(str, F)))
print('menu_columns=', len(columns), 'minimum_dual_pairing=', min(atom_pairings))
print('dual_target_pairing=', dual_value, 'edge0_determinant=', local_det)
print('finite_menu_Farkas_check=PASS; not an all-slope PSD obstruction')

def phase1(columns, rhs):
    m = len(rhs)
    n = len(columns)
    A = [[coeff(columns[j], i) for j in range(n)] for i in range(m)]
    rows = [list(row) for row in A]
    bvec = list(rhs)
    for i in range(m):
        if bvec[i] < 0:
            bvec[i] = -bvec[i]
            rows[i] = [-x for x in rows[i]]
        assert any(rows[i])
    tab = [rows[i] + [Q(int(i == j)) for j in range(m)] + [bvec[i]]
           for i in range(m)]
    basis = [n+i for i in range(m)]
    cost = [Q(0)]*n + [Q(-1)]*m
    pivots = 0
    while True:
        reduced = [
            cost[j] - sum(cost[basis[i]]*tab[i][j] for i in range(m))
            for j in range(n+m)
        ]
        enter = next((j for j, value in enumerate(reduced) if value > 0), None)
        if enter is None:
            break
        choices = [(tab[i][-1]/tab[i][enter], basis[i], i)
                   for i in range(m) if tab[i][enter] > 0]
        if not choices:
            return None, pivots
        _, _, leave = min(choices)
        pivot = tab[leave][enter]
        tab[leave] = [x/pivot for x in tab[leave]]
        for i in range(m):
            if i != leave and tab[i][enter]:
                factor = tab[i][enter]
                tab[i] = [x-factor*y for x, y in zip(tab[i], tab[leave])]
        basis[leave] = enter
        pivots += 1
    artificial_sum = sum(tab[i][-1] for i in range(m) if basis[i] >= n)
    if artificial_sum:
        return None, pivots
    solution = [Q(0)]*n
    for i, j in enumerate(basis):
        if j < n:
            solution[j] = tab[i][-1]
    assert all(x >= 0 for x in solution)
    assert [sum(solution[j]*A[i][j] for j in range(n))
            for i in range(m)] == list(rhs)
    return solution, pivots

if args.full_sweep:
    cases = [(TT, bb, dd) for TT in range(4, 21)
             for bb in range(3, 13) for dd in range(6, 11)
             if TT + 2*bb - 2*dd >= 7]
    feasible = total_pivots = 0
    first_bad = None
    print('menu_sweep_start rows=', len(cases), 'utc=',
          datetime.now(timezone.utc).isoformat(timespec='seconds'), flush=True)
    for index, case in enumerate(cases, 1):
        TT, bb, dd = case
        target = fpoly(TT, bb, dd)
        cols = menu_columns(TT, bb, dd, 10)
        solution, pivots = phase1(cols, target)
        total_pivots += pivots
        if solution is None:
            if first_bad is None:
                first_bad = case
        else:
            feasible += 1
        if index % 50 == 0:
            print('menu_progress rows=', index, '/', len(cases),
                  'feasible=', feasible, 'last=', case, 'utc=',
                  datetime.now(timezone.utc).isoformat(timespec='seconds'),
                  flush=True)
    assert len(cases) == 595 and feasible == 574 and first_bad == (4,12,10)
    assert total_pivots == 27883
    expanded, expanded_pivots = phase1(menu_columns(4,12,10,30), fpoly(4,12,10))
    assert expanded is None and expanded_pivots == 197
    print('menu_sweep_total=', len(cases), 'feasible=', feasible,
          'first_uncovered=', first_bad, 'total_pivots=', total_pivots)
    print('expanded_exponents=+-30 columns=',
          len(menu_columns(4,12,10,30)), 'feasible=', expanded is not None,
          'pivots=', expanded_pivots)
print('PASS')

