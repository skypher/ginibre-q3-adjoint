import ast
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path
from math import factorial

import numpy as np
from scipy.optimize import minimize

source = Path(
    "ginibre_q3/character_ring_iter/fm39/"
    "mech54_genfun_obstructions_repro.py"
)
lib = {}
nodes = [
    node for node in ast.parse(source.read_text()).body
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))
]
exec(
    compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"),
    lib,
)
fpoly, craw = lib["fpoly"], lib["craw"]


def mul(P, R):
    out = [Q(0)] * (len(P) + len(R) - 1)
    for i, x in enumerate(P):
        for j, y in enumerate(R):
            out[i + j] += x * y
    return out


def even(P):
    return [
        P[2 * i] if 2 * i < len(P) else Q(0)
        for i in range((len(P) + 1) // 2)
    ]


def coeff(P, i):
    return P[i] if i < len(P) else Q(0)


def gram(T, b, d, beta):
    c = craw(T - 2 * b, d)
    F = fpoly(T, b, d)
    squares = [even(mul(x, x)) for x in c]
    crosses = [even(mul(c[j], c[j + 2])) for j in range(d - 1)]

    A = [Q(0)] * (d + 1)
    for m in range(d, -1, -1):
        rhs = F[m]
        rhs -= sum(A[k] * coeff(squares[k], m)
                   for k in range(m + 1, d + 1))
        rhs -= 2 * sum(beta[j] * coeff(crosses[j], m)
                       for j in range(d - 1))
        A[m] = rhs / coeff(squares[m], m)

    reconstructed = [
        sum(A[j] * coeff(squares[j], m) for j in range(d + 1))
        + 2 * sum(beta[j] * coeff(crosses[j], m)
                  for j in range(d - 1))
        for m in range(d + 1)
    ]
    assert reconstructed == F

    G = [[Q(0) for _ in range(d + 1)] for _ in range(d + 1)]
    for j, value in enumerate(A):
        G[j][j] = value
    for j, value in enumerate(beta):
        G[j][j + 2] = G[j + 2][j] = value
    return A, G


def exact_positive_pivots(A, beta, d):
    pivots = []
    for parity in (0, 1):
        previous = None
        for j in range(parity, d + 1, 2):
            pivot = A[j]
            if previous is not None:
                if previous <= 0:
                    return False, pivots
                pivot -= beta[j - 2] ** 2 / previous
            if pivot <= 0:
                return False, pivots + [pivot]
            pivots.append(pivot)
            previous = pivot
    return True, pivots


def propose_beta(T, b, d):
    n = d - 1
    zero = [Q(0)] * n
    A0, G0 = gram(T, b, d, zero)

    def monic_matrix(G):
        return np.array([
            [float(G[i][j]) * factorial(i) * factorial(j)
             for j in range(d + 1)]
            for i in range(d + 1)
        ])

    H0 = monic_matrix(G0)
    directions = []
    for edge in range(n):
        unit = [Q(0)] * n
        unit[edge] = Q(1)
        _, G = gram(T, b, d, unit)
        delta = [
            [G[i][j] - G0[i][j] for j in range(d + 1)]
            for i in range(d + 1)
        ]
        directions.append(monic_matrix(delta))

    scale = max(
        np.linalg.norm(H0),
        *(np.linalg.norm(H) for H in directions),
        1.0,
    )
    H0 = H0 / scale
    directions = [H / scale for H in directions]
    steps = np.array([
        1.0 / max(np.linalg.norm(H), 1e-200)
        for H in directions
    ])
    directions = [H * step for H, step in zip(directions, steps)]

    even_indices = np.arange(0, d + 1, 2)
    odd_indices = np.arange(1, d + 1, 2)

    def blocks(z):
        H = H0 + sum(
            (z[j] * directions[j] for j in range(n)),
            np.zeros_like(H0),
        )
        return [
            H[np.ix_(even_indices, even_indices)],
            H[np.ix_(odd_indices, odd_indices)],
        ]

    def constraints(v):
        return np.concatenate([
            np.linalg.eigvalsh(block) - v[-1]
            for block in blocks(v[:n])
        ])

    result = minimize(
        lambda v: -v[-1],
        np.r_[np.zeros(n), 0.0],
        method="SLSQP",
        constraints=[{"type": "ineq", "fun": constraints}],
        options={"maxiter": 300, "ftol": 1e-10, "disp": False},
    )
    beta = [
        Q(str(float(result.x[j] * steps[j])))
        for j in range(n)
    ]
    return beta


rows = [
    (T, b, d)
    for T in range(4, 31)
    for b in range(3, 13)
    for d in range(6, 15)
    if T + 2 * b - 2 * d >= 7
]
assert len(rows) == 1655

certified = set()
diagonal_bad = 0
print("sweep_start_utc=", datetime.now(timezone.utc).isoformat(), flush=True)

for index, (T, b, d) in enumerate(rows, 1):
    zero = [Q(0)] * (d - 1)
    A0, _ = gram(T, b, d, zero)
    diagonal_bad += any(value < 0 for value in A0)

    beta = propose_beta(T, b, d)
    A, _ = gram(T, b, d, beta)
    ok, _ = exact_positive_pivots(A, beta, d)
    if ok:
        certified.add((T, b, d))

    if index % 100 == 0 or index == len(rows):
        print(
            "heartbeat",
            index, "/", len(rows),
            "exact_positive_pivots=", len(certified),
            "utc=", datetime.now(timezone.utc).isoformat(),
            flush=True,
        )

extra = {
    (4, 12, 10): [
        Q("-26496107/10"), Q("-2696191/8"), Q("-468475/9"),
        Q("-105898/9"), Q("-21061/6"), Q("-7208/9"),
        Q("-488/5"), Q("-44/5"), Q("-1/4"),
    ],
    (5, 11, 10): [
        Q("-1482310/3"), Q("-478074/7"), Q("-8459"),
        Q("-8432/3"), Q("-12671/8"), Q("-3338/7"),
        Q("-723/10"), Q("-7"), Q("-1/5"),
    ],
    (5, 12, 10): [
        Q("-15243737/10"), Q("-1824008/9"), Q("-129605/4"),
        Q("-59629/7"), Q("-26780/9"), Q("-2827/4"),
        Q("-82"), Q("-32/5"), Q("-2/9"),
    ],
}
for row, beta in extra.items():
    A, _ = gram(*row, beta)
    ok, pivots = exact_positive_pivots(A, beta, row[2])
    assert ok
    certified.add(row)
    print("extra_exact_witness=", row, "pivot_count=", len(pivots))

unresolved = sorted(set(rows) - certified)
print("rows=", len(rows))
print("diagonal_member_bad_rows=", diagonal_bad)
print("exact_feasible_rows=", len(certified))
print("search_unresolved_rows=", len(unresolved))
print("first_search_unresolved=", unresolved[0] if unresolved else None)
if unresolved:
    T, b, d = unresolved[0]
    print("first_unresolved_F_coefficients=", list(map(str, fpoly(T, b, d))))
print("unresolved_rows=", unresolved)
print("sweep_end_utc=", datetime.now(timezone.utc).isoformat())

