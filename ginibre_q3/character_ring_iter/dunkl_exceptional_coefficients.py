#!/usr/bin/env python3
"""Legacy Python cross-check of the exceptional Dunkl Gaussian identity.

The active load-bearing sweep is
``dunkl_exceptional_coefficients_gmp.cpp`` (GMP rationals plus OpenMP).  This
implementation is retained for independent source-level comparison and is
not invoked by the clean-room replay.

For the Coxeter Gaussian with density

    exp(-|w|^2) * prod_{alpha > 0} (alpha.w)^2 dw,

the normalized expectation of a W-invariant polynomial p is

    (exp(Delta_D / 4) p)(0),

where Delta_D is the multiplicity-one Dunkl Laplacian.  On W-invariant
polynomials this operator is

    Delta p + 2 * sum_{alpha > 0} (alpha.grad p)/(alpha.w).

For E_6, E_7, and E_8, the remaining coefficient inputs are expectations of
R(w)|w|^(2k) and S(w), where

    R(w) = sum_{alpha > 0} (alpha.w)^6,
    S(w) = sum_{alpha > 0} (alpha.w)^8.

Since the input polynomials are homogeneous, the expectation is obtained by
applying Delta_D exactly degree/2 times.

The script first checks the normalization against F_4, where the existing
Wick computation gives <S(w)> = 6265350.
"""

from __future__ import annotations

import os
import sys
from fractions import Fraction
from math import factorial
from typing import Dict, Iterable, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from c2_exact import (
    e6_simple_roots,
    e7_simple_roots,
    f4_positive_roots,
    generate_roots_from_simple,
    positive_roots,
)
from c3_exact import ROOT_DATA, compute_c3, eta_1, moment_QK

Exponent = Tuple[int, ...]
Poly = Dict[Exponent, Fraction]
Vector = Tuple[Fraction, ...]


def clean(poly: Poly) -> Poly:
    return {k: v for k, v in poly.items() if v}


def add_scaled(target: Poly, source: Poly, scale: Fraction) -> None:
    if scale == 0:
        return
    for key, coeff in source.items():
        new = target.get(key, Fraction(0)) + scale * coeff
        if new:
            target[key] = new
        elif key in target:
            del target[key]


def poly_mul(p: Poly, q: Poly) -> Poly:
    result: Poly = {}
    for ka, ca in p.items():
        for kb, cb in q.items():
            key = tuple(a + b for a, b in zip(ka, kb))
            result[key] = result.get(key, Fraction(0)) + ca * cb
    return clean(result)


def linear_power(alpha: Vector, power: int) -> Poly:
    n = len(alpha)
    lin: Poly = {}
    for i, coeff in enumerate(alpha):
        if coeff:
            key = tuple(1 if j == i else 0 for j in range(n))
            lin[key] = coeff

    result: Poly = {tuple(0 for _ in range(n)): Fraction(1)}
    for _ in range(power):
        result = poly_mul(result, lin)
    return result


def root_power_sum(roots: Iterable[Vector], power: int) -> Poly:
    total: Poly = {}
    for alpha in roots:
        add_scaled(total, linear_power(alpha, power), Fraction(1))
    return clean(total)


def ordinary_laplacian(poly: Poly, n_vars: int) -> Poly:
    result: Poly = {}
    for key, coeff in poly.items():
        for i in range(n_vars):
            if key[i] < 2:
                continue
            out_key = list(key)
            out_key[i] -= 2
            scale = coeff * key[i] * (key[i] - 1)
            out_tuple = tuple(out_key)
            result[out_tuple] = result.get(out_tuple, Fraction(0)) + scale
    return clean(result)


def directional_derivative(poly: Poly, alpha: Vector) -> Poly:
    result: Poly = {}
    for key, coeff in poly.items():
        for i, alpha_i in enumerate(alpha):
            if not alpha_i or key[i] == 0:
                continue
            out_key = list(key)
            out_key[i] -= 1
            out_tuple = tuple(out_key)
            result[out_tuple] = (
                result.get(out_tuple, Fraction(0)) + coeff * key[i] * alpha_i
            )
    return clean(result)


def divide_by_linear(poly: Poly, alpha: Vector) -> Poly:
    """Divide a homogeneous polynomial by L = alpha.x.

    The caller only invokes this on derivatives of W-invariant polynomials,
    where divisibility by every root linear form is guaranteed.  A failure
    raises an ArithmeticError instead of silently losing a remainder.
    """
    if not poly:
        return {}

    pivot = next(i for i, coeff in enumerate(alpha) if coeff)
    pivot_coeff = alpha[pivot]
    n_vars = len(alpha)
    degree_in_pivot = max(key[pivot] for key in poly)
    if degree_in_pivot == 0:
        raise ArithmeticError("nonzero polynomial cannot divide by pivot-linear form")

    buckets: List[Poly] = [{} for _ in range(degree_in_pivot + 1)]
    for key, coeff in poly.items():
        pivot_degree = key[pivot]
        y_key = list(key)
        y_key[pivot] = 0
        y_tuple = tuple(y_key)
        buckets[pivot_degree][y_tuple] = buckets[pivot_degree].get(y_tuple, Fraction(0)) + coeff
    buckets = [clean(bucket) for bucket in buckets]

    def subtract(p: Poly, q: Poly) -> Poly:
        out = dict(p)
        add_scaled(out, q, Fraction(-1))
        return out

    def scale_poly(p: Poly, scale: Fraction) -> Poly:
        return clean({key: coeff * scale for key, coeff in p.items()})

    def multiply_by_nonpivot_linear(p: Poly) -> Poly:
        out: Poly = {}
        for key, coeff in p.items():
            for i, alpha_i in enumerate(alpha):
                if i == pivot or not alpha_i:
                    continue
                out_key = list(key)
                out_key[i] += 1
                out_tuple = tuple(out_key)
                out[out_tuple] = out.get(out_tuple, Fraction(0)) + coeff * alpha_i
        return clean(out)

    quotient_buckets: List[Poly] = [{} for _ in range(degree_in_pivot)]
    next_bq: Poly = {}
    for pivot_degree in range(degree_in_pivot, 0, -1):
        rhs = subtract(buckets[pivot_degree], next_bq)
        q_part = scale_poly(rhs, Fraction(1, 1) / pivot_coeff)
        quotient_buckets[pivot_degree - 1] = q_part
        next_bq = multiply_by_nonpivot_linear(q_part)

    remainder = subtract(buckets[0], next_bq)
    if remainder:
        lead = next(iter(remainder))
        raise ArithmeticError(
            f"nonzero remainder while dividing by root linear form: {lead}"
        )

    quotient: Poly = {}
    for pivot_degree, bucket in enumerate(quotient_buckets):
        for key, coeff in bucket.items():
            out_key = list(key)
            out_key[pivot] = pivot_degree
            out_tuple = tuple(out_key)
            quotient[out_tuple] = quotient.get(out_tuple, Fraction(0)) + coeff
    return clean(quotient)


def dunkl_laplacian(poly: Poly, roots: List[Vector], progress_label: str = "") -> Poly:
    n_vars = len(roots[0])
    result = ordinary_laplacian(poly, n_vars)
    for idx, alpha in enumerate(roots, start=1):
        deriv = directional_derivative(poly, alpha)
        quotient = divide_by_linear(deriv, alpha)
        add_scaled(result, quotient, Fraction(2))
        if progress_label and (idx % 5 == 0 or idx == len(roots)):
            print(
                f"{progress_label}: root {idx}/{len(roots)}, "
                f"current terms={len(result)}",
                flush=True,
            )
    return clean(result)


def dunkl_expectation(poly: Poly, roots: List[Vector], degree: int, label: str) -> Fraction:
    if degree % 2:
        return Fraction(0)
    current = poly
    steps = degree // 2
    for step in range(steps):
        print(
            f"{label}: applying Dunkl Laplacian {step + 1}/{steps}; "
            f"input terms={len(current)}",
            flush=True,
        )
        current = dunkl_laplacian(current, roots, progress_label=f"{label} Δ{step + 1}")
    zero = tuple(0 for _ in roots[0])
    constant = current.get(zero, Fraction(0))
    return constant / (Fraction(4) ** steps * factorial(steps))


def quadratic_form_poly(matrix: List[List[Fraction]]) -> Poly:
    n_vars = len(matrix)
    square: Poly = {}
    for i in range(n_vars):
        if matrix[i][i]:
            key = tuple(2 if j == i else 0 for j in range(n_vars))
            square[key] = square.get(key, Fraction(0)) + matrix[i][i]
        for j in range(i + 1, n_vars):
            coeff = matrix[i][j] + matrix[j][i]
            if coeff:
                key = tuple((1 if k == i else 0) + (1 if k == j else 0) for k in range(n_vars))
                square[key] = square.get(key, Fraction(0)) + coeff
    return clean(square)


def quadratic_power(matrix: List[List[Fraction]], half_power: int) -> Poly:
    n_vars = len(matrix)
    square = quadratic_form_poly(matrix)
    result: Poly = {tuple(0 for _ in range(n_vars)): Fraction(1)}
    for _ in range(half_power):
        result = poly_mul(result, square)
    return result


def identity_matrix(n_vars: int) -> List[List[Fraction]]:
    return [
        [Fraction(1 if i == j else 0) for j in range(n_vars)]
        for i in range(n_vars)
    ]


def invert_matrix(matrix: List[List[Fraction]]) -> List[List[Fraction]]:
    n = len(matrix)
    aug = [
        list(row) + [Fraction(1 if i == j else 0) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col]:
                pivot = row
                break
        if pivot is None:
            raise ArithmeticError("singular matrix")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(n):
            if row == col or not aug[row][col]:
                continue
            factor = aug[row][col]
            aug[row] = [x - factor * y for x, y in zip(aug[row], aug[col])]
    return [row[n:] for row in aug]


def projection_matrix(simple_roots: List[Vector]) -> List[List[Fraction]]:
    """Orthogonal projector onto the span of the simple roots in ambient coords."""
    rank = len(simple_roots)
    n_vars = len(simple_roots[0])
    gram = [
        [sum(simple_roots[i][k] * simple_roots[j][k] for k in range(n_vars))
         for j in range(rank)]
        for i in range(rank)
    ]
    gram_inv = invert_matrix(gram)
    projector = [[Fraction(0) for _ in range(n_vars)] for _ in range(n_vars)]
    for a in range(n_vars):
        for b in range(n_vars):
            projector[a][b] = sum(
                simple_roots[i][a] * gram_inv[i][j] * simple_roots[j][b]
                for i in range(rank)
                for j in range(rank)
            )
    return projector


def positive_half(full_roots: Iterable[Vector]) -> List[Vector]:
    pos = []
    for root in full_roots:
        for coord in root:
            if coord > 0:
                pos.append(root)
                break
            if coord < 0:
                break
    return pos


def e8_positive_roots() -> List[Vector]:
    roots: List[Vector] = []
    n = 8

    for i in range(n):
        for j in range(i + 1, n):
            for si in (Fraction(1), Fraction(-1)):
                for sj in (Fraction(1), Fraction(-1)):
                    root = [Fraction(0)] * n
                    root[i] = si
                    root[j] = sj
                    roots.append(tuple(root))

    half = Fraction(1, 2)
    for mask in range(1 << n):
        signs = []
        minus_count = 0
        for i in range(n):
            if mask & (1 << i):
                signs.append(-half)
                minus_count += 1
            else:
                signs.append(half)
        if minus_count % 2 == 0:
            roots.append(tuple(signs))

    pos = positive_half(roots)
    if len(roots) != 240 or len(pos) != 120:
        raise AssertionError(f"E8 root count mismatch: full={len(roots)}, pos={len(pos)}")
    for root in pos:
        norm_sq = sum(x * x for x in root)
        if norm_sq != 2:
            raise AssertionError(f"E8 root has norm^2 {norm_sq}: {root}")
    return pos


def group_roots_and_projector(group: str) -> Tuple[List[Vector], List[List[Fraction]]]:
    if group == "F4":
        return f4_positive_roots(), identity_matrix(4)
    if group == "E6":
        simples = e6_simple_roots()
        roots = positive_roots(generate_roots_from_simple(simples, 8))
        if len(roots) != 36:
            raise AssertionError(f"E6 positive-root count mismatch: {len(roots)}")
        if {sum(x * x for x in root) for root in roots} != {Fraction(2)}:
            raise AssertionError("E6 root-norm check failed")
        if sum(sum(x * x for x in root) for root in roots) != ROOT_DATA["E6"]["A"]:
            raise AssertionError("E6 A=sum |alpha|^2 check failed")
        return roots, projection_matrix(simples)
    if group == "E7":
        simples = e7_simple_roots()
        roots = positive_roots(generate_roots_from_simple(simples, 8))
        if len(roots) != 63:
            raise AssertionError(f"E7 positive-root count mismatch: {len(roots)}")
        if {sum(x * x for x in root) for root in roots} != {Fraction(2)}:
            raise AssertionError("E7 root-norm check failed")
        if sum(sum(x * x for x in root) for root in roots) != ROOT_DATA["E7"]["A"]:
            raise AssertionError("E7 A=sum |alpha|^2 check failed")
        return roots, projection_matrix(simples)
    if group == "E8":
        return e8_positive_roots(), identity_matrix(8)
    raise ValueError(group)


def gamma_moment_from_group(group: str, k: int) -> Fraction:
    a = Fraction(ROOT_DATA[group]["d"], 2)
    out = Fraction(1)
    for i in range(k):
        out *= a + i
    return out


def root_sum_with_q_power(
    roots: List[Vector],
    projector: List[List[Fraction]],
    root_power: int,
    q_power: int,
    label: str,
) -> Fraction:
    poly = root_power_sum(roots, root_power)
    if q_power:
        poly = poly_mul(poly, quadratic_power(projector, q_power))
    return dunkl_expectation(poly, roots, root_power + 2 * q_power, label)


def c2_from_R_moment(group: str, r_w: Fraction) -> Fraction:
    rd = ROOT_DATA[group]
    d = rd["d"]
    cp = Fraction(rd["A"], rd["r"])
    a = Fraction(d, 2)

    coef_barR = -Fraction(d * d, 360 * cp**3)
    coef_w6 = Fraction(5 * d, 12 * (d + 2)) - Fraction(1, 3) + Fraction(
        (d + 12) * d, 144 * (d + 2)
    )
    coef_w4 = Fraction(d * d * (d + 1), 288 * (d + 2))
    coef_w8 = Fraction((d + 12) ** 2, 288 * (d + 2) ** 2)
    return (
        coef_barR * r_w
        + coef_w6 * gamma_moment_from_group(group, 3)
        + coef_w4 * gamma_moment_from_group(group, 2)
        + coef_w8 * gamma_moment_from_group(group, 4)
    )


def c3_from_RS_moments(
    group: str,
    r_q0_u: Fraction,
    r_q1_u: Fraction,
    r_q2_u: Fraction,
    s_q0_u: Fraction,
) -> Fraction:
    d = ROOT_DATA[group]["d"]
    d_f = Fraction(d)

    def moment_P_Qk(k: int) -> Fraction:
        return Fraction(5, d + 2) * moment_QK(group, k + 2)

    def moment_Pj_Qk(j: int, k: int) -> Fraction:
        return Fraction(5, d + 2) ** j * moment_QK(group, 2 * j + k)

    # G_3 terms.
    total = d_f**3 / 20160 * s_q0_u
    total += -d_f**2 / 360 * r_q1_u
    total += -d_f**2 / 288 * moment_Pj_Qk(2, 0)
    total += d_f / 12 * moment_P_Qk(2)
    total += -Fraction(1, 4) * moment_QK(group, 4)

    # G_1 G_2.
    total += -d_f**3 / 4320 * Fraction(5, d + 2) * r_q2_u
    total += d_f**2 / 144 * Fraction(5, d + 2) ** 2 * moment_QK(group, 5)
    total += -d_f / 36 * Fraction(5, d + 2) * moment_QK(group, 5)
    total += d_f**2 / 720 * r_q2_u
    total += -d_f / 24 * Fraction(5, d + 2) * moment_QK(group, 5)
    total += Fraction(1, 6) * moment_QK(group, 5)

    # G_1^3/6.
    total += (
        d_f**3 * Fraction(5, d + 2) ** 3 / 1728 * moment_QK(group, 6)
        - 3 * d_f**2 * Fraction(5, d + 2) ** 2 / 288 * moment_QK(group, 6)
        + 3 * d_f * Fraction(5, d + 2) / 48 * moment_QK(group, 6)
        - Fraction(1, 8) * moment_QK(group, 6)
    ) / 6

    # G_2 H_1.
    total += d_f**3 / (360 * 12) * r_q1_u
    total += -d_f**2 / 144 * moment_P_Qk(2)
    total += d_f / 36 * moment_QK(group, 4)

    # G_1^2 H_1 / 2.
    total += (
        -d_f**3 / (144 * 12) * Fraction(5, d + 2) ** 2 * moment_QK(group, 5)
        + d_f**2 / 144 * Fraction(5, d + 2) * moment_QK(group, 5)
        - d_f / 48 * moment_QK(group, 5)
    ) / 2

    # G_1 H_2.
    total += d_f**3 / (12 * 288) * moment_P_Qk(2)
    total += -d_f**3 / (12 * 1440) * Fraction(5, d + 2) ** 2 * moment_QK(group, 4)
    total += -d_f**2 / (2 * 288) * moment_QK(group, 4)
    total += d_f**2 / (2 * 1440) * moment_P_Qk(2)

    # H_3.
    total += -d_f**3 / 10368 * moment_QK(group, 3)
    total += d_f**3 / 17280 * moment_P_Qk(1)
    total += -d_f**3 / 60480 * r_q0_u

    return total


def check_f4_normalization() -> None:
    roots = f4_positive_roots()
    s_poly = root_power_sum(roots, 8)
    e_s = dunkl_expectation(s_poly, roots, 8, "F4 S")
    expected = Fraction(6_265_350)
    print(f"F4 <S(w)>_Dunkl = {e_s}", flush=True)
    print(f"F4 <S(w)>_Wick   = {expected}", flush=True)
    if e_s != expected:
        raise AssertionError(f"F4 normalization mismatch: {e_s} != {expected}")

    radial_poly = quadratic_power(identity_matrix(4), 4)
    radial = dunkl_expectation(radial_poly, roots, 8, "F4 |w|^8")
    radial_expected = gamma_moment_from_group("F4", 4)
    print(f"F4 <|w|^8>_Dunkl = {radial}", flush=True)
    print(f"F4 <|w|^8>_Gamma  = {radial_expected}", flush=True)
    if radial != radial_expected:
        raise AssertionError(f"F4 radial mismatch: {radial} != {radial_expected}")

    roots, projector = group_roots_and_projector("F4")
    r_w = root_sum_with_q_power(roots, projector, 6, 0, "F4 R")
    c2 = c2_from_R_moment("F4", r_w)
    c2_expected = Fraction(8_401_705, 486)
    print(f"F4 c_2 from Dunkl = {c2}", flush=True)
    print(f"F4 c_2 expected   = {c2_expected}", flush=True)
    if c2 != c2_expected:
        raise AssertionError(f"F4 c2 mismatch: {c2} != {c2_expected}")

    cp = Fraction(ROOT_DATA["F4"]["A"], ROOT_DATA["F4"]["r"])
    r_q0_u = r_w / cp**3
    r_q1_u = root_sum_with_q_power(roots, projector, 6, 1, "F4 RQ") / cp**3
    r_q2_u = root_sum_with_q_power(roots, projector, 6, 2, "F4 RQ2") / cp**3
    s_q0_u = e_s / cp**4
    c3 = c3_from_RS_moments("F4", r_q0_u, r_q1_u, r_q2_u, s_q0_u)
    c3_expected = Fraction(-44_918_139_383, 39_366)
    print(f"F4 c_3 from Dunkl = {c3}", flush=True)
    print(f"F4 c_3 expected   = {c3_expected}", flush=True)
    if c3 != c3_expected:
        raise AssertionError(f"F4 c3 mismatch: {c3} != {c3_expected}")


def compute_e8_c3() -> None:
    roots, projector = group_roots_and_projector("E8")

    radial_poly = quadratic_power(projector, 4)
    radial = dunkl_expectation(radial_poly, roots, 8, "E8 |w|^8")
    radial_expected = gamma_moment_from_group("E8", 4)
    print(f"E8 <|w|^8>_Dunkl = {radial}", flush=True)
    print(f"E8 <|w|^8>_Gamma  = {radial_expected}", flush=True)
    if radial != radial_expected:
        raise AssertionError(f"E8 radial mismatch: {radial} != {radial_expected}")

    s_poly = root_power_sum(roots, 8)
    e_s_w = dunkl_expectation(s_poly, roots, 8, "E8 S")

    d = ROOT_DATA["E8"]["d"]
    c_prime = Fraction(ROOT_DATA["E8"]["A"], ROOT_DATA["E8"]["r"])
    radial_s_u = eta_1("E8") * moment_QK("E8", 4) / c_prime**4
    exact_s_u = e_s_w / c_prime**4
    s_primary_u = exact_s_u - radial_s_u
    s_correction = Fraction(d**3, 20160) * s_primary_u
    c3_radial = compute_c3("E8", assume_zeros=True)
    c3_exact = c3_radial + s_correction

    print("\n=== E8 degree-8 Dunkl result ===", flush=True)
    print(f"<S(w)> = {e_s_w}", flush=True)
    print(f"<S(u)> = {exact_s_u}", flush=True)
    print(f"radial <S(u)> = {radial_s_u}", flush=True)
    print(f"<S_primary(u)> = {s_primary_u}", flush=True)
    print(f"S correction to c_3(E8) = {s_correction}", flush=True)
    print(f"c_3 radial-only(E8) = {c3_radial}", flush=True)
    print(f"c_3 exact(E8) = {c3_exact}", flush=True)
    print(f"c_3 exact(E8) approx = {float(c3_exact):.12g}", flush=True)
    threshold = (8.0 * abs(float(c3_exact))) ** (1.0 / 3.0)
    print(f"reduced threshold (8|c_3|)^(1/3) approx = {threshold:.6f}", flush=True)


def compute_group_coefficients(group: str) -> None:
    roots, projector = group_roots_and_projector(group)
    cp = Fraction(ROOT_DATA[group]["A"], ROOT_DATA[group]["r"])

    print(f"\n=== {group} Dunkl coefficient computation ===", flush=True)
    radial = dunkl_expectation(quadratic_power(projector, 3), roots, 6, f"{group} |w|^6")
    radial_expected = gamma_moment_from_group(group, 3)
    print(f"{group} <|w|^6>_Dunkl = {radial}", flush=True)
    print(f"{group} <|w|^6>_Gamma  = {radial_expected}", flush=True)
    if radial != radial_expected:
        raise AssertionError(f"{group} radial degree-6 mismatch")

    r_w = root_sum_with_q_power(roots, projector, 6, 0, f"{group} R")
    r_q1_w = root_sum_with_q_power(roots, projector, 6, 1, f"{group} RQ")
    r_q2_w = root_sum_with_q_power(roots, projector, 6, 2, f"{group} RQ2")
    s_w = root_sum_with_q_power(roots, projector, 8, 0, f"{group} S")

    c2 = c2_from_R_moment(group, r_w)
    r_q0_u = r_w / cp**3
    r_q1_u = r_q1_w / cp**3
    r_q2_u = r_q2_w / cp**3
    s_q0_u = s_w / cp**4
    c3 = c3_from_RS_moments(group, r_q0_u, r_q1_u, r_q2_u, s_q0_u)
    threshold = (8.0 * abs(float(c3))) ** (1.0 / 3.0)

    print(f"{group} <R(w)> = {r_w}", flush=True)
    print(f"{group} <R(w)Q> = {r_q1_w}", flush=True)
    print(f"{group} <R(w)Q^2> = {r_q2_w}", flush=True)
    print(f"{group} <S(w)> = {s_w}", flush=True)
    print(f"{group} c_2 exact = {c2}", flush=True)
    print(f"{group} c_2 approx = {float(c2):.12g}", flush=True)
    print(f"{group} c_3 exact = {c3}", flush=True)
    print(f"{group} c_3 approx = {float(c3):.12g}", flush=True)
    print(f"{group} reduced threshold (8|c_3|)^(1/3) approx = {threshold:.6f}", flush=True)


def main() -> None:
    check_f4_normalization()
    compute_group_coefficients("E6")
    compute_group_coefficients("E7")
    compute_e8_c3()


if __name__ == "__main__":
    main()
