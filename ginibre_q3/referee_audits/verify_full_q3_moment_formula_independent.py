#!/usr/bin/env python3
"""Independent exact audit of the full-Q3 moment and EGF formulas.

All finite hierarchy certificates consume the double moment sum in
``lem:full-moment-form``.  This author self-audit derives its coefficients by
formal polynomial multiplication, compares it against direct expectations for
several rational finite laws, reconstructs the derivative generating-function
coefficients independently, and checks the index/sign skeleton in every
load-bearing C++ implementation.  No floating-point arithmetic is used.
"""

from __future__ import annotations

import re
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "full_q3_extension.tex"
RESULT_PATTERN = re.compile(
    r"\\begin\{(theorem|proposition|lemma|corollary)\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
CPP_SOURCES = (
    "verify_full_q3_su2_gmp.cpp",
    "verify_full_q3_su3_gmp.cpp",
    "verify_full_q3_su45_gmp.cpp",
    "verify_full_q3_sun_ge6_gmp.cpp",
    "verify_full_q3_exceptional_gmp.cpp",
    "verify_full_q3_bd_residual_gmp.cpp",
    "verify_full_q3_bcd_remaining_gmp.cpp",
    "verify_full_q3_bcd_bounded_littlewood_gmp.cpp",
)
Polynomial = dict[tuple[int, int], int]


class FormulaFailure(RuntimeError):
    """A fail-closed moment-formula check failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FormulaFailure(message)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def require_fragment(text: str, fragment: str, message: str) -> None:
    require(compact(fragment) in compact(text), message)


def result_with_label(text: str, label: str) -> str:
    needle = rf"\label{{{label}}}"
    matches = [match for match in RESULT_PATTERN.finditer(text) if needle in match.group(2)]
    require(len(matches) == 1, f"expected one result labelled {label}, found {len(matches)}")
    result = matches[0]
    proof = re.match(r"\s*\\begin\{proof\}.*?\\end\{proof\}", text[result.end() :], re.DOTALL)
    require(proof is not None, f"result labelled {label} lacks an immediate proof")
    return result.group(0) + proof.group(0)


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for (left_x, left_y), left_value in left.items():
        for (right_x, right_y), right_value in right.items():
            key = (left_x + right_x, left_y + right_y)
            output[key] = output.get(key, 0) + left_value * right_value
    return {key: value for key, value in output.items() if value}


def linear_power(x_coefficient: int, y_coefficient: int, exponent: int) -> Polynomial:
    return {
        (exponent - y_degree, y_degree):
        comb(exponent, y_degree)
        * x_coefficient ** (exponent - y_degree)
        * y_coefficient**y_degree
        for y_degree in range(exponent + 1)
    }


def direct_polynomial_coefficients(a: int, n: int) -> Polynomial:
    return polynomial_multiply(linear_power(1, -1, 2 * a), linear_power(1, 1, n))


def double_sum_coefficients(a: int, n: int) -> Polynomial:
    output: Polynomial = {}
    for j in range(2 * a + 1):
        for k in range(n + 1):
            key = (2 * a - j + k, j + n - k)
            value = (-1) ** j * comb(2 * a, j) * comb(n, k)
            output[key] = output.get(key, 0) + value
    return {key: value for key, value in output.items() if value}


def moments(distribution: tuple[tuple[Fraction, Fraction], ...], maximum: int) -> list[Fraction]:
    require(sum((weight for _value, weight in distribution), Fraction()) == 1, "weights do not sum to one")
    require(all(weight >= 0 for _value, weight in distribution), "negative probability weight")
    return [
        sum((weight * value**degree for value, weight in distribution), Fraction())
        for degree in range(maximum + 1)
    ]


def moment_formula(moment: list[Fraction], a: int, n: int) -> Fraction:
    require(len(moment) > 2 * a + n, "moment list is too short")
    return sum(
        Fraction((-1) ** j * comb(2 * a, j) * comb(n, k))
        * moment[2 * a - j + k]
        * moment[j + n - k]
        for j in range(2 * a + 1)
        for k in range(n + 1)
    )


def direct_expectation(
    distribution: tuple[tuple[Fraction, Fraction], ...],
    a: int,
    n: int,
) -> Fraction:
    return sum(
        weight_x * weight_y * (x - y) ** (2 * a) * (x + y) ** n
        for x, weight_x in distribution
        for y, weight_y in distribution
    )


def derivative_series(moment: list[Fraction], derivative: int, maximum: int) -> list[Fraction]:
    require(len(moment) > derivative + maximum, "moment list is too short for derivative series")
    return [moment[derivative + degree] / factorial(degree) for degree in range(maximum + 1)]


def series_multiply(left: list[Fraction], right: list[Fraction], maximum: int) -> list[Fraction]:
    return [
        sum((left[index] * right[degree - index] for index in range(degree + 1)), Fraction())
        for degree in range(maximum + 1)
    ]


def egf_formula(moment: list[Fraction], a: int, maximum: int) -> list[Fraction]:
    output = [Fraction() for _ in range(maximum + 1)]
    for j in range(2 * a + 1):
        left = derivative_series(moment, 2 * a - j, maximum)
        right = derivative_series(moment, j, maximum)
        product = series_multiply(left, right, maximum)
        coefficient = Fraction((-1) ** j * comb(2 * a, j))
        output = [value + coefficient * term for value, term in zip(output, product, strict=True)]
    return [coefficient * factorial(degree) for degree, coefficient in enumerate(output)]


def audit_formal_coefficients() -> int:
    checked = 0
    for a in range(9):
        for n in range(17):
            require(
                direct_polynomial_coefficients(a, n) == double_sum_coefficients(a, n),
                f"formal coefficient mismatch at a={a}, n={n}",
            )
            checked += 1
    return checked


def audit_distributions() -> tuple[int, int]:
    distributions = (
        ((Fraction(0), Fraction(1)),),
        ((Fraction(-1), Fraction(1, 3)), (Fraction(2), Fraction(2, 3))),
        (
            (Fraction(-2), Fraction(1, 5)),
            (Fraction(0), Fraction(1, 2)),
            (Fraction(3), Fraction(3, 10)),
        ),
        (
            (Fraction(-3, 2), Fraction(1, 7)),
            (Fraction(-1, 3), Fraction(2, 7)),
            (Fraction(4, 5), Fraction(1, 7)),
            (Fraction(5, 2), Fraction(3, 7)),
        ),
    )
    direct_cases = 0
    egf_cases = 0
    for distribution in distributions:
        moment = moments(distribution, 30)
        for a in range(7):
            egf_values = egf_formula(moment, a, 12)
            for n in range(13):
                expected = direct_expectation(distribution, a, n)
                value = moment_formula(moment, a, n)
                require(value == expected, f"distribution formula mismatch at a={a}, n={n}")
                require(egf_values[n] == expected, f"distribution EGF mismatch at a={a}, n={n}")
                if n % 2 == 0:
                    require(expected >= 0, f"even-n pointwise sector became negative at a={a}, n={n}")
                direct_cases += 1
                egf_cases += 1
    return direct_cases, egf_cases


def audit_formal_egf() -> int:
    moment_ledgers = (
        [Fraction(2**degree + 3**degree) for degree in range(31)],
        [Fraction((-2) ** degree + 5**degree, 7) for degree in range(31)],
        [Fraction(degree**4 - 3 * degree + 11) for degree in range(31)],
    )
    cases = 0
    for moment in moment_ledgers:
        for a in range(6):
            egf_values = egf_formula(moment, a, 10)
            for n in range(11):
                require(egf_values[n] == moment_formula(moment, a, n), f"formal EGF mismatch at a={a}, n={n}")
                cases += 1
    return cases


def function_body(source: str, names: tuple[str, ...]) -> str:
    name_pattern = "|".join(re.escape(name) for name in names)
    matches = list(re.finditer(rf"\b(?:BigInt|mpz_class)\s+(?:{name_pattern})\s*\(", source))
    require(len(matches) == 1, f"expected one hierarchy function among {names}, found {len(matches)}")
    start = matches[0].start()
    brace = source.find("{", start)
    require(brace >= 0, "hierarchy function has no body")
    depth = 0
    for position in range(brace, len(source)):
        if source[position] == "{":
            depth += 1
        elif source[position] == "}":
            depth -= 1
            if depth == 0:
                return source[start : position + 1]
    raise FormulaFailure("unterminated hierarchy function")


def audit_cpp_implementations() -> int:
    directory = ROOT / "character_ring_iter"
    checked = 0
    for name in CPP_SOURCES:
        source = (directory / name).read_text(encoding="utf-8")
        body = function_body(source, ("full_q3_value", "full_q3_hierarchy_value", "hierarchy_value"))
        normalized = compact(body)
        summation_variable = "ell" if name == "verify_full_q3_exceptional_gmp.cpp" else "k"
        patterns = (
            rf"2(?:U)?\*a-j\+{summation_variable}",
            rf"j\+n-{summation_variable}",
            rf"binom(?:ial|_int)\(2(?:U)?\*a,j\)",
            rf"binom(?:ial|_int)\(n,{summation_variable}\)",
        )
        for pattern in patterns:
            require(re.search(pattern, normalized) is not None, f"C++ formula fragment missing in {name}: {pattern}")
        require("answer+=term" in normalized, f"positive-j parity branch missing in {name}")
        require("answer-=term" in normalized, f"negative-j parity branch missing in {name}")
        require(re.search(r"if\([^)]*j", normalized) is not None, f"sign branch is not controlled by j in {name}")
        checked += 1
    return checked


def validate_manuscript(text: str) -> None:
    lemma = result_with_label(text, "lem:full-moment-form")
    fragments = (
        r"\sum_{j=0}^{2a}\sum_{k=0}^{n}",
        r"(-1)^j\binom{2a}{j}\binom nk",
        r"m^G_{2a-j+k}m^G_{j+n-k}",
        r"M_G(t)=\mathbf E e^{tX}",
        r"M_G^{(2a-j)}(t)M_G^{(j)}(t)",
        r"\sum_{n\geq0}I^G_{a,n}\frac{t^n}{n!}",
        r"identity of convergent power series",
        r"The adjoint character is bounded on the compact group",
        r"all series converge absolutely and termwise integration is justified",
    )
    for fragment in fragments:
        require_fragment(lemma, fragment, f"moment-form lemma fragment missing: {fragment}")


def expect_rejected(name: str, action: Callable[[], object]) -> None:
    try:
        action()
    except FormulaFailure:
        return
    raise FormulaFailure(f"adversarial moment-formula mutation was accepted: {name}")


def mutation_self_tests(text: str) -> int:
    lemma = result_with_label(text, "lem:full-moment-form")
    mutations = (
        (r"(-1)^j", r"(-1)^k", "summation sign"),
        (r"m^G_{2a-j+k}", r"m^G_{2a+j+k}", "first moment index"),
        (r"M_G^{(2a-j)}(t)", r"M_G^{(2a+j)}(t)", "derivative index"),
        (r"\frac{t^n}{n!}", r"\frac{t^n}{(n+1)!}", "EGF normalization"),
        ("bounded", "measurable", "termwise convergence"),
    )
    for old, new, name in mutations:
        require(lemma.count(old) >= 1, f"mutation target changed: {name}")
        mutated_lemma = lemma.replace(old, new)
        mutated_text = text.replace(lemma, mutated_lemma, 1)
        expect_rejected(name, lambda value=mutated_text: validate_manuscript(value))
    return len(mutations)


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")
    formal_coefficients = audit_formal_coefficients()
    distribution_cases, distribution_egf_cases = audit_distributions()
    formal_egf_cases = audit_formal_egf()
    cpp_implementations = audit_cpp_implementations()
    validate_manuscript(text)
    mutations = mutation_self_tests(text)
    print(
        "FULL_Q3_MOMENT_FORMULA "
        f"formal_coefficient_cases={formal_coefficients} "
        f"rational_distribution_cases={distribution_cases} "
        f"distribution_egf_cases={distribution_egf_cases} formal_egf_cases={formal_egf_cases} "
        f"cpp_implementations={cpp_implementations} exact_arithmetic=PASS "
        f"sign_indices=PASS convergence_scope=PASS mutations_rejected={mutations}"
    )
    print("FULL_Q3_MOMENT_FORMULA VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FormulaFailure, OSError, UnicodeError) as error:
        raise SystemExit(f"FULL_Q3_MOMENT_FORMULA FAILURE: {error}") from error
