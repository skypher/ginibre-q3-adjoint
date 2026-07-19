#!/usr/bin/env python3
"""Independent algebraic audit of the adjoint-generator cone promotion.

For ``u=X+Y`` and ``v=X-Y``, the signed monomial factors have the exact
nonnegative expansions

    X^d + Y^d = 2^(1-d) sum_{k even} C(d,k) u^(d-k) v^k,
    X^d - Y^d = 2^(1-d) sum_{k odd } C(d,k) u^(d-k) v^k.

Thus every tuple of nonnegative polynomials expands with nonnegative
coefficients, and the parity of the total ``v`` exponent is the parity of the
number of minus signs.  This script verifies the identities over Q, audits the
paper's multiplicative factorization and full quantifiers, and checks the
uniform-closure passage.  It is a structural audit, not a replacement for the
hierarchy inequalities proved by Part III.
"""

from __future__ import annotations

import itertools
import re
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "full_q3_extension.tex"
RESULT_PATTERN = re.compile(
    r"\\begin\{(theorem|proposition|lemma|corollary)\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
DEFINITION_PATTERN = re.compile(
    r"\\begin\{definition\}(.*?)\\end\{definition\}",
    re.DOTALL,
)
Polynomial = dict[tuple[int, ...], Fraction]


class PromotionFailure(RuntimeError):
    """A fail-closed cone-promotion invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PromotionFailure(message)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def require_fragment(text: str, fragment: str, message: str) -> None:
    require(compact(fragment) in compact(text), message)


def block_with_title(text: str, pattern: re.Pattern[str], title: str) -> str:
    matches = [match.group(0) for match in pattern.finditer(text) if title in match.group(1)]
    require(len(matches) == 1, f"expected one block titled {title}, found {len(matches)}")
    return matches[0]


def result_with_label(text: str, label: str, include_proof: bool = False) -> str:
    needle = rf"\label{{{label}}}"
    matches = [match for match in RESULT_PATTERN.finditer(text) if needle in match.group(2)]
    require(len(matches) == 1, f"expected one result labelled {label}, found {len(matches)}")
    result = matches[0]
    if not include_proof:
        return result.group(0)
    proof = re.match(r"\s*\\begin\{proof\}.*?\\end\{proof\}", text[result.end() :], re.DOTALL)
    require(proof is not None, f"result labelled {label} lacks an immediate proof")
    return result.group(0) + proof.group(0)


def clean(polynomial: Polynomial) -> Polynomial:
    return {monomial: coefficient for monomial, coefficient in polynomial.items() if coefficient}


def add(*polynomials: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            output[monomial] = output.get(monomial, Fraction()) + coefficient
    return clean(output)


def scale(polynomial: Polynomial, coefficient: Fraction) -> Polynomial:
    return clean({monomial: coefficient * value for monomial, value in polynomial.items()})


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_degree + right_degree
                for left_degree, right_degree in zip(left_monomial, right_monomial, strict=True)
            )
            output[monomial] = output.get(monomial, Fraction()) + left_coefficient * right_coefficient
    return clean(output)


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    require(exponent >= 0, f"negative polynomial exponent {exponent}")
    variable_count = len(next(iter(polynomial))) if polynomial else 2
    output: Polynomial = {(0,) * variable_count: Fraction(1)}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            output = multiply(output, base)
        remaining >>= 1
        if remaining:
            base = multiply(base, base)
    return output


def variable(index: int, count: int) -> Polynomial:
    exponents = [0] * count
    exponents[index] = 1
    return {tuple(exponents): Fraction(1)}


def direct_signed_monomial(degree: int, sign: int) -> Polynomial:
    require(sign in {-1, 1}, f"invalid sign {sign}")
    x = variable(0, 2)
    y = variable(1, 2)
    return add(power(x, degree), scale(power(y, degree), Fraction(sign)))


def uv_signed_monomial(degree: int, sign: int) -> Polynomial:
    """Return the signed monomial in the formal variables (u,v)."""

    require(degree >= 0, f"negative monomial degree {degree}")
    require(sign in {-1, 1}, f"invalid sign {sign}")
    parity = 0 if sign == 1 else 1
    coefficient_scale = Fraction(2, 2**degree)
    output: Polynomial = {}
    for v_degree in range(parity, degree + 1, 2):
        output[(degree - v_degree, v_degree)] = coefficient_scale * comb(degree, v_degree)
    return clean(output)


def substitute_uv(polynomial: Polynomial) -> Polynomial:
    u = add(variable(0, 2), variable(1, 2))
    v = add(variable(0, 2), scale(variable(1, 2), Fraction(-1)))
    output: Polynomial = {}
    for (u_degree, v_degree), coefficient in polynomial.items():
        term = multiply(power(u, u_degree), power(v, v_degree))
        output = add(output, scale(term, coefficient))
    return output


def verify_monomial_expansions(maximum_degree: int = 32) -> tuple[int, int]:
    cases = 0
    coefficients = 0
    for degree in range(maximum_degree + 1):
        for sign in (1, -1):
            expansion = uv_signed_monomial(degree, sign)
            require(
                substitute_uv(expansion) == direct_signed_monomial(degree, sign),
                f"u/v expansion mismatch at degree={degree}, sign={sign}",
            )
            expected_parity = 0 if sign == 1 else 1
            for (_u_degree, v_degree), coefficient in expansion.items():
                require(coefficient > 0, f"nonpositive u/v coefficient at degree={degree}")
                require(v_degree % 2 == expected_parity, f"wrong v parity at degree={degree}")
                coefficients += 1
            if degree == 0 and sign == 1:
                require(expansion == {(0, 0): Fraction(2)}, "constant plus factor is not 2")
            if degree == 0 and sign == -1:
                require(not expansion, "constant minus factor is not zero")
            cases += 1
    return cases, coefficients


def verify_tuple_parity() -> tuple[int, int]:
    tuple_cases = 0
    terms = 0
    # Exhaust every tuple through length four with degrees 0..4.  The general
    # statement follows algebraically because multiplication adds v exponents.
    for length in range(5):
        for degrees in itertools.product(range(5), repeat=length):
            for signs in itertools.product((1, -1), repeat=length):
                expansion: Polynomial = {(0, 0): Fraction(1)}
                for degree, sign in zip(degrees, signs, strict=True):
                    expansion = multiply(expansion, uv_signed_monomial(degree, sign))
                expected_parity = sum(sign == -1 for sign in signs) % 2
                for (_u_degree, v_degree), coefficient in expansion.items():
                    require(coefficient > 0, f"nonpositive tuple coefficient: {degrees}/{signs}")
                    require(v_degree % 2 == expected_parity, f"tuple parity mismatch: {degrees}/{signs}")
                    terms += 1
                tuple_cases += 1
    return tuple_cases, terms


def verify_factorization_identity() -> int:
    fx, fy, gx, gy = (variable(index, 4) for index in range(4))
    for sign in (1, -1):
        left = add(multiply(fx, gx), scale(multiply(fy, gy), Fraction(sign)))
        right = scale(
            add(
                multiply(add(fx, fy), add(gx, scale(gy, Fraction(sign)))),
                multiply(
                    add(fx, scale(fy, Fraction(-1))),
                    add(gx, scale(gy, Fraction(-sign))),
                ),
            ),
            Fraction(1, 2),
        )
        require(left == right, f"multiplicative factorization fails for sign={sign}")
    return 2


def validate_manuscript(text: str) -> None:
    definition = block_with_title(text, DEFINITION_PATTERN, "[Adjoint-generated cone]")
    reduction = result_with_label(text, "thm:full-cone-reduction", include_proof=True)
    final = result_with_label(text, "thm:full-adjoint-generated-q3")
    fragments = (
        (definition, r"\sum_{j=0}^{d}c_j\chiad^j", "cone polynomial form"),
        (definition, r"d\geq0,\ c_j\geq0", "nonnegative cone coefficients"),
        (definition, r"^{\lVert\cdot\rVert_\infty}", "uniform cone closure"),
        (reduction, r"S_G=\{\chiad\}", "singleton generator"),
        (reduction, r"It vanishes when $r$ is odd", "odd generator-minus vanishing"),
        (reduction, r"for $r=2a$ it is $I^G_{a,n}$", "even-minus hierarchy"),
        (reduction, r"(fg)(x)\pm(fg)(y)", "multiplicative factorization left side"),
        (reduction, r"=\frac12\bigl[", "multiplicative factorization coefficient"),
        (reduction, r"(g(x)\mp g(y))", "multiplicative sign switch"),
        (reduction, r"with nonnegative coefficients", "nonnegative expansion"),
        (reduction, r"uniform limits preserve it because Haar measure is finite", "uniform closure passage"),
        (reduction, r"Proposition~1 and (1.5)", "Ginibre promotion citation"),
        (final, r"every integer $r\geq0$", "all tuple lengths"),
        (final, r"every tuple $f_1,\ldots,f_r\in\Qcone_G$", "all cone tuples"),
        (final, r"every choice of signs", "all sign choices"),
        (final, r"containing an even number of minus signs", "even-minus conclusion"),
        (final, r"If the number of minus signs is odd, the integral is zero", "odd-minus conclusion"),
    )
    for block, fragment, name in fragments:
        require_fragment(block, fragment, f"manuscript promotion fragment missing: {name}")


def replace_in_block(text: str, block: str, old: str, new: str, all_occurrences: bool = False) -> str:
    require(text.count(block) == 1, "mutation block is not unique")
    require(block.count(old) >= 1, f"mutation source absent: {old}")
    changed = block.replace(old, new) if all_occurrences else block.replace(old, new, 1)
    return text.replace(block, changed, 1)


def expect_rejected(name: str, action: Callable[[], object]) -> None:
    try:
        action()
    except PromotionFailure:
        return
    raise PromotionFailure(f"adversarial promotion mutation was accepted: {name}")


def mutation_self_tests(text: str) -> int:
    definition = block_with_title(text, DEFINITION_PATTERN, "[Adjoint-generated cone]")
    reduction = result_with_label(text, "thm:full-cone-reduction", include_proof=True)
    final = result_with_label(text, "thm:full-adjoint-generated-q3")

    mutated = replace_in_block(text, definition, r"c_j\geq0", r"c_j\in\mathbb R")
    expect_rejected("negative cone coefficients allowed", lambda: validate_manuscript(mutated))

    mutated = replace_in_block(text, definition, r"\lVert\cdot\rVert_\infty", r"\lVert\cdot\rVert_2")
    expect_rejected("nonuniform closure", lambda: validate_manuscript(mutated))

    mutated = replace_in_block(text, reduction, r"\frac12", r"-\frac12", all_occurrences=True)
    expect_rejected("factorization coefficient sign", lambda: validate_manuscript(mutated))

    mutated = replace_in_block(text, reduction, r"(g(x)\mp g(y))", r"(g(x)\pm g(y))")
    expect_rejected("factorization sign switch", lambda: validate_manuscript(mutated))

    mutated = replace_in_block(
        text,
        final,
        "containing an even number of minus",
        "containing an odd number of minus",
    )
    expect_rejected("final minus parity", lambda: validate_manuscript(mutated))
    return 5


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")
    monomial_cases, monomial_coefficients = verify_monomial_expansions()
    tuple_cases, tuple_terms = verify_tuple_parity()
    factorization_cases = verify_factorization_identity()
    validate_manuscript(text)
    mutations = mutation_self_tests(text)
    print(
        "FULL_Q3_CONE_PROMOTION "
        f"monomial_cases={monomial_cases} monomial_coefficients={monomial_coefficients} "
        f"tuple_cases={tuple_cases} tuple_terms={tuple_terms} "
        f"factorization_cases={factorization_cases} constants=PASS nonnegative=PASS "
        "minus_parity=PASS uniform_closure=PASS all_tuple_quantifiers=PASS "
        f"mutations_rejected={mutations}"
    )
    print("FULL_Q3_CONE_PROMOTION VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PromotionFailure, OSError, UnicodeError) as error:
        raise SystemExit(f"FULL_Q3_CONE_PROMOTION FAILURE: {error}") from error
