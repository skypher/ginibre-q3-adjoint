#!/usr/bin/env python3
"""Independent exact audit of the SO(3) all-positive-definite obstruction.

The Part III scope boundary rests on a concrete ``--++`` counterexample.  This
referee-side verifier reconstructs the invariant-tensor projection for the
standard SO(3) module, derives its second and fourth Haar moments, and expands
all 16 signed-product terms with exact rational arithmetic.  It also checks
that the manuscript supplies the matrix-coefficient positive-definiteness
factorization and the same vectors/sign pattern.  No floating-point arithmetic
is used.
"""

from __future__ import annotations

import itertools
import re
from fractions import Fraction
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "full_q3_extension.tex"
RESULT_PATTERN = re.compile(
    r"\\begin\{(theorem|proposition|lemma|corollary)\}(.*?)\\end\{\1\}",
    re.DOTALL,
)


class ObstructionFailure(RuntimeError):
    """A fail-closed SO(3) obstruction check failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ObstructionFailure(message)


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


def spin_one_tensor_multiplicities(power: int) -> dict[int, int]:
    multiplicities = {0: 1}
    for _ in range(power):
        following: dict[int, int] = {}
        for spin, multiplicity in multiplicities.items():
            for target in range(abs(spin - 1), spin + 2):
                following[target] = following.get(target, 0) + multiplicity
        multiplicities = following
    return multiplicities


def pairing_value(kind: int, indices: tuple[int, int, int, int]) -> int:
    i, j, k, ell = indices
    if kind == 0:
        return int(i == j and k == ell)
    if kind == 1:
        return int(i == k and j == ell)
    require(kind == 2, f"unknown pairing kind {kind}")
    return int(i == ell and j == k)


def pairing_gram(dimension: int) -> list[list[Fraction]]:
    output: list[list[Fraction]] = []
    index_tuples = itertools.product(range(dimension), repeat=4)
    cached_indices = list(index_tuples)
    for left in range(3):
        row: list[Fraction] = []
        for right in range(3):
            value = sum(
                pairing_value(left, indices) * pairing_value(right, indices)
                for indices in cached_indices
            )
            row.append(Fraction(value))
        output.append(row)
    return output


def invert(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    require(size > 0 and all(len(row) == size for row in matrix), "matrix is not square")
    augmented = [
        row[:] + [Fraction(int(column == index)) for column in range(size)]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column]), None)
        require(pivot is not None, f"singular matrix at column {column}")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [entry / pivot_value for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(augmented[row], augmented[column], strict=True)
            ]
    return [row[size:] for row in augmented]


def determinant_3(matrix: list[list[Fraction]]) -> Fraction:
    require(len(matrix) == 3 and all(len(row) == 3 for row in matrix), "expected a 3x3 matrix")
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def quadratic_form(vector: tuple[Fraction, ...], matrix: list[list[Fraction]]) -> Fraction:
    return sum(
        vector[left] * matrix[left][right] * vector[right]
        for left in range(len(vector))
        for right in range(len(vector))
    )


def exact_obstruction() -> dict[str, object]:
    # The chosen unit vectors have squared coordinate lengths 19/31 and 12/31.
    x_squared = Fraction(19, 31)
    y_squared = Fraction(12, 31)
    require(x_squared + y_squared == 1, "v1/v2 are not unit vectors")
    dot_12 = x_squared - y_squared
    dot_34 = Fraction(1)
    c = (
        dot_12 * dot_34,
        x_squared,
        x_squared,
    )
    require(c == (Fraction(7, 31), Fraction(19, 31), Fraction(19, 31)), "vector contractions changed")

    multiplicities_2 = spin_one_tensor_multiplicities(2)
    multiplicities_4 = spin_one_tensor_multiplicities(4)
    require(multiplicities_2.get(0) == 1, "standard square has wrong invariant dimension")
    require(multiplicities_4.get(0) == 3, "standard fourth power has wrong invariant dimension")

    gram = pairing_gram(3)
    expected_gram = [
        [Fraction(9 if row == column else 3) for column in range(3)]
        for row in range(3)
    ]
    require(gram == expected_gram, f"pairing-tensor Gram matrix={gram}")
    gram_determinant = determinant_3(gram)
    require(gram_determinant == 540, f"pairing Gram determinant={gram_determinant}")
    inverse = invert(gram)
    expected_inverse = [
        [Fraction(2, 15) if row == column else Fraction(-1, 30) for column in range(3)]
        for row in range(3)
    ]
    require(inverse == expected_inverse, f"pairing Gram inverse={inverse}")
    fourth_moment = quadratic_form(c, inverse)
    require(fourth_moment == Fraction(61, 961), f"fourth moment={fourth_moment}")

    squared_dots = {
        (0, 1): dot_12**2,
        (0, 2): x_squared,
        (0, 3): x_squared,
        (1, 2): x_squared,
        (1, 3): x_squared,
        (2, 3): dot_34**2,
    }
    second_moments = {pair: value / 3 for pair, value in squared_dots.items()}
    signs = (-1, -1, 1, 1)
    signed_value = Fraction()
    nonzero_terms = 0
    indices = frozenset(range(4))
    for mask in range(1 << 4):
        on_g = frozenset(index for index in range(4) if mask & (1 << index))
        on_h = indices - on_g
        coefficient = 1
        for index in on_h:
            coefficient *= signs[index]
        if len(on_g) in {1, 3}:
            contribution = Fraction()
        elif len(on_g) == 0 or len(on_g) == 4:
            contribution = fourth_moment
        else:
            pair_g = tuple(sorted(on_g))
            pair_h = tuple(sorted(on_h))
            contribution = second_moments[pair_g] * second_moments[pair_h]
        signed_value += coefficient * contribution
        nonzero_terms += int(contribution != 0)
    require(nonzero_terms == 8, f"nonzero signed-expansion terms={nonzero_terms}")
    require(signed_value == Fraction(-8, 279), f"signed obstruction={signed_value}")
    return {
        "gram": gram,
        "gram_determinant": gram_determinant,
        "inverse": inverse,
        "c": c,
        "fourth_moment": fourth_moment,
        "signed_value": signed_value,
        "nonzero_terms": nonzero_terms,
        "invariants_2": multiplicities_2[0],
        "invariants_4": multiplicities_4[0],
    }


def validate_manuscript(text: str) -> None:
    counterexample = result_with_label(text, "prop:all-pd-counterexample")
    fragments = (
        r"$G=\mathrm{SO}(3)$",
        r"four functions in $\mathcal P_G$ and the sign pattern $--++$",
        r"-\frac8{279}",
        r"f_v(g)=\langle v,gv\rangle",
        r"\left\lVert\sum_{j=1}^r z_jg_jv\right\rVert^2\geq0",
        r"diagonal entries $9$ and off-diagonal entries $3$",
        r"inverse has diagonal entries $2/15$ and off-diagonal entries $-1/30$",
        r"(c_1,c_2,c_3)=(7,19,19)/31",
        r"equals $61/961$",
        r"M_{12}M_{34}-M_{13}M_{24}-M_{14}M_{23}",
    )
    for fragment in fragments:
        require_fragment(counterexample, fragment, f"counterexample fragment missing: {fragment}")
    require_fragment(
        text,
        r"all real continuous positive-definite functions",
        "all-positive-definite scope boundary is absent",
    )
    require_fragment(
        text,
        r"Such a replacement is impossible",
        "final scope warning is absent",
    )


def expect_rejected(name: str, action: Callable[[], object]) -> None:
    try:
        action()
    except ObstructionFailure:
        return
    raise ObstructionFailure(f"adversarial obstruction mutation was accepted: {name}")


def mutation_self_tests(text: str) -> int:
    counterexample = result_with_label(text, "prop:all-pd-counterexample")
    mutations = (
        (r"$--++$", r"$-+++$", "sign pattern"),
        (r"-\frac8{279}", r"-\frac7{279}", "final fraction"),
        ("off-diagonal entries $3$", "off-diagonal entries $4$", "Gram entry"),
        (r"(c_1,c_2,c_3)=(7,19,19)/31", r"(c_1,c_2,c_3)=(8,19,19)/31", "contractions"),
        ("equals $61/961$", "equals $60/961$", "fourth moment"),
    )
    for old, new, name in mutations:
        require(counterexample.count(old) >= 1, f"mutation target changed: {name}")
        mutated_block = counterexample.replace(old, new)
        mutated_text = text.replace(counterexample, mutated_block, 1)
        expect_rejected(name, lambda value=mutated_text: validate_manuscript(value))
    return len(mutations)


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")
    values = exact_obstruction()
    validate_manuscript(text)
    mutations = mutation_self_tests(text)
    print(
        "FULL_Q3_SO3_OBSTRUCTION "
        f"invariants_V2={values['invariants_2']} invariants_V4={values['invariants_4']} "
        f"pairing_gram_det={values['gram_determinant']} contractions=7,19,19/31 "
        f"fourth_moment={values['fourth_moment']} signed_terms={values['nonzero_terms']}/16 "
        f"signed_value={values['signed_value']} pd_gram_factorization=PASS "
        f"mutations_rejected={mutations}"
    )
    print("FULL_Q3_SO3_OBSTRUCTION VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ObstructionFailure, OSError, UnicodeError) as error:
        raise SystemExit(f"FULL_Q3_SO3_OBSTRUCTION FAILURE: {error}") from error
