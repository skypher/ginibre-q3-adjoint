#!/usr/bin/env python3
"""Exact all-exponent GKS2* certificate for the SU(2)_4 fusion ring.

The script enumerates every support-disjoint sign chamber and every parity
class of its positive exponents.  After writing p_a=p_a^0+2r_a, it applies
an exact capacitated Hall transport in Q(sqrt(3)).  The fourteen regimes in
which coordinatewise transport is too strong are checked against the two
closed inequalities recorded in SU2_LEVELS_1_4_GKS2_2026_07_24.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class Quad:
    """An exact number a+b*sqrt(3), with rational a,b."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(value: int | Fraction | "Quad") -> "Quad":
        return value if isinstance(value, Quad) else Quad(Fraction(value), Fraction(0))

    def __add__(self, other: int | Fraction | "Quad") -> "Quad":
        rhs = Quad.coerce(other)
        return Quad(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> "Quad":
        return Quad(-self.a, -self.b)

    def __sub__(self, other: int | Fraction | "Quad") -> "Quad":
        return self + (-Quad.coerce(other))

    def __rsub__(self, other: int | Fraction | "Quad") -> "Quad":
        return Quad.coerce(other) - self

    def __mul__(self, other: int | Fraction | "Quad") -> "Quad":
        rhs = Quad.coerce(other)
        return Quad(
            self.a * rhs.a + 3 * self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a,
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Quad":
        if exponent < 0:
            raise ValueError("negative exponent")
        result = Quad(Fraction(1))
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power >>= 1
        return result

    def sign(self) -> int:
        """Return -1, 0, or 1 without floating-point arithmetic."""
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        rational_square = self.a * self.a
        radical_square = 3 * self.b * self.b
        if rational_square == radical_square:
            return 0
        if self.a > 0:
            return 1 if rational_square > radical_square else -1
        return -1 if rational_square > radical_square else 1


ZERO = Quad()
ONE = Quad(Fraction(1))
ROOT3 = Quad(Fraction(0), Fraction(1))

VALUES: Tuple[Tuple[Quad, ...], ...] = (
    (ROOT3, Quad(2), ROOT3, ONE),
    (ONE, ZERO, Quad(-1), Quad(-1)),
    (ZERO, Quad(-1), ZERO, ONE),
    (Quad(-1), ZERO, ONE, Quad(-1)),
    (-ROOT3, Quad(2), -ROOT3, ONE),
)
WEIGHTS: Tuple[Quad, ...] = (
    Quad(Fraction(1, 12)),
    Quad(Fraction(1, 4)),
    Quad(Fraction(1, 3)),
    Quad(Fraction(1, 4)),
    Quad(Fraction(1, 12)),
)

Status = Tuple[int, int, int, int]
Parity = Tuple[int, int, int, int]
Bases = Dict[int, Quad]
Term = Tuple[Quad, Bases, Tuple[int, int]]


def status_name(status: Status) -> str:
    return "".join("0" if value == 0 else "+" if value > 0 else "-" for value in status)


def dominates(left: Bases, right: Bases, active: Sequence[int]) -> bool:
    return all((left[index] - right[index]).sign() >= 0 for index in active)


def hall_transport(
    positive: Sequence[Term], negative: Sequence[Term], active: Sequence[int]
) -> bool:
    for subset_mask in range(1, 1 << len(negative)):
        demand = sum(
            (negative[index][0] for index in range(len(negative)) if subset_mask >> index & 1),
            ZERO,
        )
        capacity = ZERO
        for coefficient, bases, _pair in positive:
            if any(
                subset_mask >> index & 1
                and dominates(bases, negative[index][1], active)
                for index in range(len(negative))
            ):
                capacity += coefficient
        if (capacity - demand).sign() < 0:
            return False
    return True


def build_terms(status: Status, parity: Parity) -> Tuple[List[Term], List[Term], List[int]]:
    active = [index for index, value in enumerate(status) if value]
    minimum_power = {index: 1 if parity[index] == 1 else 2 for index in active}
    positive: List[Term] = []
    negative: List[Term] = []
    for left in range(5):
        for right in range(left + 1, 5):
            coefficient = Quad(2) * WEIGHTS[left] * WEIGHTS[right]
            bases: Bases = {}
            vanishes = False
            for label in active:
                base = VALUES[left][label] + status[label] * VALUES[right][label]
                if base == ZERO:
                    vanishes = True
                    break
                coefficient *= base ** minimum_power[label]
                bases[label] = base * base
            if vanishes or coefficient == ZERO:
                continue
            if coefficient.sign() > 0:
                positive.append((coefficient, bases, (left, right)))
            else:
                negative.append((-coefficient, bases, (left, right)))
    return positive, negative, active


RADIAL_EXCEPTIONS = {
    ("00-+", (-1, -1, 0, 1)),
    ("-00+", (0, -1, -1, 1)),
    ("-0-+", (1, -1, 1, 0)),
}
CROSSED_EXCEPTIONS = {
    ("+0+-", (1, -1, 1, 0)),
    ("+++-", (0, 1, 0, 0)),
    ("+++-", (1, 0, 1, 0)),
    ("+-+-", (0, 1, 0, 1)),
    ("+-+-", (1, 0, 1, 0)),
    ("-0--", (1, -1, 1, 0)),
    ("-+--", (0, 1, 0, 0)),
    ("-+--", (1, 0, 1, 0)),
    ("----", (0, 1, 0, 1)),
    ("----", (1, 0, 1, 0)),
}
COMBINED_EXCEPTIONS = {("-0-0", (1, -1, 1, -1))}


def verify_radial_template(maximum_residual: int = 40) -> None:
    for exponent in range(maximum_residual + 1):
        if 2 * 3**exponent + 12**exponent < 3 * 4**exponent:
            raise AssertionError("radial template failure")


def verify_crossed_template(maximum_left: int = 30, maximum_right: int = 30) -> None:
    large = Quad(4) + Quad(2) * ROOT3
    small = Quad(4) - Quad(2) * ROOT3
    if large * small != Quad(4):
        raise AssertionError("crossed product identity failure")
    for left in range(maximum_left + 1):
        for right in range(maximum_right + 1):
            value = large**left * small**right + small**left * large**right - Quad(2)
            if value.sign() < 0:
                raise AssertionError("crossed template failure")


def all_regimes() -> Iterable[Tuple[Status, Parity]]:
    for status_raw in product((0, 1, -1), repeat=4):
        status: Status = status_raw
        if -1 not in status:
            continue
        active = [index for index, value in enumerate(status) if value]
        for active_parity in product((0, 1), repeat=len(active)):
            if sum(
                active_parity[position]
                for position, index in enumerate(active)
                if status[index] == -1
            ) % 2:
                continue
            full = [-1, -1, -1, -1]
            for position, index in enumerate(active):
                full[index] = active_parity[position]
            yield status, tuple(full)


def main() -> None:
    verify_radial_template()
    verify_crossed_template()

    total = direct = radial = crossed = combined = 0
    unresolved: List[Tuple[str, Parity]] = []
    found_exceptions = set()

    for status, parity in all_regimes():
        total += 1
        positive, negative, active = build_terms(status, parity)
        if hall_transport(positive, negative, active):
            direct += 1
            continue
        key = (status_name(status), parity)
        found_exceptions.add(key)
        if key in RADIAL_EXCEPTIONS:
            radial += 1
        elif key in CROSSED_EXCEPTIONS:
            crossed += 1
        elif key in COMBINED_EXCEPTIONS:
            combined += 1
        else:
            unresolved.append(key)

    expected = RADIAL_EXCEPTIONS | CROSSED_EXCEPTIONS | COMBINED_EXCEPTIONS
    if found_exceptions != expected:
        raise AssertionError(
            f"exception table mismatch missing={sorted(expected-found_exceptions)} "
            f"extra={sorted(found_exceptions-expected)}"
        )
    if unresolved:
        raise AssertionError(f"unresolved regimes: {unresolved}")
    if (total, direct, radial, crossed, combined) != (272, 258, 3, 10, 1):
        raise AssertionError(
            f"unexpected coverage counts: {total=} {direct=} {radial=} "
            f"{crossed=} {combined=}"
        )

    print(
        "SU2_LEVEL4_ALL_EXPONENT_TRANSPORT PASS "
        f"regimes={total} direct_hall={direct} radial={radial} "
        f"crossed={crossed} combined={combined} unresolved=0"
    )
    print(
        "RADIAL_TEMPLATE PASS inequality=2*3^n+12^n>=3*4^n "
        "reason=n=0 equality; n>=1 uses 12^n>=3*4^n"
    )
    print(
        "CROSSED_TEMPLATE PASS L=4+2sqrt(3) S=4-2sqrt(3) LS=4 "
        "inequality=L^u*S^v+S^u*L^v>=2"
    )


if __name__ == "__main__":
    main()
