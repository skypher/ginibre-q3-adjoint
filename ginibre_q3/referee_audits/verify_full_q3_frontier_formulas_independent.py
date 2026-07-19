#!/usr/bin/env python3
"""Independently recompute the formula fields in all H25--H29 certificates.

This is a supplementary referee audit, not a replacement for the distributed
replay.  It deliberately uses Python integers/Fraction rather than the GMP C++
implementations.  The reverse-Pieri totals themselves are not recomputed; the
audit checks their initial terminal partition counts, stable values, and every
printed margin.  The separate horizontal-strip oracle checks the transition.
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


def fail(message: str) -> None:
    raise RuntimeError(message)


def expect(actual, expected, label: str) -> None:
    if actual != expected:
        fail(f"{label}: expected={expected!r}, actual={actual!r}")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        fail(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def integer_fields(line: str) -> dict[str, int]:
    return {key: int(value) for key, value in re.findall(r"(\w+)=(-?\d+)", line)}


def stable_moments(maximum: int) -> list[int]:
    values = [0] * (maximum + 1)
    values[0] = 1
    for n in range(1, maximum):
        values[n + 1] = (
            n * values[n]
            + n * values[n - 1]
            - (n * (n - 1) // 2) * (values[n - 2] if n >= 2 else 0)
        )
        if values[n + 1] < 0:
            fail(f"independent stable recurrence is negative at {n + 1}")
    return values


def stable_moments_from_closed_egf(maximum: int) -> list[int]:
    """Extract coefficients of (1-t)^(-1/2) exp(-t/2+t^2/4)."""

    inverse_square_root = [
        Fraction(comb(2 * n, n), 4**n) for n in range(maximum + 1)
    ]
    exponential: list[Fraction] = []
    for n in range(maximum + 1):
        coefficient = Fraction()
        for quadratic_degree in range(n // 2 + 1):
            linear_degree = n - 2 * quadratic_degree
            coefficient += Fraction(
                (-1) ** linear_degree,
                2**linear_degree
                * factorial(linear_degree)
                * 4**quadratic_degree
                * factorial(quadratic_degree),
            )
        exponential.append(coefficient)

    values: list[int] = []
    for n in range(maximum + 1):
        coefficient = sum(
            (
                inverse_square_root[left_degree]
                * exponential[n - left_degree]
                for left_degree in range(n + 1)
            ),
            Fraction(),
        )
        moment = coefficient * factorial(n)
        if moment.denominator != 1:
            fail(f"closed stable EGF has nonintegral coefficient at degree {n}")
        values.append(moment.numerator)
    return values


def partition_count(n: int, maximum_part: int | None = None) -> int:
    bound = n if maximum_part is None else maximum_part
    values = [0] * (n + 1)
    values[0] = 1
    for part in range(1, bound + 1):
        for total in range(part, n + 1):
            values[total] += values[total - part]
    return values[n]


def odd_double_factorial(n: int) -> int:
    value = 1
    for factor in range(1, n + 1, 2):
        value *= factor
    return value


def lower_box_ratio(h: int, quotient_floor: int) -> Fraction:
    lower_boxes = 2 * h
    one_choices = 2 * lower_boxes
    two_choices = lower_boxes * (lower_boxes + 1) // 2
    value = Fraction()
    for a in range(lower_boxes + 1):
        for c in range(lower_boxes + 1):
            if a + 2 * c > lower_boxes:
                continue
            value += Fraction(
                factorial(lower_boxes)
                * 3 ** (lower_boxes - a - c)
                * two_choices**c,
                factorial(a)
                * factorial(c)
                * quotient_floor ** (lower_boxes - a - c)
                * one_choices ** (lower_boxes - a),
            )
    return value


def pause_polynomial(h: int, q: int) -> int:
    j = 2 * q + h
    return sum(
        2 ** (h + pauses) * comb(j, pauses)
        for pauses in range(h % 2, h + 1, 2)
    )


def fpf_bound(h: int, q: int) -> int:
    j = 2 * q + h
    return sum(
        comb(j, pauses) * odd_double_factorial(j + pauses - 1)
        for pauses in range(h % 2, h + 1, 2)
    )


class Audit:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.verifier = load_module(
            root / "verify_full_q3_distributed_replay.py",
            "full_q3_formula_referee_verifier",
        )
        self.ledger = load_module(
            root / "run_full_q3_frontier_fleet_queue.py",
            "full_q3_formula_referee_ledger",
        )
        self.tasks = self.ledger.make_tasks(8)
        expect(len(self.tasks), 71, "fleet task count")
        expect(len({task.name for task in self.tasks}), 71, "unique fleet tasks")
        self.accepted_manifest = self.verifier.parse_manifest(
            root / "certificates/post_m29/post_m29_accepted_manifest.sha256"
        )
        self.block_cache: dict[str, dict[str, tuple[str, ...]]] = {}
        self.ratio_cache: dict[tuple[int, int], Fraction] = {}
        recurrence_stable = stable_moments(340)
        closed_egf_stable = stable_moments_from_closed_egf(340)
        expect(
            recurrence_stable,
            closed_egf_stable,
            "stable recurrence/closed-EGF coefficient agreement",
        )
        self.stable = closed_egf_stable
        self.counts = {
            "base": 0,
            "b": 0,
            "reverse": 0,
            "terminal": 0,
            "fpf": 0,
            "c": 0,
            "variable": 0,
            "margins": 0,
        }

    def ratio(self, h: int, quotient_floor: int) -> Fraction:
        key = (h, quotient_floor)
        if key not in self.ratio_cache:
            self.ratio_cache[key] = lower_box_ratio(h, quotient_floor)
        return self.ratio_cache[key]

    def rows(self, task) -> list[str]:
        certificate_name = self.verifier.frontier_certificate_name(task.name)
        if certificate_name not in self.block_cache:
            path = self.root / "certificates/post_m29" / certificate_name
            expect(
                self.verifier.digest(path),
                self.accepted_manifest[certificate_name],
                f"accepted identity {certificate_name}",
            )
            self.block_cache[certificate_name] = (
                self.verifier.accepted_arithmetic_blocks(path.read_text(), path)
            )
        scope = self.verifier.normalized_fleet_scope(
            self.verifier.expected_fleet_scope_line(task, 8)
        )
        blocks = self.block_cache[certificate_name]
        if scope not in blocks:
            fail(f"missing accepted scope for {task.name}: {scope}")
        return list(blocks[scope])

    def audit_standard(self, task, h: int, b_absorption_q: int, c_absorption_q: int) -> None:
        lower_boxes = 2 * h
        one_choices = 2 * lower_boxes
        ratio = self.ratio(h, 31 - h)
        ratio_ceil = (ratio.numerator + ratio.denominator - 1) // ratio.denominator
        cases: dict[tuple[int, int], dict[str, int]] = {}
        bad: dict[tuple[int, int], int] = {}
        stable: dict[tuple[int, int], int] = {}
        margins: dict[tuple[int, int], int] = {}
        fpf_bounds: dict[tuple[int, int], int] = {}
        fpf_stable: dict[tuple[int, int], int] = {}
        fpf_margins: dict[tuple[int, int], int] = {}
        pauses: dict[tuple[int, int], int] = {}
        c_bounds: dict[tuple[int, int], int] = {}
        c_stable: dict[tuple[int, int], int] = {}
        c_margins: dict[tuple[int, int], int] = {}
        reverse_rows: list[dict[str, int]] = []
        local_base_count = 0

        for line in self.rows(task):
            fields = integer_fields(line)
            if line.startswith("B_REVERSE_STEP "):
                reverse_rows.append(fields)
            elif line.startswith("B_LOWER_BOX_RATIONAL "):
                expect(
                    (fields["numerator"], fields["denominator"], fields["ceil"]),
                    (ratio.numerator, ratio.denominator, ratio_ceil),
                    f"{task.name} lower ratio",
                )
                local_base_count += 1
            elif line.startswith("B_ABSORPTION_BASE "):
                q = b_absorption_q
                left = (
                    4
                    * ratio_ceil
                    * comb(2 * q + h, lower_boxes)
                    * one_choices**lower_boxes
                    * 3 ** (2 * q - h)
                )
                right = factorial(2 * q + h - 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, right - left),
                    f"{task.name} B absorption",
                )
                local_base_count += 1
            elif line.startswith("B_ABSORPTION_RATIO_BASE "):
                q = b_absorption_q
                left = 9 * (2 * q + h + 2) * (2 * q + h + 1)
                right = 30 * (2 * q - h + 2) * (2 * q - h + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, right - left),
                    f"{task.name} B ratio",
                )
                local_base_count += 1
            elif line.startswith("C_PAUSE_ABSORPTION_BASE "):
                q = c_absorption_q
                p_value = pause_polynomial(h, q)
                central = comb(2 * q + h, q)
                expect(
                    (fields["q"], fields["P"], fields["central"], fields["margin"]),
                    (q, p_value, central, central - p_value),
                    f"{task.name} C absorption",
                )
                local_base_count += 1
            elif line.startswith("C_PAUSE_RATIO_BASE "):
                q = c_absorption_q
                left = (2 * q + h + 2) * (2 * q + h + 1)
                right = 2 * (2 * q + 2) * (2 * q + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, right - left),
                    f"{task.name} C pause ratio",
                )
                local_base_count += 1
            elif line.startswith("C_CENTRAL_RATIO_BASE "):
                q = c_absorption_q
                left = (2 * q + h + 2) * (2 * q + h + 1)
                right = 3 * (q + 1) * (q + h + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, left - right),
                    f"{task.name} C central ratio",
                )
                local_base_count += 1
            elif line.startswith("B_CASE "):
                cases[fields["q"], fields["j"]] = fields
            elif line.startswith("B_BAD_COUNT "):
                bad[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("B_STABLE "):
                stable[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("B_MARGIN "):
                margins[fields["q"], fields["j"]] = fields["stable_minus_2bad"]
            elif line.startswith("C_FPF_BOUND "):
                fpf_bounds[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_FPF_STABLE "):
                fpf_stable[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_FPF_MARGIN "):
                fpf_margins[fields["q"], fields["j"]] = fields["stable_minus_2bound"]
            elif line.startswith("C_PAUSE_POLYNOMIAL "):
                pauses[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_BOUND "):
                c_bounds[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_STABLE "):
                c_stable[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_MARGIN "):
                c_margins[fields["q"], fields["j"]] = fields["stable_minus_bound"]
            else:
                fail(f"unclassified line in {task.name}: {line}")

        expect(local_base_count, 6, f"{task.name} base-record count")
        self.counts["base"] += local_base_count
        expect(set(cases), set(bad), f"{task.name} bad-count keys")
        expect(set(cases), set(stable), f"{task.name} stable keys")
        expect(set(cases), set(margins), f"{task.name} margin keys")
        case_q = {q for q, _ in cases}
        for fields in reverse_rows:
            if (
                fields["q"] not in case_q
                or fields["states"] <= 0
                or not 0 <= fields["after_removing_to"] < 2 * fields["q"] + h
            ):
                fail(f"invalid reverse-progress record in {task.name}: {fields}")
        self.counts["reverse"] += len(reverse_rows)
        for key, fields in cases.items():
            q, j = key
            expect(j, 2 * q + h, f"{task.name} B degree")
            expect(fields["width_wall"], 2 * q, f"{task.name} width wall")
            expect(
                fields["terminal_shapes"],
                partition_count(j) - partition_count(j, 2 * q - 1),
                f"{task.name} terminal shapes",
            )
            expect(stable[key], self.stable[j], f"{task.name} B stable")
            expect(margins[key], self.stable[j] - 2 * bad[key], f"{task.name} B margin")
            if margins[key] < 0:
                fail(f"negative B margin in {task.name}: {key}")
            self.counts["b"] += 1
            self.counts["terminal"] += 1
            self.counts["margins"] += 1

        expect(set(fpf_bounds), set(fpf_stable), f"{task.name} FPF stable keys")
        expect(set(fpf_bounds), set(fpf_margins), f"{task.name} FPF margin keys")
        for key, bound in fpf_bounds.items():
            q, j = key
            expected_bound = fpf_bound(h, q)
            expect(j, 2 * q + h, f"{task.name} FPF degree")
            expect(bound, expected_bound, f"{task.name} FPF bound")
            expect(fpf_stable[key], self.stable[j], f"{task.name} FPF stable")
            expect(
                fpf_margins[key],
                self.stable[j] - 2 * expected_bound,
                f"{task.name} FPF margin",
            )
            if fpf_margins[key] < 0:
                fail(f"negative FPF margin in {task.name}: {key}")
            self.counts["fpf"] += 1
            self.counts["margins"] += 1

        expect(set(pauses), set(c_bounds), f"{task.name} C-bound keys")
        expect(set(pauses), set(c_stable), f"{task.name} C-stable keys")
        expect(set(pauses), set(c_margins), f"{task.name} C-margin keys")
        for key, p_value in pauses.items():
            q, j = key
            expected_p = pause_polynomial(h, q)
            expected_bound = 2 ** (2 * q + 1) * expected_p * self.stable[q + h]
            expect(j, 2 * q + h, f"{task.name} C degree")
            expect(p_value, expected_p, f"{task.name} pause polynomial")
            expect(c_bounds[key], expected_bound, f"{task.name} C bound")
            expect(c_stable[key], self.stable[j], f"{task.name} C stable")
            expect(c_margins[key], self.stable[j] - expected_bound, f"{task.name} C margin")
            if c_margins[key] < 0:
                fail(f"negative C margin in {task.name}: {key}")
            self.counts["c"] += 1
            self.counts["margins"] += 1

    def audit_absorption(self, task, h: int, q_min: int, q_max: int, fixed_q: int) -> None:
        lower_boxes = 2 * h
        one_choices = 2 * lower_boxes
        fixed_ratio = self.ratio(h, 31 - h)
        fixed_ceil = (
            fixed_ratio.numerator + fixed_ratio.denominator - 1
        ) // fixed_ratio.denominator
        seen_q: list[int] = []
        for line in self.rows(task):
            fields = integer_fields(line)
            if line.startswith("B_VARIABLE_ABSORPTION "):
                q = fields["q"]
                j = 2 * q + h
                ratio = self.ratio(h, 2 * q - h + 1)
                left = (
                    2
                    * ratio.numerator
                    * comb(j, lower_boxes)
                    * one_choices**lower_boxes
                    * 3 ** (2 * q - h)
                )
                right = self.stable[j] * ratio.denominator
                expect(
                    (
                        fields["j"],
                        fields["quotient_floor"],
                        fields["numerator"],
                        fields["denominator"],
                        fields["left_num"],
                        fields["right_num"],
                        fields["margin_num"],
                    ),
                    (
                        j,
                        2 * q - h + 1,
                        ratio.numerator,
                        ratio.denominator,
                        left,
                        right,
                        right - left,
                    ),
                    f"{task.name} variable absorption",
                )
                if right < left:
                    fail(f"negative variable absorption in {task.name}: q={q}")
                seen_q.append(q)
                self.counts["variable"] += 1
                self.counts["margins"] += 1
            elif line.startswith("B_FIXED_LOWER_BOX_RATIONAL "):
                expect(
                    (fields["numerator"], fields["denominator"]),
                    (fixed_ratio.numerator, fixed_ratio.denominator),
                    f"{task.name} fixed ratio",
                )
            elif line.startswith("B_FIXED_LOWER_BOX_CEIL "):
                expect(fields["value"], fixed_ceil, f"{task.name} fixed ceil")
            elif line.startswith("B_FIXED_ABSORPTION_BASE "):
                left = (
                    4
                    * fixed_ceil
                    * comb(2 * fixed_q + h, lower_boxes)
                    * one_choices**lower_boxes
                    * 3 ** (2 * fixed_q - h)
                )
                right = factorial(2 * fixed_q + h - 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (fixed_q, left, right, right - left),
                    f"{task.name} fixed absorption",
                )
                if right < left:
                    fail(f"negative fixed absorption in {task.name}")
                self.counts["margins"] += 1
            elif line.startswith("B_FIXED_ABSORPTION_RATIO_BASE "):
                left = 9 * (2 * fixed_q + h + 2) * (2 * fixed_q + h + 1)
                right = 30 * (2 * fixed_q - h + 2) * (2 * fixed_q - h + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (fixed_q, left, right, right - left),
                    f"{task.name} fixed absorption ratio",
                )
            else:
                fail(f"unclassified line in {task.name}: {line}")
        expect(seen_q, list(range(q_min, q_max + 1)), f"{task.name} q coverage")

    def audit_h29_c(self, task) -> None:
        h = 29
        c_absorption_q = 74
        fpf_bounds: dict[tuple[int, int], int] = {}
        fpf_stable: dict[tuple[int, int], int] = {}
        fpf_margins: dict[tuple[int, int], int] = {}
        pauses: dict[tuple[int, int], int] = {}
        c_bounds: dict[tuple[int, int], int] = {}
        c_stable: dict[tuple[int, int], int] = {}
        c_margins: dict[tuple[int, int], int] = {}
        base_count = 0
        for line in self.rows(task):
            fields = integer_fields(line)
            if line.startswith("C_PAUSE_ABSORPTION_BASE "):
                q = c_absorption_q
                p_value = pause_polynomial(h, q)
                central = comb(2 * q + h, q)
                expect(
                    (fields["q"], fields["P"], fields["central"], fields["margin"]),
                    (q, p_value, central, central - p_value),
                    f"{task.name} C absorption",
                )
                base_count += 1
            elif line.startswith("C_PAUSE_RATIO_BASE "):
                q = c_absorption_q
                left = (2 * q + h + 2) * (2 * q + h + 1)
                right = 2 * (2 * q + 2) * (2 * q + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, right - left),
                    f"{task.name} C pause ratio",
                )
                base_count += 1
            elif line.startswith("C_CENTRAL_RATIO_BASE "):
                q = c_absorption_q
                left = (2 * q + h + 2) * (2 * q + h + 1)
                right = 3 * (q + 1) * (q + h + 1)
                expect(
                    (fields["q"], fields["left"], fields["right"], fields["margin"]),
                    (q, left, right, left - right),
                    f"{task.name} C central ratio",
                )
                base_count += 1
            elif line.startswith("C_FPF_BOUND "):
                fpf_bounds[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_FPF_STABLE "):
                fpf_stable[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_FPF_MARGIN "):
                fpf_margins[fields["q"], fields["j"]] = fields["stable_minus_2bound"]
            elif line.startswith("C_PAUSE_POLYNOMIAL "):
                pauses[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_BOUND "):
                c_bounds[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_STABLE "):
                c_stable[fields["q"], fields["j"]] = fields["value"]
            elif line.startswith("C_MARGIN "):
                c_margins[fields["q"], fields["j"]] = fields["stable_minus_bound"]
            else:
                fail(f"unclassified line in {task.name}: {line}")
        expect(base_count, 3, f"{task.name} base count")
        expect(set(fpf_bounds), set(fpf_stable), f"{task.name} FPF stable keys")
        expect(set(fpf_bounds), set(fpf_margins), f"{task.name} FPF margin keys")
        for key, bound in fpf_bounds.items():
            q, j = key
            expected_bound = fpf_bound(h, q)
            expect(bound, expected_bound, f"{task.name} FPF bound")
            expect(fpf_stable[key], self.stable[j], f"{task.name} FPF stable")
            expect(
                fpf_margins[key],
                self.stable[j] - 2 * expected_bound,
                f"{task.name} FPF margin",
            )
            if fpf_margins[key] < 0:
                fail(f"negative FPF margin in {task.name}: {key}")
            self.counts["fpf"] += 1
            self.counts["margins"] += 1
        expect(set(pauses), set(c_bounds), f"{task.name} C-bound keys")
        expect(set(pauses), set(c_stable), f"{task.name} C-stable keys")
        expect(set(pauses), set(c_margins), f"{task.name} C-margin keys")
        for key, p_value in pauses.items():
            q, j = key
            expected_p = pause_polynomial(h, q)
            expected_bound = 2 ** (2 * q + 1) * expected_p * self.stable[q + h]
            expect(p_value, expected_p, f"{task.name} pause polynomial")
            expect(c_bounds[key], expected_bound, f"{task.name} C bound")
            expect(c_stable[key], self.stable[j], f"{task.name} C stable")
            expect(c_margins[key], self.stable[j] - expected_bound, f"{task.name} C margin")
            if c_margins[key] < 0:
                fail(f"negative C margin in {task.name}: {key}")
            self.counts["c"] += 1
            self.counts["margins"] += 1

    def run(self) -> None:
        standard = {
            "fifth": (25, 44, 64),
            "sixth": (26, 46, 66),
            "seventh": (27, 48, 69),
            "eighth": (28, 50, 71),
        }
        for task in self.tasks[:68]:
            match = re.match(r"twenty(fifth|sixth|seventh|eighth)_", task.name)
            if match is None:
                fail(f"unexpected standard task: {task.name}")
            self.audit_standard(task, *standard[match.group(1)])
        self.audit_absorption(self.tasks[68], 28, 42, 49, 50)
        self.audit_h29_c(self.tasks[69])
        self.audit_absorption(self.tasks[70], 29, 43, 52, 53)
        expected = {
            "base": 408,
            "b": 120,
            "reverse": 2344,
            "terminal": 120,
            "fpf": 60,
            "c": 179,
            "variable": 18,
            "margins": 379,
        }
        expect(self.counts, expected, "global formula-audit counts")
        print(
            "FULL_Q3_FRONTIER_FORMULA_AUDIT "
            "tasks=71 standard_tasks=68 base_records=408 B_cases=120 "
            "reverse_progress_rows=2344 terminal_partition_counts=120 "
            "C_fpf_cases=60 C_arithmetic_cases=179 "
            "variable_absorption_cases=18 exact_margins=379 "
            "stable_closed_egf_values=341 failures=0"
        )
        print("FULL_Q3_FRONTIER_FORMULA_AUDIT VERIFICATION: ALL PASS")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="ginibre_q3 source root",
    )
    args = parser.parse_args()
    Audit(args.root.resolve()).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
