#!/usr/bin/env python3
"""Generate and verify the reader-facing full-Q3 classification ledger."""

from __future__ import annotations

import argparse
import csv
import io
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
HEADER = ROOT / "character_ring_iter" / "full_q3_bcd_remaining_data.hpp"
LEDGER = ROOT / "full_q3_classification_ledger.csv"

ROW = re.compile(
    r"\{'([BCD])',\s*(\d+),\s*(\d+),\s*(\d+),\s*TailMethod::([a-z_]+)\}"
)

FIELDS = (
    "scope",
    "family",
    "rank",
    "representative",
    "stable_through",
    "finite_first_degree",
    "finite_last_degree",
    "tail_onset",
    "maximum_moment_consumed",
    "finite_pairs",
    "tail_method",
    "finite_supplier",
    "tail_supplier",
)


def odd_degrees(first: int, last: int) -> list[int]:
    return [degree for degree in range(first, last + 1) if degree % 2 == 1]


def pairs_at_degree(degree: int) -> int:
    # L=2a+n, with a>=2 and n>=1 odd.
    return max(0, (degree - 3) // 2)


def finite_pair_count(first: int, last: int) -> int:
    return sum(pairs_at_degree(degree) for degree in odd_degrees(first, last))


def record(**values: object) -> dict[str, str]:
    missing = set(FIELDS) - set(values)
    extra = set(values) - set(FIELDS)
    if missing or extra:
        raise ValueError(f"bad ledger record missing={sorted(missing)} extra={sorted(extra)}")
    return {field: str(values[field]) for field in FIELDS}


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    fixed_a = (
        # rank, representative, stable, finite first/last, tail, max moment,
        # exclusive finite count, exact checker count, supplier
        (1, "SU(2)", 2, 5, 27, 29, 27, 78, "SU2 exact recurrence"),
        (2, "SU(3)", 3, 5, 27, 29, 54, 78, "SU3 hook-length/GMP"),
        (3, "SU(4)", 4, 5, 25, 27, 94, 66, "SU4 hook-length/GMP"),
        (4, "SU(5)", 5, 7, 35, 37, 96, 135, "SU5 hook-length/GMP"),
    )
    for rank, representative, stable, first, last, tail, moment, count, supplier in fixed_a:
        assert finite_pair_count(first, last) == count
        rows.append(
            record(
                scope="fixed",
                family="A",
                rank=rank,
                representative=representative,
                stable_through=stable,
                finite_first_degree=first,
                finite_last_degree=last,
                tail_onset=tail,
                maximum_moment_consumed=moment,
                finite_pairs=count,
                tail_method="exact Haar/moment cap",
                finite_supplier=supplier,
                tail_supplier="Part III fixed-group cap tail",
            )
        )

    a_residual_total = 0
    for n_value in range(6, 29):
        rank = n_value - 1
        tail = 27 if n_value == 6 else 29 if n_value in (7, 8) else 31
        first = n_value + 1 if n_value % 2 == 0 else n_value + 2
        last = tail - 2
        count = finite_pair_count(first, last)
        a_residual_total += count
        rows.append(
            record(
                scope="fixed",
                family="A",
                rank=rank,
                representative=f"SU({n_value})",
                stable_through=n_value,
                finite_first_degree=first,
                finite_last_degree=last,
                tail_onset=tail,
                maximum_moment_consumed=32,
                finite_pairs=count,
                tail_method="uniform Paley--Zygmund cap",
                finite_supplier="SUN_GE6 exact Schur--Weyl/GMP",
                tail_supplier="Part III uniform type-A tail",
            )
        )
    assert a_residual_total == 1315
    rows.append(
        record(
            scope="parametric",
            family="A",
            rank=">=28",
            representative="SU(N), N>=29",
            stable_through="N",
            finite_first_degree="none",
            finite_last_degree="none",
            tail_onset="31 (overlaps stable range)",
            maximum_moment_consumed=32,
            finite_pairs=0,
            tail_method="stable law + uniform cap",
            finite_supplier="none",
            tail_supplier="Part III stable and uniform type-A tails",
        )
    )

    exceptional = (
        ("G2", 25, 38, 55),
        ("F4", 37, 50, 136),
        ("E6", 29, 42, 78),
        ("E7", 63, 70, 435),
        ("E8", 71, 94, 561),
    )
    for family, tail, moment, count in exceptional:
        assert finite_pair_count(5, tail - 2) == count
        rows.append(
            record(
                scope="fixed",
                family=family,
                rank="-",
                representative=family,
                stable_through="none",
                finite_first_degree=5,
                finite_last_degree=tail - 2,
                tail_onset=tail,
                maximum_moment_consumed=moment,
                finite_pairs=count,
                tail_method="exact polynomial cap",
                finite_supplier="exceptional Racah--Speiser/GMP",
                tail_supplier="Part III exceptional polynomial tail",
            )
        )

    header_text = HEADER.read_text(encoding="utf-8")
    matches = ROW.findall(header_text)
    if len(matches) != 114:
        raise ValueError(f"expected 114 classical rows, found {len(matches)}")

    classical_total = 0
    method_name = {
        "polynomial": "exact polynomial cap",
        "rational_cap": "exact rational Weyl cap",
        "directed_interval": "directed exact-MGF/Fredholm tail",
    }
    for family, rank_text, tail_text, moment_text, method in matches:
        rank = int(rank_text)
        tail = int(tail_text)
        moment = int(moment_text)
        stable = 2 * rank + 1 if family in "BC" else rank - 1
        first = stable + 1 if stable % 2 == 0 else stable + 2
        last = tail - 2
        count = finite_pair_count(first, last)
        classical_total += count
        rows.append(
            record(
                scope="fixed",
                family=family,
                rank=rank,
                representative=f"{family}_{rank}",
                stable_through=stable,
                finite_first_degree=first,
                finite_last_degree=last,
                tail_onset=tail,
                maximum_moment_consumed=moment,
                finite_pairs=count,
                tail_method=method_name[method],
                finite_supplier="bounded-Littlewood determinant/GMP",
                tail_supplier="Part III exact or directed classical tail",
            )
        )
    assert classical_total == 17862

    for family, rank, stable, onset, supplier in (
        ("B", ">=22", "2b+1", "2b+3", "finite-middle/high-rank trace tail"),
        ("C", ">=29", "2b+1", "2b+3", "finite-middle/high-rank trace tail"),
        ("D", ">=71", "b-1", "b+1 (b even), b (b odd)", "finite-middle/high-rank trace tail"),
    ):
        rows.append(
            record(
                scope="parametric",
                family=family,
                rank=rank,
                representative=f"{family}_b",
                stable_through=stable,
                finite_first_degree="none",
                finite_last_degree="none",
                tail_onset=onset,
                maximum_moment_consumed="none",
                finite_pairs=0,
                tail_method="directed Fredholm/trace tail",
                finite_supplier="none",
                tail_supplier=supplier,
            )
        )

    return rows


def render(rows: list[dict[str, str]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="fail unless the committed ledger is current")
    args = parser.parse_args()

    rows = build_rows()
    expected = render(rows)
    if args.check:
        if not LEDGER.is_file() or LEDGER.read_text(encoding="utf-8") != expected:
            raise SystemExit("FULL_Q3_CLASSIFICATION_LEDGER: FAIL: regenerate the CSV")
    else:
        LEDGER.write_text(expected, encoding="utf-8")

    fixed = sum(row["scope"] == "fixed" for row in rows)
    parametric = len(rows) - fixed
    print(
        "FULL_Q3_CLASSIFICATION_LEDGER "
        f"rows={len(rows)} fixed={fixed} parametric={parametric} "
        "classical_pairs=17862 type_a_pairs=1315 exceptional_pairs=1265"
    )
    print("FULL_Q3_CLASSIFICATION_LEDGER VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
