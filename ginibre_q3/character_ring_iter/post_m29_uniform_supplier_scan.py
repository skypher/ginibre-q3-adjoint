#!/usr/bin/env python3
"""Diagnostics for post-m=29 uniform supplier candidates.

The active target is the finite-row replacement `CMO-PostM29FiniteTail`.
This script keeps the exploratory supplier checks reproducible:

* the currently open finite residue after the direct post-m=29 rows and the
  top finite trace rows are removed;
* the corrected fixed-root-angle cap check, using an inscribed radial cap;
* the CJL/Rains Chernoff tau-source diagnostic for the trace route, replacing
  the current Rains-square additive total-variation term by the source-supported
  no-TV concentration term.

The printed inequalities are floating diagnostics.  They are meant to identify
which named estimate would shrink the finite tail; they are not proof
certificates.  The source-audited certificate is
`classical_tail_constants.py --post-m29-chernoff-trace-certificate`.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from fractions import Fraction
from typing import Iterable

import mpmath as mp

import classical_tail_constants as ctc


def pre_chernoff_rows() -> list[tuple[str, int]]:
    closed = ctc.post_m29_all_later_direct_closed_rows() | ctc.first_hit_trace_closed_rows()
    return [row for row in ctc.finite_residue_rows(296) if row not in closed]


def current_residue_rows() -> list[tuple[str, int]]:
    return ctc.post_m29_remaining_residue(296)


def row_summary(rows: Iterable[tuple[str, int]]) -> str:
    by_family: dict[str, list[int]] = defaultdict(list)
    for family, rank in rows:
        by_family[family].append(rank)
    parts = []
    for family in "BCD":
        ranks = by_family[family]
        parts.append(f"{family}: {ctc.rank_range_text(ranks)} ({len(ranks)})")
    return "; ".join(parts)


def gamma_p(shape: mp.mpf, radius: mp.mpf) -> mp.mpf:
    if radius <= 0:
        return mp.mpf("0")
    return mp.gammainc(shape, 0, radius, regularized=True)


def wronskian_fraction(shape: mp.mpf, radius: mp.mpf) -> mp.mpf:
    p0 = gamma_p(shape, radius)
    p1 = gamma_p(shape + 1, radius)
    p2 = gamma_p(shape + 2, radius)
    return (shape + 1) * p0 * p2 - shape * p1 * p1


def max_root_square(family: str) -> int:
    return 4 if family == "C" else 2


def inscribed_root_cap_slack(family: str, rank: int, delta: float) -> float:
    """Log10 slack for the corrected fixed-root-angle cap diagnostic.

    If `Q(x)/(2d) <= delta^2 kappa n/(2d Lmax)`, then the scaled radial ball
    is inside the root-angle cap `|alpha(theta)| <= delta`.
    """

    data = ctc.family_data(family, rank)
    n_value = max(63, ctc.stable_odd_reach(family, rank) + 2)
    dimension = int(data["d"])
    kappa = int(data["kappa"])
    r_plus = (dimension - int(data["r"])) // 2
    shape = mp.mpf(dimension) / 2
    radius = (
        mp.mpf(delta)
        * mp.mpf(delta)
        * kappa
        * n_value
        / (2 * dimension * max_root_square(family))
    )
    fraction = wronskian_fraction(shape, radius)
    fraction_log10 = float(mp.log10(fraction)) if fraction > 0 else -1.0e300
    sine_ratio_log10 = math.log10(2 * math.sin(delta / 2) / delta)
    required_log10 = -ctc.leading_margin(data, n_value) / math.log(10)
    return fraction_log10 + 4 * r_plus * sine_ratio_log10 - required_log10


def positive_trace_log(C: int, r_value: float, eta_value: float) -> float | None:
    L_value = 2 * r_value + 3 * eta_value
    gaussian_log = (
        -L_value * C / 2
        + 0.5 * math.log(L_value * C)
        - math.log(L_value * C + 1)
        - 0.5 * math.log(2)
    )
    if C < 132:
        return None
    tv_log = (
        math.log(5)
        + 0.25 * math.log(math.log(2 * C))
        - C * math.log(2 * C / math.e)
    )
    if tv_log >= gaussian_log:
        return None
    return gaussian_log + math.log1p(-math.exp(tv_log - gaussian_log))


def chernoff_trace_margin(
    family: str,
    rank: int,
    r_value: Fraction,
    eta_value: Fraction,
) -> tuple[str, float | None]:
    """Trace-route log10 margin using the checked Decimal certificate test."""

    reason, margin = ctc.finite_trace_pushforward_margin_chernoff_tau(
        ctc.c_value(family, rank),
        max(63, ctc.stable_odd_reach(family, rank) + 2),
        r_value,
        eta_value,
    )
    if margin is None:
        return reason, None
    return reason, float(margin / ctc.Decimal(10).ln())


def print_open_rows(rows: list[tuple[str, int]], label: str) -> None:
    print(label)
    print(f"  count: {len(rows)}")
    print(f"  {row_summary(rows)}")


def print_cap_scan(rows: list[tuple[str, int]]) -> None:
    print("corrected inscribed root-angle cap scan")
    print("  diagnostic: negative worst slack means this supplier cannot close")
    mp.mp.dps = 80
    for delta in (2.0, 2.3, 2.5, 2.8, 3.0):
        worst: tuple[float, str, int] | None = None
        bad = 0
        for family, rank in rows:
            slack = inscribed_root_cap_slack(family, rank, delta)
            if slack < 0:
                bad += 1
            if worst is None or slack < worst[0]:
                worst = (slack, family, rank)
        assert worst is not None
        print(
            f"  delta={delta:.1f}: bad={bad}, "
            f"worst={worst[1]}_{worst[2]} slack={worst[0]:.2f}"
        )


def print_sharp_tau_scan(rows: list[tuple[str, int]]) -> None:
    print("CJL/Rains Chernoff tau-source trace diagnostic")
    print(
        "  supplier: replace the Rains-square additive tau-TV term "
        "by the no-TV Chernoff concentration term"
    )
    parameter_rows = {
        "B": (
            ("even-high", Fraction(2001, 1000), Fraction(9, 20)),
            ("odd-high", Fraction(10001, 5000), Fraction(56, 125)),
            ("odd-low-219", Fraction(99, 50), Fraction(56, 125)),
        ),
        "C": (
            ("even-high", Fraction(2001, 1000), Fraction(9, 20)),
            ("odd-high", Fraction(10001, 5000), Fraction(56, 125)),
            ("odd-low-219", Fraction(99, 50), Fraction(56, 125)),
        ),
        "D": (("d-probe", Fraction(2001, 1000), Fraction(11, 20)),),
    }
    ok_rows: list[tuple[str, int, str, Fraction, Fraction, float]] = []
    reason_counts: Counter[str] = Counter()
    for family, rank in rows:
        failures: list[str] = []
        for label, r_value, eta_value in parameter_rows[family]:
            reason, margin = chernoff_trace_margin(family, rank, r_value, eta_value)
            if reason == "OK" and margin is not None and margin >= 0:
                ok_rows.append((family, rank, label, r_value, eta_value, margin))
                break
            failures.append(reason)
        else:
            reason_counts[failures[-1]] += 1
    print(f"  rows covered by this supplier: {len(ok_rows)}")
    for family in "BCD":
        family_rows = [
            (rank, label, r_value, eta_value, margin)
            for fam, rank, label, r_value, eta_value, margin in ok_rows
            if fam == family
        ]
        ranks = [rank for rank, *_rest in family_rows]
        min_margin = min((row[4] for row in family_rows), default=None)
        margin_text = "none" if min_margin is None else f"{min_margin:.2f}"
        print(
            f"  {family}: ranks {ctc.rank_range_text(ranks)} ({len(ranks)}), "
            f"min log10 margin {margin_text}"
        )
        for label in sorted({label for _rank, label, *_rest in family_rows}):
            label_rows = [row for row in family_rows if row[1] == label]
            label_ranks = [rank for rank, *_rest in label_rows]
            _rank, _label, r_value, eta_value, _margin = label_rows[0]
            print(
                f"    {label}: r={r_value}, eta={eta_value}, "
                f"ranks {ctc.rank_range_text(label_ranks)} ({len(label_ranks)})"
            )
    if reason_counts:
        print(f"  non-closed first-failure counts: {dict(reason_counts)}")
    closed = {(family, rank) for family, rank, *_rest in ok_rows}
    remaining = [row for row in rows if row not in closed]
    print(f"  remaining rows after this supplier: {len(remaining)}")
    print(f"  {row_summary(remaining)}")


def main() -> None:
    pre_rows = pre_chernoff_rows()
    print_open_rows(pre_rows, "post-m29 pre-Chernoff finite rows")
    print_sharp_tau_scan(pre_rows)
    rows = current_residue_rows()
    print_open_rows(rows, "post-m29 current residual finite rows")
    print_cap_scan(rows)


if __name__ == "__main__":
    main()
