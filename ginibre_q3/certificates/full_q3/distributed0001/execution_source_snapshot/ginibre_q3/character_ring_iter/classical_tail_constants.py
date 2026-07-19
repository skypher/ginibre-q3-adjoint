#!/usr/bin/env python3
"""Classical B/C/D positive-saddle leading constants.

This computes the Macdonald-Mehta leading coefficient used by the
positive-saddle tail comparison in CLASSICAL_M3_ONE_PROOF.md.

With long roots normalized to have squared length 2,

    A0 = z^2 * 8*a*d^2*(d/kappa)^d
         * prod_i ((degree_i - 1)!)^2 / (2*pi)^r,

where r is rank, d is dim(G), a=d/2, z=|Z(G) cap T| for the displayed
compact model, and kappa is defined by

    sum_{alpha in R_+} alpha(x)^2 = kappa * |x|^2.

The checker prints a rational lower bound using pi < 4, evaluates the
logarithmic dominance margin at the last odd n reached by the stable Chain
window, and checks the minimum of that leading-term margin over the
post-stable odd range.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial, lgamma, log, log10
from pathlib import Path
from time import time


TRACE_R = Fraction(20000001, 10000000)
TRACE_ETA = Fraction(9973, 50000)
TRACE_A = Fraction(461, 200)
TRACE_L = 2 * TRACE_R + 3 * TRACE_ETA
TRACE_MONOTONE_CUTOFF = 945
TRACE_CUTOFF = TRACE_MONOTONE_CUTOFF
TRACE_SQUARE_R = Fraction(2000001, 1000000)
TRACE_SQUARE_ETA = Fraction(359, 2000)
TRACE_SQUARE_A = Fraction(4571, 2000)
TRACE_SQUARE_L = 2 * TRACE_SQUARE_R + 3 * TRACE_SQUARE_ETA
TRACE_SQUARE_CUTOFF = 296
TRACE_CUTOFF = TRACE_SQUARE_CUTOFF
FIRST_HIT_TRACE_R = Fraction(20001, 10000)
FIRST_HIT_TRACE_ETA = Fraction(1861, 10000)
HERE = Path(__file__).resolve().parent
STABLE_MOMENTS_PATH = HERE.parent / "references" / "oeis_A002137_stable.txt"

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def decimal_fraction(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def tau_rains_correction(c_value: int) -> Decimal:
    """Rains/CJL correction N^(-5/2) exp(-N/2), N=floor(C/2)."""

    n_value = c_value // 2
    if n_value < 2:
        n_value = 2
    n_dec = Decimal(n_value)
    return (n_dec ** Decimal("-2.5")) * (-(n_dec / Decimal(2))).exp()


def sharp_tau_log(c_value: int, eta_value: Fraction = TRACE_ETA) -> Decimal:
    """Log of the retained-correction tau upper bound."""

    c = Decimal(c_value)
    eta = decimal_fraction(eta_value)
    return Decimal(2).ln() - (eta * c - 1) ** 2 / (
        Decimal(16) * (Decimal(1) + tau_rains_correction(c_value))
    )


def trace_cutoff_logs(c_value: int) -> dict[str, Decimal]:
    """High-precision log margins for the concrete trace cutoff.

    The Gaussian row uses the Mills lower bound with the conservative paper
    bound pi < 4.  The tau row retains the Rains/CJL correction instead of
    replacing it by 1.
    """

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(TRACE_R)
    eta = decimal_fraction(TRACE_ETA)
    a_const = decimal_fraction(TRACE_A)
    l_value = decimal_fraction(TRACE_L)
    log_half = Decimal("0.5").ln()
    den = Decimal(1) - (Decimal(1) / (eta * eta * c * c))
    s_log = sharp_tau_log(c_value)

    gaussian = (
        (a_const - l_value / 2) * c
        + (l_value * c).ln() / 2
        - (l_value * c + 1).ln()
        - Decimal(2).ln() / 2
        - Decimal(2).ln()
    )
    r0 = a_const * c + s_log - log_half
    r1 = (
        (Decimal(4) * (Decimal(1) + c ** Decimal(-2)) * (r / 2) ** 3 / (r * r * den)).ln()
        + a_const * c
        + s_log
        + c * (Decimal(2) / r).ln()
        - log_half
    )
    r2 = (
        (Decimal(4) * (r / eta) ** 3 / (r * r * c * c * den)).ln()
        - ((r / eta).ln() - a_const) * c
        - log_half
    )
    return {
        "gaussian_minus_required": gaussian,
        "log_R0_minus_log_half": r0,
        "log_R1_minus_log_half": r1,
        "log_R2_minus_log_half": r2,
    }


def trace_cutoff_derivative_bounds(c_value: int) -> dict[str, Decimal]:
    """Derivative signs proving the cutoff margins persist for larger C."""

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(TRACE_R)
    eta = decimal_fraction(TRACE_ETA)
    a_const = decimal_fraction(TRACE_A)
    l_value = decimal_fraction(TRACE_L)
    tau_derivative_upper = -eta * (eta * c - 1) / (
        Decimal(8) * (Decimal(1) + tau_rains_correction(c_value))
    )
    return {
        # The Mills derivative is at least this since L/(LC+1) < 1/C.
        "gaussian_derivative_lower": a_const - l_value / 2 - Decimal(1) / (2 * c),
        "log_R0_derivative": a_const + tau_derivative_upper,
        # The omitted -den'/den term only decreases this derivative.
        "log_R1_derivative_upper": (
            Decimal(2) / c
            + a_const
            + tau_derivative_upper
            + (Decimal(2) / r).ln()
        ),
        # The omitted -2/C-den'/den terms only decrease this derivative.
        "log_R2_derivative_upper": a_const - (r / eta).ln(),
    }


def trace_cutoff_certificate(cutoff: int = TRACE_MONOTONE_CUTOFF) -> None:
    """Print the concrete C_G >= cutoff trace-bridge certificate margins."""

    margins = trace_cutoff_logs(cutoff)
    derivatives = trace_cutoff_derivative_bounds(cutoff)
    print("trace cutoff certificate")
    print(f"cutoff C_G >= {cutoff}")
    print(f"r={TRACE_R} eta={TRACE_ETA} A={TRACE_A} L={TRACE_L}")
    for key, value in margins.items():
        relation = "OK" if value >= 0 and key == "gaussian_minus_required" else ""
        if key != "gaussian_minus_required":
            relation = "OK" if value <= 0 else "FAIL"
        print(f"{key} = {value} {relation}")
    for key, value in derivatives.items():
        if key == "gaussian_derivative_lower":
            relation = "OK" if value > 0 else "FAIL"
        else:
            relation = "OK" if value < 0 else "FAIL"
        print(f"{key} = {value} {relation}")
    if margins["gaussian_minus_required"] < 0:
        raise SystemExit("Gaussian margin is negative")
    for key in ("log_R0_minus_log_half", "log_R1_minus_log_half", "log_R2_minus_log_half"):
        if margins[key] > 0:
            raise SystemExit(f"{key} is positive")
    if derivatives["gaussian_derivative_lower"] <= 0:
        raise SystemExit("Gaussian derivative lower bound is not positive")
    for key in ("log_R0_derivative", "log_R1_derivative_upper", "log_R2_derivative_upper"):
        if derivatives[key] >= 0:
            raise SystemExit(f"{key} is not negative")


def decimal_one_trace_tv_log_bound(n_value: int) -> Decimal:
    """Decimal version of the CJL/Stirling log bound for E_1(n)."""

    n = Decimal(n_value)
    return Decimal(5).ln() + (Decimal(2) * n).ln().ln() / Decimal(4) - n * (
        (Decimal(2) * n).ln() - Decimal(1)
    )


def decimal_unitary_trace_tv_log_bound(n_value: int) -> Decimal:
    """Decimal CJL Theorem 1.1 upper bound for the complex unitary trace."""

    n = Decimal(n_value)
    return (
        Decimal(19).ln()
        + n.ln() / Decimal(4)
        + n.ln().ln() / Decimal(2)
        - Decimal(lgamma(n_value + 2))
    )


def decimal_logadd(first: Decimal, second: Decimal) -> Decimal:
    """Stable log(exp(first)+exp(second)) for Decimal values."""

    if second > first:
        first, second = second, first
    return first + (Decimal(1) + (second - first).exp()).ln()


def decimal_logsub(first: Decimal, second: Decimal) -> Decimal | None:
    """Stable log(exp(first)-exp(second)); return None if nonpositive."""

    if second >= first:
        return None
    return first + (Decimal(1) - (second - first).exp()).ln()


def square_tau_log(c_value: int, eta_value: Fraction = TRACE_SQUARE_ETA) -> Decimal:
    """Rains-square one-sided tau-tail log bound."""

    c = Decimal(c_value)
    eta = decimal_fraction(eta_value)
    n_half = c_value // 2
    gaussian_log = -((eta * c - Decimal(1)) ** 2) / Decimal(4)
    tv_log = Decimal(2).ln() + decimal_one_trace_tv_log_bound(n_half)
    return decimal_logadd(gaussian_log, tv_log)


def trace_square_cutoff_logs(c_value: int) -> dict[str, Decimal]:
    """High-precision log margins for the Rains-square trace cutoff."""

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(TRACE_SQUARE_R)
    eta = decimal_fraction(TRACE_SQUARE_ETA)
    a_const = decimal_fraction(TRACE_SQUARE_A)
    l_value = decimal_fraction(TRACE_SQUARE_L)
    log_half = Decimal("0.5").ln()
    den = Decimal(1) - (Decimal(1) / (eta * eta * c * c))
    s_log = square_tau_log(c_value)
    positive_tv_log = decimal_one_trace_tv_log_bound(c_value)

    gaussian = (
        (a_const - l_value / 2) * c
        + (l_value * c).ln() / 2
        - (l_value * c + 1).ln()
        - Decimal(2).ln() / 2
        - Decimal(2).ln()
    )
    r0 = a_const * c + s_log - log_half
    r1 = (
        (Decimal(4) * (Decimal(1) + c ** Decimal(-2)) * (r / 2) ** 3 / (r * r * den)).ln()
        + a_const * c
        + s_log
        + c * (Decimal(2) / r).ln()
        - log_half
    )
    r2 = (
        (Decimal(4) * (r / eta) ** 3 / (r * r * c * c * den)).ln()
        - ((r / eta).ln() - a_const) * c
        - log_half
    )
    unitary_vs_orth = (
        decimal_unitary_trace_tv_log_bound(c_value)
        - decimal_one_trace_tv_log_bound(c_value // 2)
    )
    return {
        "gaussian_minus_required": gaussian,
        "positive_tv_log_plus_A_C": positive_tv_log + a_const * c,
        "square_tau_log": s_log,
        "unitary_projection_log_minus_half_orth_log": unitary_vs_orth,
        "log_R0_minus_log_half": r0,
        "log_R1_minus_log_half": r1,
        "log_R2_minus_log_half": r2,
    }


def trace_square_derivative_bounds(c_value: int) -> dict[str, Decimal]:
    """Derivative and integer-step checks for the Rains-square cutoff."""

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(TRACE_SQUARE_R)
    eta = decimal_fraction(TRACE_SQUARE_ETA)
    a_const = decimal_fraction(TRACE_SQUARE_A)
    l_value = decimal_fraction(TRACE_SQUARE_L)
    return {
        "gaussian_derivative_lower": a_const - l_value / 2 - Decimal(1) / (2 * c),
        "positive_tv_derivative_upper": (
            a_const - (Decimal(2) * c).ln()
            + Decimal(1) / (Decimal(4) * c * (Decimal(2) * c).ln())
        ),
        "square_gaussian_R0_derivative": a_const
        - eta * (eta * c - Decimal(1)) / Decimal(2),
        "square_gaussian_R1_derivative_upper": (
            a_const
            - eta * (eta * c - Decimal(1)) / Decimal(2)
            + (Decimal(2) / r).ln()
        ),
        "log_R2_derivative_upper": a_const - (r / eta).ln(),
    }


def trace_square_cutoff_certificate(cutoff: int = TRACE_SQUARE_CUTOFF) -> None:
    """Print the Rains-square trace-bridge certificate margins."""

    margins = trace_square_cutoff_logs(cutoff)
    derivatives = trace_square_derivative_bounds(cutoff)
    print("Rains-square trace cutoff certificate")
    print(f"cutoff C_G >= {cutoff}")
    print(
        f"r={TRACE_SQUARE_R} eta={TRACE_SQUARE_ETA} "
        f"A={TRACE_SQUARE_A} L={TRACE_SQUARE_L}"
    )
    for key, value in margins.items():
        relation = "OK"
        if key == "gaussian_minus_required":
            relation = "OK" if value >= 0 else "FAIL"
        elif key == "square_tau_log":
            relation = ""
        else:
            relation = "OK" if value <= 0 else "FAIL"
        print(f"{key} = {value} {relation}")
    for key, value in derivatives.items():
        relation = "OK" if value < 0 else "FAIL"
        if key == "gaussian_derivative_lower":
            relation = "OK" if value > 0 else "FAIL"
        print(f"{key} = {value} {relation}")

    first_bad = None
    last_bad = None
    for c_value in range(248, cutoff):
        replay = trace_square_cutoff_logs(c_value)
        c_derivatives = trace_square_derivative_bounds(c_value)
        ok = (
            replay["gaussian_minus_required"] >= 0
            and replay["positive_tv_log_plus_A_C"] <= 0
            and replay["unitary_projection_log_minus_half_orth_log"] <= 0
            and replay["log_R0_minus_log_half"] <= 0
            and replay["log_R1_minus_log_half"] <= 0
            and replay["log_R2_minus_log_half"] <= 0
            and c_derivatives["gaussian_derivative_lower"] > 0
            and c_derivatives["log_R2_derivative_upper"] < 0
        )
        if not ok:
            if first_bad is None:
                first_bad = c_value
            last_bad = c_value
    if last_bad is not None:
        print(f"pre-cutoff bad range starts at {first_bad} and ends at {last_bad}")

    for c_value in range(cutoff, 10001):
        replay = trace_square_cutoff_logs(c_value)
        c_derivatives = trace_square_derivative_bounds(c_value)
        if (
            replay["gaussian_minus_required"] < 0
            or replay["positive_tv_log_plus_A_C"] > 0
            or replay["unitary_projection_log_minus_half_orth_log"] > 0
            or replay["log_R0_minus_log_half"] > 0
            or replay["log_R1_minus_log_half"] > 0
            or replay["log_R2_minus_log_half"] > 0
            or c_derivatives["gaussian_derivative_lower"] <= 0
            or c_derivatives["log_R2_derivative_upper"] >= 0
        ):
            raise SystemExit(f"Rains-square replay failed at C={c_value}")
    print(f"integer replay OK for {cutoff} <= C <= 10000")


def c_value(family: str, rank: int) -> int:
    return int(family_data(family, rank)["c"])


def stable_odd_reach(family: str, rank: int) -> int:
    """Largest odd exponent reached by Corollary 16 in this rank."""

    if family in {"B", "C"}:
        return 2 * rank - 1
    if family == "D":
        bound = rank - 3
        if bound < 1:
            return 1
        return bound if bound % 2 == 1 else bound - 1
    raise ValueError(f"unknown family {family!r}")


def rank_range_text(ranks: list[int]) -> str:
    if not ranks:
        return "none"
    ranges: list[tuple[int, int]] = []
    start = prev = ranks[0]
    for rank in ranks[1:]:
        if rank == prev + 1:
            prev = rank
        else:
            ranges.append((start, prev))
            start = prev = rank
    ranges.append((start, prev))
    return ", ".join(str(lo) if lo == hi else f"{lo}..{hi}" for lo, hi in ranges)


def parity_range_text(ranks: list[int]) -> str:
    if len(ranks) >= 3 and all(ranks[index + 1] == ranks[index] + 2 for index in range(len(ranks) - 1)):
        parity = "even" if ranks[0] % 2 == 0 else "odd"
        return f"{ranks[0]}..{ranks[-1]} {parity}"
    return rank_range_text(ranks)


def finite_residue_rows(cutoff: int = TRACE_CUTOFF) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    for family in ("B", "C"):
        for rank in range(2, cutoff):
            rows.append((family, rank))
    # D_b has C_G=b for even b and C_G=b-2 for odd b; for even cutoffs this
    # leaves one finite odd-rank row just above the cutoff.
    for rank in range(4, cutoff + 2):
        if c_value("D", rank) < cutoff:
            rows.append(("D", rank))
    return rows


def finite_residue_plan(cutoff: int = TRACE_CUTOFF, current_max_chain: int = 14) -> None:
    """Print the finite classical residue induced by the high-rank cutoff."""

    rows = finite_residue_rows(cutoff)
    current_odd = 2 * current_max_chain + 3
    max_stable = max(stable_odd_reach(family, rank) for family, rank in rows)
    print("finite classical residue plan")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(f"current exact bridge: m <= {current_max_chain}, odd n <= {current_odd}")
    print(f"finite rows total: {len(rows)}")
    for family in ("B", "C", "D"):
        family_rows = [rank for fam, rank in rows if fam == family]
        print(
            f"{family}: ranks {rank_range_text(family_rows)} "
            f"({len(family_rows)} rows)"
        )
    print(f"largest finite stable odd reach: {max_stable}")
    print("stable edge status:")
    print("  Corollary 16 already covers every finite row through its stable edge")
    print("  no exact bridge is needed merely to reach the largest stable edge")
    print("row-gated middle bridge target:")
    print("  for each finite row G, prove Chain only when")
    print("    stable_odd_reach(G) < 2m+3 < N_loc(G)")
    print("  Proposition 24 plus finite-leading dominance covers odd n >= N_loc(G)")
    print("  Corollary 16 covers odd n <= stable_odd_reach(G)")
    print("middle endpoint:")
    print("  M_mid = max { m : stable_odd_reach(G) < 2m+3 < N_loc(G)")
    print("                 for some finite row G with C_G < cutoff }")
    print("  N_loc(G) still requires calibrated Wong constants B_G and N_W(G)")


def middle_slice_rows(chain_m: int, cutoff: int = TRACE_CUTOFF) -> list[tuple[str, int]]:
    """Rows whose stable edge is before the target odd exponent 2m+3."""

    target_odd = 2 * chain_m + 3
    return [
        (family, rank)
        for family, rank in finite_residue_rows(cutoff)
        if stable_odd_reach(family, rank) < target_odd
    ]


def middle_slice_plan(
    cutoff: int = TRACE_CUTOFF,
    current_max_chain: int = 14,
    steps: int = 3,
) -> None:
    """Print the pre-Wong row gate for the next exact Chain steps."""

    print("row-gated middle slice plan")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(f"current exact bridge: m <= {current_max_chain}, odd n <= {2 * current_max_chain + 3}")
    print("rows listed still also require the inequality 2m+3 < N_loc(G)")
    for chain_m in range(current_max_chain + 1, current_max_chain + steps + 1):
        target_odd = 2 * chain_m + 3
        rows = middle_slice_rows(chain_m, cutoff)
        print(f"m={chain_m}, odd n={target_odd}: pre-Wong rows {len(rows)}")
        for family in ("B", "C", "D"):
            family_rows = [rank for fam, rank in rows if fam == family]
            print(
                f"  {family}: ranks {rank_range_text(family_rows)} "
                f"({len(family_rows)} rows)"
            )


def first_hit_onset_plan(cutoff: int = TRACE_CUTOFF, current_max_chain: int = 20) -> None:
    """Print the row-wise local-onset target equivalent to no later middle steps."""

    first_unbridged_odd = 2 * (current_max_chain + 1) + 3
    rows = finite_residue_rows(cutoff)
    print("first-hit local-onset plan")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(f"current exact bridge: m <= {current_max_chain}")
    print(f"first exact-unbridged odd n: {first_unbridged_odd}")
    print("row target: N_loc(G) <= max(first exact-unbridged odd, stable_odd_reach(G)+2)")
    for family in ("B", "C"):
        ranks = [rank for fam, rank in rows if fam == family]
        flat = [rank for rank in ranks if stable_odd_reach(family, rank) + 2 <= first_unbridged_odd]
        sloped = [rank for rank in ranks if stable_odd_reach(family, rank) + 2 > first_unbridged_odd]
        print(f"{family}: {len(ranks)} rows")
        print(f"  N_loc <= {first_unbridged_odd}: ranks {rank_range_text(flat)}")
        if sloped:
            print(
                "  N_loc <= 2b+1: "
                f"ranks {rank_range_text(sloped)} "
                f"(targets {2 * sloped[0] + 1}..{2 * sloped[-1] + 1})"
            )
    d_ranks = [rank for fam, rank in rows if fam == "D"]
    d_flat = [rank for rank in d_ranks if stable_odd_reach("D", rank) + 2 <= first_unbridged_odd]
    d_even = [
        rank
        for rank in d_ranks
        if rank % 2 == 0 and stable_odd_reach("D", rank) + 2 > first_unbridged_odd
    ]
    d_odd = [
        rank
        for rank in d_ranks
        if rank % 2 == 1 and stable_odd_reach("D", rank) + 2 > first_unbridged_odd
    ]
    print(f"D: {len(d_ranks)} rows")
    print(f"  N_loc <= {first_unbridged_odd}: ranks {rank_range_text(d_flat)}")
    if d_even:
        print(
            "  N_loc <= b-1: "
            f"even ranks {parity_range_text(d_even)} "
            f"(targets {d_even[0] - 1}..{d_even[-1] - 1})"
        )
    if d_odd:
        odd_targets = [rank - 2 for rank in d_odd]
        print(
            "  N_loc <= b-2: "
            f"odd ranks {parity_range_text(d_odd)} "
            f"(targets {parity_range_text(odd_targets)})"
        )


def finite_trace_pushforward_margin(
    c_value: int,
    odd_n: int,
    r_value: Fraction = TRACE_SQUARE_R,
    eta_value: Fraction = TRACE_SQUARE_ETA,
) -> tuple[str, Decimal | None]:
    """Test the source-audited trace-pushforward inequality at one C and odd n."""

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(r_value)
    eta = decimal_fraction(eta_value)
    l_value = Decimal(2) * r + Decimal(3) * eta
    if c_value < 132:
        return "CJL-source-range", None
    if eta * c <= 1:
        return "a<=1", None
    den = Decimal(1) - Decimal(1) / (eta * eta * c * c)
    if den <= 0:
        return "chebyshev-denominator", None

    x_square = l_value * c
    gaussian_log = (
        -x_square / 2
        + x_square.ln() / 2
        - (x_square + 1).ln()
        - Decimal(2).ln() / 2
    )
    positive_tv_log = decimal_one_trace_tv_log_bound(c_value)
    p1_log = decimal_logsub(gaussian_log, positive_tv_log)
    if p1_log is None:
        return "positive-tail-lower-nonpositive", positive_tv_log - gaussian_log

    p2_log = square_tau_log(c_value, eta_value)
    positive_gap_log = decimal_logsub(p1_log, p2_log)
    if positive_gap_log is None:
        return "tau-upper-beats-positive-tail", p2_log - p1_log

    left_log = (r * r * c * c * den).ln() + positive_gap_log
    term1_log = (
        (Decimal(2) * (c * c + Decimal(1))).ln()
        + p2_log
        + Decimal(odd_n) * (Decimal(2) / r).ln()
    )
    term2_log = Decimal(2).ln() + Decimal(odd_n) * (eta / r).ln()
    right_log = decimal_logadd(term1_log, term2_log)
    return ("OK", left_log - right_log) if left_log >= right_log else (
        "pushforward-ineq-negative",
        right_log - left_log,
    )


def cjl_chernoff_tau_log(c_value: int, eta_value: Fraction) -> Decimal:
    """CJL/Rains Chernoff tau-tail log bound without the TV approximation term."""

    c = Decimal(c_value)
    eta = decimal_fraction(eta_value)
    return Decimal(2).ln() - (eta * c - Decimal(1)) ** 2 / (
        Decimal(16) * (Decimal(1) + tau_rains_correction(c_value))
    )


def cjl_chernoff_admissible(c_value: int, eta_value: Fraction) -> bool:
    """Check the CJL power-trace Chernoff parameter window for m=2."""

    c = Decimal(c_value)
    eta = decimal_fraction(eta_value)
    n_half = max(2, c_value // 2)
    threshold = eta * c - Decimal(1)
    window = (
        Decimal(4)
        * (-(Decimal(3) / Decimal(2))).exp()
        * (Decimal(1) + tau_rains_correction(c_value))
        * Decimal(n_half)
    )
    return threshold > 0 and threshold < window


def finite_trace_pushforward_margin_chernoff_tau(
    c_value: int,
    odd_n: int,
    r_value: Fraction,
    eta_value: Fraction,
) -> tuple[str, Decimal | None]:
    """Finite-row trace-pushforward test using the CJL Chernoff tau tail."""

    getcontext().prec = 80
    c = Decimal(c_value)
    r = decimal_fraction(r_value)
    eta = decimal_fraction(eta_value)
    l_value = Decimal(2) * r + Decimal(3) * eta
    if c_value < 132:
        return "CJL-one-trace-source-range", None
    if eta * c <= 1:
        return "a<=1", None
    if not cjl_chernoff_admissible(c_value, eta_value):
        return "CJL-chernoff-window", None
    den = Decimal(1) - Decimal(1) / (eta * eta * c * c)
    if den <= 0:
        return "chebyshev-denominator", None

    x_square = l_value * c
    gaussian_log = (
        -x_square / 2
        + x_square.ln() / 2
        - (x_square + 1).ln()
        - Decimal(2).ln() / 2
    )
    positive_tv_log = decimal_one_trace_tv_log_bound(c_value)
    p1_log = decimal_logsub(gaussian_log, positive_tv_log)
    if p1_log is None:
        return "positive-tail-lower-nonpositive", positive_tv_log - gaussian_log

    p2_log = cjl_chernoff_tau_log(c_value, eta_value)
    positive_gap_log = decimal_logsub(p1_log, p2_log)
    if positive_gap_log is None:
        return "tau-upper-beats-positive-tail", p2_log - p1_log

    left_log = (r * r * c * c * den).ln() + positive_gap_log
    term1_log = (
        (Decimal(2) * (c * c + Decimal(1))).ln()
        + p2_log
        + Decimal(odd_n) * (Decimal(2) / r).ln()
    )
    term2_log = Decimal(2).ln() + Decimal(odd_n) * (eta / r).ln()
    right_log = decimal_logadd(term1_log, term2_log)
    return ("OK", left_log - right_log) if left_log >= right_log else (
        "pushforward-ineq-negative",
        right_log - left_log,
    )


def post_m29_direct_closed_rows() -> set[tuple[str, int]]:
    """Rows already closed by direct post-m=29 tail propositions."""

    return {
        ("B", 2),
        ("C", 2),
        ("B", 3),
        ("C", 3),
        ("D", 4),
        ("B", 4),
        ("C", 4),
        ("D", 5),
        ("C", 5),
        ("B", 5),
        ("C", 6),
        ("D", 6),
        ("B", 6),
        ("D", 7),
        ("B", 7),
        ("C", 7),
        ("C", 8),
        ("D", 8),
        ("B", 8),
        ("C", 9),
        ("D", 9),
        ("C", 10),
        ("B", 9),
        ("C", 11),
        ("C", 12),
        ("B", 10),
        ("D", 10),
        ("D", 11),
    }


def first_hit_trace_closed_rows() -> set[tuple[str, int]]:
    """Rows closed by Proposition 24.17D."""

    return {
        *(("B", rank) for rank in [276, *range(278, 296)]),
        *(("C", rank) for rank in [276, *range(278, 296)]),
        *(("D", rank) for rank in [276, 278, *range(280, 296), 297]),
    }


def post_m29_open_rows(cutoff: int = TRACE_CUTOFF) -> list[tuple[str, int]]:
    closed = post_m29_direct_closed_rows() | first_hit_trace_closed_rows()
    return [row for row in finite_residue_rows(cutoff) if row not in closed]


def post_m29_all_later_direct_closed_rows() -> set[tuple[str, int]]:
    """Rows removed from the all-later assumption by direct tail propositions."""

    return post_m29_direct_closed_rows() | {
        ("B", 11),
        ("B", 12),
        ("B", 13),
        ("C", 13),
        ("C", 14),
        ("C", 15),
        ("C", 16),
        ("C", 17),
        ("C", 18),
        ("C", 19),
        ("D", 12),
        ("D", 13),
        ("D", 14),
        ("D", 15),
        ("D", 16),
        ("D", 17),
        ("D", 18),
        ("D", 19),
        ("D", 20),
        ("D", 21),
        ("D", 22),
        ("D", 23),
        ("D", 24),
        ("D", 25),
        ("D", 26),
        ("D", 27),
        ("D", 28),
        ("D", 29),
        ("D", 30),
        ("D", 31),
        ("D", 32),
        ("D", 33),
        ("D", 34),
        ("D", 35),
        ("D", 36),
        ("D", 37),
        ("D", 38),
        ("D", 39),
        ("D", 40),
        ("D", 41),
        ("D", 42),
        ("D", 43),
        ("D", 44),
        ("D", 45),
        ("D", 46),
        ("D", 47),
        ("D", 48),
        ("D", 49),
        ("D", 50),
        ("D", 51),
        ("D", 52),
        ("D", 53),
        ("D", 54),
        ("D", 55),
        ("D", 56),
        ("D", 57),
        ("D", 58),
        ("D", 59),
        ("D", 60),
        ("D", 61),
        ("D", 62),
        ("D", 63),
        ("D", 64),
        ("D", 65),
        ("D", 66),
        ("D", 67),
        ("D", 68),
        ("D", 69),
        ("D", 70),
        ("D", 71),
        ("D", 72),
        ("D", 73),
        ("D", 74),
        ("D", 75),
        ("D", 76),
        ("D", 77),
        ("D", 78),
        ("D", 79),
        ("D", 80),
        ("D", 81),
        ("D", 82),
        ("D", 83),
        ("D", 84),
        ("D", 85),
        ("D", 86),
        ("D", 87),
        ("D", 88),
        ("D", 89),
        ("D", 90),
        ("D", 91),
        ("D", 92),
        ("D", 93),
        ("D", 94),
        ("D", 95),
        ("D", 96),
        ("D", 97),
    }


def post_m29_all_later_open_rows(cutoff: int = TRACE_CUTOFF) -> list[tuple[str, int]]:
    """Rows still relevant to the post-m=29 all-later assumption before trace."""

    closed = post_m29_all_later_direct_closed_rows() | first_hit_trace_closed_rows()
    return [row for row in finite_residue_rows(cutoff) if row not in closed]


def post_m29_chernoff_trace_certificate(cutoff: int = TRACE_CUTOFF) -> None:
    """Diagnose the source-supported CJL/Rains Chernoff tau supplier."""

    params = {
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
    expected = post_m29_chernoff_trace_closed_rows()
    rows = post_m29_all_later_open_rows(cutoff)
    ok_rows: list[tuple[str, int, str, Fraction, Fraction, Decimal]] = []
    reason_counts: Counter[str] = Counter()
    for family, rank in rows:
        target_odd = max(63, stable_odd_reach(family, rank) + 2)
        failures: list[str] = []
        for label, r_value, eta_value in params[family]:
            reason, margin = finite_trace_pushforward_margin_chernoff_tau(
                c_value(family, rank),
                target_odd,
                r_value,
                eta_value,
            )
            if reason == "OK":
                assert margin is not None
                ok_rows.append((family, rank, label, r_value, eta_value, margin))
                break
            failures.append(reason)
        else:
            reason_counts[failures[-1]] += 1

    ok_set = {(family, rank) for family, rank, *_rest in ok_rows}
    if ok_set != expected:
        missing = sorted(expected - ok_set)
        unexpected = sorted(ok_set - expected)
        raise SystemExit(
            "post-m29 Chernoff trace row mismatch: "
            f"missing={missing[:10]} unexpected={unexpected[:10]}"
        )
    if any(row[-1] <= 0 for row in ok_rows):
        raise SystemExit("post-m29 Chernoff trace margin is not positive")

    print("post-m29 Chernoff trace certificate")
    print(f"  open rows tested: {len(rows)}")
    print(f"  rows closed: {len(ok_rows)}")
    if not ok_rows:
        print("  source-supported no-TV Chernoff supplier closes no rows")
        if reason_counts:
            print(f"  first failure counts: {dict(reason_counts)}")
        return
    for family in "BCD":
        family_rows = [
            (rank, label, r_value, eta_value, margin)
            for fam, rank, label, r_value, eta_value, margin in ok_rows
            if fam == family
        ]
        if not family_rows:
            continue
        ranks = [rank for rank, *_rest in family_rows]
        print(f"  {family}: ranks {rank_range_text(ranks)} ({len(ranks)} rows)")
        for label in sorted({label for _rank, label, *_rest in family_rows}):
            label_rows = [row for row in family_rows if row[1] == label]
            label_ranks = [rank for rank, *_rest in label_rows]
            min_row = min(label_rows, key=lambda row: row[4])
            _rank, _label, r_value, eta_value, margin = min_row
            print(
                f"    {label}: r={r_value}, eta={eta_value}, "
                f"ranks {rank_range_text(label_ranks)} ({len(label_ranks)} rows), "
                f"minimum log margin {margin} at {family}_{min_row[0]}"
            )
    if reason_counts:
        print("  non-closed reason counts:")
        for reason, count in sorted(reason_counts.items()):
            print(f"    {reason}: {count}")


def finite_trace_tail_diagnostic(
    cutoff: int = TRACE_CUTOFF,
    current_max_chain: int = 19,
    steps: int = 3,
) -> None:
    """Print the finite-row coverage of the Rains-square trace route."""

    print("finite trace-tail diagnostic")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(
        f"constants r={TRACE_SQUARE_R} eta={TRACE_SQUARE_ETA} "
        f"L={TRACE_SQUARE_L}"
    )
    print("criterion: Proposition 24.8 + Proposition 24.10 with CJL/Rains-square tails")
    for chain_m in range(current_max_chain + 1, current_max_chain + steps + 1):
        target_odd = 2 * chain_m + 3
        rows = middle_slice_rows(chain_m, cutoff)
        reasons: Counter[str] = Counter()
        ok_rows: list[tuple[str, int, int, Decimal | None]] = []
        max_obstruction: tuple[Decimal, str, int, int, str] | None = None
        for family, rank in rows:
            c = c_value(family, rank)
            reason, margin = finite_trace_pushforward_margin(c, target_odd)
            reasons[reason] += 1
            if reason == "OK":
                ok_rows.append((family, rank, c, margin))
            elif margin is not None:
                row = (margin, family, rank, c, reason)
                if max_obstruction is None or row > max_obstruction:
                    max_obstruction = row
        print(f"m={chain_m}, odd n={target_odd}: rows {len(rows)}, trace-covered {len(ok_rows)}")
        for reason, count in sorted(reasons.items()):
            print(f"  {reason}: {count}")
        if ok_rows:
            first = ok_rows[0]
            last = ok_rows[-1]
            print(
                "  covered row range sample: "
                f"{first[0]}_{first[1]}(C={first[2]}) .. "
                f"{last[0]}_{last[1]}(C={last[2]})"
            )
        if max_obstruction is not None:
            margin, family, rank, c, reason = max_obstruction
            print(
                "  largest displayed obstruction: "
                f"{family}_{rank} C={c} {reason} margin={margin}"
            )


def first_hit_trace_diagnostic(
    cutoff: int = TRACE_CUTOFF,
    current_max_chain: int = 20,
    r_value: Fraction = FIRST_HIT_TRACE_R,
    eta_value: Fraction = FIRST_HIT_TRACE_ETA,
) -> None:
    """Test trace-pushforward closure at the first-hit exponent of each row."""

    rows = finite_residue_rows(cutoff)
    reasons: Counter[str] = Counter()
    ok_rows: list[tuple[str, int, int, int, Decimal]] = []
    min_ok: tuple[str, int, int, int, Decimal] | None = None
    max_obstruction: tuple[Decimal, str, int, int, int, str] | None = None
    first_unbridged_odd = 2 * (current_max_chain + 1) + 3
    for family, rank in rows:
        target_odd = max(first_unbridged_odd, stable_odd_reach(family, rank) + 2)
        c = c_value(family, rank)
        reason, margin = finite_trace_pushforward_margin(c, target_odd, r_value, eta_value)
        reasons[reason] += 1
        if reason == "OK":
            assert margin is not None
            row = (family, rank, c, target_odd, margin)
            ok_rows.append(row)
            if min_ok is None or margin < min_ok[-1]:
                min_ok = row
        elif margin is not None:
            row = (margin, family, rank, c, target_odd, reason)
            if max_obstruction is None or row > max_obstruction:
                max_obstruction = row

    print("first-hit trace diagnostic")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(f"current exact bridge: m <= {current_max_chain}")
    print(f"constants r={r_value} eta={eta_value}")
    print(f"rows tested: {len(rows)}, trace-covered: {len(ok_rows)}")
    for reason, count in sorted(reasons.items()):
        print(f"  {reason}: {count}")
    for family in ("B", "C", "D"):
        family_ranks = [rank for fam, rank, _, _, _ in ok_rows if fam == family]
        print(f"  {family}: ranks {rank_range_text(family_ranks)} ({len(family_ranks)} rows)")
    if min_ok is not None:
        family, rank, c, target_odd, margin = min_ok
        print(
            "  minimum OK margin: "
            f"{family}_{rank} C={c} N_hit={target_odd} margin={margin}"
        )
    if max_obstruction is not None:
        margin, family, rank, c, target_odd, reason = max_obstruction
        print(
            "  largest displayed obstruction: "
            f"{family}_{rank} C={c} N_hit={target_odd} {reason} margin={margin}"
        )


def read_stable_moments(max_index: int) -> list[int]:
    moments: dict[int, int] = {}
    with STABLE_MOMENTS_PATH.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            index_text, value_text = line.split()[:2]
            index = int(index_text)
            if index <= max_index:
                moments[index] = int(value_text)
    missing = [index for index in range(max_index + 1) if index not in moments]
    if missing:
        raise RuntimeError(f"missing stable moments: {missing}")
    return [moments[index] for index in range(max_index + 1)]


def stable_moments_formula(
    max_index: int,
    progress_label: str | None = None,
) -> list[int]:
    """Stable SO(infinity) moments from the closed generating function."""

    progress_step = max(1, max_index // 20)

    def progress(stage: str, index: int) -> None:
        if progress_label is None:
            return
        if index == 0 or index == max_index or index % progress_step == 0:
            print(
                f"{progress_label} {stage} {index}/{max_index}",
                flush=True,
            )

    moments = [0 for _ in range(max_index + 1)]
    if max_index >= 0:
        moments[0] = 1
        progress("moment", 0)
    if max_index >= 1:
        moments[1] = 0
        progress("moment", 1)
    for index in range(1, max_index):
        value = index * moments[index] + index * moments[index - 1]
        if index >= 2:
            value -= index * (index - 1) * moments[index - 2] // 2
        moments[index + 1] = value
        progress("moment", index + 1)
    if max_index >= 0 and moments[0] != 1:
        raise RuntimeError("stable moment formula produced an inconsistent m_0")
    if max_index >= 1 and moments[1] != 0:
        raise RuntimeError("stable moment formula produced an inconsistent m_1")
    return moments


def binom_integer(n_value: int, k_value: int) -> int:
    if k_value < 0 or k_value > n_value:
        return 0
    k = min(k_value, n_value - k_value)
    out = 1
    for j in range(1, k + 1):
        out = out * (n_value - k + j) // j
    return out


def binom_row(n_value: int) -> list[int]:
    row = [0 for _ in range(n_value + 1)]
    value = 1
    for k_value in range(n_value + 1):
        row[k_value] = value
        if k_value < n_value:
            value = value * (n_value - k_value) // (k_value + 1)
    return row


def integer_partitions(total: int, max_part: int | None = None):
    if total == 0:
        yield ()
        return
    if max_part is None or max_part > total:
        max_part = total
    for first in range(max_part, 0, -1):
        remaining = total - first
        next_max = min(first, remaining) if remaining else 0
        for tail in integer_partitions(remaining, next_max):
            yield (first,) + tail


def conjugate_partition(partition: tuple[int, ...]) -> tuple[int, ...]:
    if not partition:
        return ()
    return tuple(
        sum(1 for row in partition if row >= column)
        for column in range(1, partition[0] + 1)
    )


def is_vertical_two_strip(lower: tuple[int, ...], upper: tuple[int, ...]) -> bool:
    if sum(upper) - sum(lower) != 2:
        return False
    length = max(len(lower), len(upper))
    lower_rows = list(lower) + [0] * (length - len(lower))
    upper_rows = list(upper) + [0] * (length - len(upper))
    return all(
        upper_rows[index] >= lower_rows[index]
        and upper_rows[index] - lower_rows[index] in (0, 1)
        for index in range(length)
    )


@lru_cache(maxsize=None)
def pieri_e2_predecessors(partition: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    out: set[tuple[int, ...]] = set()
    for first, second in combinations(range(len(partition)), 2):
        predecessor = list(partition)
        predecessor[first] -= 1
        predecessor[second] -= 1
        predecessor_tuple = tuple(
            sorted((row for row in predecessor if row > 0), reverse=True)
        )
        if is_vertical_two_strip(predecessor_tuple, partition):
            out.add(predecessor_tuple)
    return tuple(out)


@lru_cache(maxsize=None)
def pieri_e2_coefficient(power: int, partition: tuple[int, ...]) -> int:
    if power == 0:
        return 1 if partition == () else 0
    if sum(partition) != 2 * power:
        return 0
    return sum(
        pieri_e2_coefficient(power - 1, predecessor)
        for predecessor in pieri_e2_predecessors(partition)
    )


def d_determinant_correction_pieri(rank: int, depth: int) -> int:
    if depth < 0:
        return 0
    matrix_size = 2 * rank
    power = rank + depth
    total = 0
    for excess_partition in integer_partitions(depth):
        target = tuple(
            [1 + 2 * row for row in excess_partition]
            + [1] * (matrix_size - len(excess_partition))
        )
        total += pieri_e2_coefficient(power, target)
    return total


def q3_integer(moments: list[int], n_value: int) -> int:
    row = binom_row(n_value)
    return 2 * sum(
        row[k]
        * (
            moments[k + 2] * moments[n_value - k]
            - moments[k + 1] * moments[n_value - k + 1]
        )
        for k in range(n_value + 1)
    )


def chain_diff_integer(moments: list[int], chain_m: int) -> int:
    return q3_integer(moments, 2 * chain_m + 3) - 4 * q3_integer(
        moments, 2 * chain_m + 1
    )


def chain_diff_linear_coefficients(moments: list[int], chain_m: int) -> list[int]:
    """Linear coefficients of D(chain_m) at the supplied moment vector."""

    max_index = 2 * chain_m + 5
    coeffs = [0] * (max_index + 1)
    for n_value, scale in ((2 * chain_m + 3, 1), (2 * chain_m + 1, -4)):
        row = binom_row(n_value)
        for k in range(n_value + 1):
            c = 2 * scale * row[k]
            first = k + 2
            second = n_value - k
            coeffs[first] += c * moments[second]
            coeffs[second] += c * moments[first]
            first = k + 1
            second = n_value - k + 1
            coeffs[first] -= c * moments[second]
            coeffs[second] -= c * moments[first]
    return coeffs


def chain_diff_quadratic_coefficients(chain_m: int) -> dict[tuple[int, int], int]:
    """Quadratic coefficients of D(chain_m) in commutative moment variables."""

    coeffs: dict[tuple[int, int], int] = {}
    for n_value, scale in ((2 * chain_m + 3, 1), (2 * chain_m + 1, -4)):
        row = binom_row(n_value)
        for k in range(n_value + 1):
            c = 2 * scale * row[k]
            for first, second, sign in (
                (k + 2, n_value - k, 1),
                (k + 1, n_value - k + 1, -1),
            ):
                key = tuple(sorted((first, second)))
                coeffs[key] = coeffs.get(key, 0) + sign * c
    return coeffs


BC16_REUSED_DELTAS: dict[tuple[str, int], dict[int, int]] = {
    ("B", 16): {
        34: -29671013856627,
        35: -20107792882889763,
        36: -7463034910126645010,
        37: -2010179347739912025378,
        38: -439778576182482066696518,
        39: -82991190661437015342653146,
        40: -14018475111566428408400478360,
        41: -2172791683567379592235980926420,
    },
    ("C", 16): {
        34: -6332659870762850625,
        35: -3546289527627196350001,
        36: -1085386238549398782872470,
        37: -240906540190770702448574886,
        38: -43440010478315010846896357231,
        39: -6762758539075699323052519323164,
        40: -943713224488106537386181639905578,
        41: -121053276068452153496514002390431460,
    },
}


B10_POST_M29_DELTAS: dict[int, int] = {
    22: -105089229,
    23: -27976738713,
    24: -4295265140026,
    25: -499323359850252,
    26: -48882284131436700,
    27: -4259955534370509792,
    28: -341710972133039521632,
    29: -25789649147329227629322,
    30: -1859519972814360249527589,
    31: -129522123451014892429010625,
    32: -8787934171910510215181275050,
    33: -584527223555429163503334344880,
    34: -38306540163725854247953705707090,
    35: -2483268800895316575995024378120850,
    36: -159755337735954867969182421522566380,
    37: -10226075182110393640826712744330518940,
    38: -652716453810005970170406516277872306766,
    39: -41618133724910207907296882492819875030918,
    40: -2654822364772729952773467193716489637424376,
    41: -169641269373625807217746830807893717898801320,
    42: -10870107803087217264484081401690116471482778420,
    43: -699093301667618363700703023522022909115504452380,
}


def parse_bc_delta_logs(paths: list[Path]) -> dict[tuple[str, int], dict[int, int]]:
    """Parse exact B/C boundary deltas printed by the C++ certificate."""

    pattern = re.compile(r"^\s*([BC])_(\d+) Delta_(\d+) = (-?\d+)")
    deltas: dict[tuple[str, int], dict[int, int]] = {}
    for path in paths:
        with path.open() as handle:
            for line in handle:
                match = pattern.match(line)
                if not match:
                    continue
                family = match.group(1)
                rank = int(match.group(2))
                moment = int(match.group(3))
                value = int(match.group(4))
                deltas.setdefault((family, rank), {})[moment] = value
    return deltas


def parse_d_delta_logs(paths: list[Path]) -> dict[int, dict[int, int]]:
    """Parse exact D determinant deltas printed by the C++ certificate."""

    pattern = re.compile(r"^\s*D_(\d+) Delta_(\d+) = (-?\d+)")
    deltas: dict[int, dict[int, int]] = {}
    for path in paths:
        with path.open() as handle:
            for line in handle:
                match = pattern.match(line)
                if not match:
                    continue
                rank = int(match.group(1))
                moment = int(match.group(2))
                value = int(match.group(3))
                if value < 0:
                    raise SystemExit(f"D_{rank} Delta_{moment} is negative")
                deltas.setdefault(rank, {})[moment] = value
    return deltas


def d_reused_delta_certificate(
    delta_logs: list[Path], chain_m: int, low_rank: int, high_rank: int
) -> None:
    """Check D rows from previously certified determinant deltas."""

    if low_rank < 4 or high_rank < low_rank:
        raise SystemExit("invalid D reused-delta rank range")
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
    negative_linear = [
        (index, value)
        for index, value in enumerate(linear_coeffs)
        if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    deltas = parse_d_delta_logs(delta_logs)

    print(f"m={chain_m} reused D determinant-delta certificate")
    print(f"stable Chain diff D({chain_m}) = {stable_diff}")
    print(f"certified rows: D_{low_rank}..D_{high_rank}")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")

    margins: list[tuple[int, int]] = []
    for rank in range(high_rank, low_rank - 1, -1):
        active_indices: set[int] = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        row_deltas = deltas.get(rank, {})
        missing = [
            index for index in sorted(active_indices)
            if index not in row_deltas
        ]
        if missing:
            raise SystemExit(f"D_{rank} missing active deltas: {missing}")
        negative_linear_sum = sum(
            value * row_deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_sum = sum(
            value * row_deltas[first] * row_deltas[second]
            for first, second, value in negative_quadratic
            if first >= rank and second >= rank
        )
        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}"
        )
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")
        margins.append((lower_margin, rank))

    min_margin, min_rank = min(margins)
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")
    print(f"row_count {len(margins)}")


def post_m29_chernoff_trace_closed_rows() -> set[tuple[str, int]]:
    """Rows closed by the source-supported no-TV CJL Chernoff supplier."""

    return {
        *(("B", rank) for rank in [218, *range(219, 276), 277]),
        *(("C", rank) for rank in [218, *range(219, 276), 277]),
    }


def post_m29_d_interval_closed_rows() -> set[tuple[str, int]]:
    """Rows closed by the D98-D279 interval direct-tail/frontier theorem."""

    return {*(("D", rank) for rank in [*range(98, 276), 277, 279])}


def post_m29_remaining_residue(cutoff: int = TRACE_CUTOFF) -> list[tuple[str, int]]:
    """The post-m=29 residue after direct, Chernoff, D-interval, and trace rows."""

    closed = post_m29_chernoff_trace_closed_rows() | post_m29_d_interval_closed_rows()
    return [
        row for row in post_m29_all_later_open_rows(cutoff)
        if row not in closed
    ]


def post_m29_lower_bound_supplier_certificate(
    bc_delta_logs: list[Path],
    d_delta_logs: list[Path],
    cutoff: int = TRACE_CUTOFF,
    *,
    bc_only: bool = False,
) -> None:
    """Check the first-hit Chain lower bound for the post-m=29 residue."""

    rows = post_m29_open_rows(cutoff)
    expected_bc = {
        *(("B", rank) for rank in [*range(11, 276), 277]),
        *(("C", rank) for rank in [*range(13, 276), 277]),
    }
    expected = expected_bc | {
        *(("D", rank) for rank in [*range(12, 276), 277, 279]),
    }
    if bc_only:
        rows = [row for row in rows if row[0] in {"B", "C"}]
        expected = expected_bc
    row_set = set(rows)
    if row_set != expected:
        missing = sorted(expected - row_set)
        unexpected = sorted(row_set - expected)
        raise SystemExit(
            "post-m29 lower-bound residue mismatch: "
            f"missing={missing[:10]} unexpected={unexpected[:10]}"
        )

    max_moment = max(
        2 * ((max(63, stable_odd_reach(family, rank) + 2) - 3) // 2) + 5
        for family, rank in rows
    )
    stable = stable_moments_formula(max_moment)
    table_check = read_stable_moments(min(max_moment, 100))
    if stable[: len(table_check)] != table_check:
        raise SystemExit("stable moment formula disagrees with the checked table")

    bc_deltas = parse_bc_delta_logs(bc_delta_logs)
    d_deltas = {} if bc_only else parse_d_delta_logs(d_delta_logs)
    exact_bc_rows = {
        *(("B", rank) for rank in range(11, 19)),
        *(("C", rank) for rank in range(13, 19)),
    }
    exact_d_rows = {("D", rank) for rank in range(12, 37)}
    category_margins: dict[str, list[tuple[int, str, int]]] = defaultdict(list)

    def row_chain_data(family: str, rank: int):
        target_odd = max(63, stable_odd_reach(family, rank) + 2)
        chain_m = (target_odd - 3) // 2
        moment_max = 2 * chain_m + 5
        local_stable = stable[: moment_max + 1]
        stable_diff = chain_diff_integer(local_stable, chain_m)
        if stable_diff <= 0:
            raise SystemExit(f"{family}_{rank} stable Chain diff is not positive")
        linear_coeffs = chain_diff_linear_coefficients(local_stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        return target_odd, chain_m, moment_max, stable_diff, linear_coeffs, quadratic_coeffs

    def check_bc_row(family: str, rank: int) -> None:
        (
            _target_odd,
            chain_m,
            moment_max,
            stable_diff,
            linear_coeffs,
            quadratic_coeffs,
        ) = row_chain_data(family, rank)
        first_boundary = 2 * rank + 2
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        if (family, rank) in exact_bc_rows:
            exact_through = 39
            row_deltas = {
                moment: value
                for moment, value in bc_deltas.get((family, rank), {}).items()
                if moment <= exact_through
            }
            missing = [
                moment
                for moment in range(first_boundary, exact_through + 1)
                if moment not in row_deltas
            ]
            if missing:
                raise SystemExit(
                    f"{family}_{rank} missing B/C deltas through m_{exact_through}: "
                    f"{missing}"
                )
            if any(value > 0 for value in row_deltas.values()):
                raise SystemExit(f"{family}_{rank} has a positive B/C correction")
            outside = [
                (first, second)
                for first, second, _value in negative_quadratic
                if first > exact_through or second > exact_through
            ]
            if outside:
                raise SystemExit(
                    f"{family}_{rank} has negative quadratic terms outside "
                    f"m_{exact_through}: {outside}"
                )
            tail = -sum(
                linear_coeffs[index] * stable[index]
                for index in range(exact_through + 1, moment_max + 1)
                if linear_coeffs[index] > 0
            )
            linear_exact = sum(
                linear_coeffs[index] * row_deltas[index]
                for index in row_deltas
            )
            negative_quadratic_exact = sum(
                value * row_deltas[first] * row_deltas[second]
                for first, second, value in negative_quadratic
            )
            positive_quadratic_exact = sum(
                value * row_deltas[first] * row_deltas[second]
                for (first, second), value in quadratic_coeffs.items()
                if value > 0 and first in row_deltas and second in row_deltas
            )
            lower_margin = (
                stable_diff
                + tail
                + linear_exact
                + negative_quadratic_exact
                + positive_quadratic_exact
            )
            category = "BC-exact-m30-window"
        else:
            tail = -sum(
                linear_coeffs[index] * stable[index]
                for index in range(first_boundary, moment_max + 1)
                if linear_coeffs[index] > 0
            )
            negative_quadratic_tail = sum(
                value * stable[first] * stable[second]
                for first, second, value in negative_quadratic
            )
            lower_margin = stable_diff + tail + negative_quadratic_tail
            category = (
                "BC-stable-tail-m30"
                if chain_m == 30
                else "BC-stable-tail-sloped"
            )

        if lower_margin <= 0:
            raise SystemExit(f"{family}_{rank} lower margin is not positive")
        category_margins[category].append((lower_margin, family, rank))

    def check_d_row(rank: int) -> None:
        (
            _target_odd,
            chain_m,
            _moment_max,
            stable_diff,
            linear_coeffs,
            quadratic_coeffs,
        ) = row_chain_data("D", rank)
        active_indices: set[int] = {
            index for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= rank
        }
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= rank and second >= rank
        ]
        for first, second, _value in negative_quadratic:
            active_indices.add(first)
            active_indices.add(second)

        if ("D", rank) in exact_d_rows:
            row_deltas = d_deltas.get(rank, {})
            missing = [
                index for index in sorted(active_indices)
                if index not in row_deltas
            ]
            if missing:
                raise SystemExit(f"D_{rank} missing active deltas: {missing}")
            negative_linear_sum = sum(
                value * row_deltas[index]
                for index, value in enumerate(linear_coeffs)
                if value < 0 and index >= rank
            )
            negative_quadratic_sum = sum(
                value * row_deltas[first] * row_deltas[second]
                for first, second, value in negative_quadratic
            )
            lower_margin = (
                stable_diff + negative_linear_sum + negative_quadratic_sum
            )
            category = "D-exact-m30-window"
        else:
            if active_indices:
                raise SystemExit(
                    f"D_{rank} has active negative terms without exact deltas: "
                    f"{sorted(active_indices)}"
                )
            lower_margin = stable_diff
            category = "D-stable-only"

        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")
        category_margins[category].append((lower_margin, "D", rank))

    for family, rank in rows:
        if family in {"B", "C"}:
            check_bc_row(family, rank)
        else:
            check_d_row(rank)

    total_checked = sum(len(margins) for margins in category_margins.values())
    expected_count = 530 if bc_only else 796
    if total_checked != expected_count:
        raise SystemExit(
            f"row count mismatch: checked {total_checked}, expected {expected_count}"
        )

    print(
        "post-m29 B/C first-hit lower-bound certificate"
        if bc_only
        else "post-m29 first-hit lower-bound certificate"
    )
    print(f"  finite cutoff: {cutoff}")
    print("  scope: one Chain step at N_hit^(29)(G) for each listed row")
    print(f"  rows checked: {total_checked}")
    print(f"  stable moments generated through m_{max_moment}")
    print("  delta logs:")
    for path in bc_delta_logs:
        print(f"    B/C {path}")
    if not bc_only:
        for path in d_delta_logs:
            print(f"    D {path}")
    for category in sorted(category_margins):
        margins = category_margins[category]
        min_margin, min_family, min_rank = min(margins)
        by_family: dict[str, list[int]] = defaultdict(list)
        for _margin, family, rank in margins:
            by_family[family].append(rank)
        pieces = [
            f"{family}: {rank_range_text(sorted(by_family[family]))}"
            for family in ("B", "C", "D")
            if by_family[family]
        ]
        print(
            f"  {category}: {len(margins)} rows; "
            f"{'; '.join(pieces)}; "
            f"minimum margin {min_margin} at {min_family}_{min_rank}"
        )


def d_m25_reused_delta_certificate(delta_logs: list[Path]) -> None:
    """Check m=25 D_4..D_13 rows from previously certified determinant deltas."""

    d_reused_delta_certificate(delta_logs, 25, 4, 13)


def bc_reused_window_certificate(
    delta_logs: list[Path],
    chain_m: int,
    low_rank: int,
    high_rank: int,
    *,
    stable_quadratic_tail: bool = False,
) -> None:
    """Check B/C rows from previously certified exact deltas."""

    if low_rank < 2 or high_rank < low_rank:
        raise SystemExit("invalid B/C reused-window rank range")
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]

    deltas = parse_bc_delta_logs(delta_logs)
    deltas.update(BC16_REUSED_DELTAS)

    if stable_quadratic_tail:
        print(f"m={chain_m} reused B/C stable-quadratic-tail certificate")
    else:
        print(f"m={chain_m} reused B/C exact-window certificate")
    print(f"stable Chain diff D({chain_m}) = {stable_diff}")
    print(f"certified rows: B_{low_rank}..B_{high_rank} and C_{low_rank}..C_{high_rank}")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")

    def row_margin(family: str, rank: int, exact_through: int) -> int:
        first_boundary = 2 * rank + 2
        row_deltas = {
            moment: value
            for moment, value in deltas.get((family, rank), {}).items()
            if moment <= exact_through
        }
        missing = [
            moment
            for moment in range(first_boundary, exact_through + 1)
            if moment not in row_deltas
        ]
        if missing:
            raise SystemExit(
                f"{family}_{rank} missing exact deltas through m_{exact_through}: {missing}"
            )
        positive_deltas = {
            moment: value for moment, value in row_deltas.items() if value > 0
        }
        if positive_deltas:
            raise SystemExit(
                f"{family}_{rank} boundary corrections are positive: {positive_deltas}"
            )
        active_negative_quadratic = [
            (first, second, value)
            for first, second, value in negative_quadratic
            if first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside and not stable_quadratic_tail:
            raise SystemExit(
                f"{family}_{rank} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in row_deltas.items()
        )
        negative_quadratic_exact = sum(
            value * row_deltas[first] * row_deltas[second]
            for first, second, value in active_negative_quadratic
            if first in row_deltas and second in row_deltas
        )
        negative_quadratic_tail_bound = sum(
            value * stable[first] * stable[second]
            for first, second, value in active_negative_quadratic
            if first not in row_deltas or second not in row_deltas
        )
        if not stable_quadratic_tail:
            negative_quadratic_tail_bound = 0
        positive_quadratic_exact = sum(
            value * row_deltas[first] * row_deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in row_deltas and second in row_deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + negative_quadratic_tail_bound
            + positive_quadratic_exact
        )
        if stable_quadratic_tail:
            print(
                f"{family}_{rank} exact through m_{exact_through}; "
                f"tail = {negative_linear_bound}; "
                f"linear exact = {linear_exact}; "
                f"negative quadratic exact = {negative_quadratic_exact}; "
                f"negative quadratic tail bound = {negative_quadratic_tail_bound}; "
                f"outside negative quadratic pairs = {outside}; "
                f"positive quadratic exact = {positive_quadratic_exact}; "
                f"lower margin = {lower_margin}"
            )
        else:
            print(
                f"{family}_{rank} exact through m_{exact_through}; "
                f"tail = {negative_linear_bound}; "
                f"linear exact = {linear_exact}; "
                f"negative quadratic exact = {negative_quadratic_exact}; "
                f"positive quadratic exact = {positive_quadratic_exact}; "
                f"lower margin = {lower_margin}"
            )
        return lower_margin

    margins: list[tuple[int, str, int]] = []
    for rank in range(high_rank, low_rank - 1, -1):
        for family in ("B", "C"):
            available = max(deltas.get((family, rank), {}).keys(), default=-1)
            if available < 0:
                raise SystemExit(f"missing deltas for {family}_{rank}")
            margin = row_margin(family, rank, available)
            if margin <= 0:
                raise SystemExit(f"{family}_{rank} lower margin is not positive")
            margins.append((margin, family, rank))

    min_margin, min_family, min_rank = min(margins)
    print(
        f"minimum lower margin = {min_family}_{min_rank}: {min_margin}"
    )
    print(f"row_count {len(margins)}")


def bc_m25_reused_window_certificate(delta_logs: list[Path]) -> None:
    """Check the m=25 B/C_5..B/C_16 rows from previously certified deltas."""

    bc_reused_window_certificate(delta_logs, 25, 5, 16)


def d_m21_high_middle_linear_certificate() -> None:
    """Certify the m=21 D-rank >= 27 linear monotonicity shortcut."""

    chain_m = 21
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    coeffs = chain_diff_linear_coefficients(stable, chain_m)
    negative = [index for index, value in enumerate(coeffs) if value < 0]
    negative_high = [
        index for index in range(27, moment_max + 1) if coeffs[index] < 0
    ]
    min_high_index = min(range(27, moment_max + 1), key=lambda index: coeffs[index])

    print("m=21 D high-middle linear certificate")
    print(f"stable Chain diff D(21) = {stable_diff}")
    print(f"negative linear coefficient indices = {negative}")
    print(f"first positive-only D rank = {max(negative) + 1}")
    print(
        "minimum linear coefficient on indices 27..47 = "
        f"L_{min_high_index} = {coeffs[min_high_index]}"
    )
    print("quadratic D-delta term vanishes for rank >= 27 because 27+27 > 47")
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if negative != [21, 22, 23, 24, 25, 26]:
        raise SystemExit("unexpected negative linear coefficient support")
    if negative_high:
        raise SystemExit(f"negative high linear coefficients: {negative_high}")


def d_m22_high_linear_certificate() -> None:
    """Certify the m=22 D-rank >= 29 linear monotonicity shortcut."""

    chain_m = 22
    first_rank = 29
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    coeffs = chain_diff_linear_coefficients(stable, chain_m)
    negative = [index for index, value in enumerate(coeffs) if value < 0]
    negative_high = [
        index for index in range(first_rank, moment_max + 1) if coeffs[index] < 0
    ]
    min_high_index = min(
        range(first_rank, moment_max + 1), key=lambda index: coeffs[index]
    )

    print("m=22 D high linear certificate")
    print(f"stable Chain diff D(22) = {stable_diff}")
    print(f"negative linear coefficient indices = {negative}")
    print(f"first positive-only D rank = {max(negative) + 1}")
    print(
        "minimum linear coefficient on indices 29..49 = "
        f"L_{min_high_index} = {coeffs[min_high_index]}"
    )
    print("quadratic D-delta term vanishes for rank >= 29 because 29+29 > 49")
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if negative != [22, 23, 24, 25, 26, 27, 28]:
        raise SystemExit("unexpected negative linear coefficient support")
    if negative_high:
        raise SystemExit(f"negative high linear coefficients: {negative_high}")


def d_m22_middle_lower_certificate() -> None:
    """Certify the m=22 D_18..D_28 lower-bound shortcut."""

    chain_m = 22
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)

    negative_linear = [
        (index, value) for index, value in enumerate(linear_coeffs) if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_18 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 18 and second >= 18
    ]

    print("m=22 D middle lower certificate")
    print(f"stable Chain diff D(22) = {stable_diff}")
    print(
        "negative linear coefficient indices = "
        f"{[index for index, _value in negative_linear]}"
    )
    print(
        "negative quadratic pairs active from rank >= 18 = "
        f"{active_quadratic_from_18}"
    )

    expected_pairs = [
        (18, 29),
        (19, 28),
        (20, 27),
        (22, 27),
        (23, 26),
        (24, 25),
    ]
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [index for index, _value in negative_linear] != [22, 23, 24, 25, 26, 27, 28]:
        raise SystemExit("unexpected negative linear coefficient support")
    if [(first, second) for first, second, _value in active_quadratic_from_18] != expected_pairs:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[int, int]] = []
    for rank in range(28, 17, -1):
        active_indices = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        deltas = {
            index: d_determinant_correction(rank, index)
            for index in sorted(active_indices)
        }
        negative_delta = {
            index: value for index, value in deltas.items() if value < 0
        }
        if negative_delta:
            raise SystemExit(f"negative D correction at rank {rank}: {negative_delta}")

        negative_linear_sum = sum(
            value * deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_terms: list[tuple[int, int, int, int]] = []
        negative_quadratic_sum = 0
        for first, second, value in negative_quadratic:
            if first >= rank and second >= rank:
                contribution = value * deltas[first] * deltas[second]
                negative_quadratic_sum += contribution
                if contribution != 0:
                    negative_quadratic_terms.append(
                        (first, second, value, contribution)
                    )

        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        margins.append((rank, lower_margin))
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}",
            flush=True,
        )
        print(f"  active negative quadratic terms = {negative_quadratic_terms}", flush=True)
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")

    min_rank, min_margin = min(margins, key=lambda item: item[1])
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")


def d_m22_low_middle_lower_certificate() -> None:
    """Certify the m=22 D_12..D_17 lower-bound shortcut."""

    chain_m = 22
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)

    negative_linear = [
        (index, value) for index, value in enumerate(linear_coeffs) if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_12 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 12 and second >= 12
    ]

    print("m=22 D low-middle lower certificate")
    print(f"stable Chain diff D(22) = {stable_diff}")
    print(
        "negative linear coefficient indices = "
        f"{[index for index, _value in negative_linear]}"
    )
    print(
        "negative quadratic pairs active from rank >= 12 = "
        f"{active_quadratic_from_12}"
    )

    expected_pairs = [
        (12, 35),
        (13, 34),
        (14, 33),
        (15, 32),
        (16, 31),
        (17, 30),
        (18, 29),
        (19, 28),
        (20, 27),
        (22, 27),
        (23, 26),
        (24, 25),
    ]
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [index for index, _value in negative_linear] != [22, 23, 24, 25, 26, 27, 28]:
        raise SystemExit("unexpected negative linear coefficient support")
    if [(first, second) for first, second, _value in active_quadratic_from_12] != expected_pairs:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[int, int]] = []
    for rank in range(17, 11, -1):
        print(f"begin D_{rank}", flush=True)
        active_indices = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        deltas: dict[int, int] = {}
        for index in sorted(active_indices):
            deltas[index] = d_determinant_correction_with_progress(rank, index)
        negative_delta = {
            index: value for index, value in deltas.items() if value < 0
        }
        if negative_delta:
            raise SystemExit(f"negative D correction at rank {rank}: {negative_delta}")

        negative_linear_sum = sum(
            value * deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_sum = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in negative_quadratic
            if first >= rank and second >= rank
        )
        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        margins.append((rank, lower_margin))
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}",
            flush=True,
        )
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")

    min_rank, min_margin = min(margins, key=lambda item: item[1])
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")


def d_m21_middle_lower_certificate() -> None:
    """Certify the m=21 D_20..D_26 lower-bound shortcut."""

    chain_m = 21
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)

    negative_linear = [
        (index, value) for index, value in enumerate(linear_coeffs) if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_20 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 20 and second >= 20
    ]

    print("m=21 D middle lower certificate")
    print(f"stable Chain diff D(21) = {stable_diff}")
    print(
        "negative linear coefficient indices = "
        f"{[index for index, _value in negative_linear]}"
    )
    print(
        "negative quadratic pairs active from rank >= 20 = "
        f"{active_quadratic_from_20}"
    )

    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [index for index, _value in negative_linear] != [21, 22, 23, 24, 25, 26]:
        raise SystemExit("unexpected negative linear coefficient support")
    if [(first, second) for first, second, _value in active_quadratic_from_20] != [
        (21, 26),
        (22, 25),
        (23, 24),
    ]:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[int, int]] = []
    for rank in range(20, 27):
        active_indices = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        deltas = {
            index: d_determinant_correction(rank, index)
            for index in sorted(active_indices)
        }
        negative_delta = {
            index: value for index, value in deltas.items() if value < 0
        }
        if negative_delta:
            raise SystemExit(f"negative D correction at rank {rank}: {negative_delta}")

        negative_linear_sum = sum(
            value * deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_terms: list[tuple[int, int, int, int]] = []
        negative_quadratic_sum = 0
        for first, second, value in negative_quadratic:
            if first >= rank and second >= rank:
                contribution = value * deltas[first] * deltas[second]
                negative_quadratic_sum += contribution
                if contribution != 0:
                    negative_quadratic_terms.append(
                        (first, second, value, contribution)
                    )

        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        margins.append((rank, lower_margin))
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}"
        )
        print(f"  active negative quadratic terms = {negative_quadratic_terms}")
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")

    min_rank, min_margin = min(margins, key=lambda item: item[1])
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")


def d_m21_low_middle_lower_certificate() -> None:
    """Certify the m=21 D_12..D_18 lower-bound shortcut."""

    chain_m = 21
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)

    negative_linear = [
        (index, value) for index, value in enumerate(linear_coeffs) if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_12 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 12 and second >= 12
    ]

    print("m=21 D low-middle lower certificate")
    print(f"stable Chain diff D(21) = {stable_diff}")
    print(
        "negative linear coefficient indices = "
        f"{[index for index, _value in negative_linear]}"
    )
    print(
        "negative quadratic pairs active from rank >= 12 = "
        f"{active_quadratic_from_12}"
    )

    expected_pairs = [
        (12, 33),
        (13, 32),
        (14, 31),
        (15, 30),
        (16, 29),
        (17, 28),
        (18, 27),
        (19, 26),
        (21, 26),
        (22, 25),
        (23, 24),
    ]
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [index for index, _value in negative_linear] != [21, 22, 23, 24, 25, 26]:
        raise SystemExit("unexpected negative linear coefficient support")
    if [(first, second) for first, second, _value in active_quadratic_from_12] != expected_pairs:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[int, int]] = []
    for rank in range(18, 11, -1):
        active_indices = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        deltas = {
            index: d_determinant_correction(rank, index)
            for index in sorted(active_indices)
        }
        negative_delta = {
            index: value for index, value in deltas.items() if value < 0
        }
        if negative_delta:
            raise SystemExit(f"negative D correction at rank {rank}: {negative_delta}")

        negative_linear_sum = sum(
            value * deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_sum = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in negative_quadratic
            if first >= rank and second >= rank
        )
        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        margins.append((rank, lower_margin))
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}",
            flush=True,
        )
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")

    min_rank, min_margin = min(margins, key=lambda item: item[1])
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")


def d_m21_bottom_lower_certificate() -> None:
    """Certify the m=21 D_9..D_11 lower-bound shortcut."""

    chain_m = 21
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)

    negative_linear = [
        (index, value) for index, value in enumerate(linear_coeffs) if value < 0
    ]
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_9 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 9 and second >= 9
    ]

    print("m=21 D bottom lower certificate")
    print(f"stable Chain diff D(21) = {stable_diff}")
    print(
        "negative linear coefficient indices = "
        f"{[index for index, _value in negative_linear]}"
    )
    print(
        "negative quadratic pairs active from rank >= 9 = "
        f"{active_quadratic_from_9}"
    )

    expected_pairs = [
        (9, 36),
        (10, 35),
        (11, 34),
        (12, 33),
        (13, 32),
        (14, 31),
        (15, 30),
        (16, 29),
        (17, 28),
        (18, 27),
        (19, 26),
        (21, 26),
        (22, 25),
        (23, 24),
    ]
    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [index for index, _value in negative_linear] != [21, 22, 23, 24, 25, 26]:
        raise SystemExit("unexpected negative linear coefficient support")
    if [(first, second) for first, second, _value in active_quadratic_from_9] != expected_pairs:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[int, int]] = []
    for rank in range(11, 8, -1):
        active_indices = {
            index for index, _value in negative_linear if index >= rank
        }
        for first, second, _value in negative_quadratic:
            if first >= rank and second >= rank:
                active_indices.add(first)
                active_indices.add(second)
        deltas: dict[int, int] = {}
        print(f"begin D_{rank} lower-bound row", flush=True)
        for index in sorted(active_indices):
            print(
                f"  computing determinant correction m_{index} "
                f"(depth {index - rank})",
                flush=True,
            )
            deltas[index] = d_determinant_correction_with_progress(rank, index)

        negative_delta = {
            index: value for index, value in deltas.items() if value < 0
        }
        if negative_delta:
            raise SystemExit(f"negative D correction at rank {rank}: {negative_delta}")

        negative_linear_sum = sum(
            value * deltas[index]
            for index, value in negative_linear
            if index >= rank
        )
        negative_quadratic_sum = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in negative_quadratic
            if first >= rank and second >= rank
        )
        lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum
        margins.append((rank, lower_margin))
        print(
            f"D_{rank} negative linear = {negative_linear_sum}; "
            f"negative quadratic = {negative_quadratic_sum}; "
            f"lower margin = {lower_margin}",
            flush=True,
        )
        if lower_margin <= 0:
            raise SystemExit(f"D_{rank} lower margin is not positive")

    min_rank, min_margin = min(margins, key=lambda item: item[1])
    print(f"minimum lower margin = D_{min_rank}: {min_margin}")


def b_boundary_exact_correction(rank: int, moment_index: int) -> int:
    """SO(2b+1)-minus-stable B-boundary moment correction."""

    matrix_size = 2 * rank + 1
    first_boundary = matrix_size + 1
    if moment_index < first_boundary:
        return 0
    if moment_index <= first_boundary + 1:
        return b_boundary_first_two_correction(rank, moment_index)
    total = 0
    for half_partition in integer_partitions(moment_index):
        if len(half_partition) <= matrix_size:
            continue
        target = tuple(2 * part for part in half_partition)
        total += pieri_e2_coefficient(moment_index, target)
    return -total


def c_boundary_exact_correction(rank: int, moment_index: int) -> int:
    """Sp(b)-minus-stable C-boundary moment correction."""

    first_boundary = 2 * rank + 2
    if moment_index < first_boundary:
        return 0
    if moment_index <= first_boundary + 1:
        return c_boundary_first_two_correction(rank, moment_index)
    total = 0
    for paired_partition in integer_partitions(moment_index):
        if len(paired_partition) <= rank:
            continue
        target = tuple(row for part in paired_partition for row in (part, part))
        total += pieri_e2_coefficient(moment_index, conjugate_partition(target))
    return -total


def b_boundary_exact_correction_with_progress(rank: int, moment_index: int) -> int:
    """Progress-printing version of the exact B-boundary correction."""

    matrix_size = 2 * rank + 1
    first_boundary = matrix_size + 1
    if moment_index < first_boundary:
        return 0
    if moment_index <= first_boundary + 1:
        return b_boundary_first_two_correction(rank, moment_index)

    targets = [
        tuple(2 * part for part in half_partition)
        for half_partition in integer_partitions(moment_index)
        if len(half_partition) > matrix_size
    ]
    total = 0
    start = time()
    for index, target in enumerate(targets, 1):
        total += pieri_e2_coefficient(moment_index, target)
        if index == len(targets) or index % 100 == 0:
            print(
                f"    B_{rank} m_{moment_index} target {index}/{len(targets)} "
                f"elapsed={time() - start:.1f}s",
                flush=True,
            )
    return -total


def c_boundary_exact_correction_with_progress(rank: int, moment_index: int) -> int:
    """Progress-printing version of the exact C-boundary correction."""

    first_boundary = 2 * rank + 2
    if moment_index < first_boundary:
        return 0
    if moment_index <= first_boundary + 1:
        return c_boundary_first_two_correction(rank, moment_index)

    targets = [
        conjugate_partition(tuple(row for part in paired_partition for row in (part, part)))
        for paired_partition in integer_partitions(moment_index)
        if len(paired_partition) > rank
    ]
    total = 0
    start = time()
    for index, target in enumerate(targets, 1):
        total += pieri_e2_coefficient(moment_index, target)
        if index == len(targets) or index % 100 == 0:
            print(
                f"    C_{rank} m_{moment_index} target {index}/{len(targets)} "
                f"elapsed={time() - start:.1f}s",
                flush=True,
            )
    return -total


def bc_m21_high_boundary_lower_certificate() -> None:
    """Certify the m=21 B/C rank 9..16 boundary lower bound."""

    chain_m = 21
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]

    print("m=21 B/C high-boundary lower certificate")
    print(f"stable Chain diff D(21) = {stable_diff}")
    print(
        "positive linear coefficient indices at or above 20 = "
        f"{[index for index, value in enumerate(linear_coeffs) if index >= 20 and value > 0]}"
    )
    print(
        "negative quadratic pairs active from rank-boundary >= 20 = "
        f"{[(a, b, value) for a, b, value in negative_quadratic if a >= 20 and b >= 20]}"
    )

    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")

    margins: list[tuple[str, int, int]] = []
    tail_from_30 = -sum(
        value * stable[index]
        for index, value in enumerate(linear_coeffs)
        if index >= 30 and value > 0
    )
    for family, correction in (
        ("B", b_boundary_exact_correction),
        ("C", c_boundary_exact_correction),
    ):
        for rank in range(9, 14):
            first_boundary = 2 * rank + 2
            deltas = {
                index: correction(rank, index)
                for index in range(first_boundary, 30)
            }
            active_negative_quadratic = [
                (a, b, value)
                for a, b, value in negative_quadratic
                if a >= first_boundary and b >= first_boundary
            ]
            high_negative_quadratic = [
                (a, b, value)
                for a, b, value in active_negative_quadratic
                if a >= 30 or b >= 30
            ]
            if high_negative_quadratic:
                raise SystemExit(
                    f"{family}_{rank} has an unbounded high negative quadratic term"
                )
            linear_exact = sum(
                linear_coeffs[index] * value for index, value in deltas.items()
            )
            quadratic_exact = sum(
                value * deltas[a] * deltas[b]
                for a, b, value in active_negative_quadratic
            )
            lower_margin = (
                stable_diff + tail_from_30 + linear_exact + quadratic_exact
            )
            delta_text = "; ".join(
                f"Delta_{index} = {value}" for index, value in deltas.items()
            )
            print(
                f"{family}_{rank} first boundary = {first_boundary}; "
                f"{delta_text}; tail-from-30 bound = {tail_from_30}; "
                f"negative quadratic exact = {quadratic_exact}; "
                f"lower margin = {lower_margin}"
            )
            if any(value > 0 for value in deltas.values()):
                raise SystemExit(f"{family}_{rank} boundary corrections are not nonpositive")
            if lower_margin <= 0:
                raise SystemExit(f"{family}_{rank} lower margin is not positive")
            margins.append((f"{family}_{rank}", first_boundary, lower_margin))

    for rank in range(14, 17):
        first_boundary = 2 * rank + 2
        negative_linear_bound = -sum(
            value * stable[index]
            for index, value in enumerate(linear_coeffs)
            if index >= first_boundary and value > 0
        )
        lower_margin = stable_diff + negative_linear_bound
        margins.append((f"B/C rank {rank}", first_boundary, lower_margin))
        print(
            f"B/C rank {rank} first boundary = {first_boundary}; "
            f"negative linear bound = {negative_linear_bound}; "
            f"lower margin = {lower_margin}"
        )
        if lower_margin <= 0:
            raise SystemExit(f"B/C rank {rank} lower margin is not positive")

    min_label, min_boundary, min_margin = min(margins, key=lambda item: item[2])
    print(
        "minimum lower margin = "
        f"{min_label} first boundary {min_boundary}: {min_margin}"
    )


def bc_m22_boundary_lower_certificate() -> None:
    """Certify the m=22 B/C rank 9..23 boundary lower bound."""

    chain_m = 22
    moment_max = 2 * chain_m + 5
    exact_through = 31
    tail_start_floor = exact_through + 1
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]

    print("m=22 B/C boundary lower certificate")
    print(f"stable Chain diff D(22) = {stable_diff}")
    print(
        "positive linear coefficient indices at or above 20 = "
        f"{[index for index, value in enumerate(linear_coeffs) if index >= 20 and value > 0]}"
    )
    print(
        "negative quadratic pairs active from rank-boundary >= 20 = "
        f"{[(a, b, value) for a, b, value in negative_quadratic if a >= 20 and b >= 20]}"
    )
    print(f"exact B/C boundary corrections used through index {exact_through}")

    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")

    margins: list[tuple[str, int, int]] = []
    for family, correction in (
        ("B", b_boundary_exact_correction),
        ("C", c_boundary_exact_correction),
    ):
        for rank in range(9, 24):
            first_boundary = 2 * rank + 2
            deltas = {
                index: correction(rank, index)
                for index in range(first_boundary, min(exact_through, moment_max) + 1)
            }
            active_negative_quadratic = [
                (a, b, value)
                for a, b, value in negative_quadratic
                if a >= first_boundary and b >= first_boundary
            ]
            high_negative_quadratic = [
                (a, b, value)
                for a, b, value in active_negative_quadratic
                if a > exact_through or b > exact_through
            ]
            if high_negative_quadratic:
                raise SystemExit(
                    f"{family}_{rank} has an unbounded high negative quadratic term"
                )
            if any(value > 0 for value in deltas.values()):
                raise SystemExit(f"{family}_{rank} boundary corrections are not nonpositive")

            tail_start = max(first_boundary, tail_start_floor)
            negative_linear_bound = -sum(
                value * stable[index]
                for index, value in enumerate(linear_coeffs)
                if index >= tail_start and value > 0
            )
            linear_exact = sum(
                linear_coeffs[index] * value for index, value in deltas.items()
            )
            quadratic_exact = sum(
                value * deltas[a] * deltas[b]
                for a, b, value in active_negative_quadratic
            )
            lower_margin = (
                stable_diff + negative_linear_bound + linear_exact + quadratic_exact
            )
            print(
                f"{family}_{rank} first boundary = {first_boundary}; "
                f"exact corrections = {len(deltas)}; "
                f"tail start = {tail_start}; "
                f"negative linear bound = {negative_linear_bound}; "
                f"negative quadratic exact = {quadratic_exact}; "
                f"lower margin = {lower_margin}"
            )
            if lower_margin <= 0:
                raise SystemExit(f"{family}_{rank} lower margin is not positive")
            margins.append((f"{family}_{rank}", first_boundary, lower_margin))

    min_label, min_boundary, min_margin = min(margins, key=lambda item: item[2])
    print(
        "minimum lower margin = "
        f"{min_label} first boundary {min_boundary}: {min_margin}"
    )


def bc_m22_low_boundary_lower_certificate() -> None:
    """Certify the m=22 B/C rank 6..8 boundary lower bound."""

    chain_m = 22
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
    quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
    negative_quadratic = [
        (first, second, value)
        for (first, second), value in sorted(quadratic_coeffs.items())
        if value < 0
    ]
    active_quadratic_from_14 = [
        (first, second, value)
        for first, second, value in negative_quadratic
        if first >= 14 and second >= 14
    ]
    expected_pairs = [
        (14, 33),
        (15, 32),
        (16, 31),
        (17, 30),
        (18, 29),
        (19, 28),
        (20, 27),
        (22, 27),
        (23, 26),
        (24, 25),
    ]

    print("m=22 B/C low-boundary lower certificate")
    print(f"stable Chain diff D(22) = {stable_diff}")
    print(
        "positive linear coefficient indices at or above 14 = "
        f"{[index for index, value in enumerate(linear_coeffs) if index >= 14 and value > 0]}"
    )
    print(
        "negative quadratic pairs active from rank-boundary >= 14 = "
        f"{active_quadratic_from_14}"
    )

    if stable_diff <= 0:
        raise SystemExit("stable Chain diff is not positive")
    if [(first, second) for first, second, _value in active_quadratic_from_14] != expected_pairs:
        raise SystemExit("unexpected active negative quadratic support")

    margins: list[tuple[str, int, int]] = []
    for family, correction in (
        ("B", b_boundary_exact_correction_with_progress),
        ("C", c_boundary_exact_correction_with_progress),
    ):
        for rank in range(8, 5, -1):
            first_boundary = 2 * rank + 2
            active_negative_quadratic = [
                (first, second, value)
                for first, second, value in negative_quadratic
                if first >= first_boundary and second >= first_boundary
            ]
            exact_through = max(
                [31]
                + [
                    endpoint
                    for first, second, _value in active_negative_quadratic
                    for endpoint in (first, second)
                ]
            )
            high_negative_quadratic = [
                (first, second, value)
                for first, second, value in active_negative_quadratic
                if first > exact_through or second > exact_through
            ]
            if high_negative_quadratic:
                raise SystemExit(
                    f"{family}_{rank} has an unbounded high negative quadratic term"
                )

            print(
                f"begin {family}_{rank}: first boundary = {first_boundary}; "
                f"exact through m_{exact_through}",
                flush=True,
            )
            deltas: dict[int, int] = {}
            for moment_index in range(first_boundary, exact_through + 1):
                deltas[moment_index] = correction(rank, moment_index)
                print(
                    f"  {family}_{rank} Delta_{moment_index} = "
                    f"{deltas[moment_index]}",
                    flush=True,
                )
            if any(value > 0 for value in deltas.values()):
                raise SystemExit(f"{family}_{rank} boundary corrections are not nonpositive")

            tail_start = exact_through + 1
            negative_linear_bound = -sum(
                value * stable[index]
                for index, value in enumerate(linear_coeffs)
                if index >= tail_start and value > 0
            )
            linear_exact = sum(
                linear_coeffs[index] * value for index, value in deltas.items()
            )
            quadratic_exact = sum(
                value * deltas[first] * deltas[second]
                for first, second, value in active_negative_quadratic
            )
            lower_margin = (
                stable_diff
                + negative_linear_bound
                + linear_exact
                + quadratic_exact
            )
            print(
                f"{family}_{rank} tail = {negative_linear_bound}; "
                f"linear exact = {linear_exact}; "
                f"negative quadratic exact = {quadratic_exact}; "
                f"lower margin = {lower_margin}",
                flush=True,
            )
            if lower_margin <= 0:
                raise SystemExit(f"{family}_{rank} lower margin is not positive")
            margins.append((f"{family}_{rank}", first_boundary, lower_margin))

    min_label, min_boundary, min_margin = min(margins, key=lambda item: item[2])
    print(
        "minimum lower margin = "
        f"{min_label} first boundary {min_boundary}: {min_margin}"
    )


def d_determinant_correction(rank: int, moment_index: int) -> int:
    """SO(2b)-minus-stable adjoint moment correction."""

    depth = moment_index - rank
    if depth < 0:
        return 0
    if depth == 0:
        return 1
    if depth == 1:
        return binom_integer(rank, 2)
    if depth == 2:
        return rank * (rank + 1) * (rank * rank + rank + 2) // 8
    if depth == 3:
        return (
            (rank - 1)
            * (rank + 2)
            * (rank**4 + 8 * rank**3 + 25 * rank * rank + 18 * rank - 24)
            // 48
        )
    if depth == 4:
        return (
            (rank + 3)
            * (
                rank**7
                + 17 * rank**6
                + 103 * rank**5
                + 203 * rank**4
                - 104 * rank**3
                - 460 * rank * rank
                - 48 * rank
                + 384
            )
            // 384
        )
    return d_determinant_correction_pieri(rank, depth)


def d_determinant_correction_with_progress(rank: int, moment_index: int) -> int:
    depth = moment_index - rank
    if depth <= 4:
        return d_determinant_correction(rank, moment_index)
    matrix_size = 2 * rank
    power = rank + depth
    parts = list(integer_partitions(depth))
    total = 0
    start = time()
    for index, excess_partition in enumerate(parts, 1):
        target = tuple(
            [1 + 2 * row for row in excess_partition]
            + [1] * (matrix_size - len(excess_partition))
        )
        total += pieri_e2_coefficient(power, target)
        if index == len(parts) or index % 25 == 0:
            print(
                f"    depth {depth}: partition {index}/{len(parts)} "
                f"elapsed={time() - start:.1f}s",
                flush=True,
            )
    return total


def b_boundary_first_two_correction(rank: int, moment_index: int) -> int:
    """SO(2b+1)-minus-stable first B-boundary moment correction."""

    matrix_size = 2 * rank + 1
    first_boundary = matrix_size + 1
    if moment_index < first_boundary:
        return 0
    if moment_index == first_boundary:
        return -pieri_e2_coefficient(moment_index, (2,) * first_boundary)
    if moment_index == first_boundary + 1:
        return -(
            pieri_e2_coefficient(moment_index, (4,) + (2,) * matrix_size)
            + pieri_e2_coefficient(moment_index, (2,) * (matrix_size + 2))
        )
    raise ValueError("only the first two B-boundary moments are implemented")


def b_boundary_m15_rank16_certificate() -> None:
    """Replay the B_16 odd-orthogonal first-boundary row."""

    chain_m = 15
    rank = 16
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = -40215585765779526
    print("B boundary m=15 rank-sixteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin B_{rank}=SO({2 * rank + 1})")
    for moment_index in (2 * rank + 2, 2 * rank + 3):
        correction = b_boundary_first_two_correction(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"B_{rank}=SO({2 * rank + 1})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  B_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected B_16 Chain correction")
    if row_diff <= 0:
        raise SystemExit("B_16 boundary Chain diff is not positive")


def c_boundary_first_two_correction(rank: int, moment_index: int) -> int:
    """Sp(b)-minus-stable first C-boundary moment correction."""

    first_boundary = 2 * rank + 2
    if moment_index < first_boundary:
        return 0
    if moment_index > first_boundary + 1:
        raise ValueError("only the first two C-boundary moments are implemented")
    total = 0
    for paired_partition in integer_partitions(moment_index):
        if len(paired_partition) != rank + 1:
            continue
        target = tuple(row for part in paired_partition for row in (part, part))
        total += pieri_e2_coefficient(moment_index, conjugate_partition(target))
    return -total


def c_boundary_m15_rank16_certificate() -> None:
    """Replay the C_16 symplectic first-boundary row."""

    chain_m = 15
    rank = 16
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = -7092579055254392700002
    print("C boundary m=15 rank-sixteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin C_{rank}=Sp({rank})")
    for moment_index in (2 * rank + 2, 2 * rank + 3):
        correction = c_boundary_first_two_correction(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"C_{rank}=Sp({rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  C_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected C_16 Chain correction")
    if row_diff <= 0:
        raise SystemExit("C_16 boundary Chain diff is not positive")


def d_boundary_m15_depth8_certificate() -> None:
    """Replay the D_35 through D_27 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {
        35: 2,
        34: 1122,
        33: 316200,
        32: 59175156,
        31: 8216171500,
        30: 897437738444,
        29: 79903676760500,
        28: 5935469843032484,
        27: 373756281668634552,
    }
    print("D boundary m=15 depth-eight certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in range(35, 26, -1):
        row = stable.copy()
        corrections: dict[int, int] = {}
        for moment_index in range(rank, moment_max + 1):
            correction = d_determinant_correction(rank, moment_index)
            corrections[moment_index] = correction
            row[moment_index] += correction
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        correction_text = ", ".join(
            f"m_{moment_index} += {correction}"
            for moment_index, correction in corrections.items()
        )
        print(f"D_{rank}=SO({2 * rank}): {correction_text}")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_depth11_certificate() -> None:
    """Replay the D_35 through D_24 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {
        35: 2,
        34: 1122,
        33: 316200,
        32: 59175156,
        31: 8216171500,
        30: 897437738444,
        29: 79903676760500,
        28: 5935469843032484,
        27: 373756281668634552,
        26: 20175274148084539064,
        25: 940942326400224424148,
        24: 38121899420476396012556,
    }
    print("D boundary m=15 depth-eleven certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in range(35, 23, -1):
        row = stable.copy()
        corrections: dict[int, int] = {}
        for moment_index in range(rank, moment_max + 1):
            correction = d_determinant_correction(rank, moment_index)
            corrections[moment_index] = correction
            row[moment_index] += correction
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        correction_text = ", ".join(
            f"m_{moment_index} += {correction}"
            for moment_index, correction in corrections.items()
        )
        print(f"D_{rank}=SO({2 * rank}): {correction_text}")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_depth15_certificate() -> None:
    """Replay the D_23 through D_20 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {
        23: 1346433429603175519393796,
        22: 41538310419241788865735076,
        21: 1120047581338576795332371756,
        20: 26380575894047450571141473756,
    }
    print("D boundary m=15 depth-fifteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in range(23, 19, -1):
        row = stable.copy()
        print(f"begin D_{rank}=SO({2 * rank})")
        for moment_index in range(rank, moment_max + 1):
            print(f"  computing m_{moment_index}", flush=True)
            correction = d_determinant_correction(rank, moment_index)
            row[moment_index] += correction
            print(f"  m_{moment_index} += {correction}", flush=True)
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        print(f"D_{rank}=SO({2 * rank})")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_depth16_certificate() -> None:
    """Replay the D_19 determinant-boundary row."""

    chain_m = 15
    rank = 19
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 541769025571438891945770168266
    print("D boundary m=15 depth-sixteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_19 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_19 boundary Chain diff is not positive")


def d_boundary_m15_depth17_certificate() -> None:
    """Replay the D_18 determinant-boundary row."""

    chain_m = 15
    rank = 18
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 9672183498929510843585150678186
    print("D boundary m=15 depth-seventeen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_18 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_18 boundary Chain diff is not positive")


def d_boundary_m15_depth18_certificate() -> None:
    """Replay the D_17 determinant-boundary row."""

    chain_m = 15
    rank = 17
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 149469405784662980644953679985396
    print("D boundary m=15 depth-eighteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_17 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_17 boundary Chain diff is not positive")


def d_boundary_m15_depth19_certificate() -> None:
    """Replay the D_16 determinant-boundary row."""

    chain_m = 15
    rank = 16
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 1987934151693111645655720383622436
    print("D boundary m=15 depth-nineteen certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_16 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_16 boundary Chain diff is not positive")


def d_boundary_m15_depth20_certificate() -> None:
    """Replay the D_15 determinant-boundary row."""

    chain_m = 15
    rank = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 22585579969287948920971197502904876
    print("D boundary m=15 depth-twenty certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_15 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_15 boundary Chain diff is not positive")


def d_boundary_m15_depth21_certificate() -> None:
    """Replay the D_14 determinant-boundary row."""

    chain_m = 15
    rank = 14
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 217107652771456656368979552537899036
    print("D boundary m=15 depth-twenty-one certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_14 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_14 boundary Chain diff is not positive")


def d_boundary_m15_depth22_certificate() -> None:
    """Replay the D_13 determinant-boundary row."""

    chain_m = 15
    rank = 13
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 1744119757413181392572698915333536836
    print("D boundary m=15 depth-twenty-two certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_13 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_13 boundary Chain diff is not positive")


def d_boundary_m15_depth23_certificate() -> None:
    """Replay the D_12 determinant-boundary row."""

    chain_m = 15
    rank = 12
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 11522883607577150101966986764353266596
    print("D boundary m=15 depth-twenty-three certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_12 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_12 boundary Chain diff is not positive")


def d_boundary_m15_depth24_certificate() -> None:
    """Replay the D_11 determinant-boundary row."""

    chain_m = 15
    rank = 11
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 61283860709490438492118310295552637736
    print("D boundary m=15 depth-twenty-four certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_11 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_11 boundary Chain diff is not positive")


def d_boundary_m15_depth25_certificate() -> None:
    """Replay the D_10 determinant-boundary row."""

    chain_m = 15
    rank = 10
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 254761008791285875054337542543008475880
    print("D boundary m=15 depth-twenty-five certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_10 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_10 boundary Chain diff is not positive")


def d_boundary_m15_depth26_certificate() -> None:
    """Replay the D_9 determinant-boundary row."""

    chain_m = 15
    rank = 9
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 793141899605172638000857779019212911708
    print("D boundary m=15 depth-twenty-six certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_9 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_9 boundary Chain diff is not positive")


def d_boundary_m15_depth27_certificate() -> None:
    """Replay the D_8 determinant-boundary row."""

    chain_m = 15
    rank = 8
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 1729879920827145126492482924958574042196
    print("D boundary m=15 depth-twenty-seven certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_8 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_8 boundary Chain diff is not positive")


def d_boundary_m15_depth28_certificate() -> None:
    """Replay the D_7 determinant-boundary row."""

    chain_m = 15
    rank = 7
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 2358025505083962440192612470017321728600
    print("D boundary m=15 depth-twenty-eight certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_7 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_7 boundary Chain diff is not positive")


def d_boundary_m15_depth29_certificate() -> None:
    """Replay the D_6 determinant-boundary row."""

    chain_m = 15
    rank = 6
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 1642070299122660773049079693146862223096
    print("D boundary m=15 depth-twenty-nine certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_6 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_6 boundary Chain diff is not positive")


def d_boundary_m15_depth30_certificate() -> None:
    """Replay the D_5 determinant-boundary row."""

    chain_m = 15
    rank = 5
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 510659785860825641079411590217029794160
    print("D boundary m=15 depth-thirty certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_5 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_5 boundary Chain diff is not positive")


def d_boundary_m15_depth31_certificate() -> None:
    """Replay the D_4 determinant-boundary row."""

    chain_m = 15
    rank = 4
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_chain_correction = 374319086629787669451485698346854129792
    print("D boundary m=15 depth-thirty-one certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    row = stable.copy()
    print(f"begin D_{rank}=SO({2 * rank})")
    for moment_index in range(rank, moment_max + 1):
        print(f"  computing m_{moment_index}", flush=True)
        correction = d_determinant_correction_with_progress(rank, moment_index)
        row[moment_index] += correction
        print(f"  m_{moment_index} += {correction}", flush=True)
    row_diff = chain_diff_integer(row, chain_m)
    chain_correction = row_diff - stable_diff
    print(f"D_{rank}=SO({2 * rank})")
    print(f"  boundary correction in Chain diff = {chain_correction}")
    print(f"  D_{rank} Chain diff D(15) = {row_diff}")
    if chain_correction != expected_chain_correction:
        raise SystemExit("unexpected D_4 Chain correction")
    if row_diff <= 0:
        raise SystemExit("D_4 boundary Chain diff is not positive")


def d_boundary_m15_depth4_certificate() -> None:
    """Replay the D_35 through D_31 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {
        35: 2,
        34: 1122,
        33: 316200,
        32: 59175156,
        31: 8216171500,
    }
    print("D boundary m=15 depth-four certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in (35, 34, 33, 32, 31):
        row = stable.copy()
        corrections: dict[int, int] = {}
        for moment_index in range(rank, moment_max + 1):
            correction = d_determinant_correction(rank, moment_index)
            corrections[moment_index] = correction
            row[moment_index] += correction
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        correction_text = ", ".join(
            f"m_{moment_index} += {correction}"
            for moment_index, correction in corrections.items()
        )
        print(f"D_{rank}=SO({2 * rank}): {correction_text}")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_depth3_certificate() -> None:
    """Replay the D_35 through D_32 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {35: 2, 34: 1122, 33: 316200, 32: 59175156}
    print("D boundary m=15 depth-three certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in (35, 34, 33, 32):
        row = stable.copy()
        corrections: dict[int, int] = {}
        for moment_index in range(rank, moment_max + 1):
            correction = d_determinant_correction(rank, moment_index)
            corrections[moment_index] = correction
            row[moment_index] += correction
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        correction_text = ", ".join(
            f"m_{moment_index} += {correction}"
            for moment_index, correction in corrections.items()
        )
        print(f"D_{rank}=SO({2 * rank}): {correction_text}")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_depth2_certificate() -> None:
    """Replay the D_35, D_34, and D_33 determinant-boundary rows."""

    chain_m = 15
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    stable_diff = chain_diff_integer(stable, chain_m)
    expected_corrections = {35: 2, 34: 1122, 33: 316200}
    print("D boundary m=15 depth-two certificate")
    print(f"moment range m_0..m_{moment_max}")
    print(f"stable Chain diff D(15) = {stable_diff}")
    for rank in (35, 34, 33):
        row = stable.copy()
        corrections: dict[int, int] = {}
        for moment_index in range(rank, moment_max + 1):
            correction = d_determinant_correction(rank, moment_index)
            corrections[moment_index] = correction
            row[moment_index] += correction
        row_diff = chain_diff_integer(row, chain_m)
        chain_correction = row_diff - stable_diff
        correction_text = ", ".join(
            f"m_{moment_index} += {correction}"
            for moment_index, correction in corrections.items()
        )
        print(f"D_{rank}=SO({2 * rank}): {correction_text}")
        print(f"  boundary correction in Chain diff = {chain_correction}")
        print(f"  D_{rank} Chain diff D(15) = {row_diff}")
        if chain_correction != expected_corrections[rank]:
            raise SystemExit(f"unexpected D_{rank} Chain correction")
        if row_diff <= 0:
            raise SystemExit(f"D_{rank} boundary Chain diff is not positive")


def d_boundary_m15_certificate() -> None:
    """Replay the D_35 boundary correction for the first open middle slice."""

    chain_m = 15
    rank = 35
    moment_max = 2 * chain_m + 5
    stable = read_stable_moments(moment_max)
    boundary = stable.copy()
    boundary[rank] += 1
    stable_diff = chain_diff_integer(stable, chain_m)
    boundary_diff = chain_diff_integer(boundary, chain_m)
    correction = boundary_diff - stable_diff
    print("D boundary m=15 certificate")
    print(f"group D_{rank}=SO({2 * rank})")
    print(f"moment range m_0..m_{moment_max}")
    print("moments m_0..m_34 equal the stable SO moments")
    print("Pfaffian boundary correction: m_35(D_35)=m_35(stable)+1")
    print(f"stable Chain diff D(15) = {stable_diff}")
    print(f"boundary correction in Chain diff = {correction}")
    print(f"D_35 Chain diff D(15) = {boundary_diff}")
    if correction != 2:
        raise SystemExit("unexpected boundary correction")
    if boundary_diff <= 0:
        raise SystemExit("D_35 boundary Chain diff is not positive")


def finite_leading_margin_certificate(cutoff: int = TRACE_CUTOFF) -> None:
    """Check the finite residue's leading dominance from the stable edge."""

    rows = finite_residue_rows(cutoff)
    worst: tuple[float, str, int, int, int] | None = None
    for family, rank in rows:
        data = family_data(family, rank)
        min_n, min_margin = post_stable_min_margin(data)
        margin_log10 = min_margin / log(10)
        row = (margin_log10, family, rank, min_n, stable_odd_reach(family, rank))
        if worst is None or row < worst:
            worst = row
    if worst is None:
        raise SystemExit("finite residue is empty")
    margin_log10, family, rank, min_n, stable_edge = worst
    print("finite leading dominance certificate")
    print(f"high-rank trace cutoff: C_G >= {cutoff}")
    print(f"finite rows total: {len(rows)}")
    print("A_G lower bound: Proposition 24.5 with pi < 4")
    print(
        "worst post-stable leading margin: "
        f"{family}_{rank} at odd n={min_n}, "
        f"stable edge={stable_edge}, log10 margin={margin_log10:.12f}"
    )
    if margin_log10 <= 0:
        raise SystemExit("finite leading margin is not positive")


def family_data(family: str, rank: int) -> dict[str, object]:
    b = rank
    if family == "B":
        return {
            "r": b,
            "d": b * (2 * b + 1),
            "c": b,
            "kappa": 2 * b - 1,
            "z": 1,
            "degrees": tuple(2 * i for i in range(1, b + 1)),
            "stable_last_odd": 2 * b - 1,
        }
    if family == "C":
        return {
            "r": b,
            "d": b * (2 * b + 1),
            "c": b,
            "kappa": b + 1,
            "z": 2,
            "degrees": tuple(2 * i for i in range(1, b + 1)),
            "stable_last_odd": 2 * b - 1,
        }
    if family == "D":
        return {
            "r": b,
            "d": b * (2 * b - 1),
            "c": b if b % 2 == 0 else b - 2,
            "kappa": 2 * b - 2,
            "z": 2,
            "degrees": tuple(2 * i for i in range(1, b)) + (b,),
            "stable_last_odd": max(1, b - 3),
        }
    raise ValueError(f"unknown family {family!r}")


def alpha_value(data: dict[str, object]) -> int:
    return int(data["d"]) + 2


def log_a0_lower(data: dict[str, object], pi_upper: int = 4) -> float:
    r = int(data["r"])
    d = int(data["d"])
    kappa = int(data["kappa"])
    z = int(data["z"])
    degrees = data["degrees"]
    if not isinstance(degrees, tuple):
        raise TypeError("degrees must be a tuple")
    a = Fraction(d, 2)
    out = 2 * log(z)
    out += log(8) + log(float(a)) + 2 * log(d)
    out += d * log(Fraction(d, kappa))
    out += 2 * sum(lgamma(degree) for degree in degrees)
    out -= r * log(2 * pi_upper)
    return out


def a0_lower_fraction(data: dict[str, object], pi_upper: int = 4) -> Fraction:
    """Exact rational lower bound for A_G obtained by replacing pi with pi_upper."""

    r = int(data["r"])
    d = int(data["d"])
    kappa = int(data["kappa"])
    z = int(data["z"])
    degrees = data["degrees"]
    if not isinstance(degrees, tuple):
        raise TypeError("degrees must be a tuple")
    out = Fraction(z * z, 1)
    out *= Fraction(8, 1) * Fraction(d, 2) * d * d
    out *= Fraction(d, kappa) ** d
    for degree in degrees:
        out *= factorial(degree - 1) ** 2
    out /= (2 * pi_upper) ** r
    return out


def stable_edge_margin(data: dict[str, object], pi_upper: int = 4) -> float:
    d = int(data["d"])
    c = int(data["c"])
    n = int(data["stable_last_odd"])
    alpha = alpha_value(data)
    return (
        log_a0_lower(data, pi_upper)
        + n * log(Fraction(d, c))
        - alpha * log(n)
        - log(8 * d * d)
    )


def leading_margin(data: dict[str, object], n: int, pi_upper: int = 4) -> float:
    d = int(data["d"])
    c = int(data["c"])
    alpha = alpha_value(data)
    return (
        log_a0_lower(data, pi_upper)
        + n * log(Fraction(d, c))
        - alpha * log(n)
        - log(8 * d * d)
    )


def post_stable_min_margin(data: dict[str, object]) -> tuple[int, float]:
    """Return the minimum leading margin over odd n >= stable edge.

    The continuous function n*L - alpha*log(n) is convex after its unique
    critical point alpha/L.  Hence it is enough to check the odd integers
    nearest that point, together with the stable edge.
    """

    d = int(data["d"])
    c = int(data["c"])
    alpha = alpha_value(data)
    edge = int(data["stable_last_odd"])
    critical = alpha / log(Fraction(d, c))
    candidates = {edge}
    center = int(critical)
    for n in range(max(edge, center - 5), center + 6):
        if n >= edge and n % 2 == 1:
            candidates.add(n)
    best_n = min(candidates, key=lambda n: leading_margin(data, n))
    return best_n, leading_margin(data, best_n)


def post_m29_local_onset_obstruction() -> None:
    """Replay the B2/C2 obstruction to half-leading local onset at n=63."""

    from character_ring import compute_moments
    from classical_chain_bridge import install_classical_cartan, q3

    n_value = 63
    rows = []
    for family in ("B", "C"):
        group = install_classical_cartan(family, 2)
        moments = compute_moments(group, n_value + 2, verbose=False)
        q_value = q3(moments, n_value)
        data = family_data(family, 2)
        d = int(data["d"])
        c = int(data["c"])
        alpha = alpha_value(data)
        a0_lower = a0_lower_fraction(data)
        half_leading = a0_lower * (2 * d) ** n_value / (
            2 * n_value**alpha
        )
        negative_allowance = 4 * d * d * (2 * c) ** n_value
        upper = Fraction(q_value + negative_allowance, 1)
        excess = half_leading - upper
        if excess <= 0:
            raise SystemExit(f"{family}_2 obstruction inequality did not replay")
        rows.append((family, q_value, a0_lower, half_leading, negative_allowance, excess))

    print("post-m29 local-onset obstruction")
    print("target odd n = 63")
    for family, q_value, a0_lower, half_leading, negative_allowance, excess in rows:
        print(f"{family}_2 Q3(63) = {q_value}")
        print(f"{family}_2 A_lower(pi<4) = {a0_lower}")
        print(f"{family}_2 half-leading lower at n=63 = {half_leading}")
        print(f"{family}_2 Proposition 21 negative allowance = {negative_allowance}")
        print(f"{family}_2 half-leading minus (Q3+allowance) = {excess}")


def b2_c2_post_m29_tail_certificate() -> None:
    """Replay the B2/C2 direct post-m29 pushforward tail certificate."""

    cap_width = Fraction(1, 4)
    cap_integral = cap_width**10 * (
        Fraction(1, 21) - Fraction(2, 25) + Fraction(1, 21)
    )
    cap_probability_lower = Fraction(1, 512) * Fraction(1, 4) ** 4 * cap_integral
    weighted_tail_lower = Fraction(64, 3) * cap_probability_lower
    monotone_target = 2 * Fraction(4, 7) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("B2/C2 post-m29 tail margin is not positive")

    print("B2/C2 post-m29 direct tail certificate")
    print("cap [0,1/4]^2; chi_adj >= 9 on the cap")
    print(f"cap polynomial integral = {cap_integral}")
    print(f"Pr(chi_adj >= 9) lower = {cap_probability_lower}")
    print(f"P_G(7) lower = {weighted_tail_lower}")
    print(f"2*(4/7)^63 = {monotone_target}")
    print(f"margin = {margin}")


def multiply_fraction_polynomials(
    left: dict[tuple[int, ...], Fraction],
    right: dict[tuple[int, ...], Fraction],
) -> dict[tuple[int, ...], Fraction]:
    out: dict[tuple[int, ...], Fraction] = {}
    for left_exp, left_coeff in left.items():
        for right_exp, right_coeff in right.items():
            exp = tuple(a + b for a, b in zip(left_exp, right_exp))
            out[exp] = out.get(exp, Fraction(0)) + left_coeff * right_coeff
    return {exp: coeff for exp, coeff in out.items() if coeff}


def linear_square_polynomial(coefficients: tuple[int, ...]) -> dict[tuple[int, ...], Fraction]:
    rank = len(coefficients)
    out: dict[tuple[int, ...], Fraction] = {}
    for i, coeff_i in enumerate(coefficients):
        for j, coeff_j in enumerate(coefficients):
            exp = [0] * rank
            exp[i] += 1
            exp[j] += 1
            exp_tuple = tuple(exp)
            out[exp_tuple] = out.get(exp_tuple, Fraction(0)) + coeff_i * coeff_j
    return out


def box_integral(poly: dict[tuple[int, ...], Fraction], width: Fraction) -> Fraction:
    total = Fraction(0)
    for exp, coeff in poly.items():
        term = coeff
        for power in exp:
            term *= width ** (power + 1) / Fraction(power + 1)
        total += term
    return total


def b3_c3_cap_integral(family: str) -> Fraction:
    roots: list[tuple[int, int, int]] = []
    if family == "B":
        roots.extend(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    elif family == "C":
        roots.extend(((2, 0, 0), (0, 2, 0), (0, 0, 2)))
    else:
        raise ValueError("family must be B or C")
    roots.extend(
        (
            (1, 1, 0),
            (1, -1, 0),
            (1, 0, 1),
            (1, 0, -1),
            (0, 1, 1),
            (0, 1, -1),
        )
    )
    poly: dict[tuple[int, ...], Fraction] = {(0, 0, 0): Fraction(1)}
    for root in roots:
        poly = multiply_fraction_polynomials(poly, linear_square_polynomial(root))
    return box_integral(poly, Fraction(1, 4))


def b3_c3_post_m29_tail_certificate() -> None:
    """Replay the B3/C3 direct post-m29 pushforward tail certificate."""

    print("B3/C3 post-m29 direct tail certificate")
    for family in ("B", "C"):
        cap_integral = b3_c3_cap_integral(family)
        cap_probability_lower = (
            Fraction(1, (2**3) * 6 * 8**3)
            * Fraction(1, 4) ** 9
            * cap_integral
        )
        weighted_tail_lower = Fraction(5041, 16) * cap_probability_lower * Fraction(1, 4)
        monotone_target = 2 * Fraction(2, 5) ** 63
        margin = weighted_tail_lower - monotone_target
        if margin <= 0:
            raise SystemExit(f"{family}3 post-m29 tail margin is not positive")
        print(f"{family}_3 cap polynomial integral = {cap_integral}")
        print(f"{family}_3 Pr(chi_adj >= 75/4) lower = {cap_probability_lower}")
        print(f"{family}_3 P_G(15) lower = {weighted_tail_lower}")
        print(f"{family}_3 2*(2/5)^63 = {monotone_target}")
        print(f"{family}_3 margin = {margin}")


def d4_post_m29_tail_certificate() -> None:
    """Replay the D4 Macdonald-Mehta ball-cap post-m29 tail certificate."""

    rank = 4
    r_plus = 12
    dimension = 28
    kappa = 6
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    c_value = 4
    shape = dimension // 2
    degrees = (2, 4, 4, 6)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    radius_s = Fraction(2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * Fraction(1, 4) ** r_plus
        * Fraction(1, 64)
        * Fraction(dimension, kappa) ** shape
        * radius_s**shape
        * degree_factor
        / (factorial(shape) * (2 * dimension) ** shape)
    )
    c_push = 22
    weighted_tail_lower = (
        (Fraction(dimension) - radius_s - 1) ** 2
        * ball_probability_lower
        * Fraction(1, c_value + 1)
    )
    monotone_target = 2 * Fraction(2 * c_value, c_push) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("D4 post-m29 tail margin is not positive")

    print("D4 post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 2; chi_adj >= 26")
    print(f"Pr(chi_adj >= 26) lower = {ball_probability_lower}")
    print(f"P_G(22) lower = {weighted_tail_lower}")
    print(f"2*(4/11)^63 = {monotone_target}")
    print(f"margin = {margin}")


def b4_c4_post_m29_tail_certificate() -> None:
    """Replay the root-normalized B4 Macdonald-Mehta ball-cap certificate."""

    rank = 4
    r_plus = 16
    dimension = 36
    weyl_order = 2**rank * factorial(rank)
    c_value = 4
    shape = dimension // 2
    degrees = (2, 4, 6, 8)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    radius_s = Fraction(9)
    c_push = 23
    monotone_target = 2 * Fraction(2 * c_value, c_push) ** 63

    print("B4 root-normalized post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 9; chi_adj >= 27")
    kappa = 7
    # Standard B_b Macdonald--Mehta roots have short roots sqrt(2)e_i.
    # The actual roots e_i multiply Delta(theta)^2 by 2^{-b}.
    root_normalization = Fraction(1, 2) ** rank
    ball_probability_lower = (
        root_normalization
        * Fraction(1, weyl_order)
        * Fraction(1, 4) ** r_plus
        * Fraction(1, 64)
        * Fraction(dimension, kappa) ** shape
        * radius_s**shape
        * degree_factor
        / (factorial(shape) * (2 * dimension) ** shape)
    )
    weighted_tail_lower = (
        (Fraction(dimension) - radius_s - 1) ** 2
        * ball_probability_lower
        * Fraction(1, c_value + 1)
    )
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("B4 post-m29 tail margin is not positive")
    print(f"B_4 root normalization = {root_normalization}")
    print(f"B_4 Pr(chi_adj >= 27) lower = {ball_probability_lower}")
    print(f"B_4 P_G(23) lower = {weighted_tail_lower}")
    print(f"B_4 2*(8/23)^63 = {monotone_target}")
    print(f"B_4 margin = {margin}")


def d5_post_m29_tail_certificate() -> None:
    """Replay the D5 Macdonald-Mehta ball-cap post-m29 tail certificate."""

    rank = 5
    r_plus = 20
    dimension = 45
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    c_value = 3
    degrees = (2, 4, 6, 8, 5)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)

    # The factor is a rational lower bound for
    # (2*pi)^(-5/2)/Gamma(47/2), using Gamma(47/2)=46!*sqrt(pi)/(4^23*23!),
    # pi < 4, and sqrt(2) < 3/2.
    half_gamma_factor = Fraction(4**23 * factorial(23), 384 * factorial(46))
    radius_s = Fraction(4)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * Fraction(1, 4) ** r_plus
        * degree_factor
        * half_gamma_factor
        * Fraction(1, 2) ** 45
    )
    c_push = 38
    weighted_tail_lower = (
        (Fraction(dimension) - radius_s - 1) ** 2
        * ball_probability_lower
        * Fraction(1, c_value + 1)
    )
    monotone_target = 2 * Fraction(2 * c_value, c_push) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("D5 post-m29 tail margin is not positive")

    print("D5 post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 4; chi_adj >= 41")
    print(f"Pr(chi_adj >= 41) lower = {ball_probability_lower}")
    print(f"P_G(38) lower = {weighted_tail_lower}")
    print(f"2*(3/19)^63 = {monotone_target}")
    print(f"margin = {margin}")


def c5_post_m29_tail_certificate() -> None:
    """Replay the C5 Macdonald-Mehta ball-cap post-m29 tail certificate."""

    rank = 5
    r_plus = 25
    dimension = 55
    kappa = 6
    weyl_order = 2**rank * factorial(rank)
    c_value = 5
    degrees = (2, 4, 6, 8, 10)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)

    # Rational lower bound for (2*pi)^(-5/2)/Gamma(57/2), using
    # Gamma(57/2)=56!*sqrt(pi)/(4^28*28!), pi < 4, and sqrt(2) < 3/2.
    half_gamma_factor = Fraction(4**28 * factorial(28), 384 * factorial(56))
    # (S/(2*kappa))^(55/2) at S=9 and kappa=6 is (3/4)^55 * sqrt(3/4).
    # Use sqrt(3)>5/3.
    radial_factor = Fraction(5 * 3**26, 2**55)
    radius_s = Fraction(9)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * Fraction(1, 4) ** r_plus
        * degree_factor
        * half_gamma_factor
        * radial_factor
    )
    c_push = 41
    weighted_tail_lower = (
        (Fraction(dimension) - radius_s - 1) ** 2
        * ball_probability_lower
        * Fraction(1, c_value + 1)
    )
    monotone_target = 2 * Fraction(2 * c_value, c_push) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("C5 post-m29 tail margin is not positive")

    print("C5 post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 9; chi_adj >= 46")
    print(f"Pr(chi_adj >= 46) lower = {ball_probability_lower}")
    print(f"P_G(41) lower = {weighted_tail_lower}")
    print(f"2*(10/41)^63 = {monotone_target}")
    print(f"margin = {margin}")


def b5_post_m29_tail_certificate() -> None:
    """Replay the root-normalized B5 sharpened ball-cap certificate."""

    rank = 5
    r_plus = 25
    dimension = 55
    kappa = 9
    weyl_order = 2**rank * factorial(rank)
    c_value = 5
    degrees = (2, 4, 6, 8, 10)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)

    # Rational lower bound for (2*pi)^(-5/2)/Gamma(57/2), using
    # Gamma(57/2)=56!*sqrt(pi)/(4^28*28!), pi < 22/7, and sqrt(2) < 3/2.
    half_gamma_factor = Fraction(4**28 * factorial(28), 1) / (
        Fraction(6) * Fraction(22, 7) ** 3 * factorial(56)
    )
    # At S=9 and kappa=9, (S/(2*kappa))^(55/2) = (1/2)^(55/2).
    # Use 1/sqrt(2)>2/3.
    radial_factor = Fraction(1, 2) ** 27 * Fraction(2, 3)
    sine_factor = Fraction(49, 121) ** r_plus
    radius_s = Fraction(9)
    root_normalization = Fraction(1, 2) ** rank
    ball_probability_lower = (
        root_normalization
        * Fraction(1, weyl_order)
        * sine_factor
        * degree_factor
        * half_gamma_factor
        * radial_factor
    )
    c_push = 41
    weighted_tail_lower = (
        (Fraction(dimension) - radius_s - 1) ** 2
        * ball_probability_lower
        * Fraction(1, c_value + 1)
    )
    monotone_target = 2 * Fraction(2 * c_value, c_push) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("B5 post-m29 tail margin is not positive")

    print("B5 root-normalized post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 9; chi_adj >= 46")
    print(f"B_5 root normalization = {root_normalization}")
    print(f"Pr(chi_adj >= 46) lower = {ball_probability_lower}")
    print(f"P_G(41) lower = {weighted_tail_lower}")
    print(f"2*(10/41)^63 = {monotone_target}")
    print(f"margin = {margin}")


def c6_post_m29_tail_certificate() -> None:
    """Replay the C6 central-Y post-m29 ball-cap tail certificate."""

    rank = 6
    r_plus = 36
    dimension = 78
    kappa = 7
    weyl_order = 2**rank * factorial(rank)
    degrees = (2, 4, 6, 8, 10, 12)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)

    # On |u|<=3, the alternating Taylor lower bound at |u|/2<=3/2 gives
    # 4 sin^2(u/2) >= (11/25) u^2.
    sine_factor = Fraction(11, 25) ** r_plus
    radius_s = Fraction(9)
    half_gamma_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** 3 * factorial(39)
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_factor
        * degree_factor
        * half_gamma_factor
        * Fraction(radius_s, 2 * kappa) ** 39
    )
    central_y_lower = Fraction(5, 9)
    c_push = Fraction(135, 2)
    weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
    monotone_target = 2 * Fraction(8, 45) ** 63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("C6 post-m29 tail margin is not positive")

    print("C6 post-m29 direct tail certificate")
    print("Macdonald-Mehta cap Q(theta) <= 9; chi_adj >= 69")
    print("central Y event |chi_adj| <= 3/2 has probability >= 5/9")
    print(f"Pr(chi_adj >= 69) lower = {ball_probability_lower}")
    print(f"P_G(135/2) lower = {weighted_tail_lower}")
    print(f"2*(8/45)^63 = {monotone_target}")
    print(f"margin = {margin}")


def d6_b6_post_m29_tail_certificate() -> None:
    """Replay the D6/B6 product-sine post-m29 ball-cap certificates."""

    rows = (
        {
            "family": "D",
            "rank": 6,
            "dimension": 66,
            "kappa": 10,
            "weyl_order": 2**5 * factorial(6),
            "c_value": 6,
            "degrees": (2, 4, 6, 8, 10, 6),
            "radius_s": Fraction(371, 20),
            "central_t": Fraction(7, 5),
            "c_push": Fraction(921, 20),
        },
        {
            "family": "B",
            "rank": 6,
            "dimension": 78,
            "kappa": 11,
            "weyl_order": 2**6 * factorial(6),
            "c_value": 6,
            "degrees": (2, 4, 6, 8, 10, 12),
            "radius_s": Fraction(1053, 50),
            "central_t": Fraction(3, 2),
            "c_push": Fraction(1386, 25),
        },
    )

    print("D6/B6 post-m29 product-sine direct tail certificate")
    for row in rows:
        family = str(row["family"])
        rank = int(row["rank"])
        dimension = int(row["dimension"])
        kappa = int(row["kappa"])
        weyl_order = int(row["weyl_order"])
        c_value = int(row["c_value"])
        degrees = row["degrees"]
        if not isinstance(degrees, tuple):
            raise TypeError("degrees must be a tuple")
        radius_s = Fraction(row["radius_s"])
        central_t = Fraction(row["central_t"])
        c_push = Fraction(row["c_push"])
        degree_factor = 1
        for degree in degrees:
            degree_factor *= factorial(degree)
        shape = dimension // 2
        half_gamma_factor = Fraction(1, 1) / (
            (2 * Fraction(22, 7)) ** (rank // 2) * factorial(shape)
        )
        sine_product_factor = (1 - radius_s / 24) ** 2
        ball_probability_lower = (
            Fraction(1, weyl_order)
            * sine_product_factor
            * degree_factor
            * half_gamma_factor
            * Fraction(radius_s, 2 * kappa) ** shape
        )
        if family == "B":
            ball_probability_lower *= Fraction(1, 2) ** rank
        central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)
        weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
        monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
        margin = weighted_tail_lower - monotone_target
        if margin <= 0:
            raise SystemExit(f"{family}6 post-m29 tail margin is not positive")

        print(f"{family}_6 cap Q(theta) <= {radius_s}; chi_adj >= {dimension - radius_s}")
        print(f"{family}_6 central Y event |chi_adj| <= {central_t}")
        print(f"{family}_6 sine product factor = {sine_product_factor}")
        if family == "B":
            print(f"{family}_6 root normalization = {Fraction(1, 2) ** rank}")
        print(f"{family}_6 Pr cap lower = {ball_probability_lower}")
        print(f"{family}_6 P_G({c_push}) lower = {weighted_tail_lower}")
        print(f"{family}_6 2*(2C/c)^63 = {monotone_target}")
        print(f"{family}_6 margin = {margin}")


def d7_b7_c7_post_m29_tail_certificate() -> None:
    """Replay the D7 product-sine post-m29 ball-cap certificate."""

    rows = (
        {
            "family": "D",
            "rank": 7,
            "dimension": 91,
            "kappa": 12,
            "weyl_order": 2**6 * factorial(7),
            "c_value": 5,
            "degrees": (2, 4, 6, 8, 10, 12, 7),
            "sqrt_q_lower": Fraction(433, 500),
            "c_push": Fraction(357, 5),
        },
    )

    radius_s = Fraction(18)
    central_t = Fraction(8, 5)
    sine_product_factor = (1 - radius_s / 24) ** 2
    central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)

    print("D7 post-m29 product-sine direct tail certificate")
    for row in rows:
        family = str(row["family"])
        rank = int(row["rank"])
        dimension = int(row["dimension"])
        kappa = int(row["kappa"])
        weyl_order = int(row["weyl_order"])
        c_value = int(row["c_value"])
        degrees = row["degrees"]
        if not isinstance(degrees, tuple):
            raise TypeError("degrees must be a tuple")
        sqrt_q_lower = Fraction(row["sqrt_q_lower"])
        c_push = Fraction(row["c_push"])
        degree_factor = 1
        for degree in degrees:
            degree_factor *= factorial(degree)
        half_index = (dimension + 1) // 2
        half_gamma_factor = (
            Fraction(4**half_index * factorial(half_index), 1)
            / (Fraction(12) * Fraction(22, 7) ** 4 * factorial(2 * half_index))
        )
        q_value = Fraction(radius_s, 2 * kappa)
        ball_probability_lower = (
            Fraction(1, weyl_order)
            * sine_product_factor
            * degree_factor
            * half_gamma_factor
            * q_value ** ((dimension - 1) // 2)
            * sqrt_q_lower
        )
        weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
        monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
        margin = weighted_tail_lower - monotone_target
        if margin <= 0:
            raise SystemExit(f"{family}7 post-m29 tail margin is not positive")

        print(f"{family}_7 cap Q(theta) <= {radius_s}; chi_adj >= {dimension - radius_s}")
        print(f"{family}_7 sqrt(S/(2*kappa)) lower = {sqrt_q_lower}")
        print(f"{family}_7 sine product factor = {sine_product_factor}")
        print(f"{family}_7 central Y event |chi_adj| <= {central_t}")
        print(f"{family}_7 Pr cap lower = {ball_probability_lower}")
        print(f"{family}_7 P_G({c_push}) lower = {weighted_tail_lower}")
        print(f"{family}_7 2*(2C/c)^63 = {monotone_target}")
        print(f"{family}_7 margin = {margin}")


def c8_post_m29_tail_certificate() -> None:
    """Replay the C8 product-sine post-m29 ball-cap certificate."""

    rank = 8
    dimension = 136
    kappa = 9
    weyl_order = 2**rank * factorial(rank)
    c_value = 8
    degrees = (2, 4, 6, 8, 10, 12, 14, 16)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    radius_s = Fraction(463, 20)
    central_t = Fraction(17, 10)
    c_push = Fraction(2223, 20)
    sine_product_factor = (1 - radius_s / 24) ** 2
    central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)
    half_gamma_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * half_gamma_factor
        * Fraction(radius_s, 2 * kappa) ** (dimension // 2)
    )
    weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
    monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("C8 post-m29 tail margin is not positive")

    print("C8 post-m29 product-sine direct tail certificate")
    print(f"C_8 cap Q(theta) <= {radius_s}; chi_adj >= {dimension - radius_s}")
    print(f"C_8 sine product factor = {sine_product_factor}")
    print(f"C_8 central Y event |chi_adj| <= {central_t}")
    print(f"C_8 Pr cap lower = {ball_probability_lower}")
    print(f"C_8 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_8 2*(2C/c)^63 = {monotone_target}")
    print(f"C_8 margin = {margin}")


def d8_post_m29_tail_certificate() -> None:
    """Replay the D8 radial-sine post-m29 ball-cap certificate."""

    rank = 8
    dimension = 120
    kappa = 14
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    c_value = 8
    degrees = (2, 4, 6, 8, 10, 12, 14, 8)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    radius_s = Fraction(36)
    central_t = Fraction(8, 5)
    c_push = Fraction(412, 5)
    # From pi > 333/106, sin(3)=sin(pi-3)>=sin(15/106), and
    # sin(x)>=x-x^3/6 for x>=0.
    sin_three_lower = Fraction(15, 106) - Fraction(15, 106) ** 3 / 6
    radial_sine_factor = sin_three_lower**2 / 9
    central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)
    half_gamma_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * radial_sine_factor
        * degree_factor
        * half_gamma_factor
        * Fraction(radius_s, 2 * kappa) ** (dimension // 2)
    )
    weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
    monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("D8 post-m29 tail margin is not positive")

    print("D8 post-m29 radial-sine direct tail certificate")
    print(f"D_8 cap Q(theta) <= {radius_s}; chi_adj >= {dimension - radius_s}")
    print(f"D_8 sin(3) lower = {sin_three_lower}")
    print(f"D_8 radial sine factor = {radial_sine_factor}")
    print(f"D_8 central Y event |chi_adj| <= {central_t}")
    print(f"D_8 Pr cap lower = {ball_probability_lower}")
    print(f"D_8 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_8 2*(2C/c)^63 = {monotone_target}")
    print(f"D_8 margin = {margin}")


def half_integer_density_factor(rank: int, dimension: int) -> Fraction:
    """Lower bound for (2*pi)^(-rank/2)/Gamma(dimension/2+1), odd parity."""

    if rank % 2 != 1 or dimension % 2 != 1:
        raise ValueError("half-integer density factor requires odd rank and dimension")
    half_index = (dimension + 1) // 2
    pi_power = (rank + 1) // 2
    return (
        Fraction(7, 5)
        * Fraction(4**half_index * factorial(half_index), 1)
        / ((2 * Fraction(22, 7)) ** pi_power * factorial(2 * half_index))
    )


def b8_c9_d9_post_m29_tail_certificate() -> None:
    """Replay the root-normalized B8 and type-D D9 cap certificates."""

    rows = (
        {
            "family": "B",
            "rank": 8,
            "dimension": 136,
            "kappa": 15,
            "root_length_sq": 2,
            "weyl_order": 2**8 * factorial(8),
            "c_value": 8,
            "degrees": (2, 4, 6, 8, 10, 12, 14, 16),
            "radius_s": Fraction(1609, 25),
            "central_t": Fraction(3, 2),
            "sqrt_q_lower": None,
        },
        {
            "family": "D",
            "rank": 9,
            "dimension": 153,
            "kappa": 16,
            "root_length_sq": 2,
            "weyl_order": 2**8 * factorial(9),
            "c_value": 7,
            "degrees": (2, 4, 6, 8, 10, 12, 14, 16, 9),
            "radius_s": Fraction(763, 10),
            "central_t": Fraction(8, 5),
            "sqrt_q_lower": Fraction(193, 125),
        },
    )

    print("B8/D9 post-m29 root-length direct tail certificate")
    for row in rows:
        family = str(row["family"])
        rank = int(row["rank"])
        dimension = int(row["dimension"])
        kappa = int(row["kappa"])
        root_length_sq = int(row["root_length_sq"])
        weyl_order = int(row["weyl_order"])
        c_value = int(row["c_value"])
        degrees = row["degrees"]
        if not isinstance(degrees, tuple):
            raise TypeError("degrees must be a tuple")
        radius_s = Fraction(row["radius_s"])
        central_t = Fraction(row["central_t"])
        c_push = Fraction(dimension) - radius_s - central_t
        exponent = Fraction(2 * kappa, root_length_sq)
        if exponent.denominator != 1:
            raise SystemExit(f"{family}{rank} has nonintegral root-length exponent")
        sine_product_factor = (
            1 - Fraction(root_length_sq, 1) * radius_s / (24 * kappa)
        ) ** exponent.numerator
        degree_factor = 1
        for degree in degrees:
            degree_factor *= factorial(degree)
        if dimension % 2 == 0:
            density_factor = Fraction(1, 1) / (
                (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
            )
            radial_factor = Fraction(radius_s, 2 * kappa) ** (dimension // 2)
        else:
            sqrt_q_lower = Fraction(row["sqrt_q_lower"])
            density_factor = half_integer_density_factor(rank, dimension)
            radial_factor = (
                Fraction(radius_s, 2 * kappa) ** ((dimension - 1) // 2)
                * sqrt_q_lower
            )
        ball_probability_lower = (
            Fraction(1, weyl_order)
            * sine_product_factor
            * degree_factor
            * density_factor
            * radial_factor
        )
        if family == "B":
            ball_probability_lower *= Fraction(1, 2) ** rank
        central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)
        weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
        monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
        margin = weighted_tail_lower - monotone_target
        if margin <= 0:
            raise SystemExit(f"{family}{rank} post-m29 tail margin is not positive")

        print(
            f"{family}_{rank} cap Q(theta) <= {radius_s}; "
            f"chi_adj >= {dimension - radius_s}"
        )
        print(f"{family}_{rank} root-length sine factor = {sine_product_factor}")
        if family == "B":
            print(
                f"{family}_{rank} root normalization = "
                f"{Fraction(1, 2) ** rank}"
            )
        if dimension % 2 == 1:
            print(f"{family}_{rank} sqrt(S/(2*kappa)) lower = {row['sqrt_q_lower']}")
        print(f"{family}_{rank} central Y event |chi_adj| <= {central_t}")
        print(f"{family}_{rank} Pr cap lower = {ball_probability_lower}")
        print(f"{family}_{rank} P_G({c_push}) lower = {weighted_tail_lower}")
        print(f"{family}_{rank} 2*(2C/c)^63 = {monotone_target}")
        print(f"{family}_{rank} margin = {margin}")


def c10_post_m29_tail_certificate() -> None:
    """Replay the C10 half-integer root-length post-m29 ball-cap certificate."""

    rank = 10
    dimension = 210
    kappa = 11
    weyl_order = 2**10 * factorial(10)
    c_value = 10
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(308, 5)
    central_t = Fraction(37, 20)
    sqrt_base_lower = Fraction(1, 4)
    c_push = Fraction(dimension) - radius_s - central_t
    root_length_base = 1 - radius_s / (6 * kappa)
    if sqrt_base_lower * sqrt_base_lower > root_length_base:
        raise SystemExit("C10 root-length square-root lower bound is invalid")
    sine_product_factor = root_length_base ** ((kappa - 1) // 2) * sqrt_base_lower
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    central_y_lower = 1 - Fraction(1, 1) / (central_t * central_t)
    weighted_tail_lower = c_push**2 * ball_probability_lower * central_y_lower
    monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("C10 post-m29 tail margin is not positive")

    print("C10 post-m29 half-integer root-length direct tail certificate")
    print(f"C_10 cap Q(theta) <= {radius_s}; chi_adj >= {dimension - radius_s}")
    print(f"C_10 root-length base = {root_length_base}; sqrt(base) lower = 1/4")
    print(f"C_10 root-length sine factor = {sine_product_factor}")
    print(f"C_10 central Y event |chi_adj| <= {central_t}")
    print(f"C_10 Pr cap lower = {ball_probability_lower}")
    print(f"C_10 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_10 2*(2C/c)^63 = {monotone_target}")
    print(f"C_10 margin = {margin}")


def exp_upper_bound(
    value: Fraction,
    terms: int,
    progress_label: str | None = None,
) -> Fraction:
    """Rational upper bound for exp(value) from a Taylor sum and tail ratio."""

    if terms < 0:
        raise ValueError("terms must be nonnegative")
    progress_step = max(1, terms // 20) if terms else 1
    term = Fraction(1, 1)
    total = term
    if progress_label is not None:
        print(f"{progress_label} exp_upper 0/{terms}", flush=True)
    for index in range(1, terms + 1):
        term *= value / index
        total += term
        if progress_label is not None and (
            index == terms or index % progress_step == 0
        ):
            print(f"{progress_label} exp_upper {index}/{terms}", flush=True)
    next_term = term * value / (terms + 1)
    ratio_bound = value / (terms + 2)
    if ratio_bound >= 1:
        raise ValueError("tail ratio bound must be below one")
    return total + next_term / (1 - ratio_bound)


def b9_post_m29_tail_certificate() -> None:
    """Replay the B9 exponential-sine/quartic-tail post-m29 certificate."""

    rank = 9
    dimension = 171
    kappa = 17
    weyl_order = 2**9 * factorial(9)
    c_value = 9
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(92, 1)
    central_t = Fraction(59, 50)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    # Euler-product sine lower bound with pi^2 > 39/4 and exp(A) < 3155.
    q_bound = Fraction(2) * radius_s / (39 * kappa)
    if q_bound >= 1:
        raise SystemExit(f"{label} sine-product cap radius is outside the checked range")
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 80)
    if exp_upper >= 3155:
        raise SystemExit("B9 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 3155)

    # Quartic majorant for E[(X-Y)^2 1_{Y < -T}] at X >= x_floor.
    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    if pushforward_weight_lower <= 0:
        raise SystemExit("B9 quartic pushforward weight lower bound is not positive")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(16449, 10000)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("B9 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (radius_s / (2 * kappa)) ** ((dimension - 1) // 2) * sqrt_radial_lower
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower
    monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("B9 post-m29 tail margin is not positive")

    print("B9 post-m29 exponential-sine quartic-tail certificate")
    print(f"B_9 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"B_9 sine exponent upper A = {sine_exponent_bound}")
    print(f"B_9 exp(A) upper = {exp_upper} < 3155")
    print(f"B_9 sine product factor = {sine_product_factor}")
    print(f"B_9 radial sqrt lower = {sqrt_radial_lower}")
    print(f"B_9 quartic moment = {quartic_moment}")
    print(f"B_9 quartic tail upper = {tail_weight_upper}")
    print(f"B_9 pushforward weight lower = {pushforward_weight_lower}")
    print(f"B_9 Pr cap lower = {ball_probability_lower}")
    print(f"B_9 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"B_9 2*(2C/c)^63 = {monotone_target}")
    print(f"B_9 margin = {margin}")


def c11_post_m29_tail_certificate() -> None:
    """Replay the C11 exponential-sine/quartic-tail post-m29 certificate."""

    rank = 11
    dimension = 253
    kappa = 12
    weyl_order = 2**11 * factorial(11)
    c_value = 11
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(96, 1)
    central_t = Fraction(36, 25)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    # Euler-product sine lower bound with pi^2 > 39/4 and exp(A) < 210000.
    q_bound = radius_s / 117
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (4 * 1440)
        + Fraction(7) * radius_s**3 / (144 * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 90)
    if exp_upper >= 210000:
        raise SystemExit("C11 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 210000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    if pushforward_weight_lower <= 0:
        raise SystemExit("C11 quartic pushforward weight lower bound is not positive")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(2, 1)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("C11 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (radius_s / (2 * kappa)) ** ((dimension - 1) // 2) * sqrt_radial_lower
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower
    monotone_target = 2 * Fraction(2 * c_value, 1) ** 63 / c_push**63
    margin = weighted_tail_lower - monotone_target
    if margin <= 0:
        raise SystemExit("C11 post-m29 tail margin is not positive")

    print("C11 post-m29 exponential-sine quartic-tail certificate")
    print(f"C_11 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_11 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_11 exp(A) upper = {exp_upper} < 210000")
    print(f"C_11 sine product factor = {sine_product_factor}")
    print(f"C_11 radial sqrt lower = {sqrt_radial_lower}")
    print(f"C_11 quartic moment = {quartic_moment}")
    print(f"C_11 quartic tail upper = {tail_weight_upper}")
    print(f"C_11 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_11 Pr cap lower = {ball_probability_lower}")
    print(f"C_11 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_11 2*(2C/c)^63 = {monotone_target}")
    print(f"C_11 margin = {margin}")


def c12_post_m29_tail_certificate() -> None:
    """Replay the C12 sharpened negative-tail post-m29 certificate."""

    rank = 12
    dimension = 300
    kappa = 13
    weyl_order = 2**rank * factorial(rank)
    c_value = 12
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(461, 4)
    central_t = Fraction(31, 20)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    # Euler-product sine lower bound with pi^2 > 39/4 and exp(A) < 172000000.
    q_bound = 4 * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + Fraction(38) * radius_s * radius_s / (1440 * kappa * kappa)
        + Fraction(86)
        * radius_s**3
        / (90720 * kappa**3 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 130)
    if exp_upper >= 172000000:
        raise SystemExit("C12 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 172000000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C12 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C12 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    stable_moments = read_stable_moments(18)
    q3_even_16 = q3_integer(stable_moments, 16)
    if q3_even_16 != 1960653454347808:
        raise SystemExit("unexpected stable Q3(16) value")
    negative_cutoff = Fraction(96, 5)
    negative_tail_mass_upper = Fraction(q3_even_16, 1) / negative_cutoff**16
    sharpened_target = (
        negative_tail_mass_upper * Fraction(2 * c_value, 1) ** 63 / c_push**63
        + 2 * negative_cutoff**63 / c_push**63
    )
    margin = weighted_tail_lower - sharpened_target
    if margin <= 0:
        raise SystemExit("C12 post-m29 sharpened tail margin is not positive")

    print("C12 post-m29 sharpened negative-tail certificate")
    print(f"C_12 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_12 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_12 exp(A) upper = {exp_upper} < 172000000")
    print(f"C_12 sine product factor = {sine_product_factor}")
    print(f"C_12 quartic moment = {quartic_moment}")
    print(f"C_12 quartic tail upper = {tail_weight_upper}")
    print(f"C_12 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_12 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_12 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_12 Pr cap lower = {ball_probability_lower}")
    print(f"C_12 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_12 Q3(16) = {q3_even_16}")
    print(f"C_12 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_12 sharpened target = {sharpened_target}")
    print(f"C_12 margin = {margin}")


def b10_post_m29_tail_certificate(delta_logs: list[Path] | None = None) -> None:
    """Replay the B10 post-m29 Chain bridge and sharpened tail certificate."""

    family = "B"
    rank = 10
    first_boundary = 2 * rank + 2
    exact_through = max(B10_POST_M29_DELTAS)
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [
        index for index in expected_indices if index not in B10_POST_M29_DELTAS
    ]
    if missing_indices:
        raise SystemExit(f"B10 embedded deltas have gaps: {missing_indices}")
    positive_deltas = {
        index: value for index, value in B10_POST_M29_DELTAS.items() if value > 0
    }
    if positive_deltas:
        raise SystemExit(f"B10 embedded deltas are positive: {positive_deltas}")
    if delta_logs:
        parsed = parse_bc_delta_logs(delta_logs)
        source_deltas = parsed.get((family, rank), {})
        if source_deltas != B10_POST_M29_DELTAS:
            raise SystemExit("B10 embedded deltas differ from the supplied C++ logs")

    print("B10 post-m29 Chain bridge and sharpened negative-tail certificate")
    print(f"B_10 embedded exact deltas through m_{exact_through}")
    if delta_logs:
        print("B_10 embedded deltas matched supplied C++ logs:")
        for path in delta_logs:
            print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m, stable_quadratic_tail in ((30, False), (31, False), (32, True)):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0
        ]
        active_negative_quadratic = [
            (first, second, value)
            for first, second, value in negative_quadratic
            if first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside and not stable_quadratic_tail:
            raise SystemExit(
                f"B10 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in B10_POST_M29_DELTAS.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * B10_POST_M29_DELTAS[first] * B10_POST_M29_DELTAS[second]
            for first, second, value in active_negative_quadratic
            if first in B10_POST_M29_DELTAS and second in B10_POST_M29_DELTAS
        )
        negative_quadratic_tail_bound = sum(
            value * stable[first] * stable[second]
            for first, second, value in active_negative_quadratic
            if first not in B10_POST_M29_DELTAS or second not in B10_POST_M29_DELTAS
        )
        if not stable_quadratic_tail:
            negative_quadratic_tail_bound = 0
        positive_quadratic_exact = sum(
            value * B10_POST_M29_DELTAS[first] * B10_POST_M29_DELTAS[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0
            and first in B10_POST_M29_DELTAS
            and second in B10_POST_M29_DELTAS
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + negative_quadratic_tail_bound
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"B10 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        if stable_quadratic_tail:
            print(
                f"B_10 m={chain_m}: stable D={stable_diff}; "
                f"tail={negative_linear_bound}; linear exact={linear_exact}; "
                f"negative quadratic exact={negative_quadratic_exact}; "
                f"negative quadratic tail bound={negative_quadratic_tail_bound}; "
                f"outside negative quadratic pairs={outside}; "
                f"positive quadratic exact={positive_quadratic_exact}; "
                f"lower margin={lower_margin}"
            )
        else:
            print(
                f"B_10 m={chain_m}: stable D={stable_diff}; "
                f"tail={negative_linear_bound}; linear exact={linear_exact}; "
                f"negative quadratic exact={negative_quadratic_exact}; "
                f"positive quadratic exact={positive_quadratic_exact}; "
                f"lower margin={lower_margin}"
            )

    min_chain_m, min_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"B_10 minimum Chain bridge margin m={min_chain_m}: {min_chain_margin}")

    dimension = 210
    kappa = 19
    weyl_order = 2**rank * factorial(rank)
    c_value = 10
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(481, 4)
    central_t = Fraction(6, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 100)
    if exp_upper >= 55000:
        raise SystemExit("B10 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 55000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("B10 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("B10 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    stable_moments = read_stable_moments(16)
    q3_even_14 = q3_integer(stable_moments, 14)
    if q3_even_14 != 6424471800784:
        raise SystemExit("unexpected stable Q3(14) value")
    negative_cutoff = Fraction(417, 25)
    negative_tail_mass_upper = Fraction(q3_even_14, 1) / negative_cutoff**14
    target_exponent = 67
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    margin = weighted_tail_lower - sharpened_target
    if margin <= 0:
        raise SystemExit("B10 post-m29 sharpened tail margin is not positive")

    print(f"B_10 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"B_10 sine exponent upper A = {sine_exponent_bound}")
    print(f"B_10 exp(A) upper = {exp_upper} < 55000")
    print(f"B_10 sine product factor = {sine_product_factor}")
    print(f"B_10 quartic moment = {quartic_moment}")
    print(f"B_10 quartic tail upper = {tail_weight_upper}")
    print(f"B_10 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"B_10 pushforward second derivative = {pushforward_second_derivative}")
    print(f"B_10 pushforward weight lower = {pushforward_weight_lower}")
    print(f"B_10 Pr cap lower = {ball_probability_lower}")
    print(f"B_10 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"B_10 Q3(14) = {q3_even_14}")
    print(f"B_10 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"B_10 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"B_10 margin = {margin}")


def b11_post_m29_direct_tail_certificate() -> None:
    """Replay the B11 sharpened negative-tail certificate from n=77 onward."""

    rank = 11
    dimension = 253
    kappa = 21
    weyl_order = 2**rank * factorial(rank)
    c_value = 11
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(1402939, 10000)
    central_t = Fraction(6, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 160)
    if exp_upper >= 254491:
        raise SystemExit("B11 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 254491)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("B11 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("B11 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(4569, 2500)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("B11 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    stable_moments = read_stable_moments(16)
    q3_even_14 = q3_integer(stable_moments, 14)
    if q3_even_14 != 6424471800784:
        raise SystemExit("unexpected stable Q3(14) value")
    negative_cutoff = Fraction(37, 2)
    negative_tail_mass_upper = Fraction(q3_even_14, 1) / negative_cutoff**14
    target_exponent = 77
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    margin = weighted_tail_lower - sharpened_target
    if margin <= 0:
        raise SystemExit("B11 post-m29 sharpened tail margin is not positive")

    print("B11 post-m29 direct tail certificate")
    print(f"B_11 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"B_11 sine exponent upper A = {sine_exponent_bound}")
    print(f"B_11 exp(A) upper = {exp_upper} < 254491")
    print(f"B_11 sine product factor = {sine_product_factor}")
    print(f"B_11 quartic moment = {quartic_moment}")
    print(f"B_11 quartic tail upper = {tail_weight_upper}")
    print(f"B_11 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"B_11 pushforward second derivative = {pushforward_second_derivative}")
    print(f"B_11 pushforward weight lower = {pushforward_weight_lower}")
    print(f"B_11 radial square-root lower = {sqrt_radial_lower}")
    print(f"B_11 Pr cap lower = {ball_probability_lower}")
    print(f"B_11 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"B_11 Q3(14) = {q3_even_14}")
    print(f"B_11 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"B_11 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"B_11 margin = {margin}")


def low_bc_exact_onset_bridge_certificate(delta_logs: list[Path]) -> None:
    """Replay the finite Chain gaps preceding the fixed exact-MGF onsets."""

    if not delta_logs:
        raise SystemExit(
            "low-B/C exact-onset bridge replay needs --bc-window-delta-log inputs"
        )
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    parsed = parse_bc_delta_logs(delta_logs)
    # The exact-MGF tail begins at the odd onset in the final column.  Starting
    # from the independently checked odd exponent 63, Chain steps m=31,...,
    # (onset-5)/2 cover precisely the intervening odd exponents.
    rows = (
        ("B", 11, 51, 75),
        ("B", 12, 53, 81),
        ("B", 13, 57, 83),
        ("C", 13, 43, 71),
        ("C", 14, 51, 75),
        ("C", 17, 45, 73),
    )
    total_steps = 0
    global_minimum: tuple[str, int, int] | None = None

    print("Low-rank B/C exact-MGF onset bridge certificate")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")

    for family, rank, exact_through, onset in rows:
        first_boundary = 2 * rank + 2
        deltas = parsed.get((family, rank), {})
        expected_indices = list(range(first_boundary, exact_through + 1))
        missing_indices = [
            index for index in expected_indices if index not in deltas
        ]
        if missing_indices:
            raise SystemExit(
                f"{family}{rank} supplied deltas have gaps: {missing_indices}"
            )
        positive_deltas = {
            index: value for index, value in deltas.items() if value > 0
        }
        if positive_deltas:
            raise SystemExit(
                f"{family}{rank} supplied deltas are positive: {positive_deltas}"
            )

        last_chain_m = (onset - 5) // 2
        row_margins: list[tuple[int, int]] = []
        print(
            f"{family}_{rank} exact deltas through m_{exact_through}; "
            f"exact-MGF onset n={onset}"
        )
        for chain_m in range(31, last_chain_m + 1):
            moment_max = 2 * chain_m + 5
            stable = read_stable_moments(moment_max)
            stable_diff = chain_diff_integer(stable, chain_m)
            linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
            quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
            tail_start = exact_through + 1
            negative_linear_bound = -sum(
                linear_coeffs[index] * stable[index]
                for index in range(tail_start, moment_max + 1)
                if linear_coeffs[index] > 0
            )
            linear_exact = sum(
                linear_coeffs[index] * value
                for index, value in deltas.items()
                if index <= moment_max
            )

            negative_quadratic_exact = 0
            negative_quadratic_stable_bound = 0
            stable_bound_pairs: list[tuple[int, int]] = []
            for (first, second), coefficient in sorted(quadratic_coeffs.items()):
                if coefficient >= 0 or first < first_boundary or second < first_boundary:
                    continue
                if first in deltas and second in deltas:
                    negative_quadratic_exact += (
                        coefficient * deltas[first] * deltas[second]
                    )
                else:
                    first_bound = -deltas[first] if first in deltas else stable[first]
                    second_bound = (
                        -deltas[second] if second in deltas else stable[second]
                    )
                    negative_quadratic_stable_bound += (
                        coefficient * first_bound * second_bound
                    )
                    stable_bound_pairs.append((first, second))

            positive_quadratic_exact = sum(
                coefficient * deltas[first] * deltas[second]
                for (first, second), coefficient in quadratic_coeffs.items()
                if coefficient > 0 and first in deltas and second in deltas
            )
            lower_margin = (
                stable_diff
                + negative_linear_bound
                + linear_exact
                + negative_quadratic_exact
                + negative_quadratic_stable_bound
                + positive_quadratic_exact
            )
            if lower_margin <= 0:
                raise SystemExit(
                    f"{family}{rank} m={chain_m} lower margin is not positive"
                )
            row_margins.append((chain_m, lower_margin))
            total_steps += 1
            print(
                f"{family}_{rank} m={chain_m}: stable D={stable_diff}; "
                f"tail={negative_linear_bound}; linear exact={linear_exact}; "
                f"negative quadratic exact={negative_quadratic_exact}; "
                f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
                f"positive quadratic exact={positive_quadratic_exact}; "
                f"stable-bound pairs={stable_bound_pairs}; "
                f"lower margin={lower_margin}"
            )

        if row_margins:
            row_minimum_m, row_minimum = min(row_margins, key=lambda item: item[1])
            print(
                f"{family}_{rank} minimum bridge margin "
                f"m={row_minimum_m}: {row_minimum}"
            )
            candidate = (f"{family}_{rank}", row_minimum_m, row_minimum)
            if global_minimum is None or candidate[2] < global_minimum[2]:
                global_minimum = candidate

    if global_minimum is None:
        raise SystemExit("low-B/C exact-onset bridge checked no Chain steps")
    print(
        f"minimum bridge row={global_minimum[0]} m={global_minimum[1]} "
        f"lower_margin={global_minimum[2]}"
    )
    print(
        f"SUMMARY rows_checked={len(rows)} chain_steps_checked={total_steps} "
        "failures=0"
    )
    print("__EXIT_STATUS=0")


def b11_post_m29_bridge_range_certificate(delta_logs: list[Path]) -> None:
    """Replay the B11 finite Chain bridge from n=65 through n=75."""

    if not delta_logs:
        raise SystemExit("B11 bridge replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "B"
    rank = 11
    first_boundary = 2 * rank + 2
    exact_through = 51
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"B11 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"B11 supplied deltas are positive: {positive_deltas}")

    print("B11 post-m29 finite Chain bridge certificate")
    print(f"B_11 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 37):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0
        ]
        active_negative_quadratic = [
            (first, second, value)
            for first, second, value in negative_quadratic
            if first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"B11 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"B11 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"B_11 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_margin = min(chain_margins, key=lambda item: item[1])
    print(f"B_11 minimum bridge margin m={minimum_chain_m}: {minimum_margin}")


def b12_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the B12 finite bridge and Chebyshev tail from n=81 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("B12 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "B"
    rank = 12
    first_boundary = 2 * rank + 2
    exact_through = 53
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"B12 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"B12 supplied deltas are positive: {positive_deltas}")

    print("B12 post-m29 finite-tail certificate")
    print(f"B_12 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 39):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"B12 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"B12 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"B_12 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"B_12 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 300
    kappa = 23
    c_value = 12
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(367, 2)
    central_t = Fraction(6, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 120)
    if exp_upper >= 18000000:
        raise SystemExit("B12 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 18000000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("B12 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("B12 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    stable_moments = read_stable_moments(42)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_degree = 10
    chebyshev_scale = Fraction(8500, 1)
    negative_cutoff = Fraction(35, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("B12 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 81
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("B12 post-m29 sharpened tail margin is not positive")

    print(f"B_12 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"B_12 sine exponent upper A = {sine_exponent_bound}")
    print(f"B_12 exp(A) upper = {exp_upper} < 18000000")
    print(f"B_12 sine product factor = {sine_product_factor}")
    print(f"B_12 quartic moment = {quartic_moment}")
    print(f"B_12 quartic tail upper = {tail_weight_upper}")
    print(f"B_12 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"B_12 pushforward second derivative = {pushforward_second_derivative}")
    print(f"B_12 pushforward weight lower = {pushforward_weight_lower}")
    print(f"B_12 Pr cap lower = {ball_probability_lower}")
    print(f"B_12 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"B_12 Chebyshev degree = {chebyshev_degree}")
    print(f"B_12 Chebyshev scale = {chebyshev_scale}")
    print(f"B_12 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"B_12 Chebyshev denominator = {denominator}")
    print(f"B_12 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"B_12 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"B_12 tail margin = {tail_margin}")


def b13_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the B13 finite bridge and Chebyshev tail from n=91 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("B13 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "B"
    rank = 13
    first_boundary = 2 * rank + 2
    exact_through = 57
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"B13 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"B13 supplied deltas are positive: {positive_deltas}")

    print("B13 post-m29 finite-tail certificate")
    print(f"B_13 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 44):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in active_negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_abs_bound = -deltas[first] if first in deltas else stable[first]
                second_abs_bound = -deltas[second] if second in deltas else stable[second]
                negative_quadratic_stable_bound += (
                    value * first_abs_bound * second_abs_bound
                )
                negative_quadratic_stable_pairs.append((first, second))
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"B13 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"B_13 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"B_13 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 351
    kappa = 25
    c_value = 13
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(1089, 5)
    central_t = Fraction(5, 4)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(391, 100)
    quartic_s = Fraction(37, 100)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("B13 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("B13 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(417, 200)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("B13 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 14
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(14000, 1)
    negative_cutoff = Fraction(37, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("B13 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 91
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("B13 post-m29 sharpened tail margin is not positive")

    print(f"B_13 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"B_13 sine exponent upper A = {sine_exponent_bound}")
    print(f"B_13 exp(A) upper = {exp_upper}")
    print(f"B_13 sine product factor = {sine_product_factor}")
    print(f"B_13 quartic moment = {quartic_moment}")
    print(f"B_13 quartic tail upper = {tail_weight_upper}")
    print(f"B_13 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"B_13 pushforward second derivative = {pushforward_second_derivative}")
    print(f"B_13 pushforward weight lower = {pushforward_weight_lower}")
    print(f"B_13 radial square-root lower = {sqrt_radial_lower}")
    print(f"B_13 Pr cap lower = {ball_probability_lower}")
    print(f"B_13 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"B_13 Chebyshev degree = {chebyshev_degree}")
    print(f"B_13 Chebyshev scale = {chebyshev_scale}")
    print(f"B_13 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"B_13 Chebyshev denominator = {denominator}")
    print(f"B_13 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"B_13 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"B_13 tail margin = {tail_margin}")


def c13_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C13 finite bridge and Chebyshev tail from n=71 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C13 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 13
    first_boundary = 2 * rank + 2
    exact_through = 43
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C13 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C13 supplied deltas are positive: {positive_deltas}")

    print("C13 post-m29 finite-tail certificate")
    print(f"C_13 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 34):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"C13 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"C13 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"C_13 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"C_13 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 351
    kappa = 14
    c_value = 13
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(125, 1)
    central_t = Fraction(8, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(4) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + Fraction(38) * radius_s * radius_s / (1440 * kappa * kappa)
        + Fraction(86)
        * radius_s**3
        / (90720 * kappa**3 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    if exp_upper >= 823713432:
        raise SystemExit("C13 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 823713432)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C13 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C13 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(2641, 1250)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("C13 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(9000, 1)
    negative_cutoff = Fraction(37, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C13 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 71
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C13 post-m29 sharpened tail margin is not positive")

    print(f"C_13 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_13 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_13 exp(A) upper = {exp_upper} < 823713432")
    print(f"C_13 sine product factor = {sine_product_factor}")
    print(f"C_13 quartic moment = {quartic_moment}")
    print(f"C_13 quartic tail upper = {tail_weight_upper}")
    print(f"C_13 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_13 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_13 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_13 radial square-root lower = {sqrt_radial_lower}")
    print(f"C_13 Pr cap lower = {ball_probability_lower}")
    print(f"C_13 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_13 Chebyshev degree = {chebyshev_degree}")
    print(f"C_13 Chebyshev scale = {chebyshev_scale}")
    print(f"C_13 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_13 Chebyshev denominator = {denominator}")
    print(f"C_13 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_13 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_13 tail margin = {tail_margin}")


def c14_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C14 finite bridge and Chebyshev tail from n=83 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C14 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 14
    first_boundary = 2 * rank + 2
    exact_through = 51
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C14 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C14 supplied deltas are positive: {positive_deltas}")

    print("C14 post-m29 finite-tail certificate")
    print(f"C_14 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 40):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"C14 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"C14 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"C_14 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"C_14 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 406
    kappa = 15
    c_value = 14
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(135, 1)
    central_t = Fraction(9, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(4) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + Fraction(38) * radius_s * radius_s / (1440 * kappa * kappa)
        + Fraction(86)
        * radius_s**3
        / (90720 * kappa**3 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    if exp_upper >= 5197422006:
        raise SystemExit("C14 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 5197422006)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C14 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C14 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(9000, 1)
    negative_cutoff = Fraction(41, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C14 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 83
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C14 post-m29 sharpened tail margin is not positive")

    print(f"C_14 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_14 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_14 exp(A) upper = {exp_upper} < 5197422006")
    print(f"C_14 sine product factor = {sine_product_factor}")
    print(f"C_14 quartic moment = {quartic_moment}")
    print(f"C_14 quartic tail upper = {tail_weight_upper}")
    print(f"C_14 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_14 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_14 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_14 Pr cap lower = {ball_probability_lower}")
    print(f"C_14 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_14 Chebyshev degree = {chebyshev_degree}")
    print(f"C_14 Chebyshev scale = {chebyshev_scale}")
    print(f"C_14 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_14 Chebyshev denominator = {denominator}")
    print(f"C_14 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_14 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_14 tail margin = {tail_margin}")


def c15_post_m29_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C15 Chebyshev negative-tail post-m29 certificate."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C15 tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 15
    first_boundary = 2 * rank + 2
    exact_through = 39
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C15 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C15 supplied deltas are positive: {positive_deltas}")

    dimension = 465
    kappa = 16
    c_value = 15
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(280, 1)
    central_t = Fraction(8, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 200)
    if exp_upper >= 10_000_000_000_000_000:
        raise SystemExit("C15 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 10_000_000_000_000_000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C15 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C15 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(1479, 500)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("C15 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(11000, 1)
    negative_cutoff = Fraction(41, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C15 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 63
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C15 post-m29 sharpened tail margin is not positive")

    print("C15 post-m29 Chebyshev negative-tail certificate")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    print(f"C_15 supplied exact deltas through m_{exact_through}")
    print(f"C_15 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_15 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_15 exp(A) upper = {exp_upper} < 10000000000000000")
    print(f"C_15 sine product factor = {sine_product_factor}")
    print(f"C_15 quartic moment = {quartic_moment}")
    print(f"C_15 quartic tail upper = {tail_weight_upper}")
    print(f"C_15 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_15 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_15 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_15 radial square-root lower = {sqrt_radial_lower}")
    print(f"C_15 Pr cap lower = {ball_probability_lower}")
    print(f"C_15 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_15 Chebyshev degree = {chebyshev_degree}")
    print(f"C_15 Chebyshev scale = {chebyshev_scale}")
    print(f"C_15 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_15 Chebyshev denominator = {denominator}")
    print(f"C_15 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_15 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_15 tail margin = {tail_margin}")


def c16_post_m29_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C16 Chebyshev negative-tail post-m29 certificate."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C16 tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 16
    first_boundary = 2 * rank + 2
    exact_through = 39
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C16 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C16 supplied deltas are positive: {positive_deltas}")

    dimension = 528
    kappa = 17
    c_value = 16
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(302, 1)
    central_t = Fraction(8, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 120)
    if exp_upper >= 1_000_000_000_000_000_000:
        raise SystemExit("C16 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 1_000_000_000_000_000_000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C16 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C16 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(13000, 1)
    negative_cutoff = Fraction(41, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C16 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 63
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C16 post-m29 sharpened tail margin is not positive")

    print("C16 post-m29 Chebyshev negative-tail certificate")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    print(f"C_16 supplied exact deltas through m_{exact_through}")
    print(f"C_16 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_16 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_16 exp(A) upper = {exp_upper} < 1000000000000000000")
    print(f"C_16 sine product factor = {sine_product_factor}")
    print(f"C_16 quartic moment = {quartic_moment}")
    print(f"C_16 quartic tail upper = {tail_weight_upper}")
    print(f"C_16 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_16 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_16 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_16 Pr cap lower = {ball_probability_lower}")
    print(f"C_16 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_16 Chebyshev degree = {chebyshev_degree}")
    print(f"C_16 Chebyshev scale = {chebyshev_scale}")
    print(f"C_16 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_16 Chebyshev denominator = {denominator}")
    print(f"C_16 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_16 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_16 tail margin = {tail_margin}")


def c17_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C17 finite bridge and Chebyshev tail from n=73 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C17 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 17
    first_boundary = 2 * rank + 2
    exact_through = 45
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C17 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C17 supplied deltas are positive: {positive_deltas}")

    print("C17 post-m29 finite-tail certificate")
    print(f"C_17 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 35):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"C17 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"C17 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"C_17 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"C_17 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 595
    kappa = 18
    c_value = 17
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(322, 1)
    central_t = Fraction(17, 10)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    if exp_upper >= 10**20:
        raise SystemExit("C17 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 10**20)

    quartic_r = Fraction(1, 4)
    quartic_s = Fraction(3, 1)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C17 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C17 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(299, 100)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("C17 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 9
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(14500, 1)
    negative_cutoff = Fraction(45, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C17 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 73
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C17 post-m29 sharpened tail margin is not positive")

    print(f"C_17 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_17 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_17 exp(A) upper = {exp_upper} < 100000000000000000000")
    print(f"C_17 sine product factor = {sine_product_factor}")
    print(f"C_17 quartic moment = {quartic_moment}")
    print(f"C_17 quartic tail upper = {tail_weight_upper}")
    print(f"C_17 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_17 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_17 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_17 radial square-root lower = {sqrt_radial_lower}")
    print(f"C_17 Pr cap lower = {ball_probability_lower}")
    print(f"C_17 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_17 Chebyshev degree = {chebyshev_degree}")
    print(f"C_17 Chebyshev scale = {chebyshev_scale}")
    print(f"C_17 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_17 Chebyshev denominator = {denominator}")
    print(f"C_17 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_17 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_17 tail margin = {tail_margin}")


def c18_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C18 finite bridge and Chebyshev tail from n=85 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C18 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 18
    first_boundary = 2 * rank + 2
    exact_through = 50
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C18 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C18 supplied deltas are positive: {positive_deltas}")

    print("C18 post-m29 finite-tail certificate")
    print(f"C_18 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 41):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"C18 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"C18 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"C_18 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"C_18 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 666
    kappa = 19
    c_value = 18
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(340, 1)
    central_t = Fraction(17, 10)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(14, 5)
    quartic_s = Fraction(1, 20)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C18 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C18 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 9
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(18000, 1)
    negative_cutoff = Fraction(25, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C18 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 85
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C18 post-m29 sharpened tail margin is not positive")

    print(f"C_18 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_18 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_18 exp(A) upper = {exp_upper}")
    print(f"C_18 sine product factor = {sine_product_factor}")
    print(f"C_18 quartic moment = {quartic_moment}")
    print(f"C_18 quartic tail upper = {tail_weight_upper}")
    print(f"C_18 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_18 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_18 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_18 Pr cap lower = {ball_probability_lower}")
    print(f"C_18 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_18 Chebyshev degree = {chebyshev_degree}")
    print(f"C_18 Chebyshev scale = {chebyshev_scale}")
    print(f"C_18 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_18 Chebyshev denominator = {denominator}")
    print(f"C_18 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_18 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_18 tail margin = {tail_margin}")


def c19_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the C19 finite bridge and Chebyshev tail from n=99 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("C19 finite-tail replay needs --bc-window-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    family = "C"
    rank = 19
    first_boundary = 2 * rank + 2
    exact_through = 58
    parsed = parse_bc_delta_logs(delta_logs)
    deltas = parsed.get((family, rank), {})
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"C19 supplied deltas have gaps: {missing_indices}")
    positive_deltas = {index: value for index, value in deltas.items() if value > 0}
    if positive_deltas:
        raise SystemExit(f"C19 supplied deltas are positive: {positive_deltas}")

    print("C19 post-m29 finite-tail certificate")
    print(f"C_19 supplied exact deltas through m_{exact_through}")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 48):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        active_negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]
        outside = [
            (first, second)
            for first, second, _value in active_negative_quadratic
            if first > exact_through or second > exact_through
        ]
        if outside:
            raise SystemExit(
                f"C19 m={chain_m} has active negative quadratic terms outside "
                f"m_{exact_through}: {outside}"
            )
        tail_start = exact_through + 1
        negative_linear_bound = -sum(
            linear_coeffs[index] * stable[index]
            for index in range(tail_start, moment_max + 1)
            if linear_coeffs[index] > 0
        )
        linear_exact = sum(
            linear_coeffs[index] * value
            for index, value in deltas.items()
            if index <= moment_max
        )
        negative_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for first, second, value in active_negative_quadratic
        )
        positive_quadratic_exact = sum(
            value * deltas[first] * deltas[second]
            for (first, second), value in quadratic_coeffs.items()
            if value > 0 and first in deltas and second in deltas
        )
        lower_margin = (
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact
        )
        if lower_margin <= 0:
            raise SystemExit(f"C19 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"C_19 m={chain_m}: stable D={stable_diff}; "
            f"tail={negative_linear_bound}; linear exact={linear_exact}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"positive quadratic exact={positive_quadratic_exact}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"C_19 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 741
    kappa = 20
    c_value = 19
    weyl_order = 2**rank * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank + 1))
    radius_s = Fraction(359, 1)
    central_t = Fraction(17, 10)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(27, 10)
    quartic_s = Fraction(1, 8)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("C19 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("C19 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(299, 100)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("C19 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 10
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(21000, 1)
    negative_cutoff = Fraction(55, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("C19 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 99
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("C19 post-m29 sharpened tail margin is not positive")

    print(f"C_19 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"C_19 sine exponent upper A = {sine_exponent_bound}")
    print(f"C_19 exp(A) upper = {exp_upper}")
    print(f"C_19 sine product factor = {sine_product_factor}")
    print(f"C_19 quartic moment = {quartic_moment}")
    print(f"C_19 quartic tail upper = {tail_weight_upper}")
    print(f"C_19 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"C_19 pushforward second derivative = {pushforward_second_derivative}")
    print(f"C_19 pushforward weight lower = {pushforward_weight_lower}")
    print(f"C_19 radial square-root lower = {sqrt_radial_lower}")
    print(f"C_19 Pr cap lower = {ball_probability_lower}")
    print(f"C_19 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"C_19 Chebyshev degree = {chebyshev_degree}")
    print(f"C_19 Chebyshev scale = {chebyshev_scale}")
    print(f"C_19 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"C_19 Chebyshev denominator = {denominator}")
    print(f"C_19 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"C_19 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"C_19 tail margin = {tail_margin}")


def d10_post_m29_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D10 Chebyshev negative-tail post-m29 certificate."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    rank = 10
    dimension = 190
    kappa = 18
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    c_value = 10
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(107, 1)
    central_t = Fraction(6, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    required_deltas = list(range(rank, 51))
    missing = [index for index in required_deltas if index not in deltas]
    if missing:
        raise SystemExit(f"D10 missing exact deltas for tail certificate: {missing}")
    stable_moments = read_stable_moments(50)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value

    # Euler-product sine lower bound with pi^2 > 39/4 and exp(A) < 12400.
    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 90)
    if exp_upper >= 12400:
        raise SystemExit("D10 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 12400)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D10 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D10 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    upper_support = Fraction(2 * dimension, 1)
    chebyshev_degree = 12
    chebyshev_scale = Fraction(5925, 1)
    negative_cutoff = Fraction(27, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(49)]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D10 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 63
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    margin = weighted_tail_lower - sharpened_target
    if margin <= 0:
        raise SystemExit("D10 post-m29 sharpened tail margin is not positive")

    print("D10 post-m29 Chebyshev negative-tail certificate")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    print(f"D_10 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_10 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_10 exp(A) upper = {exp_upper} < 12400")
    print(f"D_10 sine product factor = {sine_product_factor}")
    print(f"D_10 quartic moment = {quartic_moment}")
    print(f"D_10 quartic tail upper = {tail_weight_upper}")
    print(f"D_10 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_10 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_10 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_10 Pr cap lower = {ball_probability_lower}")
    print(f"D_10 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_10 Chebyshev degree = {chebyshev_degree}")
    print(f"D_10 Chebyshev scale = {chebyshev_scale}")
    print(f"D_10 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_10 Chebyshev denominator = {denominator}")
    print(f"D_10 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_10 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_10 margin = {margin}")


def d11_post_m29_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D11 Chebyshev negative-tail post-m29 certificate."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    rank = 11
    dimension = 231
    kappa = 20
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    c_value = 9
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(140, 1)
    central_t = Fraction(6, 5)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    required_deltas = list(range(rank, 51))
    missing = [index for index in required_deltas if index not in deltas]
    if missing:
        raise SystemExit(f"D11 missing exact deltas for tail certificate: {missing}")
    stable_moments = read_stable_moments(52)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value

    # Euler-product sine lower bound with pi^2 > 39/4 and exp(A) < 260000.
    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 140)
    if exp_upper >= 260000:
        raise SystemExit("D11 exponential sine lower bound is not strong enough")
    sine_product_factor = Fraction(1, 260000)

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(15, 4)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D11 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D11 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(187, 100)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("D11 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    upper_support = Fraction(2 * dimension, 1)
    chebyshev_degree = 12
    chebyshev_scale = Fraction(7400, 1)
    negative_cutoff = Fraction(25, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(49)]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D11 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 63
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    margin = weighted_tail_lower - sharpened_target
    if margin <= 0:
        raise SystemExit("D11 post-m29 sharpened tail margin is not positive")

    print("D11 post-m29 Chebyshev negative-tail certificate")
    print("delta logs:")
    for path in delta_logs:
        print(f"  {path}")
    print(f"D_11 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_11 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_11 exp(A) upper = {exp_upper} < 260000")
    print(f"D_11 sine product factor = {sine_product_factor}")
    print(f"D_11 quartic moment = {quartic_moment}")
    print(f"D_11 quartic tail upper = {tail_weight_upper}")
    print(f"D_11 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_11 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_11 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_11 radial sqrt lower = {sqrt_radial_lower}")
    print(f"D_11 Pr cap lower = {ball_probability_lower}")
    print(f"D_11 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_11 Chebyshev degree = {chebyshev_degree}")
    print(f"D_11 Chebyshev scale = {chebyshev_scale}")
    print(f"D_11 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_11 Chebyshev denominator = {denominator}")
    print(f"D_11 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_11 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_11 margin = {margin}")


def d12_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D12 finite bridge and Chebyshev tail from n=79 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D12 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 12
    first_boundary = rank
    exact_through = 51
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D12 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D12 supplied determinant deltas are negative: {negative_deltas}")

    print("D12 post-m29 finite-tail certificate")
    print(f"D_12 supplied exact deltas through m_{exact_through}")
    print("D_12 bridge uses stable absolute-moment envelopes beyond m_51")
    for path in delta_logs:
        print(f"  {path}")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, 38):
        moment_max = 2 * chain_m + 5
        stable = read_stable_moments(moment_max)
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D12 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_12 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_12 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 276
    kappa = 22
    c_value = 12
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(164, 1)
    central_t = Fraction(4, 3)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 180)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(6, 5)
    quartic_s = Fraction(3, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D12 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D12 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 16
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(33, 2)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D12 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 79
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D12 post-m29 sharpened tail margin is not positive")

    print(f"D_12 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_12 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_12 exp(A) upper = {exp_upper}")
    print(f"D_12 sine product factor = {sine_product_factor}")
    print(f"D_12 quartic moment = {quartic_moment}")
    print(f"D_12 quartic tail upper = {tail_weight_upper}")
    print(f"D_12 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_12 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_12 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_12 Pr cap lower = {ball_probability_lower}")
    print(f"D_12 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_12 Chebyshev degree = {chebyshev_degree}")
    print(f"D_12 Chebyshev scale = {chebyshev_scale}")
    print(f"D_12 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_12 Chebyshev denominator = {denominator}")
    print(f"D_12 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_12 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_12 tail margin = {tail_margin}")


def d13_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D13 finite bridge and Chebyshev tail from n=201 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D13 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 13
    first_boundary = rank
    exact_through = 50
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D13 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D13 supplied determinant deltas are negative: {negative_deltas}")

    print("D13 post-m29 finite-tail certificate")
    print(f"D_13 supplied exact deltas through m_{exact_through}")
    print("D_13 bridge uses stable absolute-moment envelopes beyond m_50")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 98
    stable_all = stable_moments_formula(
        2 * bridge_max_chain + 5,
        progress_label=f"{row_label} stable moments",
    )
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D13 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D13 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_13 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_13 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 325
    kappa = 24
    c_value = 11
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(141, 1)
    central_t = Fraction(11, 10)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 220)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(19, 32)
    quartic_s = Fraction(7, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D13 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D13 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    sqrt_radial_lower = Fraction(17, 10)
    if sqrt_radial_lower * sqrt_radial_lower > radius_s / (2 * kappa):
        raise SystemExit("D13 radial square-root lower bound is invalid")
    density_factor = half_integer_density_factor(rank, dimension)
    radial_factor = (
        (radius_s / (2 * kappa)) ** ((dimension - 1) // 2)
        * sqrt_radial_lower
    )
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 12
    stable_moments = read_stable_moments(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(10000, 1)
    negative_cutoff = Fraction(25, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D13 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 201
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D13 post-m29 sharpened tail margin is not positive")

    print(f"D_13 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_13 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_13 exp(A) upper = {exp_upper}")
    print(f"D_13 sine product factor = {sine_product_factor}")
    print(f"D_13 quartic moment = {quartic_moment}")
    print(f"D_13 quartic tail upper = {tail_weight_upper}")
    print(f"D_13 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_13 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_13 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_13 radial sqrt lower = {sqrt_radial_lower}")
    print(f"D_13 Pr cap lower = {ball_probability_lower}")
    print(f"D_13 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_13 Chebyshev degree = {chebyshev_degree}")
    print(f"D_13 Chebyshev scale = {chebyshev_scale}")
    print(f"D_13 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_13 Chebyshev denominator = {denominator}")
    print(f"D_13 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_13 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_13 tail margin = {tail_margin}")


def d14_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D14 finite bridge and Chebyshev tail from n=325 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D14 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 14
    first_boundary = rank
    exact_through = 42
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D14 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D14 supplied determinant deltas are negative: {negative_deltas}")

    print("D14 post-m29 finite-tail certificate")
    print(f"D_14 supplied exact deltas through m_{exact_through}")
    print("D_14 bridge uses stable absolute-moment envelopes beyond m_42")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 160
    stable_all = stable_moments_formula(
        2 * bridge_max_chain + 5,
        progress_label=f"{row_label} stable moments",
    )
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D14 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D14 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_14 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_14 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 378
    kappa = 26
    c_value = 14
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(133, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 220)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(7, 10)
    quartic_s = Fraction(9, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D14 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D14 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 10
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(60, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D14 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 325
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D14 post-m29 sharpened tail margin is not positive")

    print(f"D_14 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_14 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_14 exp(A) upper = {exp_upper}")
    print(f"D_14 sine product factor = {sine_product_factor}")
    print(f"D_14 quartic moment = {quartic_moment}")
    print(f"D_14 quartic tail upper = {tail_weight_upper}")
    print(f"D_14 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_14 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_14 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_14 Pr cap lower = {ball_probability_lower}")
    print(f"D_14 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_14 Chebyshev degree = {chebyshev_degree}")
    print(f"D_14 Chebyshev scale = {chebyshev_scale}")
    print(f"D_14 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_14 Chebyshev denominator = {denominator}")
    print(f"D_14 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_14 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_14 tail margin = {tail_margin}")


def d15_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D15 finite bridge and Chebyshev tail from n=525 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D15 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 15
    first_boundary = rank
    exact_through = 48
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D15 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D15 supplied determinant deltas are negative: {negative_deltas}")

    print("D15 post-m29 finite-tail certificate")
    print(f"D_15 supplied exact deltas through m_{exact_through}")
    print("D_15 bridge uses stable absolute-moment envelopes beyond m_48")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 260
    stable_all = stable_moments_formula(
        2 * bridge_max_chain + 5,
        progress_label=f"{row_label} stable moments",
    )
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D15 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D15 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_15 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_15 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 435
    kappa = 28
    c_value = 15
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(131, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 240)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(1, 1)
    quartic_s = Fraction(13, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D15 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D15 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    half_index = (dimension + 1) // 2
    pi_power = (rank + 1) // 2
    density_factor = (
        Fraction(7, 5)
        * 4**half_index
        * factorial(half_index)
        / ((Fraction(44, 7) ** pi_power) * factorial(2 * half_index))
    )
    radial_factor = (radius_s / (2 * kappa)) ** half_index
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(10000, 1)
    negative_cutoff = Fraction(70, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D15 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 525
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D15 post-m29 sharpened tail margin is not positive")

    print(f"D_15 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_15 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_15 exp(A) upper = {exp_upper}")
    print(f"D_15 sine product factor = {sine_product_factor}")
    print(f"D_15 quartic moment = {quartic_moment}")
    print(f"D_15 quartic tail upper = {tail_weight_upper}")
    print(f"D_15 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_15 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_15 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_15 Pr cap lower = {ball_probability_lower}")
    print(f"D_15 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_15 Chebyshev degree = {chebyshev_degree}")
    print(f"D_15 Chebyshev scale = {chebyshev_scale}")
    print(f"D_15 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_15 Chebyshev denominator = {denominator}")
    print(f"D_15 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_15 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_15 tail margin = {tail_margin}")


def d16_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D16 finite bridge and Chebyshev tail from n=275 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D16 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 16
    first_boundary = rank
    exact_through = 45
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D16 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D16 supplied determinant deltas are negative: {negative_deltas}")

    print("D16 post-m29 finite-tail certificate")
    print(f"D_16 supplied exact deltas through m_{exact_through}")
    print("D_16 bridge uses stable absolute-moment envelopes beyond m_45")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 135
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D16 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D16 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_16 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_16 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 496
    kappa = 30
    c_value = 16
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(223, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 260)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(7, 10)
    quartic_s = Fraction(9, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D16 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D16 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(10000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D16 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 275
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D16 post-m29 sharpened tail margin is not positive")

    print(f"D_16 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_16 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_16 exp(A) upper = {exp_upper}")
    print(f"D_16 sine product factor = {sine_product_factor}")
    print(f"D_16 quartic moment = {quartic_moment}")
    print(f"D_16 quartic tail upper = {tail_weight_upper}")
    print(f"D_16 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_16 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_16 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_16 Pr cap lower = {ball_probability_lower}")
    print(f"D_16 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_16 Chebyshev degree = {chebyshev_degree}")
    print(f"D_16 Chebyshev scale = {chebyshev_scale}")
    print(f"D_16 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_16 Chebyshev denominator = {denominator}")
    print(f"D_16 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_16 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_16 tail margin = {tail_margin}")


def d17_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D17 finite bridge and Chebyshev tail from n=287 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D17 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 17
    first_boundary = rank
    exact_through = 44
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D17 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D17 supplied determinant deltas are negative: {negative_deltas}")

    print("D17 post-m29 finite-tail certificate")
    print(f"D_17 supplied exact deltas through m_{exact_through}")
    print("D_17 bridge uses stable absolute-moment envelopes beyond m_44")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 141
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D17 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D17 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_17 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_17 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 561
    kappa = 32
    c_value = 17
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(262, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 280)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(37, 10)
    quartic_s = Fraction(1, 2)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D17 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D17 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    half_index = (dimension + 1) // 2
    pi_power = (rank + 1) // 2
    density_factor = (
        Fraction(7, 5)
        * 4**half_index
        * factorial(half_index)
        / ((Fraction(44, 7) ** pi_power) * factorial(2 * half_index))
    )
    radial_factor = (radius_s / (2 * kappa)) ** half_index
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D17 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 287
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D17 post-m29 sharpened tail margin is not positive")

    print(f"D_17 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_17 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_17 exp(A) upper = {exp_upper}")
    print(f"D_17 sine product factor = {sine_product_factor}")
    print(f"D_17 quartic moment = {quartic_moment}")
    print(f"D_17 quartic tail upper = {tail_weight_upper}")
    print(f"D_17 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_17 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_17 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_17 Pr cap lower = {ball_probability_lower}")
    print(f"D_17 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_17 Chebyshev degree = {chebyshev_degree}")
    print(f"D_17 Chebyshev scale = {chebyshev_scale}")
    print(f"D_17 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_17 Chebyshev denominator = {denominator}")
    print(f"D_17 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_17 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_17 tail margin = {tail_margin}")


def d18_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D18 finite bridge and Chebyshev tail from n=365 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D18 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 18
    first_boundary = rank
    exact_through = 43
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D18 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D18 supplied determinant deltas are negative: {negative_deltas}")

    print("D18 post-m29 finite-tail certificate")
    print(f"D_18 supplied exact deltas through m_{exact_through}")
    print("D_18 bridge uses stable absolute-moment envelopes beyond m_43")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 180
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D18 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D18 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_18 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_18 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 630
    kappa = 34
    c_value = 18
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(291, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 320)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(4, 1)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D18 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D18 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (2 * Fraction(22, 7)) ** (rank // 2) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D18 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 365
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D18 post-m29 sharpened tail margin is not positive")

    print(f"D_18 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_18 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_18 exp(A) upper = {exp_upper}")
    print(f"D_18 sine product factor = {sine_product_factor}")
    print(f"D_18 quartic moment = {quartic_moment}")
    print(f"D_18 quartic tail upper = {tail_weight_upper}")
    print(f"D_18 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_18 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_18 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_18 Pr cap lower = {ball_probability_lower}")
    print(f"D_18 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_18 Chebyshev degree = {chebyshev_degree}")
    print(f"D_18 Chebyshev scale = {chebyshev_scale}")
    print(f"D_18 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_18 Chebyshev denominator = {denominator}")
    print(f"D_18 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_18 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_18 tail margin = {tail_margin}")


def d19_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D19 finite bridge and Chebyshev tail from n=405 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D19 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 19
    first_boundary = rank
    exact_through = 44
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D19 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D19 supplied determinant deltas are negative: {negative_deltas}")

    print("D19 post-m29 finite-tail certificate")
    print(f"D_19 supplied exact deltas through m_{exact_through}")
    print("D_19 bridge uses stable absolute-moment envelopes beyond m_44")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 200
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D19 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D19 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_19 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_19 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 703
    kappa = 36
    c_value = 19
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(321, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 360)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(18, 5)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D19 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D19 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    half_index = (dimension + 1) // 2
    pi_power = (rank + 1) // 2
    density_factor = (
        Fraction(7, 5)
        * 4**half_index
        * factorial(half_index)
        / ((Fraction(44, 7) ** pi_power) * factorial(2 * half_index))
    )
    radial_factor = (radius_s / (2 * kappa)) ** half_index
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D19 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 405
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D19 post-m29 sharpened tail margin is not positive")

    print(f"D_19 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_19 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_19 exp(A) upper = {exp_upper}")
    print(f"D_19 sine product factor = {sine_product_factor}")
    print(f"D_19 quartic moment = {quartic_moment}")
    print(f"D_19 quartic tail upper = {tail_weight_upper}")
    print(f"D_19 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_19 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_19 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_19 Pr cap lower = {ball_probability_lower}")
    print(f"D_19 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_19 Chebyshev degree = {chebyshev_degree}")
    print(f"D_19 Chebyshev scale = {chebyshev_scale}")
    print(f"D_19 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_19 Chebyshev denominator = {denominator}")
    print(f"D_19 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_19 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_19 tail margin = {tail_margin}")


def d20_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D20 finite bridge and Chebyshev tail from n=405 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D20 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 20
    first_boundary = rank
    exact_through = 41
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D20 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D20 supplied determinant deltas are negative: {negative_deltas}")

    print("D20 post-m29 finite-tail certificate")
    print(f"D_20 supplied exact deltas through m_{exact_through}")
    print("D_20 bridge uses stable absolute-moment envelopes beyond m_41")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 200
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D20 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D20 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_20 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_20 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 780
    kappa = 38
    c_value = 20
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(362, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 380)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(4, 1)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D20 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D20 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    density_factor = Fraction(1, 1) / (
        (Fraction(44, 7) ** (rank // 2)) * factorial(dimension // 2)
    )
    radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D20 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 405
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D20 post-m29 sharpened tail margin is not positive")

    print(f"D_20 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_20 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_20 exp(A) upper = {exp_upper}")
    print(f"D_20 sine product factor = {sine_product_factor}")
    print(f"D_20 quartic moment = {quartic_moment}")
    print(f"D_20 quartic tail upper = {tail_weight_upper}")
    print(f"D_20 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_20 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_20 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_20 Pr cap lower = {ball_probability_lower}")
    print(f"D_20 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_20 Chebyshev degree = {chebyshev_degree}")
    print(f"D_20 Chebyshev scale = {chebyshev_scale}")
    print(f"D_20 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_20 Chebyshev denominator = {denominator}")
    print(f"D_20 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_20 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_20 tail margin = {tail_margin}")


def d21_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D21 finite bridge and Chebyshev tail from n=405 onward."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    if not delta_logs:
        raise SystemExit("D21 finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    rank = 21
    first_boundary = rank
    exact_through = 40
    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {})
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"D21 supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(f"D21 supplied determinant deltas are negative: {negative_deltas}")

    print("D21 post-m29 finite-tail certificate")
    print(f"D_21 supplied exact deltas through m_{exact_through}")
    print("D_21 bridge uses stable absolute-moment envelopes beyond m_40")
    for path in delta_logs:
        print(f"  {path}")

    bridge_max_chain = 200
    stable_all = stable_moments_formula(2 * bridge_max_chain + 5)
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit("D21 stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"D21 m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        print(
            f"D_21 m={chain_m}: stable D={stable_diff}; "
            f"negative linear exact={negative_linear_exact}; "
            f"negative linear stable bound={negative_linear_stable_bound}; "
            f"negative quadratic exact={negative_quadratic_exact}; "
            f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
            f"stable-bound linear indices={negative_linear_stable_indices}; "
            f"stable-bound pairs={negative_quadratic_stable_pairs}; "
            f"lower margin={lower_margin}"
        )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"D_21 minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    dimension = 861
    kappa = 40
    c_value = 21
    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    radius_s = Fraction(360, 1)
    central_t = Fraction(1, 1)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    exp_upper = exp_upper_bound(sine_exponent_bound, 420)
    sine_product_factor = Fraction(1, 1) / exp_upper

    quartic_r = Fraction(1, 2)
    quartic_s = Fraction(1, 1)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s + 4 * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit("D21 quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit("D21 quartic pushforward monotonicity check failed")

    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    half_index = (dimension + 1) // 2
    pi_power = (rank + 1) // 2
    density_factor = (
        Fraction(7, 5)
        * 4**half_index
        * factorial(half_index)
        / ((Fraction(44, 7) ** pi_power) * factorial(2 * half_index))
    )
    radial_factor = (radius_s / (2 * kappa)) ** half_index
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit("D21 Chebyshev negative-tail upper bound is not positive")

    target_exponent = 405
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit("D21 post-m29 sharpened tail margin is not positive")

    print(f"D_21 cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"D_21 sine exponent upper A = {sine_exponent_bound}")
    print(f"D_21 exp(A) upper = {exp_upper}")
    print(f"D_21 sine product factor = {sine_product_factor}")
    print(f"D_21 quartic moment = {quartic_moment}")
    print(f"D_21 quartic tail upper = {tail_weight_upper}")
    print(f"D_21 pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"D_21 pushforward second derivative = {pushforward_second_derivative}")
    print(f"D_21 pushforward weight lower = {pushforward_weight_lower}")
    print(f"D_21 Pr cap lower = {ball_probability_lower}")
    print(f"D_21 P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"D_21 Chebyshev degree = {chebyshev_degree}")
    print(f"D_21 Chebyshev scale = {chebyshev_scale}")
    print(f"D_21 Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"D_21 Chebyshev denominator = {denominator}")
    print(f"D_21 N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"D_21 sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"D_21 tail margin = {tail_margin}")


def d_high_rank_post_m29_finite_tail_certificate(
    delta_logs: list[Path],
    label: str,
    row_label: str,
    rank: int,
    exact_through: int,
    dimension: int,
    kappa: int,
    c_value: int,
    radius_s: Fraction,
    exp_terms: int,
    quartic_s: Fraction,
    density_kind: str,
    central_t: Fraction = Fraction(1, 1),
    quartic_r: Fraction = Fraction(1, 2),
    target_exponent: int = 405,
    bridge_max_chain: int = 200,
    verbose_bridge_terms: bool = True,
) -> None:
    """Replay a high-rank D finite bridge and Chebyshev tail."""

    def poly_add(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        length = max(len(first), len(second))
        out = [Fraction(0) for _ in range(length)]
        for index, value in enumerate(first):
            out[index] += value
        for index, value in enumerate(second):
            out[index] += value
        return out

    def poly_scale(poly: list[Fraction], factor: Fraction) -> list[Fraction]:
        return [factor * value for value in poly]

    def poly_mul(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(first) + len(second) - 1)]
        for first_index, first_value in enumerate(first):
            for second_index, second_value in enumerate(second):
                out[first_index + second_index] += first_value * second_value
        return out

    def chebyshev_t(degree: int) -> list[Fraction]:
        if degree == 0:
            return [Fraction(1)]
        if degree == 1:
            return [Fraction(0), Fraction(1)]
        previous = [Fraction(1)]
        current = [Fraction(0), Fraction(1)]
        for _index in range(2, degree + 1):
            previous, current = current, poly_add(
                poly_scale([Fraction(0)] + current, Fraction(2)),
                poly_scale(previous, Fraction(-1)),
            )
        return current

    def compose(poly: list[Fraction], argument: list[Fraction]) -> list[Fraction]:
        out = [Fraction(0)]
        argument_power = [Fraction(1)]
        for coefficient in poly:
            out = poly_add(out, poly_scale(argument_power, coefficient))
            argument_power = poly_mul(argument_power, argument)
        return out

    def eval_poly(poly: list[Fraction], value: Fraction) -> Fraction:
        out = Fraction(0)
        value_power = Fraction(1)
        for coefficient in poly:
            out += coefficient * value_power
            value_power *= value
        return out

    first_boundary = rank
    if not delta_logs and exact_through >= first_boundary:
        raise SystemExit(f"{label} finite-tail replay needs --d-delta-log inputs")
    for path in delta_logs:
        text = path.read_text()
        if "__EXIT_STATUS=0" not in text:
            raise SystemExit(f"{path} does not record __EXIT_STATUS=0")

    parsed_deltas = parse_d_delta_logs(delta_logs).get(rank, {}) if delta_logs else {}
    deltas = {
        index: value
        for index, value in parsed_deltas.items()
        if index <= exact_through
    }
    expected_indices = list(range(first_boundary, exact_through + 1))
    missing_indices = [index for index in expected_indices if index not in deltas]
    if missing_indices:
        raise SystemExit(f"{label} supplied deltas have gaps: {missing_indices}")
    negative_deltas = {
        index: value for index, value in deltas.items() if value < 0
    }
    if negative_deltas:
        raise SystemExit(
            f"{label} supplied determinant deltas are negative: {negative_deltas}"
        )

    print(f"{label} post-m29 finite-tail certificate")
    if exact_through < first_boundary:
        print(f"{row_label} uses no exact determinant deltas before m_{first_boundary}")
    else:
        print(f"{row_label} supplied exact deltas through m_{exact_through}")
    print(
        f"{row_label} bridge uses stable absolute-moment envelopes beyond "
        f"m_{exact_through}"
    )
    for path in delta_logs:
        print(f"  {path}")

    stable_all = stable_moments_formula(
        2 * bridge_max_chain + 5,
        progress_label=f"{row_label} stable moments",
    )
    table_check = read_stable_moments(100)
    if stable_all[: len(table_check)] != table_check:
        raise SystemExit(f"{label} stable moment formula disagrees with the checked table")

    chain_margins: list[tuple[int, int]] = []
    for chain_m in range(31, bridge_max_chain + 1):
        moment_max = 2 * chain_m + 5
        stable = stable_all[: moment_max + 1]
        stable_diff = chain_diff_integer(stable, chain_m)
        linear_coeffs = chain_diff_linear_coefficients(stable, chain_m)
        quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m)
        negative_linear_indices = [
            index
            for index, value in enumerate(linear_coeffs)
            if value < 0 and index >= first_boundary
        ]
        negative_quadratic = [
            (first, second, value)
            for (first, second), value in sorted(quadratic_coeffs.items())
            if value < 0 and first >= first_boundary and second >= first_boundary
        ]

        negative_linear_exact = 0
        negative_linear_stable_bound = 0
        negative_linear_stable_indices: list[int] = []
        negative_linear_stable_count = 0
        for index in negative_linear_indices:
            value = linear_coeffs[index]
            if index in deltas:
                negative_linear_exact += value * deltas[index]
            else:
                negative_linear_stable_bound += value * stable[index]
                negative_linear_stable_count += 1
                if verbose_bridge_terms:
                    negative_linear_stable_indices.append(index)

        negative_quadratic_exact = 0
        negative_quadratic_stable_bound = 0
        negative_quadratic_stable_pairs: list[tuple[int, int]] = []
        negative_quadratic_stable_count = 0
        for first, second, value in negative_quadratic:
            if first in deltas and second in deltas:
                negative_quadratic_exact += value * deltas[first] * deltas[second]
            else:
                first_bound = abs(deltas[first]) if first in deltas else stable[first]
                second_bound = abs(deltas[second]) if second in deltas else stable[second]
                negative_quadratic_stable_bound += value * first_bound * second_bound
                negative_quadratic_stable_count += 1
                if verbose_bridge_terms:
                    negative_quadratic_stable_pairs.append((first, second))

        lower_margin = (
            stable_diff
            + negative_linear_exact
            + negative_linear_stable_bound
            + negative_quadratic_exact
            + negative_quadratic_stable_bound
        )
        if lower_margin <= 0:
            raise SystemExit(f"{label} m={chain_m} lower margin is not positive")
        chain_margins.append((chain_m, lower_margin))
        if verbose_bridge_terms:
            print(
                f"{row_label} m={chain_m}: stable D={stable_diff}; "
                f"negative linear exact={negative_linear_exact}; "
                f"negative linear stable bound={negative_linear_stable_bound}; "
                f"negative quadratic exact={negative_quadratic_exact}; "
                f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
                f"stable-bound linear indices={negative_linear_stable_indices}; "
                f"stable-bound pairs={negative_quadratic_stable_pairs}; "
                f"lower margin={lower_margin}"
            )
        else:
            print(
                f"{row_label} m={chain_m}: stable D={stable_diff}; "
                f"negative linear exact={negative_linear_exact}; "
                f"negative linear stable bound={negative_linear_stable_bound}; "
                f"negative linear stable count={negative_linear_stable_count}; "
                f"negative quadratic exact={negative_quadratic_exact}; "
                f"negative quadratic stable bound={negative_quadratic_stable_bound}; "
                f"negative quadratic stable count={negative_quadratic_stable_count}; "
                f"lower margin={lower_margin}"
            )

    minimum_chain_m, minimum_chain_margin = min(chain_margins, key=lambda item: item[1])
    print(f"{row_label} minimum bridge margin m={minimum_chain_m}: {minimum_chain_margin}")

    weyl_order = 2 ** (rank - 1) * factorial(rank)
    degrees = tuple(2 * index for index in range(1, rank)) + (rank,)
    x_floor = Fraction(dimension) - radius_s
    c_push = x_floor - central_t

    q_bound = Fraction(2) * radius_s / (39 * kappa)
    sine_exponent_bound = (
        radius_s / 12
        + radius_s * radius_s / (1440 * kappa)
        + radius_s**3 / (kappa * kappa * 90720 * (1 - q_bound))
    )
    print(f"{row_label} direct tail exp bound start", flush=True)
    exp_upper = exp_upper_bound(
        sine_exponent_bound,
        exp_terms,
        progress_label=row_label,
    )
    print(f"{row_label} direct tail exp bound done", flush=True)
    sine_product_factor = Fraction(1, 1) / exp_upper

    print(f"{row_label} direct tail quartic setup", flush=True)
    quartic_moment = (
        Fraction(6)
        - 2 * (quartic_r + quartic_s)
        + (
            quartic_r * quartic_r
            + quartic_s * quartic_s
            + 4 * quartic_r * quartic_s
        )
        + quartic_r * quartic_r * quartic_s * quartic_s
    )
    quartic_denominator = (central_t + quartic_r) ** 2 * (central_t + quartic_s) ** 2
    tail_weight_upper = (x_floor + central_t) ** 2 * quartic_moment / quartic_denominator
    pushforward_weight_lower = x_floor * x_floor + 1 - tail_weight_upper
    pushforward_derivative_lower = (
        2 * x_floor
        - 2 * quartic_moment * (x_floor + central_t) / quartic_denominator
    )
    pushforward_second_derivative = 2 - 2 * quartic_moment / quartic_denominator
    if pushforward_weight_lower <= 0:
        raise SystemExit(f"{label} quartic pushforward weight lower bound is not positive")
    if pushforward_derivative_lower <= 0 or pushforward_second_derivative <= 0:
        raise SystemExit(f"{label} quartic pushforward monotonicity check failed")

    print(f"{row_label} direct tail density setup", flush=True)
    degree_factor = 1
    for degree in degrees:
        degree_factor *= factorial(degree)
    if density_kind == "even":
        density_factor = Fraction(1, 1) / (
            (Fraction(44, 7) ** (rank // 2)) * factorial(dimension // 2)
        )
        radial_factor = (radius_s / (2 * kappa)) ** (dimension // 2)
    elif density_kind == "odd":
        half_index = (dimension + 1) // 2
        pi_power = (rank + 1) // 2
        density_factor = (
            Fraction(7, 5)
            * 4**half_index
            * factorial(half_index)
            / ((Fraction(44, 7) ** pi_power) * factorial(2 * half_index))
        )
        radial_factor = (radius_s / (2 * kappa)) ** half_index
    else:
        raise ValueError(f"unknown density kind {density_kind}")
    ball_probability_lower = (
        Fraction(1, weyl_order)
        * sine_product_factor
        * degree_factor
        * density_factor
        * radial_factor
    )
    weighted_tail_lower = pushforward_weight_lower * ball_probability_lower

    print(f"{row_label} direct tail Chebyshev setup", flush=True)
    chebyshev_degree = 8
    stable_moments = stable_moments_formula(4 * chebyshev_degree + 2)
    moments = stable_moments[:]
    for index, value in deltas.items():
        if index < len(moments):
            moments[index] += value
    upper_support = Fraction(2 * dimension, 1)
    chebyshev_scale = Fraction(12000, 1)
    negative_cutoff = Fraction(80, 1)
    q0 = negative_cutoff * (negative_cutoff + upper_support)
    chebyshev_poly = chebyshev_t(chebyshev_degree)
    argument_poly = [
        Fraction(1),
        -2 * upper_support / chebyshev_scale,
        Fraction(2, 1) / chebyshev_scale,
    ]
    tail_majorant = poly_mul(
        compose(chebyshev_poly, argument_poly),
        compose(chebyshev_poly, argument_poly),
    )
    denominator = eval_poly(
        chebyshev_poly,
        1 + 2 * q0 / chebyshev_scale,
    ) ** 2
    q3_moments = [q3_integer(moments, index) for index in range(len(tail_majorant))]
    negative_tail_mass_upper = (
        sum(
            coefficient * q3_moments[index]
            for index, coefficient in enumerate(tail_majorant)
        )
        / denominator
    )
    if negative_tail_mass_upper <= 0:
        raise SystemExit(f"{label} Chebyshev negative-tail upper bound is not positive")

    print(f"{row_label} direct tail margin setup", flush=True)
    sharpened_target = (
        negative_tail_mass_upper
        * Fraction(2 * c_value, 1) ** target_exponent
        / c_push**target_exponent
        + 2 * negative_cutoff**target_exponent / c_push**target_exponent
    )
    tail_margin = weighted_tail_lower - sharpened_target
    if tail_margin <= 0:
        raise SystemExit(f"{label} post-m29 sharpened tail margin is not positive")

    print(f"{row_label} cap Q(theta) <= {radius_s}; chi_adj >= {x_floor}")
    print(f"{row_label} sine exponent upper A = {sine_exponent_bound}")
    print(f"{row_label} exp(A) upper = {exp_upper}")
    print(f"{row_label} sine product factor = {sine_product_factor}")
    print(f"{row_label} quartic moment = {quartic_moment}")
    print(f"{row_label} quartic tail upper = {tail_weight_upper}")
    print(f"{row_label} pushforward derivative at x_floor = {pushforward_derivative_lower}")
    print(f"{row_label} pushforward second derivative = {pushforward_second_derivative}")
    print(f"{row_label} pushforward weight lower = {pushforward_weight_lower}")
    print(f"{row_label} Pr cap lower = {ball_probability_lower}")
    print(f"{row_label} P_G({c_push}) lower = {weighted_tail_lower}")
    print(f"{row_label} Chebyshev degree = {chebyshev_degree}")
    print(f"{row_label} Chebyshev scale = {chebyshev_scale}")
    print(f"{row_label} Chebyshev z0 = {1 + 2 * q0 / chebyshev_scale}")
    print(f"{row_label} Chebyshev denominator = {denominator}")
    print(f"{row_label} N_G({negative_cutoff}) upper = {negative_tail_mass_upper}")
    print(f"{row_label} sharpened target at n={target_exponent} = {sharpened_target}")
    print(f"{row_label} tail margin = {tail_margin}")


def d22_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D22 finite bridge and Chebyshev tail from n=405 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D22",
        row_label="D_22",
        rank=22,
        exact_through=39,
        dimension=946,
        kappa=42,
        c_value=22,
        radius_s=Fraction(360, 1),
        exp_terms=460,
        quartic_s=Fraction(1, 1),
        density_kind="even",
    )


def d23_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D23 finite bridge and Chebyshev tail from n=405 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D23",
        row_label="D_23",
        rank=23,
        exact_through=38,
        dimension=1035,
        kappa=44,
        c_value=23,
        radius_s=Fraction(360, 1),
        exp_terms=500,
        quartic_s=Fraction(1, 1),
        density_kind="odd",
    )


def d24_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D24 finite bridge and Chebyshev tail from n=405 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D24",
        row_label="D_24",
        rank=24,
        exact_through=37,
        dimension=1128,
        kappa=46,
        c_value=24,
        radius_s=Fraction(410, 1),
        exp_terms=540,
        quartic_s=Fraction(1, 1),
        density_kind="even",
    )


def d25_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D25 finite bridge and Chebyshev tail from n=405 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D25",
        row_label="D_25",
        rank=25,
        exact_through=36,
        dimension=1225,
        kappa=48,
        c_value=25,
        radius_s=Fraction(500, 1),
        exp_terms=760,
        quartic_s=Fraction(1, 1),
        density_kind="odd",
    )


def d26_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D26 finite bridge and Chebyshev tail from n=405 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D26",
        row_label="D_26",
        rank=26,
        exact_through=35,
        dimension=1326,
        kappa=50,
        c_value=26,
        radius_s=Fraction(640, 1),
        exp_terms=900,
        quartic_s=Fraction(1, 1),
        density_kind="even",
    )


def d27_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D27 finite bridge and Chebyshev tail from n=413 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D27",
        row_label="D_27",
        rank=27,
        exact_through=35,
        dimension=1431,
        kappa=52,
        c_value=27,
        radius_s=Fraction(820, 1),
        exp_terms=1200,
        quartic_s=Fraction(1, 1),
        density_kind="odd",
        target_exponent=413,
        bridge_max_chain=204,
    )


def d28_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D28 finite bridge and Chebyshev tail from n=433 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D28",
        row_label="D_28",
        rank=28,
        exact_through=35,
        dimension=1540,
        kappa=54,
        c_value=28,
        radius_s=Fraction(870, 1),
        exp_terms=1400,
        quartic_s=Fraction(1, 3),
        density_kind="even",
        central_t=Fraction(5, 4),
        quartic_r=Fraction(16, 5),
        target_exponent=433,
        bridge_max_chain=214,
    )


def d29_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D29 finite bridge and Chebyshev tail from n=449 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D29",
        row_label="D_29",
        rank=29,
        exact_through=35,
        dimension=1653,
        kappa=56,
        c_value=29,
        radius_s=Fraction(925, 1),
        exp_terms=1400,
        quartic_s=Fraction(7, 2),
        density_kind="odd",
        central_t=Fraction(5, 4),
        quartic_r=Fraction(1, 3),
        target_exponent=449,
        bridge_max_chain=222,
    )


def d30_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D30 finite bridge and Chebyshev tail from n=481 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D30",
        row_label="D_30",
        rank=30,
        exact_through=35,
        dimension=1770,
        kappa=58,
        c_value=30,
        radius_s=Fraction(970, 1),
        exp_terms=1600,
        quartic_s=Fraction(1, 4),
        density_kind="even",
        central_t=Fraction(4, 3),
        quartic_r=Fraction(10, 3),
        target_exponent=481,
        bridge_max_chain=238,
    )


def d31_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D31 finite bridge and Chebyshev tail from n=497 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D31",
        row_label="D_31",
        rank=31,
        exact_through=35,
        dimension=1891,
        kappa=60,
        c_value=29,
        radius_s=Fraction(1020, 1),
        exp_terms=1800,
        quartic_s=Fraction(1, 4),
        density_kind="odd",
        central_t=Fraction(4, 3),
        quartic_r=Fraction(10, 3),
        target_exponent=497,
        bridge_max_chain=246,
    )


def d32_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D32 finite bridge and Chebyshev tail from n=513 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D32",
        row_label="D_32",
        rank=32,
        exact_through=35,
        dimension=2016,
        kappa=62,
        c_value=32,
        radius_s=Fraction(1067, 1),
        exp_terms=2000,
        quartic_s=Fraction(1, 4),
        density_kind="even",
        central_t=Fraction(4, 3),
        quartic_r=Fraction(10, 3),
        target_exponent=513,
        bridge_max_chain=254,
    )


def d33_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D33 finite bridge and Chebyshev tail from n=545 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D33",
        row_label="D_33",
        rank=33,
        exact_through=35,
        dimension=2145,
        kappa=64,
        c_value=31,
        radius_s=Fraction(1110, 1),
        exp_terms=2200,
        quartic_s=Fraction(1, 4),
        density_kind="odd",
        central_t=Fraction(4, 3),
        quartic_r=Fraction(10, 3),
        target_exponent=545,
        bridge_max_chain=270,
    )


def d34_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D34 finite bridge and Chebyshev tail from n=577 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D34",
        row_label="D_34",
        rank=34,
        exact_through=35,
        dimension=2278,
        kappa=66,
        c_value=34,
        radius_s=Fraction(1152, 1),
        exp_terms=2400,
        quartic_s=Fraction(1, 5),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(3, 1),
        target_exponent=577,
        bridge_max_chain=286,
    )


def d35_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D35 finite bridge and Chebyshev tail from n=593 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D35",
        row_label="D_35",
        rank=35,
        exact_through=35,
        dimension=2415,
        kappa=68,
        c_value=33,
        radius_s=Fraction(1194, 1),
        exp_terms=2600,
        quartic_s=Fraction(1, 4),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(3, 1),
        target_exponent=593,
        bridge_max_chain=294,
    )


def d36_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D36 finite bridge and Chebyshev tail from n=625 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D36",
        row_label="D_36",
        rank=36,
        exact_through=35,
        dimension=2556,
        kappa=70,
        c_value=36,
        radius_s=Fraction(1235, 1),
        exp_terms=2800,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=625,
        bridge_max_chain=310,
    )


def d37_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D37 finite bridge and Chebyshev tail from n=657 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D37",
        row_label="D_37",
        rank=37,
        exact_through=35,
        dimension=2701,
        kappa=72,
        c_value=35,
        radius_s=Fraction(1275, 1),
        exp_terms=3200,
        quartic_s=Fraction(13, 4),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=657,
        bridge_max_chain=327,
    )


def d38_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D38 finite bridge and Chebyshev tail from n=673 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D38",
        row_label="D_38",
        rank=38,
        exact_through=35,
        dimension=2850,
        kappa=74,
        c_value=38,
        radius_s=Fraction(1316, 1),
        exp_terms=3400,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=673,
        bridge_max_chain=335,
    )


def d39_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D39 finite bridge and Chebyshev tail from n=705 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D39",
        row_label="D_39",
        rank=39,
        exact_through=35,
        dimension=3003,
        kappa=76,
        c_value=37,
        radius_s=Fraction(1356, 1),
        exp_terms=3600,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=705,
        bridge_max_chain=351,
    )


def d40_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D40 finite bridge and Chebyshev tail from n=737 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D40",
        row_label="D_40",
        rank=40,
        exact_through=35,
        dimension=3160,
        kappa=78,
        c_value=40,
        radius_s=Fraction(1396, 1),
        exp_terms=3800,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=737,
        bridge_max_chain=367,
    )


def d41_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D41 finite bridge and Chebyshev tail from n=769 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D41",
        row_label="D_41",
        rank=41,
        exact_through=35,
        dimension=3321,
        kappa=80,
        c_value=39,
        radius_s=Fraction(1434, 1),
        exp_terms=4000,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=769,
        bridge_max_chain=383,
    )


def d42_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D42 finite bridge and Chebyshev tail from n=801 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D42",
        row_label="D_42",
        rank=42,
        exact_through=35,
        dimension=3486,
        kappa=82,
        c_value=42,
        radius_s=Fraction(1474, 1),
        exp_terms=4200,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=801,
        bridge_max_chain=399,
    )


def d43_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D43 finite bridge and Chebyshev tail from n=833 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D43",
        row_label="D_43",
        rank=43,
        exact_through=35,
        dimension=3655,
        kappa=84,
        c_value=41,
        radius_s=Fraction(1513, 1),
        exp_terms=4400,
        quartic_s=Fraction(7, 2),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=833,
        bridge_max_chain=415,
    )


def d44_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D44 finite bridge and Chebyshev tail from n=881 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D44",
        row_label="D_44",
        rank=44,
        exact_through=35,
        dimension=3828,
        kappa=86,
        c_value=44,
        radius_s=Fraction(1552, 1),
        exp_terms=4600,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=881,
        bridge_max_chain=439,
    )


def d45_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D45 finite bridge and Chebyshev tail from n=913 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D45",
        row_label="D_45",
        rank=45,
        exact_through=35,
        dimension=4005,
        kappa=88,
        c_value=43,
        radius_s=Fraction(1590, 1),
        exp_terms=4800,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(7, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=913,
        bridge_max_chain=455,
    )


def d46_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D46 finite bridge and Chebyshev tail from n=977 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D46",
        row_label="D_46",
        rank=46,
        exact_through=35,
        dimension=4186,
        kappa=90,
        c_value=46,
        radius_s=Fraction(1628, 1),
        exp_terms=5000,
        quartic_s=Fraction(7, 2),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=977,
        bridge_max_chain=487,
    )


def d47_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D47 finite bridge and Chebyshev tail from n=1009 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D47",
        row_label="D_47",
        rank=47,
        exact_through=35,
        dimension=4371,
        kappa=92,
        c_value=45,
        radius_s=Fraction(1667, 1),
        exp_terms=5200,
        quartic_s=Fraction(13, 4),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1009,
        bridge_max_chain=503,
    )


def d48_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D48 finite bridge and Chebyshev tail from n=1073 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D48",
        row_label="D_48",
        rank=48,
        exact_through=35,
        dimension=4560,
        kappa=94,
        c_value=48,
        radius_s=Fraction(1705, 1),
        exp_terms=5400,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1073,
        bridge_max_chain=535,
    )


def d49_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D49 finite bridge and Chebyshev tail from n=1093 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D49",
        row_label="D_49",
        rank=49,
        exact_through=35,
        dimension=4753,
        kappa=96,
        c_value=47,
        radius_s=Fraction(1744, 1),
        exp_terms=5600,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1093,
        bridge_max_chain=545,
    )


def d50_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D50 finite bridge and Chebyshev tail from n=1155 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D50",
        row_label="D_50",
        rank=50,
        exact_through=35,
        dimension=4950,
        kappa=98,
        c_value=50,
        radius_s=Fraction(1782, 1),
        exp_terms=5800,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1155,
        bridge_max_chain=576,
    )


def d51_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D51 finite bridge and Chebyshev tail from n=1191 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D51",
        row_label="D_51",
        rank=51,
        exact_through=35,
        dimension=5151,
        kappa=100,
        c_value=49,
        radius_s=Fraction(1821, 1),
        exp_terms=6000,
        quartic_s=Fraction(7, 2),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1191,
        bridge_max_chain=594,
    )


def d52_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D52 finite bridge and Chebyshev tail from n=1255 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D52",
        row_label="D_52",
        rank=52,
        exact_through=35,
        dimension=5356,
        kappa=102,
        c_value=52,
        radius_s=Fraction(1859, 1),
        exp_terms=6200,
        quartic_s=Fraction(7, 2),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1255,
        bridge_max_chain=626,
    )


def d53_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D53 finite bridge and Chebyshev tail from n=1293 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D53",
        row_label="D_53",
        rank=53,
        exact_through=35,
        dimension=5565,
        kappa=104,
        c_value=51,
        radius_s=Fraction(1897, 1),
        exp_terms=6400,
        quartic_s=Fraction(31, 10),
        density_kind="odd",
        central_t=Fraction(8, 5),
        quartic_r=Fraction(1, 5),
        target_exponent=1293,
        bridge_max_chain=645,
    )


def d54_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D54 finite bridge and Chebyshev tail from n=1359 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D54",
        row_label="D_54",
        rank=54,
        exact_through=35,
        dimension=5778,
        kappa=106,
        c_value=54,
        radius_s=Fraction(1936, 1),
        exp_terms=6600,
        quartic_s=Fraction(25, 8),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1359,
        bridge_max_chain=678,
    )


def d55_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D55 finite bridge and Chebyshev tail from n=1399 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D55",
        row_label="D_55",
        rank=55,
        exact_through=35,
        dimension=5995,
        kappa=108,
        c_value=53,
        radius_s=Fraction(1974, 1),
        exp_terms=6800,
        quartic_s=Fraction(29, 10),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(3, 10),
        target_exponent=1399,
        bridge_max_chain=698,
    )


def d56_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D56 finite bridge and Chebyshev tail from n=1469 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D56",
        row_label="D_56",
        rank=56,
        exact_through=35,
        dimension=6216,
        kappa=110,
        c_value=56,
        radius_s=Fraction(2012, 1),
        exp_terms=7200,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1469,
        bridge_max_chain=733,
    )


def d57_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D57 finite bridge and Chebyshev tail from n=1511 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D57",
        row_label="D_57",
        rank=57,
        exact_through=35,
        dimension=6441,
        kappa=112,
        c_value=55,
        radius_s=Fraction(2051, 1),
        exp_terms=7400,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1511,
        bridge_max_chain=754,
    )


def d58_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D58 finite bridge and Chebyshev tail from n=1585 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D58",
        row_label="D_58",
        rank=58,
        exact_through=35,
        dimension=6670,
        kappa=114,
        c_value=58,
        radius_s=Fraction(2089, 1),
        exp_terms=7600,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=1585,
        bridge_max_chain=791,
    )


def d59_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D59 finite bridge and Chebyshev tail from n=1627 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D59",
        row_label="D_59",
        rank=59,
        exact_through=35,
        dimension=6903,
        kappa=116,
        c_value=57,
        radius_s=Fraction(2127, 1),
        exp_terms=7800,
        quartic_s=Fraction(63, 20),
        density_kind="odd",
        central_t=Fraction(47, 30),
        quartic_r=Fraction(4, 15),
        target_exponent=1627,
        bridge_max_chain=812,
    )


def d60_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D60 finite bridge and Chebyshev tail from n=1703 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D60",
        row_label="D_60",
        rank=60,
        exact_through=35,
        dimension=7140,
        kappa=118,
        c_value=60,
        radius_s=Fraction(2165, 1),
        exp_terms=8000,
        quartic_s=Fraction(1, 5),
        density_kind="even",
        central_t=Fraction(8, 5),
        quartic_r=Fraction(16, 5),
        target_exponent=1703,
        bridge_max_chain=850,
    )


def d61_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D61 finite bridge and Chebyshev tail from n=1753 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D61",
        row_label="D_61",
        rank=61,
        exact_through=35,
        dimension=7381,
        kappa=120,
        c_value=59,
        radius_s=Fraction(2204, 1),
        exp_terms=8200,
        quartic_s=Fraction(47, 200),
        density_kind="odd",
        central_t=Fraction(38, 25),
        quartic_r=Fraction(141, 50),
        target_exponent=1753,
        bridge_max_chain=875,
    )


def d62_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D62 finite bridge and Chebyshev tail from n=1833 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D62",
        row_label="D_62",
        rank=62,
        exact_through=35,
        dimension=7626,
        kappa=122,
        c_value=62,
        radius_s=Fraction(2242, 1),
        exp_terms=8600,
        quartic_s=Fraction(21, 100),
        density_kind="even",
        central_t=Fraction(77, 50),
        quartic_r=Fraction(13, 5),
        target_exponent=1833,
        bridge_max_chain=915,
    )


def d63_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D63 finite bridge and Chebyshev tail from n=1891 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D63",
        row_label="D_63",
        rank=63,
        exact_through=35,
        dimension=7875,
        kappa=124,
        c_value=61,
        radius_s=Fraction(2281, 1),
        exp_terms=9000,
        quartic_s=Fraction(21, 100),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(313, 100),
        target_exponent=1891,
        bridge_max_chain=944,
    )


def d64_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D64 finite bridge and Chebyshev tail from n=1965 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D64",
        row_label="D_64",
        rank=64,
        exact_through=35,
        dimension=8128,
        kappa=126,
        c_value=64,
        radius_s=Fraction(2319, 1),
        exp_terms=9400,
        quartic_s=Fraction(8, 25),
        density_kind="even",
        central_t=Fraction(39, 25),
        quartic_r=Fraction(73, 25),
        target_exponent=1965,
        bridge_max_chain=981,
    )


def d65_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D65 finite bridge and Chebyshev tail from n=2061 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D65",
        row_label="D_65",
        rank=65,
        exact_through=35,
        dimension=8385,
        kappa=128,
        c_value=63,
        radius_s=Fraction(2357, 1),
        exp_terms=10000,
        quartic_s=Fraction(21, 100),
        density_kind="odd",
        central_t=Fraction(157, 100),
        quartic_r=Fraction(283, 100),
        target_exponent=2061,
        bridge_max_chain=1029,
    )


def d66_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D66 finite bridge and Chebyshev tail from n=2101 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D66",
        row_label="D_66",
        rank=66,
        exact_through=35,
        dimension=8646,
        kappa=130,
        c_value=66,
        radius_s=Fraction(2395, 1),
        exp_terms=10400,
        quartic_s=Fraction(21, 100),
        density_kind="even",
        central_t=Fraction(157, 100),
        quartic_r=Fraction(3, 1),
        target_exponent=2101,
        bridge_max_chain=1050,
    )


def d67_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D67 finite bridge and Chebyshev tail from n=2187 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D67",
        row_label="D_67",
        rank=67,
        exact_through=35,
        dimension=8911,
        kappa=132,
        c_value=65,
        radius_s=Fraction(2433, 1),
        exp_terms=6000,
        quartic_s=Fraction(67, 20),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=2187,
        bridge_max_chain=1092,
    )


def d68_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D68 finite bridge and Chebyshev tail from n=2251 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D68",
        row_label="D_68",
        rank=68,
        exact_through=35,
        dimension=9180,
        kappa=134,
        c_value=68,
        radius_s=Fraction(2471, 1),
        exp_terms=6000,
        quartic_s=Fraction(3, 1),
        density_kind="even",
        central_t=Fraction(8, 5),
        quartic_r=Fraction(1, 5),
        target_exponent=2251,
        bridge_max_chain=1124,
    )


def d69_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D69 finite bridge and Chebyshev tail from n=2281 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D69",
        row_label="D_69",
        rank=69,
        exact_through=35,
        dimension=9453,
        kappa=136,
        c_value=67,
        radius_s=Fraction(2510, 1),
        exp_terms=6000,
        quartic_s=Fraction(3, 1),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(1, 4),
        target_exponent=2281,
        bridge_max_chain=1139,
    )


def d70_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D70 finite bridge and Chebyshev tail from n=2369 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D70",
        row_label="D_70",
        rank=70,
        exact_through=35,
        dimension=9730,
        kappa=138,
        c_value=70,
        radius_s=Fraction(2548, 1),
        exp_terms=6000,
        quartic_s=Fraction(13, 5),
        density_kind="even",
        central_t=Fraction(8, 5),
        quartic_r=Fraction(1, 4),
        target_exponent=2369,
        bridge_max_chain=1183,
    )


def d71_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D71 finite bridge and Chebyshev tail from n=2439 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D71",
        row_label="D_71",
        rank=71,
        exact_through=35,
        dimension=10011,
        kappa=140,
        c_value=69,
        radius_s=Fraction(2587, 1),
        exp_terms=6000,
        quartic_s=Fraction(2, 9),
        density_kind="odd",
        central_t=Fraction(5, 3),
        quartic_r=Fraction(3, 1),
        target_exponent=2439,
        bridge_max_chain=1218,
    )


def d72_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D72 finite bridge and Chebyshev tail from n=2525 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D72",
        row_label="D_72",
        rank=72,
        exact_through=35,
        dimension=10296,
        kappa=142,
        c_value=72,
        radius_s=Fraction(2625, 1),
        exp_terms=6000,
        quartic_s=Fraction(1, 5),
        density_kind="even",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(3, 1),
        target_exponent=2525,
        bridge_max_chain=1261,
    )


def d73_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D73 finite bridge and Chebyshev tail from n=2575 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D73",
        row_label="D_73",
        rank=73,
        exact_through=35,
        dimension=10585,
        kappa=144,
        c_value=71,
        radius_s=Fraction(2662, 1),
        exp_terms=6000,
        quartic_s=Fraction(1, 5),
        density_kind="odd",
        central_t=Fraction(3, 2),
        quartic_r=Fraction(3, 1),
        target_exponent=2575,
        bridge_max_chain=1286,
        verbose_bridge_terms=False,
    )


def d74_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D74 finite bridge and Chebyshev tail from n=2669 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D74",
        row_label="D_74",
        rank=74,
        exact_through=35,
        dimension=10878,
        kappa=146,
        c_value=74,
        radius_s=Fraction(2701, 1),
        exp_terms=6000,
        quartic_s=Fraction(31, 100),
        density_kind="even",
        central_t=Fraction(31, 20),
        quartic_r=Fraction(103, 40),
        target_exponent=2669,
        bridge_max_chain=1334,
        verbose_bridge_terms=False,
    )


def d75_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D75 finite bridge and Chebyshev tail from n=2729 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D75",
        row_label="D_75",
        rank=75,
        exact_through=35,
        dimension=11175,
        kappa=148,
        c_value=73,
        radius_s=Fraction(2740, 1),
        exp_terms=6000,
        quartic_s=Fraction(69, 25),
        density_kind="odd",
        central_t=Fraction(157, 100),
        quartic_r=Fraction(8, 25),
        target_exponent=2729,
        bridge_max_chain=1363,
        verbose_bridge_terms=False,
    )


def d76_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D76 finite bridge and Chebyshev tail from n=2827 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D76",
        row_label="D_76",
        rank=76,
        exact_through=35,
        dimension=11476,
        kappa=150,
        c_value=76,
        radius_s=Fraction(2778, 1),
        exp_terms=6000,
        quartic_s=Fraction(7, 25),
        density_kind="even",
        central_t=Fraction(33, 20),
        quartic_r=Fraction(11, 4),
        target_exponent=2827,
        bridge_max_chain=1412,
        verbose_bridge_terms=False,
    )


def d77_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D77 finite bridge and Chebyshev tail from n=2889 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D77",
        row_label="D_77",
        rank=77,
        exact_through=35,
        dimension=11781,
        kappa=152,
        c_value=75,
        radius_s=Fraction(2816, 1),
        exp_terms=6000,
        quartic_s=Fraction(1, 5),
        density_kind="odd",
        central_t=Fraction(41, 25),
        quartic_r=Fraction(7, 2),
        target_exponent=2889,
        bridge_max_chain=1443,
        verbose_bridge_terms=False,
    )


def d78_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D78 finite bridge and Chebyshev tail from n=2989 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D78",
        row_label="D_78",
        rank=78,
        exact_through=35,
        dimension=12090,
        kappa=154,
        c_value=78,
        radius_s=Fraction(28537, 10),
        exp_terms=6000,
        quartic_s=Fraction(31, 10),
        density_kind="even",
        central_t=Fraction(8, 5),
        quartic_r=Fraction(11, 50),
        target_exponent=2989,
        bridge_max_chain=1493,
        verbose_bridge_terms=False,
    )


def d79_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D79 finite bridge and Chebyshev tail from n=3055 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D79",
        row_label="D_79",
        rank=79,
        exact_through=35,
        dimension=12403,
        kappa=156,
        c_value=77,
        radius_s=Fraction(144601, 50),
        exp_terms=6000,
        quartic_s=Fraction(281, 100),
        density_kind="odd",
        central_t=Fraction(19, 12),
        quartic_r=Fraction(13, 50),
        target_exponent=3055,
        bridge_max_chain=1526,
        verbose_bridge_terms=False,
    )


def d80_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D80 finite bridge and Chebyshev tail from n=3155 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D80",
        row_label="D_80",
        rank=80,
        exact_through=35,
        dimension=12720,
        kappa=158,
        c_value=80,
        radius_s=Fraction(146521, 50),
        exp_terms=6000,
        quartic_s=Fraction(109, 500),
        density_kind="even",
        central_t=Fraction(749, 500),
        quartic_r=Fraction(141, 50),
        target_exponent=3155,
        bridge_max_chain=1576,
        verbose_bridge_terms=False,
    )


def d81_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D81 finite bridge and Chebyshev tail from n=3221 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D81",
        row_label="D_81",
        rank=81,
        exact_through=35,
        dimension=13041,
        kappa=160,
        c_value=79,
        radius_s=Fraction(148439, 50),
        exp_terms=6000,
        quartic_s=Fraction(311, 100),
        density_kind="odd",
        central_t=Fraction(202, 125),
        quartic_r=Fraction(137, 500),
        target_exponent=3221,
        bridge_max_chain=1609,
        verbose_bridge_terms=False,
    )


def d82_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D82 finite bridge and Chebyshev tail from n=3327 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D82",
        row_label="D_82",
        rank=82,
        exact_through=35,
        dimension=13366,
        kappa=162,
        c_value=82,
        radius_s=Fraction(300697, 100),
        exp_terms=6000,
        quartic_s=Fraction(249, 1000),
        density_kind="even",
        central_t=Fraction(787, 500),
        quartic_r=Fraction(1699, 500),
        target_exponent=3327,
        bridge_max_chain=1662,
        verbose_bridge_terms=False,
    )


def d83_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D83 finite bridge and Chebyshev tail from n=3395 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D83",
        row_label="D_83",
        rank=83,
        exact_through=35,
        dimension=13695,
        kappa=164,
        c_value=81,
        radius_s=Fraction(304533, 100),
        exp_terms=6000,
        quartic_s=Fraction(49, 200),
        density_kind="odd",
        central_t=Fraction(81, 50),
        quartic_r=Fraction(79, 25),
        target_exponent=3395,
        bridge_max_chain=1696,
        verbose_bridge_terms=False,
    )


def d84_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D84 finite bridge and Chebyshev tail from n=3507 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D84",
        row_label="D_84",
        rank=84,
        exact_through=35,
        dimension=14028,
        kappa=166,
        c_value=84,
        radius_s=Fraction(30833575, 10000),
        exp_terms=6000,
        quartic_s=Fraction(180203, 50000),
        density_kind="even",
        central_t=Fraction(7887, 5000),
        quartic_r=Fraction(791, 2500),
        target_exponent=3507,
        bridge_max_chain=1752,
        verbose_bridge_terms=False,
    )


def d85_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D85 finite bridge and Chebyshev tail from n=3583 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D85",
        row_label="D_85",
        rank=85,
        exact_through=35,
        dimension=14365,
        kappa=168,
        c_value=83,
        radius_s=Fraction(31214533, 10000),
        exp_terms=6000,
        quartic_s=Fraction(13905, 50000),
        density_kind="odd",
        central_t=Fraction(2073, 1250),
        quartic_r=Fraction(3682, 1250),
        target_exponent=3583,
        bridge_max_chain=1790,
        verbose_bridge_terms=False,
    )


def d86_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D86 finite bridge and Chebyshev tail from n=3693 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D86",
        row_label="D_86",
        rank=86,
        exact_through=35,
        dimension=14706,
        kappa=170,
        c_value=86,
        radius_s=Fraction(31596951, 10000),
        exp_terms=6000,
        quartic_s=Fraction(138612, 50000),
        density_kind="even",
        central_t=Fraction(313189, 200000),
        quartic_r=Fraction(20015, 100000),
        target_exponent=3693,
        bridge_max_chain=1845,
        verbose_bridge_terms=False,
    )


def d87_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D87 finite bridge and Chebyshev tail from n=3769 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D87",
        row_label="D_87",
        rank=87,
        exact_through=35,
        dimension=15051,
        kappa=172,
        c_value=85,
        radius_s=Fraction(31980022, 10000),
        exp_terms=6000,
        quartic_s=Fraction(2577915, 1000000),
        density_kind="odd",
        central_t=Fraction(1657167, 1000000),
        quartic_r=Fraction(229201, 1000000),
        target_exponent=3769,
        bridge_max_chain=1883,
        verbose_bridge_terms=False,
    )


def d88_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D88 finite bridge and Chebyshev tail from n=3887 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D88",
        row_label="D_88",
        rank=88,
        exact_through=35,
        dimension=15400,
        kappa=174,
        c_value=88,
        radius_s=Fraction(32364259, 10000),
        exp_terms=6000,
        quartic_s=Fraction(3566194, 1000000),
        density_kind="even",
        central_t=Fraction(7669, 5000),
        quartic_r=Fraction(249111, 1000000),
        target_exponent=3887,
        bridge_max_chain=1942,
        verbose_bridge_terms=False,
    )


def d89_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D89 finite bridge and Chebyshev tail from n=3961 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D89",
        row_label="D_89",
        rank=89,
        exact_through=35,
        dimension=15753,
        kappa=176,
        c_value=87,
        radius_s=Fraction(327407696, 100000),
        exp_terms=6000,
        quartic_s=Fraction(1640615, 500000),
        density_kind="odd",
        central_t=Fraction(791943, 500000),
        quartic_r=Fraction(111003, 500000),
        target_exponent=3961,
        bridge_max_chain=1979,
        verbose_bridge_terms=False,
    )


def d90_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D90 finite bridge and Chebyshev tail from n=4085 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D90",
        row_label="D_90",
        rank=90,
        exact_through=35,
        dimension=16110,
        kappa=178,
        c_value=90,
        radius_s=Fraction(165648409, 50000),
        exp_terms=6000,
        quartic_s=Fraction(724017, 200000),
        density_kind="even",
        central_t=Fraction(317271, 200000),
        quartic_r=Fraction(37406, 125000),
        target_exponent=4085,
        bridge_max_chain=2041,
        verbose_bridge_terms=False,
    )


def d91_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D91 finite bridge and Chebyshev tail from n=4159 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D91",
        row_label="D_91",
        rank=91,
        exact_through=35,
        dimension=16471,
        kappa=180,
        c_value=89,
        radius_s=Fraction(335100941, 100000),
        exp_terms=6000,
        quartic_s=Fraction(227327, 1000000),
        density_kind="odd",
        central_t=Fraction(1582931, 1000000),
        quartic_r=Fraction(3139354, 1000000),
        target_exponent=4159,
        bridge_max_chain=2078,
        verbose_bridge_terms=False,
    )


def d92_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D92 finite bridge and Chebyshev tail from n=4289 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D92",
        row_label="D_92",
        rank=92,
        exact_through=35,
        dimension=16836,
        kappa=182,
        c_value=92,
        radius_s=Fraction(338929755, 100000),
        exp_terms=6000,
        quartic_s=Fraction(511, 2500),
        density_kind="even",
        central_t=Fraction(1507867, 1000000),
        quartic_r=Fraction(3284591, 1000000),
        target_exponent=4289,
        bridge_max_chain=2143,
        verbose_bridge_terms=False,
    )


def d93_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D93 finite bridge and Chebyshev tail from n=4363 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D93",
        row_label="D_93",
        rank=93,
        exact_through=35,
        dimension=17205,
        kappa=184,
        c_value=91,
        radius_s=Fraction(342807625, 100000),
        exp_terms=6000,
        quartic_s=Fraction(3085228, 1000000),
        density_kind="odd",
        central_t=Fraction(1636008, 1000000),
        quartic_r=Fraction(212251, 1000000),
        target_exponent=4363,
        bridge_max_chain=2180,
        verbose_bridge_terms=False,
    )


def d94_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D94 finite bridge and Chebyshev tail from n=4497 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D94",
        row_label="D_94",
        rank=94,
        exact_through=35,
        dimension=17578,
        kappa=186,
        c_value=94,
        radius_s=Fraction(346590543, 100000),
        exp_terms=6000,
        quartic_s=Fraction(266427, 1000000),
        density_kind="even",
        central_t=Fraction(1550968, 1000000),
        quartic_r=Fraction(3100350, 1000000),
        target_exponent=4497,
        bridge_max_chain=2247,
        verbose_bridge_terms=False,
    )


def d95_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D95 finite bridge and Chebyshev tail from n=4567 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D95",
        row_label="D_95",
        rank=95,
        exact_through=35,
        dimension=17955,
        kappa=188,
        c_value=93,
        radius_s=Fraction(350386012, 100000),
        exp_terms=6000,
        quartic_s=Fraction(2711718, 1000000),
        density_kind="odd",
        central_t=Fraction(1576362, 1000000),
        quartic_r=Fraction(392194, 1000000),
        target_exponent=4567,
        bridge_max_chain=2282,
        verbose_bridge_terms=False,
    )


def d96_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D96 finite bridge and Chebyshev tail from n=4693 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D96",
        row_label="D_96",
        rank=96,
        exact_through=35,
        dimension=18336,
        kappa=190,
        c_value=96,
        radius_s=Fraction(354288486, 100000),
        exp_terms=6000,
        quartic_s=Fraction(274661, 1000000),
        density_kind="even",
        central_t=Fraction(1568734, 1000000),
        quartic_r=Fraction(2884871, 1000000),
        target_exponent=4693,
        bridge_max_chain=2345,
        verbose_bridge_terms=False,
    )


def d97_post_m29_finite_tail_certificate(delta_logs: list[Path]) -> None:
    """Replay the D97 finite bridge and Chebyshev tail from n=4779 onward."""

    d_high_rank_post_m29_finite_tail_certificate(
        delta_logs=delta_logs,
        label="D97",
        row_label="D_97",
        rank=97,
        exact_through=35,
        dimension=18721,
        kappa=192,
        c_value=95,
        radius_s=Fraction(358036174, 100000),
        exp_terms=6000,
        quartic_s=Fraction(2742965, 1000000),
        density_kind="odd",
        central_t=Fraction(1631155, 1000000),
        quartic_r=Fraction(201778, 1000000),
        target_exponent=4779,
        bridge_max_chain=2388,
        verbose_bridge_terms=False,
    )


def stable_edge_cap_ratio(data: dict[str, object], delta: float) -> float:
    """Return S/a for the scaled root-angle cap at the stable edge.

    In the identity-saddle scaling v=x/sqrt(n), the fixed cap
    |alpha(v)| <= delta for every positive root implies

        S = Q(x)/(2d) <= |R_+| * delta^2 * n/(2d).

    The Macdonald-Mehta radial variable S has gamma shape a=d/2.  Thus this
    S/a upper envelope measures how far into the radial mass any fixed cap can
    reach.
    """

    r = int(data["r"])
    d = int(data["d"])
    n = int(data["stable_last_odd"])
    r_plus = (d - r) // 2
    a = d / 2
    s_cap = r_plus * delta * delta * n / (2 * d)
    return s_cap / a


def stable_edge_cap_log10_upper(data: dict[str, object], delta: float) -> float:
    """Chernoff upper bound for fixed-cap Wronskian mass at the stable edge."""

    d = int(data["d"])
    r = int(data["r"])
    n = int(data["stable_last_odd"])
    r_plus = (d - r) // 2
    a = d / 2
    s_cap = r_plus * delta * delta * n / (2 * d)
    if s_cap <= 0:
        return float("-inf")
    # If X,Y are independent Gamma(a,1), then the capped Wronskian fraction is
    # at most S^2/a * P(X<=S)^2.  For S<a, Chernoff gives
    # P(X<=S) <= (e*S/a)^a.
    if s_cap >= a:
        return 0.0
    log_bound = (2 * log10(s_cap) - log10(a)) + 2 * a * (
        1 / log(10) + log10(s_cap / a)
    )
    return min(0.0, log_bound)


def normal_tail_log_lower(c_value: int) -> float:
    """Displayed Mills Gaussian lower bound from Proposition 24.14."""

    c = float(c_value)
    l_value = float(TRACE_L)
    lc_value = l_value * c
    return (
        -lc_value / 2
        + 0.5 * log(lc_value)
        - log(lc_value + 1)
        - 0.5 * log(2)
    )


def one_trace_tv_log_bound(n_value: int) -> float:
    """Displayed CJL/Stirling log bound for E_1(n)."""

    n = float(n_value)
    return log(5) + 0.25 * log(log(2 * n)) - n * log(2 * n / 2.718281828459045)


def s_eta_log(c_value: int) -> float:
    return float(sharp_tau_log(c_value))


def quad_ratio_logs(c_value: int) -> tuple[float, float, float]:
    c = float(c_value)
    r = float(TRACE_R)
    eta = float(TRACE_ETA)
    a_const = float(TRACE_A)
    den = 1 - eta ** -2 * c ** -2
    s_log = s_eta_log(c_value)
    r0 = a_const * c + s_log
    r1 = (
        log(144 * c * c * (r / 2) ** 3 / (r * r * den))
        + a_const * c
        + s_log
        + c * log(2 / r)
    )
    r2 = (
        log(4 * (r / eta) ** 3 / (r * r * c * c * den))
        - (log(r / eta) - a_const) * c
    )
    return r0, r1, r2


def first_threshold(predicate, start: int = 1, limit: int = 100_000) -> int:
    for value in range(start, limit + 1):
        if predicate(value):
            return value
    raise RuntimeError(f"threshold not found up to {limit}")


def trace_bridge_diagnostic() -> None:
    """Print the displayed high-rank trace-bridge cutoff diagnostics."""

    a_const = float(TRACE_A)
    c_gauss = first_threshold(
        lambda c: normal_tail_log_lower(c) >= log(2) - a_const * c
    )
    c_tv1 = first_threshold(
        lambda n: n >= 124 and one_trace_tv_log_bound(n) <= -a_const * n,
        start=124,
    )
    c_quad = first_threshold(
        lambda c: c > 1 / float(TRACE_ETA)
        and max(quad_ratio_logs(c)) <= log(0.5),
        start=124,
    )
    print("trace bridge constants")
    print(f"r={TRACE_R} eta={TRACE_ETA} A={TRACE_A} L={TRACE_L}")
    print(f"C_gauss <= {c_gauss}")
    print(f"C_TV1   <= {c_tv1}")
    print(f"C_quad  <= {c_quad}")
    print(f"monotone displayed cutoff <= {max(c_gauss, c_tv1, c_quad)}")
    print(f"retained-correction effective cutoff <= {TRACE_CUTOFF}")


def row(
    family: str, rank: int
) -> tuple[str, int, int, int, float, float, int, float]:
    data = family_data(family, rank)
    min_n, min_margin = post_stable_min_margin(data)
    return (
        family,
        rank,
        int(data["d"]),
        int(data["stable_last_odd"]),
        log_a0_lower(data) / log(10),
        stable_edge_margin(data) / log(10),
        min_n,
        min_margin / log(10),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rank", type=int, default=20)
    parser.add_argument(
        "--cap-diagnostic",
        action="store_true",
        help="also print fixed root-angle cap obstruction data at the stable edge",
    )
    parser.add_argument(
        "--delta",
        type=float,
        default=3.141592653589793,
        help="root-angle cap used by --cap-diagnostic",
    )
    parser.add_argument(
        "--trace-bridge-diagnostic",
        action="store_true",
        help="print the displayed high-rank trace-bridge cutoff diagnostics",
    )
    parser.add_argument(
        "--trace-cutoff-certificate",
        action="store_true",
        help="check the concrete displayed monotone trace cutoff margins",
    )
    parser.add_argument(
        "--trace-square-cutoff-certificate",
        action="store_true",
        help="check the Rains-square trace cutoff margins",
    )
    parser.add_argument(
        "--post-m29-local-onset-obstruction",
        action="store_true",
        help="check the B2/C2 half-leading obstruction at odd n=63",
    )
    parser.add_argument(
        "--b2-c2-post-m29-tail-certificate",
        action="store_true",
        help="check the B2/C2 direct post-m29 pushforward tail certificate",
    )
    parser.add_argument(
        "--b3-c3-post-m29-tail-certificate",
        action="store_true",
        help="check the B3/C3 direct post-m29 pushforward tail certificate",
    )
    parser.add_argument(
        "--d4-post-m29-tail-certificate",
        action="store_true",
        help="check the D4 direct post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--b4-c4-post-m29-tail-certificate",
        action="store_true",
        help="check the B4/C4 direct post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--d5-post-m29-tail-certificate",
        action="store_true",
        help="check the D5 direct post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--c5-post-m29-tail-certificate",
        action="store_true",
        help="check the C5 direct post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--b5-post-m29-tail-certificate",
        action="store_true",
        help="check the B5 direct post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--c6-post-m29-tail-certificate",
        action="store_true",
        help="check the C6 central-Y post-m29 ball-cap tail certificate",
    )
    parser.add_argument(
        "--d6-b6-post-m29-tail-certificate",
        action="store_true",
        help="check the D6/B6 product-sine post-m29 ball-cap certificates",
    )
    parser.add_argument(
        "--d7-b7-c7-post-m29-tail-certificate",
        action="store_true",
        help="check the D7/B7/C7 product-sine post-m29 ball-cap certificates",
    )
    parser.add_argument(
        "--c8-post-m29-tail-certificate",
        action="store_true",
        help="check the C8 product-sine post-m29 ball-cap certificate",
    )
    parser.add_argument(
        "--d8-post-m29-tail-certificate",
        action="store_true",
        help="check the D8 radial-sine post-m29 ball-cap certificate",
    )
    parser.add_argument(
        "--b8-c9-d9-post-m29-tail-certificate",
        action="store_true",
        help="check the B8/C9/D9 root-length post-m29 ball-cap certificates",
    )
    parser.add_argument(
        "--c10-post-m29-tail-certificate",
        action="store_true",
        help="check the C10 half-integer root-length post-m29 ball-cap certificate",
    )
    parser.add_argument(
        "--b9-post-m29-tail-certificate",
        action="store_true",
        help="check the B9 exponential-sine quartic-tail post-m29 certificate",
    )
    parser.add_argument(
        "--c11-post-m29-tail-certificate",
        action="store_true",
        help="check the C11 exponential-sine quartic-tail post-m29 certificate",
    )
    parser.add_argument(
        "--c12-post-m29-tail-certificate",
        action="store_true",
        help="check the C12 sharpened negative-tail post-m29 certificate",
    )
    parser.add_argument(
        "--b10-post-m29-tail-certificate",
        action="store_true",
        help="check the B10 Chain bridge and sharpened post-m29 certificate",
    )
    parser.add_argument(
        "--b11-post-m29-direct-tail-certificate",
        action="store_true",
        help="check the B11 sharpened direct tail from n=77 onward",
    )
    parser.add_argument(
        "--b11-post-m29-bridge-range-certificate",
        action="store_true",
        help="check the B11 finite Chain bridge from n=65 through n=75",
    )
    parser.add_argument(
        "--low-bc-exact-onset-bridge-certificate",
        action="store_true",
        help="check the finite Chain gaps preceding the fixed low-B/C exact-MGF onsets",
    )
    parser.add_argument(
        "--b12-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the B12 finite Chain bridge through n=79 and direct tail from n=81 onward",
    )
    parser.add_argument(
        "--b13-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the B13 finite Chain bridge through n=89 and direct tail from n=91 onward",
    )
    parser.add_argument(
        "--c13-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the C13 finite Chain bridge through n=69 and direct tail from n=71 onward",
    )
    parser.add_argument(
        "--c14-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the C14 finite Chain bridge through n=81 and direct tail from n=83 onward",
    )
    parser.add_argument(
        "--c15-post-m29-tail-certificate",
        action="store_true",
        help="check the C15 direct post-m29 Chebyshev tail from n=63 onward",
    )
    parser.add_argument(
        "--c16-post-m29-tail-certificate",
        action="store_true",
        help="check the C16 direct post-m29 Chebyshev tail from n=63 onward",
    )
    parser.add_argument(
        "--c17-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the C17 finite Chain bridge through n=71 and direct tail from n=73 onward",
    )
    parser.add_argument(
        "--c18-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the C18 finite Chain bridge through n=83 and direct tail from n=85 onward",
    )
    parser.add_argument(
        "--c19-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the C19 finite Chain bridge through n=97 and direct tail from n=99 onward",
    )
    parser.add_argument(
        "--d10-post-m29-tail-certificate",
        action="store_true",
        help="check the D10 Chebyshev sharpened post-m29 certificate",
    )
    parser.add_argument(
        "--d11-post-m29-tail-certificate",
        action="store_true",
        help="check the D11 Chebyshev sharpened post-m29 certificate",
    )
    parser.add_argument(
        "--d12-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D12 finite Chain bridge through n=77 and direct tail from n=79 onward",
    )
    parser.add_argument(
        "--d13-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D13 finite Chain bridge through n=199 and direct tail from n=201 onward",
    )
    parser.add_argument(
        "--d14-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D14 finite Chain bridge through n=323 and direct tail from n=325 onward",
    )
    parser.add_argument(
        "--d15-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D15 finite Chain bridge through n=523 and direct tail from n=525 onward",
    )
    parser.add_argument(
        "--d16-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D16 finite Chain bridge through n=273 and direct tail from n=275 onward",
    )
    parser.add_argument(
        "--d17-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D17 finite Chain bridge through n=285 and direct tail from n=287 onward",
    )
    parser.add_argument(
        "--d18-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D18 finite Chain bridge through n=363 and direct tail from n=365 onward",
    )
    parser.add_argument(
        "--d19-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D19 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d20-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D20 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d21-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D21 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d22-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D22 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d23-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D23 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d24-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D24 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d25-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D25 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d26-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D26 finite Chain bridge through n=403 and direct tail from n=405 onward",
    )
    parser.add_argument(
        "--d27-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D27 finite Chain bridge through n=411 and direct tail from n=413 onward",
    )
    parser.add_argument(
        "--d28-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D28 finite Chain bridge through n=431 and direct tail from n=433 onward",
    )
    parser.add_argument(
        "--d29-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D29 finite Chain bridge through n=447 and direct tail from n=449 onward",
    )
    parser.add_argument(
        "--d30-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D30 finite Chain bridge through n=479 and direct tail from n=481 onward",
    )
    parser.add_argument(
        "--d31-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D31 finite Chain bridge through n=495 and direct tail from n=497 onward",
    )
    parser.add_argument(
        "--d32-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D32 finite Chain bridge through n=511 and direct tail from n=513 onward",
    )
    parser.add_argument(
        "--d33-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D33 finite Chain bridge through n=543 and direct tail from n=545 onward",
    )
    parser.add_argument(
        "--d34-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D34 finite Chain bridge through n=575 and direct tail from n=577 onward",
    )
    parser.add_argument(
        "--d35-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D35 finite Chain bridge through n=591 and direct tail from n=593 onward",
    )
    parser.add_argument(
        "--d36-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D36 finite Chain bridge through n=623 and direct tail from n=625 onward",
    )
    parser.add_argument(
        "--d37-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D37 finite Chain bridge through n=655 and direct tail from n=657 onward",
    )
    parser.add_argument(
        "--d38-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D38 finite Chain bridge through n=671 and direct tail from n=673 onward",
    )
    parser.add_argument(
        "--d39-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D39 finite Chain bridge through n=703 and direct tail from n=705 onward",
    )
    parser.add_argument(
        "--d40-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D40 finite Chain bridge through n=735 and direct tail from n=737 onward",
    )
    parser.add_argument(
        "--d41-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D41 finite Chain bridge through n=767 and direct tail from n=769 onward",
    )
    parser.add_argument(
        "--d42-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D42 finite Chain bridge through n=799 and direct tail from n=801 onward",
    )
    parser.add_argument(
        "--d43-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D43 finite Chain bridge through n=831 and direct tail from n=833 onward",
    )
    parser.add_argument(
        "--d44-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D44 finite Chain bridge through n=879 and direct tail from n=881 onward",
    )
    parser.add_argument(
        "--d45-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D45 finite Chain bridge through n=911 and direct tail from n=913 onward",
    )
    parser.add_argument(
        "--d46-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D46 finite Chain bridge through n=975 and direct tail from n=977 onward",
    )
    parser.add_argument(
        "--d47-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D47 finite Chain bridge through n=1007 and direct tail from n=1009 onward",
    )
    parser.add_argument(
        "--d48-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D48 finite Chain bridge through n=1071 and direct tail from n=1073 onward",
    )
    parser.add_argument(
        "--d49-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D49 finite Chain bridge through n=1091 and direct tail from n=1093 onward",
    )
    parser.add_argument(
        "--d50-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D50 finite Chain bridge through n=1153 and direct tail from n=1155 onward",
    )
    parser.add_argument(
        "--d51-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D51 finite Chain bridge through n=1189 and direct tail from n=1191 onward",
    )
    parser.add_argument(
        "--d52-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D52 finite Chain bridge through n=1253 and direct tail from n=1255 onward",
    )
    parser.add_argument(
        "--d53-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D53 finite Chain bridge through n=1291 and direct tail from n=1293 onward",
    )
    parser.add_argument(
        "--d54-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D54 finite Chain bridge through n=1357 and direct tail from n=1359 onward",
    )
    parser.add_argument(
        "--d55-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D55 finite Chain bridge through n=1397 and direct tail from n=1399 onward",
    )
    parser.add_argument(
        "--d56-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D56 finite Chain bridge through n=1467 and direct tail from n=1469 onward",
    )
    parser.add_argument(
        "--d57-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D57 finite Chain bridge through n=1509 and direct tail from n=1511 onward",
    )
    parser.add_argument(
        "--d58-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D58 finite Chain bridge through n=1583 and direct tail from n=1585 onward",
    )
    parser.add_argument(
        "--d59-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D59 finite Chain bridge through n=1625 and direct tail from n=1627 onward",
    )
    parser.add_argument(
        "--d60-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D60 finite Chain bridge through n=1701 and direct tail from n=1703 onward",
    )
    parser.add_argument(
        "--d61-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D61 finite Chain bridge through n=1751 and direct tail from n=1753 onward",
    )
    parser.add_argument(
        "--d62-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D62 finite Chain bridge through n=1831 and direct tail from n=1833 onward",
    )
    parser.add_argument(
        "--d63-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D63 finite Chain bridge through n=1889 and direct tail from n=1891 onward",
    )
    parser.add_argument(
        "--d64-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D64 finite Chain bridge through n=1963 and direct tail from n=1965 onward",
    )
    parser.add_argument(
        "--d65-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D65 finite Chain bridge through n=2059 and direct tail from n=2061 onward",
    )
    parser.add_argument(
        "--d66-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D66 finite Chain bridge through m=1050 and direct tail from n=2101 onward",
    )
    parser.add_argument(
        "--d67-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D67 finite Chain bridge through n=2185 and direct tail from n=2187 onward",
    )
    parser.add_argument(
        "--d68-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D68 finite Chain bridge through n=2249 and direct tail from n=2251 onward",
    )
    parser.add_argument(
        "--d69-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D69 finite Chain bridge through n=2279 and direct tail from n=2281 onward",
    )
    parser.add_argument(
        "--d70-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D70 finite Chain bridge through n=2367 and direct tail from n=2369 onward",
    )
    parser.add_argument(
        "--d71-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D71 finite Chain bridge through n=2437 and direct tail from n=2439 onward",
    )
    parser.add_argument(
        "--d72-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D72 finite Chain bridge through n=2523 and direct tail from n=2525 onward",
    )
    parser.add_argument(
        "--d73-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D73 finite Chain bridge through n=2573 and direct tail from n=2575 onward",
    )
    parser.add_argument(
        "--d74-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D74 finite Chain bridge through n=2667 and direct tail from n=2669 onward",
    )
    parser.add_argument(
        "--d75-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D75 finite Chain bridge through n=2727 and direct tail from n=2729 onward",
    )
    parser.add_argument(
        "--d76-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D76 finite Chain bridge through n=2825 and direct tail from n=2827 onward",
    )
    parser.add_argument(
        "--d77-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D77 finite Chain bridge through n=2887 and direct tail from n=2889 onward",
    )
    parser.add_argument(
        "--d78-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D78 finite Chain bridge through n=2987 and direct tail from n=2989 onward",
    )
    parser.add_argument(
        "--d79-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D79 finite Chain bridge through n=3053 and direct tail from n=3055 onward",
    )
    parser.add_argument(
        "--d80-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D80 finite Chain bridge through n=3153 and direct tail from n=3155 onward",
    )
    parser.add_argument(
        "--d81-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D81 finite Chain bridge through n=3219 and direct tail from n=3221 onward",
    )
    parser.add_argument(
        "--d82-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D82 finite Chain bridge through n=3325 and direct tail from n=3327 onward",
    )
    parser.add_argument(
        "--d83-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D83 finite Chain bridge through n=3393 and direct tail from n=3395 onward",
    )
    parser.add_argument(
        "--d84-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D84 finite Chain bridge through n=3505 and direct tail from n=3507 onward",
    )
    parser.add_argument(
        "--d85-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D85 finite Chain bridge through n=3581 and direct tail from n=3583 onward",
    )
    parser.add_argument(
        "--d86-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D86 finite Chain bridge through n=3691 and direct tail from n=3693 onward",
    )
    parser.add_argument(
        "--d87-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D87 finite Chain bridge through n=3767 and direct tail from n=3769 onward",
    )
    parser.add_argument(
        "--d88-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D88 finite Chain bridge through n=3885 and direct tail from n=3887 onward",
    )
    parser.add_argument(
        "--d89-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D89 finite Chain bridge through n=3959 and direct tail from n=3961 onward",
    )
    parser.add_argument(
        "--d90-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D90 finite Chain bridge through n=4083 and direct tail from n=4085 onward",
    )
    parser.add_argument(
        "--d91-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D91 finite Chain bridge through n=4157 and direct tail from n=4159 onward",
    )
    parser.add_argument(
        "--d92-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D92 finite Chain bridge through n=4287 and direct tail from n=4289 onward",
    )
    parser.add_argument(
        "--d93-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D93 finite Chain bridge through n=4361 and direct tail from n=4363 onward",
    )
    parser.add_argument(
        "--d94-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D94 finite Chain bridge through n=4495 and direct tail from n=4497 onward",
    )
    parser.add_argument(
        "--d95-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D95 finite Chain bridge through n=4565 and direct tail from n=4567 onward",
    )
    parser.add_argument(
        "--d96-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D96 finite Chain bridge through n=4691 and direct tail from n=4693 onward",
    )
    parser.add_argument(
        "--d97-post-m29-finite-tail-certificate",
        action="store_true",
        help="check the D97 finite Chain bridge through n=4777 and direct tail from n=4779 onward",
    )
    parser.add_argument(
        "--d-m21-high-middle-linear-certificate",
        action="store_true",
        help="check the m=21 D-rank >= 27 linear monotonicity certificate",
    )
    parser.add_argument(
        "--d-m22-high-linear-certificate",
        action="store_true",
        help="check the m=22 D_29..D_49 linear monotonicity certificate",
    )
    parser.add_argument(
        "--d-m22-middle-lower-certificate",
        action="store_true",
        help="check the m=22 D_18..D_28 lower-bound certificate",
    )
    parser.add_argument(
        "--d-m22-low-middle-lower-certificate",
        action="store_true",
        help="check the m=22 D_12..D_17 lower-bound certificate",
    )
    parser.add_argument(
        "--d-m21-middle-lower-certificate",
        action="store_true",
        help="check the m=21 D_20..D_26 lower-bound certificate",
    )
    parser.add_argument(
        "--d-m21-low-middle-lower-certificate",
        action="store_true",
        help="check the m=21 D_12..D_18 lower-bound certificate",
    )
    parser.add_argument(
        "--d-m21-bottom-lower-certificate",
        action="store_true",
        help="check the m=21 D_9..D_11 lower-bound certificate",
    )
    parser.add_argument(
        "--bc-m21-high-boundary-lower-certificate",
        action="store_true",
        help="check the m=21 B/C rank 9..16 boundary lower bound",
    )
    parser.add_argument(
        "--bc-m22-boundary-lower-certificate",
        action="store_true",
        help="check the m=22 B/C rank 9..23 boundary lower bound",
    )
    parser.add_argument(
        "--bc-m22-low-boundary-lower-certificate",
        action="store_true",
        help="check the m=22 B/C rank 6..8 boundary lower bound",
    )
    parser.add_argument(
        "--bc-m25-reused-window-certificate",
        action="store_true",
        help="check m=25 B/C_5..B/C_16 using previously certified exact deltas",
    )
    parser.add_argument(
        "--bc-reused-window-certificate",
        action="store_true",
        help="check a B/C rank range using previously certified exact deltas",
    )
    parser.add_argument(
        "--bc-reused-window-stable-quadratic-certificate",
        action="store_true",
        help=(
            "check a B/C rank range using exact deltas and stable bounds for "
            "outside negative quadratic terms"
        ),
    )
    parser.add_argument(
        "--bc-window-delta-log",
        type=Path,
        action="append",
        default=[],
        help="C++ exact-window log supplying B/C Delta lines",
    )
    parser.add_argument(
        "--d-m25-reused-delta-certificate",
        action="store_true",
        help="check m=25 D_4..D_13 using previously certified determinant deltas",
    )
    parser.add_argument(
        "--d-reused-delta-certificate",
        action="store_true",
        help="check a D rank range using previously certified determinant deltas",
    )
    parser.add_argument(
        "--d-delta-log",
        type=Path,
        action="append",
        default=[],
        help="C++ determinant-delta log supplying D Delta lines",
    )
    parser.add_argument(
        "--reused-chain-m",
        type=int,
        default=25,
        help="Chain step used by the generic reused-delta/window certificates",
    )
    parser.add_argument(
        "--bc-reused-low-rank",
        type=int,
        default=5,
        help="first B/C rank for --bc-reused-window-certificate",
    )
    parser.add_argument(
        "--bc-reused-high-rank",
        type=int,
        default=16,
        help="last B/C rank for --bc-reused-window-certificate",
    )
    parser.add_argument(
        "--d-reused-low-rank",
        type=int,
        default=4,
        help="first D rank for --d-reused-delta-certificate",
    )
    parser.add_argument(
        "--d-reused-high-rank",
        type=int,
        default=13,
        help="last D rank for --d-reused-delta-certificate",
    )
    parser.add_argument(
        "--finite-residue-plan",
        action="store_true",
        help="print the finite B/C/D residue induced by the trace cutoff",
    )
    parser.add_argument(
        "--finite-leading-certificate",
        action="store_true",
        help="check leading-term dominance over the finite residue from the stable edge",
    )
    parser.add_argument(
        "--middle-slice-plan",
        action="store_true",
        help="print the pre-Wong row gate for the next exact Chain steps",
    )
    parser.add_argument(
        "--finite-trace-tail-diagnostic",
        action="store_true",
        help="test finite-row coverage by the Rains-square trace route",
    )
    parser.add_argument(
        "--first-hit-onset-plan",
        action="store_true",
        help="print the row-wise N_loc target equivalent to no later middle steps",
    )
    parser.add_argument(
        "--first-hit-trace-diagnostic",
        action="store_true",
        help="test trace-pushforward closure at each first-hit exponent",
    )
    parser.add_argument(
        "--post-m29-chernoff-trace-certificate",
        action="store_true",
        help="check the post-m=29 finite rows closed by the CJL Chernoff tau supplier",
    )
    parser.add_argument(
        "--post-m29-lower-bound-supplier-certificate",
        action="store_true",
        help="check the first-hit lower-bound certificate on the 796-row residue",
    )
    parser.add_argument(
        "--post-m29-bc-lower-bound-supplier-certificate",
        action="store_true",
        help="check only the accepted 530-row B/C first-hit lower-bound certificate",
    )
    parser.add_argument(
        "--d-boundary-m15-certificate",
        action="store_true",
        help="replay the D_35 Pfaffian boundary correction at Chain step m=15",
    )
    parser.add_argument(
        "--d-boundary-m15-depth2-certificate",
        action="store_true",
        help="replay the D_35, D_34, and D_33 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth3-certificate",
        action="store_true",
        help="replay the D_35 through D_32 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth4-certificate",
        action="store_true",
        help="replay the D_35 through D_31 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth8-certificate",
        action="store_true",
        help="replay the D_35 through D_27 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth11-certificate",
        action="store_true",
        help="replay the D_35 through D_24 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth15-certificate",
        action="store_true",
        help="replay the D_23 through D_20 determinant-boundary corrections",
    )
    parser.add_argument(
        "--d-boundary-m15-depth16-certificate",
        action="store_true",
        help="replay the D_19 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth17-certificate",
        action="store_true",
        help="replay the D_18 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth18-certificate",
        action="store_true",
        help="replay the D_17 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth19-certificate",
        action="store_true",
        help="replay the D_16 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth20-certificate",
        action="store_true",
        help="replay the D_15 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth21-certificate",
        action="store_true",
        help="replay the D_14 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth22-certificate",
        action="store_true",
        help="replay the D_13 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth23-certificate",
        action="store_true",
        help="replay the D_12 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth24-certificate",
        action="store_true",
        help="replay the D_11 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth25-certificate",
        action="store_true",
        help="replay the D_10 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth26-certificate",
        action="store_true",
        help="replay the D_9 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth27-certificate",
        action="store_true",
        help="replay the D_8 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth28-certificate",
        action="store_true",
        help="replay the D_7 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth29-certificate",
        action="store_true",
        help="replay the D_6 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth30-certificate",
        action="store_true",
        help="replay the D_5 determinant-boundary correction",
    )
    parser.add_argument(
        "--d-boundary-m15-depth31-certificate",
        action="store_true",
        help="replay the D_4 determinant-boundary correction",
    )
    parser.add_argument(
        "--b-boundary-m15-rank16-certificate",
        action="store_true",
        help="replay the B_16 odd-orthogonal first-boundary correction",
    )
    parser.add_argument(
        "--c-boundary-m15-rank16-certificate",
        action="store_true",
        help="replay the C_16 symplectic first-boundary correction",
    )
    parser.add_argument(
        "--trace-cutoff",
        type=int,
        default=TRACE_MONOTONE_CUTOFF,
        help="C_G cutoff used by --trace-cutoff-certificate",
    )
    parser.add_argument(
        "--finite-cutoff",
        type=int,
        default=TRACE_CUTOFF,
        help="effective C_G cutoff used by finite-residue diagnostics",
    )
    parser.add_argument(
        "--current-max-chain",
        type=int,
        default=14,
        help="current exact bridge depth used by --finite-residue-plan",
    )
    parser.add_argument(
        "--slice-steps",
        type=int,
        default=3,
        help="number of next Chain steps shown by --middle-slice-plan",
    )
    args = parser.parse_args()

    if args.trace_bridge_diagnostic:
        trace_bridge_diagnostic()
        return
    if args.trace_cutoff_certificate:
        trace_cutoff_certificate(args.trace_cutoff)
        return
    if args.trace_square_cutoff_certificate:
        trace_square_cutoff_certificate()
        return
    if args.post_m29_local_onset_obstruction:
        post_m29_local_onset_obstruction()
        return
    if args.b2_c2_post_m29_tail_certificate:
        b2_c2_post_m29_tail_certificate()
        return
    if args.b3_c3_post_m29_tail_certificate:
        b3_c3_post_m29_tail_certificate()
        return
    if args.d4_post_m29_tail_certificate:
        d4_post_m29_tail_certificate()
        return
    if args.b4_c4_post_m29_tail_certificate:
        b4_c4_post_m29_tail_certificate()
        return
    if args.d5_post_m29_tail_certificate:
        d5_post_m29_tail_certificate()
        return
    if args.c5_post_m29_tail_certificate:
        c5_post_m29_tail_certificate()
        return
    if args.b5_post_m29_tail_certificate:
        b5_post_m29_tail_certificate()
        return
    if args.c6_post_m29_tail_certificate:
        c6_post_m29_tail_certificate()
        return
    if args.d6_b6_post_m29_tail_certificate:
        d6_b6_post_m29_tail_certificate()
        return
    if args.d7_b7_c7_post_m29_tail_certificate:
        d7_b7_c7_post_m29_tail_certificate()
        return
    if args.c8_post_m29_tail_certificate:
        c8_post_m29_tail_certificate()
        return
    if args.d8_post_m29_tail_certificate:
        d8_post_m29_tail_certificate()
        return
    if args.b8_c9_d9_post_m29_tail_certificate:
        b8_c9_d9_post_m29_tail_certificate()
        return
    if args.c10_post_m29_tail_certificate:
        c10_post_m29_tail_certificate()
        return
    if args.b9_post_m29_tail_certificate:
        b9_post_m29_tail_certificate()
        return
    if args.c11_post_m29_tail_certificate:
        c11_post_m29_tail_certificate()
        return
    if args.c12_post_m29_tail_certificate:
        c12_post_m29_tail_certificate()
        return
    if args.b10_post_m29_tail_certificate:
        b10_post_m29_tail_certificate(args.bc_window_delta_log)
        return
    if args.b11_post_m29_direct_tail_certificate:
        b11_post_m29_direct_tail_certificate()
        return
    if args.b11_post_m29_bridge_range_certificate:
        b11_post_m29_bridge_range_certificate(args.bc_window_delta_log)
        return
    if args.low_bc_exact_onset_bridge_certificate:
        low_bc_exact_onset_bridge_certificate(args.bc_window_delta_log)
        return
    if args.b12_post_m29_finite_tail_certificate:
        b12_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.b13_post_m29_finite_tail_certificate:
        b13_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.c13_post_m29_finite_tail_certificate:
        c13_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.c14_post_m29_finite_tail_certificate:
        c14_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.c15_post_m29_tail_certificate:
        c15_post_m29_tail_certificate(args.bc_window_delta_log)
        return
    if args.c16_post_m29_tail_certificate:
        c16_post_m29_tail_certificate(args.bc_window_delta_log)
        return
    if args.c17_post_m29_finite_tail_certificate:
        c17_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.c18_post_m29_finite_tail_certificate:
        c18_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.c19_post_m29_finite_tail_certificate:
        c19_post_m29_finite_tail_certificate(args.bc_window_delta_log)
        return
    if args.d10_post_m29_tail_certificate:
        d10_post_m29_tail_certificate(args.d_delta_log)
        return
    if args.d11_post_m29_tail_certificate:
        d11_post_m29_tail_certificate(args.d_delta_log)
        return
    if args.d12_post_m29_finite_tail_certificate:
        d12_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d13_post_m29_finite_tail_certificate:
        d13_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d14_post_m29_finite_tail_certificate:
        d14_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d15_post_m29_finite_tail_certificate:
        d15_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d16_post_m29_finite_tail_certificate:
        d16_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d17_post_m29_finite_tail_certificate:
        d17_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d18_post_m29_finite_tail_certificate:
        d18_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d19_post_m29_finite_tail_certificate:
        d19_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d20_post_m29_finite_tail_certificate:
        d20_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d21_post_m29_finite_tail_certificate:
        d21_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d22_post_m29_finite_tail_certificate:
        d22_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d23_post_m29_finite_tail_certificate:
        d23_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d24_post_m29_finite_tail_certificate:
        d24_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d25_post_m29_finite_tail_certificate:
        d25_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d26_post_m29_finite_tail_certificate:
        d26_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d27_post_m29_finite_tail_certificate:
        d27_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d28_post_m29_finite_tail_certificate:
        d28_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d29_post_m29_finite_tail_certificate:
        d29_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d30_post_m29_finite_tail_certificate:
        d30_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d31_post_m29_finite_tail_certificate:
        d31_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d32_post_m29_finite_tail_certificate:
        d32_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d33_post_m29_finite_tail_certificate:
        d33_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d34_post_m29_finite_tail_certificate:
        d34_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d35_post_m29_finite_tail_certificate:
        d35_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d36_post_m29_finite_tail_certificate:
        d36_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d37_post_m29_finite_tail_certificate:
        d37_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d38_post_m29_finite_tail_certificate:
        d38_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d39_post_m29_finite_tail_certificate:
        d39_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d40_post_m29_finite_tail_certificate:
        d40_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d41_post_m29_finite_tail_certificate:
        d41_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d42_post_m29_finite_tail_certificate:
        d42_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d43_post_m29_finite_tail_certificate:
        d43_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d44_post_m29_finite_tail_certificate:
        d44_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d45_post_m29_finite_tail_certificate:
        d45_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d46_post_m29_finite_tail_certificate:
        d46_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d47_post_m29_finite_tail_certificate:
        d47_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d48_post_m29_finite_tail_certificate:
        d48_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d49_post_m29_finite_tail_certificate:
        d49_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d50_post_m29_finite_tail_certificate:
        d50_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d51_post_m29_finite_tail_certificate:
        d51_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d52_post_m29_finite_tail_certificate:
        d52_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d53_post_m29_finite_tail_certificate:
        d53_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d54_post_m29_finite_tail_certificate:
        d54_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d55_post_m29_finite_tail_certificate:
        d55_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d56_post_m29_finite_tail_certificate:
        d56_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d57_post_m29_finite_tail_certificate:
        d57_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d58_post_m29_finite_tail_certificate:
        d58_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d59_post_m29_finite_tail_certificate:
        d59_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d60_post_m29_finite_tail_certificate:
        d60_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d61_post_m29_finite_tail_certificate:
        d61_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d62_post_m29_finite_tail_certificate:
        d62_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d63_post_m29_finite_tail_certificate:
        d63_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d64_post_m29_finite_tail_certificate:
        d64_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d65_post_m29_finite_tail_certificate:
        d65_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d66_post_m29_finite_tail_certificate:
        d66_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d67_post_m29_finite_tail_certificate:
        d67_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d68_post_m29_finite_tail_certificate:
        d68_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d69_post_m29_finite_tail_certificate:
        d69_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d70_post_m29_finite_tail_certificate:
        d70_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d71_post_m29_finite_tail_certificate:
        d71_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d72_post_m29_finite_tail_certificate:
        d72_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d73_post_m29_finite_tail_certificate:
        d73_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d74_post_m29_finite_tail_certificate:
        d74_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d75_post_m29_finite_tail_certificate:
        d75_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d76_post_m29_finite_tail_certificate:
        d76_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d77_post_m29_finite_tail_certificate:
        d77_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d78_post_m29_finite_tail_certificate:
        d78_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d79_post_m29_finite_tail_certificate:
        d79_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d80_post_m29_finite_tail_certificate:
        d80_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d81_post_m29_finite_tail_certificate:
        d81_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d82_post_m29_finite_tail_certificate:
        d82_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d83_post_m29_finite_tail_certificate:
        d83_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d84_post_m29_finite_tail_certificate:
        d84_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d85_post_m29_finite_tail_certificate:
        d85_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d86_post_m29_finite_tail_certificate:
        d86_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d87_post_m29_finite_tail_certificate:
        d87_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d88_post_m29_finite_tail_certificate:
        d88_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d89_post_m29_finite_tail_certificate:
        d89_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d90_post_m29_finite_tail_certificate:
        d90_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d91_post_m29_finite_tail_certificate:
        d91_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d92_post_m29_finite_tail_certificate:
        d92_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d93_post_m29_finite_tail_certificate:
        d93_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d94_post_m29_finite_tail_certificate:
        d94_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d95_post_m29_finite_tail_certificate:
        d95_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d96_post_m29_finite_tail_certificate:
        d96_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d97_post_m29_finite_tail_certificate:
        d97_post_m29_finite_tail_certificate(args.d_delta_log)
        return
    if args.d_m21_high_middle_linear_certificate:
        d_m21_high_middle_linear_certificate()
        return
    if args.d_m22_high_linear_certificate:
        d_m22_high_linear_certificate()
        return
    if args.d_m22_middle_lower_certificate:
        d_m22_middle_lower_certificate()
        return
    if args.d_m22_low_middle_lower_certificate:
        d_m22_low_middle_lower_certificate()
        return
    if args.d_m21_middle_lower_certificate:
        d_m21_middle_lower_certificate()
        return
    if args.d_m21_low_middle_lower_certificate:
        d_m21_low_middle_lower_certificate()
        return
    if args.d_m21_bottom_lower_certificate:
        d_m21_bottom_lower_certificate()
        return
    if args.bc_m21_high_boundary_lower_certificate:
        bc_m21_high_boundary_lower_certificate()
        return
    if args.bc_m22_boundary_lower_certificate:
        bc_m22_boundary_lower_certificate()
        return
    if args.bc_m22_low_boundary_lower_certificate:
        bc_m22_low_boundary_lower_certificate()
        return
    if args.bc_m25_reused_window_certificate:
        bc_m25_reused_window_certificate(args.bc_window_delta_log)
        return
    if args.bc_reused_window_certificate:
        bc_reused_window_certificate(
            args.bc_window_delta_log,
            args.reused_chain_m,
            args.bc_reused_low_rank,
            args.bc_reused_high_rank,
        )
        return
    if args.bc_reused_window_stable_quadratic_certificate:
        bc_reused_window_certificate(
            args.bc_window_delta_log,
            args.reused_chain_m,
            args.bc_reused_low_rank,
            args.bc_reused_high_rank,
            stable_quadratic_tail=True,
        )
        return
    if args.d_m25_reused_delta_certificate:
        d_m25_reused_delta_certificate(args.d_delta_log)
        return
    if args.d_reused_delta_certificate:
        d_reused_delta_certificate(
            args.d_delta_log,
            args.reused_chain_m,
            args.d_reused_low_rank,
            args.d_reused_high_rank,
        )
        return
    if args.finite_residue_plan:
        finite_residue_plan(args.finite_cutoff, args.current_max_chain)
        return
    if args.middle_slice_plan:
        middle_slice_plan(args.finite_cutoff, args.current_max_chain, args.slice_steps)
        return
    if args.first_hit_onset_plan:
        first_hit_onset_plan(args.finite_cutoff, args.current_max_chain)
        return
    if args.first_hit_trace_diagnostic:
        first_hit_trace_diagnostic(args.finite_cutoff, args.current_max_chain)
        return
    if args.post_m29_chernoff_trace_certificate:
        post_m29_chernoff_trace_certificate(args.finite_cutoff)
        return
    if args.post_m29_bc_lower_bound_supplier_certificate:
        post_m29_lower_bound_supplier_certificate(
            args.bc_window_delta_log,
            [],
            args.finite_cutoff,
            bc_only=True,
        )
        return
    if args.post_m29_lower_bound_supplier_certificate:
        post_m29_lower_bound_supplier_certificate(
            args.bc_window_delta_log,
            args.d_delta_log,
            args.finite_cutoff,
        )
        return
    if args.finite_trace_tail_diagnostic:
        finite_trace_tail_diagnostic(args.finite_cutoff, args.current_max_chain, args.slice_steps)
        return
    if args.d_boundary_m15_certificate:
        d_boundary_m15_certificate()
        return
    if args.d_boundary_m15_depth2_certificate:
        d_boundary_m15_depth2_certificate()
        return
    if args.d_boundary_m15_depth3_certificate:
        d_boundary_m15_depth3_certificate()
        return
    if args.d_boundary_m15_depth4_certificate:
        d_boundary_m15_depth4_certificate()
        return
    if args.d_boundary_m15_depth8_certificate:
        d_boundary_m15_depth8_certificate()
        return
    if args.d_boundary_m15_depth11_certificate:
        d_boundary_m15_depth11_certificate()
        return
    if args.d_boundary_m15_depth15_certificate:
        d_boundary_m15_depth15_certificate()
        return
    if args.d_boundary_m15_depth16_certificate:
        d_boundary_m15_depth16_certificate()
        return
    if args.d_boundary_m15_depth17_certificate:
        d_boundary_m15_depth17_certificate()
        return
    if args.d_boundary_m15_depth18_certificate:
        d_boundary_m15_depth18_certificate()
        return
    if args.d_boundary_m15_depth19_certificate:
        d_boundary_m15_depth19_certificate()
        return
    if args.d_boundary_m15_depth20_certificate:
        d_boundary_m15_depth20_certificate()
        return
    if args.d_boundary_m15_depth21_certificate:
        d_boundary_m15_depth21_certificate()
        return
    if args.d_boundary_m15_depth22_certificate:
        d_boundary_m15_depth22_certificate()
        return
    if args.d_boundary_m15_depth23_certificate:
        d_boundary_m15_depth23_certificate()
        return
    if args.d_boundary_m15_depth24_certificate:
        d_boundary_m15_depth24_certificate()
        return
    if args.d_boundary_m15_depth25_certificate:
        d_boundary_m15_depth25_certificate()
        return
    if args.d_boundary_m15_depth26_certificate:
        d_boundary_m15_depth26_certificate()
        return
    if args.d_boundary_m15_depth27_certificate:
        d_boundary_m15_depth27_certificate()
        return
    if args.d_boundary_m15_depth28_certificate:
        d_boundary_m15_depth28_certificate()
        return
    if args.d_boundary_m15_depth29_certificate:
        d_boundary_m15_depth29_certificate()
        return
    if args.d_boundary_m15_depth30_certificate:
        d_boundary_m15_depth30_certificate()
        return
    if args.d_boundary_m15_depth31_certificate:
        d_boundary_m15_depth31_certificate()
        return
    if args.b_boundary_m15_rank16_certificate:
        b_boundary_m15_rank16_certificate()
        return
    if args.c_boundary_m15_rank16_certificate:
        c_boundary_m15_rank16_certificate()
        return
    if args.finite_leading_certificate:
        finite_leading_margin_certificate(args.finite_cutoff)
        return

    print(
        f"{'fam':<4}{'rank':>6}{'d':>8}{'edge n':>8}"
        f"{'log10 A0 lower':>18}{'edge margin':>16}"
        f"{'min n':>8}{'min margin':>14}"
    )
    for family in ("B", "C", "D"):
        start = 4 if family == "D" else 2
        for rank in range(start, args.max_rank + 1):
            fam, b, d, n, log_a0, margin, min_n, min_margin = row(family, rank)
            print(
                f"{fam:<4}{b:>6}{d:>8}{n:>8}"
                f"{log_a0:>18.2f}{margin:>16.2f}"
                f"{min_n:>8}{min_margin:>14.2f}"
            )

    if args.cap_diagnostic:
        print()
        print(
            f"fixed root-angle cap diagnostic at delta={args.delta:.6g}"
        )
        print(
            f"{'fam':<4}{'rank':>6}{'edge n':>8}"
            f"{'S_root/a':>14}{'log10 frac upper':>20}"
        )
        for family in ("B", "C", "D"):
            start = 4 if family == "D" else 2
            for rank in range(start, args.max_rank + 1):
                data = family_data(family, rank)
                print(
                    f"{family:<4}{rank:>6}{int(data['stable_last_odd']):>8}"
                    f"{stable_edge_cap_ratio(data, args.delta):>14.4e}"
                    f"{stable_edge_cap_log10_upper(data, args.delta):>20.2f}"
                )


if __name__ == "__main__":
    main()
