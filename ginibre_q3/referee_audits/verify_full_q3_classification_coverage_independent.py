#!/usr/bin/env python3
"""Independent family/rank/degree coverage audit for the Part III theorem.

This is an author structural self-audit.  It does not verify the arithmetic
inside any certificate.  Instead it parses the theorem statements and the
independent low-row ledger, reconstructs every case split used by the final
classification proof, and proves that no nonautomatic odd-total-degree case
is omitted.  Infinite tails are represented symbolically by their upward
closed onset; finite residual cardinalities are recomputed from scratch.
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "full_q3_extension.tex"
LOW_ROW_LEDGER = ROOT / "character_ring_iter" / "full_q3_bcd_remaining_data.hpp"
RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")
RESULT_PATTERN = re.compile(
    r"\\begin\{(" + "|".join(RESULT_ENVS) + r")\}(.*?)\\end\{\1\}",
    re.DOTALL,
)


class AuditFailure(RuntimeError):
    """A fail-closed audit rejection."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditFailure(message)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def require_compact(text: str, fragment: str, message: str) -> None:
    require(compact(fragment) in compact(text), message)


def statement_with_label(text: str, label: str) -> str:
    needle = rf"\label{{{label}}}"
    matches = [match.group(2) for match in RESULT_PATTERN.finditer(text) if needle in match.group(2)]
    require(len(matches) == 1, f"expected one statement labelled {label}, found {len(matches)}")
    return matches[0]


def result_and_proof_with_label(text: str, label: str) -> str:
    needle = rf"\label{{{label}}}"
    matches = [match for match in RESULT_PATTERN.finditer(text) if needle in match.group(2)]
    require(len(matches) == 1, f"expected one result labelled {label}, found {len(matches)}")
    result = matches[0]
    proof = re.match(r"\s*\\begin\{proof\}.*?\\end\{proof\}", text[result.end() :], re.DOTALL)
    require(proof is not None, f"result labelled {label} lacks an immediate proof")
    return result.group(0) + proof.group(0)


def odd_degrees(stop: int, start: int = 5) -> range:
    require(start >= 5, f"invalid odd-degree start {start}")
    first = start if start % 2 == 1 else start + 1
    return range(first, stop + 1, 2)


def pair_count_at_degree(total_degree: int) -> int:
    """Count (a,n) with a>=2, n>=1 odd, and 2a+n=L."""

    require(total_degree >= 5 and total_degree % 2 == 1, f"invalid total degree {total_degree}")
    return (total_degree - 3) // 2


def residual_pair_count(stable_endpoint: int, tail_onset: int) -> int:
    return sum(
        pair_count_at_degree(degree)
        for degree in odd_degrees(tail_onset - 2)
        if degree > stable_endpoint
    )


def below_tail_pair_count(tail_onset: int) -> int:
    return sum(pair_count_at_degree(degree) for degree in odd_degrees(tail_onset - 2))


@dataclass(frozen=True)
class LowRow:
    family: str
    rank: int
    tail_onset: int
    moment_through: int
    method: str

    @property
    def stable_endpoint(self) -> int:
        return 2 * self.rank + 1 if self.family in {"B", "C"} else self.rank - 1


def parse_low_row_ledger(text: str) -> dict[tuple[str, int], LowRow]:
    pattern = re.compile(
        r"\{'([BCD])',\s*(\d+),\s*(\d+),\s*(\d+),\s*"
        r"TailMethod::([a-z_]+)\}"
    )
    rows: dict[tuple[str, int], LowRow] = {}
    for family, rank_text, onset_text, moment_text, method in pattern.findall(text):
        row = LowRow(family, int(rank_text), int(onset_text), int(moment_text), method)
        key = (row.family, row.rank)
        require(key not in rows, f"duplicate low-row ledger entry {key}")
        rows[key] = row
    return rows


def parse_exact_low_tail_rows(statement: str) -> dict[tuple[str, int], tuple[int, str]]:
    rows: dict[tuple[str, int], tuple[int, str]] = {}
    for raw_line in statement.splitlines():
        line = raw_line.strip()
        families = re.findall(r"([BCD])_(\d+)", line)
        if not families or "&" not in line:
            continue
        table_part = line.split(r"\\", 1)[0]
        fields = [field.strip() for field in table_part.split("&")]
        if r"\text{polynomial}" in table_part:
            require(len(fields) >= 6, f"malformed polynomial row: {line}")
            onset_text = fields[-2]
            method = "polynomial"
        else:
            require(len(fields) >= 5, f"malformed rational-cap row: {line}")
            onset_text = fields[-1]
            method = "rational_cap"
        require(re.fullmatch(r"\d+", onset_text) is not None, f"invalid tail onset in {line}")
        onset = int(onset_text)
        for family, rank_text in families:
            key = (family, int(rank_text))
            require(key not in rows, f"duplicate exact-tail statement row {key}")
            rows[key] = (onset, method)
    return rows


def parse_directed_low_tail_rows(statement: str) -> dict[tuple[str, int], tuple[int, str]]:
    pattern = re.compile(
        r"^b&(?P<ranks>[0-9&]+)\\\\[ \t]*\\hline[ \t]*\n"
        r"^K_b\((?P<family>[BCD])\)&(?P<onsets>[0-9&]+)\.?[ \t]*$",
        re.MULTILINE,
    )
    rows: dict[tuple[str, int], tuple[int, str]] = {}
    tables = list(pattern.finditer(statement))
    require(len(tables) == 5, f"expected five directed-tail tables, found {len(tables)}")
    for table in tables:
        family = table.group("family")
        ranks = [int(value) for value in table.group("ranks").split("&")]
        onsets = [int(value) for value in table.group("onsets").split("&")]
        require(len(ranks) == len(onsets), f"directed-tail row length mismatch in type {family}")
        for rank, onset in zip(ranks, onsets, strict=True):
            key = (family, rank)
            require(key not in rows, f"duplicate directed-tail statement row {key}")
            rows[key] = (onset, "directed_interval")
    return rows


def expected_low_rank_keys() -> set[tuple[str, int]]:
    return (
        {("B", rank) for rank in range(2, 18)}
        | {("C", rank) for rank in range(2, 29)}
        | {("D", rank) for rank in range(4, 31)}
    )


def validate_low_rows(
    rows: dict[tuple[str, int], LowRow],
    paper_tails: dict[tuple[str, int], tuple[int, str]],
) -> int:
    require(set(rows) == expected_low_rank_keys(), "low-row rank ledger is not exactly B2--17/C2--28/D4--30")
    require(set(paper_tails) == set(rows), "paper low-tail tables and row ledger have different scopes")
    method_counts = {method: 0 for method in ("polynomial", "rational_cap", "directed_interval")}
    residual_total = 0
    for key, row in rows.items():
        require(row.family in {"B", "C", "D"}, f"unknown family in {key}")
        require(row.tail_onset % 2 == 1, f"even low-row tail onset in {key}")
        require(row.tail_onset > row.stable_endpoint, f"tail begins before stable endpoint in {key}")
        require(
            row.moment_through >= row.tail_onset - 2,
            f"moment ledger stops before the residual box in {key}",
        )
        require(row.method in method_counts, f"unknown low-row tail method in {key}")
        require(
            paper_tails[key] == (row.tail_onset, row.method),
            f"paper/header tail mismatch in {key}: paper={paper_tails[key]}, header={row}",
        )
        method_counts[row.method] += 1
        residual_total += residual_pair_count(row.stable_endpoint, row.tail_onset)

        # The three theorem intervals are symbolic: stable is downward closed,
        # tail is upward closed, and the finite strict interval is residual.
        for degree in odd_degrees(row.tail_onset + 20):
            owners = {
                name
                for name, applies in (
                    ("stable", degree <= row.stable_endpoint),
                    ("residual", row.stable_endpoint < degree < row.tail_onset),
                    ("tail", degree >= row.tail_onset),
                )
                if applies
            }
            require(
                len(owners) == 1,
                f"low-row degree partition {key}, L={degree}: {owners}",
            )

    require(
        method_counts
        == {"polynomial": 6, "rational_cap": 6, "directed_interval": 58},
        f"low-row method counts={method_counts}",
    )
    require(residual_total == 12993, f"low-row residual pair count={residual_total}")
    return residual_total


def parse_exceptional_rows(text: str) -> dict[str, tuple[int, int, int]]:
    expected_groups = {"G_2", "F_4", "E_6", "E_7", "E_8"}
    tail_rows: dict[str, int] = {}
    tail_pattern = re.compile(
        r"^\s*\$(G_2|F_4|E_6|E_7|E_8)\$&([^\n\\]+?)(?:\\\\)?\s*$",
        re.MULTILINE,
    )
    for group, fields_text in tail_pattern.findall(text):
        fields = [field.strip() for field in fields_text.split("&")]
        require(len(fields) == 6 and fields[-1].isdigit(), f"malformed exceptional tail row {group}")
        require(group not in tail_rows, f"duplicate exceptional tail row {group}")
        tail_rows[group] = int(fields[-1])
    require(set(tail_rows) == expected_groups, f"exceptional tail groups={sorted(tail_rows)}")

    residual_statement = statement_with_label(text, "prop:exceptional-full-residual")
    residual_pattern = re.compile(
        r"^\s*(G_2|F_4|E_6|E_7|E_8)&(\d+)&(\d+)&(\d+)(?:\\\\)?\s*$",
        re.MULTILINE,
    )
    residual_rows = {
        group: (int(largest), int(count), int(minimum))
        for group, largest, count, minimum in residual_pattern.findall(residual_statement)
    }
    require(set(residual_rows) == expected_groups, f"exceptional residual groups={sorted(residual_rows)}")
    expected_onsets = {"G_2": 25, "F_4": 37, "E_6": 29, "E_7": 63, "E_8": 71}
    require(tail_rows == expected_onsets, f"exceptional tail onsets={tail_rows}")

    output: dict[str, tuple[int, int, int]] = {}
    for group in sorted(expected_groups):
        onset = tail_rows[group]
        largest, count, minimum = residual_rows[group]
        require(onset % 2 == 1 and largest == onset - 2, f"exceptional degree boundary mismatch in {group}")
        require(count == below_tail_pair_count(onset), f"exceptional residual count mismatch in {group}")
        require(minimum > 0, f"nonpositive stated exceptional minimum in {group}")
        output[group] = (onset, largest, count)
    require(sum(value[2] for value in output.values()) == 1265, "exceptional aggregate residual count mismatch")
    return output


def validate_explicit_type_a(text: str) -> int:
    cases = {
        "A1": ("lem:su2-effective-rectangle", "prop:su2-residual-certificate", 29, 78),
        "A2": ("prop:su3-cap-residual-certificate", "prop:su3-cap-residual-certificate", 29, 78),
        "A3": ("thm:su45-full-cone", "prop:su45-cap-residual-certificate", 27, 66),
        "A4": ("thm:su45-full-cone", "prop:su45-cap-residual-certificate", 37, 136),
    }
    residual_total = 0
    for family, (tail_label, residual_label, onset, count) in cases.items():
        tail_statement = statement_with_label(text, tail_label)
        tail_material = result_and_proof_with_label(text, tail_label)
        residual_statement = statement_with_label(text, residual_label)
        if family == "A1":
            require_compact(tail_statement, r"2a+n\geq29", "A1 tail boundary missing")
            require_compact(residual_statement, r"2a+n\leq27", "A1 residual boundary missing")
        elif family == "A2":
            require_compact(residual_statement, r"2a+n\leq27", "A2 residual boundary missing")
            theorem = result_and_proof_with_label(text, "thm:su3-full-cone")
            require_compact(theorem, r"totaldegree$2a+n\geq29$", "A2 tail boundary missing")
        elif family == "A3":
            require_compact(residual_statement, r"2a+n\leq25", "A3 residual boundary missing")
            require_compact(tail_material, r"atleast$27$", "A3 tail boundary missing")
        else:
            require_compact(residual_statement, r"2a+n\leq35", "A4 residual boundary missing")
            require_compact(tail_material, r"$2a+n\geq37$", "A4 tail boundary missing")
        computed = below_tail_pair_count(onset)
        require(computed == count, f"{family} residual count={computed}, expected={count}")
        residual_total += computed
        for degree in odd_degrees(onset + 20):
            owners = {"residual" if degree < onset else "tail"}
            require(len(owners) == 1, f"{family} degree gap at L={degree}")
    return residual_total


def type_a_tail_onset(n_value: int) -> int:
    if n_value == 6:
        return 27
    if n_value in {7, 8}:
        return 29
    return 31


def validate_general_type_a(text: str, residual_rank_maximum: int = 28) -> int:
    tail_statement = statement_with_label(text, "prop:sun-ge6-tail")
    residual_statement = statement_with_label(text, "prop:sun-ge6-residual")
    require_compact(tail_statement, r"27,&N=6", "type-A N=6 tail onset missing")
    require_compact(tail_statement, r"29,&N=7,8", "type-A N=7,8 tail onset missing")
    require_compact(tail_statement, r"31,&N\geq9", "type-A stable tail onset missing")
    require_compact(
        residual_statement,
        rf"6\leqN\leq{residual_rank_maximum}",
        "type-A residual rank interval mismatch",
    )
    require_compact(residual_statement, r"N<2a+n<K_N", "type-A residual degree interval missing")

    residual_total = 0
    gaps: list[tuple[int, int]] = []
    for n_value in range(6, 513):
        onset = type_a_tail_onset(n_value)
        for degree in odd_degrees(max(onset + 20, n_value + 20)):
            owners = set()
            if degree <= n_value:
                owners.add("stable")
            if degree >= onset:
                owners.add("tail")
            if 6 <= n_value <= residual_rank_maximum and n_value < degree < onset:
                owners.add("residual")
            if not owners:
                gaps.append((n_value, degree))
        if 6 <= n_value <= residual_rank_maximum:
            residual_total += sum(
                pair_count_at_degree(degree)
                for degree in odd_degrees(onset - 2)
                if degree > n_value
            )
    require(not gaps, f"type-A rank/degree gaps begin {gaps[:5]}")
    require(residual_total == 1315, f"type-A N>=6 residual count={residual_total}")
    # These are the only two boundary mechanisms not represented by a finite
    # residual row: parity closes N=29, and stable/tail adjacency closes N>=30.
    require(not [degree for degree in odd_degrees(29) if 29 < degree < 31], "N=29 parity boundary failed")
    for n_value in range(30, 513):
        require(
            not [degree for degree in odd_degrees(type_a_tail_onset(n_value) - 2) if degree > n_value],
            f"unexpected post-stable pre-tail degree at N={n_value}",
        )

    # This is the symbolic infinite-rank step behind the finite sanity loop:
    # for N>=30, an odd L>N has L>=31 and is therefore in the constant tail.
    require(type_a_tail_onset(29) == 31, "N=29 tail boundary changed")
    require(type_a_tail_onset(30) == 31, "N>=30 tail boundary changed")
    require(30 + 1 >= type_a_tail_onset(30), "N>=30 symbolic stable/tail adjacency failed")
    return residual_total


def parse_additional_tail_onsets(text: str) -> tuple[dict[int, int], dict[int, int]]:
    statement = statement_with_label(text, "prop:bcd-low-tail")
    normalized = compact(statement)
    b_match = re.search(
        r"\(K_\{18\},K_\{19\},K_\{20\},K_\{21\}\)=\((\d+),(\d+),(\d+),(\d+)\)",
        normalized,
    )
    require(b_match is not None, "cannot parse B18--B21 tail-onset tuple")
    b_values = [int(value) for value in b_match.groups()]
    b_onsets = dict(zip(range(18, 22), b_values, strict=True))

    d_match = re.search(
        r"\\begin\{split\}\s*\(\s*&?(.*?)\)\.\s*\\end\{split\}",
        statement,
        re.DOTALL,
    )
    require(d_match is not None, "cannot parse D31--D70 tail-onset array")
    d_values = [int(value) for value in re.findall(r"\d+", d_match.group(1))]
    require(len(d_values) == 40, f"D31--D70 tail-onset count={len(d_values)}")
    d_onsets = dict(zip(range(31, 71), d_values, strict=True))
    return b_onsets, d_onsets


def classical_negative_endpoint(family: str, rank: int) -> int:
    if family in {"B", "C"}:
        return rank
    require(family == "D", f"unknown classical family {family}")
    return rank if rank % 2 == 0 else rank - 2


def validate_finite_bd_rows(text: str) -> tuple[int, int]:
    b_onsets, d_onsets = parse_additional_tail_onsets(text)
    require(set(b_onsets) == set(range(18, 22)), "B finite-tail rank scope mismatch")
    require(set(d_onsets) == set(range(31, 71)), "D finite-tail rank scope mismatch")
    require(all(onset % 2 == 1 for onset in b_onsets.values()), "even B finite-tail onset")
    require(all(onset % 2 == 1 for onset in d_onsets.values()), "even D finite-tail onset")
    b_count = sum(residual_pair_count(2 * rank + 1, onset) for rank, onset in b_onsets.items())
    d_count = sum(residual_pair_count(rank - 1, onset) for rank, onset in d_onsets.items())
    require(b_count == 137, f"B18--B21 residual pair count={b_count}")
    require(d_count == 4732, f"D31--D70 residual pair count={d_count}")

    statement = statement_with_label(text, "prop:bd-finite-residual-certificate")
    require_compact(statement, r"18\leqb\leq21", "B finite residual rank scope missing")
    require_compact(statement, r"31\leqb\leq70", "D finite residual rank scope missing")
    proof_start = text.find(
        r"\begin{proof}",
        text.find(r"\label{prop:bd-finite-residual-certificate}"),
    )
    proof_end = text.find(r"\end{proof}", proof_start)
    proof = text[proof_start:proof_end]
    require_compact(
        proof,
        r"exactly$137$type-$B$and$4732$type-$D$residualpairs",
        "finite residual scope counts missing",
    )
    return b_count, d_count


def validate_classical_rank_coverage(
    middle_ranges: dict[str, tuple[int, int]] | None = None,
) -> tuple[int, int]:
    middle = middle_ranges or {"B": (22, 295), "C": (29, 295), "D": (71, 297)}
    require(middle["B"] == (22, 295), f"type-B middle interval={middle['B']}")
    require(middle["C"] == (29, 295), f"type-C middle interval={middle['C']}")
    require(middle["D"] == (71, 297), f"type-D middle interval={middle['D']}")

    # Symbolic infinite-rank closure.  Types B/C have C_G=b, hence their
    # high theorem begins at b=296.  In type D it begins at even b=296 and
    # odd b=299.  Thus every b>=298 is high; D296 is the sole overlap with
    # the stated middle interval, while D297 is middle only.
    require(classical_negative_endpoint("B", 296) == 296, "type-B high boundary failed")
    require(classical_negative_endpoint("C", 296) == 296, "type-C high boundary failed")
    require(classical_negative_endpoint("D", 296) == 296, "even type-D high boundary failed")
    require(classical_negative_endpoint("D", 297) == 295, "odd type-D middle boundary failed")
    require(classical_negative_endpoint("D", 298) >= 296, "D298 high coverage failed")
    require(classical_negative_endpoint("D", 299) >= 296, "D299 high coverage failed")

    lower_rank = {"B": 2, "C": 2, "D": 4}
    overlaps = 0
    checked = 0
    for family in ("B", "C", "D"):
        for rank in range(lower_rank[family], 2001):
            owners = set()
            if family == "B" and 2 <= rank <= 17:
                owners.add("low")
            if family == "C" and 2 <= rank <= 28:
                owners.add("low")
            if family == "D" and 4 <= rank <= 30:
                owners.add("low")
            if family == "B" and 18 <= rank <= 21:
                owners.add("finite")
            if family == "D" and 31 <= rank <= 70:
                owners.add("finite")
            middle_start, middle_end = middle[family]
            if middle_start <= rank <= middle_end:
                owners.add("middle")
            if classical_negative_endpoint(family, rank) >= 296:
                owners.add("high")
            require(owners, f"classical rank gap at {family}{rank}")
            if len(owners) > 1:
                overlaps += 1
                require(
                    (family, rank, owners) == ("D", 296, {"middle", "high"}),
                    f"unexpected rank overlap {family}{rank}: {owners}",
                )
            checked += 1
    require(overlaps == 1, f"classical rank overlap count={overlaps}")
    return checked, overlaps


def validate_classical_statements(text: str) -> None:
    middle = statement_with_label(text, "prop:bcd-middle-rank-full-tail")
    require_compact(middle, r"B_b\(22\leqb\leq295)", "middle B interval missing")
    require_compact(middle, r"C_b\(29\leqb\leq295)", "middle C interval missing")
    require_compact(middle, r"D_b\(71\leqb\leq297)", "middle D interval missing")
    high = statement_with_label(text, "prop:bcd-high-rank-full-tail")
    require_compact(high, r"C_G\geq296", "high-rank endpoint hypothesis missing")
    endpoints = statement_with_label(text, "lem:full-classical-negative-endpoints")
    require_compact(endpoints, r"C_{B_b}=b", "type-B negative endpoint missing")
    require_compact(endpoints, r"C_{C_b}=b", "type-C negative endpoint missing")
    require_compact(endpoints, r"b,&b\text{even}", "even type-D endpoint missing")
    require_compact(endpoints, r"b-2,&b\text{odd}", "odd type-D endpoint missing")


def validate_final_scope(statement: str) -> None:
    required = (
        r"anycompactconnectedLiegroupwithsimpleLiealgebra",
        r"everyinteger$r\geq0$",
        r"everytuple$f_1,\ldots,f_r\in\Qcone_G$",
        r"everychoiceofsigns$\varepsilon_i\in\{+1,-1\}$",
        r"containinganevennumberofminussigns",
        r"Ifthenumberofminussignsisodd,theintegraliszero",
        r"allofitsoriginalfinite-tupleandsign-patternquantifiers",
    )
    normalized = compact(statement)
    for fragment in required:
        require(fragment in normalized, f"final theorem scope fragment missing: {fragment}")


def expect_rejected(name: str, action: Callable[[], object]) -> None:
    try:
        action()
    except AuditFailure:
        return
    raise AuditFailure(f"adversarial mutation was accepted: {name}")


def mutation_self_tests(
    text: str,
    rows: dict[tuple[str, int], LowRow],
    paper_tails: dict[tuple[str, int], tuple[int, str]],
) -> int:
    mutations = 0

    missing_row = copy.deepcopy(rows)
    del missing_row[("B", 17)]
    expect_rejected("omit B17 low row", lambda: validate_low_rows(missing_row, paper_tails))
    mutations += 1

    changed_tail = copy.deepcopy(rows)
    original = changed_tail[("C", 28)]
    changed_tail[("C", 28)] = LowRow(
        original.family,
        original.rank,
        original.tail_onset + 2,
        original.moment_through,
        original.method,
    )
    expect_rejected("change C28 onset", lambda: validate_low_rows(changed_tail, paper_tails))
    mutations += 1

    expect_rejected(
        "remove B295 from middle ranks",
        lambda: validate_classical_rank_coverage({"B": (22, 294), "C": (29, 295), "D": (71, 297)}),
    )
    mutations += 1

    expect_rejected("truncate type-A residual ranks", lambda: validate_general_type_a(text, 27))
    mutations += 1

    final_statement = statement_with_label(text, "thm:full-adjoint-generated-q3")
    mutated_scope = final_statement.replace("an even number of minus", "an odd number of minus", 1)
    expect_rejected("reverse final sign quantifier", lambda: validate_final_scope(mutated_scope))
    mutations += 1
    return mutations


def main() -> int:
    text = PAPER.read_text(encoding="utf-8")
    ledger_text = LOW_ROW_LEDGER.read_text(encoding="utf-8")

    exact_statement = statement_with_label(text, "prop:bcd-remaining-exact-tails")
    directed_statement = statement_with_label(text, "prop:bcd-remaining-directed-tails")
    exact_tails = parse_exact_low_tail_rows(exact_statement)
    directed_tails = parse_directed_low_tail_rows(directed_statement)
    require(not (set(exact_tails) & set(directed_tails)), "a low row has two tail methods")
    paper_tails = exact_tails | directed_tails
    low_rows = parse_low_row_ledger(ledger_text)

    low_residual = validate_low_rows(low_rows, paper_tails)
    explicit_a_residual = validate_explicit_type_a(text)
    general_a_residual = validate_general_type_a(text)
    exceptional = parse_exceptional_rows(text)
    b_finite_residual, d_finite_residual = validate_finite_bd_rows(text)
    validate_classical_statements(text)
    classical_ranks, rank_overlaps = validate_classical_rank_coverage()

    final_statement = statement_with_label(text, "thm:full-adjoint-generated-q3")
    validate_final_scope(final_statement)
    require_compact(text, r"-\frac8{279}", "SO(3) all-PD obstruction absent")
    require_compact(
        text,
        r"Itdoesnotreplace$\Qcone_G$bytheconeofallrealcontinuouspositive-definitefunctions",
        "final scope boundary absent",
    )
    require_compact(
        text,
        r"Classificationleavestype$A$,types$B,C,D$,andthefiveexceptionaltypes",
        "final classification exhaustion sentence absent",
    )

    mutation_count = mutation_self_tests(text, low_rows, paper_tails)
    exceptional_residual = sum(row[2] for row in exceptional.values())
    print(
        "FULL_Q3_CLASSIFICATION_COVERAGE "
        f"low_rows={len(low_rows)} low_residual_pairs={low_residual} "
        f"explicit_A_residual_pairs={explicit_a_residual} "
        f"general_A_residual_pairs={general_a_residual} "
        f"exceptional_groups={len(exceptional)} exceptional_residual_pairs={exceptional_residual} "
        f"finite_B_residual_pairs={b_finite_residual} "
        f"finite_D_residual_pairs={d_finite_residual} "
        f"sampled_classical_ranks={classical_ranks} rank_overlaps={rank_overlaps} "
        "rank_gaps=0 degree_gaps=0 scope_gaps=0 "
        f"mutations_rejected={mutation_count}"
    )
    print("FULL_Q3_CLASSIFICATION_COVERAGE VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AuditFailure as error:
        raise SystemExit(f"FULL_Q3_CLASSIFICATION_COVERAGE FAILURE: {error}") from error
