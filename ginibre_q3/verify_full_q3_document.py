#!/usr/bin/env python3
"""Fail-closed structural audit for the Part III full-Q3 manuscript."""

from __future__ import annotations

import re
import csv
from fractions import Fraction
from pathlib import Path


SOURCE = Path(__file__).with_name("full_q3_extension.tex")
RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FULL_Q3_DOCUMENT FAILURE: {message}")


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")

    result_pattern = re.compile(
        r"\\begin\{(" + "|".join(RESULT_ENVS) + r")\}(.*?)\\end\{\1\}",
        re.DOTALL,
    )
    results = list(result_pattern.finditer(text))
    proofs = list(re.finditer(r"\\begin\{proof\}.*?\\end\{proof\}", text, re.DOTALL))
    require(len(results) == 51, f"expected 51 numbered results, found {len(results)}")
    require(len(proofs) == 51, f"expected 51 proofs, found {len(proofs)}")

    for index, result in enumerate(results, 1):
        following = text[result.end() :]
        require(
            re.match(r"\s*\\begin\{proof\}", following) is not None,
            f"numbered result {index} is not followed immediately by a proof",
        )

    # A theorem may of course be stated before a later theorem, but its proof
    # must not depend on that later theorem.  Such a forward edge can conceal a
    # simultaneous-computation packaging cycle: theorem A cites certificate B,
    # while B cites A to interpret its other output.  Associate every label in
    # a result/proof interval with that result and reject every proof reference
    # whose owner occurs later in the manuscript.  References to equations in
    # the same proof and to unnumbered background labels are harmless.
    result_records: list[tuple[re.Match[str], re.Match[str]]] = []
    for index, result in enumerate(results, 1):
        proof = re.match(
            r"\s*(\\begin\{proof\}.*?\\end\{proof\})",
            text[result.end() :],
            re.DOTALL,
        )
        require(proof is not None, f"cannot parse proof following result {index}")
        result_records.append((result, proof))

    label_owner: dict[str, int] = {}
    for index, (result, _proof) in enumerate(result_records):
        interval_end = (
            result_records[index + 1][0].start()
            if index + 1 < len(result_records)
            else len(text)
        )
        for label in re.findall(r"\\label\{([^}]+)\}", text[result.start() : interval_end]):
            require(label not in label_owner, f"duplicate owned label {label}")
            label_owner[label] = index

    forward_dependencies: list[tuple[int, str, int]] = []
    for index, (result, proof_match) in enumerate(result_records):
        proof_text = proof_match.group(1)
        for target in re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", proof_text):
            owner = label_owner.get(target)
            if owner is not None and owner > index:
                forward_dependencies.append((index + 1, target, owner + 1))
    require(
        not forward_dependencies,
        f"forward proof dependencies: {forward_dependencies}",
    )

    labels = re.findall(r"\\label\{([^}]+)\}", text)
    require(len(labels) == len(set(labels)), "duplicate LaTeX label")
    references = re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", text)
    missing_references = sorted(set(references) - set(labels))
    require(not missing_references, f"undefined references: {missing_references}")

    citation_keys: list[str] = []
    for group in re.findall(r"\\cite(?:\[[^]]*\])?\{([^}]+)\}", text):
        citation_keys.extend(key.strip() for key in group.split(","))
    bibliography_keys = re.findall(r"\\bibitem\{([^}]+)\}", text)
    require(
        set(citation_keys) == set(bibliography_keys),
        "citation keys and bibliography items do not agree exactly",
    )
    require(len(bibliography_keys) == len(set(bibliography_keys)), "duplicate bibliography key")

    source_table = SOURCE.with_name("references") / "active_paper_sources.tsv"
    with source_table.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    mapped_keys = {row["bibkey"] for row in source_rows}
    external_keys = set(citation_keys) - {"AdjointTwoMinus"}
    require(external_keys <= mapped_keys, "an external citation lacks a primary-source mapping")
    for row in source_rows:
        if row["bibkey"] not in external_keys:
            continue
        require(bool(row["verified_locus"].strip()), f"empty source locus for {row['bibkey']}")
        require(
            (source_table.parent / row["path"]).is_file(),
            f"missing source artifact for {row['bibkey']}: {row['path']}",
        )

    required_labels = {
        "lem:adjoint-cone-positive-definite",
        "thm:full-cone-reduction",
        "prop:all-pd-counterexample",
        "lem:finite-rank-fredholm-mgf",
        "prop:bcd-middle-rank-full-tail",
        "prop:bcd-bounded-littlewood-residual-certificate",
        "thm:bcd-full-hierarchy",
        "thm:full-adjoint-generated-q3",
    }
    require(required_labels <= set(labels), "a required proof-spine label is absent")
    require("-\\frac8{279}" in text, "the exact SO(3) obstruction is absent")
    require(
        "the fixed $+1$ eigenvalue cancels that shift exactly" in text,
        "the type-B deterministic/random trace cancellation is absent",
    )
    require("no symmetry of the trace law is being assumed" in text, "two-tail justification is absent")
    require("all real continuous positive-definite functions" in text, "scope boundary is absent")
    require(
        {"Sylvester1980", "Abdesselam2022", "ChevyrevGarban2025"} <= set(citation_keys),
        "the adjacent nonabelian literature is not cited",
    )
    require(
        "We therefore make no claim that this is the first nonabelian" in text,
        "Sylvester priority boundary is absent",
    )
    require(
        "does not settle the invariant-observable GKS2" in text,
        "finite-O(N) GKS2 scope boundary is absent",
    )
    require(
        "It does not by itself prove correlation inequalities for quantum" in text,
        "quantum-application limitation is absent",
    )
    require(
        "the occurrence of an adjoint representation in a gauge model neither" in text,
        "gauge-theory non-implication is absent",
    )
    require(
        "reflection-positive transfer-matrix or" in text,
        "quantum reconstruction prerequisites are absent",
    )
    require("Leslie P. Polzer" in text, "publication author identity is absent")
    require("\\date{July 19, 2026}" in text, "publication date is stale")
    require("Code and data availability" in text, "availability section is absent")
    require(
        "github.com/skypher/ginibre-q3-adjoint" in text,
        "official publication repository is absent",
    )
    require(
        "intermediate cone of all real continuous" in text
        and "central positive-definite functions" in text,
        "central positive-definite cone scope boundary is absent",
    )
    require(
        "two inseparable mathematical components" in text
        and "Part~III is not offered as" in text,
        "formal Parts I--III submission boundary is absent",
    )
    require("tab:certificate-map" in text, "trusted-computation map is absent")
    require(
        all(
            identifier in text
            for identifier in (
                "distributed0001",
                "bcdanalytic0003",
                "bcdlowtail0001",
                "bcd0010",
                "bcd0002",
            )
        ),
        "a load-bearing certificate family is absent from the map",
    )
    require(
        "exact integer/rational or directed-rounding" in text
        and "outward-rounded 384-bit MPFR" in text,
        "certificate arithmetic boundaries are absent",
    )
    require(
        "the Ginibre monotonicity consequence established there is abelian" in text,
        "Chevyrev--Garban gauge/Ginibre scope distinction is absent",
    )

    c1, c2, c3 = (Fraction(7, 31), Fraction(19, 31), Fraction(19, 31))
    fourth_moment = Fraction(2, 15) * (c1 * c1 + c2 * c2 + c3 * c3)
    fourth_moment -= Fraction(1, 15) * (c1 * c2 + c1 * c3 + c2 * c3)
    signed_value = 2 * (
        fourth_moment + Fraction(1, 9) * (c1 * c1 - c2 * c2 - c3 * c3)
    )
    require(fourth_moment == Fraction(61, 961), "SO(3) fourth moment mismatch")
    require(signed_value == Fraction(-8, 279), "SO(3) signed obstruction mismatch")

    print(
        "FULL_Q3_DOCUMENT results=51 proofs=51 "
        f"labels={len(labels)} references={len(references)} "
        f"citations={len(set(citation_keys))} source_mappings={len(external_keys)} "
        f"forward_dependencies={len(forward_dependencies)} "
        f"so3_obstruction={signed_value}"
    )
    print("FULL_Q3_DOCUMENT VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
