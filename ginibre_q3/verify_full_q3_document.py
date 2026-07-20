#!/usr/bin/env python3
"""Fail-closed structural audit for the Part III full-Q3 manuscript."""

from __future__ import annotations

import re
import csv
from fractions import Fraction
from pathlib import Path


SOURCE = Path(__file__).with_name("full_q3_extension.tex")
COMPACT_SOURCE = SOURCE.with_name("paper.tex")
READER_PART_III = SOURCE.with_name("full_q3_main.tex")
SUPPLEMENT_WRAPPER = SOURCE.with_name("paper_full.tex")
UNIFIED_WRAPPER = SOURCE.with_name("submission.tex")
ENVIRONMENT = SOURCE.with_name("ENVIRONMENT.md")
ACTIVE_PROOF = SOURCE.with_name("ACTIVE_PROOF_SUPPLEMENT.md")
FINAL_REPLAY = SOURCE.with_name("run_final_publication_replay.sh")
FINAL_ARCHIVE = SOURCE.with_name("build_final_release_archive.py")
ARTIFACT_VERIFIER = SOURCE.with_name("verify_publication_artifacts.py")
ARTIFACT_MANIFEST = SOURCE.with_name("publication_artifacts.sha256")
MODULAR_CHECKER = SOURCE.parent / "character_ring_iter" / (
    "verify_full_q3_bcd_modular_moment_checker.cpp"
)
MODULAR_TRANSCRIPT = SOURCE.parent / "certificates" / "full_q3" / (
    "fullq3bcdmodularfinal0001_current_source.log"
)
RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FULL_Q3_DOCUMENT FAILURE: {message}")


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")
    compact_text = COMPACT_SOURCE.read_text(encoding="utf-8")
    reader_part_iii = READER_PART_III.read_text(encoding="utf-8")
    supplement_wrapper = SUPPLEMENT_WRAPPER.read_text(encoding="utf-8")
    unified_wrapper = UNIFIED_WRAPPER.read_text(encoding="utf-8")
    require(ENVIRONMENT.is_file(), "validated environment record is absent")
    require(ACTIVE_PROOF.is_file(), "active proof supplement is absent")
    require(FINAL_REPLAY.is_file(), "final-source replay driver is absent")
    require(FINAL_ARCHIVE.is_file(), "final release archive builder is absent")
    require(ARTIFACT_VERIFIER.is_file(), "publication artifact verifier is absent")
    require(ARTIFACT_MANIFEST.is_file(), "publication artifact manifest is absent")
    require(MODULAR_CHECKER.is_file(), "compressed modular checker is absent")
    require(MODULAR_TRANSCRIPT.is_file(), "compressed modular transcript is absent")
    environment = ENVIRONMENT.read_text(encoding="utf-8")
    require(
        "No theorem imports a result" in compact_text
        and "optional derivation archive" in compact_text,
        "formal Parts I--II do not exclude the optional derivation archive",
    )
    require(
        r"\label{prop:bc-active-correction-prefix-contract}" in compact_text
        and compact_text.count(r"0\le r\le27") >= 2
        and "exactly the correction" in compact_text
        and "prefix consumed by the main theorem" in compact_text,
        "compact Parts I--II do not state the active B/C contract and consumed offset",
    )
    require(
        r"\cref{prop:post29-bc-local-half-bridge}" in compact_text
        and "post_m29_bc_interval_bridge_frontier_gmp.cpp" in compact_text
        and r"D_G(2m+1)\ge\mathcal L_m" in compact_text,
        "compact Parts I--II omit the proved half-stable bridge",
    )
    require(
        r"\label{lem:bc-active-finite-adjoint-moments}" in compact_text
        and r"\label{lem:bc-active-bounded-littlewood-determinants}" in compact_text
        and r"\label{lem:bc-active-two-power-functional}" in compact_text
        and r"m_j^{C_b}=[m_{(2^j)}]\mathcal C_b" in compact_text,
        "compact Parts I--II omit the exact finite B/C moment specification",
    )
    require(
        "$B_1=C_1=A_1$" in compact_text
        and "$D_3=A_3$" in compact_text
        and "$D_2$ is not simple" in compact_text,
        "compact Parts I--II omit the low-rank classification conventions",
    )
    require(
        r"\section{Conclusion and limitations}" in compact_text,
        "compact Parts I--II omit the conclusion and scope limitations",
    )
    require(
        "Optional derivation archive" in supplement_wrapper
        and r"\def\GinibreFullProof{1}" in supplement_wrapper,
        "optional derivation archive wrapper is absent or inactive",
    )
    require(
        "submission.pdf" in text
        and "load-bearing theorem dependency" in text
        and "GinibreDevelopmentArchive" in text,
        "Part III supplement does not state the reader/supplement architecture",
    )
    require(
        "paper.pdf" in unified_wrapper
        and "full_q3_main.pdf" in unified_wrapper
        and "full_q3_extension.pdf" in unified_wrapper
        and r"paper\_full.pdf" in unified_wrapper
        and "every load-bearing proof component" in unified_wrapper,
        "unified wrapper omits a load-bearing manuscript component",
    )
    require(
        r"\label{prop:reduction}" in reader_part_iii
        and r"\label{lem:bl}" in reader_part_iii
        and r"\label{prop:contract}" in reader_part_iii
        and r"\label{thm:main}" in reader_part_iii
        and "17{,}862" in reader_part_iii
        and "included as the final component" in reader_part_iii,
        "compact Part III omits a load-bearing reader interface",
    )
    require(
        "certificates/full_q3/full_q3_source_manifest.sha256" in text,
        "Part III does not name the live final-source manifest path",
    )
    require(
        "publication_artifacts.sha256" in text
        and "publication-artifact-audit" in FINAL_REPLAY.read_text(encoding="utf-8"),
        "final publication flow does not rebuild and bind the reader PDFs",
    )
    proof_spine_text_coverage = ACTIVE_PROOF.with_name("PUBLICATION_PROOF_SPINE.md").read_text(
        encoding="utf-8"
    )
    modular_transcript = MODULAR_TRANSCRIPT.read_text(encoding="utf-8")
    require(
        "12,993-case low-row subledger" in proof_spine_text_coverage
        and "other 4,869 cases" in proof_spine_text_coverage
        and "no single secondary checker covers the complete box" in proof_spine_text_coverage
        and "prefix_pairs=7484 high_pairs=5509" in modular_transcript,
        "secondary finite-box coverage is overstated or inconsistent",
    )
    require(
        "run_final_publication_replay.sh" in text
        and "build_final_release_archive.py" in text
        and "verify_full_q3_bcd_modular_moment_checker.cpp" in text
        and "$538$" in text
        and "$5{,}509$" in text
        and "$628$" in text
        and "$607$" in text,
        "Part III omits the final-source replay/archive or compressed cross-check",
    )
    require(
        r"\path{ENVIRONMENT.md}" in text
        and "Ubuntu 24.04.4 LTS" in environment
        and "GNU g++ 13.3.0" in environment
        and "MPFR 4.2.1" in environment,
        "validated replay environment is incomplete or not cited",
    )
    require(
        "maximal cone to which" not in text
        and "largest function class that" not in compact_text,
        "unsupported maximality wording remains in the manuscript",
    )

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
    require("v1.0.0" in text, "historical snapshot tag is absent")
    require(
        "91363c3cebb37a44349f613ab0a5aa6dcd412af3" in text,
        "historical execution commit is absent",
    )
    require(
        "intermediate cone of all real continuous" in text
        and "central positive-definite functions" in text,
        "central positive-definite cone scope boundary is absent",
    )
    require(
        "submission.pdf" in text
        and "optional derivation archive" in text
        and "Part~III is not offered as" in text,
        "unified Parts I--III boundary is absent",
    )
    require("tab:certificate-map" in text, "trusted-computation map is absent")
    require(
        all(
            identifier in text
            for identifier in (
                "distributed0001",
                "bcdanalytic0003",
                "bcdlowtail0001",
                "bcdboundedfinal0001",
                "bcdmodularfinal0001",
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
    require(
        "author-controlled self-audits" in text
        and "not independent peer review" in text,
        "author self-audit boundary is absent",
    )
    require(
        "three reproducibility tiers" in text
        and "five-minute ceiling" in text
        and "proved analytic reduction" in text,
        "reproduction-tier resource disclosure is absent",
    )
    require(
        "C_\\ell=\\mathrm{Sp}(2\\ell)" in text,
        "compact symplectic-group notation is ambiguous",
    )
    proof_spine = SOURCE.with_name("PUBLICATION_PROOF_SPINE.md")
    require(proof_spine.is_file(), "publication proof-spine guide is absent")
    proof_spine_text = proof_spine.read_text(encoding="utf-8")
    require(
        all(
            marker in proof_spine_text
            for marker in (
                "thm:main",
                "thm:full-cone-reduction",
                "thm:bcd-full-hierarchy",
                "thm:full-adjoint-generated-q3",
                "author-controlled",
            )
        ),
        "publication proof-spine guide is incomplete",
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
        "FULL_Q3_DOCUMENT results=51 proofs=51 proof_spine=PASS "
        f"labels={len(labels)} references={len(references)} "
        f"citations={len(set(citation_keys))} source_mappings={len(external_keys)} "
        f"forward_dependencies={len(forward_dependencies)} "
        f"so3_obstruction={signed_value}"
    )
    print("FULL_Q3_DOCUMENT VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
