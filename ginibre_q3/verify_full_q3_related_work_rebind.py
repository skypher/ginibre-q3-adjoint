#!/usr/bin/env python3
"""Fail-closed audit of the publication-only metadata and scope rebind.

The distributed computation and the first paper rebind were accepted before
this edit.  This verifier proves that the subsequent changes are fixed
authorship metadata, related-work and application-scope remarks, bibliography
and availability additions, their source audit, and fail-closed document
checks.  Reversing exactly those changes must recover the accepted pre-edit
hashes.  It also verifies every final file in the reference, replay, and Part
III manifests.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent

FILES = {
    "extension": ROOT / "full_q3_extension.tex",
    "document_verifier": ROOT / "verify_full_q3_document.py",
    "source_table": ROOT / "references" / "active_paper_sources.tsv",
    "source_audit": ROOT / "references" / "GINIBRE_RELATED_WORK_SOURCE_AUDIT.md",
    "reference_manifest": ROOT / "references" / "references_manifest.sha256",
    "full_manifest": ROOT / "certificates" / "full_q3" / "full_q3_source_manifest.sha256",
    "paper": ROOT / "paper.tex",
    "readme": ROOT / "README.md",
    "replay_manifest": ROOT / "replay_sources.sha256",
}

BASELINE = {
    "extension": "48d269b8d0ab67d15e49f99256fd0e36dd559691bf2039d2c373ca99fc483c02",
    "document_verifier": "a92f75a0cfdd5ea7f7bdea8eafc11a3481ff4227d2f29b7de6f156bcc1cdc26e",
    "source_table": "323d268e12bb285d61038e28e80f71a3920f2af36ff929e8019846c35908bba1",
    "reference_manifest": "d03b7c4b2450239a26156fcbfa51daad54b022c08f7cd08550327b0af1667713",
    "full_manifest": "1cbd33b571e990c8b8bce47a1bd549129dea6cc89d965ec5951ec6f8f265fb59",
    "paper": "d2fb8944ec081f47a1bfcf558bb65adb01936b0ac1d0dbbb1f561d63d7197010",
    "readme": "b78ad3a817f7a966003fa8147808d2ffea8bfa67c83c1a27b86f727695470464",
    "replay_manifest": "bfe68eaf0ee4a5fe8263902b1145804a630381a54bfd4a9b4bb3135517196874",
}

FINAL = {
    "extension": "b7b872f7d382f895bbf1ce5e16f123642e36d94745194e5b08e20c14fe48f5b8",
    "document_verifier": "9d07c924434f41a5395331171c7a593d92052a5064ab3b4cb69fcbc40d9ca4ba",
    "source_table": "b5cf73cae7dfa1499962d03d8c65128b09ed9b4c87950b4f64117923674bb03e",
    "source_audit": "0302d465759c5aa02bdd91d83271edd6645b27131062adf2a6c23a55a527f68b",
    "reference_manifest": "c226c70b6305409665cb4a4f740cbb5cdb86a5879bb4aba995bb2fcd80103d52",
    "full_manifest": "433b6073e40ea393f86e9c3df2a7019852f29cbf3208fd171992f1479143909a",
    "paper": "8b58f9015dc4738c0e48e3d340f512a74b54a6832f83a3ae378b910d75e6060f",
    "readme": "8de77a263088bb9e46ce93340777ea22f37586f08e4d01774307df339f926d6a",
    "replay_manifest": "04442536fdfaf8505d56432a2b28fc957063eb70be2edb8766cb0decbc835923",
}

FINAL_EXTENSION_METADATA = """\\author{Leslie P. Polzer}
\\address{}
\\email{polzer@fastmail.com}"""
BASELINE_EXTENSION_METADATA = r"\author{Anonymous}"
FINAL_PAPER_METADATA = """\\author{Leslie P. Polzer}
\\address{}
\\email{polzer@fastmail.com}"""
BASELINE_PAPER_METADATA = """\\author{Anonymous}
\\address{}
\\email{}"""

RELATED_REMARK = r"""
\begin{remark}[Relation to earlier nonabelian correlation inequalities]
\label{rem:related-nonabelian-inequalities}
Sylvester~\cite[Secs.~1, 3--4]{Sylvester1980} had already disproved the
general Ginibre inequality for higher-dimensional $O(N)$ spin systems by
constructing graph-cycle counterexamples.  His variables are spins on a
product of spheres and his signed factors are indexed by graph edges.  The
present Proposition~\ref{prop:all-pd-counterexample} supplies a different,
exact witness directly in the compact-group formulation proposed at the end
of Ginibre's 1970 paper: four explicit real positive-definite functions on
$\mathrm{SO}(3)$, integrated over two normalized Haar variables, give
$-8/279$.  We therefore make no claim that this is the first nonabelian
failure of a Ginibre inequality.

Abdesselam~\cite[Sec.~1 and Thms.~2.1--2.3]{Abdesselam2022} subsequently
formulated general and padded Ginibre inequalities for invariant observables
of $O(N)$ and $\mathbb{CP}^{N-1}$ models.  He proves large-power asymptotics
in terms of Kirchhoff polynomials and proves all padded Ginibre inequalities
for inverse half-integer powers of stable determinantal polynomials.  Those
results concern an asymptotic regime and a determinantal-polynomial
substitute; they do not establish the finite-rank compact-group Haar
hierarchy proved here.  Conversely, the present theorem is restricted to the
adjoint-generated cone and does not settle the invariant-observable GKS2
problem for finite $O(N)$ spin systems, which Abdesselam records as open for
$N\geq3$.
\end{remark}

\begin{remark}[Scope of potential quantum applications]
\label{rem:quantum-application-scope}
Theorem~\ref{thm:full-adjoint-generated-q3} is an algebraic Haar-positivity
theorem.  Through Ginibre's general mechanism it can serve as an input for
group-valued lattice models whose interactions genuinely lie in
$\Qcone_G$.  It does not by itself prove correlation inequalities for quantum
Heisenberg or general $O(N)$ vector-spin systems, nonlinear sigma models,
lattice gauge theories, or continuum quantum field theories.  In particular,
the occurrence of an adjoint representation in a gauge model neither places
its plaquette interactions and observables in $\Qcone_G$ nor removes the
shared-link constraints.  This distinction is visible in the recent work of
Chevyrev and Garban~\cite[Thm.~1.5, Cor.~1.6, and Rem.~1.7]{ChevyrevGarban2025}:
their Villain-action limit applies to nonabelian lattice gauge theories, while
the Ginibre monotonicity consequence established there is abelian.

Such applications require additional model-specific results: Gibbs-measure
promotion on the relevant observable cone, a thermodynamic-limit argument,
and, for a quantum reconstruction, a reflection-positive transfer-matrix or
equivalent construction.  A direct route toward the conventional finite-$N$
$O(N)$ problem would instead require the two-minus or padded inequality for
defining-representation matrix coefficients on products of spheres.  That is
a separate theorem and is not asserted here.
\end{remark}
"""

BIBLIOGRAPHY_ADDITIONS = r"""
\bibitem{Sylvester1980}
G.~S. Sylvester,
The Ginibre inequality,
\emph{Comm. Math. Phys.} \textbf{73} (1980), no.~2, 105--114.

\bibitem{Abdesselam2022}
A.~Abdesselam,
Non-Abelian correlation inequalities and stable determinantal polynomials,
arXiv:2207.07603, 2022.

\bibitem{ChevyrevGarban2025}
I.~Chevyrev and C.~Garban,
Villain action in lattice gauge theory,
\emph{J. Stat. Phys.} \textbf{192} (2025), Paper~38.
"""

AVAILABILITY_START = "\n\\section*{Code and data availability}\n"
AVAILABILITY_END = "\\begin{thebibliography}{99}"

SCOPE_CORRECTION = r"""The matrix coefficients in Proposition~\ref{prop:all-pd-counterexample} are
not central.  Consequently neither that counterexample nor the affirmative
theorem decides $Q_3$ for the intermediate cone of all real continuous
central positive-definite functions on a nonabelian compact group.
"""

FINAL_EXTENSION_DATE = r"\date{July 19, 2026}"
BASELINE_EXTENSION_DATE = r"\date{July 15, 2026}"
FINAL_PAPER_DATE = r"\date{July 19, 2026}"
BASELINE_PAPER_DATE = r"\date{July 14, 2026}"

README_ADDITIONS = r"""
The related-work boundary is explicit.  Sylvester (Comm. Math. Phys. 73
(1980), 105--114) already gave graph-cycle counterexamples for
higher-dimensional `O(N)` spins, so the `SO(3)` calculation is not claimed as
the first nonabelian failure of a Ginibre inequality.  Abdesselam
(arXiv:2207.07603) proves large-power asymptotics and padded inequalities for
stable-determinantal substitutes; those results neither imply nor are implied
by the fixed finite-rank adjoint-character Haar hierarchy proved here.  The
source-level comparison is recorded in
`references/GINIBRE_RELATED_WORK_SOURCE_AUDIT.md`.

The quantum-application boundary is equally explicit.  Part III is an
algebraic Haar-positivity theorem and can feed Ginibre's mechanism for
group-valued lattice models only when their interactions actually lie in the
adjoint-generated cone.  It does not establish correlation inequalities for
quantum Heisenberg or general `O(N)` vector spins, nonlinear sigma models,
lattice gauge theories, or continuum QFT.  Those uses need separate
model-specific Gibbs promotion, thermodynamic limits, and quantum
reconstruction machinery.

The audited computation package is archived in the official publication
repository, `https://github.com/skypher/ginibre-q3-adjoint`. Parts I--II
(`paper.tex`, `paper_full.tex`), Part III
(`full_q3_extension.tex`), all replay sources, accepted transcripts, and
SHA-256 manifests are archived together there. Parts I--II and Part III are
inseparable formal submission components.  Part III contains a
publication-facing code and data availability statement, exact build commands,
software requirements, and a result-to-certificate map.  It also states
explicitly that the intermediate cone of all real continuous central
positive-definite functions remains undecided.
"""

FINAL_README_PROOF_LIST = """- `CLASSICAL_M3_ONE_PROOF.md` supplies the numbered classical proof chain.
- `character_ring_iter/CLOSURE_STATUS.md` records the exceptional-family
  prefix and tail split."""
BASELINE_README_PROOF_LIST = """- `CLASSICAL_M3_ONE_PROOF.md` supplies the numbered classical proof chain.
- `POST_M29_UNCONDITIONAL_BRIEF.md` records the final post-`m=29` assembly;
  only its opening current-status section is authoritative, while its later
  route diary is historical.
- `character_ring_iter/CLOSURE_STATUS.md` records the exceptional-family
  prefix and tail split."""

FINAL_README_HISTORY = """Historical route diaries, superseded referee drafts, and unrelated research
notes are intentionally excluded from this publication repository."""
BASELINE_README_HISTORY = """`NONABELIAN.md`, the long post-`m=29` brief, old referee reports, and files
whose titles contain `TODO` are retained as research history.  They do not
override `paper.tex`, this status page, or `REPLAY.md`."""

FINAL_COMPANION_BIB = r"""\bibitem{AdjointTwoMinus}
L.~P. Polzer,
The two-minus adjoint-character case of Ginibre's $Q_3$ condition,
Parts I--II: classification, exact certificates, and expanded proof
companion, companion manuscript and source archive, 2026,
\url{https://github.com/skypher/ginibre-q3-adjoint}."""

BASELINE_COMPANION_BIB = r"""\bibitem{AdjointTwoMinus}
Anonymous,
The two-minus adjoint-character case of Ginibre's $Q_3$ condition,
Parts I--II: classification, exact certificates, and expanded proof
companion, companion manuscript, 2026."""

DOCUMENT_CHECKS = '''    require(
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
    require("\\\\date{July 19, 2026}" in text, "publication date is stale")
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
'''

SOURCE_ROWS = (
    "Sylvester1980\tGINIBRE_RELATED_WORK_SOURCE_AUDIT.md\tSections 1, 3--4 and appendix; graph-cycle formulation, higher-dimensional counterexamples, and exact graph-15 value; direct author-uploaded primary text access recorded in audit\n"
    "Abdesselam2022\tGINIBRE_RELATED_WORK_SOURCE_AUDIT.md\tSection 1 general/padded Ginibre definitions, Sylvester comparison, and finite-N GKS2 open status; Theorems 2.1--2.3 large-power asymptotics and stable-determinantal PGG theorem; arXiv PDF digest recorded in audit\n"
    "ChevyrevGarban2025\tGINIBRE_RELATED_WORK_SOURCE_AUDIT.md\tTheorem 1.5 nonabelian Villain-action limit; Corollary 1.6 and Remark 1.7 abelian Ginibre monotonicity scope; direct open-access publisher reading recorded in audit\n"
)


class RebindFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RebindFailure(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read(key: str) -> bytes:
    path = FILES[key]
    require(path.is_file(), f"missing {path.relative_to(ROOT)}")
    return path.read_bytes()


def remove_once(text: str, passage: str, description: str) -> str:
    require(text.count(passage) == 1, f"expected one {description}")
    return text.replace(passage, "", 1)


def remove_between_once(text: str, start: str, end: str, description: str) -> str:
    require(text.count(start) == 1, f"expected one start marker for {description}")
    require(text.count(end) == 1, f"expected one end marker for {description}")
    left, remainder = text.split(start, 1)
    _removed, right = remainder.split(end, 1)
    return left + "\n" + end + right


def replace_once(text: str, new: str, old: str, description: str) -> str:
    require(text.count(new) == 1, f"expected one final {description}")
    require(old not in text, f"baseline {description} remains in final file")
    return text.replace(new, old, 1)


def parse_manifest(text: str, path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for number, line in enumerate(text.splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"malformed manifest row {path}:{number}")
        value, recorded = match.groups()
        require(recorded not in records, f"duplicate manifest path {recorded}")
        records[recorded] = value
    require(records, f"empty manifest {path}")
    return records


def verify_manifest(path: Path) -> int:
    records = parse_manifest(path.read_text(encoding="utf-8"), path)
    root = ROOT.resolve()
    for recorded, expected in records.items():
        source = (path.parent / recorded).resolve()
        require(root in source.parents, f"manifest path escapes source root: {recorded}")
        require(source.is_file(), f"missing manifest input: {recorded}")
        require(digest(source.read_bytes()) == expected, f"manifest mismatch: {recorded}")
    return len(records)


def main() -> int:
    for key, expected in FINAL.items():
        require(digest(read(key)) == expected, f"unexpected final hash for {key}")

    extension = read("extension").decode("utf-8")
    extension = replace_once(
        extension,
        FINAL_EXTENSION_DATE,
        BASELINE_EXTENSION_DATE,
        "Part III publication date",
    )
    extension = remove_once(extension, SCOPE_CORRECTION, "central-cone scope correction")
    extension = remove_once(extension, RELATED_REMARK, "related-work remark")
    extension = remove_once(extension, BIBLIOGRAPHY_ADDITIONS, "bibliography additions")
    extension = remove_between_once(
        extension,
        AVAILABILITY_START,
        AVAILABILITY_END,
        "availability section",
    )
    extension = replace_once(
        extension,
        FINAL_COMPANION_BIB,
        BASELINE_COMPANION_BIB,
        "companion bibliography entry",
    )
    extension = replace_once(
        extension,
        FINAL_EXTENSION_METADATA,
        BASELINE_EXTENSION_METADATA,
        "Part III author metadata",
    )
    require(digest(extension.encode()) == BASELINE["extension"], "extension reversal mismatch")

    paper = read("paper").decode("utf-8")
    paper = replace_once(
        paper,
        FINAL_PAPER_DATE,
        BASELINE_PAPER_DATE,
        "Parts I--II publication date",
    )
    paper = replace_once(
        paper,
        FINAL_PAPER_METADATA,
        BASELINE_PAPER_METADATA,
        "Parts I--II author metadata",
    )
    require(digest(paper.encode()) == BASELINE["paper"], "companion reversal mismatch")

    readme = read("readme").decode("utf-8")
    readme = remove_once(readme, README_ADDITIONS, "README publication additions")
    readme = replace_once(
        readme,
        FINAL_README_PROOF_LIST,
        BASELINE_README_PROOF_LIST,
        "README proof-file list",
    )
    readme = replace_once(
        readme,
        FINAL_README_HISTORY,
        BASELINE_README_HISTORY,
        "README history boundary",
    )
    require(digest(readme.encode()) == BASELINE["readme"], "README reversal mismatch")

    verifier = read("document_verifier").decode("utf-8")
    verifier = remove_once(verifier, DOCUMENT_CHECKS, "document-verifier checks")
    require(
        digest(verifier.encode()) == BASELINE["document_verifier"],
        "document-verifier reversal mismatch",
    )

    table = read("source_table").decode("utf-8")
    table = remove_once(table, SOURCE_ROWS, "source-table rows")
    require(digest(table.encode()) == BASELINE["source_table"], "source-table reversal mismatch")

    references = read("reference_manifest").decode("utf-8")
    references = remove_once(
        references,
        FINAL["source_audit"] + "  GINIBRE_RELATED_WORK_SOURCE_AUDIT.md\n",
        "reference-manifest source-audit row",
    )
    references = replace_once(
        references,
        FINAL["source_table"] + "  active_paper_sources.tsv",
        BASELINE["source_table"] + "  active_paper_sources.tsv",
        "reference-manifest source-table row",
    )
    require(
        digest(references.encode()) == BASELINE["reference_manifest"],
        "reference-manifest reversal mismatch",
    )

    replay = read("replay_manifest").decode("utf-8")
    replay = replace_once(
        replay,
        FINAL["readme"] + "  README.md",
        BASELINE["readme"] + "  README.md",
        "replay-manifest README row",
    )
    replay = replace_once(
        replay,
        FINAL["paper"] + "  paper.tex",
        BASELINE["paper"] + "  paper.tex",
        "replay-manifest paper row",
    )
    require(
        digest(replay.encode()) == BASELINE["replay_manifest"],
        "replay-manifest reversal mismatch",
    )

    full = read("full_manifest").decode("utf-8")
    for key, recorded in (
        ("extension", "../../full_q3_extension.tex"),
        ("document_verifier", "../../verify_full_q3_document.py"),
        ("paper", "../../paper.tex"),
        ("replay_manifest", "../../replay_sources.sha256"),
        ("reference_manifest", "../../references/references_manifest.sha256"),
    ):
        full = replace_once(
            full,
            FINAL[key] + "  " + recorded,
            BASELINE[key] + "  " + recorded,
            f"full-manifest {recorded} row",
        )
    require(digest(full.encode()) == BASELINE["full_manifest"], "full-manifest reversal mismatch")

    reference_records = verify_manifest(FILES["reference_manifest"])
    full_records = verify_manifest(FILES["full_manifest"])
    print(
        "FULL_Q3_RELATED_WORK_REBIND "
        f"extension={FINAL['extension']} citations_added=3 "
        f"reference_records={reference_records} full_records={full_records}"
    )
    print("FULL_Q3_RELATED_WORK_REBIND VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RebindFailure as error:
        raise SystemExit(f"FULL_Q3_RELATED_WORK_REBIND FAILURE: {error}")
