#!/usr/bin/env python3
"""Fail-closed binding of the final publication sources.

The historical filename is retained for manifest stability.  The verifier now
binds the exact final sources directly, checks their publication and self-audit
boundaries, validates the live manifests, and proves that every arithmetic
source in the Part III manifest is byte-identical to the immutable execution
snapshot.  Historical prose-reversal constants remain below as provenance but
are no longer the publication acceptance condition.
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
    "replay_doc": ROOT / "REPLAY.md",
    "proof_spine": ROOT / "PUBLICATION_PROOF_SPINE.md",
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
    "extension": "9362671cfb5fbd4750a11a072b13f07a64073362e43e84da198ef73f12bd9439",
    "document_verifier": "c0b09862c12ea951f1a532a55848deabb6f222a6c4b9a77ba7b23be524b80997",
    "source_table": "a1556141b7032fc256da94a6fe72eaf80e2d1023e29265df361e7a1acc0eb39e",
    "source_audit": "0302d465759c5aa02bdd91d83271edd6645b27131062adf2a6c23a55a527f68b",
    "reference_manifest": "925ce91059884fe78f91a750513ecd085b7cfcd15564b83b00c2189486bac700",
    "full_manifest": "62a7592804a090191d5103120908ceeac10005a7649e7fba898481a4a244eefe",
    "paper": "dba5091d17cd85aeeee73deaafda86df1fe05e395625489a391201898910e325",
    "readme": "b9b48665b1513f426fa549e2d0dd7e6a00e2a4b1d923f7f66dbc53c4fe0c6597",
    "replay_doc": "65572ec18410f7bcf31417fdf276765d9fd50271bb7db2ca7b1f84bf1409881d",
    "proof_spine": "810a546718d3506008f79baf319c7e98755d8fb9f9625472cb1682b72381eaff",
    "replay_manifest": "83376b78f1dd8effce4bea92ede23d04bc0522b2d9c6feddae6819f637a5c2b3",
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

The audited computation package is rooted at immutable publication import
commit `91363c3cebb37a44349f613ab0a5aa6dcd412af3` in
`https://github.com/skypher/ginibre-q3-adjoint`; tag `v1.0.0` identifies the
publication release containing this binding. Parts I--II (`paper.tex`,
`paper_full.tex`), Part III
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
\url{https://github.com/skypher/ginibre-q3-adjoint/tree/91363c3cebb37a44349f613ab0a5aa6dcd412af3/ginibre_q3}."""

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
    require("v1.0.0" in text, "publication release tag is absent")
    require(
        "91363c3cebb37a44349f613ab0a5aa6dcd412af3" in text,
        "immutable publication import commit is absent",
    )
    require(
        "intermediate cone of all real continuous" in text
        and "central positive-definite functions" in text,
        "central positive-definite cone scope boundary is absent",
    )
    require(
        "three inseparable manuscript files" in text
        and "formal detailed" in text
        and "Part~III is not offered as" in text,
        "formal Parts I--III and supplement boundary is absent",
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
    paper = read("paper").decode("utf-8")
    readme = read("readme").decode("utf-8")
    replay_doc = read("replay_doc").decode("utf-8")
    proof_spine = read("proof_spine").decode("utf-8")
    require("pdfauthor={Leslie P. Polzer}" in extension, "Part III PDF metadata is absent")
    require("pdfauthor={Leslie P. Polzer}" in paper, "companion PDF metadata is absent")
    require(
        "author-controlled self-audits" in extension
        and "not independent peer review" in extension,
        "Part III self-audit boundary is absent",
    )
    require(
        "author-controlled historical self-audits" in readme,
        "README self-audit boundary is absent",
    )
    require(
        "Reproduction tiers and resource envelope" in replay_doc
        and "44 min 55 s" in replay_doc
        and "eight-hour guards" in replay_doc,
        "resource-tier disclosure is incomplete",
    )
    require(
        "# Publication proof spine" in proof_spine
        and "thm:full-adjoint-generated-q3" in proof_spine
        and "author-controlled" in proof_spine,
        "publication proof spine is incomplete",
    )

    full_records_live = parse_manifest(
        FILES["full_manifest"].read_text(encoding="utf-8"), FILES["full_manifest"]
    )
    snapshot_manifest = (
        ROOT
        / "certificates/full_q3/distributed0001/execution_source_snapshot/"
          "ginibre_q3/certificates/full_q3/full_q3_source_manifest.sha256"
    )
    require(snapshot_manifest.is_file(), "immutable execution manifest is absent")
    snapshot_records = parse_manifest(
        snapshot_manifest.read_text(encoding="utf-8"), snapshot_manifest
    )
    arithmetic_paths = [
        path
        for path in full_records_live
        if path == "../../Makefile"
        or path.startswith("../../character_ring_iter/")
        or path.startswith("../../run_full")
    ]
    require(arithmetic_paths, "no arithmetic sources found in live manifest")
    for path in arithmetic_paths:
        require(path in snapshot_records, f"arithmetic source absent from snapshot: {path}")
        require(
            full_records_live[path] == snapshot_records[path],
            f"arithmetic source changed after execution snapshot: {path}",
        )

    reference_records = verify_manifest(FILES["reference_manifest"])
    full_records = verify_manifest(FILES["full_manifest"])
    print(
        "FULL_Q3_FINAL_SOURCE_BINDING "
        f"extension={FINAL['extension']} arithmetic_sources={len(arithmetic_paths)} "
        f"reference_records={reference_records} full_records={full_records}"
    )
    print("FULL_Q3_FINAL_SOURCE_BINDING VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RebindFailure as error:
        raise SystemExit(f"FULL_Q3_RELATED_WORK_REBIND FAILURE: {error}")
