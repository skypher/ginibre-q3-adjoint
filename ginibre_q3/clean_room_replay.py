#!/usr/bin/env python3
"""Isolated one-command replay for the correction-aware Ginibre Q3 manuscript.

Python is used only for orchestration, hashing, log parsing, and lightweight
exact audits.  Load-bearing sweeps run in C++ with GMP/OpenMP or directed MPFR.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


HASH_LINE = re.compile(r"^([0-9a-f]{64})\s+[*]?(.+?)\s*$")
HASH_SPLIT = re.compile(r"\\hashsplit\{([0-9a-f]+)\}\{([0-9a-f]+)\}")
PRINTED_LOG = re.compile(
    r"^\s+(?:(?:B/C|D)\s+)?((?:ginibre_q3/|certificates/)\S+\.log)\s*$"
)
CLASSICAL_DELTA_CLAIM = re.compile(r"([BCD])_(\d+).*?Delta_(\d+)\s*=\s*(-?\d+)")
CLASSICAL_MOMENT_CLAIM = re.compile(
    r"D_(\d+)\s+m_(\d+).*?moment_\d+\s*=\s*(-?\d+)"
)
CLASSICAL_DIRECT_MOMENT_CLAIM = re.compile(
    r"D_(\d+).*?m_(\d+)\s*=\s*(-?\d+)"
)
LATEX_FAILURE = re.compile(
    r"LaTeX Warning: (?:There were undefined references|Citation .+ undefined)|"
    r"Package rerunfilecheck Warning: File .+ has changed"
)
LATEX_PAGES = re.compile(r"Output written on .+? \((\d+) pages?,")


class ReplayFailure(RuntimeError):
    pass


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def resolve_manifest_path(root: Path, manifest: Path, text: str) -> Path:
    root = root.resolve()
    manifest = manifest.resolve()
    if Path(text).is_absolute():
        raise ReplayFailure(f"absolute path in {manifest}: {text}")
    candidates = []
    if text.startswith("ginibre_q3/"):
        candidates.append(root.parent / text)
    else:
        candidates.extend((manifest.parent / text, root / text))
    existing = [candidate.resolve() for candidate in candidates if candidate.exists()]
    if not existing:
        raise ReplayFailure(f"missing artifact named by {manifest}: {text}")
    outside = [path for path in existing if path != root and root not in path.parents]
    if outside:
        raise ReplayFailure(
            f"manifest path escapes the isolated tree in {manifest}: {text} -> {outside[0]}"
        )
    return existing[0]


def verify_manifests(root: Path) -> tuple[int, int]:
    root = root.resolve()
    symlinks = sorted(path for path in root.rglob("*") if path.is_symlink())
    if symlinks:
        raise ReplayFailure(
            "symlinks are forbidden in the isolated replay tree: "
            + ", ".join(str(path.relative_to(root)) for path in symlinks)
        )
    manifests = sorted(root.rglob("*.sha256"))
    checked = 0
    for manifest in manifests:
        entries = 0
        seen: set[Path] = set()
        for line_number, line in enumerate(manifest.read_text(errors="strict").splitlines(), 1):
            match = HASH_LINE.match(line)
            if not match:
                continue
            expected, relative = match.groups()
            artifact = resolve_manifest_path(root, manifest, relative)
            if artifact in seen:
                raise ReplayFailure(
                    f"duplicate artifact in {manifest}:{line_number}: {relative}"
                )
            seen.add(artifact)
            actual = digest(artifact)
            if actual != expected:
                raise ReplayFailure(
                    f"hash mismatch in {manifest}:{line_number}: {relative}\n"
                    f"expected {expected}\nactual   {actual}"
                )
            entries += 1
            checked += 1
        if entries == 0:
            raise ReplayFailure(f"manifest has no SHA-256 entries: {manifest}")
    return len(manifests), checked


def audit_certificate_identity_ledger(root: Path) -> tuple[int, int, int]:
    """Require every certified-output hash printed in the paper to be manifested.

    The expanded identity table records a historical generator identity first
    and one or more accepted output identities afterward.  Historical source
    snapshots are provenance; every output, however, must resolve through the
    verified archive manifests.
    """

    paper = (root / "paper.tex").read_text(errors="strict")
    start_marker = "The following table is an identity ledger"
    end_marker = r"\section{Post-"
    try:
        start = paper.index(start_marker)
        end = paper.index(end_marker, start)
    except ValueError as error:
        raise ReplayFailure("could not locate the certificate identity ledger") from error

    manifested: set[str] = set()
    for manifest in root.rglob("*.sha256"):
        for line in manifest.read_text(errors="strict").splitlines():
            match = HASH_LINE.match(line)
            if match is not None:
                manifested.add(match.group(1))

    rows = 0
    generators: list[str] = []
    outputs: list[str] = []
    for chunk in paper[start:end].split(r"\addlinespace"):
        hashes = [left + right for left, right in HASH_SPLIT.findall(chunk)]
        if not hashes:
            continue
        if any(len(value) != 64 for value in hashes):
            raise ReplayFailure("malformed SHA-256 in the certificate identity ledger")
        rows += 1
        generators.append(hashes[0])
        outputs.extend(hashes[1:])

    if (rows, len(outputs)) != (115, 119):
        raise ReplayFailure(
            f"certificate identity ledger shape changed: rows={rows}, outputs={len(outputs)}"
        )
    missing = sorted({value for value in outputs if value not in manifested})
    if missing:
        raise ReplayFailure(
            "certified output hashes are absent from verified manifests: " + ", ".join(missing)
        )
    outside_ledger = paper[:start] + paper[end:]
    inline_raw = {
        value
        for value in re.findall(
            r"(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])", outside_ledger
        )
        if re.search(r"[a-f]", value)
    }
    inline_split = {left + right for left, right in HASH_SPLIT.findall(outside_ledger)}
    if any(len(value) != 64 for value in inline_split):
        raise ReplayFailure("malformed split SHA-256 outside the identity ledger")
    inline = inline_raw | inline_split
    if len(inline) != 96:
        raise ReplayFailure(f"inline SHA-256 ledger shape changed: hashes={len(inline)}")
    # Historical generator identities are provenance fields in the identity
    # ledger and need not equal a current manifested source snapshot.  Hashes
    # of manifest/classification containers authenticate those containers
    # directly, rather than appearing as entries inside themselves.  No other
    # unmanifested hash is allowed elsewhere in the manuscript.
    container_digests = {
        digest(path)
        for path in root.rglob("*")
        if path.is_file() and path.suffix in {".sha256", ".tsv"}
    }
    missing_inline = sorted(
        inline - manifested - set(generators) - container_digests
    )
    if missing_inline:
        raise ReplayFailure(
            "inline paper hashes are absent from verified manifests: "
            + ", ".join(missing_inline)
        )
    return rows, len(outputs), len(inline)


def audit_reference_archive(root: Path) -> tuple[int, int, int]:
    """Pin the complete reference directory and map every bibliography key."""

    references = root / "references"
    manifest = references / "references_manifest.sha256"
    if not manifest.is_file():
        raise ReplayFailure("missing complete reference manifest")
    actual = {
        path.resolve()
        for path in references.rglob("*")
        if path.is_file() and path != manifest
    }
    listed: set[Path] = set()
    for line_number, line in enumerate(manifest.read_text(errors="strict").splitlines(), 1):
        match = HASH_LINE.match(line)
        if match is None:
            raise ReplayFailure(f"malformed reference manifest at line {line_number}")
        listed.add(resolve_manifest_path(root, manifest, match.group(2)))
    if actual != listed:
        missing = sorted(str(path.relative_to(root)) for path in actual - listed)
        extra = sorted(str(path.relative_to(root)) for path in listed - actual)
        raise ReplayFailure(
            f"reference manifest is not directory-complete: missing={missing}, extra={extra}"
        )

    source_map = references / "active_paper_sources.tsv"
    with source_map.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"bibkey", "path", "verified_locus"}:
        raise ReplayFailure("malformed active paper source map")
    if any(not all(row.get(field, "").strip() for field in row) for row in rows):
        raise ReplayFailure("active paper source map contains an empty field")
    paper = (root / "paper.tex").read_text(errors="strict")
    bibliography_keys = set(re.findall(r"\\bibitem\{([^}]+)\}", paper))
    cited_keys: set[str] = set()
    for group in re.findall(r"\\cite(?:\[[^]]*\])?\{([^}]+)\}", paper):
        cited_keys.update(key.strip() for key in group.split(","))
    mapped_keys = {row["bibkey"] for row in rows}
    mapped_files = {(references / row["path"]).resolve() for row in rows}
    if bibliography_keys != mapped_keys:
        raise ReplayFailure(
            f"bibliography/source-map mismatch: bibliography={sorted(bibliography_keys)}, "
            f"mapped={sorted(mapped_keys)}"
        )
    if cited_keys != bibliography_keys:
        raise ReplayFailure(
            f"citation/bibliography mismatch: cited={sorted(cited_keys)}, "
            f"bibliography={sorted(bibliography_keys)}"
        )
    if not mapped_files <= actual:
        raise ReplayFailure("active paper source map names an unarchived file")
    if (len(actual), len(bibliography_keys), len(rows)) != (45, 14, 15):
        raise ReplayFailure(
            "reference archive shape changed: "
            f"files={len(actual)}, keys={len(bibliography_keys)}, mappings={len(rows)}"
        )
    return len(actual), len(bibliography_keys), len(rows)


def proof_graph(
    root: Path,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, dict[str, object]],
    set[str],
]:
    """Pair numbered results with proofs and walk dependencies from the main theorem."""

    paper = (root / "paper.tex").read_text(errors="strict")
    token = re.compile(
        r"\\begin\{(theorem|lemma|proposition|corollary|proof)\}"
        r"(?:\[([^]]*)\])?|"
        r"\\end\{(theorem|lemma|proposition|corollary|proof)\}"
    )
    stack: list[tuple[str, str, int]] = []
    results: list[dict[str, object]] = []
    proofs: list[dict[str, object]] = []
    for match in token.finditer(paper):
        begin, title, end = match.groups()
        if begin is not None:
            stack.append((begin, title or "", match.start()))
            continue
        if not stack or stack[-1][0] != end:
            raise ReplayFailure(f"theorem/proof environment mismatch near byte {match.start()}")
        kind, saved_title, start = stack.pop()
        record = {"kind": kind, "title": saved_title, "start": start, "end": match.end()}
        (proofs if kind == "proof" else results).append(record)
    if stack:
        raise ReplayFailure("unterminated theorem/proof environment")
    if len(results) != len(proofs):
        raise ReplayFailure(
            f"theorem/proof count mismatch: results={len(results)}, proofs={len(proofs)}"
        )

    nodes: dict[str, dict[str, object]] = {}
    attached: set[int] = set()
    main_proofs = [proof for proof in proofs if "thm:main" in str(proof["title"])]
    if len(main_proofs) != 1:
        raise ReplayFailure(f"expected one deferred main proof, found {len(main_proofs)}")
    for index, result in enumerate(results):
        statement = paper[int(result["start"]):int(result["end"])]
        labels = re.findall(r"\\label\{([^}]+)\}", statement)
        if len(labels) != 1 or labels[0] in nodes:
            raise ReplayFailure(
                f"numbered result has missing or duplicate label near byte {result['start']}"
            )
        label = labels[0]
        next_start = int(results[index + 1]["start"]) if index + 1 < len(results) else len(paper)
        if label == "thm:main":
            candidates = main_proofs
        else:
            candidates = [
                proof
                for proof in proofs
                if int(result["end"]) <= int(proof["start"]) < next_start
                and "thm:main" not in str(proof["title"])
            ]
        if len(candidates) != 1:
            raise ReplayFailure(
                f"result {label} has {len(candidates)} attached proofs instead of one"
            )
        proof = candidates[0]
        proof_start = int(proof["start"])
        if proof_start in attached:
            raise ReplayFailure(f"proof attached more than once at byte {proof_start}")
        attached.add(proof_start)
        body = statement + paper[proof_start:int(proof["end"])]
        edge_groups = re.findall(r"\\(?:c|C)ref\{([^}]+)\}", body)
        edges = {value.strip() for group in edge_groups for value in group.split(",")}
        nodes[label] = {
            "title": result["title"],
            "body": body,
            "edges": edges,
        }
    if len(attached) != len(proofs):
        raise ReplayFailure(f"unattached proofs remain: {len(proofs) - len(attached)}")

    reachable: set[str] = set()
    frontier = ["thm:main"]
    while frontier:
        label = frontier.pop()
        if label in reachable:
            continue
        reachable.add(label)
        frontier.extend(
            edge for edge in nodes[label]["edges"] if edge in nodes and edge not in reachable
        )
    if "thm:main" not in reachable or len(reachable) < 3:
        raise ReplayFailure("main-theorem dependency walk is unexpectedly empty")
    return results, proofs, nodes, reachable


def audit_proof_structure(root: Path) -> tuple[int, int, int]:
    """Audit proof attachment and unresolved-status leaks on the main proof spine.

    This is deliberately a source-structure check, not a theorem prover: the
    mathematical discharge of hypotheses is stated and proved in the paper.
    """

    results, proofs, nodes, reachable = proof_graph(root)
    forbidden = ("not proved", "unproved", "would have to", "conditional arithmetic", "superseded")
    leaks = [
        (label, phrase)
        for label in reachable
        for phrase in forbidden
        if phrase in str(nodes[label]["body"]).lower()
    ]
    if leaks:
        raise ReplayFailure(
            "explicit unresolved-status phrase is reachable from the main theorem: "
            f"phrases={leaks}"
        )
    return len(results), len(proofs), len(reachable)


def reachable_certificate_logs(root: Path) -> set[str]:
    """Return certificate paths cited on the main-theorem proof spine."""

    _results, _proofs, nodes, reachable = proof_graph(root)
    certificate = re.compile(
        r"certificates/\s*"
        r"((?:[A-Za-z0-9_.-]+\s*/\s*)*[A-Za-z0-9_.-]+\.log)"
    )
    logs = {
        "certificates/" + re.sub(r"\s+", "", match.group(1))
        for label in reachable
        for match in certificate.finditer(
            str(nodes[label]["body"]).replace(r"\_", "_")
        )
    }
    if not logs:
        raise ReplayFailure("main-theorem dependency walk exposes no certificate logs")
    return logs


def audit_bc_caller_ranges(root: Path) -> tuple[int, int, int]:
    """Check the quantified B/C correction domains used by the final bridge.

    The citation graph above cannot detect a theorem whose label is present but
    whose rank interval is too short for its caller.  This audit parses the
    correction propositions themselves, checks their moment offsets and rank
    domains, checks the two analytic tail suppliers, and independently counts
    the odd exponents in the former rank-62-and-above coverage gap.
    """

    paper = (root / "paper.tex").read_text(errors="strict")
    result_pattern = re.compile(
        r"\\begin\{(?:theorem|lemma|proposition|corollary)\}"
        r"(?:\[[^]]*\])?(.*?)"
        r"\\end\{(?:theorem|lemma|proposition|corollary)\}",
        re.DOTALL,
    )
    statements: dict[str, str] = {}
    for match in result_pattern.finditer(paper):
        statement = match.group(0)
        labels = re.findall(r"\\label\{([^}]+)\}", statement)
        if len(labels) == 1:
            statements[labels[0]] = statement

    def require_statement(label: str) -> str:
        try:
            return statements[label]
        except KeyError as error:
            raise ReplayFailure(f"missing B/C caller-range result: {label}") from error

    def domains(label: str) -> list[tuple[int, int]]:
        return [
            (int(low), int(high))
            for low, high in re.findall(
                r"(\d+)\\le\s*b\\le(\d+)", require_statement(label)
            )
        ]

    def correction_offset(label: str) -> int:
        matches = re.findall(
            r"-\\frac12\s+s_\{2q(?:\+(\d+))?\}\\le\s+\\delta",
            require_statement(label),
        )
        if len(matches) != 1:
            raise ReplayFailure(
                f"could not read the unique correction offset in {label}: {matches}"
            )
        return int(matches[0]) if matches[0] else 0

    b_labels = sorted(
        label
        for label in statements
        if re.fullmatch(
            r"prop:bc-b-(?:boundary|[a-z]+-nonboundary)-lower-correction", label
        )
    )
    c_labels = sorted(
        label
        for label in statements
        if re.fullmatch(
            r"prop:bc-c-(?:boundary|[a-z]+-nonboundary)-lower-correction", label
        )
    )
    b_high = "prop:bc-b-twentyninth-high-nonboundary-lower-correction"
    c_overlap = "prop:bc-c-twentyninth-nonboundary-lower-correction"
    if len(b_labels) != 29 or len(c_labels) != 30:
        raise ReplayFailure(
            "B/C correction proposition count changed: "
            f"B_prefix={len(b_labels)}, C_prefix={len(c_labels)}"
        )

    b_offsets = {correction_offset(label) for label in b_labels}
    c_offsets = {correction_offset(label) for label in c_labels}
    if b_offsets != set(range(29)):
        raise ReplayFailure(f"B correction offsets are not 0..28: {sorted(b_offsets)}")
    if c_offsets != set(range(30)):
        raise ReplayFailure(f"C correction offsets are not 0..29: {sorted(c_offsets)}")
    for label in b_labels:
        if domains(label) != [(14, 123)]:
            raise ReplayFailure(f"short B caller range in {label}: {domains(label)}")
    for label in c_labels:
        if domains(label) != [(20, 217)]:
            raise ReplayFailure(f"short C caller range in {label}: {domains(label)}")
    if domains(b_high) != [(42, 123)] or correction_offset(b_high) != 29:
        raise ReplayFailure(
            f"bad high-B twenty-ninth correction domain/offset: "
            f"domains={domains(b_high)}, offset={correction_offset(b_high)}"
        )
    if correction_offset(c_overlap) != 29:
        raise ReplayFailure("bad C offset-29 overlap result")

    active_contract = "prop:bc-active-correction-prefix-contract"
    active_statement = require_statement(active_contract)
    if domains(active_contract) != [(14, 123), (20, 217)]:
        raise ReplayFailure(
            f"bad active B/C contract domains: {domains(active_contract)}"
        )
    if active_statement.count(r"0\le r\le28") != 2:
        raise ReplayFailure("active B/C contract does not state both offset-0..28 ranges")
    active_b_labels = set(b_labels)
    active_c_labels = {
        label for label in c_labels if correction_offset(label) <= 28
    }
    for label in sorted(active_b_labels | active_c_labels):
        if active_statement.count(label) != 1:
            raise ReplayFailure(
                f"active B/C contract table must contain {label} exactly once"
            )
    for label in (b_high, c_overlap):
        if label in active_statement:
            raise ReplayFailure(f"non-load-bearing offset-29 result entered active contract: {label}")

    residual = "prop:post29-bc-residual-closure"
    if domains(residual) != [(14, 123), (20, 217)]:
        raise ReplayFailure(f"bad residual B/C caller domains: {domains(residual)}")
    residual_proof_match = re.search(
        rf"\\label\{{{residual}\}}.*?\\begin\{{proof\}}(.*?)\\end\{{proof\}}",
        paper,
        re.DOTALL,
    )
    if residual_proof_match is None:
        raise ReplayFailure("cannot parse the residual B/C proof")
    residual_proof = residual_proof_match.group(1)
    half_bridge = "prop:post29-bc-half-bridge"
    if rf"\contractref{{{half_bridge}}}" not in residual_proof:
        raise ReplayFailure("residual B/C proof omits its half-stable bridge import")
    if (
        "post_m29_bc_interval_bridge_frontier_gmp.cpp" not in residual_proof
        or r"\mathcal L_m" not in residual_proof
        or r"D_G(2m+1)\ge\mathcal L_m" not in residual_proof
    ):
        raise ReplayFailure("residual B/C proof omits the exact bridge predicate/source")
    for supplier in (
        "prop:post29",
        "prop:post29-bc-layered-mgf",
        "prop:post29-b-tilted-mgf",
        "prop:post29-c-tilted-mgf",
        "prop:post29-bc-interval-direct",
    ):
        if supplier not in residual_proof:
            raise ReplayFailure(f"residual B/C proof omits supplier {supplier}")
    proof_spine = (root / "PUBLICATION_PROOF_SPINE.md").read_text(errors="strict")
    if half_bridge not in proof_spine or "explicit incoming node" not in proof_spine:
        raise ReplayFailure("publication proof spine omits the half-stable bridge edge")
    tail_specs = (
        ("prop:post29-b-tilted-mgf", [(62, 123)]),
        ("prop:post29-c-tilted-mgf", [(62, 217)]),
    )
    for label, expected_domain in tail_specs:
        statement = require_statement(label)
        if domains(label) != expected_domain or not re.search(
            r"n\\ge\s*2b\+29", statement
        ):
            raise ReplayFailure(
                f"bad tilted-tail caller contract in {label}: {domains(label)}"
            )

    # The first-hit value is supplied separately.  To propagate from it to the
    # tilted-tail onset 2b+29, a target odd n uses correction moments only
    # through n+2, hence correction offset n-2b relative to 2b+2.
    newly_closed = 0
    max_required_offset = 0
    for low, high, available_offsets in (
        (62, 123, b_offsets),
        (62, 217, c_offsets),
    ):
        for rank in range(low, high + 1):
            first_hit = max(63, 2 * rank + 1)
            for odd_target in range(first_hit + 2, 2 * rank + 29, 2):
                required = odd_target - 2 * rank
                max_required_offset = max(max_required_offset, required)
                if not set(range(required + 1)) <= available_offsets:
                    raise ReplayFailure(
                        "B/C local bridge lacks a correction prefix at "
                        f"rank={rank}, target={odd_target}, offset={required}"
                    )
                newly_closed += 1
    if (newly_closed, max_required_offset) != (2834, 27):
        raise ReplayFailure(
            "B/C former-gap enumeration changed: "
            f"odd_targets={newly_closed}, max_offset={max_required_offset}"
        )

    arithmetic_log = (
        root / "certificates/post_m29/bc_pieri_powerloss_2q_gmp_certificate.log"
    ).read_text(errors="strict")
    required_log_fragments = (
        "rows=402 total_checks=3904626 max_index=30899",
        "SUMMARY failures=0 rows=402 total_checks=3904626 max_index=30899",
        "__EXIT_STATUS=0",
    )
    if any(fragment not in arithmetic_log for fragment in required_log_fragments):
        raise ReplayFailure("B/C all-row GMP arithmetic transcript has the wrong scope/status")
    arithmetic_source = (
        root / "character_ring_iter/post_m29_bc_pieri_powerloss_gmp.cpp"
    ).read_text(errors="strict")
    for declaration in (
        "int b_rank_lo = 14;",
        "int b_rank_hi = 217;",
        "int c_rank_lo = 20;",
        "int c_rank_hi = 217;",
    ):
        if declaration not in arithmetic_source:
            raise ReplayFailure(
                f"B/C GMP arithmetic generator range changed: missing {declaration}"
            )

    return len(active_b_labels) + len(active_c_labels), newly_closed, 402


def audit_accepted_replay_logs(root: Path) -> int:
    classification = root / "certificates/post_m29/post_m29_artifact_classification.tsv"
    accepted = 0
    with classification.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        for row in rows:
            if row["status"] != "accepted" or row["role"] != "replay":
                continue
            path = classification.parent / row["path"]
            if not path.exists():
                raise ReplayFailure(f"accepted replay is missing: {path}")
            text = path.read_text(errors="strict")
            if "__EXIT_STATUS=0" not in text:
                raise ReplayFailure(f"accepted replay lacks __EXIT_STATUS=0: {path}")
            accepted += 1
    return accepted


def audit_post_m29_artifact_boundary(root: Path) -> int:
    """Require the accepted manifest to equal the accepted classification rows."""

    certificate = root / "certificates/post_m29"
    classification = certificate / "post_m29_artifact_classification.tsv"
    accepted_manifest = certificate / "post_m29_accepted_manifest.sha256"

    classified: set[Path] = set()
    with classification.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if rows.fieldnames != ["path", "role", "status", "proof_use"]:
            raise ReplayFailure("malformed post-m29 artifact-classification header")
        for row in rows:
            if row["status"] != "accepted":
                continue
            path = (certificate / row["path"]).resolve()
            if path in classified:
                raise ReplayFailure(f"duplicate accepted classification row: {row['path']}")
            if not path.exists():
                raise ReplayFailure(f"accepted classified artifact is missing: {path}")
            classified.add(path)

    manifested: set[Path] = set()
    for line_number, line in enumerate(
        accepted_manifest.read_text(errors="strict").splitlines(), 1
    ):
        match = HASH_LINE.match(line)
        if match is None:
            continue
        path = resolve_manifest_path(root, accepted_manifest, match.group(2))
        if path in manifested:
            raise ReplayFailure(
                f"duplicate accepted-manifest artifact at line {line_number}: {path}"
            )
        manifested.add(path)

    missing = sorted(str(path) for path in classified - manifested)
    extra = sorted(str(path) for path in manifested - classified)
    if missing or extra:
        raise ReplayFailure(
            "post-m29 accepted boundary mismatch: "
            f"missing_from_manifest={missing}, unclassified_manifest_entries={extra}"
        )
    return len(classified)


def purge_elf_binaries(root: Path) -> int:
    removed = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        with path.open("rb") as handle:
            magic = handle.read(4)
        if magic == b"\x7fELF":
            path.unlink()
            removed += 1
    for cache in root.rglob("__pycache__"):
        if cache.is_dir():
            shutil.rmtree(cache)
    return removed


def purge_latex_outputs(root: Path) -> int:
    removed = 0
    for stem in ("paper", "paper_full"):
        for suffix in (".aux", ".log", ".out", ".pdf", ".toc"):
            path = root / f"{stem}{suffix}"
            if path.exists():
                path.unlink()
                removed += 1
    return removed


def tail(path: Path, lines: int = 80) -> str:
    content = path.read_text(errors="replace").splitlines()
    return "\n".join(content[-lines:])


class Runner:
    def __init__(self, log_dir: Path, environment: dict[str, str]):
        self.log_dir = log_dir
        self.environment = environment
        self.results: list[tuple[str, float]] = []
        self.replayed_certificates: set[str] = set()

    def run(
        self,
        name: str,
        command: list[str],
        cwd: Path,
        required_text: str | tuple[str, ...] | None = None,
        environment_overrides: dict[str, str] | None = None,
        certificates: tuple[str, ...] = (),
        timeout: int = 7200,
    ) -> Path:
        log_path = self.log_dir / f"{len(self.results) + 1:02d}_{name}.log"
        print(f"[replay] {name} ...", flush=True)
        started = time.monotonic()
        environment = self.environment.copy()
        if environment_overrides:
            environment.update(environment_overrides)
        with log_path.open("w") as output:
            output.write("COMMAND: " + " ".join(command) + "\n")
            if environment_overrides:
                output.write(
                    "ENV: "
                    + " ".join(
                        f"{key}={environment_overrides[key]}"
                        for key in sorted(environment_overrides)
                    )
                    + "\n"
                )
            output.flush()
            try:
                completed = subprocess.run(
                    command,
                    cwd=cwd,
                    env=environment,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired as error:
                raise ReplayFailure(f"{name} timed out after {timeout}s; log: {log_path}") from error
        elapsed = time.monotonic() - started
        if completed.returncode != 0:
            raise ReplayFailure(
                f"{name} exited {completed.returncode}; log: {log_path}\n{tail(log_path)}"
            )
        required = () if required_text is None else (
            (required_text,) if isinstance(required_text, str) else required_text
        )
        log_text = log_path.read_text(errors="replace")
        missing = tuple(fragment for fragment in required if fragment not in log_text)
        if missing:
            raise ReplayFailure(
                f"{name} did not emit {missing!r}; log: {log_path}\n{tail(log_path)}"
            )
        for certificate in certificates:
            relative = Path(certificate)
            if (
                relative.is_absolute()
                or ".." in relative.parts
                or len(relative.parts) not in {1, 2}
                or relative.suffix != ".log"
            ):
                raise ReplayFailure(f"invalid replay certificate basename: {certificate}")
            normalized = (
                Path("certificates/post_m29") / relative
                if len(relative.parts) == 1
                else Path("certificates") / relative
            )
            self.replayed_certificates.add(normalized.as_posix())
        self.results.append((name, elapsed))
        print(f"[replay] {name}: PASS ({elapsed:.2f}s)", flush=True)
        return log_path


def printed_inputs(root: Path, replay_logs: list[Path]) -> list[Path]:
    paths: set[Path] = set()
    for replay_log in replay_logs:
        for line in replay_log.read_text(errors="strict").splitlines():
            match = PRINTED_LOG.match(line)
            if not match:
                continue
            printed = match.group(1)
            path = ((root.parent if printed.startswith("ginibre_q3/") else root) / printed).resolve()
            if not path.exists():
                raise ReplayFailure(f"input printed by {replay_log} is missing: {path}")
            paths.add(path)
    return sorted(paths)


def manifested_inputs(root: Path, manifests: list[Path]) -> list[Path]:
    """Resolve authenticated input logs without consulting superseded replays."""

    paths: set[Path] = set()
    for manifest in manifests:
        entries = 0
        for line_number, line in enumerate(
            manifest.read_text(errors="strict").splitlines(), 1
        ):
            match = HASH_LINE.match(line)
            if match is None:
                raise ReplayFailure(
                    f"malformed input manifest at {manifest}:{line_number}"
                )
            path = resolve_manifest_path(root, manifest, match.group(2))
            if path.suffix != ".log":
                raise ReplayFailure(f"non-log bridge input in {manifest}: {path}")
            paths.add(path)
            entries += 1
        if entries == 0:
            raise ReplayFailure(f"bridge input manifest is empty: {manifest}")
    return sorted(paths)


def classical_claim_groups(
    paths: list[Path],
    allowed: set[str] | None = None,
) -> dict[str, tuple[int, list[Path]]]:
    """Group exact classical source claims by family/rank and needed degree."""

    grouped: dict[str, tuple[int, set[Path]]] = {}
    for path in paths:
        text_value = path.read_text(errors="strict")
        claims: list[tuple[str, int]] = []
        for match in CLASSICAL_DELTA_CLAIM.finditer(text_value):
            claims.append((match.group(1) + match.group(2), int(match.group(3))))
        for pattern in (CLASSICAL_MOMENT_CLAIM, CLASSICAL_DIRECT_MOMENT_CLAIM):
            for match in pattern.finditer(text_value):
                claims.append(("D" + match.group(1), int(match.group(2))))
        for group, degree in claims:
            if allowed is not None and group not in allowed:
                continue
            maximum, sources = grouped.setdefault(group, (0, set()))
            sources.add(path)
            grouped[group] = (max(maximum, degree), sources)
    return {
        group: (maximum, sorted(sources))
        for group, (maximum, sources) in grouped.items()
    }


def replay_classical_source_pairings(
    runner: Runner,
    character: Path,
    stable: str,
    paths: list[Path],
    name_prefix: str,
    allowed: set[str] | None = None,
) -> None:
    """Regenerate claimed finite-rank moments by an independent root-datum path."""

    groups = classical_claim_groups(paths, allowed)
    if not groups:
        raise ReplayFailure(f"{name_prefix} exposed no classical moment claims")
    for group in sorted(groups, key=lambda value: (value[0], int(value[1:]))):
        maximum, sources = groups[group]
        command = [
            "./verify_character_ring_moment_sources_gmp",
            group,
            str(maximum),
            "--stable-moments",
            stable,
        ]
        command.extend(str(path) for path in sources)
        runner.run(
            f"{name_prefix}_{group.lower()}",
            command,
            character,
            required_text=f"SUMMARY group={group[0]}_{group[1:]} ",
        )


def stable_moment_recurrence_replay(runner: Runner, root: Path) -> None:
    """Validate the archived stable table from the proved recurrence first."""

    runner.run(
        "stable_moment_recurrence_comparison",
        [
            "./verify_character_ring_moment_sources_gmp",
            "D4",
            "59",
            "--stable-moments",
            str(root / "references/oeis_A002137_stable.txt"),
            str(
                root
                / "certificates/post_m29/"
                "d4_modular_weyl_exact_moments_certificate.log"
            ),
        ],
        root / "character_ring_iter",
        required_text="stable_recurrence_values_checked=101",
    )


def classical_replays(runner: Runner, root: Path) -> None:
    character = root / "character_ring_iter"
    script = "classical_tail_constants.py"
    simple_flags = [
        "--b2-c2-post-m29-tail-certificate",
        "--b3-c3-post-m29-tail-certificate",
        "--d4-post-m29-tail-certificate",
        "--b4-c4-post-m29-tail-certificate",
        "--d5-post-m29-tail-certificate",
        "--b5-post-m29-tail-certificate",
        "--d6-b6-post-m29-tail-certificate",
        "--d7-b7-c7-post-m29-tail-certificate",
        "--d8-post-m29-tail-certificate",
        "--b8-c9-d9-post-m29-tail-certificate",
    ]
    for flag in simple_flags:
        replay_name = flag.removeprefix("--").replace("-", "_")
        runner.run(
            f"classical_low_rank_exact_{replay_name}",
            [sys.executable, script, flag],
            character,
        )

    runner.run(
        "classical_low_bc_exact_mgf",
        ["./post_m29_bc_layered_mgf_mpfr", "--low-bc-exact-certificate"],
        character,
        required_text=(
            "SUMMARY rows_checked=22 B_rows_checked=6 C_rows_checked=16 "
            "D_rows_checked=0 failures=0"
        ),
        certificates=("post_m29_low_bc_exact_mgf_mpfr_certificate.log",),
    )

    certificate = root / "certificates/post_m29"
    exact_source = certificate / "d4_24_exact_adjoint_sources"
    d10_sources = sorted(exact_source.glob("d4_11_exact_adjoint_*.log"))
    exact_tail_command = ["./verify_d10_exact_tail_gmp"]
    for path in d10_sources:
        exact_tail_command.extend(("--source-log", str(path)))
    runner.run(
        "classical_d10_exact_tail",
        exact_tail_command,
        character,
        required_text="SUMMARY rows_checked=1 failures=0",
        certificates=("post_m29_d10_exact_tail_gmp_certificate.log",),
    )

    d11_mgf_sources = [
        exact_source / "d4_11_exact_adjoint_m4_20.log",
        exact_source / "d4_11_exact_adjoint_m21_26.log",
    ]
    d11_mgf_command = [
        "./post_m29_bc_layered_mgf_mpfr",
        "--d-only",
        "--d-rank-low", "11",
        "--d-rank-high", "11",
        "--allow-d-shallow-wall",
        "--exact-d-square-mgf",
        "--exact-d-defining-mgf",
        "--d-polynomial-k", "3",
        "--d-polynomial-m", "7",
        "--d-square-alpha", "24", "5",
        "--d-square-lambda-fraction", "1", "1",
    ]
    for path in d11_mgf_sources:
        d11_mgf_command.extend(("--d-adjoint-moment-source-log", str(path)))
    runner.run(
        "classical_d11_exact_mgf_onset75",
        d11_mgf_command,
        character,
        required_text=(
            "row=D_11 n=75 polynomial_k=3 polynomial_m=7 "
            "push_moments=reconstructed_from_exact_adjoint_source"
        ),
        certificates=("post_m29_d11_exact_mgf_mpfr_machine_c_certificate.log",),
    )

    d11_bridge_command = [
        "./post_m29_d_stable_window_interval_gmp",
        "--rank", "11",
        "--onset", "75",
    ]
    d11_bridge_sources = sorted(exact_source.glob("d4_11_exact_adjoint_*.log"))
    d11_bridge_sources.extend(sorted(exact_source.glob("d11_exact_adjoint_*.log")))
    for path in d11_bridge_sources:
        d11_bridge_command.extend(("--exact-correction-log", str(path)))
    runner.run(
        "classical_d11_onset75_bridge",
        d11_bridge_command,
        character,
        required_text="SUMMARY ranks_checked=1 chain_steps_checked=34 failures=0",
        certificates=("post_m29_d11_onset75_bridge_gmp_machine_c_certificate.log",),
    )

    bc_manifests = [
        # The shared first-hit sources supply the early consecutive prefixes
        # (B12 Delta_26..39 and the corresponding B13/C13/C14/C17 data).  The
        # row-specific manifests below begin only at their later continuation
        # points, so omitting this accepted input manifest makes the focused
        # onset replay correctly reject a gapped source sequence.
        certificate / "post_m29_input_manifest.sha256",
        certificate / "b11_delta_logs_bonly.sha256",
        certificate / "b12_delta_logs_bonly.sha256",
        certificate / "b13_delta_logs_bonly.sha256",
        certificate / "c13_delta_logs_conly.sha256",
        certificate / "c14_delta_logs_conly.sha256",
        certificate / "c17_delta_logs_conly.sha256",
    ]
    supplier = certificate / "post_m29_lower_bound_supplier_certificate_replay.log"
    bc_inputs = manifested_inputs(root, bc_manifests)
    # The accepted m=28 bridge source supplies C17 Delta_40 and Delta_41,
    # between the shared prefix through Delta_39 and the row-specific shard
    # beginning at Delta_42.
    bc_inputs.extend(
        (
            root
            / "certificates/classical_bridge/raw_logs/ginibre_m28_logs/bc/"
            / "m28_BC17_rank_pair_exact41.rerun2.log",
        )
    )
    bc_inputs = sorted(set(bc_inputs))
    is_d_log = lambda path: re.search(r"(?:^|_)D\d", path.name) is not None
    bc_inputs = [path for path in bc_inputs if not is_d_log(path)]

    replay_classical_source_pairings(
        runner,
        character,
        str(root / "references/oeis_A002137_stable.txt"),
        bc_inputs,
        "classical_low_bc_source_pairing",
    )

    command = [
        sys.executable,
        script,
        "--low-bc-exact-onset-bridge-certificate",
    ]
    for path in bc_inputs:
        command.extend(("--bc-window-delta-log", str(path)))
    runner.run(
        "classical_low_bc_exact_onset_bridge",
        command,
        character,
        required_text="SUMMARY rows_checked=6 chain_steps_checked=34 failures=0",
        certificates=("post_m29_low_bc_exact_onset_bridge_certificate.log",),
    )

    supplier_bc = printed_inputs(root, [supplier])
    supplier_bc = [path for path in supplier_bc if not is_d_log(path)]
    command = [
        sys.executable,
        script,
        "--post-m29-bc-lower-bound-supplier-certificate",
    ]
    for path in supplier_bc:
        command.extend(("--bc-window-delta-log", str(path)))
    runner.run(
        "classical_post_m29_bc_first_hit",
        command,
        character,
        required_text="rows checked: 530",
    )


def exceptional_replays(runner: Runner, root: Path, threads: int) -> None:
    character = root / "character_ring_iter"
    pairing_sources = (
        ("G2", 200, [character / "logs/g2_200.log"]),
        ("F4", 220, [character / "logs/f4_220c.log"]),
        ("E6", 80, [character / "logs/e6_80.log"]),
        ("E7", 70, [character / "logs/e7_70.log"]),
        (
            "E8",
            100,
            [
                character / "logs/e8_70.log",
                root / "references/arxiv_2412_21189_e8_m71_m100.txt",
            ],
        ),
    )
    for group, maximum, sources in pairing_sources:
        runner.run(
            f"exceptional_source_pairing_{group.lower()}",
            [
                "./verify_character_ring_moment_sources_gmp",
                group,
                str(maximum),
                *(str(path) for path in sources),
            ],
            character,
            required_text=(
                f"SUMMARY group={group} moments={maximum + 1} "
                f"source_values={maximum + 1} failures=0"
            ),
        )
    for group, log_name in (
        ("G2", "g2_200.log"),
        ("F4", "f4_220c.log"),
        ("E6", "e6_80.log"),
        ("E7", "e7_70.log"),
        ("E8", "e8_70.log"),
    ):
        runner.run(
            f"exceptional_prefix_{group.lower()}",
            [
                sys.executable,
                "verify_chain.py",
                group,
                str(character / "logs" / log_name),
                str(root / "references"),
            ],
            character,
            required_text="Chain verification: ALL PASS",
        )
    runner.run(
        "exceptional_dunkl_coefficients",
        ["./dunkl_exceptional_coefficients_gmp", "--threads", str(threads)],
        character,
        required_text="All exceptional Dunkl coefficient certificates verified",
    )
    runner.run(
        "exceptional_g2_f4_tails",
        ["./direct_chain_rect_g2_f4_certificate"],
        character,
        required_text="RESULT: G2/F4 RECTANGULAR TAILS PASS",
    )
    runner.run(
        "exceptional_e6_e7_rank_tails",
        [
            "./direct_chain_rect_e6_e7_rank_certificate",
            "--e7-moments",
            str(character / "logs/e7_70.log"),
        ],
        character,
        required_text="RESULT: ALL PASS",
    )
    runner.run(
        "exceptional_e8_negative_bounds",
        [
            "./verify_e8_negative_bounds_gmp",
            "--table",
            str(root / "certificates/exceptional_rect/e8_rect_negative_bounds.tsv"),
            "--factors",
            str(root / "certificates/exceptional_rect/e8_negative_majorant_factors.tsv"),
            "--local-moments",
            str(character / "logs/e8_70.log"),
            "--ancillary-moments",
            str(root / "references/arxiv_2412_21189_e8_m71_m100.txt"),
            "--oeis-moments",
            str(root / "references/oeis_A179663.txt"),
        ],
        character,
        required_text="E8_NEGATIVE_BOUND_AUDIT_GMP: ALL PASS",
    )
    runner.run(
        "exceptional_e8_rectangles",
        [
            "./direct_chain_rect_e8_parallel_certificate",
            "--negative-bounds",
            str(root / "certificates/exceptional_rect/e8_rect_negative_bounds.tsv"),
        ],
        character,
        required_text="RESULT: ALL PASS",
    )
    runner.run(
        "exceptional_e8_finite_bridge",
        [sys.executable, "direct_chain_e8_arxiv_m100_finite_bridge.py"],
        character,
        required_text="extended_chain_diff_67_95_positive: OK",
    )


def cpp_endpoint_replays(runner: Runner, root: Path, threads: int) -> None:
    character = root / "character_ring_iter"
    certificate = root / "certificates/post_m29"
    runner.run(
        "sun_recurrence_derivation_exact",
        [sys.executable, "verify_sun_recurrence_derivation_exact.py"],
        character,
        required_text=(
            "SUMMARY SU(3..5) recurrence derivations verified by exact integer "
            "Bessel/Weyl/differential algebra: ALL PASS"
        ),
    )
    runner.run(
        "sun_repair_gmp",
        ["./verify_sun_repair_certificates_gmp"],
        character,
        required_text="All SU(N) repair certificates verified with exact GMP arithmetic",
    )
    stable = str(root / "references/oeis_A002137_stable.txt")
    runner.run(
        "classical_bc_row_gated_source_bridge",
        [
            "./verify_bc_row_gated_bridge_gmp",
            "--threads",
            str(min(threads, 3)),
            "--stable",
            stable,
        ],
        character,
        required_text=(
            "families=2 ranks=58 row_gated_steps=870 "
            "exact_steps=416 interval_steps=454 failures=0 status=PASS"
        ),
    )
    runner.run(
        "classical_d53_63_row_gated_source_bridge",
        [
            "./post_m29_d_stable_window_interval_gmp",
            "--certificate-d53-63-row-gated",
        ],
        character,
        required_text=(
            "SUMMARY ranks_checked=11 chain_steps_checked=36 failures=0"
        ),
    )
    runner.run(
        "classical_source_m29_weyl",
        [
            "./classical_boundary_source_replay_gmp",
            "--b2-c2-weyl-certificate",
            "--b3-c3-weyl-certificate",
            "--chain-m",
            "29",
            "--stable-moments",
            stable,
        ],
        character,
        required_text="C_3 Chain diff D(29) =",
    )
    runner.run(
        "classical_source_m29_high_rows",
        [
            "./classical_boundary_source_replay_gmp",
            "--bc-lower-range",
            "19",
            "30",
            "--d-lower-range",
            "35",
            "63",
            "--chain-m",
            "29",
            "--stable-moments",
            stable,
        ],
        character,
        required_text="minimum lower margin = B_19",
    )
    runner.run(
        "trace_square_cutoff_mpfr",
        ["./trace_square_cutoff_mpfr"],
        character,
        required_text="first_hit_rows=57 minimum_row=D_281",
        certificates=("classical_trace/first_hit_trace_replay.log",),
    )
    runner.run(
        "d4_j18_denv_counterexample_gmp",
        [
            "./classical_boundary_source_replay_gmp",
            "--d-exact-moment-range",
            "4",
            "4",
            "18",
            "18",
            "--stable-moments",
            stable,
        ],
        character,
        required_text=(
            "Delta_18 = -80781084484224; moment_18 = 578125687744444"
        ),
    )
    layered_mgf_log: Path | None = None
    for name, binary, marker, certificate_name in (
        (
            "b_unitary_square_mpfr",
            "post_m29_b_unitary_square_trace_mpfr",
            "failures=0",
            "post_m29_b_unitary_square_trace_mpfr_certificate.log",
        ),
        (
            "bc_tilted_mgf_mpfr",
            "post_m29_bc_tilted_mgf_mpfr",
            "failures=0",
            "post_m29_bc_tilted_mgf_mpfr_certificate.log",
        ),
        (
            "bc_layered_mgf_mpfr",
            "post_m29_bc_layered_mgf_mpfr",
            "failures=0",
            "post_m29_bc_layered_mgf_mpfr_certificate.log",
        ),
    ):
        replay_log = runner.run(
            name,
            [f"./{binary}"],
            character,
            required_text=marker,
            certificates=(certificate_name,),
        )
        if name == "bc_layered_mgf_mpfr":
            layered_mgf_log = replay_log
    exact_source = certificate / "d4_24_exact_adjoint_sources"
    classical_bridge = root / "certificates/classical_bridge/raw_logs"
    rank4 = classical_bridge / "ginibre_m27_logs/rank4"
    d4_crt_command = [sys.executable, "verify_d4_modular_weyl_moments.py"]
    for prime in (5433, 5507, 5549, 5643, 5783):
        d4_crt_command.extend(
            (
                "--mod-log",
                str(rank4 / f"m27_D4_rank4_mod_prime_922337203685477{prime}.log"),
            )
        )
    d4_crt_command.extend(
        (
            "--reconstruction-log",
            str(rank4 / "m27_D4_rank4_modular_weyl.log"),
            "--reconstruction-check-log",
            str(rank4 / "m27_D4_rank4_modular_weyl.rerun_check.log"),
            "--stable-moments",
            stable,
        )
    )
    d4_cross_checks = sorted(exact_source.glob("d4_11_exact_adjoint_*.log"))
    d4_cross_checks.extend(sorted(exact_source.glob("d4_7_exact_adjoint_m39_*.log")))
    d4_cross_checks.extend(sorted(exact_source.glob("d4_7_exact_adjoint_m40_*.log")))
    for path in d4_cross_checks:
        d4_crt_command.extend(("--source-cross-check-log", str(path)))
    d4_crt_log = runner.run(
        "d4_modular_weyl_exact_moments",
        d4_crt_command,
        character,
        required_text="SUMMARY primes=5 moments=60 cross_checks=37 failures=0",
    )
    runner.run(
        "d4_character_pairing_m0_m59",
        [
            "./verify_character_ring_moment_sources_gmp",
            "D4",
            "59",
            "--stable-moments",
            stable,
            str(d4_crt_log),
        ],
        character,
        required_text="SUMMARY group=D_4 moments=60 source_values=120 failures=0",
    )
    exact_source_logs = sorted(exact_source.glob("*exact_adjoint*.log"))
    replay_classical_source_pairings(
        runner,
        character,
        stable,
        exact_source_logs,
        "classical_d4_d24_source_pairing",
        {f"D{rank}" for rank in range(4, 25)},
    )
    d4_24_command = [
        "./post_m29_d_stable_window_interval_gmp",
        "--certificate-d4-24",
        "--exact-correction-log",
        str(d4_crt_log),
    ]
    for path in exact_source_logs:
        d4_24_command.extend(("--exact-correction-log", str(path)))
    runner.run(
        "d4_d24_exact_prefix_interval_gmp",
        d4_24_command,
        character,
        required_text="SUMMARY ranks_checked=21 chain_steps_checked=530 failures=0",
        certificates=("post_m29_d4_24_exact_prefix_interval_gmp_certificate.log",),
    )
    d12_24_tail_command = [
        "./post_m29_bc_layered_mgf_mpfr",
        "--d12-24-exact-moment-certificate",
    ]
    d12_24_sources = sorted(exact_source.glob("d12_*exact_adjoint*.log"))
    d12_24_sources.extend(sorted(exact_source.glob("d18*exact_adjoint*.log")))
    for path in d12_24_sources:
        d12_24_tail_command.extend(
            ("--d-adjoint-moment-source-log", str(path))
        )
    runner.run(
        "d12_d24_exact_mgf_onset63_mpfr",
        d12_24_tail_command,
        character,
        required_text=(
            "SUMMARY rows_checked=13 B_rows_checked=0 C_rows_checked=0 "
            "D_rows_checked=13 failures=0"
        ),
        certificates=("post_m29_d12_24_exact_mgf_mpfr_certificate.log",),
    )
    d19_artifacts = certificate / "d19_exact_mgf"
    d19_moment_logs = [
        d19_artifacts / f"exact_D19_m{moment}_machine_c.log"
        for moment in range(19, 33)
    ] + [
        d19_artifacts / "exact_D15_23_m33.log",
        d19_artifacts / "exact_D15_23_m34.log",
    ]
    d19_bridge_logs = [
        certificate / "d19_delta_logs/d19_delta_19_42_from_m29_D16_63_lower.log",
        certificate / "d19_delta_logs/d19_delta_27.log",
        d19_artifacts / "exact_D19_m39.log",
    ]
    replay_classical_source_pairings(
        runner,
        character,
        stable,
        sorted(set(d19_moment_logs + d19_bridge_logs)),
        "classical_d19_source_pairing",
        {"D19"},
    )
    d19_tail_command = [
        "./post_m29_bc_layered_mgf_mpfr",
        "--d-only",
        "--d-rank-low",
        "19",
        "--d-rank-high",
        "19",
        "--allow-d-shallow-wall",
        "--exact-d-square-mgf",
        "--exact-d-defining-mgf",
        "--d-polynomial-k",
        "6",
        "--d-polynomial-m",
        "10",
        "--d-square-alpha",
        "9",
        "2",
        "--d-square-lambda-fraction",
        "1",
        "1",
    ]
    for source_log in d19_moment_logs:
        d19_tail_command.extend(
            ["--d-adjoint-moment-source-log", str(source_log)]
        )
    runner.run(
        "d19_exact_mgf_onset63_mpfr",
        d19_tail_command,
        character,
        required_text=(
            "row=D_19 n=63 polynomial_k=6 polynomial_m=10 "
            "push_moments=reconstructed_from_exact_adjoint_source"
        ),
    )
    runner.run(
        "d19_correction_bridge_gmp",
        [
            "./post_m29_d_stable_window_interval_gmp",
            "--rank",
            "19",
            "--onset",
            "63",
            "--exact-correction-log",
            str(d19_bridge_logs[0]),
            "--exact-correction-log",
            str(d19_bridge_logs[1]),
            "--exact-correction-log",
            str(d19_bridge_logs[2]),
        ],
        character,
        required_text="SUMMARY ranks_checked=1 chain_steps_checked=24 failures=0",
    )
    runner.run(
        "d25_d52_exact_mgf_mpfr",
        ["./post_m29_bc_layered_mgf_mpfr", "--d-low-certificate"],
        character,
        required_text=(
            "SUMMARY rows_checked=28 B_rows_checked=0 C_rows_checked=0 "
            "D_rows_checked=28 failures=0"
        ),
        certificates=("post_m29_d25_52_exact_mgf_mpfr_machine_c_certificate.log",),
    )
    runner.run(
        "d25_d52_correction_interval_gmp",
        ["./post_m29_d_stable_window_interval_gmp", "--certificate-d25-52"],
        character,
        required_text="SUMMARY ranks_checked=28 chain_steps_checked=707 failures=0",
        certificates=("post_m29_d25_52_correction_interval_gmp_machine_c_certificate.log",),
    )
    runner.run(
        "d53_d295_a_free_frontier_gmp",
        [
            "./post_m29_d_interval_bridge_frontier_gmp",
            "--rank-lo",
            "53",
            "--rank-hi",
            "295",
            "--chain-lo",
            "30",
            "--onset-log",
            str(layered_mgf_log),
        ],
        character,
        required_text=(
            "minimum_a_free_slack=1"
        ),
        certificates=("post_m29_d53_295_short_frontier_gmp_machine_c_certificate.log",),
    )


def bc_residual_certificate_replays(runner: Runner, root: Path, threads: int) -> None:
    """Recompute every main-proof B/C residual certificate omitted by the old replay."""

    character = root / "character_ring_iter"
    certificate = root / "certificates/post_m29"
    runner.run(
        "bc_cheb80_negative_tail_gmp",
        [
            "./post_m29_bc_cheb80_negative_tail_gmp",
            "--bc-delta-log",
            str(certificate / "m30_BC13_18_pair_range_exact39.log"),
        ],
        character,
        required_text=(
            "B/C degree-8 cutoff-80 Chebyshev negative-tail GMP verifier",
            "SUMMARY failures=0 worst_row=B_217",
        ),
        certificates=("bc_cheb80_negative_tail_gmp_certificate.log",),
    )
    interval_log = runner.run(
        "bc_interval_formula_mpfr",
        ["./post_m29_bc_interval_formula_mpfr"],
        character,
        required_text=(
            "B/C interval fixed-template direct-tail MPFR onset verifier",
            "SUMMARY failures=0 worst_onset=B_20 onset=309",
        ),
        certificates=("bc_interval_formula_mpfr_certificate.log",),
    )
    runner.run(
        "bc_lambda_half_bridge_gmp",
        [
            "./post_m29_bc_interval_bridge_frontier_gmp",
            "--onset-log",
            str(interval_log),
            "--lambda-num",
            "1",
            "--lambda-den",
            "2",
            "--progress",
            "--progress-step",
            "100",
        ],
        character,
        required_text=(
            "lambda=1/2",
            "SUMMARY failures=0 worst_m=15447 worst_odd_n=30897",
        ),
        # The verifier keeps one large exact state map per OpenMP worker.
        # Letting a higher-core, lower-memory host use every advertised CPU
        # can exhaust RAM (48 workers reached about 35 GiB before m=13300 on
        # a 47-GiB host).  The documented primary host has 32 workers, so this
        # cap preserves its replay while making the contract portable.
        environment_overrides={"OMP_NUM_THREADS": str(min(threads, 32))},
        certificates=("bc_lambda_half_bridge_gmp_certificate.log",),
        timeout=28800,
    )
    runner.run(
        "bc_pieri_powerloss_2q_gmp",
        [
            "./post_m29_bc_pieri_powerloss_gmp",
            "--onset-log",
            str(interval_log),
            "--loss-q-mult",
            "2",
            "--loss-add",
            "0",
            "--progress",
        ],
        character,
        required_text=(
            "VERIFIES exact inequality: 2^(A*q+B+1)*binom(j,q)*s_{j-q} <= s_j",
            "SUMMARY failures=0 rows=402 total_checks=3904626 max_index=30899 "
            "loss_q_mult=2 loss_add=0",
        ),
        certificates=("bc_pieri_powerloss_2q_gmp_certificate.log",),
        timeout=28800,
    )

    simple_frontiers = (
        (
            "b14_j37_badshape",
            "post_m29_bc_b14_j37_badshape_gmp",
            "B14_J37_BADSHAPE_GMP q=15 j=37 width_wall=30",
            "bc_b14_j37_badshape_gmp_certificate.log",
        ),
        (
            "b_eighth_frontier",
            "post_m29_bc_b_eighth_frontier_gmp",
            "B_EIGHTH_FRONTIER_GMP cases=3",
            "bc_b_eighth_frontier_gmp_certificate.log",
        ),
        (
            "ninth_frontier",
            "post_m29_bc_ninth_frontier_gmp",
            "BC_NINTH_FRONTIER_GMP b_cases=4 c_arithmetic_cases=2",
            "bc_ninth_frontier_gmp_certificate.log",
        ),
        (
            "tenth_frontier",
            "post_m29_bc_tenth_frontier_gmp",
            "BC_TENTH_FRONTIER_GMP b_cases=5 c_arithmetic_cases=5",
            "bc_tenth_frontier_gmp_certificate.log",
        ),
        (
            "eleventh_frontier",
            "post_m29_bc_eleventh_frontier_gmp",
            "BC_ELEVENTH_FRONTIER_GMP b_cases=7 c_arithmetic_cases=7",
            "bc_eleventh_frontier_gmp_certificate.log",
        ),
        (
            "twelfth_frontier",
            "post_m29_bc_twelfth_frontier_gmp",
            "BC_TWELFTH_FRONTIER_GMP b_cases=8 c_arithmetic_cases=10",
            "bc_twelfth_frontier_gmp_certificate.log",
        ),
        (
            "thirteenth_frontier",
            "post_m29_bc_thirteenth_frontier_gmp",
            "BC_THIRTEENTH_FRONTIER_GMP b_cases=10 c_arithmetic_cases=12",
            "bc_thirteenth_frontier_gmp_certificate.log",
        ),
        (
            "fourteenth_frontier",
            "post_m29_bc_fourteenth_frontier_gmp",
            "BC_FOURTEENTH_FRONTIER_GMP b_cases=11 c_arithmetic_cases=15",
            "bc_fourteenth_frontier_gmp_certificate.log",
        ),
        (
            "fifteenth_frontier",
            "post_m29_bc_fifteenth_frontier_gmp",
            "BC_FIFTEENTH_FRONTIER_GMP b_cases=13 c_fpf_cases=1 c_arithmetic_cases=16",
            "bc_fifteenth_frontier_gmp_certificate.log",
        ),
        (
            "sixteenth_frontier",
            "post_m29_bc_sixteenth_frontier_gmp",
            "BC_SIXTEENTH_FRONTIER_GMP b_cases=14 c_fpf_cases=2 c_arithmetic_cases=18",
            "bc_sixteenth_frontier_gmp_certificate.log",
        ),
        (
            "seventeenth_frontier",
            "post_m29_bc_seventeenth_frontier_gmp",
            "BC_SEVENTEENTH_FRONTIER_GMP b_cases=16 c_fpf_cases=3 c_arithmetic_cases=20",
            "bc_seventeenth_frontier_gmp_certificate.log",
        ),
        (
            "eighteenth_frontier",
            "post_m29_bc_eighteenth_frontier_gmp",
            "BC_EIGHTEENTH_FRONTIER_GMP b_cases=17 c_fpf_cases=3 c_arithmetic_cases=22",
            "bc_eighteenth_frontier_gmp_certificate.log",
        ),
        (
            "nineteenth_frontier",
            "post_m29_bc_nineteenth_frontier_gmp",
            "BC_NINETEENTH_FRONTIER_GMP b_cases=19 c_fpf_cases=4 c_arithmetic_cases=24",
            "bc_nineteenth_frontier_gmp_certificate.log",
        ),
        (
            "twentieth_frontier",
            "post_m29_bc_twentieth_frontier_gmp",
            "BC_TWENTIETH_FRONTIER_GMP b_cases=20 c_fpf_cases=5 c_arithmetic_cases=25",
            "bc_twentieth_frontier_gmp_certificate.log",
        ),
        (
            "twentyfirst_frontier",
            "post_m29_bc_twentyfirst_frontier_gmp",
            "BC_TWENTYFIRST_FRONTIER_GMP b_cases=22 c_fpf_cases=6 c_arithmetic_cases=27",
            "bc_twentyfirst_frontier_gmp_certificate.log",
        ),
        (
            "twentysecond_frontier",
            "post_m29_bc_twentysecond_frontier_gmp",
            "BC_TWENTYSECOND_FRONTIER_GMP b_cases=23 c_fpf_cases=7 c_arithmetic_cases=28",
            "bc_twentysecond_frontier_gmp_certificate.log",
        ),
    )
    for name, binary, scope, certificate_name in simple_frontiers:
        runner.run(
            f"bc_{name}_gmp",
            [f"./{binary}"],
            character,
            required_text=(scope, "SUMMARY failures=0"),
            certificates=(certificate_name,),
            # H21 takes about 6,600 seconds on the documented 32-thread
            # host, leaving too little headroom under the generic two-hour
            # guard; H22 is larger still.  These are exact finite replays,
            # so give the unsharded H8--H22 family the same explicit budget
            # as the other long GMP stages.
            timeout=28800,
        )

    sharded_frontiers = (
        (
            "twentythird",
            15,
            39,
            ((15, 17), (18, 20), (21, 23), (24, 26), (27, 29), (30, 32),
             (33, 35), (36, 39)),
            8,
            30,
        ),
        (
            "twentyfourth",
            15,
            41,
            ((15, 17), (18, 20), (21, 23), (24, 26), (27, 29), (30, 32),
             (33, 35), (36, 38), (39, 41)),
            9,
            31,
        ),
        (
            "twentyfifth",
            15,
            43,
            ((15, 17), (18, 20), (21, 23), (24, 26), (27, 29), (30, 32),
             (33, 35), (36, 38), (39, 41), (42, 43)),
            10,
            33,
        ),
        (
            "twentysixth",
            15,
            45,
            ((15, 17), (18, 20), (21, 23), (24, 26), (27, 29), (30, 32),
             (33, 35), (36, 38), (39, 41), (42, 43), (44, 44), (45, 45)),
            11,
            34,
        ),
        (
            "twentyseventh",
            15,
            47,
            ((15, 17), (18, 20), (21, 23), (24, 26), (27, 29), (30, 32),
             *((rank, rank) for rank in range(33, 48))),
            12,
            36,
        ),
        (
            "twentyeighth",
            15,
            41,
            ((15, 17), (18, 20), (21, 23), *((rank, rank) for rank in range(24, 42))),
            13,
            37,
        ),
    )
    for ordinal, q_low, q_high, shards, c_fpf_cases, c_arithmetic_cases in sharded_frontiers:
        covered = [rank for low, high in shards for rank in range(low, high + 1)]
        if covered != list(range(q_low, q_high + 1)):
            raise ReplayFailure(f"noncontiguous B-frontier shard plan: {ordinal}")
        binary = f"post_m29_bc_{ordinal}_frontier_gmp"
        header = f"BC_{ordinal.upper()}_FRONTIER_GMP"
        for low, high in shards:
            case_count = high - low + 1
            runner.run(
                f"bc_{ordinal}_frontier_b_{low}_{high}",
                [f"./{binary}"],
                character,
                required_text=(
                    f"{header} b_cases={case_count} c_fpf_cases={c_fpf_cases} "
                    f"c_arithmetic_cases={c_arithmetic_cases} run_b=1 run_c=0 "
                    f"b_q_min={low} b_q_max={high}",
                    "SUMMARY failures=0",
                ),
                environment_overrides={
                    "RUN_B": "1",
                    "RUN_C": "0",
                    "B_Q_MIN": str(low),
                    "B_Q_MAX": str(high),
                    "OMP_NUM_THREADS": str(min(threads, 8)),
                },
            )
        runner.run(
            f"bc_{ordinal}_frontier_c",
            [f"./{binary}"],
            character,
            required_text=(
                f"{header} b_cases=0 c_fpf_cases={c_fpf_cases} "
                f"c_arithmetic_cases={c_arithmetic_cases} run_b=0 run_c=1",
                "SUMMARY failures=0",
            ),
            environment_overrides={
                "RUN_B": "0",
                "RUN_C": "1",
                "OMP_NUM_THREADS": str(min(threads, 64)),
            },
            certificates=(f"bc_{ordinal}_frontier_gmp_certificate.log",),
        )

    for name, binary, scope, certificate_name in (
        (
            "twentyeighth_b_absorption",
            "post_m29_bc_twentyeighth_b_absorption_gmp",
            "BC_TWENTYEIGHTH_B_ABSORPTION_GMP h=28 lower_boxes=56 "
            "one_lower_choices=112 two_lower_choices=1596",
            "bc_twentyeighth_b_absorption_gmp_certificate.log",
        ),
        (
            "twentyninth_c_frontier",
            "post_m29_bc_twentyninth_c_frontier_gmp",
            "BC_TWENTYNINTH_C_FRONTIER_GMP c_fpf_cases=14 c_arithmetic_cases=39",
            "bc_twentyninth_c_frontier_gmp_certificate.log",
        ),
        (
            "twentyninth_b_absorption",
            "post_m29_bc_twentyninth_b_absorption_gmp",
            "BC_TWENTYNINTH_B_ABSORPTION_GMP h=29 lower_boxes=58 "
            "one_lower_choices=116 two_lower_choices=1711",
            "bc_twentyninth_b_absorption_gmp_certificate.log",
        ),
    ):
        runner.run(
            f"bc_{name}_gmp",
            [f"./{binary}"],
            character,
            required_text=(scope, "SUMMARY failures=0"),
            certificates=(certificate_name,),
        )


def audit_replayed_certificate_coverage(root: Path, replayed: set[str]) -> int:
    """Require exact equality between theorem-reachable and recomputed certificates."""

    reachable = reachable_certificate_logs(root)
    classification = root / "certificates/post_m29/post_m29_artifact_classification.tsv"
    statuses: dict[str, str] = {}
    with classification.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (Path("certificates/post_m29") / row["path"]).as_posix()
            if name in statuses:
                raise ReplayFailure(f"duplicate post-m29 classification path: {name}")
            statuses[name] = row["status"]
    post_m29 = {
        name for name in reachable
        if name.startswith("certificates/post_m29/")
    }
    unclassified = sorted(post_m29 - statuses.keys())
    nonaccepted = sorted(
        f"{name}:{statuses[name]}"
        for name in post_m29 & statuses.keys()
        if statuses[name] != "accepted"
    )
    manifested: set[Path] = set()
    for manifest in root.rglob("*.sha256"):
        for line in manifest.read_text(errors="strict").splitlines():
            match = HASH_LINE.match(line)
            if match is not None:
                manifested.add(resolve_manifest_path(root, manifest, match.group(2)))
    unauthenticated = sorted(
        name for name in reachable
        if (root / name).resolve() not in manifested
    )
    missing = sorted(reachable - replayed)
    extra = sorted(replayed - reachable)
    if unclassified or nonaccepted or unauthenticated or missing or extra:
        raise ReplayFailure(
            "main-theorem certificate replay mismatch: "
            f"unclassified={unclassified}, nonaccepted={nonaccepted}, "
            f"unauthenticated={unauthenticated}, "
            f"not_replayed={missing}, replayed_but_not_reachable={extra}"
        )
    return len(reachable)


def audit_replay_input_coverage(root: Path) -> int:
    """Require every in-tree file passed to a proof-spine stage to be manifested."""

    class Collector:
        def __init__(self) -> None:
            self.calls: list[tuple[list[str], Path]] = []

        def run(
            self,
            _name: str,
            command: list[str],
            cwd: Path,
            **_kwargs: object,
        ) -> None:
            self.calls.append((command, Path(cwd)))

    collector = Collector()
    stable_moment_recurrence_replay(collector, root)  # type: ignore[arg-type]
    cpp_endpoint_replays(collector, root, 1)  # type: ignore[arg-type]
    bc_residual_certificate_replays(collector, root, 1)  # type: ignore[arg-type]
    classical_replays(collector, root)  # type: ignore[arg-type]
    exceptional_replays(collector, root, 1)  # type: ignore[arg-type]

    manifested: set[Path] = set()
    for manifest in root.rglob("*.sha256"):
        for line in manifest.read_text(errors="strict").splitlines():
            match = HASH_LINE.match(line)
            if match is not None:
                manifested.add(resolve_manifest_path(root, manifest, match.group(2)))

    input_suffixes = {".cpp", ".log", ".md", ".py", ".sha256", ".tex", ".tsv", ".txt"}
    consumed: set[Path] = {
        (root / "paper.tex").resolve(),
        (root / "paper_full.tex").resolve(),
        (root / "character_ring_iter/Makefile").resolve(),
    }
    for command, cwd in collector.calls:
        for argument in command:
            if not isinstance(argument, str) or argument.startswith("-"):
                continue
            raw = Path(argument)
            candidates = [raw] if raw.is_absolute() else [cwd / raw, root / raw]
            for candidate in candidates:
                if candidate.is_file() and root in candidate.resolve().parents:
                    resolved = candidate.resolve()
                    if resolved.suffix in input_suffixes:
                        consumed.add(resolved)
                    break

    missing = sorted(path.relative_to(root) for path in consumed - manifested)
    if missing:
        raise ReplayFailure(
            "proof-spine inputs are absent from verified manifests: "
            + ", ".join(map(str, missing))
        )
    return len(consumed)


def build_latex_document(runner: Runner, root: Path, stem: str) -> int:
    for pass_number in range(1, 4):
        runner.run(
            f"{stem}_pdflatex_pass_{pass_number}",
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"],
            root,
        )
    paper_log = root / f"{stem}.log"
    paper_pdf = root / f"{stem}.pdf"
    if not paper_pdf.is_file() or paper_pdf.stat().st_size == 0:
        raise ReplayFailure(f"pdflatex did not produce a nonempty {stem}.pdf")
    log_text = paper_log.read_text(errors="replace")
    if LATEX_FAILURE.search(log_text):
        raise ReplayFailure(f"{stem} build retained unresolved references; see {paper_log}")
    page_matches = LATEX_PAGES.findall(log_text)
    if not page_matches:
        raise ReplayFailure(f"could not read the page count from {paper_log}")
    return int(page_matches[-1])


def build_papers(runner: Runner, root: Path) -> None:
    main_pages = build_latex_document(runner, root, "paper")
    full_pages = build_latex_document(runner, root, "paper_full")
    if full_pages <= main_pages:
        raise ReplayFailure(
            f"detailed computational supplement is incomplete: main={main_pages}, full={full_pages}"
        )
    print(
        f"[replay] document split: PASS (main={main_pages} pages, "
        f"supplement={full_pages} pages)",
        flush=True,
    )


def preflight(skip_paper: bool = False) -> None:
    required = ["g++", "make"]
    if not skip_paper:
        required.append("pdflatex")
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise ReplayFailure("missing required tools: " + ", ".join(missing))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--keep", action="store_true", help="keep the isolated tree after success")
    parser.add_argument("--skip-paper", action="store_true", help="skip only the final PDF rebuild")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")

    source = Path(__file__).resolve().parent
    workspace = Path(tempfile.mkdtemp(prefix="ginibre-q3-clean-room-"))
    isolated = workspace / "ginibre_q3"
    logs = workspace / "replay-logs"
    logs.mkdir()
    succeeded = False
    started = time.monotonic()
    try:
        preflight(args.skip_paper)
        print(f"[replay] isolated workspace: {workspace}", flush=True)
        shutil.copytree(source, isolated, symlinks=True)
        manifest_count, entry_count = verify_manifests(isolated)
        print(
            f"[replay] artifact integrity: PASS ({manifest_count} manifests, "
            f"{entry_count} entries)",
            flush=True,
        )
        reference_files, bibliography_keys, source_mappings = audit_reference_archive(isolated)
        print(
            f"[replay] reference archive: PASS ({reference_files} files, "
            f"{bibliography_keys} bibliography keys, {source_mappings} mappings)",
            flush=True,
        )
        result_count, proof_count, reachable_count = audit_proof_structure(isolated)
        print(
            f"[replay] proof structure: PASS ({result_count} results, {proof_count} proofs, "
            f"{reachable_count} main-theorem-reachable results)",
            flush=True,
        )
        bc_corrections, bc_closed_odds, bc_arithmetic_rows = audit_bc_caller_ranges(isolated)
        print(
            f"[replay] B/C caller-range closure: PASS ({bc_corrections} active correction "
            f"propositions, {bc_closed_odds} formerly uncovered odd targets, "
            f"{bc_arithmetic_rows} exact-arithmetic rows)",
            flush=True,
        )
        ledger_rows, ledger_outputs, inline_hashes = audit_certificate_identity_ledger(isolated)
        print(
            f"[replay] certificate identity ledger: PASS ({ledger_rows} rows, "
            f"{ledger_outputs} table outputs, {inline_hashes} inline hashes)",
            flush=True,
        )
        accepted_artifacts = audit_post_m29_artifact_boundary(isolated)
        print(
            f"[replay] post-m29 artifact boundary: PASS ({accepted_artifacts} accepted artifacts)",
            flush=True,
        )
        replay_count = audit_accepted_replay_logs(isolated)
        print(f"[replay] accepted replay status: PASS ({replay_count} logs)", flush=True)
        replay_inputs = audit_replay_input_coverage(isolated)
        print(
            f"[replay] proof-spine input coverage: PASS ({replay_inputs} files)",
            flush=True,
        )
        elf_count = purge_elf_binaries(isolated)
        print(f"[replay] purged {elf_count} pre-existing ELF binaries", flush=True)
        latex_count = purge_latex_outputs(isolated)
        print(f"[replay] purged {latex_count} pre-existing LaTeX outputs", flush=True)

        environment = os.environ.copy()
        environment.update(
            {
                "LC_ALL": "C",
                "PYTHONHASHSEED": "0",
                "OMP_NUM_THREADS": str(args.threads),
                "SOURCE_DATE_EPOCH": "946684800",
            }
        )
        runner = Runner(logs, environment)
        runner.run(
            "build_cpp_replayers",
            ["make", "-j", str(min(args.threads, 8)), "replay-build"],
            isolated / "character_ring_iter",
        )
        stable_moment_recurrence_replay(runner, isolated)
        cpp_endpoint_replays(runner, isolated, args.threads)
        bc_residual_certificate_replays(runner, isolated, args.threads)
        classical_replays(runner, isolated)
        exceptional_replays(runner, isolated, args.threads)
        reachable_certificates = audit_replayed_certificate_coverage(
            isolated, runner.replayed_certificates
        )
        print(
            f"[replay] main-theorem certificate coverage: PASS "
            f"({reachable_certificates} reachable logs recomputed)",
            flush=True,
        )
        if not args.skip_paper:
            build_papers(runner, isolated)

        elapsed = time.monotonic() - started
        print(
            f"[replay] ALL GINIBRE Q3 CHECKS PASSED: {len(runner.results)} executable "
            f"stages, {manifest_count} manifests, {entry_count} hashes, {elapsed:.2f}s",
            flush=True,
        )
        succeeded = True
        return 0
    except (ReplayFailure, OSError) as error:
        print(f"[replay] FAILURE: {error}", file=sys.stderr, flush=True)
        print(f"[replay] preserved failed workspace: {workspace}", file=sys.stderr, flush=True)
        return 1
    finally:
        if succeeded and not args.keep:
            shutil.rmtree(workspace)
        elif succeeded:
            print(f"[replay] preserved workspace: {workspace}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
