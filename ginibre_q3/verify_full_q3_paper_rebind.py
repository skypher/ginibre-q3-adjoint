#!/usr/bin/env python3
"""Verify the historical rebind and the live manifested ``paper.tex``.

The 200-stage replay must be collected before the edit is made.  Its archived
execution-source snapshot supplies ``--execution-root``.  This verifier then
proves that the first publication paper differs from those executed bytes by
one fixed prose replacement and removes exactly one cyclic theorem edge.  For
the live revision it additionally permits only the explicitly audited B/C
correction-contract result and residual-closure supplier edges, requires the
half-stable bridge import, and validates the live manifests.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path


BASELINE_SHA256 = "6b2f273aa082f105df4bf60f8ba7047c62cafac88a3bc2c0a5cffdca95607d27"
FINAL_SHA256 = "d2fb8944ec081f47a1bfcf558bb65adb01936b0ac1d0dbbb1f561d63d7197010"
BASELINE_REPLAY_MANIFEST_SHA256 = (
    "030229492cda95a0bb9c1d3183c3880880c2e29e8ee9e2cec717d32733b72ed5"
)
FINAL_REPLAY_MANIFEST_SHA256 = (
    "bfe68eaf0ee4a5fe8263902b1145804a630381a54bfd4a9b4bb3135517196874"
)
BASELINE_FULL_MANIFEST_SHA256 = (
    "6270eabf52062b8478238c5128a914ae0f07d5b31292e09d4c03561d49e278a8"
)
FINAL_FULL_MANIFEST_SHA256 = (
    "1cbd33b571e990c8b8bce47a1bd549129dea6cc89d965ec5951ec6f8f265fb59"
)
REPLAY_MANIFEST = Path("replay_sources.sha256")
FULL_MANIFEST = Path("certificates/full_q3/full_q3_source_manifest.sha256")

OLD_PASSAGE = """The
$B/C$ first-hit replay supplies only this one Chain step; its all-later tails
are the separate propositions cited in \\cref{thm:classical}."""

NEW_PASSAGE = """The $B/C$ first-hit replay supplies only this one Chain step and is not used as
an all-later tail. The all-later claims are proved separately below by the
interval, layered-MGF, tilted-MGF, and correction-bridge propositions."""

RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")
RESULT_PATTERN = re.compile(
    r"\\begin\{(" + "|".join(RESULT_ENVS) + r")\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
IMMEDIATE_PROOF_PATTERN = re.compile(
    r"\s*(\\begin\{proof\}(?:\[[^]]*\])?.*?\\end\{proof\})",
    re.DOTALL,
)
MAIN_PROOF_PATTERN = re.compile(
    r"\\begin\{proof\}\[Proof of \\cref\{thm:main\}\].*?\\end\{proof\}",
    re.DOTALL,
)
RESULT_REFERENCE_PATTERN = re.compile(
    r"\\(?:[Cc]ref|ref|contractref)\{([^}]+)\}"
)
RESULT_LABEL_PREFIXES = ("thm:", "prop:", "lem:", "cor:")


class RebindFailure(RuntimeError):
    """A fail-closed rebind invariant was violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RebindFailure(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_manifest_record(
    text: str, old_hash: str, new_hash: str, recorded_path: str
) -> str:
    old_line = f"{old_hash}  {recorded_path}"
    new_line = f"{new_hash}  {recorded_path}"
    require(
        text.splitlines().count(old_line) == 1,
        f"manifest lacks unique baseline record for {recorded_path}",
    )
    require(
        new_line not in text.splitlines(),
        f"manifest already has final record for {recorded_path}",
    )
    return text.replace(old_line, new_line)


def parse_manifest(text: str, name: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for number, line in enumerate(text.splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"malformed manifest row {name}:{number}")
        value, recorded = match.groups()
        require(recorded not in records, f"duplicate manifest path {recorded} in {name}")
        require(not Path(recorded).is_absolute(), f"absolute manifest path {recorded} in {name}")
        records[recorded] = value
    require(records, f"empty manifest {name}")
    return records


def verify_manifest_files(root: Path, relative: Path, text: str) -> int:
    manifest_directory = (root / relative).parent
    records = parse_manifest(text, relative)
    resolved_root = root.resolve()
    symlinks = sorted(path for path in resolved_root.rglob("*") if path.is_symlink())
    require(
        not symlinks,
        "symlinks are forbidden in the publication source: "
        + ", ".join(str(path.relative_to(resolved_root)) for path in symlinks),
    )
    seen_sources: set[Path] = set()
    for recorded, expected in records.items():
        source = (manifest_directory / recorded).resolve()
        require(
            source == resolved_root or resolved_root in source.parents,
            f"manifest path escapes source root: {relative}:{recorded}",
        )
        require(source.is_file(), f"missing manifest input {relative}:{recorded}")
        require(source not in seen_sources, f"duplicate manifest artifact {relative}:{recorded}")
        seen_sources.add(source)
        require(digest(source.read_bytes()) == expected, f"manifest mismatch {relative}:{recorded}")
    return len(records)


@dataclass(frozen=True)
class GraphAudit:
    results: int
    proofs: int
    labels: tuple[str, ...]
    edges: frozenset[tuple[str, str]]
    nontrivial_components: tuple[frozenset[str], ...]
    reachable: frozenset[str]


def strongly_connected_components(
    labels: tuple[str, ...], edges: frozenset[tuple[str, str]]
) -> tuple[frozenset[str], ...]:
    adjacency = {label: set() for label in labels}
    for source, target in edges:
        adjacency[source].add(target)

    next_index = 0
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    stack: list[str] = []
    on_stack: set[str] = set()
    components: list[frozenset[str]] = []

    def visit(vertex: str) -> None:
        nonlocal next_index
        indices[vertex] = next_index
        lowlinks[vertex] = next_index
        next_index += 1
        stack.append(vertex)
        on_stack.add(vertex)

        for target in sorted(adjacency[vertex]):
            if target not in indices:
                visit(target)
                lowlinks[vertex] = min(lowlinks[vertex], lowlinks[target])
            elif target in on_stack:
                lowlinks[vertex] = min(lowlinks[vertex], indices[target])

        if lowlinks[vertex] != indices[vertex]:
            return
        component: set[str] = set()
        while True:
            target = stack.pop()
            on_stack.remove(target)
            component.add(target)
            if target == vertex:
                break
        if len(component) > 1:
            components.append(frozenset(component))

    for label in labels:
        if label not in indices:
            visit(label)
    return tuple(sorted(components, key=lambda item: tuple(sorted(item))))


def reachable_labels(
    labels: tuple[str, ...], edges: frozenset[tuple[str, str]], start: str
) -> frozenset[str]:
    require(start in labels, f"missing reachability root {start}")
    adjacency = {label: set() for label in labels}
    for source, target in edges:
        adjacency[source].add(target)
    reached: set[str] = set()
    frontier = [start]
    while frontier:
        label = frontier.pop()
        if label in reached:
            continue
        reached.add(label)
        frontier.extend(sorted(adjacency[label] - reached))
    return frozenset(reached)


def audit_graph(text: str, expected_results: int = 850) -> GraphAudit:
    results = list(RESULT_PATTERN.finditer(text))
    require(
        len(results) == expected_results,
        f"expected {expected_results} numbered results, found {len(results)}",
    )

    main_proofs = list(MAIN_PROOF_PATTERN.finditer(text))
    require(len(main_proofs) == 1, f"expected one deferred main proof, found {len(main_proofs)}")
    main_proof = main_proofs[0]

    labels: list[str] = []
    proofs: list[str] = []
    proof_spans: list[tuple[int, int]] = []
    for index, result in enumerate(results):
        result_labels = re.findall(r"\\label\{([^}]+)\}", result.group(0))
        require(
            len(result_labels) == 1,
            f"result {index + 1} has {len(result_labels)} statement labels",
        )
        label = result_labels[0]
        require(label not in labels, f"duplicate result label {label}")
        labels.append(label)

        if index == 0:
            require(label == "thm:main", f"first result is {label}, not thm:main")
            proofs.append(main_proof.group(0))
            proof_spans.append(main_proof.span())
            continue
        following = text[result.end() :]
        proof = IMMEDIATE_PROOF_PATTERN.match(following)
        require(proof is not None, f"result {label} lacks its immediate proof")
        proofs.append(proof.group(1))
        proof_spans.append(
            (result.end() + proof.start(1), result.end() + proof.end(1))
        )

    global_proofs = list(re.finditer(r"\\begin\{proof\}(?:\[[^]]*\])?", text))
    require(
        len(global_proofs) == expected_results,
        f"expected {expected_results} proofs, found {len(global_proofs)}",
    )
    require(len(set(proof_spans)) == expected_results, "a proof was paired more than once")
    require(
        {start for start, _end in proof_spans} == {proof.start() for proof in global_proofs},
        "the result/proof pairing omits or invents a proof block",
    )

    owner = {label: index for index, label in enumerate(labels)}
    edge_indices: set[tuple[int, int]] = set()
    for source, proof in enumerate(proofs):
        for group in RESULT_REFERENCE_PATTERN.findall(proof):
            for target in map(str.strip, group.split(",")):
                if target.startswith(RESULT_LABEL_PREFIXES):
                    require(target in owner, f"undefined result reference {target}")
                if target in owner and owner[target] != source:
                    edge_indices.add((source, owner[target]))

    label_tuple = tuple(labels)
    edges = frozenset((label_tuple[source], label_tuple[target]) for source, target in edge_indices)
    components = strongly_connected_components(label_tuple, edges)
    reachable = reachable_labels(label_tuple, edges, "thm:main")
    return GraphAudit(
        results=len(results),
        proofs=len(proofs),
        labels=label_tuple,
        edges=edges,
        nontrivial_components=components,
        reachable=reachable,
    )


def verify_baseline(data: bytes) -> tuple[str, GraphAudit]:
    require(digest(data) == BASELINE_SHA256, "execution paper is not the replayed baseline")
    text = data.decode("utf-8")
    require(text.count(OLD_PASSAGE) == 1, "baseline does not contain exactly one old passage")
    require(NEW_PASSAGE not in text, "baseline already contains the publication rebind")
    graph = audit_graph(text)
    expected_component = frozenset(
        {
            "prop:post29",
            "prop:post29-bc-half-bridge",
            "prop:post29-direct",
            "thm:classical",
        }
    )
    require(len(graph.edges) == 2099, f"baseline edge count changed: {len(graph.edges)}")
    require(
        len(graph.reachable) == 348,
        f"baseline reachable count changed: {len(graph.reachable)}",
    )
    require(
        graph.nontrivial_components == (expected_component,),
        f"baseline dependency component changed: {graph.nontrivial_components}",
    )
    return text, graph


def verify_candidate(
    baseline_text: str, baseline_graph: GraphAudit, candidate_data: bytes
) -> GraphAudit:
    expected_text = baseline_text.replace(OLD_PASSAGE, NEW_PASSAGE)
    require(
        candidate_data == expected_text.encode("utf-8"),
        "publication paper has changes beyond the one authorized prose rebind",
    )
    require(digest(candidate_data) == FINAL_SHA256, "publication paper final digest changed")
    candidate_text = candidate_data.decode("utf-8")
    require(OLD_PASSAGE not in candidate_text, "old cyclic passage remains")
    require(candidate_text.count(NEW_PASSAGE) == 1, "new rebind passage is not unique")
    graph = audit_graph(candidate_text)
    require(graph.labels == baseline_graph.labels, "numbered-result label order changed")
    require(len(graph.edges) == 2098, f"final edge count changed: {len(graph.edges)}")
    removed = baseline_graph.edges - graph.edges
    added = graph.edges - baseline_graph.edges
    require(
        removed == {("prop:post29", "thm:classical")},
        f"wrong dependency edge removed: {sorted(removed)}",
    )
    require(not added, f"new dependency edges introduced: {sorted(added)}")
    require(
        graph.reachable == baseline_graph.reachable,
        "publication rebind changes the main-theorem reachable result set",
    )
    require(
        not graph.nontrivial_components,
        f"final proof graph remains cyclic: {graph.nontrivial_components}",
    )
    return graph


def expected_final_manifests(
    baseline_replay: str, baseline_full: str
) -> tuple[str, str]:
    final_replay = replace_manifest_record(
        baseline_replay, BASELINE_SHA256, FINAL_SHA256, "paper.tex"
    )
    require(
        digest(final_replay.encode("utf-8")) == FINAL_REPLAY_MANIFEST_SHA256,
        "computed final replay manifest digest changed",
    )
    final_full = replace_manifest_record(
        baseline_full, BASELINE_SHA256, FINAL_SHA256, "../../paper.tex"
    )
    final_full = replace_manifest_record(
        final_full,
        BASELINE_REPLAY_MANIFEST_SHA256,
        FINAL_REPLAY_MANIFEST_SHA256,
        "../../replay_sources.sha256",
    )
    require(
        digest(final_full.encode("utf-8")) == FINAL_FULL_MANIFEST_SHA256,
        "computed final Part III manifest digest changed",
    )
    return final_replay, final_full


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execution-root",
        type=Path,
        required=True,
        help="executed ginibre_q3 source root, normally the archived snapshot",
    )
    parser.add_argument("--final-root", type=Path, help="publication ginibre_q3 source root")
    parser.add_argument(
        "--preflight",
        action="store_true",
        help="audit the exact in-memory candidate without changing paper.tex",
    )
    args = parser.parse_args()
    require(
        args.preflight != (args.final_root is not None),
        "choose exactly one of --preflight and --final-root",
    )

    execution_root = args.execution_root.resolve()
    baseline_data = (execution_root / "paper.tex").read_bytes()
    baseline_text, baseline_graph = verify_baseline(baseline_data)
    baseline_replay = (execution_root / REPLAY_MANIFEST).read_text(encoding="utf-8")
    baseline_full = (execution_root / FULL_MANIFEST).read_text(encoding="utf-8")
    require(
        digest(baseline_replay.encode("utf-8")) == BASELINE_REPLAY_MANIFEST_SHA256,
        "execution replay manifest identity changed",
    )
    require(
        digest(baseline_full.encode("utf-8")) == BASELINE_FULL_MANIFEST_SHA256,
        "execution Part III manifest identity changed",
    )
    replay_entries = verify_manifest_files(execution_root, REPLAY_MANIFEST, baseline_replay)
    full_entries = verify_manifest_files(execution_root, FULL_MANIFEST, baseline_full)
    require(replay_entries == 85, f"execution replay manifest has {replay_entries} entries")
    require(full_entries == 36, f"execution Part III manifest has {full_entries} entries")
    final_replay, final_full = expected_final_manifests(baseline_replay, baseline_full)

    intermediate_data = baseline_text.replace(OLD_PASSAGE, NEW_PASSAGE).encode("utf-8")
    intermediate_graph = verify_candidate(baseline_text, baseline_graph, intermediate_data)

    if args.preflight:
        candidate_data = intermediate_data
        final_graph = intermediate_graph
        publication_replay_hash = FINAL_REPLAY_MANIFEST_SHA256
        publication_full_hash = FINAL_FULL_MANIFEST_SHA256
        publication_entries = (85, 36)
        mode = "preflight"
    else:
        require(args.final_root is not None, "missing final source root")
        final_root = args.final_root.resolve()
        candidate_data = (final_root / "paper.tex").read_bytes()
        candidate_text = candidate_data.decode("utf-8")
        require(OLD_PASSAGE not in candidate_text, "old cyclic passage remains in final revision")
        require(candidate_text.count(NEW_PASSAGE) == 1, "new rebind passage is not unique")
        final_graph = audit_graph(candidate_text, 852)
        added_labels = {
            "prop:bc-active-correction-prefix-contract",
            "prop:post29-bc-local-half-bridge",
        }
        require(
            tuple(label for label in final_graph.labels if label not in added_labels)
            == intermediate_graph.labels,
            "final revision changes historical result labels or their order",
        )
        require(
            set(final_graph.labels) - set(intermediate_graph.labels) == added_labels,
            "final revision adds an unauthorized numbered result",
        )
        expected_added_edges = {
            (
                "prop:bc-active-correction-prefix-contract",
                "prop:bc-pieri-correction-sign",
            ),
            ("prop:post29-bc-local-half-bridge", "lem:moment-formula"),
            ("prop:post29-bc-residual-closure", "lem:moment-formula"),
            (
                "prop:post29-bc-residual-closure",
                "prop:bc-active-correction-prefix-contract",
            ),
            ("prop:post29-bc-residual-closure", "prop:post29"),
            ("prop:post29-bc-residual-closure", "prop:post29-b-tilted-mgf"),
            (
                "prop:post29-bc-residual-closure",
                "prop:post29-bc-local-half-bridge",
            ),
            ("prop:post29-bc-residual-closure", "prop:post29-bc-interval-direct"),
            ("prop:post29-bc-residual-closure", "prop:post29-bc-layered-mgf"),
            ("prop:post29-bc-residual-closure", "prop:post29-c-tilted-mgf"),
        }
        require(
            final_graph.edges - intermediate_graph.edges == expected_added_edges,
            "final revision adds unauthorized proof dependencies",
        )
        removed_edges = intermediate_graph.edges - final_graph.edges
        require(
            len(removed_edges) == 67
            and all(source == "thm:classical" for source, _target in removed_edges),
            "final revision removes dependencies outside the redundant classical list",
        )
        removed_targets = {target for _source, target in removed_edges}
        removed_corrections = {
            target
            for target in removed_targets
            if target.startswith(("prop:bc-b-", "prop:bc-c-"))
            and target.endswith("lower-correction")
        }
        require(
            len(removed_corrections) == 60
            and removed_targets - removed_corrections
            == {
                "prop:bc-pieri-correction-sign",
                "prop:post29",
                "prop:post29-b-tilted-mgf",
                "prop:post29-bc-half-bridge",
                "prop:post29-bc-interval-direct",
                "prop:post29-bc-layered-mgf",
                "prop:post29-c-tilted-mgf",
            },
            "final revision removes the wrong classical direct suppliers",
        )
        require(
            r"\contractref{prop:post29-bc-local-half-bridge}" in candidate_text
            and "post_m29_bc_interval_bridge_frontier_gmp.cpp" in candidate_text
            and r"\mathcal L_m>0" in candidate_text,
            "final revision omits the explicit half-stable bridge predicate",
        )
        require(
            len(final_graph.reachable) == 90
            and added_labels <= final_graph.reachable,
            "final revision has the wrong compact proof-spine reachability",
        )
        require(not final_graph.nontrivial_components, "final revision introduces a proof cycle")
        actual_replay = (final_root / REPLAY_MANIFEST).read_text(encoding="utf-8")
        actual_full = (final_root / FULL_MANIFEST).read_text(encoding="utf-8")
        replay_records = parse_manifest(actual_replay, REPLAY_MANIFEST)
        full_records = parse_manifest(actual_full, FULL_MANIFEST)
        require(
            replay_records.get("paper.tex") == digest(candidate_data),
            "final replay manifest does not bind revised paper.tex",
        )
        require(
            full_records.get("../../paper.tex") == digest(candidate_data),
            "final Part III manifest does not bind revised paper.tex",
        )
        live_replay_entries = verify_manifest_files(final_root, REPLAY_MANIFEST, actual_replay)
        live_full_entries = verify_manifest_files(final_root, FULL_MANIFEST, actual_full)
        require(live_replay_entries >= replay_entries, "live replay manifest lost inputs")
        require(live_full_entries >= full_entries, "live Part III manifest lost inputs")
        publication_replay_hash = digest(actual_replay.encode("utf-8"))
        publication_full_hash = digest(actual_full.encode("utf-8"))
        publication_entries = (live_replay_entries, live_full_entries)
        mode = "final-revision"

    print(
        f"FULL_Q3_PAPER_REBIND mode={mode} results={final_graph.results} "
        f"proofs={final_graph.proofs} edges={len(baseline_graph.edges)}->{len(final_graph.edges)} "
        "historical_removed_edge=prop:post29->thm:classical "
        "current_redundant_edges_removed=67 components=1->0 "
        f"reachable={len(baseline_graph.reachable)}->{len(final_graph.reachable)} "
        f"baseline_sha256={BASELINE_SHA256} final_sha256={digest(candidate_data)} "
        f"replay_manifest={BASELINE_REPLAY_MANIFEST_SHA256}->{publication_replay_hash} "
        f"full_manifest={BASELINE_FULL_MANIFEST_SHA256}->{publication_full_hash} "
        f"manifest_entries={replay_entries}+{full_entries}"
        f"->{publication_entries[0]}+{publication_entries[1]} "
        f"verifier_sha256={digest(Path(__file__).resolve().read_bytes())}"
    )
    print("FULL_Q3_PAPER_REBIND VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, UnicodeError, RebindFailure) as error:
        print(f"FULL_Q3_PAPER_REBIND VERIFICATION: FAIL ({error})")
        raise SystemExit(1)
