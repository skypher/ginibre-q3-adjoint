#!/usr/bin/env python3
"""Independent all-reference-macro dependency audit for Part III."""

from __future__ import annotations

import re
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "full_q3_extension.tex"
RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")
RESULT_PATTERN = re.compile(
    r"\\begin\{(" + "|".join(RESULT_ENVS) + r")\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
REFERENCE_PATTERN = re.compile(
    r"\\(eqref|ref|cref|Cref|autoref)\{([^}]+)\}"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FULL_Q3_DEPENDENCY_AUDIT FAILURE: {message}")


def references(text: str) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    for macro, group in REFERENCE_PATTERN.findall(text):
        output.extend((macro, label.strip()) for label in group.split(","))
    return output


def nontrivial_sccs(
    vertex_count: int, edges: set[tuple[int, int]]
) -> list[list[int]]:
    adjacency = {vertex: [] for vertex in range(vertex_count)}
    for source, target in edges:
        if source != target:
            adjacency[source].append(target)

    next_index = 0
    indices: dict[int, int] = {}
    lowlinks: dict[int, int] = {}
    stack: list[int] = []
    on_stack: set[int] = set()
    components: list[list[int]] = []

    def visit(vertex: int) -> None:
        nonlocal next_index
        indices[vertex] = next_index
        lowlinks[vertex] = next_index
        next_index += 1
        stack.append(vertex)
        on_stack.add(vertex)
        for target in adjacency[vertex]:
            if target not in indices:
                visit(target)
                lowlinks[vertex] = min(lowlinks[vertex], lowlinks[target])
            elif target in on_stack:
                lowlinks[vertex] = min(lowlinks[vertex], indices[target])
        if lowlinks[vertex] != indices[vertex]:
            return
        component: list[int] = []
        while True:
            target = stack.pop()
            on_stack.remove(target)
            component.append(target)
            if target == vertex:
                break
        if len(component) > 1:
            components.append(sorted(component))

    for vertex in range(vertex_count):
        if vertex not in indices:
            visit(vertex)
    return components


def main() -> int:
    text = SOURCE.read_text(encoding="utf-8")
    result_matches = list(RESULT_PATTERN.finditer(text))
    proof_matches = list(
        re.finditer(r"\\begin\{proof\}.*?\\end\{proof\}", text, re.DOTALL)
    )
    require(len(result_matches) == 51, f"results={len(result_matches)}")
    require(len(proof_matches) == 51, f"proofs={len(proof_matches)}")

    records: list[tuple[re.Match[str], str, list[str]]] = []
    for number, result in enumerate(result_matches, 1):
        proof = re.match(
            r"\s*(\\begin\{proof\}(.*?)\\end\{proof\})",
            text[result.end() :],
            re.DOTALL,
        )
        require(proof is not None, f"result {number} lacks an immediate proof")
        labels = re.findall(r"\\label\{([^}]+)\}", result.group(2))
        require(labels, f"result {number} has no statement label")
        records.append((result, proof.group(1), labels))

    owner: dict[str, int] = {}
    for index, (_result, _proof, labels) in enumerate(records):
        for label in labels:
            require(label not in owner, f"duplicate owned label {label}")
            owner[label] = index
    require(len(owner) == 68, f"owned labels={len(owner)}")

    all_labels = re.findall(r"\\label\{([^}]+)\}", text)
    require(len(all_labels) == len(set(all_labels)), "duplicate document label")
    all_references = references(text)
    missing = sorted({label for _macro, label in all_references} - set(all_labels))
    require(not missing, f"undefined references={missing}")
    require("\\begin{assumption}" not in text, "assumption environment present")
    require(
        not any(label.startswith("ass:") for label in all_labels),
        "assumption-labelled object present",
    )

    proof_edges: set[tuple[int, int]] = set()
    whole_edges: set[tuple[int, int]] = set()
    proof_macro_counts = {name: 0 for name in ("eqref", "ref", "cref", "Cref", "autoref")}
    for index, (result, proof, _labels) in enumerate(records):
        for macro, label in references(proof):
            proof_macro_counts[macro] += 1
            target = owner.get(label)
            if target is not None and target != index:
                proof_edges.add((index, target))
        for _macro, label in references(result.group(2) + proof):
            target = owner.get(label)
            if target is not None and target != index:
                whole_edges.add((index, target))

    require(len(proof_edges) == 96, f"proof edges={len(proof_edges)}")
    require(len(whole_edges) == 103, f"whole-result edges={len(whole_edges)}")
    proof_forward = sorted(edge for edge in proof_edges if edge[1] > edge[0])
    whole_forward = sorted(edge for edge in whole_edges if edge[1] > edge[0])
    require(not proof_forward, f"forward proof edges={proof_forward}")
    require(
        whole_forward == [(0, 3), (1, 2)],
        f"unexpected statement-level forward edges={whole_forward}",
    )
    require(not nontrivial_sccs(len(records), proof_edges), "proof graph cycle")
    require(not nontrivial_sccs(len(records), whole_edges), "whole-result graph cycle")

    root = owner.get("thm:full-adjoint-generated-q3")
    require(root is not None, "final theorem label absent")
    reachable: set[int] = set()
    pending = [root]
    while pending:
        vertex = pending.pop()
        if vertex in reachable:
            continue
        reachable.add(vertex)
        pending.extend(target for source, target in proof_edges if source == vertex)
    require(len(reachable) == 47, f"reachable results={len(reachable)}")
    unreachable_primary = {
        records[index][2][0]
        for index in range(len(records))
        if index not in reachable
    }
    require(
        unreachable_primary
        == {
            "lem:adjoint-cone-positive-definite",
            "prop:all-pd-counterexample",
            "cor:fixed-group-finite-residual",
            "lem:full-classical-negative-endpoints",
        },
        f"unexpected unreachable results={sorted(unreachable_primary)}",
    )

    # The pinned document checker currently scans ref/eqref.  Confirm that no
    # proof dependency escapes it through a cleveref-style macro, while this
    # supplementary audit remains future-proof if one is introduced later.
    require(proof_macro_counts["cref"] == 0, "cref proof reference introduced")
    require(proof_macro_counts["Cref"] == 0, "Cref proof reference introduced")
    require(proof_macro_counts["autoref"] == 0, "autoref proof reference introduced")

    print(
        "FULL_Q3_DEPENDENCY_AUDIT "
        f"results={len(records)} proofs={len(proof_matches)} "
        f"owned_labels={len(owner)} proof_edges={len(proof_edges)} "
        f"whole_edges={len(whole_edges)} proof_forward={len(proof_forward)} "
        "proof_scc=0 whole_scc=0 "
        f"final_reachable={len(reachable)} assumptions=0 "
        f"proof_ref={proof_macro_counts['ref']} "
        f"proof_eqref={proof_macro_counts['eqref']} "
        "proof_cref=0"
    )
    print("FULL_Q3_DEPENDENCY_AUDIT VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
