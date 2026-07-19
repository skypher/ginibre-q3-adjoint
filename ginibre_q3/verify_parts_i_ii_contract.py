#!/usr/bin/env python3
"""Verify the compact Parts I--II residual B/C proof interface."""

from __future__ import annotations

from pathlib import Path
import re

from clean_room_replay import ReplayFailure, audit_bc_caller_ranges


HALF_BRIDGE = "prop:post29-bc-half-bridge"
RESIDUAL = "prop:post29-bc-residual-closure"
TAIL_SUPPLIERS = (
    "prop:post29",
    "prop:post29-bc-layered-mgf",
    "prop:post29-b-tilted-mgf",
    "prop:post29-c-tilted-mgf",
    "prop:post29-bc-interval-direct",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReplayFailure(message)


def audit_explicit_half_bridge(root: Path) -> None:
    """Require the formerly implicit supplement dependency to stay explicit."""

    paper = (root / "paper.tex").read_text(encoding="utf-8")
    spine = (root / "PUBLICATION_PROOF_SPINE.md").read_text(encoding="utf-8")
    extension = (root / "full_q3_extension.tex").read_text(encoding="utf-8")

    require(
        rf"\label{{{HALF_BRIDGE}}}" in paper,
        "the numbered half-stable bridge reduction is absent",
    )
    residual_match = re.search(
        rf"\\label\{{{RESIDUAL}\}}.*?\\begin\{{proof\}}(.*?)\\end\{{proof\}}",
        paper,
        re.DOTALL,
    )
    require(residual_match is not None, "cannot parse the residual B/C proof")
    residual_proof = residual_match.group(1)
    require(
        rf"\contractref{{{HALF_BRIDGE}}}" in residual_proof,
        "residual B/C proof does not explicitly import the half-stable bridge",
    )
    require(
        "post_m29_bc_interval_bridge_frontier_gmp.cpp" in residual_proof
        and r"\mathcal L_m" in residual_proof
        and r"D_G(2m+1)\ge\mathcal L_m" in residual_proof,
        "residual B/C proof does not state the exact bridge predicate and source",
    )
    for label in TAIL_SUPPLIERS:
        require(
            label in residual_proof,
            f"residual B/C proof does not name its supplier {label}",
        )
    require(
        HALF_BRIDGE in spine
        and "explicit incoming node" in spine,
        "publication proof spine omits the half-stable bridge incoming edge",
    )
    require(
        "half-stable bridge reduction" in extension,
        "Part III submission boundary omits the half-stable bridge dependency",
    )


def main() -> int:
    root = Path(__file__).resolve().parent
    try:
        audit_explicit_half_bridge(root)
        active, targets, arithmetic_rows = audit_bc_caller_ranges(root)
    except ReplayFailure as error:
        raise SystemExit(f"PARTS_I_II_CONTRACT FAILURE: {error}") from error
    if (active, targets, arithmetic_rows) != (58, 2834, 402):
        raise SystemExit(
            "PARTS_I_II_CONTRACT FAILURE: unexpected summary "
            f"active={active} targets={targets} arithmetic_rows={arithmetic_rows}"
        )
    print(
        "PARTS_I_II_CONTRACT "
        f"active_corrections={active} formerly_uncovered_targets={targets} "
        "max_consumed_offset=27 checked_overlap_offset=28 "
        f"arithmetic_rows={arithmetic_rows} half_bridge_edge=PASS "
        f"named_suppliers={len(TAIL_SUPPLIERS)}"
    )
    print("PARTS_I_II_CONTRACT VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
