#!/usr/bin/env python3
"""Run the complete short, non-arithmetic publication audit suite."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
AUDITS = (
    "verify_full_q3_document.py",
    "verify_parts_i_ii_contract.py",
    "referee_audits/verify_full_q3_moment_formula_independent.py",
    "referee_audits/verify_full_q3_cone_promotion_independent.py",
    "referee_audits/verify_full_q3_classification_coverage_independent.py",
    "referee_audits/verify_full_q3_so3_obstruction_independent.py",
    "referee_audits/verify_full_q3_dependency_graph_independent.py",
    "referee_audits/verify_full_q3_two_minus_interface_independent.py",
    "referee_audits/verify_full_q3_frontier_formulas_independent.py",
)


def main() -> int:
    for relative in AUDITS:
        source = ROOT / relative
        if not source.is_file():
            raise SystemExit(f"PUBLICATION_SHORT_AUDITS FAILURE: missing {relative}")
        completed = subprocess.run(
            [sys.executable, str(source)],
            cwd=ROOT.parent,
            check=False,
        )
        if completed.returncode != 0:
            raise SystemExit(
                f"PUBLICATION_SHORT_AUDITS FAILURE: {relative} exited "
                f"{completed.returncode}"
            )
    print(f"PUBLICATION_SHORT_AUDITS count={len(AUDITS)}")
    print("PUBLICATION_SHORT_AUDITS VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
