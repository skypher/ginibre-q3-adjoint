#!/usr/bin/env python3
"""Verify the compact Parts I--II residual B/C proof interface."""

from __future__ import annotations

from pathlib import Path

from clean_room_replay import ReplayFailure, audit_bc_caller_ranges


def main() -> int:
    root = Path(__file__).resolve().parent
    try:
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
        f"arithmetic_rows={arithmetic_rows}"
    )
    print("PARTS_I_II_CONTRACT VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
