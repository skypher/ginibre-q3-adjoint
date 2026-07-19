#!/usr/bin/env python3
"""Run a source-pinned classical or exceptional post-frontier replay."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--section", choices=("classical", "exceptional"), required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads must be positive")

    root = args.root.resolve()
    sys.path.insert(0, str(root))
    import clean_room_replay as replay

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment.update(
        {
            "LC_ALL": "C",
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": str(args.threads),
            "SOURCE_DATE_EPOCH": "946684800",
        }
    )
    runner = replay.Runner(output, environment)

    class Collector:
        def __init__(self) -> None:
            self.names: list[str] = []

        def run(self, name: str, _command: list[str], _cwd: Path, **_kwargs: object) -> Path:
            self.names.append(name)
            return root.parent / "replay-logs" / f"{len(self.names):02d}_{name}.log"

    expected = Collector()
    if args.section == "classical":
        replay.classical_replays(expected, root)
        replay.classical_replays(runner, root)
        prefix = "FULL_Q3_CLASSICAL_POST"
        expected_count = 30
    else:
        replay.exceptional_replays(expected, root, args.threads)
        replay.exceptional_replays(runner, root, args.threads)
        prefix = "FULL_Q3_EXCEPTIONAL"
        expected_count = 16

    actual_names = [name for name, _elapsed in runner.results]
    if (
        len(expected.names) != expected_count
        or len(set(expected.names)) != expected_count
        or actual_names != expected.names
    ):
        raise replay.ReplayFailure(
            f"{args.section} stage coverage mismatch: "
            f"expected={expected.names}, actual={actual_names}"
        )
    print(f"__RUN_ID={args.run_id}")
    print(f"__HOST={os.uname().nodename}")
    print(f"__WRAPPER_SHA256={replay.digest(Path(__file__).resolve())}")
    print(f"__CLEAN_DRIVER_SHA256={replay.digest(root / 'clean_room_replay.py')}")
    print(f"__SOURCE_MANIFEST_SHA256={replay.digest(root / 'replay_sources.sha256')}")
    print(
        f"{prefix} tasks_passed={expected_count} "
        f"tasks_expected={expected_count} failures=0"
    )
    print(f"{prefix} VERIFICATION: ALL PASS")
    print("__EXIT_STATUS=0")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FULL_Q3_POST_REPLAY VERIFICATION: FAIL ({error})", file=sys.stderr)
        print("__EXIT_STATUS=1", file=sys.stderr)
        raise SystemExit(1)
