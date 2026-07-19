#!/usr/bin/env python3
"""Run one exact, disjoint A/B/C partition of the H25--H29 fleet ledger.

The arithmetic task definitions and the per-task fail-closed logger are loaded
from ``run_full_q3_frontier_fleet_queue.py``.  This wrapper only selects one of
three fixed partitions and schedules its independent C++/GMP processes.  The
three fixed selections have cardinalities 24, 7, and 40 and their union is the
original 71-task ledger.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import os
import sys
from pathlib import Path


PARTITIONS = ("machine_a_low", "machine_b_mid", "machine_c_heavy")


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def load_ledger(path: Path):
    spec = importlib.util.spec_from_file_location("full_q3_frontier_ledger", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import task ledger {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def partition_name(task_name: str) -> str:
    if task_name.startswith(("twentyfifth_", "twentysixth_")):
        return "machine_a_low"
    if task_name == "twentyseventh_c":
        return "machine_b_mid"
    if task_name.startswith("twentyseventh_b_"):
        high = int(task_name.rsplit("_", 1)[1])
        return "machine_b_mid" if high <= 32 else "machine_c_heavy"
    return "machine_c_heavy"


def partition_tasks(ledger, threads: int, partition: str):
    all_tasks = ledger.make_tasks(threads)
    grouped = {
        name: [task for task in all_tasks if partition_name(task.name) == name]
        for name in PARTITIONS
    }
    expected_counts = {
        "machine_a_low": 24,
        "machine_b_mid": 7,
        "machine_c_heavy": 40,
    }
    actual_counts = {name: len(tasks) for name, tasks in grouped.items()}
    if actual_counts != expected_counts:
        raise RuntimeError(
            f"partition cardinality changed: expected={expected_counts}, "
            f"actual={actual_counts}"
        )
    names = [task.name for tasks in grouped.values() for task in tasks]
    if len(names) != 71 or len(set(names)) != 71:
        raise RuntimeError("partition union is not the original 71-task ledger")
    return grouped[partition]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger-script", type=Path, required=True)
    parser.add_argument("--character-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--partition", choices=PARTITIONS, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--timeout", type=int, default=28800)
    parser.add_argument("--host-label", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if min(args.workers, args.threads, args.timeout) < 1:
        parser.error("workers, threads, and timeout must be positive")
    if not args.run_id or any(
        character not in "abcdefghijklmnopqrstuvwxyz0123456789_-"
        for character in args.run_id
    ):
        parser.error("run-id must contain only lowercase ASCII letters, digits, _ or -")

    ledger_path = args.ledger_script.resolve()
    if not ledger_path.is_file():
        parser.error(f"missing ledger script: {ledger_path}")
    ledger = load_ledger(ledger_path)
    tasks = partition_tasks(ledger, args.threads, args.partition)
    character_dir = args.character_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    wrapper_hash = digest(Path(__file__).resolve())
    ledger_hash = digest(ledger_path)
    print(
        "FULL_Q3_FLEET_PARTITION start "
        f"run_id={args.run_id} partition={args.partition} tasks={len(tasks)} "
        f"workers={args.workers} threads_per_task={args.threads} "
        f"wrapper_sha256={wrapper_hash} ledger_sha256={ledger_hash}",
        flush=True,
    )

    results: list[tuple[str, float, str]] = []
    failure: BaseException | None = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                ledger.run_task,
                task,
                character_dir,
                output_dir,
                args.timeout,
                args.host_label,
                args.run_id,
            ): task
            for task in tasks
        }
        for future in concurrent.futures.as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(
                    f"FULL_Q3_FLEET_PARTITION task={result[0]} status=PASS "
                    f"elapsed={result[1]:.2f}s sha256={result[2]}",
                    flush=True,
                )
            except BaseException as error:
                failure = error
                print(
                    f"FULL_Q3_FLEET_PARTITION task={task.name} status=FAIL "
                    f"error={error}",
                    flush=True,
                )
                for pending in futures:
                    pending.cancel()
                break

    if failure is not None:
        print(
            f"FULL_Q3_FLEET_PARTITION VERIFICATION: FAIL ({failure})",
            flush=True,
        )
        return 1
    results.sort()
    expected_names = {task.name for task in tasks}
    if len(results) != len(tasks) or {name for name, _, _ in results} != expected_names:
        print(
            "FULL_Q3_FLEET_PARTITION VERIFICATION: FAIL "
            f"(coverage={len(results)}/{len(tasks)})",
            flush=True,
        )
        return 1
    manifest = output_dir / f"{args.run_id}_logs.sha256"
    with manifest.open("w") as output:
        for name, _elapsed, log_hash in results:
            output.write(f"{log_hash}  {name}.log\n")
    print(
        f"FULL_Q3_FLEET_PARTITION partition={args.partition} "
        f"tasks_passed={len(tasks)} tasks_expected={len(tasks)} failures=0",
        flush=True,
    )
    print("FULL_Q3_FLEET_PARTITION VERIFICATION: ALL PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
