#!/usr/bin/env python3
"""Parallel exact replay of the long H25--H29 B/C frontier certificates.

Python only orchestrates independent C++/GMP binaries.  Every arithmetic task
records its source identity, exact scope marker, terminal failure count, and
exit status in a separate immutable log.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class Task:
    name: str
    binary: str
    source: str
    environment: tuple[tuple[str, str], ...]
    required: tuple[str, ...]


ORDINALS = {
    "twentyfifth": {
        "shards": (
            (15, 17), (18, 20), (21, 23), (24, 26), (27, 29),
            (30, 32), (33, 35), (36, 38), (39, 41), (42, 43),
        ),
        "c_fpf": 10,
        "c_arithmetic": 33,
    },
    "twentysixth": {
        "shards": (
            (15, 17), (18, 20), (21, 23), (24, 26), (27, 29),
            (30, 32), (33, 35), (36, 38), (39, 41), (42, 43),
            (44, 44), (45, 45),
        ),
        "c_fpf": 11,
        "c_arithmetic": 34,
    },
    "twentyseventh": {
        "shards": (
            (15, 17), (18, 20), (21, 23), (24, 26), (27, 29),
            (30, 32), *((rank, rank) for rank in range(33, 48)),
        ),
        "c_fpf": 12,
        "c_arithmetic": 36,
    },
    "twentyeighth": {
        "shards": (
            (15, 17), (18, 20), (21, 23),
            *((rank, rank) for rank in range(24, 42)),
        ),
        "c_fpf": 13,
        "c_arithmetic": 37,
    },
}


FINAL_TASKS = (
    (
        "twentyeighth_b_absorption",
        "post_m29_bc_twentyeighth_b_absorption_gmp",
        "BC_TWENTYEIGHTH_B_ABSORPTION_GMP h=28 lower_boxes=56 "
        "one_lower_choices=112 two_lower_choices=1596",
    ),
    (
        "twentyninth_c_frontier",
        "post_m29_bc_twentyninth_c_frontier_gmp",
        "BC_TWENTYNINTH_C_FRONTIER_GMP c_fpf_cases=14 c_arithmetic_cases=39",
    ),
    (
        "twentyninth_b_absorption",
        "post_m29_bc_twentyninth_b_absorption_gmp",
        "BC_TWENTYNINTH_B_ABSORPTION_GMP h=29 lower_boxes=58 "
        "one_lower_choices=116 two_lower_choices=1711",
    ),
)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def make_tasks(threads: int) -> list[Task]:
    tasks: list[Task] = []
    for ordinal, spec in ORDINALS.items():
        shards = spec["shards"]
        if not isinstance(shards, tuple):
            raise RuntimeError(f"invalid shard ledger type for {ordinal}")
        covered = [rank for low, high in shards for rank in range(low, high + 1)]
        if covered != list(range(15, max(covered) + 1)):
            raise RuntimeError(f"noncontiguous shard ledger for {ordinal}")
        binary = f"post_m29_bc_{ordinal}_frontier_gmp"
        source = f"{binary}.cpp"
        header = f"BC_{ordinal.upper()}_FRONTIER_GMP"
        c_fpf = int(spec["c_fpf"])
        c_arithmetic = int(spec["c_arithmetic"])
        for low, high in shards:
            cases = high - low + 1
            required = (
                f"{header} b_cases={cases} c_fpf_cases={c_fpf} "
                f"c_arithmetic_cases={c_arithmetic} run_b=1 run_c=0 "
                f"b_q_min={low} b_q_max={high}",
                "SUMMARY failures=0",
            )
            tasks.append(
                Task(
                    name=f"{ordinal}_b_{low}_{high}",
                    binary=binary,
                    source=source,
                    environment=(
                        ("RUN_B", "1"),
                        ("RUN_C", "0"),
                        ("B_Q_MIN", str(low)),
                        ("B_Q_MAX", str(high)),
                        ("OMP_NUM_THREADS", str(threads)),
                    ),
                    required=required,
                )
            )
        tasks.append(
            Task(
                name=f"{ordinal}_c",
                binary=binary,
                source=source,
                environment=(
                    ("RUN_B", "0"),
                    ("RUN_C", "1"),
                    ("OMP_NUM_THREADS", str(threads)),
                ),
                required=(
                    f"{header} b_cases=0 c_fpf_cases={c_fpf} "
                    f"c_arithmetic_cases={c_arithmetic} run_b=0 run_c=1",
                    "SUMMARY failures=0",
                ),
            )
        )
    for name, binary, scope in FINAL_TASKS:
        tasks.append(
            Task(
                name=name,
                binary=binary,
                source=f"{binary}.cpp",
                environment=(("OMP_NUM_THREADS", str(threads)),),
                required=(scope, "SUMMARY failures=0"),
            )
        )
    if len(tasks) != 71 or len({task.name for task in tasks}) != 71:
        raise RuntimeError(f"fleet task ledger changed: tasks={len(tasks)}")
    return tasks


def run_task(
    task: Task,
    character_dir: Path,
    output_dir: Path,
    timeout: int,
    host_label: str,
    run_id: str,
) -> tuple[str, float, str]:
    binary = character_dir / task.binary
    source = character_dir / task.source
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise RuntimeError(f"missing executable for {task.name}: {binary}")
    if not source.is_file():
        raise RuntimeError(f"missing source for {task.name}: {source}")
    source_hash = digest(source)
    binary_hash = digest(binary)
    log = output_dir / f"{task.name}.log"
    environment = os.environ.copy()
    environment.update({"LC_ALL": "C", "PYTHONHASHSEED": "0"})
    environment.update(dict(task.environment))
    started = time.monotonic()
    with log.open("w") as output:
        output.write(f"__RUN_ID={run_id}\n")
        output.write(f"__HOST_LABEL={host_label}\n")
        output.write(f"__TASK={task.name}\n")
        output.write(f"__START_UTC={utc_now()}\n")
        output.write(f"__SOURCE_SHA256={source_hash}  {task.source}\n")
        output.write(f"__BINARY_SHA256={binary_hash}  {task.binary}\n")
        output.write("__ENV=" + " ".join(f"{k}={v}" for k, v in task.environment) + "\n")
        output.write(f"__COMMAND=./{task.binary}\n")
        output.flush()
        try:
            completed = subprocess.run(
                [str(binary)],
                cwd=character_dir,
                env=environment,
                stdout=output,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            status = completed.returncode
        except subprocess.TimeoutExpired:
            status = 124
            output.write(f"__TIMEOUT_SECONDS={timeout}\n")
        elapsed = time.monotonic() - started
        output.write(f"__EXIT_STATUS={status}\n")
        output.write(f"__ELAPSED_SECONDS={elapsed:.6f}\n")
        output.write(f"__END_UTC={utc_now()}\n")
    text = log.read_text(errors="strict")
    missing = [fragment for fragment in task.required if fragment not in text]
    if status != 0 or missing:
        raise RuntimeError(
            f"{task.name} failed: status={status}, missing={missing}, log={log}"
        )
    return task.name, elapsed, digest(log)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--character-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=28800)
    parser.add_argument("--host-label", default=os.uname().nodename)
    parser.add_argument("--run-id", default="fullq3cleanfleet0002")
    args = parser.parse_args()
    if min(args.workers, args.threads, args.timeout) < 1:
        parser.error("workers, threads, and timeout must be positive")
    if not args.run_id or any(
        character not in "abcdefghijklmnopqrstuvwxyz0123456789_-"
        for character in args.run_id
    ):
        parser.error("run-id must contain only lowercase ASCII letters, digits, _ or -")
    character_dir = args.character_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    tasks = make_tasks(args.threads)
    script_hash = digest(Path(__file__).resolve())
    print(
        f"FULL_Q3_FLEET start run_id={args.run_id} tasks={len(tasks)} workers={args.workers} "
        f"threads_per_task={args.threads} script_sha256={script_hash}",
        flush=True,
    )
    results: list[tuple[str, float, str]] = []
    failure: BaseException | None = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                run_task,
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
                    f"FULL_Q3_FLEET task={result[0]} status=PASS "
                    f"elapsed={result[1]:.2f}s sha256={result[2]}",
                    flush=True,
                )
            except BaseException as error:  # fail closed after stopping queued work
                failure = error
                print(f"FULL_Q3_FLEET task={task.name} status=FAIL error={error}", flush=True)
                for pending in futures:
                    pending.cancel()
                break
    if failure is not None:
        print(f"FULL_Q3_FLEET VERIFICATION: FAIL ({failure})", flush=True)
        return 1
    results.sort()
    if len(results) != 71 or {name for name, _, _ in results} != {task.name for task in tasks}:
        print(f"FULL_Q3_FLEET VERIFICATION: FAIL (coverage={len(results)}/71)", flush=True)
        return 1
    manifest = output_dir / f"{args.run_id}_logs.sha256"
    with manifest.open("w") as output:
        for name, _elapsed, log_hash in results:
            output.write(f"{log_hash}  {name}.log\n")
    print("FULL_Q3_FLEET tasks_passed=71 tasks_expected=71 failures=0", flush=True)
    print("FULL_Q3_FLEET VERIFICATION: ALL PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
