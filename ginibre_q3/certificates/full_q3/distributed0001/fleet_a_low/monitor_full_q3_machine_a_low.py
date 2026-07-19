#!/usr/bin/env python3
"""Freeze the accepted machine-A H25--H26 subset of a live 71-task queue.

The original machine-A run predates the disjoint partition wrapper and has all
71 tasks queued.  Its first 24 ledger entries are exactly H25--H26.  This
watcher authenticates those completed logs against the live top transcript,
stops the process group as soon as all 24 have passed, and writes the exact
subset manifest needed by the distributed verifier.  Partial later-task logs
are deliberately excluded.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import re
import signal
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


FLEET_SCRIPT_SHA256 = (
    "28b51756e053cb160e9cf5dde4f97a27ae416f084f7ef23ffedd02409098464a"
)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_ledger(path: Path):
    spec = importlib.util.spec_from_file_location("full_q3_machine_a_ledger", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import fleet ledger: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def process_group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    return True


def stop_process_group(pgid: int) -> str:
    """Freeze, terminate, and if necessary kill the authenticated run group."""

    if not process_group_exists(pgid):
        return "already_exited"
    os.killpg(pgid, signal.SIGSTOP)
    os.killpg(pgid, signal.SIGTERM)
    os.killpg(pgid, signal.SIGCONT)
    deadline = time.monotonic() + 10.0
    while process_group_exists(pgid) and time.monotonic() < deadline:
        time.sleep(0.2)
    if process_group_exists(pgid):
        os.killpg(pgid, signal.SIGKILL)
        return "sigkill_after_sigterm"
    return "sigterm"


def atomic_write(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--orchestrator-pid", type=int, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--poll-seconds", type=float, default=30.0)
    args = parser.parse_args()
    if args.orchestrator_pid < 2 or args.poll_seconds <= 0:
        parser.error("invalid PID or polling interval")

    run_root = args.run_root.resolve()
    ledger_path = run_root / "run_full_q3_frontier_fleet_queue.py"
    top_log = run_root / f"{args.run_id}_a.log"
    results = run_root / "results"
    receipt = run_root / "machine_a_low_freeze.log"
    manifest = results / "machine_a_low_logs.sha256"
    if digest(ledger_path) != FLEET_SCRIPT_SHA256:
        raise RuntimeError("machine-A staged ledger identity mismatch")
    ledger = load_ledger(ledger_path)
    tasks = ledger.make_tasks(8)
    selected = [
        task
        for task in tasks
        if task.name.startswith(("twentyfifth_", "twentysixth_"))
    ]
    selected_names = {task.name for task in selected}
    if len(tasks) != 71 or len(selected) != 24 or len(selected_names) != 24:
        raise RuntimeError("machine-A partition ledger changed")

    pid = args.orchestrator_pid
    command_path = Path(f"/proc/{pid}/cmdline")
    if not command_path.is_file():
        raise RuntimeError(f"orchestrator PID is not live: {pid}")
    command = command_path.read_bytes().replace(b"\0", b" ").decode(errors="strict")
    if (
        "run_full_q3_frontier_fleet_queue.py" not in command
        or args.run_id not in command
    ):
        raise RuntimeError(f"PID {pid} is not the authenticated machine-A run")
    pgid = os.getpgid(pid)
    if pgid != pid:
        raise RuntimeError(f"machine-A orchestrator is not its group leader: {pid}/{pgid}")

    pass_pattern = re.compile(
        r"^FULL_Q3_FLEET task=([a-z0-9_]+) status=PASS "
        r"elapsed=([0-9]+(?:\.[0-9]+)?)s sha256=([0-9a-f]{64})$",
        re.MULTILINE,
    )
    while True:
        top_text = top_log.read_text(errors="strict")
        if " status=FAIL " in top_text or "VERIFICATION: FAIL" in top_text:
            stopped = stop_process_group(pgid)
            atomic_write(
                receipt,
                f"__UTC={utc_now()}\n__STATUS=FAIL\n__STOP_METHOD={stopped}\n",
            )
            raise RuntimeError("machine-A queue reported a task failure")
        records = pass_pattern.findall(top_text)
        record_by_name: dict[str, tuple[str, str]] = {}
        for name, elapsed, value in records:
            if name in record_by_name:
                raise RuntimeError(f"duplicate machine-A PASS record: {name}")
            record_by_name[name] = (elapsed, value)
        passed = selected_names & set(record_by_name)
        if passed == selected_names:
            # Freeze first so neither the transcript nor any selected log can
            # change while the exact subset is authenticated.
            os.killpg(pgid, signal.SIGSTOP)
            try:
                top_text = top_log.read_text(errors="strict")
                records = pass_pattern.findall(top_text)
                record_by_name = {}
                for name, elapsed, value in records:
                    if name in record_by_name:
                        raise RuntimeError(f"duplicate machine-A PASS record: {name}")
                    record_by_name[name] = (elapsed, value)
                if not selected_names <= set(record_by_name):
                    raise RuntimeError("machine-A transcript changed during freeze")
                rows: list[tuple[str, str]] = []
                for task in selected:
                    log = results / f"{task.name}.log"
                    text = log.read_text(errors="strict")
                    lines = text.splitlines()
                    statuses = [
                        line for line in lines if line.startswith("__EXIT_STATUS=")
                    ]
                    if (
                        f"__RUN_ID={args.run_id}" not in lines
                        or f"__TASK={task.name}" not in lines
                        or "SUMMARY failures=0" not in lines
                        or statuses != ["__EXIT_STATUS=0", "__EXIT_STATUS=0"]
                        or len(lines) < 3
                        or lines[-3] != "__EXIT_STATUS=0"
                        or not lines[-2].startswith("__ELAPSED_SECONDS=")
                        or not lines[-1].startswith("__END_UTC=")
                    ):
                        raise RuntimeError(f"invalid completed machine-A log: {log}")
                    value = digest(log)
                    if record_by_name[task.name][1] != value:
                        raise RuntimeError(f"machine-A transcript hash mismatch: {log}")
                    rows.append((task.name, value))
            except Exception as error:
                stopped = stop_process_group(pgid)
                atomic_write(
                    receipt,
                    "\n".join(
                        (
                            f"__UTC={utc_now()}",
                            "__STATUS=FAIL",
                            f"__STOP_METHOD={stopped}",
                            f"__ERROR={error}",
                            "",
                        )
                    ),
                )
                raise
            stopped = stop_process_group(pgid)
            rows.sort()
            atomic_write(
                manifest,
                "".join(f"{value}  {name}.log\n" for name, value in rows),
            )
            atomic_write(
                receipt,
                "\n".join(
                    (
                        f"__UTC={utc_now()}",
                        "__STATUS=PASS",
                        f"__RUN_ID={args.run_id}",
                        f"__ORCHESTRATOR_PID={pid}",
                        f"__PROCESS_GROUP={pgid}",
                        f"__STOP_METHOD={stopped}",
                        f"__LEDGER_SHA256={FLEET_SCRIPT_SHA256}",
                        f"__TOP_LOG_SHA256={digest(top_log)}",
                        "MACHINE_A_LOW tasks_passed=24 tasks_expected=24 failures=0",
                        "MACHINE_A_LOW FREEZE: ALL PASS",
                        "",
                    )
                ),
            )
            print("MACHINE_A_LOW FREEZE: ALL PASS", flush=True)
            return 0
        if not process_group_exists(pgid):
            raise RuntimeError(
                f"machine-A run exited with only {len(passed)}/24 selected tasks"
            )
        print(
            f"MACHINE_A_LOW waiting passed={len(passed)}/24 utc={utc_now()}",
            flush=True,
        )
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"MACHINE_A_LOW FREEZE: FAIL ({error})", file=sys.stderr, flush=True)
        raise SystemExit(1)
