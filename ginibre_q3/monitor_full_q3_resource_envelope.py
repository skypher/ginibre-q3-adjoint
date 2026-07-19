#!/usr/bin/env python3
"""Fail closed when a detached full-Q3 run exceeds a host memory envelope."""

from __future__ import annotations

import argparse
import os
import signal
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def memory_kib() -> tuple[int, int]:
    values: dict[str, int] = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, raw = line.split(":", 1)
        values[key] = int(raw.split()[0])
    available = values["MemAvailable"]
    swap_used = values["SwapTotal"] - values["SwapFree"]
    return available, swap_used


def group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    return True


def terminate_group(pgid: int) -> str:
    if not group_exists(pgid):
        return "already_exited"
    os.killpg(pgid, signal.SIGTERM)
    deadline = time.monotonic() + 10.0
    while group_exists(pgid) and time.monotonic() < deadline:
        time.sleep(0.2)
    if group_exists(pgid):
        os.killpg(pgid, signal.SIGKILL)
        return "sigkill_after_sigterm"
    return "sigterm"


def atomic_write(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orchestrator-pid", type=int, required=True)
    parser.add_argument("--run-token", required=True)
    parser.add_argument("--min-available-kib", type=int, required=True)
    parser.add_argument("--max-swap-used-kib", type=int, required=True)
    parser.add_argument("--poll-seconds", type=float, default=30.0)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if min(
        args.orchestrator_pid,
        args.min_available_kib,
        args.max_swap_used_kib,
    ) < 1 or args.poll_seconds <= 0:
        parser.error("PID, thresholds, and polling interval must be positive")

    pid = args.orchestrator_pid
    command_path = Path(f"/proc/{pid}/cmdline")
    if not command_path.is_file():
        raise RuntimeError(f"orchestrator PID is not live: {pid}")
    command = command_path.read_bytes().replace(b"\0", b" ").decode(errors="strict")
    if args.run_token not in command:
        raise RuntimeError(f"PID {pid} does not contain run token {args.run_token}")
    pgid = os.getpgid(pid)
    if pgid != pid:
        raise RuntimeError(f"orchestrator is not its process-group leader: {pid}/{pgid}")

    while group_exists(pgid):
        available, swap_used = memory_kib()
        print(
            f"RESOURCE_ENVELOPE utc={utc_now()} available_kib={available} "
            f"swap_used_kib={swap_used}",
            flush=True,
        )
        if (
            available < args.min_available_kib
            or swap_used > args.max_swap_used_kib
        ):
            method = terminate_group(pgid)
            atomic_write(
                args.receipt,
                "\n".join(
                    (
                        f"__UTC={utc_now()}",
                        "__STATUS=ABORTED_RESOURCE_ENVELOPE",
                        f"__RUN_TOKEN={args.run_token}",
                        f"__PROCESS_GROUP={pgid}",
                        f"__AVAILABLE_KIB={available}",
                        f"__SWAP_USED_KIB={swap_used}",
                        f"__MIN_AVAILABLE_KIB={args.min_available_kib}",
                        f"__MAX_SWAP_USED_KIB={args.max_swap_used_kib}",
                        f"__STOP_METHOD={method}",
                        "",
                    )
                ),
            )
            print("RESOURCE_ENVELOPE VERIFICATION: ABORTED", flush=True)
            return 2
        time.sleep(args.poll_seconds)

    atomic_write(
        args.receipt,
        "\n".join(
            (
                f"__UTC={utc_now()}",
                "__STATUS=PROCESS_GROUP_EXITED_WITHIN_ENVELOPE",
                f"__RUN_TOKEN={args.run_token}",
                f"__PROCESS_GROUP={pgid}",
                "RESOURCE_ENVELOPE VERIFICATION: PASS",
                "",
            )
        ),
    )
    print("RESOURCE_ENVELOPE VERIFICATION: PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
