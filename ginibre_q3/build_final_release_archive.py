#!/usr/bin/env python3
"""Build a deterministic archive of the exact tracked publication package."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
from pathlib import Path
import subprocess
import tarfile


FIXED_MTIME = 1_784_419_200


def git(repo: Path, *arguments: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(repo), *arguments])


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FINAL_RELEASE_ARCHIVE: FAIL: {message}")


def add_bytes(archive: tarfile.TarFile, name: str, data: bytes, mode: int) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mode = mode
    info.mtime = FIXED_MTIME
    info.uid = 0
    info.gid = 0
    info.uname = "root"
    info.gname = "root"
    archive.addfile(info, io.BytesIO(data))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-log", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    project = Path(__file__).resolve().parent
    repo = Path(git(project, "rev-parse", "--show-toplevel").decode().strip())
    status = git(repo, "status", "--porcelain=v1", "--untracked-files=all")
    require(not status, "worktree must be clean")
    require(args.replay_log.is_absolute(), "--replay-log must be absolute")
    require(args.output.is_absolute(), "--output must be absolute")
    require(repo not in args.output.parents, "output archive must be outside the repository")
    require(args.replay_log.is_file(), "final replay log is absent")
    replay = args.replay_log.read_bytes()
    require(
        b"FINAL_PUBLICATION_REPLAY: ALL PASS" in replay,
        "replay log lacks the terminal success marker",
    )
    require(not args.output.exists(), "refusing to overwrite output archive")

    tracked_raw = git(repo, "ls-files", "-z", "--", "ginibre_q3")
    tracked = sorted(
        Path(item.decode()) for item in tracked_raw.split(b"\0") if item
    )
    require(tracked, "no tracked publication files found")
    for relative in tracked:
        require((repo / relative).is_file(), f"tracked file is absent: {relative}")
        require(not (repo / relative).is_symlink(), f"symlink rejected: {relative}")

    commit = git(repo, "rev-parse", "HEAD").decode().strip()
    manifest = project / "certificates/full_q3/full_q3_source_manifest.sha256"
    artifact_manifest = project / "publication_artifacts.sha256"
    require(artifact_manifest.is_file(), "publication artifact manifest is absent")
    metadata = (
        "GINIBRE_Q3_FINAL_RELEASE\n"
        f"commit={commit}\n"
        f"source_manifest_sha256={hashlib.sha256(manifest.read_bytes()).hexdigest()}\n"
        f"publication_artifact_manifest_sha256="
        f"{hashlib.sha256(artifact_manifest.read_bytes()).hexdigest()}\n"
        f"replay_log_sha256={hashlib.sha256(replay).hexdigest()}\n"
        f"source_date_epoch={FIXED_MTIME}\n"
    ).encode()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w") as archive:
                add_bytes(archive, "GINIBRE_Q3_FINAL_RELEASE/ARCHIVE_METADATA.txt", metadata, 0o644)
                add_bytes(archive, "GINIBRE_Q3_FINAL_RELEASE/FINAL_REPLAY_LOG.txt", replay, 0o644)
                for relative in tracked:
                    source = repo / relative
                    mode = 0o755 if os.access(source, os.X_OK) else 0o644
                    add_bytes(
                        archive,
                        f"GINIBRE_Q3_FINAL_RELEASE/{relative.as_posix()}",
                        source.read_bytes(),
                        mode,
                    )

    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(f"FINAL_RELEASE_ARCHIVE files={len(tracked) + 2}")
    print(f"FINAL_RELEASE_ARCHIVE sha256={digest} path={args.output}")
    print("FINAL_RELEASE_ARCHIVE: ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
