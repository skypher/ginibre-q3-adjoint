#!/usr/bin/env python3
"""Fail-closed composition audit for the complete Ginibre-Q3 replay.

The serial clean-room driver has 200 arithmetic/build stages before its PDF
passes.  This verifier accepts the same stage ledger split as follows:

* stages 1--49 from the final-source machine-C clean replay;
* stages 50--83 from a byte-equivalent older machine-C clean snapshot;
* stages 84--154 from three disjoint current-source fleet partitions;
* stages 155--184 from the current-source classical post-frontier replay; and
* stages 185--200 from the current-source exceptional replay.

It reconstructs the authoritative stage list from ``clean_room_replay.py``,
checks every command, environment override, required output marker, source
identity, exit marker, and log hash, and then invokes the driver's exact
main-theorem certificate-coverage audit on the composed ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


CURRENT_DRIVER_SHA256 = (
    "98cdbb03930cf2bc139364d59f703a3352b0007803c5d47c15282ce696f480eb"
)
CURRENT_REPLAY_MANIFEST_SHA256 = (
    "030229492cda95a0bb9c1d3183c3880880c2e29e8ee9e2cec717d32733b72ed5"
)
LEGACY_DRIVER_SHA256 = (
    "51e028eac01c923c935266171f4a27bd25bec49f39b30fea511f406222d8d55a"
)
LEGACY_REPLAY_MANIFEST_SHA256 = (
    "6f7ff746cb3331500d25972d65341f64a9ff2250b142fd8499c652a392ba62f8"
)
FLEET_SCRIPT_SHA256 = (
    "28b51756e053cb160e9cf5dde4f97a27ae416f084f7ef23ffedd02409098464a"
)
# Machine A's build transcript names this earlier orchestration-file digest.
# That Python file is not a compiler input: the transcript separately records
# every compile command and binds every C++ source and resulting binary.  No
# result is inferred from the earlier script.  The actual run top log and every
# task are checked against FLEET_SCRIPT_SHA256 below.  Pinning the historical
# line, rather than rewriting it or pretending that it is the run ledger,
# preserves provenance without granting it a proof role.
MACHINE_A_BUILD_LEDGER_SHA256 = (
    "5883f7d8961d34adde5ca78335480b8f412b30de35a4eadcd1a47bcfe3bb1922"
)
FLEET_PARTITION_SCRIPT_SHA256 = (
    "784c47f797a61e805d509c8e6262ccf5abd7f3bf7757a3b8fb4dad6d84b80adf"
)
CURRENT_BC_CODE_MANIFEST_SHA256 = (
    "dd8879a180b19702bfa562c79eab5483a6d4cfe9d013038042b8d06947f5cffa"
)
LEGACY_BC_CODE_MANIFEST_SHA256 = (
    "cfa293c681de81fec8cd35e96db531df3b1f7ce7bd75f29009a9ad7b50c69d0f"
)
GEN_HEADER_SHA256 = (
    "3c3e6c8d3b53cf725e1ac1c6d772a97e9ddeb9dfc9e4925d5888eb0a6a42ae83"
)
LIE_DATA_PY_SHA256 = (
    "ca8425b7c07edc5602210c4a8d6f6ad01f95ca34012c2c5b4729f745bdbf6ad3"
)
LIE_DATA_HEADER_SHA256 = (
    "ba30755944b7afa31ff8d7de00ff1b3384274e3bacd59edadb5d3c6f1025f576"
)
FRONTIER_KERNEL_AUDIT_SHA256 = (
    "c48157f690d0f23fca66f267183e54d8d69a9d39853794dc824439936237baf8"
)
POST_WRAPPER_SHA256 = (
    "77ad5695f7bafd0d9e1569902d9f303773c59c3192c88a0338a34b8d4933b8d9"
)
PART_III_TEX_SHA256 = (
    "48d269b8d0ab67d15e49f99256fd0e36dd559691bf2039d2c373ca99fc483c02"
)

CHANGED_SNAPSHOT_PATHS = {
    "Makefile",
    "README.md",
    "REPLAY.md",
    "clean_room_replay.py",
    "paper.tex",
    "character_ring_iter/Makefile",
    "character_ring_iter/post_m29_bc_layered_mgf_mpfr.cpp",
    "character_ring_iter/verify_chain.py",
}

PREFLIGHT_MARKERS = (
    "[replay] artifact integrity: PASS",
    "[replay] reference archive: PASS",
    "[replay] proof structure: PASS",
    "[replay] B/C caller-range closure: PASS",
    "[replay] certificate identity ledger: PASS",
    "[replay] post-m29 artifact boundary: PASS",
    "[replay] accepted replay status: PASS",
    "[replay] proof-spine input coverage: PASS",
    "[replay] purged ",
)


class VerificationFailure(RuntimeError):
    pass


@dataclass(frozen=True)
class StageSpec:
    number: int
    name: str
    command: tuple[str, ...]
    environment: tuple[tuple[str, str], ...]
    required: tuple[str, ...]
    certificates: tuple[str, ...]


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise VerificationFailure(f"missing file: {path}")
    return path


def require_fragment(text: str, fragment: str, locus: Path) -> None:
    if fragment not in text:
        raise VerificationFailure(f"missing {fragment!r} in {locus}")


def require_line_once(text: str, expected: str, locus: Path) -> None:
    """Require one complete output line, not merely an embedded fragment."""

    count = sum(line == expected for line in text.splitlines())
    if count != 1:
        raise VerificationFailure(
            f"expected one exact line {expected!r} in {locus}, found {count}"
        )


def require_replay_pass_line(text: str, stage_name: str, locus: Path) -> None:
    """Require one canonical timed PASS line for a clean-driver stage."""

    pattern = re.compile(
        rf"^\[replay\] {re.escape(stage_name)}: PASS "
        rf"\(([0-9]+(?:\.[0-9]+)?)s\)$",
        re.MULTILINE,
    )
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise VerificationFailure(
            f"expected one exact replay PASS line for {stage_name} in {locus}, "
            f"found {len(matches)}"
        )


def require_hash_record(
    text: str,
    expected_hash: str,
    path_suffix: str,
    locus: Path,
) -> None:
    """Require one sha256sum-style record, independently of archive relocation.

    Fleet builds are performed in a fresh temporary source tree and their raw
    logs intentionally retain that original absolute path.  The authenticated
    evidence is subsequently copied into the repository certificate archive.
    Requiring the archived absolute path to occur in the immutable build log
    would therefore make a valid replay unverifiable after relocation.  We
    instead pin the digest and the complete path suffix, while callers also
    hash the corresponding archived file itself.
    """

    pattern = re.compile(
        rf"^{re.escape(expected_hash)}\s+[*]?(?:.*/)?"
        rf"{re.escape(path_suffix)}$",
        re.MULTILINE,
    )
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise VerificationFailure(
            f"expected one hash record for {path_suffix} in {locus}, "
            f"found {len(matches)}"
        )


def parse_manifest(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    pattern = re.compile(r"^([0-9a-f]{64})\s+[*]?(.+?)\s*$")
    for line_number, line in enumerate(
        require_file(path).read_text(errors="strict").splitlines(), 1
    ):
        match = pattern.fullmatch(line)
        if match is None:
            raise VerificationFailure(f"malformed manifest line {path}:{line_number}")
        value, name = match.groups()
        if name in rows:
            raise VerificationFailure(f"duplicate manifest path in {path}: {name}")
        rows[name] = value
    if not rows:
        raise VerificationFailure(f"empty manifest: {path}")
    return rows


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise VerificationFailure(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def collect_stage_specs(driver, root: Path, threads: int) -> list[StageSpec]:
    class Collector:
        def __init__(self) -> None:
            self.specs: list[StageSpec] = []

        def run(
            self,
            name: str,
            command: list[str],
            _cwd: Path,
            required_text: str | tuple[str, ...] | None = None,
            environment_overrides: dict[str, str] | None = None,
            certificates: tuple[str, ...] = (),
            **_kwargs: object,
        ) -> Path:
            required = (
                ()
                if required_text is None
                else ((required_text,) if isinstance(required_text, str) else required_text)
            )
            number = len(self.specs) + 2  # stage 1 is replay-build
            stage = StageSpec(
                number=number,
                name=name,
                command=tuple(map(str, command)),
                environment=tuple(sorted((environment_overrides or {}).items())),
                required=tuple(required),
                certificates=tuple(certificates),
            )
            self.specs.append(stage)
            return root.parent / "replay-logs" / f"{number:02d}_{name}.log"

    collector = Collector()
    driver.stable_moment_recurrence_replay(collector, root)
    driver.cpp_endpoint_replays(collector, root, threads)
    driver.bc_residual_certificate_replays(collector, root, threads)
    driver.classical_replays(collector, root)
    driver.exceptional_replays(collector, root, threads)
    specs = collector.specs
    if len(specs) != 199 or specs[0].number != 2 or specs[-1].number != 200:
        raise VerificationFailure(
            f"clean replay stage ledger changed: count={len(specs)}"
        )
    boundaries = {
        2: "stable_moment_recurrence_comparison",
        49: "bc_b14_j37_badshape_gmp",
        50: "bc_b_eighth_frontier_gmp",
        83: "bc_twentyfourth_frontier_c",
        84: "bc_twentyfifth_frontier_b_15_17",
        154: "bc_twentyninth_b_absorption_gmp",
        155: "classical_low_rank_exact_b2_c2_post_m29_tail_certificate",
        184: "classical_post_m29_bc_first_hit",
        185: "exceptional_source_pairing_g2",
        200: "exceptional_e8_finite_bridge",
    }
    by_number = {stage.number: stage.name for stage in specs}
    if any(by_number.get(number) != name for number, name in boundaries.items()):
        raise VerificationFailure(f"clean replay partition boundaries changed: {by_number}")
    if len({stage.name for stage in specs}) != len(specs):
        raise VerificationFailure("duplicate clean replay stage name")
    return specs


def normalize_argument(argument: str) -> str:
    if "/replay-logs/" in argument:
        return "<REPLAY_LOG>/" + Path(argument).name
    marker = "/ginibre_q3/"
    if marker in argument:
        # An execution snapshot is stored below the live project tree, so its
        # absolute path contains more than one ``/ginibre_q3/`` component.
        # The command is rooted at the innermost project tree; normalizing the
        # first component incorrectly retains the archive prefix and rejects
        # a byte-identical archived replay.
        return "<ROOT>/" + argument.rsplit(marker, 1)[1]
    if argument.endswith("/ginibre_q3"):
        return "<ROOT>"
    return argument


def normalized_command(line: str) -> tuple[str, ...]:
    try:
        return tuple(normalize_argument(value) for value in shlex.split(line))
    except ValueError as error:
        raise VerificationFailure(f"cannot parse command line: {line}") from error


def find_stage_log(stage_dir: Path, name: str) -> Path:
    matches = sorted(stage_dir.glob(f"*_{name}.log"))
    if len(matches) != 1:
        raise VerificationFailure(
            f"expected one stage log for {name} in {stage_dir}, found {len(matches)}"
        )
    return matches[0]


def verify_stage_log(
    stage: StageSpec,
    stage_dir: Path,
    top_log_text: str,
    top_log_path: Path,
    require_exit: bool = False,
) -> str:
    log = find_stage_log(stage_dir, stage.name)
    text = log.read_text(errors="strict")
    lines = text.splitlines()
    if not lines or not lines[0].startswith("COMMAND: "):
        raise VerificationFailure(f"missing command header in {log}")
    actual_command = normalized_command(lines[0].removeprefix("COMMAND: "))
    expected_command = tuple(normalize_argument(value) for value in stage.command)
    if actual_command != expected_command:
        raise VerificationFailure(
            f"command mismatch for {stage.name}:\n"
            f"expected={expected_command}\nactual={actual_command}"
        )
    expected_environment = " ".join(f"{key}={value}" for key, value in stage.environment)
    environment_lines = [line.removeprefix("ENV: ") for line in lines if line.startswith("ENV: ")]
    if stage.environment:
        if environment_lines != [expected_environment]:
            raise VerificationFailure(
                f"environment mismatch for {stage.name}: {environment_lines}"
            )
    elif environment_lines:
        raise VerificationFailure(f"unexpected environment override for {stage.name}")
    for marker in stage.required:
        require_fragment(text, marker, log)
    if require_exit:
        require_line_once(text, "__EXIT_STATUS=0", log)
    require_replay_pass_line(top_log_text, stage.name, top_log_path)
    return digest(log)


def verify_current_source_identity(root: Path, legacy_root: Path) -> tuple[dict[str, str], dict[str, str]]:
    current_driver = root / "clean_room_replay.py"
    current_manifest_path = root / "replay_sources.sha256"
    legacy_driver = legacy_root / "clean_room_replay.py"
    legacy_manifest_path = legacy_root / "replay_sources.sha256"
    if not legacy_manifest_path.is_file():
        legacy_manifest_path = legacy_root / "replay_sources_snapshot.txt"
    expected_hashes = (
        (current_driver, CURRENT_DRIVER_SHA256),
        (current_manifest_path, CURRENT_REPLAY_MANIFEST_SHA256),
        (legacy_driver, LEGACY_DRIVER_SHA256),
        (legacy_manifest_path, LEGACY_REPLAY_MANIFEST_SHA256),
    )
    for path, expected in expected_hashes:
        actual = digest(require_file(path))
        if actual != expected:
            raise VerificationFailure(
                f"identity mismatch for {path}: expected={expected}, actual={actual}"
            )
    current = parse_manifest(current_manifest_path)
    legacy = parse_manifest(legacy_manifest_path)
    if len(current) != 85 or set(current) != set(legacy):
        raise VerificationFailure(
            f"snapshot manifest shape mismatch: current={len(current)}, legacy={len(legacy)}"
        )
    changed = {name for name in current if current[name] != legacy[name]}
    if changed != CHANGED_SNAPSHOT_PATHS:
        raise VerificationFailure(
            f"unexpected snapshot source delta: {sorted(changed)}"
        )
    for name, expected in current.items():
        path = root / name
        if digest(require_file(path)) != expected:
            raise VerificationFailure(f"current manifest mismatch: {name}")
    return current, legacy


def verify_part_iii_scope(root: Path) -> Fraction:
    source = require_file(root / "full_q3_extension.tex")
    if digest(source) != PART_III_TEX_SHA256:
        raise VerificationFailure("Part III source identity mismatch")
    text = source.read_text(errors="strict")
    fragments = (
        r"\title[Full adjoint-generated Ginibre $Q_3$]",
        r"Every element of $\Qcone_G$ is a real continuous positive-definite function",
        r"Then $\mathcal P_G$ does not satisfy Ginibre's",
        r"\begin{theorem}[Full adjoint-generated Ginibre $Q_3$]",
        r"f_1,\ldots,f_r\in\Qcone_G",
        r"It does not replace $\Qcone_G$ by the cone of",
        r"impossible: Proposition~\ref{prop:all-pd-counterexample}",
    )
    for fragment in fragments:
        require_fragment(text, fragment, source)
    require_fragment(text, r"-\frac8{279}", source)

    c1, c2, c3 = Fraction(7, 31), Fraction(19, 31), Fraction(19, 31)
    fourth = Fraction(2, 15) * (c1 * c1 + c2 * c2 + c3 * c3)
    fourth -= Fraction(1, 15) * (c1 * c2 + c1 * c3 + c2 * c3)
    signed = 2 * (fourth + Fraction(1, 9) * (c1 * c1 - c2 * c2 - c3 * c3))
    if fourth != Fraction(61, 961) or signed != Fraction(-8, 279):
        raise VerificationFailure("exact SO(3) full-cone obstruction mismatch")
    return signed


def make_rule(makefile: Path, target: str) -> str:
    lines = makefile.read_text(errors="strict").splitlines()
    for index, line in enumerate(lines):
        if line.startswith(target + ":"):
            block = [line]
            cursor = index + 1
            while cursor < len(lines) and (
                lines[cursor].startswith("\t") or not lines[cursor].strip()
            ):
                block.append(lines[cursor])
                cursor += 1
            return "\n".join(block).strip()
    raise VerificationFailure(f"missing Makefile rule {target} in {makefile}")


def verify_header_generator(root: Path) -> None:
    character = root / "character_ring_iter"
    identities = (
        (character / "gen_header.py", GEN_HEADER_SHA256),
        (character / "lie_data.py", LIE_DATA_PY_SHA256),
        (character / "lie_data.h", LIE_DATA_HEADER_SHA256),
    )
    for path, expected in identities:
        actual = digest(require_file(path))
        if actual != expected:
            raise VerificationFailure(
                f"Lie-data generator identity mismatch for {path}: {actual}"
            )
    with tempfile.TemporaryDirectory(prefix="full-q3-lie-data-") as raw:
        temporary = Path(raw)
        shutil.copy2(character / "gen_header.py", temporary / "gen_header.py")
        shutil.copy2(character / "lie_data.py", temporary / "lie_data.py")
        completed = subprocess.run(
            [sys.executable, str(temporary / "gen_header.py")],
            cwd=temporary,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            raise VerificationFailure(
                f"Lie-data header regeneration failed: {completed.stdout}"
            )
        regenerated = temporary / "lie_data.h"
        if digest(require_file(regenerated)) != LIE_DATA_HEADER_SHA256:
            raise VerificationFailure("Lie-data generator does not reproduce pinned header")


def verify_frontier_kernel(root: Path) -> int:
    """Compile the actual reverse-Pieri kernel against an independent oracle."""

    character = root / "character_ring_iter"
    audit_source = require_file(
        character / "verify_frontier_horizontal_strip_kernel.cpp"
    )
    if digest(audit_source) != FRONTIER_KERNEL_AUDIT_SHA256:
        raise VerificationFailure("frontier horizontal-strip audit identity mismatch")

    kernel_blocks: list[str] = []
    for ordinal in ("twentyfifth", "twentysixth", "twentyseventh", "twentyeighth"):
        source = require_file(character / f"post_m29_bc_{ordinal}_frontier_gmp.cpp")
        text = source.read_text(errors="strict")
        start = text.find("std::string encode(")
        end = text.find("mpz_class compute_b_width_bad_count(")
        if start < 0 or end <= start:
            raise VerificationFailure(f"cannot isolate reverse-Pieri kernel in {source}")
        kernel_blocks.append(text[start:end])
    if len(set(kernel_blocks)) != 1:
        raise VerificationFailure("H25--H28 reverse-Pieri kernels are not identical")

    with tempfile.TemporaryDirectory(prefix="full-q3-frontier-kernel-") as raw:
        binary = Path(raw) / "verify_frontier_horizontal_strip_kernel"
        compile_run = subprocess.run(
            [
                "g++",
                "-O2",
                "-std=c++20",
                "-Wall",
                "-Werror",
                "-fopenmp",
                str(audit_source),
                "-lgmpxx",
                "-lgmp",
                "-o",
                str(binary),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if compile_run.returncode != 0:
            raise VerificationFailure(
                f"frontier horizontal-strip audit build failed: {compile_run.stdout}"
            )
        replay = subprocess.run(
            [str(binary)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        required = (
            "FRONTIER_HORIZONTAL_STRIP_KERNEL partitions_checked=1595 "
            "predecessor_sets_checked=1595 maximum_live_part=121 "
            "encoding_limit=255 failures=0"
        )
        if (
            replay.returncode != 0
            or replay.stdout.count(required) != 1
            or replay.stdout.count(
                "FRONTIER_HORIZONTAL_STRIP_KERNEL VERIFICATION: ALL PASS"
            )
            != 1
            or replay.stdout.splitlines()[-1:] != ["__EXIT_STATUS=0"]
        ):
            raise VerificationFailure(
                f"frontier horizontal-strip audit failed: {replay.stdout}"
            )
    return 1595


def verify_replay_build_dependencies(root: Path, driver) -> int:
    """Require every direct or generated Make dependency to be manifested."""

    character = root / "character_ring_iter"
    raw_lines = (character / "Makefile").read_text(errors="strict").splitlines()
    logical: list[str] = []
    current = ""
    for line in raw_lines:
        current = current + line.lstrip() if current else line
        if current.endswith("\\"):
            current = current[:-1] + " "
        else:
            logical.append(current)
            current = ""
    if current:
        raise VerificationFailure("unterminated Makefile continuation")
    variables: dict[str, list[str]] = {}
    for line in logical:
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", line)
        if match is not None:
            variables[match.group(1)] = match.group(2).split()

    def expand(tokens: list[str]) -> list[str]:
        expanded: list[str] = []
        for token in tokens:
            match = re.fullmatch(r"\$\(([^)]+)\)", token)
            expanded.extend(variables.get(match.group(1), [token]) if match else [token])
        return expanded

    rules: dict[str, list[str]] = {}
    for line in logical:
        if not line or line[0].isspace() or ":" not in line or line.startswith(".PHONY"):
            continue
        left, right = line.split(":", 1)
        if "%" in left:
            continue
        for target in left.split():
            rules[target] = expand(right.split())
    replay_targets = rules.get("replay-build", [])
    reachable: set[str] = set()
    frontier = list(replay_targets)
    while frontier:
        target = frontier.pop()
        if target in reachable:
            continue
        reachable.add(target)
        frontier.extend(dependency for dependency in rules.get(target, []) if dependency in rules)
    inputs = {
        (character / dependency).resolve()
        for target in reachable
        for dependency in rules.get(target, [])
        if (character / dependency).is_file()
    }
    manifested: set[Path] = set()
    for manifest in root.rglob("*.sha256"):
        for line in manifest.read_text(errors="strict").splitlines():
            match = driver.HASH_LINE.match(line)
            if match is not None:
                manifested.add(driver.resolve_manifest_path(root, manifest, match.group(2)))
    missing = sorted(path.relative_to(root) for path in inputs - manifested)
    if (
        len(replay_targets) != 47
        or len(reachable) != 48
        or len(inputs) != 52
        or missing
    ):
        raise VerificationFailure(
            "replay-build dependency audit changed: "
            f"targets={len(replay_targets)}, rules={len(reachable)}, "
            f"inputs={len(inputs)}, missing={missing}"
        )
    return len(inputs)


def verify_legacy_sources(
    root: Path,
    legacy_root: Path,
    _current_manifest: dict[str, str],
    _legacy_manifest: dict[str, str],
    stages: list[StageSpec],
) -> None:
    binaries = {Path(stage.command[0]).name for stage in stages}
    expected_binaries = {
        "post_m29_bc_b_eighth_frontier_gmp",
        *(f"post_m29_bc_{name}_frontier_gmp" for name in (
            "ninth", "tenth", "eleventh", "twelfth", "thirteenth",
            "fourteenth", "fifteenth", "sixteenth", "seventeenth",
            "eighteenth", "nineteenth", "twentieth", "twentyfirst",
            "twentysecond", "twentythird", "twentyfourth",
        )),
    }
    if binaries != expected_binaries:
        raise VerificationFailure(f"legacy binary ledger changed: {sorted(binaries)}")
    current_makefile = root / "character_ring_iter/Makefile"
    legacy_makefile = legacy_root / "character_ring_iter/Makefile"
    current_code_path = root / "certificates/post_m29/bc_interval_tail_code.sha256"
    legacy_code_path = legacy_root / "certificates/post_m29/bc_interval_tail_code.sha256"
    if not legacy_code_path.is_file():
        legacy_code_path = (
            legacy_root
            / "certificates/post_m29/bc_interval_tail_code_snapshot.txt"
        )
    if digest(require_file(current_code_path)) != CURRENT_BC_CODE_MANIFEST_SHA256:
        raise VerificationFailure("current B/C code-manifest identity mismatch")
    if digest(require_file(legacy_code_path)) != LEGACY_BC_CODE_MANIFEST_SHA256:
        raise VerificationFailure("legacy B/C code-manifest identity mismatch")
    current_code = parse_manifest(current_code_path)
    legacy_code = parse_manifest(legacy_code_path)
    code_changed = {
        name for name in set(current_code) & set(legacy_code)
        if current_code[name] != legacy_code[name]
    }
    if set(current_code) != set(legacy_code) or code_changed != {
        "../../character_ring_iter/Makefile",
        "../../character_ring_iter/post_m29_bc_layered_mgf_mpfr.cpp",
    }:
        raise VerificationFailure(
            f"unexpected B/C code-manifest delta: {sorted(code_changed)}"
        )
    if digest(legacy_makefile) != legacy_code["../../character_ring_iter/Makefile"]:
        raise VerificationFailure("staged legacy Makefile does not match its code manifest")
    for binary in sorted(binaries):
        source_name = f"character_ring_iter/{binary}.cpp"
        manifest_name = f"../../{source_name}"
        current_source = root / source_name
        legacy_source = legacy_root / source_name
        current_hash = digest(require_file(current_source))
        legacy_hash = digest(require_file(legacy_source))
        if not (
            current_hash
            == legacy_hash
            == current_code[manifest_name]
            == legacy_code[manifest_name]
        ):
            raise VerificationFailure(f"legacy/current source mismatch: {source_name}")
        if make_rule(current_makefile, binary) != make_rule(legacy_makefile, binary):
            raise VerificationFailure(f"legacy/current build-rule mismatch: {binary}")


def verify_top_log(path: Path, markers: tuple[str, ...] = ()) -> str:
    text = require_file(path).read_text(errors="strict")
    for marker in markers:
        require_fragment(text, marker, path)
    if "[replay] FAILURE:" in text or "Traceback (most recent call last)" in text:
        raise VerificationFailure(f"failure marker in {path}")
    return text


def fleet_stage_name(task_name: str) -> str:
    finals = {
        "twentyeighth_b_absorption": "bc_twentyeighth_b_absorption_gmp",
        "twentyninth_c_frontier": "bc_twentyninth_c_frontier_gmp",
        "twentyninth_b_absorption": "bc_twentyninth_b_absorption_gmp",
    }
    if task_name in finals:
        return finals[task_name]
    match = re.fullmatch(r"(twenty(?:fifth|sixth|seventh|eighth))_(b_\d+_\d+|c)", task_name)
    if match is None:
        raise VerificationFailure(f"unrecognized fleet task name: {task_name}")
    ordinal, suffix = match.groups()
    return f"bc_{ordinal}_frontier_{suffix}"


def frontier_partition_name(task_name: str) -> str:
    if task_name.startswith(("twentyfifth_", "twentysixth_")):
        return "machine_a_low"
    if task_name == "twentyseventh_c":
        return "machine_b_mid"
    if task_name.startswith("twentyseventh_b_"):
        high = int(task_name.rsplit("_", 1)[1])
        return "machine_b_mid" if high <= 32 else "machine_c_heavy"
    return "machine_c_heavy"


def frontier_partition_tasks(ledger, threads: int, partition: str):
    tasks = ledger.make_tasks(threads)
    if len(tasks) != 71 or len({task.name for task in tasks}) != 71:
        raise VerificationFailure("fleet task ledger is not 71 unique tasks")
    selected = [task for task in tasks if frontier_partition_name(task.name) == partition]
    expected_counts = {
        "machine_a_low": 24,
        "machine_b_mid": 7,
        "machine_c_heavy": 40,
    }
    if partition not in expected_counts or len(selected) != expected_counts[partition]:
        raise VerificationFailure(
            f"fleet partition cardinality mismatch for {partition}: {len(selected)}"
        )
    return selected


def verify_fleet_task_semantics(
    tasks,
    expected_stages: list[StageSpec],
    threads: int,
    partition: str,
) -> None:
    """Compare the fleet ledger independently with the clean-driver ledger.

    OpenMP thread count is intentionally host-dependent.  It is the only
    permitted environment difference: commands, finite rank shards, scope
    markers, and every other environment override must agree exactly.
    """

    stages_by_name = {stage.name: stage for stage in expected_stages}
    if len(stages_by_name) != len(expected_stages):
        raise VerificationFailure(f"duplicate expected clean stage in {partition}")
    for task in tasks:
        stage_name = fleet_stage_name(task.name)
        if stage_name not in stages_by_name:
            raise VerificationFailure(
                f"fleet task {task.name} has no clean-driver stage in {partition}"
            )
        stage = stages_by_name[stage_name]
        if stage.command != (f"./{task.binary}",):
            raise VerificationFailure(
                f"fleet command differs from clean driver for {task.name}: "
                f"fleet=./{task.binary}, clean={stage.command}"
            )
        if task.source != f"{task.binary}.cpp":
            raise VerificationFailure(
                f"fleet source/binary naming mismatch for {task.name}"
            )
        if tuple(task.required) != stage.required:
            raise VerificationFailure(
                f"fleet scope markers differ from clean driver for {task.name}"
            )
        fleet_environment = dict(task.environment)
        clean_environment = dict(stage.environment)
        fleet_threads = fleet_environment.pop("OMP_NUM_THREADS", None)
        clean_environment.pop("OMP_NUM_THREADS", None)
        if fleet_threads != str(threads):
            raise VerificationFailure(
                f"fleet OpenMP ledger mismatch for {task.name}: {fleet_threads}"
            )
        if fleet_environment != clean_environment:
            raise VerificationFailure(
                f"fleet environment differs from clean driver for {task.name}: "
                f"fleet={fleet_environment}, clean={clean_environment}"
            )


def expected_fleet_scope_line(task, threads: int) -> str:
    """Return the complete deterministic scope-header line for a fleet task."""

    scope = task.required[0]
    environment = dict(task.environment)
    if "RUN_B" not in environment:
        return scope
    if environment == {
        "RUN_B": "1",
        "RUN_C": "0",
        "B_Q_MIN": environment.get("B_Q_MIN"),
        "B_Q_MAX": environment.get("B_Q_MAX"),
        "OMP_NUM_THREADS": str(threads),
    }:
        return f"{scope} omp_max_threads={threads}"
    if environment == {
        "RUN_B": "0",
        "RUN_C": "1",
        "OMP_NUM_THREADS": str(threads),
    }:
        return (
            f"{scope} b_q_min=-1000000 b_q_max=1000000 "
            f"omp_max_threads={threads}"
        )
    raise VerificationFailure(f"unrecognized fleet environment for {task.name}")


def frontier_certificate_name(task_name: str) -> str:
    finals = {
        "twentyeighth_b_absorption": (
            "bc_twentyeighth_b_absorption_gmp_certificate.log"
        ),
        "twentyninth_c_frontier": "bc_twentyninth_c_frontier_gmp_certificate.log",
        "twentyninth_b_absorption": (
            "bc_twentyninth_b_absorption_gmp_certificate.log"
        ),
    }
    if task_name in finals:
        return finals[task_name]
    for ordinal in ("twentyfifth", "twentysixth", "twentyseventh", "twentyeighth"):
        if task_name.startswith(ordinal + "_"):
            return f"bc_{ordinal}_frontier_gmp_certificate.log"
    raise VerificationFailure(f"no accepted certificate for fleet task {task_name}")


def arithmetic_output_lines(text: str) -> list[str]:
    """Remove orchestration/scope lines, retaining every arithmetic datum."""

    return [
        line
        for line in text.splitlines()
        if line
        and not line.startswith("__")
        and not line.startswith("BC_")
        and not line.startswith("SUMMARY ")
    ]


def normalized_fleet_scope(scope: str) -> str:
    """Remove only the host-dependent OpenMP suffix from a scope line."""

    return re.sub(r" omp_max_threads=\d+$", "", scope)


def accepted_arithmetic_blocks(text: str, locus: Path) -> dict[str, tuple[str, ...]]:
    """Index the paper-cited transcript by its unique complete run scopes.

    The accepted H25--H28 certificates concatenate one C run and all B shards.
    Thread counts differ between the historical transcript and the clean fleet,
    so the terminal ``omp_max_threads`` field is normalized; every other scope
    field and every arithmetic output line remain byte-for-byte significant.
    """

    raw_lines = text.splitlines()
    starts = [index for index, line in enumerate(raw_lines) if line.startswith("BC_")]
    if not starts:
        raise VerificationFailure(f"accepted certificate has no scope blocks: {locus}")
    blocks: dict[str, tuple[str, ...]] = {}
    for block_index, start in enumerate(starts):
        end = starts[block_index + 1] if block_index + 1 < len(starts) else len(raw_lines)
        scope = normalized_fleet_scope(raw_lines[start])
        if scope in blocks:
            raise VerificationFailure(
                f"duplicate accepted certificate scope in {locus}: {scope}"
            )
        arithmetic = tuple(arithmetic_output_lines("\n".join(raw_lines[start + 1 : end])))
        if not arithmetic:
            raise VerificationFailure(
                f"empty accepted arithmetic block in {locus}: {scope}"
            )
        blocks[scope] = arithmetic
    return blocks


def verify_frontier_certificate_scope_coverage(root: Path, ledger) -> tuple[int, int]:
    """Require the cited frontier blocks to equal, not merely contain, the ledger.

    Per-partition replay checks below prove that every expected task regenerates
    its accepted arithmetic block.  This global check closes the converse:
    every arithmetic block in each cited H25--H29 certificate must be replayed
    by exactly one of the 71 tasks.  Thus an extra, stale, or unpartitioned
    certificate block cannot pass unnoticed.
    """

    tasks = ledger.make_tasks(8)
    if len(tasks) != 71 or len({task.name for task in tasks}) != 71:
        raise VerificationFailure("frontier scope ledger is not 71 unique tasks")
    grouped: dict[str, list] = {}
    for task in tasks:
        grouped.setdefault(frontier_certificate_name(task.name), []).append(task)
    if len(grouped) != 7:
        raise VerificationFailure(
            f"frontier certificate ledger changed: files={len(grouped)}"
        )

    accepted_manifest = parse_manifest(
        root / "certificates/post_m29/post_m29_accepted_manifest.sha256"
    )
    seen_scopes: set[str] = set()
    block_count = 0
    for certificate_name, certificate_tasks in sorted(grouped.items()):
        certificate = require_file(root / "certificates/post_m29" / certificate_name)
        if accepted_manifest.get(certificate_name) != digest(certificate):
            raise VerificationFailure(
                f"accepted frontier certificate identity mismatch: {certificate_name}"
            )
        blocks = accepted_arithmetic_blocks(
            certificate.read_text(errors="strict"), certificate
        )
        expected_scopes = {
            normalized_fleet_scope(expected_fleet_scope_line(task, 8))
            for task in certificate_tasks
        }
        actual_scopes = set(blocks)
        if actual_scopes != expected_scopes:
            raise VerificationFailure(
                f"accepted frontier scope coverage mismatch for {certificate_name}: "
                f"missing={sorted(expected_scopes - actual_scopes)}, "
                f"extra={sorted(actual_scopes - expected_scopes)}"
            )
        overlap = seen_scopes & actual_scopes
        if overlap:
            raise VerificationFailure(
                f"frontier scopes occur in multiple certificates: {sorted(overlap)}"
            )
        seen_scopes.update(actual_scopes)
        block_count += len(actual_scopes)
    if block_count != 71 or len(seen_scopes) != 71:
        raise VerificationFailure(
            "frontier certificate block coverage changed: "
            f"blocks={block_count}, unique_scopes={len(seen_scopes)}"
        )
    return len(grouped), block_count


def resolve_stage_certificate(root: Path, certificate: str) -> Path:
    relative = Path(certificate)
    return require_file(
        (root / "certificates/post_m29" / relative)
        if len(relative.parts) == 1
        else (root / "certificates" / relative)
    )


def normalized_direct_certificate_line(stage: StageSpec, line: str) -> str:
    """Normalize only reproducibility metadata, never arithmetic fields."""

    line = re.sub(r"omp_max_threads=\d+", "omp_max_threads=<THREADS>", line)
    line = re.sub(r"OpenMP threads=\d+", "OpenMP threads=<THREADS>", line)
    if stage.number in {44, 47, 48}:
        line = re.sub(r"onset_log=\S+", "onset_log=<PATH>", line)
    if stage.number in {38, 166, 167}:
        for prefix in (
            "exact_D_adjoint_moment_source_log=",
            "source_log=",
            "stable_moments=",
        ):
            if line.startswith(prefix):
                line = prefix + Path(line[len(prefix) :]).name
    stripped = line.lstrip()
    leading = line[: len(line) - len(stripped)]
    if stage.number == 183 and "ginibre_q3/" in stripped:
        stripped = stripped[stripped.index("ginibre_q3/") :]
        line = leading + stripped
    if stage.number in {45, 46} and "certificates/" in stripped:
        stripped = stripped[stripped.index("certificates/") :]
        line = leading + stripped
    return line


def direct_certificate_payload(stage: StageSpec, text: str) -> tuple[str, ...]:
    """Return the deterministic, load-bearing output of one direct replay."""

    ignored_progress: dict[int, tuple[str, ...]] = {
        45: ("row_done ",),
        46: ("factorial_progress ", "row_done "),
        47: ("progress completed=",),
        48: ("progress completed=",),
    }
    ignored = ignored_progress.get(stage.number, ())
    return tuple(
        normalized_direct_certificate_line(stage, line)
        for line in text.splitlines()
        if line
        and not line.startswith(("COMMAND: ", "ENV: ", "__"))
        and not line.startswith(ignored)
    )


def clean_aggregate_frontier_block(lines: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    """Remove aggregate-file separators and clean-run command metadata."""

    return tuple(
        line
        for line in lines
        if not line.startswith(("COMMAND: ", "ENV: ", "===== ", "=== "))
    )


def verify_nonfleet_certificate_projections(
    root: Path,
    specs: list[StageSpec],
    current_prefix_stages: Path,
    legacy_stages: Path,
    classical_stages: Path,
) -> int:
    """Bind every frozen direct/H23/H24 replay to its cited certificate data."""

    direct = [
        stage
        for stage in specs
        if stage.certificates
        and not 65 <= stage.number <= 154
    ]
    if len(direct) != 34 or any(len(stage.certificates) != 1 for stage in direct):
        raise VerificationFailure(
            f"direct certificate-projection ledger changed: {len(direct)}"
        )

    def stage_directory(stage: StageSpec) -> Path:
        if stage.number <= 49:
            return current_prefix_stages
        if stage.number <= 64:
            return legacy_stages
        if 155 <= stage.number <= 184:
            return classical_stages
        raise VerificationFailure(
            f"unrouted direct certificate stage: {stage.number} {stage.name}"
        )

    projections = 0
    for stage in direct:
        replay = find_stage_log(stage_directory(stage), stage.name)
        certificate = resolve_stage_certificate(root, stage.certificates[0])
        actual = direct_certificate_payload(stage, replay.read_text(errors="strict"))
        expected = direct_certificate_payload(
            stage, certificate.read_text(errors="strict")
        )
        if actual != expected:
            first_difference = next(
                (
                    index
                    for index in range(max(len(actual), len(expected)))
                    if index >= len(actual)
                    or index >= len(expected)
                    or actual[index] != expected[index]
                ),
                None,
            )
            raise VerificationFailure(
                f"direct certificate payload mismatch for stage {stage.number} "
                f"{stage.name}: expected_lines={len(expected)}, "
                f"actual_lines={len(actual)}, first_difference={first_difference}"
            )
        projections += 1

    by_number = {stage.number: stage for stage in specs}
    aggregate_specs = (
        (65, 73, 73),
        (74, 83, 83),
    )
    for low, high, certificate_stage_number in aggregate_specs:
        certificate_stage = by_number[certificate_stage_number]
        if len(certificate_stage.certificates) != 1:
            raise VerificationFailure(
                f"aggregate frontier certificate ledger changed at stage "
                f"{certificate_stage_number}"
            )
        certificate = resolve_stage_certificate(
            root, certificate_stage.certificates[0]
        )
        blocks = accepted_arithmetic_blocks(
            certificate.read_text(errors="strict"), certificate
        )
        for number in range(low, high + 1):
            stage = by_number[number]
            replay = find_stage_log(legacy_stages, stage.name)
            text = replay.read_text(errors="strict")
            scope_lines = [line for line in text.splitlines() if line.startswith("BC_")]
            if len(scope_lines) != 1:
                raise VerificationFailure(
                    f"expected one aggregate frontier scope in {replay}, "
                    f"found {len(scope_lines)}"
                )
            scope = normalized_fleet_scope(scope_lines[0])
            if scope not in blocks:
                raise VerificationFailure(
                    f"accepted aggregate frontier block missing for {stage.name}: {scope}"
                )
            actual = clean_aggregate_frontier_block(arithmetic_output_lines(text))
            expected = clean_aggregate_frontier_block(blocks[scope])
            if actual != expected:
                raise VerificationFailure(
                    f"aggregate frontier payload mismatch for stage {number} "
                    f"{stage.name}: expected_lines={len(expected)}, "
                    f"actual_lines={len(actual)}"
                )
            projections += 1

    if projections != 53:
        raise VerificationFailure(
            f"nonfleet certificate-projection count changed: {projections}"
        )
    return projections


def output_q_values(lines: list[str], prefix: str, locus: Path) -> list[int]:
    values: list[int] = []
    pattern = re.compile(rf"^{re.escape(prefix)} q=(\d+)\b")
    for line in lines:
        match = pattern.match(line)
        if match is not None:
            values.append(int(match.group(1)))
    if len(values) != len(set(values)):
        raise VerificationFailure(f"duplicate {prefix} rank in {locus}: {values}")
    return values


def verify_fleet_arithmetic_projection(
    task,
    text: str,
    log: Path,
    expected_lines: tuple[str, ...],
) -> int:
    """Check complete ordered equality with the paper-cited task block."""

    lines = arithmetic_output_lines(text)
    if not lines:
        raise VerificationFailure(f"empty arithmetic output in {log}")
    if tuple(lines) != expected_lines:
        first_difference = next(
            (
                index
                for index in range(max(len(lines), len(expected_lines)))
                if index >= len(lines)
                or index >= len(expected_lines)
                or lines[index] != expected_lines[index]
            ),
            None,
        )
        raise VerificationFailure(
            "fleet arithmetic is not the complete accepted task projection in "
            f"{log}: expected_lines={len(expected_lines)}, actual_lines={len(lines)}, "
            f"first_difference={first_difference}"
        )

    environment = dict(task.environment)
    if environment.get("RUN_B") == "1":
        expected_q = list(
            range(int(environment["B_Q_MIN"]), int(environment["B_Q_MAX"]) + 1)
        )
        for prefix in ("B_CASE", "B_BAD_COUNT", "B_STABLE", "B_MARGIN"):
            if output_q_values(lines, prefix, log) != expected_q:
                raise VerificationFailure(
                    f"incomplete {prefix} rank coverage in {log}"
                )
    elif environment.get("RUN_C") == "1":
        scope_match = re.search(
            r"c_fpf_cases=(\d+) c_arithmetic_cases=(\d+)", task.required[0]
        )
        if scope_match is None:
            raise VerificationFailure(f"cannot parse C scope ledger for {task.name}")
        fpf_count, arithmetic_count = map(int, scope_match.groups())
        fpf_q = list(range(21, 21 + fpf_count))
        arithmetic_q = list(
            range(21 + fpf_count, 21 + fpf_count + arithmetic_count)
        )
        for prefix in ("C_FPF_BOUND", "C_FPF_STABLE", "C_FPF_MARGIN"):
            if output_q_values(lines, prefix, log) != fpf_q:
                raise VerificationFailure(
                    f"incomplete {prefix} rank coverage in {log}"
                )
        for prefix in ("C_PAUSE_POLYNOMIAL", "C_BOUND", "C_STABLE", "C_MARGIN"):
            if output_q_values(lines, prefix, log) != arithmetic_q:
                raise VerificationFailure(
                    f"incomplete {prefix} rank coverage in {log}"
                )
    elif task.name == "twentyeighth_b_absorption":
        if output_q_values(lines, "B_VARIABLE_ABSORPTION", log) != list(range(42, 50)):
            raise VerificationFailure(f"incomplete H28 absorption coverage in {log}")
    elif task.name == "twentyninth_c_frontier":
        fpf_q = list(range(21, 35))
        arithmetic_q = list(range(35, 74))
        for prefix in ("C_FPF_BOUND", "C_FPF_STABLE", "C_FPF_MARGIN"):
            if output_q_values(lines, prefix, log) != fpf_q:
                raise VerificationFailure(f"incomplete H29 {prefix} coverage in {log}")
        for prefix in ("C_PAUSE_POLYNOMIAL", "C_BOUND", "C_STABLE", "C_MARGIN"):
            if output_q_values(lines, prefix, log) != arithmetic_q:
                raise VerificationFailure(f"incomplete H29 {prefix} coverage in {log}")
    elif task.name == "twentyninth_b_absorption":
        if output_q_values(lines, "B_VARIABLE_ABSORPTION", log) != list(range(43, 53)):
            raise VerificationFailure(f"incomplete H29 absorption coverage in {log}")
    else:
        raise VerificationFailure(f"unrecognized arithmetic fleet task {task.name}")
    return len(lines)


def verify_fleet_partition(
    root: Path,
    ledger,
    partition_wrapper,
    partition: str,
    workers: int,
    threads: int,
    run_id: str,
    host_label: str,
    build_log: Path,
    top_log: Path,
    results: Path,
    manifest_path: Path,
    expected_stages: list[StageSpec],
    original_full_queue: bool,
) -> list[str]:
    tasks = frontier_partition_tasks(ledger, threads, partition)
    if not original_full_queue:
        wrapper_tasks = partition_wrapper.partition_tasks(ledger, threads, partition)
        if [task.name for task in wrapper_tasks] != [task.name for task in tasks]:
            raise VerificationFailure(
                f"partition wrapper selection mismatch for {partition}"
            )
    mapped = {fleet_stage_name(task.name) for task in tasks}
    expected = {stage.name for stage in expected_stages}
    if mapped != expected or len(mapped) != len(tasks):
        raise VerificationFailure(
            f"fleet/clean mismatch for {partition}: "
            f"missing={sorted(expected-mapped)}, extra={sorted(mapped-expected)}"
        )
    verify_fleet_task_semantics(tasks, expected_stages, threads, partition)

    build_text = require_file(build_log).read_text(errors="strict")
    require_line_once(build_text, "__BUILD_EXIT_STATUS=0", build_log)
    require_line_once(build_text, "__SOURCE_BINARY_IDENTITIES_BEGIN", build_log)
    require_line_once(build_text, "__SOURCE_BINARY_IDENTITIES_END", build_log)
    if not original_full_queue:
        staged_wrapper = require_file(
            results.parent / "run_full_q3_frontier_partition.py"
        )
        staged_ledger = require_file(
            results.parent / "run_full_q3_frontier_fleet_queue.py"
        )
        if digest(staged_wrapper) != FLEET_PARTITION_SCRIPT_SHA256:
            raise VerificationFailure(
                f"staged fleet partition wrapper mismatch for {partition}"
            )
        if digest(staged_ledger) != FLEET_SCRIPT_SHA256:
            raise VerificationFailure(
                f"staged fleet task ledger mismatch for {partition}"
            )
        require_hash_record(
            build_text,
            FLEET_PARTITION_SCRIPT_SHA256,
            "run_full_q3_frontier_partition.py",
            build_log,
        )
    else:
        staged_ledger = require_file(
            results.parent / "run_full_q3_frontier_fleet_queue.py"
        )
        if digest(staged_ledger) != FLEET_SCRIPT_SHA256:
            raise VerificationFailure("staged machine-A fleet ledger mismatch")
        require_hash_record(
            build_text,
            MACHINE_A_BUILD_LEDGER_SHA256,
            "run_full_q3_frontier_fleet_queue.py",
            build_log,
        )
    staged_character = results.parent / "source/character_ring_iter"
    unique_binaries = {task.binary: task.source for task in tasks}
    for binary_name, source_name in sorted(unique_binaries.items()):
        source = require_file(staged_character / source_name)
        current_source = require_file(root / "character_ring_iter" / source_name)
        binary = require_file(staged_character / binary_name)
        if digest(source) != digest(current_source):
            raise VerificationFailure(f"fleet staged source mismatch: {source_name}")
        require_hash_record(
            build_text,
            digest(source),
            f"source/character_ring_iter/{source_name}",
            build_log,
        )
        require_hash_record(
            build_text,
            digest(binary),
            f"source/character_ring_iter/{binary_name}",
            build_log,
        )
        compile_pattern = re.compile(
            rf"^g\+\+ .*\b{re.escape(source_name)}\b .* -o {re.escape(binary_name)}$",
            re.MULTILINE,
        )
        if compile_pattern.search(build_text) is None:
            raise VerificationFailure(
                f"missing fleet compile command for {binary_name} in {build_log}"
            )

    top = verify_top_log(top_log)
    if " status=FAIL " in top or "VERIFICATION: FAIL" in top:
        raise VerificationFailure(f"fleet failure marker in {top_log}")
    if original_full_queue:
        require_line_once(
            top,
            f"FULL_Q3_FLEET start run_id={run_id} tasks=71 workers=4 "
            f"threads_per_task={threads} script_sha256={FLEET_SCRIPT_SHA256}",
            top_log,
        )
        pass_prefix = "FULL_Q3_FLEET"
    else:
        require_line_once(
            top,
            f"FULL_Q3_FLEET_PARTITION start run_id={run_id} "
            f"partition={partition} tasks={len(tasks)} workers={workers} "
            f"threads_per_task={threads} "
            f"wrapper_sha256={FLEET_PARTITION_SCRIPT_SHA256} "
            f"ledger_sha256={FLEET_SCRIPT_SHA256}",
            top_log,
        )
        require_line_once(
            top,
            f"FULL_Q3_FLEET_PARTITION partition={partition} "
            f"tasks_passed={len(tasks)} tasks_expected={len(tasks)} failures=0",
            top_log,
        )
        require_line_once(
            top, "FULL_Q3_FLEET_PARTITION VERIFICATION: ALL PASS", top_log
        )
        require_line_once(top, "__ORCHESTRATOR_EXIT_STATUS=0", top_log)
        pass_prefix = "FULL_Q3_FLEET_PARTITION"

    manifest = parse_manifest(manifest_path)
    expected_log_names = {f"{task.name}.log" for task in tasks}
    if len(manifest) != len(tasks) or set(manifest) != expected_log_names:
        raise VerificationFailure(
            f"fleet manifest does not exactly cover {partition}: {manifest_path}"
        )
    accepted_manifest = parse_manifest(
        root / "certificates/post_m29/post_m29_accepted_manifest.sha256"
    )
    accepted_data: dict[str, dict[str, tuple[str, ...]]] = {}
    for task in tasks:
        certificate_name = frontier_certificate_name(task.name)
        if certificate_name in accepted_data:
            continue
        certificate = require_file(root / "certificates/post_m29" / certificate_name)
        if accepted_manifest.get(certificate_name) != digest(certificate):
            raise VerificationFailure(
                f"accepted frontier certificate identity mismatch: {certificate_name}"
            )
        accepted_data[certificate_name] = accepted_arithmetic_blocks(
            certificate.read_text(errors="strict"), certificate
        )
    hashes: list[str] = []
    for task in tasks:
        log = require_file(results / f"{task.name}.log")
        text = log.read_text(errors="strict")
        lines = text.splitlines()
        require_line_once(text, f"__RUN_ID={run_id}", log)
        require_line_once(text, f"__HOST_LABEL={host_label}", log)
        require_line_once(text, f"__TASK={task.name}", log)
        require_line_once(text, f"__COMMAND=./{task.binary}", log)
        status_lines = [line for line in lines if line.startswith("__EXIT_STATUS=")]
        if status_lines != ["__EXIT_STATUS=0", "__EXIT_STATUS=0"]:
            raise VerificationFailure(
                f"expected successful program and orchestrator statuses in {log}: "
                f"{status_lines}"
            )
        require_line_once(text, "SUMMARY failures=0", log)
        require_line_once(text, expected_fleet_scope_line(task, threads), log)
        start_matches = re.findall(
            r"^__START_UTC=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)$",
            text,
            re.MULTILINE,
        )
        end_matches = re.findall(
            r"^__END_UTC=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)$",
            text,
            re.MULTILINE,
        )
        elapsed_matches = re.findall(
            r"^__ELAPSED_SECONDS=([0-9]+(?:\.[0-9]+)?)$",
            text,
            re.MULTILINE,
        )
        if (
            len(start_matches) != 1
            or len(end_matches) != 1
            or len(elapsed_matches) != 1
            or end_matches[0] < start_matches[0]
            or float(elapsed_matches[0]) <= 0.0
            or len(lines) < 3
            or lines[-3] != "__EXIT_STATUS=0"
            or lines[-2] != f"__ELAPSED_SECONDS={elapsed_matches[0]}"
            or lines[-1] != f"__END_UTC={end_matches[0]}"
        ):
            raise VerificationFailure(f"invalid timing ledger in {log}")
        expected_environment = " ".join(
            f"{key}={value}" for key, value in task.environment
        )
        require_line_once(text, f"__ENV={expected_environment}", log)
        source = root / "character_ring_iter" / task.source
        source_hash = digest(require_file(source))
        require_line_once(text, f"__SOURCE_SHA256={source_hash}  {task.source}", log)
        binary_matches = re.findall(
            rf"^__BINARY_SHA256=([0-9a-f]{{64}})  {re.escape(task.binary)}$",
            text,
            re.MULTILINE,
        )
        if len(binary_matches) != 1:
            raise VerificationFailure(f"missing binary identity in {log}")
        staged_binary = staged_character / task.binary
        if digest(require_file(staged_binary)) != binary_matches[0]:
            raise VerificationFailure(
                f"staged binary changed after fleet run: {task.binary}"
            )
        actual_hash = digest(log)
        if manifest[f"{task.name}.log"] != actual_hash:
            raise VerificationFailure(f"fleet log hash mismatch: {log}")
        pass_pattern = re.compile(
            rf"^{re.escape(pass_prefix)} task={re.escape(task.name)} status=PASS "
            rf"elapsed=([0-9]+(?:\.[0-9]+)?)s sha256=([0-9a-f]{{64}})$",
            re.MULTILINE,
        )
        pass_matches = pass_pattern.findall(top)
        if (
            len(pass_matches) != 1
            or float(pass_matches[0][0]) <= 0.0
            or pass_matches[0][1] != actual_hash
        ):
            raise VerificationFailure(
                f"fleet PASS record does not bind the accepted log for {task.name}"
            )
        certificate_name = frontier_certificate_name(task.name)
        accepted_scope = normalized_fleet_scope(
            expected_fleet_scope_line(task, threads)
        )
        accepted_blocks = accepted_data[certificate_name]
        if accepted_scope not in accepted_blocks:
            raise VerificationFailure(
                f"no unique accepted arithmetic block for {task.name}: "
                f"{accepted_scope}"
            )
        verify_fleet_arithmetic_projection(
            task,
            text,
            log,
            accepted_blocks[accepted_scope],
        )
        hashes.append(actual_hash)
    return hashes


def verify_wrapper(
    label: str,
    top_log: Path,
    stage_dir: Path,
    stages: list[StageSpec],
    terminal_count: int,
) -> list[str]:
    text = verify_top_log(top_log)
    require_line_once(text, f"__RUN_ID={label}", top_log)
    require_line_once(text, "__HOST=optimus", top_log)
    require_line_once(text, f"__WRAPPER_SHA256={POST_WRAPPER_SHA256}", top_log)
    require_line_once(text, f"__CLEAN_DRIVER_SHA256={CURRENT_DRIVER_SHA256}", top_log)
    require_line_once(
        text,
        f"__SOURCE_MANIFEST_SHA256={CURRENT_REPLAY_MANIFEST_SHA256}",
        top_log,
    )
    prefix = "FULL_Q3_CLASSICAL_POST" if label == "fullq3classicalpost0001" else "FULL_Q3_EXCEPTIONAL"
    require_line_once(
        text,
        f"{prefix} tasks_passed={terminal_count} tasks_expected={terminal_count} failures=0",
        top_log,
    )
    require_line_once(text, f"{prefix} VERIFICATION: ALL PASS", top_log)
    require_line_once(text, "__EXIT_STATUS=0", top_log)
    return [verify_stage_log(stage, stage_dir, text, top_log) for stage in stages]


def verify_post_build(
    build_log: Path,
    root: Path,
    rebuilt: tuple[tuple[str, str], ...],
    hashed_only: tuple[str, ...] = (),
) -> None:
    text = require_file(build_log).read_text(errors="strict")
    require_line_once(text, "__EXIT_STATUS=0", build_log)
    character = root / "character_ring_iter"
    for source_name, binary_name in rebuilt:
        source = require_file(character / source_name)
        binary = require_file(character / binary_name)
        require_line_once(text, f"{digest(source)}  {source_name}", build_log)
        require_line_once(text, f"{digest(binary)}  {binary_name}", build_log)
        compile_pattern = re.compile(
            rf"^g\+\+ .*\b{re.escape(source_name)}\b .* -o {re.escape(binary_name)}$",
            re.MULTILINE,
        )
        if compile_pattern.search(text) is None:
            raise VerificationFailure(
                f"missing explicit compile command for {binary_name} in {build_log}"
            )
    for name in hashed_only:
        path = require_file(character / name)
        require_line_once(text, f"{digest(path)}  {name}", build_log)


def replayed_certificate_set(stages: list[StageSpec]) -> set[str]:
    replayed: set[str] = set()
    for stage in stages:
        for certificate in stage.certificates:
            relative = Path(certificate)
            normalized = (
                Path("certificates/post_m29") / relative
                if len(relative.parts) == 1
                else Path("certificates") / relative
            )
            replayed.add(normalized.as_posix())
    return replayed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--current-prefix-log", type=Path, required=True)
    parser.add_argument("--current-prefix-stages", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--legacy-log", type=Path, required=True)
    parser.add_argument("--legacy-stages", type=Path, required=True)
    parser.add_argument("--fleet-ledger-script", type=Path, required=True)
    parser.add_argument("--fleet-partition-script", type=Path, required=True)
    for label in ("a", "b", "c"):
        parser.add_argument(f"--fleet-{label}-build-log", type=Path, required=True)
        parser.add_argument(f"--fleet-{label}-log", type=Path, required=True)
        parser.add_argument(f"--fleet-{label}-results", type=Path, required=True)
        parser.add_argument(f"--fleet-{label}-manifest", type=Path, required=True)
    parser.add_argument("--classical-build-log", type=Path, required=True)
    parser.add_argument("--classical-log", type=Path, required=True)
    parser.add_argument("--classical-stages", type=Path, required=True)
    parser.add_argument("--exceptional-build-log", type=Path, required=True)
    parser.add_argument("--exceptional-log", type=Path, required=True)
    parser.add_argument("--exceptional-stages", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    legacy_root = args.legacy_root.resolve()
    current_manifest, legacy_manifest = verify_current_source_identity(root, legacy_root)
    all_pd_obstruction = verify_part_iii_scope(root)
    verify_header_generator(root)
    frontier_kernel_partitions = verify_frontier_kernel(root)
    driver = load_module(root / "clean_room_replay.py", "full_q3_clean_driver_audit")
    build_dependencies = verify_replay_build_dependencies(root, driver)
    specs = collect_stage_specs(driver, root, 32)
    prefix_stages = [stage for stage in specs if 2 <= stage.number <= 49]
    legacy_stages = [stage for stage in specs if 50 <= stage.number <= 83]
    fleet_stages = [stage for stage in specs if 84 <= stage.number <= 154]
    classical_stages = [stage for stage in specs if 155 <= stage.number <= 184]
    exceptional_stages_32 = [stage for stage in specs if 185 <= stage.number <= 200]
    if tuple(map(len, (prefix_stages, legacy_stages, fleet_stages, classical_stages, exceptional_stages_32))) != (48, 34, 71, 30, 16):
        raise VerificationFailure("distributed replay partition cardinality mismatch")

    current_top = verify_top_log(args.current_prefix_log, PREFLIGHT_MARKERS)
    require_replay_pass_line(
        current_top, "build_cpp_replayers", args.current_prefix_log
    )
    build_log = find_stage_log(args.current_prefix_stages, "build_cpp_replayers")
    build_lines = build_log.read_text(errors="strict").splitlines()
    if not build_lines or normalized_command(build_lines[0].removeprefix("COMMAND: ")) != (
        "make", "-j", "8", "replay-build"
    ):
        raise VerificationFailure(f"current replay build command mismatch in {build_log}")
    all_hashes = [digest(build_log)]
    all_hashes.extend(
        verify_stage_log(stage, args.current_prefix_stages, current_top, args.current_prefix_log)
        for stage in prefix_stages
    )

    verify_legacy_sources(root, legacy_root, current_manifest, legacy_manifest, legacy_stages)
    legacy_top = verify_top_log(args.legacy_log)
    all_hashes.extend(
        verify_stage_log(
            stage,
            args.legacy_stages,
            legacy_top,
            args.legacy_log,
            require_exit=True,
        )
        for stage in legacy_stages
    )
    if digest(require_file(args.fleet_ledger_script)) != FLEET_SCRIPT_SHA256:
        raise VerificationFailure("fleet task-ledger identity mismatch")
    if digest(require_file(args.fleet_partition_script)) != FLEET_PARTITION_SCRIPT_SHA256:
        raise VerificationFailure("fleet partition-wrapper identity mismatch")
    fleet_ledger = load_module(
        args.fleet_ledger_script,
        "full_q3_fleet_task_ledger_audit",
    )
    fleet_partition_wrapper = load_module(
        args.fleet_partition_script,
        "full_q3_fleet_partition_wrapper_audit",
    )
    frontier_certificate_files, frontier_certificate_blocks = (
        verify_frontier_certificate_scope_coverage(root, fleet_ledger)
    )
    fleet_stage_by_name = {stage.name: stage for stage in fleet_stages}

    def partition_expected(partition: str, threads: int) -> list[StageSpec]:
        tasks = frontier_partition_tasks(fleet_ledger, threads, partition)
        return [fleet_stage_by_name[fleet_stage_name(task.name)] for task in tasks]

    all_hashes.extend(
        verify_fleet_partition(
            root,
            fleet_ledger,
            fleet_partition_wrapper,
            "machine_a_low",
            4,
            8,
            "fullq3cleanfleet0002",
            "optimus",
            args.fleet_a_build_log,
            args.fleet_a_log,
            args.fleet_a_results,
            args.fleet_a_manifest,
            partition_expected("machine_a_low", 8),
            True,
        )
    )
    all_hashes.extend(
        verify_fleet_partition(
            root,
            fleet_ledger,
            fleet_partition_wrapper,
            "machine_b_mid",
            1,
            24,
            "fullq3fleetbmid0001",
            "machine_b",
            args.fleet_b_build_log,
            args.fleet_b_log,
            args.fleet_b_results,
            args.fleet_b_manifest,
            partition_expected("machine_b_mid", 24),
            False,
        )
    )
    all_hashes.extend(
        verify_fleet_partition(
            root,
            fleet_ledger,
            fleet_partition_wrapper,
            "machine_c_heavy",
            1,
            8,
            "fullq3fleetcheavy0004",
            "machine_c",
            args.fleet_c_build_log,
            args.fleet_c_log,
            args.fleet_c_results,
            args.fleet_c_manifest,
            partition_expected("machine_c_heavy", 8),
            False,
        )
    )
    verify_post_build(
        args.classical_build_log,
        root,
        (
            ("post_m29_bc_layered_mgf_mpfr.cpp", "post_m29_bc_layered_mgf_mpfr"),
            ("verify_d10_exact_tail_gmp.cpp", "verify_d10_exact_tail_gmp"),
            (
                "post_m29_d_stable_window_interval_gmp.cpp",
                "post_m29_d_stable_window_interval_gmp",
            ),
        ),
        hashed_only=(
            "verify_character_ring_moment_sources_gmp.cpp",
            "verify_character_ring_moment_sources_gmp",
            "lie_data.h",
        ),
    )
    all_hashes.extend(
        verify_wrapper(
            "fullq3classicalpost0001",
            args.classical_log,
            args.classical_stages,
            classical_stages,
            30,
        )
    )
    exceptional_specs_8 = collect_stage_specs(driver, root, 8)
    exceptional_stages_8 = [
        stage for stage in exceptional_specs_8 if 185 <= stage.number <= 200
    ]
    if {stage.name for stage in exceptional_stages_8} != {
        stage.name for stage in exceptional_stages_32
    }:
        raise VerificationFailure("exceptional thread-variant stage ledger mismatch")
    verify_post_build(
        args.exceptional_build_log,
        root,
        (
            (
                "verify_character_ring_moment_sources_gmp.cpp",
                "verify_character_ring_moment_sources_gmp",
            ),
            (
                "dunkl_exceptional_coefficients_gmp.cpp",
                "dunkl_exceptional_coefficients_gmp",
            ),
            (
                "direct_chain_rect_g2_f4_certificate.cpp",
                "direct_chain_rect_g2_f4_certificate",
            ),
            (
                "direct_chain_rect_e6_e7_rank_certificate.cpp",
                "direct_chain_rect_e6_e7_rank_certificate",
            ),
            ("verify_e8_negative_bounds_gmp.cpp", "verify_e8_negative_bounds_gmp"),
            (
                "direct_chain_rect_e8_parallel_certificate.cpp",
                "direct_chain_rect_e8_parallel_certificate",
            ),
        ),
        hashed_only=("gen_header.py", "lie_data.py", "lie_data.h"),
    )
    all_hashes.extend(
        verify_wrapper(
            "fullq3exceptional0001",
            args.exceptional_log,
            args.exceptional_stages,
            exceptional_stages_8,
            16,
        )
    )

    nonfleet_certificate_projections = verify_nonfleet_certificate_projections(
        root,
        specs,
        args.current_prefix_stages,
        args.legacy_stages,
        args.classical_stages,
    )
    certificate_data_projections = (
        nonfleet_certificate_projections + len(fleet_stages)
    )
    if certificate_data_projections != 124:
        raise VerificationFailure(
            "certificate data-projection count changed: "
            f"{certificate_data_projections}"
        )
    if len(all_hashes) != 200 or len(set(all_hashes)) != 200:
        raise VerificationFailure(
            f"composed log ledger is not 200 distinct stages: {len(all_hashes)}"
        )
    replayed = replayed_certificate_set(specs)
    reachable = driver.audit_replayed_certificate_coverage(root, replayed)
    if reachable != 43 or len(replayed) != 43:
        raise VerificationFailure(
            f"main-theorem certificate ledger changed: reachable={reachable}, replayed={len(replayed)}"
        )
    manifest_count, manifest_entries = driver.verify_manifests(root)
    print("__RUN_ID=fullq3distributed0001")
    print(f"__HOST={os.uname().nodename}")
    print(f"__VERIFIER_SHA256={digest(Path(__file__).resolve())}")
    print(f"__CURRENT_DRIVER_SHA256={CURRENT_DRIVER_SHA256}")
    print(f"__CURRENT_REPLAY_MANIFEST_SHA256={CURRENT_REPLAY_MANIFEST_SHA256}")
    print(
        "FULL_Q3_DISTRIBUTED source_manifest_inputs=85 unchanged_inputs=77 "
        "changed_inputs=8"
    )
    print(
        "FULL_Q3_DISTRIBUTED theorem_scope=adjoint_generated "
        f"all_positive_definite_obstruction={all_pd_obstruction}"
    )
    print(
        "FULL_Q3_DISTRIBUTED stages=200 prefix=49 legacy_h8_h24=34 "
        "fleet_a_h25_h26=24 fleet_b_h27_low=7 fleet_c_h27_h29_heavy=40 "
        "classical_post=30 exceptional=16 failures=0"
    )
    print(
        f"FULL_Q3_DISTRIBUTED reachable_certificates={reachable} "
        f"build_dependencies={build_dependencies} manifests={manifest_count} "
        f"manifest_entries={manifest_entries} "
        f"frontier_kernel_partitions={frontier_kernel_partitions} "
        f"frontier_certificate_files={frontier_certificate_files} "
        f"frontier_certificate_blocks={frontier_certificate_blocks} "
        f"certificate_data_projections={certificate_data_projections}"
    )
    print(
        "FULL_Q3_DISTRIBUTED machine_a_build_ledger="
        f"{MACHINE_A_BUILD_LEDGER_SHA256} "
        f"machine_a_run_ledger={FLEET_SCRIPT_SHA256} "
        "build_ledger_non_compilation_input=1 run_ledger_authenticated=1"
    )
    print("FULL_Q3_DISTRIBUTED VERIFICATION: ALL PASS")
    print("__EXIT_STATUS=0")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FULL_Q3_DISTRIBUTED VERIFICATION: FAIL ({error})", file=sys.stderr)
        print("__EXIT_STATUS=1", file=sys.stderr)
        raise SystemExit(1)
