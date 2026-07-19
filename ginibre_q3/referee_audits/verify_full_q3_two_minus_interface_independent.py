#!/usr/bin/env python3
"""Fail-closed cross-document audit of the Part III two-minus import.

Part III uses the Parts I--II theorem as its complete ``a=1`` row.  This
author self-audit parses both theorem statements and both defining Haar
integrals, substitutes ``a=1`` in the Part III signature, and checks that the
group scope, exponent range, normalization, theorem number, bibliography
entry, and manifest bindings agree exactly.  It is structural and does not
replace the arithmetic replay proving the imported theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


DEFAULT_ROOT = Path(__file__).resolve().parents[1]
ALLOWED_COMPANION_HASHES = {
    "6b2f273aa082f105df4bf60f8ba7047c62cafac88a3bc2c0a5cffdca95607d27",
    "d2fb8944ec081f47a1bfcf558bb65adb01936b0ac1d0dbbb1f561d63d7197010",
    "eb91422f3a32840d4a2649f4b60f09ccdd4599da42ad02c7c31276a9c7c535f4",
    "8b58f9015dc4738c0e48e3d340f512a74b54a6832f83a3ae378b910d75e6060f",
    "bd2a53c86a00caa7f7e4d9a7d61f6a02f2059fadde313e080b3360d49da03f93",
    "dba5091d17cd85aeeee73deaafda86df1fe05e395625489a391201898910e325",
    "14881fb1fdc1a1e4d0c80faaaf710e880f768b4e1b820624298345e60fcf1cbd",
}
EXPECTED_PART_THREE_HASH = (
    "9362671cfb5fbd4750a11a072b13f07a64073362e43e84da198ef73f12bd9439"
)
RESULT_ENVS = ("theorem", "proposition", "lemma", "corollary")
RESULT_PATTERN = re.compile(
    r"\\begin\{(" + "|".join(RESULT_ENVS) + r")\}(.*?)\\end\{\1\}",
    re.DOTALL,
)
DEFINITION_PATTERN = re.compile(
    r"\\begin\{definition\}(.*?)\\end\{definition\}",
    re.DOTALL,
)
INTEGRAL_PATTERN = re.compile(
    r"\\int_\{G\\timesG\}"
    r"\(\\chiad\(g\)-\\chiad\(h\)\)\^"
    r"(?:\{(?P<difference_braced>[^{}]+)\}|(?P<difference_plain>[A-Za-z0-9]+))"
    r"\(\\chiad\(g\)\+\\chiad\(h\)\)\^"
    r"(?:\{(?P<plus_braced>[^{}]+)\}|(?P<plus_plain>[A-Za-z0-9]+))"
    r"(?:\\,)?d(?P<measure_g>\\Haar|\\mu_G)\(g\)"
    r"(?:\\,)?d(?P<measure_h>\\Haar|\\mu_G)\(h\)"
)


class InterfaceFailure(RuntimeError):
    """A required cross-document invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise InterfaceFailure(message)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block_with_title(text: str, pattern: re.Pattern[str], title: str) -> str:
    matches = [match.group(0) for match in pattern.finditer(text) if title in match.group(1)]
    require(len(matches) == 1, f"expected one block titled {title}, found {len(matches)}")
    return matches[0]


def result_with_label(text: str, label: str) -> str:
    needle = rf"\label{{{label}}}"
    matches = [match.group(0) for match in RESULT_PATTERN.finditer(text) if needle in match.group(2)]
    require(len(matches) == 1, f"expected one result labelled {label}, found {len(matches)}")
    return matches[0]


@dataclass(frozen=True)
class IntegralSignature:
    difference_exponent: str
    plus_exponent: str
    measure: str


def parse_integral(block: str, name: str) -> IntegralSignature:
    matches = list(INTEGRAL_PATTERN.finditer(compact(block)))
    require(len(matches) == 1, f"expected one canonical adjoint integral in {name}, found {len(matches)}")
    match = matches[0]
    difference = match.group("difference_braced") or match.group("difference_plain")
    plus = match.group("plus_braced") or match.group("plus_plain")
    measure_g = match.group("measure_g")
    measure_h = match.group("measure_h")
    require(measure_g == measure_h, f"different measures in the two variables of {name}")
    return IntegralSignature(difference, plus, measure_g)


def canonical_signature(signature: IntegralSignature, a_value: int | None) -> tuple[int, str]:
    if a_value is None:
        require(signature.difference_exponent == "2", "companion minus exponent is not 2")
        difference = 2
    else:
        require(signature.difference_exponent == "2a", "Part III minus exponent is not 2a")
        difference = 2 * a_value
    require(signature.plus_exponent == "n", "plus exponent is not n")
    return difference, signature.plus_exponent


def require_fragment(text: str, fragment: str, message: str) -> None:
    require(compact(fragment) in compact(text), message)


def theorem_number(companion: str) -> tuple[int, int]:
    main = result_with_label(companion, "thm:main")
    main_start = companion.index(main)
    sections = list(re.finditer(r"\\section\{([^}]+)\}", companion[:main_start]))
    require(len(sections) == 2, f"main theorem section number changed: {len(sections)}")
    require(
        sections[-1].group(1) == "Statement and General Reductions",
        "main theorem is not in Statement and General Reductions",
    )
    section_body = companion[sections[-1].end() : main_start]
    preceding_results = list(
        re.finditer(
            r"\\begin\{(?:definition|theorem|proposition|lemma|corollary|remark)\}",
            section_body,
        )
    )
    require(len(preceding_results) == 1, f"objects before main theorem={len(preceding_results)}")
    require(
        section_body[preceding_results[0].start() :].startswith(r"\begin{definition}"),
        "the object before the main theorem is not its functional definition",
    )
    return len(sections), len(preceding_results) + 1


def validate_semantics(companion: str, part_three: str) -> tuple[IntegralSignature, IntegralSignature]:
    companion_definition = block_with_title(
        companion,
        DEFINITION_PATTERN,
        r"[Adjoint $Q_3$ functional]",
    )
    part_three_definition = block_with_title(
        part_three,
        DEFINITION_PATTERN,
        "[Adjoint-generated cone]",
    )
    companion_signature = parse_integral(companion_definition, "Parts I--II definition")
    part_three_signature = parse_integral(part_three_definition, "Part III definition")
    require(
        canonical_signature(companion_signature, None)
        == canonical_signature(part_three_signature, 1),
        "a=1 does not identify the two integral signatures",
    )
    require(companion_signature.measure == r"\Haar", "companion does not use its Haar measure")
    require(part_three_signature.measure == r"\mu_G", "Part III does not use mu_G")

    require_fragment(
        companion_definition,
        r"For a compact connected Lie group $G$ with simple Lie algebra and normalized Haar probability measure",
        "companion functional scope or Haar normalization changed",
    )
    target_section = part_three[: part_three.index(part_three_definition)]
    require_fragment(
        target_section,
        r"let $\mu_G$ be normalized Haar measure",
        "Part III Haar normalization changed",
    )

    companion_main = result_with_label(companion, "thm:main")
    part_three_import = result_with_label(part_three, "thm:two-minus-import")
    common_scope = r"For every compact connected Lie group $G$ with simple Lie algebra and every integer"
    require_fragment(companion_main, common_scope, "companion group scope changed")
    require_fragment(part_three_import, common_scope, "Part III import group scope changed")
    require_fragment(companion_main, r"$n\ge0$", "companion exponent range changed")
    require_fragment(part_three_import, r"$n\geq0$", "Part III exponent range changed")
    require_fragment(companion_main, r"\Qthree^G(n)\ge0", "companion conclusion changed")
    require_fragment(part_three_import, r"I^G_{1,n}\geq0", "Part III import is not the a=1 row")

    require(theorem_number(companion) == (2, 2), "companion main result is not Theorem 2.2")
    require_fragment(
        part_three,
        r"\cite[Theorem~2.2]{AdjointTwoMinus}",
        "Part III cites the wrong companion theorem number",
    )
    require_fragment(
        part_three,
        r"its functional $Q_3^G(n)$ is exactly $I^G_{1,n}$ above",
        "Part III exact-interface assertion is absent",
    )
    require_fragment(
        part_three,
        r"Theorem~\ref{thm:two-minus-import} below is exactly the row $a=1$",
        "Part III hierarchy does not identify the imported row",
    )

    require_fragment(
        companion,
        r"it is not a proof of the full cone-valued Q3 condition",
        "companion scope disclaimer is absent",
    )
    require_fragment(
        companion,
        r"We do not assert Q3 for a cone of functions, for unequal characters, or for four or more minus factors",
        "companion introduction overstates its scope",
    )
    require_fragment(
        part_three,
        r"Parts I--II: classification and exact certificates, journal manuscript with formal detailed supplement, 2026",
        "Part III companion bibliography identity changed",
    )

    companion_cover = result_with_label(companion, "lem:central-quotient")
    part_three_cover = result_with_label(part_three, "lem:full-central-cover")
    require_fragment(companion_cover, r"finite central cover", "companion central-cover interface absent")
    require_fragment(part_three_cover, r"finite central covering", "Part III central-cover interface absent")
    require_fragment(
        companion_cover,
        r"\Qthree^{\widetilde G}(n)=\Qthree^G(n)",
        "companion central-cover equality changed",
    )
    require_fragment(
        part_three_cover,
        r"I^{\widetilde G}_{a,n}=I^G_{a,n}",
        "Part III central-cover equality changed",
    )
    return companion_signature, part_three_signature


def parse_manifest(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"malformed manifest row {path}:{number}")
        value, recorded = match.groups()
        require(recorded not in records, f"duplicate manifest path {path}:{recorded}")
        records[recorded] = value
    return records


def validate_manifest_bindings(root: Path) -> tuple[str, str, str]:
    companion_path = root / "paper.tex"
    part_three_path = root / "full_q3_extension.tex"
    replay_manifest = root / "replay_sources.sha256"
    part_three_manifest = root / "certificates/full_q3/full_q3_source_manifest.sha256"
    companion_hash = digest(companion_path)
    part_three_hash = digest(part_three_path)
    replay_hash = digest(replay_manifest)
    require(companion_hash in ALLOWED_COMPANION_HASHES, "companion source is not an authorized baseline/final version")
    require(part_three_hash == EXPECTED_PART_THREE_HASH, "Part III source identity changed")
    replay = parse_manifest(replay_manifest)
    full = parse_manifest(part_three_manifest)
    require(replay.get("paper.tex") == companion_hash, "replay manifest does not bind the companion")
    require(full.get("../../paper.tex") == companion_hash, "Part III manifest does not bind the companion")
    require(
        full.get("../../full_q3_extension.tex") == part_three_hash,
        "Part III manifest does not bind the importing manuscript",
    )
    require(
        full.get("../../replay_sources.sha256") == replay_hash,
        "Part III manifest does not bind the live replay manifest",
    )
    return companion_hash, part_three_hash, replay_hash


def replace_in_block(text: str, block: str, old: str, new: str) -> str:
    require(text.count(block) == 1, "mutation target block is not unique")
    require(block.count(old) == 1, f"mutation source is not unique in block: {old}")
    return text.replace(block, block.replace(old, new, 1), 1)


def expect_rejected(name: str, action: Callable[[], object]) -> None:
    try:
        action()
    except InterfaceFailure:
        return
    raise InterfaceFailure(f"adversarial interface mutation was accepted: {name}")


def mutation_self_tests(companion: str, part_three: str) -> int:
    companion_definition = block_with_title(
        companion,
        DEFINITION_PATTERN,
        r"[Adjoint $Q_3$ functional]",
    )
    mutated = replace_in_block(companion, companion_definition, r")^2", r")^4")
    expect_rejected("companion minus exponent", lambda: validate_semantics(mutated, part_three))

    part_three_definition = block_with_title(
        part_three,
        DEFINITION_PATTERN,
        "[Adjoint-generated cone]",
    )
    mutated = replace_in_block(part_three, part_three_definition, r"d\mu_G(h)", r"d\mu_G(g)")
    expect_rejected("Part III second Haar variable", lambda: validate_semantics(companion, mutated))

    companion_main = result_with_label(companion, "thm:main")
    mutated = replace_in_block(companion, companion_main, "simple Lie algebra", "semisimple Lie algebra")
    expect_rejected("companion group scope", lambda: validate_semantics(mutated, part_three))

    part_three_import = result_with_label(part_three, "thm:two-minus-import")
    mutated = replace_in_block(part_three, part_three_import, r"I^G_{1,n}", r"I^G_{2,n}")
    expect_rejected("imported hierarchy row", lambda: validate_semantics(companion, mutated))

    require(part_three.count(r"\cite[Theorem~2.2]{AdjointTwoMinus}") == 1, "citation mutation target changed")
    mutated = part_three.replace(
        r"\cite[Theorem~2.2]{AdjointTwoMinus}",
        r"\cite[Theorem~2.3]{AdjointTwoMinus}",
        1,
    )
    expect_rejected("companion theorem number", lambda: validate_semantics(companion, mutated))
    return 5


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="ginibre_q3 source root (default: root containing this audit)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    companion = (root / "paper.tex").read_text(encoding="utf-8")
    part_three = (root / "full_q3_extension.tex").read_text(encoding="utf-8")
    companion_signature, part_three_signature = validate_semantics(companion, part_three)
    companion_hash, part_three_hash, replay_hash = validate_manifest_bindings(root)
    mutations = mutation_self_tests(companion, part_three)
    print(
        "FULL_Q3_TWO_MINUS_INTERFACE "
        f"companion_signature=minus^{companion_signature.difference_exponent},"
        f"plus^{companion_signature.plus_exponent} "
        f"part3_signature=minus^{part_three_signature.difference_exponent},plus^{part_three_signature.plus_exponent} "
        "substitution=a1 theorem=2.2 group_scope=compact_connected_simple "
        "range=n_ge_0 Haar_normalizations=2 central_cover_interfaces=2 "
        f"companion_sha256={companion_hash} part3_sha256={part_three_hash} "
        f"replay_manifest_sha256={replay_hash} manifest_bindings=4 "
        f"mutations_rejected={mutations}"
    )
    print("FULL_Q3_TWO_MINUS_INTERFACE VERIFICATION: ALL PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (InterfaceFailure, OSError, UnicodeError) as error:
        raise SystemExit(f"FULL_Q3_TWO_MINUS_INTERFACE FAILURE: {error}") from error
