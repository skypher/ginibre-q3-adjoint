#!/usr/bin/env python3
"""Independently compare every consumed exceptional moment with BPV.

The standard replay regenerates the complete theorem-consumed ranges from
the Cartan data by exact Racah--Speiser iteration and character pairing.  This
separate audit compares those local ledgers with the archived BPV table; the
table is a control, not a supplier for a theorem-consumed value.
"""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANCILLARY = ROOT / "references/bourjaily_plesser_vergu_2024_adjoint_tensor_ranks.tex"
MOMENT = re.compile(r"\bm_(\d+)\s*=\s*(-?\d+)")
GROUPS = {
    "G2": ("g[2]", ROOT / "character_ring_iter/logs/g2_200.log", 38),
    "F4": ("f[4]", ROOT / "character_ring_iter/logs/f4_220c.log", 65),
    "E6": ("e[6]", ROOT / "character_ring_iter/logs/e6_80.log", 80),
    "E7": ("e[7]", ROOT / "character_ring_iter/logs/e7_70.log", 70),
    "E8": ("e[8]", ROOT / "character_ring_iter/logs/e8_70.log", 70),
}


def moment_rows(path: Path) -> dict[int, int]:
    rows: dict[int, int] = {}
    for degree_text, value_text in MOMENT.findall(path.read_text(encoding="utf-8")):
        degree, value = int(degree_text), int(value_text)
        if degree in rows and rows[degree] != value:
            raise AssertionError(f"conflicting m_{degree} in {path}")
        rows[degree] = value
    return rows


def ancillary_rows(text: str, key: str) -> dict[int, int]:
    match = re.search(
        r"\{" + re.escape(key) + r",\s*\{(.*?)\}\}", text, re.DOTALL
    )
    if match is None:
        raise AssertionError(f"missing ancillary block {key}")
    # The ancillary lists begin with tensor power one: 0,1,1,5,... .
    # TeX source line-wraps very long integers as ``123\`` newline ``456``.
    # Remove only that explicit continuation before tokenizing values.
    payload = re.sub(r"\\\s*", "", match.group(1))
    values = [int(value) for value in re.findall(r"\d+", payload)]
    return {degree: value for degree, value in enumerate(values, 1)}


def main() -> None:
    text = ANCILLARY.read_text(encoding="utf-8")
    checked = 0
    for group, (key, local_path, maximum) in GROUPS.items():
        source = ancillary_rows(text, key)
        local = moment_rows(local_path)
        expected_degrees = set(range(1, maximum + 1))
        missing_source = sorted(expected_degrees - source.keys())
        missing_local = sorted(set(range(maximum + 1)) - local.keys())
        if missing_source or missing_local:
            raise AssertionError(
                f"{group} incomplete scope: source={missing_source} local={missing_local}"
            )
        mismatches = [
            degree for degree in expected_degrees
            if local[degree] != source[degree]
        ]
        if mismatches:
            raise AssertionError(f"{group} ancillary mismatches: {mismatches}")
        if local[0] != 1:
            raise AssertionError(f"{group} has m_0={local[0]}, expected 1")
        checked += maximum + 1
        print(
            f"EXCEPTIONAL_ANCILLARY group={group} moments={maximum + 1} "
            "mismatches=0"
        )

    # Part III consumes E8 through degree 100 from the same primary ancillary
    # artifact.  Check that the maintained transcription is byte-for-value
    # exact in every one of its 30 rows.
    e8_source = ancillary_rows(text, "e[8]")
    e8_transcription = moment_rows(
        ROOT / "references/arxiv_2412_21189_e8_m71_m100.txt"
    )
    if set(e8_transcription) != set(range(71, 101)):
        raise AssertionError("E8 m_71..m_100 transcription scope mismatch")
    mismatches = [
        degree for degree in range(71, 101)
        if e8_transcription[degree] != e8_source[degree]
    ]
    if mismatches:
        raise AssertionError(f"E8 transcription mismatches: {mismatches}")
    checked += 30
    print("EXCEPTIONAL_ANCILLARY e8_transcription=71..100 mismatches=0")
    print(f"EXCEPTIONAL_ANCILLARY checked_values={checked} failures=0")
    print("EXCEPTIONAL_ANCILLARY VERIFICATION: ALL PASS")


if __name__ == "__main__":
    main()
