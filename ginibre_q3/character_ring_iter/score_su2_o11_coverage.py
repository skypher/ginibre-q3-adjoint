#!/usr/bin/env python3
"""Score sector theorems against the reconstructed O_11 residual key index.

Reads the census transcript produced by analyze_su2_o11_frontier.cpp and
reports how many residual keys any given sector theorem can reach.

The point is to answer a question that was unanswerable while the level-eleven
package was lost: given a new sector theorem, what fraction of the hard
regimes does it actually close?

Key encoding, matching the analyzer:
  support   base-3 over the five orbit labels B_1..B_5,
            digit 0 absent, 1 plus, 2 minus (least significant = B_1);
  parity    bitmask over active labels, set = odd minimum exponent;
  residual  bitmask over active labels, set = that coordinate is free.

Usage:
  score_su2_o11_coverage.py [census.log]
"""

from __future__ import annotations

import collections
import re
import sys
from pathlib import Path

LABELS = 5
FAIL_RE = re.compile(
    r"FAIL support=(\d+) parity=(\d+) residual=(\d+) vars=([\d,]*) "
    r"neg=(\d+) pos=(\d+)"
)


def decode_support(support: int) -> tuple[list[int], list[int]]:
    """Return (signs, active) with signs[l] in {0,+1,-1} for label B_(l+1)."""
    signs, z = [], support
    for _ in range(LABELS):
        d = z % 3
        z //= 3
        signs.append(0 if d == 0 else (1 if d == 1 else -1))
    return signs, [l for l in range(LABELS) if signs[l] != 0]


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parents[1] / "certificates" / "su2_o11_frontier_census.log"
    )
    if not src.is_file():
        print(f"census transcript not found: {src}", file=sys.stderr)
        return 1

    keys = []
    header = ""
    for line in src.read_text().splitlines():
        if line.startswith("O11_DIRECT_CENSUS"):
            header = line
        m = FAIL_RE.match(line)
        if m:
            support, parity, residual = int(m[1]), int(m[2]), int(m[3])
            signs, active = decode_support(support)
            varlist = [int(v) - 1 for v in m[4].split(",") if v]
            if not set(varlist).issubset(set(active)):
                print(f"inconsistent key: {line}", file=sys.stderr)
                return 1
            keys.append((support, parity, residual, tuple(active), tuple(signs)))

    if not keys:
        print("no residual keys parsed", file=sys.stderr)
        return 1

    total = len(keys)
    print(header)
    print(f"\nresidual keys parsed: {total}\n")

    by_n = collections.Counter(len(k[3]) for k in keys)
    print("by number of active labels:")
    cum = 0
    for n in sorted(by_n):
        cum += by_n[n]
        print(f"  {n} labels {by_n[n]:>6}   cumulative {cum:>6}  ({100*cum/total:5.1f}%)")

    two = [k for k in keys if len(k[3]) == 2]
    print(f"\ntwo-label keys, by pair ({len(two)} total):")
    for pair, c in sorted(collections.Counter(k[3] for k in two).items()):
        print(f"  B{pair[0]+1},B{pair[1]+1}  {c}")

    # Sector scoring.  A sector is a predicate on the active label set.
    sectors = {
        "all two-label pairs (fan verifier scope)": lambda a: len(a) <= 2,
        "B1,B5 only (top-ray reduction)": lambda a: set(a) <= {0, 4},
        "at most three labels": lambda a: len(a) <= 3,
    }
    print("\nsector coverage:")
    for name, pred in sectors.items():
        hit = sum(1 for k in keys if pred(k[3]))
        print(f"  {name:<42} {hit:>5} / {total}  ({100*hit/total:5.1f}%)")

    high = sum(1 for k in keys if len(k[3]) >= 4)
    print(
        f"\nfour- and five-label keys: {high} / {total} ({100*high/total:.1f}%)"
        "\nThat mass is what the bulk AM-GM machinery covered, and no"
        "\ntwo-label theorem reaches it."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
