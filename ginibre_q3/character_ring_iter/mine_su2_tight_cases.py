#!/usr/bin/env python3
"""Mine the engine's boundary regimes for the factor-axis induction.

The regimes that need the heavy certificate classes -- Farey merges, Abel
chains, exact leaves -- are the tight and zero cases of a level: exactly the
boundary structure any uniform-in-k induction must reproduce.  The engine
emits them with `--dump-certified`; this script aggregates that output into
the shape the factor-axis program needs.

For each boundary regime it decodes the chamber: which labels are active,
with which signs and minimum-exponent parities.  The output groups boundary
regimes by their signed support pattern, so recurring patterns across levels
become visible -- a pattern that needs a Farey merge at level 6 and again at
level 8 is a family, and families are what an induction quantifies over.

Usage:
  mine_su2_tight_cases.py LEVEL < dump.txt
  ./verify_su2_orbit_amgm_stack --level 6 --decompose --dump-certified \
      | mine_su2_tight_cases.py 6
"""

from __future__ import annotations

import collections
import re
import sys

CERT = re.compile(
    r"CERT support=(\d+) parity=(\d+) res=(\d+) k=(\d+) class=(\w+) nodes=(\d+)"
)


def decode(support: int, parity: int, labels: int):
    """Signed word skeleton from the chamber encoding."""
    signs, z = [], support
    for _ in range(labels):
        d = z % 3
        z //= 3
        signs.append(0 if d == 0 else (1 if d == 1 else -1))
    active = [l for l in range(labels) if signs[l]]
    word = []
    for i, l in enumerate(active):
        p0 = 1 if (parity >> i) & 1 else 2
        word.append((l + 1, signs[l], p0))       # (label value, sign, min exponent)
    return tuple(word)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    labels = int(sys.argv[1])
    by_class = collections.Counter()
    families = collections.defaultdict(lambda: collections.Counter())
    for line in sys.stdin:
        m = CERT.search(line)
        if not m:
            continue
        support, parity, _res, _k = (int(m[i]) for i in range(1, 5))
        cls = m[5]
        by_class[cls] += 1
        if cls in ("farey", "abel", "leaves"):
            families[decode(support, parity, labels)][cls] += 1

    print("certified regimes by closing class:")
    for cls, n in by_class.most_common():
        print(f"  {cls:<8} {n}")

    print(f"\nboundary families (signed word skeletons needing farey/abel/leaves): "
          f"{len(families)}")
    ranked = sorted(families.items(), key=lambda kv: -sum(kv[1].values()))
    for word, cnt in ranked[:30]:
        pretty = " ".join(
            f"V{a}{'+' if s > 0 else '-'}^{p}{'+2r' if True else ''}" for a, s, p in word
        )
        print(f"  {sum(cnt.values()):>4}  {dict(cnt)}  {pretty}")
    if len(ranked) > 30:
        print(f"  ... and {len(ranked) - 30} more families")
    return 0


if __name__ == "__main__":
    sys.exit(main())
