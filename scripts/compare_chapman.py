#!/usr/bin/env python3
"""Compare Chapman-benchmark instruction counts from two profile runs.

Reads two TSV files produced by ``scripts/profile_chapman.sh`` (columns:
kind, instructions) and prints a side-by-side comparison. Exits non-zero if
any matrix ordering has a higher instruction count in the PR than in the
base, subject to an absolute-count tolerance (default 0). This is intended
as the check step in the perf-regression CI workflow.

Usage:
    scripts/compare_chapman.py base.txt pr.txt [--tolerance N]

Instruction counts are deterministic; any real hot-path regression will
show up as a strictly positive delta. A small tolerance is available for
noise-free but harmless attribution shifts (e.g. inlining boundary changes
between compiler versions).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse(path: Path) -> dict[str, int]:
    """Parse a profile_chapman.sh TSV file into {kind: instructions}."""
    result: dict[str, int] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2 or parts[0] == "kind":
            continue
        try:
            result[parts[0]] = int(parts[1])
        except ValueError:
            continue
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, help="baseline profile TSV")
    parser.add_argument("pr", type=Path, help="candidate profile TSV")
    parser.add_argument(
        "--tolerance",
        type=int,
        default=0,
        help="allowed absolute instruction-count increase per ordering (default: 0)",
    )
    args = parser.parse_args()

    base = parse(args.base)
    pr = parse(args.pr)

    if not base or not pr:
        print(f"error: could not parse counts from {args.base} or {args.pr}", file=sys.stderr)
        return 2

    kinds = sorted(set(base) & set(pr))
    if not kinds:
        print("error: no matrix orderings in common between base and pr", file=sys.stderr)
        return 2

    print(f"{'kind':<9} {'base':>15} {'pr':>15} {'delta':>15} {'delta%':>8}")
    regressed: list[str] = []
    for kind in kinds:
        b = base[kind]
        p = pr[kind]
        d = p - b
        pct = (d / b * 100.0) if b > 0 else 0.0
        marker = ""
        if d > args.tolerance:
            regressed.append(kind)
            marker = "  <-- regression"
        print(f"{kind:<9} {b:>15,} {p:>15,} {d:>+15,} {pct:>+7.2f}%{marker}")

    if regressed:
        print(
            f"\nFAIL: PR adds instructions to the hot path for: {', '.join(regressed)}",
            file=sys.stderr,
        )
        return 1
    print("\nOK: PR does not increase Chapman hot-path instruction count.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
