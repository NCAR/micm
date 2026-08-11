#!/usr/bin/env python3
"""Compare solver-benchmark instruction counts from two profile runs.

Reads two TSV files produced by ``scripts/profile_micm.sh`` (columns:
kind, instructions) and prints a side-by-side comparison. Exits non-zero if
any matrix ordering has a higher instruction count in the PR than in the
base, subject to an absolute-count tolerance (default 0). This is intended
as the check step in the perf-regression CI workflow.

Both files must come from the same mechanism. Each file carries a
"# mechanism=<name>" marker line, and this script exits non-zero when the two
markers disagree. A file written before the marker existed carries none, and
the check passes.

Usage:
    scripts/compare_micm.py base.txt pr.txt [--tolerance N]

Instruction counts are deterministic; any real hot-path regression will
show up as a strictly positive delta. A small tolerance is available for
noise-free but harmless attribution shifts (e.g. inlining boundary changes
between compiler versions).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


MECHANISM_MARKER = "# mechanism="


def parse(path: Path) -> tuple[str | None, dict[str, int]]:
    """Parse a profile_micm.sh TSV file into (mechanism, {kind: instructions}).

    The mechanism is None when the file carries no marker line.
    """
    mechanism: str | None = None
    result: dict[str, int] = {}
    for line in path.read_text().splitlines():
        if line.startswith(MECHANISM_MARKER):
            mechanism = line[len(MECHANISM_MARKER):].strip()
            continue
        parts = line.split()
        if len(parts) != 2 or parts[0] == "kind":
            continue
        try:
            result[parts[0]] = int(parts[1])
        except ValueError:
            continue
    return mechanism, result


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

    base_mechanism, base = parse(args.base)
    pr_mechanism, pr = parse(args.pr)

    if not base or not pr:
        print(f"error: could not parse counts from {args.base} or {args.pr}", file=sys.stderr)
        return 2

    # A count from one mechanism against a count from another is meaningless,
    # and the delta looks like an enormous regression.
    if base_mechanism is not None and pr_mechanism is not None and base_mechanism != pr_mechanism:
        print(
            f"error: {args.base} measures '{base_mechanism}' and {args.pr} measures "
            f"'{pr_mechanism}'; profile both with the same mechanism",
            file=sys.stderr,
        )
        return 2

    mechanism = base_mechanism or pr_mechanism
    if mechanism is not None:
        print(f"mechanism: {mechanism}\n")

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
    print("\nOK: PR does not increase the hot-path instruction count.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
