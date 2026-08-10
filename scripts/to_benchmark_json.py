#!/usr/bin/env python3
"""Convert bench_micm.sh or profile_micm.sh TSV output to github-action-benchmark JSON."""
import sys, json, argparse

parser = argparse.ArgumentParser()
parser.add_argument("--unit", default="instructions")
parser.add_argument("--value-type", choices=["int", "float"], default="int")
args = parser.parse_args()

rows = []
for line in sys.stdin:
    parts = line.split()
    if len(parts) != 2 or parts[0] in ("kind", "backend", "best_ms"):
        continue
    try:
        val = int(parts[1]) if args.value_type == "int" else float(parts[1])
        rows.append({"name": parts[0], "unit": args.unit, "value": val})
    except ValueError:
        continue

json.dump(rows, sys.stdout, indent=2)
print()
