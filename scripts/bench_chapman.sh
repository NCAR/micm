#!/usr/bin/env bash
#
# Chapman benchmark — wall-clock timing across matrix orderings.
#
# Runs the standalone chapman_bench binary (built by the normal cmake build)
# once per matrix ordering (best-of-3) and prints a compact table. Wall-clock
# time is noisy; if you need a deterministic measurement (e.g. for regression
# gates), use profile_chapman.sh instead.
#
# Usage:
#   scripts/bench_chapman.sh [BUILD_DIR] [CELLS] [STEPS] [BACKEND] [LU_TYPE] [LU_ALGORITHM] [MATRIX]
#
# Defaults:
#   BUILD_DIR    = build
#   CELLS        = 10000
#   STEPS        = 30
#   BACKEND      = cpu (other option: gpu)
#   LU_TYPE      = in-place (other option: separate)
#   LU_ALGORITHM = mozart (other option: doolittle)
#   MATRIX       = standard vector1 vector2 vector4 vector8 vector128
#
# Requires: chapman_bench built inside BUILD_DIR (i.e. cmake --build BUILD_DIR
# --target chapman_bench).

set -euo pipefail

build="${1:-build}"
cells="${2:-10000}"
steps="${3:-30}"
backend="${4:-cpu}"
lu_type="${5:-in-place}"
lu_algorithm="${6:-mozart}"
shift $(( $# > 6 ? 6 : $# )) || true
kinds=("$@")
if [[ ${#kinds[@]} -eq 0 ]]; then
  kinds=(standard vector1 vector2 vector4 vector8 vector128)
fi

bin="${build}/chapman_bench"
if [[ ! -x "$bin" ]]; then
  echo "chapman_bench not found at $bin — build it with:" >&2
  echo "  cmake -S . -B $build -D CMAKE_BUILD_TYPE=Release" >&2
  echo "  cmake --build $build --target chapman_bench --parallel" >&2
  exit 1
fi

echo "backend = $backend; LU = $lu_algorithm / $lu_type"
printf '%-9s %10s\n' "kind" "best_ms"
for kind in "${kinds[@]}"; do
  best=""
  for _ in 1 2 3; do
    line=$("$bin" "$cells" "$steps" 30.0 "$backend" "$kind" "$lu_type" "$lu_algorithm")
    ms=$(printf '%s' "$line" | sed -E 's/.*elapsed_ms=([0-9.]+).*/\1/')
    if [[ -z "$best" ]] || awk "BEGIN{ exit !($ms < $best) }"; then
      best="$ms"
    fi
  done
  printf '%-9s %10.2f\n' "$kind" "$best"
done
