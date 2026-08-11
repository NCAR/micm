#!/usr/bin/env bash
#
# micm solver benchmark — wall-clock timing across matrix orderings.
#
# Runs the standalone micm_bench binary once per matrix ordering (best-of-3)
# and prints a compact table. Wall-clock
# time is noisy; if you need a deterministic measurement (e.g. for regression
# gates), use profile_micm.sh instead.
#
# Usage:
#   scripts/bench_micm.sh [BUILD_DIR] [CELLS] [STEPS] [BACKEND] [LU_TYPE] \
#                         [LU_ALGORITHM] [MECHANISM] [MATRIX...]
#
# Defaults:
#   BUILD_DIR    = build
#   CELLS        = 10000
#   STEPS        = 30
#   BACKEND      = cpu (other option: gpu)
#   LU_TYPE      = in-place (other option: separate)
#   LU_ALGORITHM = mozart (other option: doolittle)
#   MECHANISM    = chapman (other option: ts1)
#   MATRIX       = standard vector1 vector2 vector4 vector8 vector128
#
# ts1 has 547 reactions against Chapman's 7, so start with far fewer cells.
# Keep CELLS a multiple of 128, or the vector128 ordering pads its last group.
#
# The gpu backend supports only the vector orderings, so name them explicitly:
#   scripts/bench_micm.sh build 256 30 gpu in-place mozart ts1 vector1 vector4
#
# Requires: micm_bench built inside BUILD_DIR, i.e. configured with
# -D MICM_ENABLE_BENCHMARK=ON.

set -euo pipefail

build="${1:-build}"
cells="${2:-10000}"
steps="${3:-30}"
backend="${4:-cpu}"
lu_type="${5:-in-place}"
lu_algorithm="${6:-mozart}"
mechanism="${7:-chapman}"
shift $(( $# > 7 ? 7 : $# )) || true
kinds=("$@")
if [[ ${#kinds[@]} -eq 0 ]]; then
  kinds=(standard vector1 vector2 vector4 vector8 vector128)
fi

bin="${build}/micm_bench"
if [[ ! -x "$bin" ]]; then
  echo "micm_bench not found at $bin — build it with:" >&2
  echo "  cmake -S . -B $build -D CMAKE_BUILD_TYPE=Release -D MICM_ENABLE_BENCHMARK=ON" >&2
  echo "  cmake --build $build --target micm_bench --parallel" >&2
  exit 1
fi

echo "mechanism = $mechanism; backend = $backend; LU = $lu_algorithm / $lu_type"
# A machine-readable copy of the mechanism, to match profile_micm.sh.
echo "# mechanism=$mechanism"
printf '%-9s %10s\n' "kind" "best_ms"
for kind in "${kinds[@]}"; do
  best=""
  for _ in 1 2 3; do
    line=$("$bin" "$cells" "$steps" 30.0 "$backend" "$kind" "$lu_type" "$lu_algorithm" "$mechanism")
    ms=$(printf '%s' "$line" | sed -E 's/.*elapsed_ms=([0-9.]+).*/\1/')
    if [[ -z "$best" ]] || awk "BEGIN{ exit !($ms < $best) }"; then
      best="$ms"
    fi
  done
  printf '%-9s %10.2f\n' "$kind" "$best"
done
