#!/usr/bin/env bash
#
# Chapman benchmark: deterministic instruction counts via callgrind.
#
# Runs the standalone chapman_bench binary once per matrix ordering under
# valgrind/callgrind. Instrumentation and collection are both disabled at
# startup; the benchmark itself flips them on (via CALLGRIND_START_INSTRUMENTATION
# / CALLGRIND_TOGGLE_COLLECT) only around the timed Solve() loop. That means
# the reported instruction count covers exactly Rosenbrock Solve() and
# everything it calls: not solver construction, not state init, not rate
# constants, and not the warm-up call.
#
# Instruction counts are deterministic; they don't vary with CPU load,
# thermals, or turbo, which makes them suitable for CI regression gates on
# the hot path.
#
# Because callgrind has ~50x overhead, this uses smaller CELLS/STEPS defaults
# than bench_chapman.sh. That's fine for regression checks: any real hot-path
# change is visible even at small counts.
#
# Usage:
#   scripts/profile_chapman.sh [BUILD_DIR] [CELLS] [STEPS] [KIND...]
#
# Defaults:
#   BUILD_DIR = build
#   CELLS     = 2000
#   STEPS     = 5
#   KINDS     = standard vector1 vector2 vector4 vector8
#
# Requires: chapman_bench built inside BUILD_DIR, and valgrind on PATH.
# For the tightest scoping the benchmark should also see <valgrind/callgrind.h>
# at compile time (valgrind-devel / valgrind package on most distros).
#
# Output: a TSV table on stdout with columns "kind" and "instructions"
# (space-padded). Raw callgrind files are written to $OUT (default /tmp).

set -euo pipefail

build="${1:-build}"
cells="${2:-2000}"
steps="${3:-5}"
shift $(( $# > 3 ? 3 : $# )) || true
kinds=("$@")
if [[ ${#kinds[@]} -eq 0 ]]; then
  kinds=(standard vector1 vector2 vector4 vector8)
fi

bin="${build}/chapman_bench"
if [[ ! -x "$bin" ]]; then
  echo "chapman_bench not found at $bin — build it with:" >&2
  echo "  cmake -S . -B $build -D CMAKE_BUILD_TYPE=Release" >&2
  echo "  cmake --build $build --target chapman_bench --parallel" >&2
  exit 1
fi
if ! command -v valgrind >/dev/null; then
  echo "valgrind not found — required for callgrind profiling" >&2
  exit 1
fi

out="${OUT:-/tmp}"
mkdir -p "$out"

printf '%-9s %20s\n' "kind" "instructions"
for kind in "${kinds[@]}"; do
  cg="${out}/cg_${kind}.out"
  # --instr-atstart=no: skip JIT-instrumentation of setup entirely (fast).
  # --collect-atstart=no: even after the benchmark turns instrumentation on,
  #     don't count until CALLGRIND_TOGGLE_COLLECT runs.
  valgrind --tool=callgrind --callgrind-out-file="$cg" \
      --instr-atstart=no --collect-atstart=no \
      "$bin" "$cells" "$steps" 30.0 "$kind" >/dev/null 2>&1
  # PROGRAM TOTALS line from callgrind_annotate: deterministic instruction count.
  # NB: grep -m1 exits early on match, which SIGPIPEs callgrind_annotate; that's
  # fine under pipefail because we redirect stderr and the count is captured
  # before the terminating signal reaches us.
  totals=$(callgrind_annotate --threshold=100 --auto=no "$cg" 2>/dev/null \
           | grep -m1 "PROGRAM TOTALS" || true)
  ir=${totals%% *}
  ir=${ir//,/}
  printf '%-9s %20s\n' "$kind" "$ir"
done
