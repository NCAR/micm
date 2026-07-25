# Performance checks

MICM ships a small Chapman micro-benchmark used to catch hot-path regressions
in the Rosenbrock `Solve()` function. This document explains how to run the
same checks the CI uses.

## Contents

- `test/benchmark/chapman_bench.cpp` — the benchmark itself. Runs the same
  7-reaction Chapman mechanism as `test_chapman_integration` but scaled up so
  per-`Solve()` work dominates over per-call overhead.
- `scripts/bench_chapman.sh` — wall-clock timing driver (noisy; use for quick
  local development iteration).
- `scripts/profile_chapman.sh` — callgrind driver that produces deterministic
  instruction counts (used by CI regression gates).
- `scripts/compare_chapman.py` — diffs two profile outputs and fails when the
  PR adds instructions for any matrix ordering (standard, vector1,2,4,8).
- `.github/workflows/perf-regression.yml` — CI workflow that runs the
  profile against `main` on every PR.

## Building the benchmark

The benchmark is part of the normal CMake build and doesn't add ctest
targets, so nothing runs by default.

```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release
cmake --build build --target chapman_bench --parallel
```

The binary lands at `build/chapman_bench`.

## Wall-clock timing

Good for quick local iteration; **not** for accept/reject decisions because
wall-clock is affected by CPU load, thermals, and turbo.

```bash
scripts/bench_chapman.sh                         # defaults: 10000 cells, 30 steps
scripts/bench_chapman.sh build 20000 50          # 20k cells, 50 steps
scripts/bench_chapman.sh build 10000 30 vector4  # just vector4
```

Sample output:

```
kind         best_ms
standard      635.24
vector1       664.13
vector2       578.55
vector4       411.94
vector8       368.51
```

## Instruction counts (deterministic)

`profile_chapman.sh` runs the benchmark under `valgrind --tool=callgrind`.
The benchmark itself flips callgrind instrumentation on only around the
timed `Solve()` loop (via `CALLGRIND_START_INSTRUMENTATION` /
`CALLGRIND_TOGGLE_COLLECT`), so the reported instruction count covers
**exactly Rosenbrock `Solve()` and everything it calls** — not solver
construction, not state initialization, not the rate-constant update, and
not the warm-up call. Instruction counts are deterministic across runs,
machines, and CPU generations, which is what makes them suitable for
regression gates.

```bash
sudo apt-get install valgrind                    # ships the client-request header too
scripts/profile_chapman.sh                       # defaults: 2000 cells, 5 steps
scripts/profile_chapman.sh build 2000 5 vector4  # just vector4
```

Sample output (TSV — copy into a spreadsheet if you like):

```
kind              instructions
standard             634826135
vector1              622410872
vector2              550138401
vector4              418247669
vector8              384791122
```

Because callgrind has ~50× overhead, the defaults use smaller `CELLS` and
`STEPS`. Any real hot-path change is still visible at that size; cost is
linear in cells x steps.

## Comparing two runs

```bash
scripts/profile_chapman.sh build 2000 5 > base.txt
# ... make changes and rebuild ...
scripts/profile_chapman.sh build 2000 5 > pr.txt
scripts/compare_chapman.py base.txt pr.txt
```

Output:

```
kind                 base              pr           delta   delta%
standard      634,826,135     634,900,000         +73,865   +0.01%  <-- regression
vector1       622,410,872     620,000,000      -2,410,872   -0.39%
...

FAIL: PR adds instructions to the hot path for: standard
```

Exit code is 0 when nothing regressed and 1 otherwise. Use `--tolerance N`
to allow a small absolute instruction-count increase per ordering; **the
default is 0**, matching what CI enforces.

## Typical workflow when refactoring a hot-path function

1. Before you start, capture a baseline against `main`:
   ```bash
   git checkout main && cmake --build build --target chapman_bench --parallel
   scripts/profile_chapman.sh build > /tmp/base.txt
   ```
2. Refactor.
3. Rebuild and compare:
   ```bash
   cmake --build build --target chapman_bench --parallel
   scripts/profile_chapman.sh build > /tmp/pr.txt
   scripts/compare_chapman.py /tmp/base.txt /tmp/pr.txt
   ```
4. If any ordering regressed, `callgrind_annotate /tmp/cg_<kind>.out` on the
   file the profile script left behind gives you a per-source-line
   breakdown. Diff the "hot lines" between base and PR — the delta is
   always attributable to specific lines because instruction counts are
   deterministic.

## CI: `.github/workflows/perf-regression.yml`

The workflow runs on every PR to `main`:

1. Checks out the PR head and the base branch side-by-side.
2. Overlays the PR's benchmark harness onto the base checkout (so both
   revisions run the exact same benchmark code — otherwise, changing the
   benchmark itself could mask regressions).
3. Builds `chapman_bench` from both.
4. Runs `profile_chapman.sh` against both.
5. Fails if any matrix ordering has strictly more instructions on the PR.

The gate is intentionally strict (tolerance = 0). If a change adds
instructions on purpose (e.g. new required work for a feature), bump the
`--tolerance` value in the workflow, or add a `[skip perf]` label check as
appropriate.
