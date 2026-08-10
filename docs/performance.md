# Performance checks

MICM ships a small micro-benchmark used to catch hot-path regressions
in the Rosenbrock `Solve()` function. This document explains how to run the
same checks the CI uses.

## Contents

- `test/benchmark/micm_bench.cpp` — the benchmark itself. Selects one
  mechanism and one solver configuration, then scales the mechanism up so
  per-`Solve()` work dominates over per-call overhead. Run
  `build/micm_bench list` to print every valid combination.
- `test/benchmark/bench_harness.hpp` — the shared timing loop, the callgrind
  scoping, and the table that maps a combination to a run. Each solver
  configuration is generated, not written out by hand.
- `test/benchmark/chapman_mechanism.hpp` — the 7-reaction Chapman mechanism,
  the same system as `test_chapman_integration`. This is the default and the
  one the CI charts track.
- `test/benchmark/ts1_mechanism.hpp` — MOZART-TS1, 210 species and 547
  reactions, every one written out as ordinary micm types. Use it to see how
  the solver scales with mechanism size. **This file is generated. Do not edit
  it.**
- `test/benchmark/import_ts1.py` — generates the header above. micm has no
  mechanism-configuration parser, so the script uses the musica Python API to
  read the mechanism and emits the micm calls that build it. The header is
  checked in, so a normal build needs neither musica nor a network. Re-run the
  script only when the mechanism changes:

  ```bash
  pip install musica
  test/benchmark/import_ts1.py          # musica ships the TS1 configuration
  test/benchmark/import_ts1.py --source path/to/musica/configs/v1/ts1
  ```

  `CreateGasPhase` and `CreateProcesses` are plain functions rather than
  templates, so the mechanism compiles once instead of once per solver
  configuration.
- `scripts/bench_micm.sh` — wall-clock timing driver (noisy; use for quick
  local development iteration).
- `scripts/profile_micm.sh` — callgrind driver that produces deterministic
  instruction counts (used by CI regression gates).
- `scripts/compare_micm.py` — diffs two profile outputs and exits non-zero when
  the PR adds instructions for any matrix ordering (standard, vector1,2,4,8,128).
- `.github/workflows/perf-regression.yml` — CI workflow that runs the
  profile against `main` on every PR.

## Building the benchmark

`MICM_ENABLE_BENCHMARK` defaults to **OFF**, so a normal build never compiles
the benchmark. It registers no ctest, and only the jobs that measure
performance turn it on. The option depends on `MICM_ENABLE_TESTS`: with the
tests off, the benchmark stays off.

```bash
cmake -S . -B build -D CMAKE_BUILD_TYPE=Release -D MICM_ENABLE_BENCHMARK=ON
cmake --build build --target micm_bench --parallel
```

The binary lands at `build/micm_bench`. It always holds both mechanisms, and
TS1 dominates the compile time.

Three workflows pass the option: `benchmark-charts.yml`, `perf-regression.yml`,
and `runner.yml`. No other CI job compiles the benchmark.

## Wall-clock timing

Good for quick local iteration; **not** for accept/reject decisions because
wall-clock is affected by CPU load, thermals, and turbo.

Both drivers take the same positional arguments:

    BUILD_DIR  CELLS  STEPS  BACKEND  LU_TYPE  LU_ALGORITHM  MECHANISM  [MATRIX...]

Every argument has a default, but a later argument needs every earlier one.
The trailing `MATRIX` list is variadic, so it comes last.

    scripts/bench_micm.sh                                            # defaults: 10000 cells, 30 steps
    scripts/bench_micm.sh build 20000 50                             # 20k cells, 50 steps
    scripts/bench_micm.sh build 10000 30 cpu in-place mozart chapman vector4   # just vector4
    scripts/bench_micm.sh build 256 30 cpu in-place mozart ts1        # the TS1 mechanism

TS1 has 547 reactions against Chapman's 7, so reduce the cell count by a long
way. Keep `CELLS` a multiple of 128: a `vector128` group holds 128 cells, so
any other count makes that one ordering pad its last group and solve more
cells than the other five.

Sample output:

```
mechanism = chapman; backend = cpu; LU = mozart / in-place
kind         best_ms
standard      676.85
vector1       676.51
vector2       428.13
vector4       354.28
vector8       292.94
vector128     213.87
```

## Instruction counts (deterministic)

`profile_micm.sh` runs the benchmark under `valgrind --tool=callgrind`.
The benchmark itself flips callgrind instrumentation on only around the
timed `Solve()` loop (via `CALLGRIND_START_INSTRUMENTATION` /
`CALLGRIND_TOGGLE_COLLECT`), so the reported instruction count covers
**exactly Rosenbrock `Solve()` and everything it calls**.
Instruction counts are deterministic across runs, machines, and CPU
generations.

```bash
sudo apt-get install valgrind            # ships the client-request header too
scripts/profile_micm.sh                  # defaults: 2000 cells, 5 steps
scripts/profile_micm.sh build 2000 5 cpu in-place mozart chapman vector4   # just vector4
```

Sample output (TSV — copy into a spreadsheet if you like):

```
mechanism = chapman; backend = cpu; LU = mozart / in-place
kind              instructions
standard             611346438
vector1              613011430
vector2              406257254
vector4              314970571
vector8              263484576
vector128            220045464
```

Callgrind is much slower than a native run, so the defaults use smaller
`CELLS` and `STEPS`. Any real hot-path change is still visible at that size.

TS1 works here too, but measure its wall-clock cost first. TS1 costs much more
per cell, and callgrind multiplies that again, so pick `CELLS` and `STEPS` to
match:

```bash
scripts/profile_micm.sh build 256 5 cpu in-place mozart ts1
```

The script leaves the raw callgrind files in `$OUT` (default `/tmp`), named
`cg_<mechanism>_<kind>.out`. Set `OUT` to keep two runs apart.

## Comparing two runs

Both files must come from the same mechanism. The TSV holds no mechanism
column, so `compare_micm.py` cannot detect a mismatch.

```bash
scripts/profile_micm.sh build 2000 5 > base.txt
# ... make changes and rebuild ...
scripts/profile_micm.sh build 2000 5 > pr.txt
scripts/compare_micm.py base.txt pr.txt
```

Output:

```
kind                 base              pr           delta   delta%
standard      510,009,297     611,346,438    +101,337,141  +19.87%  <-- regression
vector1       627,328,869     613,011,430     -14,317,439   -2.28%
vector128     227,174,334     220,045,464      -7,128,870   -3.14%
vector2       457,832,087     406,257,254     -51,574,833  -11.27%
vector4       344,510,309     314,970,571     -29,539,738   -8.57%
vector8       276,683,770     263,484,576     -13,199,194   -4.77%
...

FAIL: PR adds instructions to the hot path for: standard
```

Exit code is 0 when nothing regressed and 1 otherwise. Use `--tolerance N`
to allow a small absolute instruction-count increase per ordering; **the
default is 0**, which is the value CI passes.

## Typical workflow when refactoring a hot-path function

1. Before you start, capture a baseline against `main`:
   ```bash
   git checkout main
   cmake -S . -B build -D CMAKE_BUILD_TYPE=Release -D MICM_ENABLE_BENCHMARK=ON
   cmake --build build --target micm_bench --parallel
   OUT=/tmp/cg-base scripts/profile_micm.sh build > /tmp/base.txt
   ```
2. Refactor.
3. Rebuild and compare:
   ```bash
   cmake --build build --target micm_bench --parallel
   OUT=/tmp/cg-pr scripts/profile_micm.sh build > /tmp/pr.txt
   scripts/compare_micm.py /tmp/base.txt /tmp/pr.txt
   ```
4. If any ordering regressed, run
   `callgrind_annotate /tmp/cg-pr/cg_chapman_<kind>.out` on the file the
   profile script left behind. It gives you a per-source-line breakdown. Diff
   the "hot lines" between base and PR — the delta is always attributable to
   specific lines because instruction counts are deterministic. Separate `OUT`
   directories are what keep both sides' files available.

## CI: `.github/workflows/perf-regression.yml`

The workflow runs on every PR to `main`:

1. Checks out the PR head and the base branch side-by-side.
2. Overlays the PR's benchmark harness onto the base checkout (so both
   revisions run the exact same benchmark code — otherwise, changing the
   benchmark itself could mask regressions).
3. Builds `micm_bench` from both.
4. Runs `profile_micm.sh` against both, for Chapman only.
5. Runs `compare_micm.py`, which marks any matrix ordering that has strictly
   more instructions on the PR.

The comparison reports; it does not gate. The step carries
`continue-on-error: true`, so a marked ordering shows as a failed step inside
a passing workflow, and the PR stays mergeable. Read the step's log before you
merge a change to a hot path.

To make the comparison block a merge instead, remove `continue-on-error` from
the "Compare instruction counts" step. The comparison already runs at
tolerance 0, so expect it to catch deliberate additions too; raise
`--tolerance` in the workflow when a change adds instructions on purpose.
