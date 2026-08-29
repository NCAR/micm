# Baseline — DAE constraint-init line-search work

- **Date:** 2026-08-29
- **Branch:** `dae-init-damping`, cut from `fix/dae-constraint-tolerance-measure` (NCAR/micm PR #1083)
- **SHA:** `f72c25df4551754f3eeacefa8883c9639b494f06`
- **Plan:** `docs/superpowers/plans/2026-08-29-dae-init-cold-start.md`
- **Spec:** `docs/superpowers/specs/2026-08-28-dae-init-cold-start-spec.md`

Run before any source change, in three fresh build directories. The pre-existing
`build/`, `build-float/` and `build-kokkos/` directories were **not** reused: they
were configured on the reference branch `dae-rosenbrock-benchmark`, hold stale
executables, and `build/Testing/Temporary/LastTestsFailed.log` recorded
`constraint_initialization` as its last failure set. Nothing in this tree evidenced
a green baseline until this run.

## Configure

```bash
cmake -S . -B bld-double -DCMAKE_BUILD_TYPE=Release -DMICM_ENABLE_TESTS=ON \
  -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake -S . -B bld-float  -DCMAKE_BUILD_TYPE=Release -DMICM_ENABLE_TESTS=ON \
  -DMICM_USE_DOUBLE=OFF -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake -S . -B bld-kokkos -DCMAKE_BUILD_TYPE=Release -DMICM_ENABLE_TESTS=ON \
  -DMICM_ENABLE_KOKKOS=ON -DKokkos_ENABLE_SERIAL=ON \
  -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
```

`FETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER` is required: Homebrew GoogleTest 1.17.0
is installed and would otherwise be found, and tests built against it abort at
startup in this environment.

## Results — all three green

| Configuration | ctest result | Test time |
|---|---|---|
| double (`bld-double`) | `100% tests passed, 0 tests failed out of 65` | 16.51 s |
| Kokkos serial (`bld-kokkos`) | `100% tests passed, 0 tests failed out of 90` | 44.71 s |
| float (`bld-float`) | `100% tests passed, 0 tests failed out of 6` | 1.05 s |

The float configuration builds and runs only the six constraint/DAE targets that
spec R4 refers to:

```bash
cmake --build bld-float -j10 --target test_constraint_initialization \
  test_equilibrium_constraint test_linear_constraint test_external_model_constraints \
  test_dae_constraint_overshoot test_dae_algebraic_error_insensitivity
(cd bld-float && ctest -R 'constraint_initialization|equilibrium_constraint|linear_constraint|external_model_constraints|dae_constraint_overshoot|dae_algebraic_error_insensitivity')
```

Note `-R 'constraint|dae'` would select a seventh test (`constraint_set`) and is
**not** the R4 set.

This baseline confirms the spec's claim that PR #1083 is green in all three
configurations, and is the reference every later "still green" gate compares against.
