# Active Plan — DAE Rosenbrock Performance and Rigor

The July 15 constraint tolerance and convergence experiment plan is complete:
the weighted-correction initialization shipped in `d7bf2087` and its results
are recorded in
`docs/superpowers/notes/2026-07-15-dae-constraint-convergence-experiments.md`
(the previous version of this file is preserved in git history at that
commit).

The current plan is
`docs/superpowers/plans/2026-07-20-dae-performance-and-rigor.md`: a
performance thread (step-size persistence, work–precision measurement rig,
constraint-evaluation cost, Schur reduction of algebraic rows, norm policy and
diagnostics, RODAS-P tableaus) and a rigor thread (Van der Pol ε-limit, IVP
test-set problems, diurnal photolysis, conservation audit, structural
adversarial tests, mechanism scale-up), each on its own dev branch off
`dae-rosenbrock-benchmark`.

## Branch status

| Phase | Branch | Status |
|---|---|---|
| 1 — step-size persistence | `dae-step-persistence` | done (`c9ba020f`): 9000→1008 segmented steps |
| 2 — work–precision rig | `dae-work-precision` | done (`ab9231d7`): external Radau refs, ODE tracks rtol, DAE floors at model error |
| 3 — constraint evaluation cost | `dae-constraint-cost` | done (`cdf71882`): mass-coupling inlined; DAE per-step at parity, faster than ODE on both chemistry sweeps |
| 4 — Schur reduction | `dae-schur-reduction`, `dae-schur-core` | done: `SchurStageSolver` implemented in full (exact, order-preserving, cached, tested) + fill-guard early-out; both measured regimes unprofitable — kept as guarded opt-in, see design note follow-ons |
| 5 — norm policy + diagnostics | `dae-norms-diagnostics` | done (`73b9acda`): cellwise-max WRMS batch-invariant; pivot ratio tracks conditioning |
| 6 — RODAS-P tableaus | `dae-rodas-p` | done (`e9753b66`): RODAS4P holds stiff order 3.03 where RODAS4 drops to 1.02 |
| 7a — Van der Pol ε-sweep | `dae-vdp-epsilon` | done (`e0f076a5`): uniform first-order ODE→DAE convergence |
| 7b — IVP test-set problems | `dae-ivp-testset` | done (`1a1c0e96`): conservation-Robertson t=1e11, Akzo Nobel DAE, HIRES — all match published references |
| 7c — diurnal photolysis | `dae-diurnal` | done (`cff52dbb`): guarded recipe (ODE at night); twilight error envelope 0.87, day 3.4e-2 |
| 7d — conservation audit | `dae-conservation` | done (`aa996a48`): NOy conserved to ULP by both formulations |
| 7e — structural/adversarial | `dae-structural` | done (`a33ff0a4`): index-2 clean failure; 200-mechanism random sweep green |
| 7f — mechanism scale-up | `dae-mechanism-scale` | done (`f02f4692`): cell scaling linear, ratio batch-invariant |
| 8 — paper | (ODE repo `paper_dae/`) | done — v3 (`35427af`) current through batching + TS1 |
| follow-on — batched constraint evaluation | `dae-constraint-batch` | done (`775426f3`): type-packed loops; equilibrium family at ODE parity (N=256 ratio 2.59→1.28) |
| follow-on — TS1 mechanism scale-up | `dae-constraint-batch` | done (`34deacbc`): real MOZART-TS1 imported (230 species/501 reactions); no step advantage at 4/9/17 radicals — boundary-of-applicability result; Schur fill guard added (overhead 2.9×→6%) |
