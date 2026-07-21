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
| 1 — step-size persistence | `dae-step-persistence` | implemented on branch (`c9ba020f`) |
| 2 — work–precision rig | `dae-work-precision` | implemented on branch (`ab9231d7`) |
| 3 — constraint evaluation cost | `dae-constraint-cost` | planned |
| 4 — Schur reduction | `dae-schur-reduction` | planned |
| 5 — norm policy + diagnostics | `dae-norms-diagnostics` | planned |
| 6 — RODAS-P tableaus | `dae-rodas-p` | planned |
| 7a — Van der Pol ε-sweep | `dae-vdp-epsilon` | implemented on branch (`e0f076a5`) |
| 7b — IVP test-set problems | `dae-ivp-testset` | planned |
| 7c — diurnal photolysis | `dae-diurnal` | planned |
| 7d — conservation audit | `dae-conservation` | planned |
| 7e — structural/adversarial | `dae-structural` | planned |
| 7f — mechanism scale-up | `dae-mechanism-scale` | planned |
| 8 — paper v2 | (DAE analysis repo) | planned |
