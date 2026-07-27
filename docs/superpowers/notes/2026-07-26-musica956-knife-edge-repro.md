# musica#956 Case-2 reproduction: the t = 0 temperature knife edge

- **Date:** 2026-07-26
- **Program:** `benchmark/musica956_knife_edge.cpp` (standalone, header-only
  link against any MICM checkout; not wired into CMake — compile directly:
  `clang++ -std=c++20 -O1 -I <micm>/include musica956_knife_edge.cpp`)
- **Context:** NCAR/musica#956 reports cloud-chemistry DAE failures at t = 0
  with ~1 K sensitivity (Case 1: 280 K ok; Case 2: 286 K fails; Case 3:
  285 K ok). The issue ships no config or mechanism, so this is a
  mechanism-faithful reproduction, not an exact replication: a Henry's-law
  dissolution (Van't Hoff K1(T)) feeding an aqueous dissociation quadratic
  in the solved ion (K2·[AQ] = [XM]²), both as built-in
  `EquilibriumConstraint`s, one slow kinetic process, T swept 278–292 K.

## Result matrix (identical program, three MICM trees)

| Tree | 278–292 K sweep |
|---|---|
| `main` (`72d41e16`) | erratic 1 K flips: fails at 284, 287, 290 K; converges at neighbors |
| `fix-rejection-alpha` (`a826c8d0`) | identical flips — the alpha fix does not touch t = 0 |
| `init-weighted-correction` (`03dc6693`) | converges at all 15 temperatures, same manifold values |

The `AQ after init` column is identical across all three trees at every
converging temperature — and, tellingly, at the *failing* ones on `main`:
the projection reaches the manifold and reports failure anyway (the
842-class false failure of the 5400-case sweep, observed live).

## Mechanism, refined by the experiment

The raw-residual acceptance rule (`max_residual < 1e-10`, absolute) compares
a quantity of magnitude ~K_eq(T)·[AQ] against a fixed tolerance. For a
constraint quadratic in the solved variable, the best representable root
leaves a residual of order ulp(K2·[AQ])/2. Two regimes:

- Scale well below ~9e5: floor < 1e-10, always passes. Scale well above:
  floor > 1e-10, always fails (with the state actually converged).
- **In the crossing band, pass/fail depends on where the nearest
  representable root lands relative to the tolerance — effectively
  pseudo-random under 1 K parameter changes.** This matches the field
  report's erratic pattern (280 ok / 286 fail / 285 ok) better than the
  naive monotone-threshold story: Van't Hoff moves the scale ~9%/K here,
  and each temperature re-rolls the ulp dice.

Two adjacent failure modes surfaced while tuning, both plausible Case-2
stories and both fixed by the port:

1. **Cold-start budget exhaustion:** starting the ion at XM = 1 (far below
   its ~600–800 root), main's undamped Newton overshoots to ~[AQ]/2 on the
   first sweep and then halves its way down, exhausting the 10-iteration
   budget at every temperature. The ported damped line search converges
   from the same start.
2. **All-fail with converged state:** at larger gas loading (1e7), the
   floor exceeds the tolerance at every T — main reports failure across the
   board while holding the correct manifold state.

The weighted-correction rule is immune to all three by construction: it
measures the Newton correction in state units against the integrator's own
atol/rtol weights, so constraint-row scale cancels (affine covariance) and
K_eq(T) cannot move the accept/reject decision.

## What this does and does not establish

- Establishes: main's initialization acceptance reproduces the *signature*
  of musica#956 Case 2 (t = 0 failure, ~1 K sensitivity, non-monotone in
  T) on a faithful aqueous-equilibrium system, and the
  `init-weighted-correction` stack removes it; the alpha fix alone does
  not (and was not expected to — its target is the mid-run Case 4).
- Does not establish: that the reporter's exact mechanism fails by this
  route. Confirming that requires their config — request it in the issue,
  along with the returned `SolverState` per failing case.
