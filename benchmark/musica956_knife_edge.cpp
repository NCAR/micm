// Mechanism-faithful reproducer for NCAR/musica#956 Case 2 (t = 0 failure
// with ~1 K sensitivity): a Henry's-law dissolution feeding an aqueous
// dissociation, both as MICM EquilibriumConstraints with a Van't Hoff
// K_eq(T). On main's raw-residual initialization rule, the converged
// residual sits at a floating-point cancellation floor proportional to
// K_eq(T)*[GAS]; Van't Hoff moves that floor exponentially in T across the
// fixed absolute tolerance (1e-10), so initialization flips from Converged
// to failed within a ~1 K window. The weighted-correction rule measures the
// Newton correction in state units against the integrator's own tolerances,
// which removes the scale dependence entirely.
//
// Build against different MICM checkouts to compare:
//   clang++ -std=c++20 -O1 <sysroot flags> -I <micm>/include knife_edge.cpp
#include <micm/CPU.hpp>

#include <cstdio>
#include <utility>
#include <vector>

int main()
{
  const auto gas = micm::Species("GAS");
  const auto aq = micm::Species("AQ");
  const auto hp = micm::Species("HP");
  const auto xm = micm::Species("XM");
  const auto sink = micm::Species("SINK");
  const micm::Phase phase{ "cloud", std::vector<micm::PhaseSpecies>{ gas, aq, hp, xm, sink } };

  // One slow kinetic process so the system has a nontrivial ODE part.
  micm::Process loss = micm::ChemicalReactionBuilder()
                           .SetReactants({ gas })
                           .SetProducts({ micm::StoichSpecies(sink, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1.0e-12 })
                           .SetPhase(phase)
                           .Build();

  auto run = [&](double temperature)
  {
    // Henry's law: K1(T)*[GAS] = [AQ]; exothermic in MICM's sign convention
    // so K1 grows with T (warmer -> larger residual scale -> Case-2 analog).
    // Dissociation: K2*[AQ] = [HP]*[XM], T-independent.
    std::vector<micm::Constraint> constraints;
    constraints.emplace_back(micm::EquilibriumConstraint(
        "henry",
        aq,
        std::vector<micm::StoichSpecies>{ { gas, 1.0 } },
        std::vector<micm::StoichSpecies>{ { aq, 1.0 } },
        micm::VantHoffParam{ .K_HLC_ref_ = 0.5, .delta_H_ = -60000.0 }));
    // Quadratic in the solved ion (X^2, symmetric dissociation), as in real
    // aqueous chemistry — the converged residual then sits at a genuine
    // cancellation floor of a few ulps of K2*[AQ] instead of cancelling
    // exactly.
    constraints.emplace_back(micm::EquilibriumConstraint(
        "dissoc",
        xm,
        std::vector<micm::StoichSpecies>{ { aq, 1.0 } },
        std::vector<micm::StoichSpecies>{ { xm, 2.0 } },
        micm::VantHoffParam{ .K_HLC_ref_ = 1.0, .delta_H_ = 0.0 }));

    auto params = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params)
                      .SetSystem(micm::System(phase))
                      .SetReactions({ loss })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("GAS")] = 4.0e6;
    state.variables_[0][state.variable_map_.at("AQ")] = 1.0;    // off-manifold start
    state.variables_[0][state.variable_map_.at("HP")] = 1.0e3;  // fixed differential
    state.variables_[0][state.variable_map_.at("XM")] = 900.0;  // near-manifold start: Newton
                                                                // converges in a few sweeps at
                                                                // every T, so only the residual
                                                                // floor vs. tolerance decides
    state.variables_[0][state.variable_map_.at("SINK")] = 0.0;
    state.conditions_[0].temperature_ = temperature;
    state.conditions_[0].pressure_ = 85000.0;
    solver.UpdateStateParameters(state);

    auto result = solver.Solve(1.0e-6, state);
    return std::make_pair(result.state_, state.variables_[0][state.variable_map_.at("AQ")]);
  };

  std::printf("T (K)   solver state                       AQ after init\n");
  for (double T = 278.0; T <= 292.0 + 1e-9; T += 1.0)
  {
    const auto [status, aq_value] = run(T);
    std::printf("%5.1f   %-32s   %.6e\n", T, micm::SolverStateToString(status).c_str(), aq_value);
  }
  return 0;
}
