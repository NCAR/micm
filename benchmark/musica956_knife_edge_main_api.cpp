// musica#956 Case-2 knife-edge reproducer, ported to the templated constraint
// API present on NCAR/micm main and on PR #1083's branch.
// Mechanism is unchanged from benchmark/musica956_knife_edge.cpp:
//   Henry's law   K1(T)*[GAS] = [AQ]      (Van't Hoff, exothermic)
//   Dissociation  K2*[AQ]     = [XM]^2    (quadratic in the solved ion)
#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <cstdio>
#include <utility>
#include <vector>

using namespace micm;

using StandardBuilder =
    CpuSolverBuilder<RosenbrockSolverParameters, Matrix<micm::Real>, SparseMatrix<micm::Real, SparseMatrixStandardOrdering>>;
using DenseMatrix = StandardBuilder::DenseMatrixPolicyType;
using SparseMatrixT = StandardBuilder::SparseMatrixPolicyType;

int main()
{
  const auto gas = Species("GAS");
  const auto aq = Species("AQ");
  const auto hp = Species("HP");
  const auto xm = Species("XM");
  const auto sink = Species("SINK");
  const Phase phase{ "cloud", std::vector<PhaseSpecies>{ gas, aq, hp, xm, sink } };

  Process loss = ChemicalReactionBuilder()
                     .SetReactants({ gas })
                     .SetProducts({ StoichSpecies(sink, 1) })
                     .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0e-12 })
                     .SetPhase(phase)
                     .Build();

  auto run = [&](double temperature)
  {
    std::vector<Constraint<DenseMatrix, SparseMatrixT>> constraints;
    constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
        "henry",
        aq,
        std::vector<StoichSpecies>{ { gas, 1.0 } },
        std::vector<StoichSpecies>{ { aq, 1.0 } },
        { .K_HLC_ref_ = 0.5, .delta_H_ = -60000.0 }));
    constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
        "dissoc",
        xm,
        std::vector<StoichSpecies>{ { aq, 1.0 } },
        std::vector<StoichSpecies>{ { xm, 2.0 } },
        { .K_HLC_ref_ = 1.0, .delta_H_ = 0.0 }));

    auto params = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = StandardBuilder(params)
                      .SetSystem(System(phase))
                      .SetReactions({ loss })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("GAS")] = 4.0e6;
    state.variables_[0][state.variable_map_.at("AQ")] = 1.0;
    state.variables_[0][state.variable_map_.at("HP")] = 1.0e3;
    state.variables_[0][state.variable_map_.at("XM")] = 900.0;
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
    std::printf("%5.1f   %-32s   %.6e\n", T, SolverStateToString(status).c_str(), aq_value);
  }
  return 0;
}
