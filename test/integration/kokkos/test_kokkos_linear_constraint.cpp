// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../precision_matchers.hpp"

#include <micm/Kokkos.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <type_traits>

using namespace micm;
using DenseMatrix = KokkosDenseMatrix<Real, MICM_KOKKOS_DEFAULT_VECTOR_SIZE>;
using StdSparseMatrix = KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<MICM_KOKKOS_DEFAULT_VECTOR_SIZE>>;

TEST(DAESolveWithConstraint, TerminatorAndRobertson)
{
  auto Cl2 = micm::Species("Cl2");
  auto Cl = micm::Species("Cl");
  Cl2.SetProperty("absolute tolerance", micm::Real{ 1.0e-20 });
  Cl.SetProperty("absolute tolerance", micm::Real{ 1.0e-20 });

  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", { Cl2, Cl, A, B, C } };

  micm::Process terminator_r1 = micm::ChemicalReactionBuilder()
                                    .SetReactants({ Cl2 })
                                    .SetProducts({ micm::StoichSpecies(Cl, 2.0) })
                                    .SetPhase(gas_phase)
                                    .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "terminator_k1" })
                                    .Build();

  micm::Process terminator_r2 = micm::ChemicalReactionBuilder()
                                    .SetReactants({ Cl, Cl })
                                    .SetProducts({ micm::StoichSpecies(Cl2, 1.0) })
                                    .SetPhase(gas_phase)
                                    .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 1.0 })
                                    .Build();

  micm::Process robertson_r1 = micm::ChemicalReactionBuilder()
                                   .SetReactants({ A })
                                   .SetProducts({ micm::StoichSpecies(B, 1) })
                                   .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "robertson_r1" })
                                   .SetPhase(gas_phase)
                                   .Build();

  micm::Process robertson_r2 = micm::ChemicalReactionBuilder()
                                   .SetReactants({ B, B })
                                   .SetProducts({ micm::StoichSpecies(B, 1), micm::StoichSpecies(C, 1) })
                                   .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "robertson_r2" })
                                   .SetPhase(gas_phase)
                                   .Build();

  micm::Process robertson_r3 = micm::ChemicalReactionBuilder()
                                   .SetReactants({ B, C })
                                   .SetProducts({ micm::StoichSpecies(A, 1), micm::StoichSpecies(C, 1) })
                                   .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "robertson_r3" })
                                   .SetPhase(gas_phase)
                                   .Build();

  std::vector<micm::Process> processes{ terminator_r1, terminator_r2, robertson_r1, robertson_r2, robertson_r3 };

  // ---------------------------------------------------------------------------
  // Constraint: A + B + C = 1
  // ---------------------------------------------------------------------------

  micm::Real sum_initial_conc = 1.0;

  std::vector<Constraint<DenseMatrix, StdSparseMatrix>> constraints;
  constraints.emplace_back(LinearConstraint<DenseMatrix, StdSparseMatrix>(
      "mass_conservation", C, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, sum_initial_conc));

  // ---------------------------------------------------------------------------
  // Solver
  // ---------------------------------------------------------------------------

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

  auto solver = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(gas_phase))
                    .SetReactions(processes)
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  // 1e-8 is below float epsilon -- see the same override in terminator.hpp.
  state.SetRelativeTolerance(micm::Real{ std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-5 });

  // Robertson rates
  state.SetCustomRateParameter("robertson_r1", 0.04);
  state.SetCustomRateParameter("robertson_r2", 3.0e7);
  state.SetCustomRateParameter("robertson_r3", 1.0e4);

  // Initial conditions
  state[Cl] = 1.2e-6;
  state[Cl2] = 1.8e-10;
  state[A] = 1.0;
  state[B] = 0.0;
  state[C] = 0.0;

  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101300.0;
  state.conditions_[0].air_density_ = 42.0;

  state.variables_.CopyToDevice();
  state.conditions_.CopyToDevice();

  solver.UpdateStateParameters(state);

  state.custom_rate_parameters_.CopyToHost();

  constexpr micm::Index N = 12;
  micm::Real time_step = 1.0;

  for (micm::Index i = 0; i < N; ++i)
  {
    micm::Real advanced = 0.0;

    while (advanced < time_step)
    {
      auto result = solver.Solve(time_step - advanced, state);
      advanced += result.stats_.final_time_;
    }

    state.variables_.CopyToHost();

    // 1. Mass conservation enforced by DAE constraint
    EXPECT_REAL_CLOSE(state[A] + state[B] + state[C], sum_initial_conc, 1e-10);

    time_step *= 10.0;
  }
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
