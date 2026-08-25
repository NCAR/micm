// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "../precision_matchers.hpp"

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>

/// @brief Integrate the Terminator and Robertson mechanisms side by side under a linear
///        mass-conservation constraint on the Robertson species.
///
/// The two sub-systems share a phase but not a reaction, so each one's invariants must hold
/// independently:
///
///   Terminator:  Cl2 --terminator_k1--> 2 Cl,  Cl + Cl --k=1--> Cl2
///                total chlorine 2*[Cl2] + [Cl] is invariant under both reactions.
///                terminator_k1 is deliberately left at its default of zero, so this
///                sub-system only ever runs forward, Cl -> Cl2.
///   Robertson:   A -> B (0.04), B + B -> B + C (3e7), B + C -> A + C (1e4)
///                driven entirely by user-defined rate constants, and closed by the DAE
///                constraint A + B + C = 1 with C as the algebraic variable.
///
/// The Robertson rates are host inputs, so they have to be pushed to the device *before*
/// UpdateStateParameters() computes rate constants from them. Getting that wrong is not a
/// crash: the device simply reads a freshly allocated, zero-filled rate array, no Robertson
/// reaction fires, and A + B + C = 1 stays trivially true at its initial values forever.
/// That is why this test asserts the trajectory actually moved, not just that the constraint
/// held -- see the EvolutionCheck assertions below.
template<class BuilderPolicy>
void TestTerminatorAndRobertson(BuilderPolicy builder)
{
  using DenseMatrix = typename BuilderPolicy::DenseMatrixPolicyType;
  using StdSparseMatrix = typename BuilderPolicy::SparseMatrixPolicyType;

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

  std::vector<micm::Constraint<DenseMatrix, StdSparseMatrix>> constraints;
  constraints.emplace_back(micm::LinearConstraint<DenseMatrix, StdSparseMatrix>(
      "mass_conservation", C, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, sum_initial_conc));

  // ---------------------------------------------------------------------------
  // Solver
  // ---------------------------------------------------------------------------

  auto solver = builder.SetSystem(micm::System(gas_phase))
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

  // Total chlorine, conserved by both terminator reactions: Cl2 -> 2 Cl adds two Cl for
  // every Cl2 removed, and Cl + Cl -> Cl2 does the reverse.
  const micm::Real total_chlorine_initial = 2.0 * micm::Real{ 1.8e-10 } + micm::Real{ 1.2e-6 };

  state.variables_.CopyToDevice();
  state.conditions_.CopyToDevice();
  // The user-defined Robertson rates must reach the device before rate constants are
  // computed from them.
  state.custom_rate_parameters_.CopyToDevice();

  solver.UpdateStateParameters(state);

  constexpr micm::Index N = 12;
  micm::Real time_step = 1.0;

  for (micm::Index i = 0; i < N; ++i)
  {
    micm::Real advanced = 0.0;

    while (advanced < time_step)
    {
      auto result = solver.Solve(time_step - advanced, state);
      // A solve that ends short returns only part of the interval it was given, so once what it
      // returns falls under half an ulp of `advanced` the accumulator stops moving and this loop
      // spins forever. See the same guard in TestAnalyticalOregonator's sub-step loop.
      if (advanced + result.stats_.final_time_ == advanced)
      {
        break;
      }
      advanced += result.stats_.final_time_;
    }

    state.variables_.CopyToHost();

    // VariableProxy converts to micm::Real but has no operator<<, so pull the values into
    // named locals: gtest then reports the number rather than a byte dump on failure.
    const micm::Real cl2 = state[Cl2];
    const micm::Real cl = state[Cl];
    const micm::Real a = state[A];
    const micm::Real b = state[B];
    const micm::Real c = state[C];

    // 1. Mass conservation enforced by DAE constraint
    EXPECT_REAL_CLOSE(a + b + c, sum_initial_conc, 1e-10);

    // 2. Total chlorine is conserved by the kinetics rather than by a constraint, so it is
    //    only held to the accuracy of the integration.
    EXPECT_REAL_SOLVE_REL(2.0 * cl2 + cl, total_chlorine_initial, 1.0e-6)
        << "chlorine not conserved after step " << i;

    // 3. Nothing may go negative or non-finite.
    for (const auto& [name, value] :
         { std::pair{ "Cl2", cl2 }, std::pair{ "Cl", cl }, std::pair{ "A", a }, std::pair{ "B", b }, std::pair{ "C", c } })
    {
      EXPECT_TRUE(std::isfinite(static_cast<double>(value))) << name << " not finite at step " << i;
      EXPECT_GE(value, micm::Real{ 0.0 }) << name << " went negative at step " << i;
    }

    // 4. EvolutionCheck: the Robertson sub-system is driven only by the user-defined rate
    //    constants, so these are the assertions that fail if those rates never reach the
    //    device. After the first second A has been measurably consumed and C produced
    //    (the reference solution is A ~ 0.966, B ~ 3.1e-5, C ~ 0.033 at t = 1 s); the bounds
    //    are one-sided and loose so they hold in single precision and under any reasonable
    //    step-size history, but they are nowhere near satisfied by A = 1, B = C = 0.
    if (i == 0)
    {
      EXPECT_LT(a, micm::Real{ 0.99 }) << "A was not consumed -- did the rate constants reach the device?";
      EXPECT_GT(b, micm::Real{ 0.0 }) << "B was never produced";
      EXPECT_GT(c, micm::Real{ 1.0e-3 }) << "C was never produced";
    }

    time_step *= 10.0;
  }

  // 5. By the end of the run (t ~ 1e11 s) Robertson has run to completion, C -> 1, and the
  //    terminator has converted essentially all Cl into Cl2.
  const micm::Real a_final = state[A];
  const micm::Real c_final = state[C];
  const micm::Real cl_final = state[Cl];
  const micm::Real cl2_final = state[Cl2];
  EXPECT_GT(c_final, micm::Real{ 0.95 }) << "Robertson did not run to completion";
  EXPECT_LT(a_final, micm::Real{ 0.05 });
  EXPECT_GT(cl2_final, micm::Real{ 0.4 } * total_chlorine_initial) << "Cl was not converted to Cl2";
  EXPECT_LT(cl_final, micm::Real{ 0.1 } * total_chlorine_initial);
}
