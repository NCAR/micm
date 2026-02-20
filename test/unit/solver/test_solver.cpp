#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

namespace
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b } };

  micm::LambdaRateConstant makeLambdaRateConstant()
  {
    return micm::LambdaRateConstant(
        micm::LambdaRateConstantParameters{ .label_ = "lambda_rc",
                                            .lambda_function_ = [](const micm::Conditions& conditions)
                                            { return 1.0e-3 * conditions.temperature_; } });
  }
}  // namespace

TEST(Solver, GetRateConstantByNameCanOverrideLambda)
{
  auto lambda_rate_constant = makeLambdaRateConstant();
  micm::Process reaction = micm::ChemicalReactionBuilder()
                               .SetReactants({ a })
                               .SetProducts({ micm::StoichSpecies(b, 1) })
                               .SetRateConstant(lambda_rate_constant)
                               .SetPhase(gas_phase)
                               .Build();

  auto system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });

  auto solver = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
                    .SetSystem(system)
                    .SetReactions({ reaction })
                    .Build();

  auto& lambda_ref = solver.GetRateConstantByName("lambda_rc");
  lambda_ref.parameters_.lambda_function_ = [](const micm::Conditions& conditions)
  { return 2.0e-3 * conditions.temperature_; };

  auto state = solver.GetState(1);
  state.SetAbsoluteTolerances({ 1.0e-6, 1.0e-6 });
  state.variables_[0] = { 1.0, 0.0 };
  state.conditions_[0].temperature_ = 300.0;
  state.conditions_[0].pressure_ = 101253.3;
  state.conditions_[0].air_density_ = 1.0e6;

  solver.CalculateRateConstants(state);
  EXPECT_NEAR(state.rate_constants_[0][0], 0.6, 1.0e-12);

  EXPECT_NO_THROW({
    auto result = solver.Solve(1.0, state);
    EXPECT_EQ(result.stats_.final_time_, 1.0);
  });
}
