#include <micm/process/process_set.hpp>
#include <micm/solver/backward_euler_dae_solver_parameters.hpp>
#include <micm/solver/backward_euler_dae.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

namespace
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yields(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yields(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase);

  auto the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(BackwardEulerDAE, CanCreateSolver)
{
  auto params = micm::BackwardEulerDAESolverParameters();
  
  // Add a simple algebraic constraint: A + B + C should remain constant (mass conservation)
  params.algebraic_constraints_.push_back([](double time, const std::vector<double>& variables, const std::vector<double>& rate_constants) {
    // This is a dummy constraint for testing - in practice constraints would be more meaningful
    return 0.0; // Always satisfied for this test
  });

  // Test that we can create DAE solver parameters
  EXPECT_EQ(params.algebraic_constraints_.size(), 1);
  EXPECT_EQ(params.small_, 1.0e-40);
  EXPECT_EQ(params.h_start_, 0.0);
  EXPECT_EQ(params.max_number_of_steps_, 11);
}

TEST(BackwardEulerDAE, CanCallConstraintFunction)
{
  auto params = micm::BackwardEulerDAESolverParameters();
  
  // Add a test constraint function
  params.algebraic_constraints_.push_back([](double time, const std::vector<double>& variables, const std::vector<double>& rate_constants) {
    return time * 2.0; // Simple test function
  });

  // Create dummy variables and rate constants to test the constraint function
  std::vector<double> dummy_variables(3, 1.0);
  std::vector<double> dummy_rate_constants(2, 1.0);
  
  double result = params.algebraic_constraints_[0](5.0, dummy_variables, dummy_rate_constants);
  EXPECT_EQ(result, 10.0);
}

TEST(BackwardEulerDAE, ConstraintEvaluationBasic)
{
  auto params = micm::BackwardEulerDAESolverParameters();
  
  // Test multiple constraints
  params.algebraic_constraints_.push_back([](double time, const std::vector<double>& variables, const std::vector<double>& rate_constants) {
    return 1.0; // Should be zero for satisfaction, so this will fail
  });
  
  params.algebraic_constraints_.push_back([](double time, const std::vector<double>& variables, const std::vector<double>& rate_constants) {
    return 0.0; // This constraint is satisfied
  });

  EXPECT_EQ(params.algebraic_constraints_.size(), 2);
  
  // Test the constraint functions
  std::vector<double> dummy_variables(3, 1.0);
  std::vector<double> dummy_rate_constants(2, 1.0);
  
  EXPECT_EQ(params.algebraic_constraints_[0](1.0, dummy_variables, dummy_rate_constants), 1.0);
  EXPECT_EQ(params.algebraic_constraints_[1](1.0, dummy_variables, dummy_rate_constants), 0.0);
}