#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

#include <chrono>
#include <iomanip>
#include <iostream>

// Use our namespace so that this example is easier to read
using namespace micm;

int main()
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ std::vector<Species>{ a, b, c } };

  Process r1 = Process::Create()
                   .SetReactants({ a })
                   .SetProducts({ Yields(b, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r1" }))
                   .SetPhase(gas_phase);

  Process r2 = Process::Create()
                   .SetReactants({ b, b })
                   .SetProducts({ Yields(b, 1), Yields(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r2" }))
                   .SetPhase(gas_phase);

  Process r3 = Process::Create()
                   .SetReactants({ b, c })
                   .SetProducts({ Yields(a, 1), Yields(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r3" }))
                   .SetPhase(gas_phase);

  RosenbrockSolver<> solver{ System(SystemParameters{ .gas_phase_ = gas_phase }),
                             std::vector<Process>{ r1, r2, r3 },
                             RosenbrockSolverParameters::ThreeStageRosenbrockParameters(3, false) };

  auto state = solver.GetState();

  // mol m-3
  state.SetConcentration(a, std::vector<double>{ 1, 2, 0.5 });
  state.SetConcentration(b, std::vector<double>(3, 0));
  state.SetConcentration(c, std::vector<double>(3, 0));

  double k1 = 0.04;
  double k2 = 3e7;
  double k3 = 1e4;
  state.SetCustomRateParameter("r1", std::vector<double>(3, k1));
  state.SetCustomRateParameter("r2", std::vector<double>(3, k2));
  state.SetCustomRateParameter("r3", std::vector<double>(3, k3));

  double temperature = 272.5;  // [K]
  double pressure = 101253.3;  // [Pa]
  double air_density = 1e6;    // [mol m-3]

  for (size_t cell = 0; cell < solver.parameters_.number_of_grid_cells_; ++cell)
  {
    state.conditions_[cell].temperature_ = temperature;
    state.conditions_[cell].pressure_ = pressure;
    state.conditions_[cell].air_density_ = air_density;
  }

  // choose a timestep and print the initial state
  double time_step = 200;  // s

  auto result = solver.Solve(time_step, state);
  std::cout << "Solver state: " << StateToString(result.state_) << std::endl;
  std::cout << "accepted: " << result.stats_.accepted << std::endl;
  std::cout << "function_calls: " << result.stats_.function_calls << std::endl;
  std::cout << "jacobian_updates: " << result.stats_.jacobian_updates << std::endl;
  std::cout << "number_of_steps: " << result.stats_.number_of_steps << std::endl;
  std::cout << "accepted: " << result.stats_.accepted << std::endl;
  std::cout << "rejected: " << result.stats_.rejected << std::endl;
  std::cout << "decompositions: " << result.stats_.decompositions << std::endl;
  std::cout << "solves: " << result.stats_.solves << std::endl;
  std::cout << "singular: " << result.stats_.singular << std::endl;
  std::cout << "final simulation time: " << result.final_time_ << std::endl;
}