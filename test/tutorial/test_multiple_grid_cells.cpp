#include <iomanip>
#include <iostream>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

int main()
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b, b })
                         .products({ yields(b, 1), yields(c, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b, c })
                         .products({ yields(a, 1), yields(c, 1) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .phase(gas_phase);

  micm::RosenbrockSolver<> solver{ micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
                                   std::vector<micm::Process>{ r1, r2, r3 },
                                   micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(3, false) };

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

  state.PrintHeader();
  state.PrintState(0);

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }
    state.PrintState(time_step * (i + 1));
  }
}