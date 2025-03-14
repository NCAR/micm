// Each rate constant is in its own header file
#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

#include <iomanip>
#include <iostream>
#include <map>

// Use our namespace so that this example is easier to read
using namespace micm;

int main(const int argc, const char* argv[])
{
  SolverConfig solverConfig;

  std::string config_path = "./configs/rate_constants_no_user_defined";
  try
  {
    solverConfig.ReadAndParse(config_path);
  }
  catch (const std::system_error& e)
  {
    std::cerr << "Error reading and parsing config file: " << e.what() << std::endl;
    return 1;
  }

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();

  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  std::unordered_map<std::string, std::vector<double>> intial_concentration = {
    { "A", { 1.0 } },  // mol m-3
    { "B", { 0.0 } },  // mol m-3
    { "C", { 0.0 } },  // mol m-3
    { "D", { 0.0 } },  // mol m-3
    { "E", { 0.0 } },  // mol m-3
    { "F", { 0.0 } },  // mol m-3
    { "G", { 0.0 } },  // mol m-3
  };

  state.SetConcentrations(intial_concentration);

  state.SetCustomRateParameter("SURF.C surface.effective radius [m]", 1e-7);
  state.SetCustomRateParameter("SURF.C surface.particle number concentration [# m-3]", 2.5e6);

  // choose a timestep and print the initial state
  double time_step = 500;  // s

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
    solver.CalculateRateConstants(state);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time += result.final_time_;
    }

    state.PrintState(time_step * (i + 1));
  }

  return 0;
}