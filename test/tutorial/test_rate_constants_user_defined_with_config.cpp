#include <iomanip>
#include <iostream>
#include <map>

// Each rate constant is in its own header file
#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

int main(const int argc, const char* argv[])
{
  SolverConfig solverConfig;

  std::string config_path = "./configs/rate_constants_user_defined";
  ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa

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

  // choose a timestep a print the initial state
  double time_step = 500;  // s

  state.PrintHeader();
  state.PrintState(0);

  double photo_rate = 1e-10;
  double emission_rate = 1e-20;
  double loss = emission_rate * 1e-3;
  // these rates are constant through the simulation
  state.SetCustomRateParameter("EMIS.my emission rate", emission_rate);
  state.SetCustomRateParameter("LOSS.my loss rate", loss);

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;
    // This rate is updated at each time step and would typically vary with time
    state.SetCustomRateParameter("PHOTO.my photolysis rate", photo_rate);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }

    state.PrintState(time_step * (i + 1));
    photo_rate *= 1.5;
  }

  return 0;
}