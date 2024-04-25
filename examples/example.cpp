#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <iomanip>
#include <map>
#include <stdexcept>


// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/configure/solver_config.hpp>
#include "util_example.cpp"

namespace fs = std::filesystem;
using namespace micm;

int main(const int argc, const char *argv[])
{
    if (argc < 3)
    {
        std::cout << "  path                  Path to the folder holding the configuration data" << std::endl;
        std::cout << "  initial condition     csv file contanining initial condtions" << std::endl;
        return 1;
    }

    fs::path config_path {argv[1]};
    std::string initial_conditions_file {argv[2]};

    // Read csv file containing initial conditions 
    auto dataMap = readCSVToMap(initial_conditions_file);
    if (dataMap.empty())
    {
        std::cout << "Error: Failed in reading CSV file. Please verify file path and/or file name of \"" << initial_conditions_file << "\""<< std::endl; 
        return 1;
    }
    else if (dataMap.is_incomplete())
    {
        std::cout << "Error: Parsing csv file is incomplete . Please verify content formatting of \"" << initial_conditions_file << "\""<< std::endl; 
        return 1;
    }

    // Read configure files
    SolverConfig solverConfig;
    solverConfig.ReadAndParse(config_path);

    SolverParameters solver_params = solverConfig.GetSolverParams();

    // add third-body species parameterizaton on air density
    for (auto& species : solver_params.system_.gas_phase_.species_)
      if (species.name_ == "M")
        species.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    for (auto& process : solver_params.processes_)
    {
      for (auto& reactant : process.reactants_)
        if (reactant.name_ == "M")
          reactant.parameterize_ = [](const Conditions& c) { return c.air_density_; };
      for (auto& product : process.products_)
        if (product.first.name_ == "M")
          product.first.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    }

    auto chemical_system = solver_params.system_;
    auto reactions = solver_params.processes_;

    auto params = solver_params.parameters_;
    RosenbrockSolver<> solver{ chemical_system, reactions, params};

    State state = solver.GetState();

    state.conditions_[0].temperature_ = dataMap.environments["temperature"];  // K
    state.conditions_[0].pressure_ = dataMap.environments["pressure"];        // Pa
    state.conditions_[0].air_density_ = dataMap.environments["air_density"];  // mol m-3

    state.SetConcentrations(dataMap.concentrations);

    if (!dataMap.custom_rate_params.empty())
    {   
        for(auto& [rate_constant_key, rate_constant_val] : dataMap.custom_rate_params)
        {
            state.SetCustomRateParameter(rate_constant_key,  rate_constant_val[0]);
        }
    }

    // choose a timestep and print the initial state
    double time_step = dataMap.environments["time_step"];  // s

    state.PrintHeader();
    state.PrintState(0);

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
      if (result.state_ != SolverState::Converged)
      {
        std::cout << "solver failed to converge" << std::endl;
        return 1;
      }
    }
    state.PrintState(time_step);

    return 0;
}