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

// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/configure/solver_config.hpp>

namespace fs = std::filesystem;
using namespace micm;

// TODO(jiwon): type check - double vs vector of double 
struct InitialConditions
{
    std::unordered_map<std::string, double> environments;
    std::unordered_map<std::string, std::vector<double>> concentrations;
    std::unordered_map<std::string, std::vector<double>> custom_rate_params;

    bool empty()
    {
        if (!environments.empty() && !concentrations.empty()) return false;
        else return true; 
    }
};

/// @brief Read CSV file and create maps for 
/// @param filename
/// @return unordered map
InitialConditions readCSVToMap(const std::string& filename)
{
    const std::string CONC_PREFIX = "CONC.";
    const std::string ENV_PREFIX = "ENV.";
    const std::string PHOTO_PREFIX = "PHOTO.";
    const std::string EMIS_PREFIX = "EMIS.";
    const std::string LOSS_PREFIX = "LOSS.";
    constexpr int CONC_POS = 5;
    constexpr int ENV_POS = 4;

    InitialConditions dataMap;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return dataMap;
    }

    std::string line;
    std::string key, value_string;
    size_t delimiter_pos;  
    double value;
    while (std::getline(file, line))
    {   
        if (line.empty()) continue;
        
        // Find concentrations 
        else if (line.find(CONC_PREFIX) != std::string::npos)
        {
            delimiter_pos = line.find(",");
            key = line.substr(CONC_POS, delimiter_pos - CONC_POS);
            value_string = line.substr(delimiter_pos + 1);
            value = std::stod(value_string);
            dataMap.concentrations[key].emplace_back(value);
        }
        // Find environment parameters (temp, pressure, density)
        else if (line.find(ENV_PREFIX) != std::string::npos)
        {
            delimiter_pos = line.find(",");
            key = line.substr(ENV_POS, delimiter_pos - ENV_POS);
            value_string = line.substr(delimiter_pos + 1);
            value = std::stod(value_string);
            dataMap.environments[key] = value;
        }
        // Find custom rate constants
        else if (line.find(PHOTO_PREFIX) != std::string::npos 
               || line.find(EMIS_PREFIX) != std::string::npos 
               || line.find(LOSS_PREFIX) != std::string::npos)
        {
            delimiter_pos = line.find(",");
            key = line.substr(0, delimiter_pos);
            value_string = line.substr(delimiter_pos + 1);
            value = std::stod(value_string);
            dataMap.custom_rate_params[key].emplace_back(value);
        }
        else
        {
            std::cerr << "Error parsing value for key. " << line << std::endl;
        }
    }

    return dataMap;
}


int main(const int argc, const char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: performance_test_micm_performance <path> <matrix-ordering>" << std::endl;
        std::cout << std::endl;
        std::cout << "  path                  Path to the folder holding the configuration data" << std::endl;
        std::cout << "  initial condition     csv file contanining initial condtions" << std::endl;
        return 1;
    }

    fs::path config_path {argv[1]};
    std::string initial_conditions_file {argv[2]};

    // Read configure files
    SolverConfig solverConfig;
    ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
    if (status != ConfigParseStatus::Success) 
    {
        std::cout << "Parsing failed" << std::endl; 
        return 1;
    }

    // Read csv file containing initial conditions 
    auto dataMap = readCSVToMap(initial_conditions_file);
    if (dataMap.empty())
    {
        std::cout << "Please verify initial condition file path, file name and/or formatting" << std::endl; 
        return 1;
    }

    SolverParameters solver_params = solverConfig.GetSolverParams();
    auto chemical_system = solver_params.system_;
    auto reactions = solver_params.processes_;

    RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

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
      // std::cout << "solver state: " << StateToString(result.state_) << std::endl;
      state.variables_ = result.result_;
      
    }
    state.PrintState(time_step);


    return 0;
}