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

namespace fs = std::filesystem;
using namespace micm;

/// @brief Struct of mapping keys in certain type to its initial condition values 
struct InitialConditions
{
    std::unordered_map<std::string, double> environments;
    std::unordered_map<std::string, std::vector<double>> concentrations;
    std::unordered_map<std::string, std::vector<double>> custom_rate_params;

    bool empty()
    {
        if (!environments.empty())    
            return false;
        else 
            return true; 
    }
    
    bool is_incomplete()
    {
        return !status_;
    }

    void incomplete_parsing()
    {   
        status_ = false;
    }

    private:
        bool status_ = true;
};

/// @brief Reads CSV file and creates InitialConditions object 
///        holding initial values for input data
/// @param filename
/// @return InitialCondtions object
InitialConditions readCSVToMap(const std::string& filename)
{
    const std::string CONC_PREFIX = "CONC.";
    const std::string ENV_PREFIX = "ENV.";
    const std::string PHOTO_PREFIX = "PHOTO.";
    const std::string EMIS_PREFIX = "EMIS.";
    const std::string LOSS_PREFIX = "LOSS.";
    const std::string USER_PREFIX = "USER.";
    const std::string SURF_PREFIX = "SURF.";
    constexpr int CONC_POS = 5;
    constexpr int ENV_POS = 4;

    InitialConditions dataMap;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file \"" << filename << "\"" << std::endl;
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
            delimiter_pos = line.find_last_of(',');
            if (delimiter_pos == std::string::npos)
            {
                std::cerr << "Error: Unable to find the delimiter \',' in \"" << line << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
            
            key = line.substr(CONC_POS, delimiter_pos - CONC_POS);
            value_string = line.substr(delimiter_pos + 1);
            
            try 
            {
                value = std::stod(value_string);
                dataMap.concentrations[key].emplace_back(value);
            }
            catch(std::invalid_argument)
            {
                std::cerr << "Parsing Error: Unable to convert string to double for the value of "<< "\""<< value_string << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
        }

        // Find environment parameters (temp, pressure, density)
        else if (line.find(ENV_PREFIX) != std::string::npos)
        {
            delimiter_pos = line.find_last_of(',');
            if (delimiter_pos == std::string::npos)
            {
                std::cerr << "Error: Unable to find the delimiter \',' in \"" << line << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }

            key = line.substr(ENV_POS, delimiter_pos - ENV_POS);
            value_string = line.substr(delimiter_pos + 1);

            try 
            {
                value = std::stod(value_string);
                dataMap.environments[key] = value;
            }
            catch(std::invalid_argument)
            {
                std::cerr << "Parsing Error: Unable to convert string to double for the value of "<< "\""<< value_string << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
        }
        // Find custom rate constants that use UserDefinedRateConstant class
        else if (line.find(PHOTO_PREFIX) != std::string::npos 
               || line.find(EMIS_PREFIX) != std::string::npos 
               || line.find(LOSS_PREFIX) != std::string::npos
               || line.find(USER_PREFIX) != std::string::npos)
        {
            delimiter_pos = line.find_last_of(',');
            if (delimiter_pos == std::string::npos)
            {
                std::cerr << "Error: Unable to find the delimiter ',' in \"" << line << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }

            key = line.substr(0, delimiter_pos);
            value_string = line.substr(delimiter_pos + 1);
            
            try 
            {
                value = std::stod(value_string);
                dataMap.custom_rate_params[key].emplace_back(value);
            }
            catch(std::invalid_argument)
            {
                std::cerr << "Parsing Error: Unable to convert string to double for the value of "<< "\""<< value_string << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
        }
        else if (line.find(SURF_PREFIX) != std::string::npos)
        {
            auto last_delimiter_pos = line.find_last_of(',');
            auto second_last_delimiter_pos = line.substr(0,last_delimiter_pos-1).find_last_of(',');
            if (last_delimiter_pos == std::string::npos || second_last_delimiter_pos == std::string::npos)
            {
                std::cerr << "Error: Unable to find both delimiters ',' in \"" << line << "\"" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
            key = line.substr(0, second_last_delimiter_pos);
            value_string = line.substr(second_last_delimiter_pos + 1, last_delimiter_pos - second_last_delimiter_pos);
            try
            {
                value = std::stod(value_string);
                dataMap.custom_rate_params[key + ".effective radius [m]"].emplace_back(value);
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << ": Parsing Error: Unable to convert string to double for the value of '" << value_string << "'" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }
            value_string = line.substr(last_delimiter_pos + 1);
            try
            {
                value = std::stod(value_string);
                dataMap.custom_rate_params[key + ".particle number concentration [# m-3]"].emplace_back(value);
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << ": Parsing Error: Unable to convert string to double for the value of '" << value_string << "'" << std::endl;
                dataMap.incomplete_parsing();
                return dataMap;
            }            
        }
        else
        {
            std::cerr << "Error: Unable to parse string \"" << line << "\"" << std::endl;
            dataMap.incomplete_parsing();
            return dataMap;
        }
    }

    return dataMap;
}


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

    // Read configure files
    SolverConfig solverConfig;
    ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
    if (status != ConfigParseStatus::Success) 
    {
        std::cout << "Error: Parsing configuration data failed" << std::endl; 
        return 1;
    }

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

    auto params = RosenbrockSolverParameters::three_stage_rosenbrock_parameters();
    params.relative_tolerance_ = 0.1;
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