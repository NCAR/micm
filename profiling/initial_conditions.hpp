#pragma once

#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <micm/solver/rosenbrock.hpp>

void repeat_concentrations_for_each_grid_cell(std::unordered_map<std::string, std::vector<double>> &concentrations, int grid_cells)
{
  for (auto &pair : concentrations)
  {
    pair.second = std::vector<double>(grid_cells, pair.second[0]);
  }
}

struct InitialConditions
{
  std::unordered_map<std::string, double> environments;
  std::unordered_map<std::string, std::vector<double>> concentrations;
  std::unordered_map<std::string, std::vector<double>> custom_rate_params;

  bool empty()
  {
    if (!environments.empty() && !concentrations.empty())
      return false;
    else
      return true;
  }
};

/// @brief Read CSV file and create maps for
/// @param filename
/// @return unordered map
InitialConditions read_initial_conditions(const std::string &filename)
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

  InitialConditions data_map;
  std::ifstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return data_map;
  }

  std::string line;
  std::string key, value_string;
  size_t delimiter_pos;
  double value;
  while (std::getline(file, line))
  {
    if (line.empty())
      continue;

    // Find concentrations
    else if (line.find(CONC_PREFIX) != std::string::npos)
    {
      delimiter_pos = line.find(",");
      key = line.substr(CONC_POS, delimiter_pos - CONC_POS);
      value_string = line.substr(delimiter_pos + 1);
      value = std::stod(value_string);
      data_map.concentrations[key].emplace_back(value);
    }
    // Find environment parameters (temp, pressure, density)
    else if (line.find(ENV_PREFIX) != std::string::npos)
    {
      delimiter_pos = line.find(",");
      key = line.substr(ENV_POS, delimiter_pos - ENV_POS);
      value_string = line.substr(delimiter_pos + 1);
      value = std::stod(value_string);
      data_map.environments[key] = value;
    }
    // Find custom rate constants
    else if (line.find(PHOTO_PREFIX) != std::string::npos || line.find(EMIS_PREFIX) != std::string::npos || line.find(USER_PREFIX) != std::string::npos || line.find(LOSS_PREFIX) != std::string::npos)
    {
      key = line.substr(0, line.find(","));
      value_string = line.substr(line.find(",") + 1);
      value = std::stod(value_string);
      data_map.custom_rate_params[key].emplace_back(value);
    }
    else if (line.find(SURF_PREFIX) != std::string::npos)
    {
      std::size_t first_comma = line.find(",");
      std::size_t second_comma = line.find(",", first_comma + 1);
      key = line.substr(0, first_comma);
      double eff_rad = std::stod(line.substr(first_comma + 1, second_comma - first_comma));
      double num_conc = std::stod(line.substr(second_comma + 1));
      data_map.custom_rate_params[key + ".effective radius [m]"].emplace_back(eff_rad);
      data_map.custom_rate_params[key + ".particle number concentration [# m-3]"].emplace_back(num_conc);
    }
    else
    {
      std::cerr << "Error parsing value for key. " << line << std::endl;
    }
  }

  return data_map;
}

template <template <class> class MatrixType, template <class> class SparseMatrixType>
void SetState(InitialConditions &initial_conditions, const int grid_cells, micm::State<MatrixType, SparseMatrixType> &state)
{
  if (initial_conditions.concentrations.size() != grid_cells)
  {
    repeat_concentrations_for_each_grid_cell(initial_conditions.concentrations, grid_cells);
    repeat_concentrations_for_each_grid_cell(initial_conditions.custom_rate_params, grid_cells);
    for (size_t i = 0; i < grid_cells; ++i)
    {
      state.conditions_[i].pressure_ = initial_conditions.environments["pressure"];
      state.conditions_[i].temperature_ = initial_conditions.environments["temperature"];
      state.conditions_[i].air_density_ = initial_conditions.environments["air_density"];
    }

    for (auto &[rate_key, val] : initial_conditions.custom_rate_params)
    {
      state.SetCustomRateParameter(rate_key, val);
    }

    double p = initial_conditions.environments["pressure"];
    double t = initial_conditions.environments["temperature"];
    double m = initial_conditions.environments["air_density"] * micm::MolesM3ToMoleculesCm3;

    // k_co_oh_jpl19 non-standard algorithm
    {
      double k0 = 6.9e-33 * std::pow(298.0 / t, 2.1);
      double kinf = 1.1e-12 * std::pow(298.0 / t, -1.3);
      double term2 = 1.0 / (1.0 + std::pow(std::log10(k0 * m / kinf), 2));
      double term1 = (kinf * k0 * m / (kinf + k0 * m)) * std::pow(0.6, term2);
      double k_co_oh_jpl19 = term1 + 1.85e-13 * std::exp(-65.0 / t) * (1.0 - term1 / kinf);
      k_co_oh_jpl19 *= micm::MolesM3ToMoleculesCm3;
      std::vector<double> usr_CO_OH(grid_cells, k_co_oh_jpl19);
      state.SetCustomRateParameter("USER.usr_CO_OH", usr_CO_OH);
    }

    // k_dms_oh_jpl15 non-standard algorithm
    {
      double k_dms_oh_jpl15 = 8.2e-39 * std::exp(5376.0 / t) * m * 0.21 /
                            (1.0 + 1.05e-5 * std::exp(3644.0 / t) * 0.21);
      k_dms_oh_jpl15 *= micm::MolesM3ToMoleculesCm3;
      std::vector<double> usr_DMS_OH(grid_cells, k_dms_oh_jpl15);
      state.SetCustomRateParameter("USER.usr_DMS_OH", usr_DMS_OH);
    }
  }
  state.SetConcentrations(initial_conditions.concentrations);
}