#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

/// @brief Struct of mapping keys in certain type to its initial condition values
struct InitialConditions
{
  std::unordered_map<std::string, double> environments;
  std::unordered_map<std::string, std::vector<double>> concentrations;
  std::unordered_map<std::string, std::vector<double>> custom_rate_params;

  bool empty()
  {
    return environments.empty();
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
    if (line.empty())
      continue;

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
      catch (std::invalid_argument)
      {
        std::cerr << "Parsing Error: Unable to convert string to double for the value of "
                  << "\"" << value_string << "\"" << std::endl;
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
      catch (std::invalid_argument)
      {
        std::cerr << "Parsing Error: Unable to convert string to double for the value of "
                  << "\"" << value_string << "\"" << std::endl;
        dataMap.incomplete_parsing();
        return dataMap;
      }
    }
    // Find custom rate constants that use UserDefinedRateConstant class
    else if (
        line.find(PHOTO_PREFIX) != std::string::npos || line.find(EMIS_PREFIX) != std::string::npos ||
        line.find(LOSS_PREFIX) != std::string::npos || line.find(USER_PREFIX) != std::string::npos)
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
      catch (std::invalid_argument)
      {
        std::cerr << "Parsing Error: Unable to convert string to double for the value of "
                  << "\"" << value_string << "\"" << std::endl;
        dataMap.incomplete_parsing();
        return dataMap;
      }
    }
    else if (line.find(SURF_PREFIX) != std::string::npos)
    {
      auto last_delimiter_pos = line.find_last_of(',');
      auto second_last_delimiter_pos = line.substr(0, last_delimiter_pos - 1).find_last_of(',');
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
      catch (std::invalid_argument)
      {
        std::cerr << "Parsing Error: Unable to convert string to double for the value of "
                  << "\"" << value_string << "\"" << std::endl;
        dataMap.incomplete_parsing();
        return dataMap;
      }
      value_string = line.substr(last_delimiter_pos + 1);
      try
      {
        value = std::stod(value_string);
        dataMap.custom_rate_params[key + ".particle number concentration [# m-3]"].emplace_back(value);
      }
      catch (std::invalid_argument)
      {
        std::cerr << "Parsing Error: Unable to convert string to double for the value of "
                  << "\"" << value_string << "\"" << std::endl;
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