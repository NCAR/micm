// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <unordered_map>
#include <string>
#include <vector>

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
InitialConditions readCSVToMap(const std::string& filename);