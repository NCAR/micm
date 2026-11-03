// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <random>
#include <sstream>
#include <string>

namespace micm
{
  /// @brief Generates a version 4 random Universally Unique Identifier
  ///        Used to build parameter names that are guaranteed unique across models
  /// @return UUID string in the form xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx, where y is one of 8, 9, a, b
  inline std::string GenerateUuid()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 15);
    std::uniform_int_distribution<> dis2(8, 11);

    std::stringstream ss;
    ss << std::hex;
    for (int i = 0; i < 8; i++)
    {
      ss << dis(gen);
    }
    ss << "-";
    for (int i = 0; i < 4; i++)
    {
      ss << dis(gen);
    }
    ss << "-4";  // UUID version 4
    for (int i = 0; i < 3; i++)
    {
      ss << dis(gen);
    }
    ss << "-";
    ss << dis2(gen);  // Variant bits
    for (int i = 0; i < 3; i++)
    {
      ss << dis(gen);
    }
    ss << "-";
    for (int i = 0; i < 12; i++)
    {
      ss << dis(gen);
    }
    return ss.str();
  }
}  // namespace micm
