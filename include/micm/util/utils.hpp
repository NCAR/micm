// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>
#include <vector>
#include <format>

namespace micm
{
  inline std::string JoinStrings(const std::vector<std::string>& names)
  {
    std::string result;
    for (size_t i = 0; i < names.size(); ++i)
    {
      if (!names[i].empty())
      {
        if (!result.empty())
          result += ".";
        result += names[i];
      }
    }
    return result;
  }

}  // namespace micm