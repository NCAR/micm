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
    bool first = true;

    for (const auto& name : names)
    {
        if (name.empty())
            continue;

        if (!first)
            std::format_to(std::back_inserter(result), ".");

        std::format_to(std::back_inserter(result), "{}", name);
        first = false;
    }

    return result;
}

}  // namespace micm